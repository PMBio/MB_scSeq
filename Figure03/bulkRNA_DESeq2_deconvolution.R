#############################################################################################################################
##                                                                                                                      
##  PERFORM AN ASSESSMENT OF THE SIMILARITY BETWEEN THE TRANSCRIPTOME OF THE PRIMARY TUMOR DATA AND CELL TYPES IN ALDINGER ET AL
##                                                                                                                      
##  Date: 06 OCTOBER 2020                                                                                                                   
##  
##  Author: Moritz Przybilla
##
##                                                                                                                      
############################################################################################################################
# clear workspace
rm(list=ls())
set.seed(16011985) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "optparse", "MuSiC", "RColorBrewer", "ggplot2", "biomaRt", "httr", "data.table", 
                      "dplyr", "Biobase","tidyverse", "ComplexHeatmap", "SummarizedExperiment", "DESeq2", "EnhancedVolcano")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))
httr::set_config(httr::config(ssl_verifypeer=0L))

`%nin%` = Negate(`%in%`)
############################################################################
##                          READ IN THE DATA
############################################################################

# read in the bulk expression data from RNA-seq
lfs.mb.matrix <- read.table("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/MB_TP53_FFPE_counts_full.gn.txt.gz")
colnames(lfs.mb.matrix) <- gsub("X", "ID_", colnames(lfs.mb.matrix))

# read metadata
lfs.metadata <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/metadata_SHH-TP53.txt", header = T)
lfs.metadata <- lfs.metadata[,c(1:9)]
lfs.metadata <- lfs.metadata[-c(6, 12 , 16:nrow(lfs.metadata))]
lfs.metadata <- as.data.frame(lfs.metadata)
rownames(lfs.metadata) <- paste0("ID_", lfs.metadata$RNA_seq)

# order metadata
lfs.metadata <- lfs.metadata[colnames(lfs.mb.matrix),]
lfs.metadata <- lfs.metadata[complete.cases(lfs.metadata),]

# create metadata
lfs.phenoData <- new("AnnotatedDataFrame", data = lfs.metadata)

# remove two samples
lfs.mb.matrix[,c("ID_65390", "ID_105768")] <- NULL

# create lfs mb expression set
lfs.mb.expSet <- ExpressionSet(as.matrix(lfs.mb.matrix), phenoData = lfs.phenoData)

# read in the scanpy anndata object
nuclei.data <- LoadH5Seurat("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/adata_nuclei.h5seurat")

# get the matrix and metadata
pData<- read.csv('/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/Nuclei_metadata_final.csv')
rownames(pData) <- pData$X
pData$X <- NULL
pData$Cell_types <- as.factor(pData$Cell_types)

# add metadata
phenoData <- new("AnnotatedDataFrame", data = pData)

# create nuclei expression dataframe
celltypes.lfs.scrna <- ExpressionSet(as.matrix(nuclei.data@assays$RNA@counts), phenoData = phenoData)
levels(celltypes.lfs.scrna$Cell_types)

#####################################################################################
# READ IN THE DATA
#####################################################################################

# read in the object
deseq2.obj <- readRDS("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/ICGC_MB_SHH_DESEQ_rmbad.rds")

# read in the old gene list from John
diff.genes <- read.table("/omics/groups/OE0540/internal/projects/przybilm/MB_SHH_DEG_BY_CT.tsv", header = T)
diff.genes <- diff.genes[,c("ens_fullid", "ens_id", "gene_name")]

# read in the metadata
metadata <- read.table("/omics/groups/OE0540/internal/projects/przybilm/mb_ssh_groups.txt", header = T)
metadata <- metadata[, c("id", "prediction")]

# tcc content
tcc.content <- read.table("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/MB_CT_STATUS_TCC.tsv", header = T)
tcc.content$ID <- str_split_fixed(tcc.content$ID, "_", 2)[,2]
tcc.content[tcc.content$ID == "LFS_MB1", "ID"] <- "LFS_MB_P"
tcc.content[tcc.content$ID == "LFS_MB2", "ID"] <- "LFS_MB_1R"

# wrangle into shape
metadata <- metadata[order(metadata$prediction),]

# adapt id information 
metadata$id <- str_split_fixed(metadata$id, "_", 2)[,2]
metadata[metadata$id == "72688_diag", "id"] <- "LFS_MB_P"
metadata[metadata$id == "81216_diag", "id"] <- "LFS_MB_1R"

# combine tumor cell content and metadata
tcc.metadata <- merge(metadata, tcc.content, by.x = "id", by.y = "ID")
tcc.metadata <- tcc.metadata[,c("id", "prediction", "TCC")]

# create new metadata
deseq.metadata <- colData(deseq2.obj)
new.deseq.metadata <- merge(deseq.metadata, tcc.metadata, by.x = "rownames", by.y = "id")
rownames(new.deseq.metadata) <- new.deseq.metadata$rownames

# add the type
mcols(new.deseq.metadata) <- DataFrame("type" = c("input", "input", "intermediate", "intermediate", "input", "input"), "description" = "")
colData(deseq2.obj) <- new.deseq.metadata

table(new.deseq.metadata$CT, new.deseq.metadata$prediction)

# convert to factor
deseq2.obj$prediction <- factor(deseq2.obj$prediction, levels = c("MB_SHH_1", "MB_SHH_2", "MB_SHH_3", "MB_SHH_4"))

# perform some filtering to only incorporate well captured genes
keep <- rowSums(counts(deseq2.obj) >= 50) >= 3
deseq2.obj <- deseq2.obj[keep,]

#####################################################################################
# DETERMINE DIFFERENTIAL GENE EXPRESSION
#####################################################################################

# read in differential expressed genes for the nuclei celltypes
nuclei.celltypes <- read.csv("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/Nuclei_Celltypes_DEG_wilcoxon.csv")
nuclei.celltypes$X <- NULL
celltypes <- unique(str_split_fixed(colnames(nuclei.celltypes), "_", 2)[,1])

# load Riemondy et al. cell types
riemondy.celltypes <- read.csv("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/Riemondy_Coarse_Celltypes_DEG_wilcoxon.csv")
riemondy.celltypes$X <- NULL
celltypes.to.include <- unique(str_split_fixed(colnames(riemondy.celltypes), "_", 2)[,1])

# combine data
celltypes <- c(celltypes, celltypes.to.include)
nuclei.celltypes <- cbind(nuclei.celltypes[1:20000, ], riemondy.celltypes[1:20000, ])

# remove the malignant cells
celltypes <- celltypes[-grep("Malignant|malignant", celltypes)]

primary.markers <- list()
i <- 7
for (i in 1:length(celltypes)){
  
  celltype.tmp <- celltypes[i]
  data.tmp <- nuclei.celltypes[,grep(celltype.tmp, colnames(nuclei.celltypes))]
  
  # subset to significant markers
  data.tmp <- data.tmp[data.tmp[,4] < 0.00001, ]
  data.tmp <- data.tmp[order(data.tmp[,3], decreasing = T),]
  genes <- data.tmp[1:100, 1]
  
  primary.markers[[celltype.tmp]] <- genes
}

# 
primary.markers <- as.data.frame(as.vector(unlist(primary.markers)))
colnames(primary.markers) <- "gene_name"

# additional markers
add.markers <- diff.genes[grep("HLA|IGF|IGH|IGM|IGS", diff.genes$gene_name),"gene_name"]
add.markers <- data.frame("gene_name" = add.markers)

# merge
primary.markers <- rbind(primary.markers, add.markers)

# deduplicate
primary.markers <- as.data.frame(primary.markers[!duplicated(primary.markers), ])
colnames(primary.markers) <- "gene_name"

# merge the results with the gene names
primary.markers <- merge(primary.markers, diff.genes, by = "gene_name", all.x = T)
genes.to.exclude <- unique(primary.markers$ens_id)
genes.to.exclude <- genes.to.exclude[!is.na(genes.to.exclude)]

# remove primary marker genes from analysis
matrix <- deseq2.obj@assays@data$counts
matrix <- as.data.frame(matrix)
dim(matrix)

matrix$genes <- str_split_fixed(rownames(matrix), "\\.", 2)[,1]
matrix <- matrix[matrix$genes %nin% genes.to.exclude, ]
dim(matrix)
matrix$genes <- NULL

# assess differentially expressed genes for chromothriptic and non-chromothriptic tumors accounting for SHH group 
deseq2.ct.obj <- DESeqDataSetFromMatrix(countData = matrix,
                                        colData = new.deseq.metadata[colnames(matrix),],
                                        design= ~ CT)


# perform analysis
ddsMF <- DESeq(deseq2.ct.obj)

# investigate the results
resMF <- results(ddsMF)
resMF$ens_id <- str_split_fixed(rownames(resMF), "\\.", 2)[,1]

# merge the results with the gene names
resMF <- merge(DataFrame(resMF), diff.genes, by = "ens_id")

# remove NAs
resMF <- resMF[complete.cases(resMF),]

# make dataframe and order according to log2FC and pvalue
results <- as.data.frame(resMF)
results <- results[order(results$pvalue, decreasing =  F),]

# write table to file
# write.table(results, "//omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/MB_CT_vs_NCT_DESeq2.txt", col.names = T, row.names = F, quote = F, sep = "\t")

results$color <- "NS"
results[abs(results$log2FoldChange) > 1, "color"] <- "Log2FC"
results[results$pvalue < 0.00001, "color"] <- "p-value"
results[abs(results$log2FoldChange) > 1 & results$pvalue < 0.00001, "color"] <- "Log2FC & p-value"
results$color <- factor(results$color, levels = c("NS", "Log2FC", "p-value", "Log2FC & p-value"))

col.values <- c("lightgrey", pal_lancet("lanonc")(3))
options(ggrepel.max.overlaps = Inf)

ggplot(results, aes(x=log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = color), size = 1) +
  geom_vline(xintercept=c(-1, 1), linetype="dotted", color = "black") +
  geom_hline(yintercept = -log10(0.00001), linetype="dotted", color = "black") +
  theme_classic() +
  # geom_smooth(method=lm, color="darkblue") +
  labs(x="Log2FC",
       y="-Log10(p-value)") + 
  xlim(-8, 8) +
  theme(legend.position="bottom") +
  scale_color_manual(values = col.values) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold")) +
  geom_text_repel(data = results %>% filter(pvalue < 0.00001 , log2FoldChange <= -1),
                  aes(label = paste0(gene_name)) ,
                  hjust = -.35,
                  nudge_x = -0.75,
                  direction = "y",
                  fontface = "bold",
                  size = 2) +
  geom_text_repel(data = results %>% filter(pvalue < 0.00001 , log2FoldChange >= 1),
                  aes(label = paste0(gene_name)) ,
                  hjust = -.35,
                  nudge_x = 0.5,
                  direction = "y",
                  fontface = "bold",
                  size = 2)
ggsave(paste0("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/bulkRNA_volcano_without_immune.pdf"), width = 6, height = 6, dpi = 600)

#####################################################################################
# DETERMINE DIFFERENTIAL GENE EXPRESSION WITH TCC AND SUBGROUP
#####################################################################################

# assess differentially expressed genes for chromothriptic and non-chromothriptic tumors accounting for SHH group 
deseq2.ct.tcc.subgroup.obj <- DESeqDataSetFromMatrix(countData = matrix,
                                                     colData = new.deseq.metadata[colnames(matrix),],
                                                     design= ~ prediction + TCC + CT)
ddsMF <- DESeq(deseq2.ct.tcc.subgroup.obj)

# investigate the results
resMF <- results(ddsMF)
resMF$ens_id <- str_split_fixed(rownames(resMF), "\\.", 2)[,1]

# merge the results with the gene names
resMF <- merge(DataFrame(resMF), diff.genes, by = "ens_id")

# remove NAs
resMF <- resMF[complete.cases(resMF),]

# make dataframe and order according to log2FC and pvalue
results <- as.data.frame(resMF)
results <- results[order(results$pvalue, decreasing =  F),]

# write table to file
write.table(results, "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/MB_CT_vs_NCT_DESeq2_TCC_GROUP_wo_Immune.txt", col.names = T, row.names = F, quote = F, sep = "\t")

results$color <- "NS"
results[abs(results$log2FoldChange) > 1, "color"] <- "Log2FC"
results[results$pvalue < 0.00001, "color"] <- "p-value"
results[abs(results$log2FoldChange) > 1 & results$pvalue < 0.00001, "color"] <- "Log2FC & p-value"
results$color <- factor(results$color, levels = c("NS", "Log2FC", "p-value", "Log2FC & p-value"))

col.values <- c("lightgrey", pal_lancet("lanonc")(3))
options(ggrepel.max.overlaps = Inf)

up.data <- results %>% filter(pvalue < 0.00001 , log2FoldChange >= 1)
up.data <- up.data[-grep("RP", up.data$gene_name),]
down.data <- results %>% filter(pvalue < 0.00001 , log2FoldChange <= -1)
down.data <- down.data[-grep("RP", down.data$gene_name),]

ggplot(results, aes(x=log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = color), size = 1) +
  geom_vline(xintercept=c(-1, 1), linetype="dotted", color = "black") +
  geom_hline(yintercept = -log10(0.00001), linetype="dotted", color = "black") +
  theme_classic() +
  # geom_smooth(method=lm, color="darkblue") +
  labs(x="Log2FC",
       y="-Log10(p-value)") + 
  xlim(-12, 12) +
  theme(legend.position="bottom") +
  scale_color_manual(values = col.values) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold")) +
  geom_text_repel(data = down.data,
                  aes(label = paste0(gene_name)) ,
                  hjust = -.35,
                  nudge_x = -0.75,
                  direction = "y",
                  fontface = "bold",
                  size = 2) +
  geom_text_repel(data = up.data,
                  aes(label = paste0(gene_name)) ,
                  hjust = -.35,
                  nudge_x = 0.5,
                  direction = "y",
                  fontface = "bold",
                  size = 2)
ggsave(paste0("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/bulkRNA_multiVariate_TCC_GROUP_volcano_wo_Immune.pdf"), width = 8, height = 6, dpi = 600)
