#############################################################################################################################
##                                                                                                                      
##  VISUALIZE VOLCANO PLOT FOR THE FFPE DATA
##                                                                                                                      
##  Date: 07 MAY 2022                                                                                                                   
##  
##  Author: Moritz Przybilla
##
##                                                                                                                      
############################################################################################################################
# clear workspace
rm(list=ls())
set.seed(16011985) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "RColorBrewer", "ggplot2", "biomaRt", "httr", "data.table", 
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

# read in the FFPE data
shh.mb.matrix <- read.table("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/MB_TP53_FFPE_counts.RPKM.tsv.gz")
matrix <- sapply(shh.mb.matrix[,c(1:ncol(shh.mb.matrix))], round)
matrix <- sapply(shh.mb.matrix[,c(1:ncol(shh.mb.matrix))], as.integer)
rownames(matrix) <- rownames(shh.mb.matrix)

shh.mb.metadata <- read.table("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/MB_TP53_FFPE_ann.tsv")
shh.mb.metadata$Group <- str_split_fixed(shh.mb.metadata$Type.MET, "_",3)[,3]
shh.mb.metadata[shh.mb.metadata$Group != "TP53", "Group"] <- "WT"
shh.mb.metadata$Group <- factor(shh.mb.metadata$Group, levels = c("TP53", "WT"))

# create metadata
shh.phenoData <- new("AnnotatedDataFrame", data = shh.mb.metadata)

# create lfs mb expression set
shh.mb.expSet <- ExpressionSet(as.matrix(shh.mb.matrix), phenoData = shh.phenoData)

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
genes.to.exclude <- unique(primary.markers$gene_name)

# remove primary marker genes from analysis
matrix <- as.data.frame(matrix)
matrix$genes <- rownames(matrix)
matrix <- matrix[matrix$genes %nin% genes.to.exclude, ]
dim(matrix)
matrix$genes <- NULL

# assess differentially expressed genes for chromothriptic and non-chromothriptic tumors accounting for SHH group 
deseq2.obj <- DESeqDataSetFromMatrix(countData = matrix,
                                     colData = shh.mb.metadata[colnames(matrix),],
                                     design= ~ Group)

# perform some filtering to only incorporate well captured genes
keep <- rowSums(counts(deseq2.obj) >= 50) >= 3
deseq2.obj <- deseq2.obj[keep,]
deseq2.obj$Group <- relevel(deseq2.obj$Group, ref = "WT")

# perform analysis
ddsMF <- DESeq(deseq2.obj)

# investigate the results
resMF <- results(ddsMF)
resMF$gene_name <- rownames(resMF)

# make dataframe and order according to log2FC and pvalue
results <- as.data.frame(resMF)
results <- results[order(results$pvalue, decreasing =  F),]

# write table to file
write.table(results, "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/MB_TP53_vs_WT_DESeq2_wo_Immune.txt", col.names = T, row.names = F, quote = F, sep = "\t")

results$color <- "NS"
results[abs(results$log2FoldChange) > 1, "color"] <- "Log2FC"
results[results$pvalue < 0.00001, "color"] <- "p-value"
results[abs(results$log2FoldChange) > 1 & results$pvalue < 0.00001, "color"] <- "Log2FC & p-value"
results$color <- factor(results$color, levels = c("NS", "Log2FC", "p-value", "Log2FC & p-value"))

# calculate the variance for each gene
vsd <- vst(ddsMF, blind=FALSE)
rv <- rowVars(assay(vsd))
ntop=2000

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))

loadings <- as.data.frame(pca$rotation)

p.variance.explained = pca$sdev^2 / sum(pca$sdev^2)

# plot percentage of variance explained for each principal component    
pdf("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/bulkRNA_FFPE_without_immune_barplot_PCA.pdf")
barplot(100*p.variance.explained[1:50], las=2, xlab='', ylab='% Variance Explained')
dev.off()

# merge the results with the gene names
loadings$hgnc_symbol <- rownames(loadings)
write.table(loadings, "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/MB_TP53_vs_WT_DESeq2_PCA_component.txt", col.names = T, row.names = F, quote = F, sep = "\t")

col.values <- c("lightgrey", pal_lancet("lanonc")(3))
options(ggrepel.max.overlaps = Inf)

up.data <- results %>% filter(pvalue < 0.0001 , log2FoldChange >= 1)
# up.data <- up.data[-grep("RP", up.data$gene_name),]
up.data <- up.data[grep("MKI67IP|TSN|CLASP1|IWS1|DDX1|HMGB2|CCNB1|FKBP4|NTRK3|FOXS1|TP53BP1|TERC|KDM6B|ROBO2|MYCN|PTCH2", up.data$gene_name),]
down.data <- results %>% filter(pvalue < 0.0001 , log2FoldChange <= -1)
# down.data <- down.data[-grep("RP", down.data$gene_name),]
down.data <- down.data[grep("MKI67IP|TSN|CLASP1|IWS1|DDX1|HMGB2|CCNB1|FKBP4|NTRK3|FOXS1|TP53BP1|TERC|KDM6B|ROBO2|MYCN|PTCH2", down.data$gene_name),]

results[grep("TBX5|NBL1|HEY1|MKI67P|CLASP1|GLI2|CDKN2A", results$gene_name),]

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
ggsave(paste0("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/bulkRNA_FFPE_without_immune.pdf"), width = 6, height = 6, dpi = 600)





