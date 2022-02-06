#############################################################################################################################
##                                                                                                                      
##  Analysis of individual subclones in scRNA-seq with ECB information from Sequencing 18
##                                                                                                                      
##  Date: 23 July 2020                                                                                                                    
##  
##  Author: Moritz Przybilla and Kasper Karlsson                                                                                                                 
##                                                                                                                      
#############################################################################################################################

# clear workspace
rm(list = ls())
set.seed(14) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("tidyverse", "patchwork", "Seurat", "Matrix", "biomaRt", "viridis",
                      "GO.db", "HTSanalyzeR2", "org.Hs.eg.db", "KEGGREST", "igraph", "tidyr", "stats", "reshape2", "ggplot2",
                      "forcats","heatmap3", "gplots", "ggpubr", "dplyr", "ComplexHeatmap", "ggplotify", "karyoploteR", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                      "regioneR", "Repitools", "mclust")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

###########################################################################
#                               FUNCTIONS
###########################################################################

# set up functions which are used 
`%notin%` <- Negate(`%in%`)

### FUNCTION TO GET GSEA ENRICHMENT FOR SEURAT OUTPUT FROM FIND MARKERS

geneSetAnalysis <- function(dfile) {
  dfile_subs <- dfile ### SUBSET DEGS FILE TO PVAL <= 0.05
  phenotype <- as.vector(as.numeric(dfile_subs$avg_logFC)) ### USE LOG FOLD CHANGE AS PHENOTYPE VECTOR
  names(phenotype) <- rownames(dfile_subs)
  
  ## specify the gene sets type you want to analyze
  
  # HALLMARK
  MSig_H <- MSigDBGeneSets(species = "Hs", collection = "H", subcategory = NULL) # Hallmarks!
  ListGSC <- list(MSig_H=MSig_H)
  
  ## iniate a *GSCA* object
  gsca <- GSCA(listOfGeneSetCollections=ListGSC, 
               geneList=phenotype)
  
  ## preprocess
  gsca1 <- preprocess(gsca, species="Hs", initialIDs="SYMBOL",
                      keepMultipleMappings=TRUE, duplicateRemoverMethod="max",
                      orderAbsValue=FALSE)
  
  
  ## analysis
  if (requireNamespace("doParallel", quietly=TRUE)) {
    doParallel::registerDoParallel(cores=4)
  }  ## support parallel calculation using multiple cores
  
  
  gsca2 <- analyze(gsca1, 
                   para=list(pValueCutoff=0.05, pAdjustMethod="BH",
                             nPermutations=100, minGeneSetSize=1,
                             exponent=1), 
                   doGSOA = FALSE)
  # return(getResult(gsca2)$GSEA.results$MSig_H)
  return(getResult(gsca2)$GSEA.results$MSig_H)
  message("Successfully ran GSEA.")
}

###########################################################################
#                       READ IN THE DATA OF INTEREST
###########################################################################

# define output directory
o.dir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/scRNA_analysis"

# create o.dir
dir.create(o.dir)
setwd(o.dir)

# get sample.ids
sample.tmp <- "LFSMBP-PDX"

# which sample? 
message(sample.tmp)
dir.create(paste0(o.dir, "/", sample.tmp, "/DEG"))

# read in files of interest
marker.files <- list.files(paste0("/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/scRNA_analysis/scanpy"), pattern="scDNA_clones_only", full.names=T, recursive = T)
marker.files <- marker.files[grep("LFSMBP-PDX", marker.files)]
marker.files <- read.table(marker.files[1], sep = ",", stringsAsFactors = F)
colnames(marker.files) <- marker.files[1,]
marker.files <- marker.files[-1,]
colnames(marker.files)[1] <- "index"

# split into distinct clones
clone1.file <- marker.files[,grep("Clone1", colnames(marker.files))]
clone2.file <- marker.files[,grep("Clone2", colnames(marker.files))]
clone3.file <- marker.files[,grep("Clone3", colnames(marker.files))]
clone4.file <- marker.files[,grep("Clone4", colnames(marker.files))]
clone5.file <- marker.files[,grep("Clone5", colnames(marker.files))]
clone6.file <- marker.files[,grep("Clone6", colnames(marker.files))]

# iterate over marker file and select 250 genes
processed.marker.files <- list(clone1.file, clone2.file, clone3.file, clone4.file, clone5.file, clone6.file)
marker.files <- list()
final.marker.files <- list()
final.clone.ids <- c()

# subclones present
subclone.id <- c("clone1", "clone2", "clone3", "clone4", "clone5", "clone6")

i <- 1
for (i in 1:length(processed.marker.files)){
  
  # which clone?
  clone.tmp <- subclone.id[i]
  
  # get the data
  data <- processed.marker.files[[i]]
  
  # change colnames
  colnames(data) <- c("hgnc_symbol", "avg_logFC", "z_score", "p_val")
  rownames(data) <- data$hgnc_symbol
  
  # make numeric 
  data$z_score <- as.numeric(data$z_score)
  data$p_val <- as.numeric(data$p_val)
  
  # save a list of the files
  marker.files[[i]] <- data
  
  if(nrow(data[data$p_val < 0.05,]) < 1){
    
    cat(paste0(clone.tmp, " does not have any significant genes and will thus not be tested.\n"))
    next
    
  } else {
    
    # filter the top250 up and top250 downregulated genes
    up.marker.files <- data[order(data$p_val, decreasing = F),]
    
    # combine
    # new.data <- new.data[new.data$p_val < 0.05,]
    final.marker.files[[clone.tmp]] <- up.marker.files
    final.clone.ids <- c(final.clone.ids, clone.tmp)
    
  }
  
}

###########################################################################
#           GET THE NUMBER OF GENES WHICH IS EXPRESSED PER ARM
###########################################################################

# GET GENE SET ENRICHMENT FOR EACH DEG FILE
GSEA.results <- lapply(final.marker.files, try(geneSetAnalysis))

### PRINT TO FILE
for (i in 1:length(GSEA.results)){
  
  # add the hallmark column to the file
  results <- GSEA.results[[i]]
  results$hallmark <- rownames(results)
  
  # add a subclone column
  results$subclone <- final.clone.ids[i]
  
  # exchange the files
  GSEA.results[[i]] <- results
  
  write.table(results, paste0(sample.tmp, "/DEG/", final.clone.ids[i], "_Hallmark_GSEA_vs_clones.txt"),sep="\t", quote=FALSE, row.names = TRUE, col.names=TRUE)
}

GSEA.results <- list.files(paste0(sample.tmp, "/DEG"), pattern = "_Hallmark_GSEA_vs_clones", recursive = T, full.names = T)
GSEA.results <- lapply(GSEA.results, read.table, header = T)

# combine all the files together
all <- bind_rows(GSEA.results)
all$Leading.Edge <- NULL
all[all$Adjusted.Pvalue < 0.05,]

# replace score according to down or up-regulation
all[all$Observed.score < 0, "Observed.score"] <- -1
all[all$Observed.score > 0, "Observed.score"] <- 1

# replace 0 is in the pvalue column by 0.00000000001
all[all$Adjusted.Pvalue > 0.05, "Adjusted.Pvalue"] <- 1
all[all$Adjusted.Pvalue <= 0.05 & all$Adjusted.Pvalue > 0.001, "Adjusted.Pvalue"] <- 0.05
all[all$Adjusted.Pvalue <= 0.001 & all$Adjusted.Pvalue > 0.00001, "Adjusted.Pvalue"] <- 0.001
all[all$Adjusted.Pvalue == 0, "Adjusted.Pvalue"] <- 0.00001

# reorder file to make it consistent
all <- all[,c(5,4,1:3)]
colnames(all) <- c("sample","hallmark","direction","pVal","adjp")
unique(all$sample)

# make a pvalue column with - and + according to up and down-regulation
all$pVal_dir <- all$adjp*all$direction

# Remove for plot non-essential columns
all[,c("direction","pVal","adjp")] <- NULL

# tranpose with Hallmarks as columns instead
all.wide <- all %>% spread(hallmark, pVal_dir,fill=1) # Go from tall to wide format
rownames(all.wide) <- all.wide$sample # Make samples name as row name

# remove sample column
all.wide$sample <- NULL

# Transpose
all.wide.t <- as.data.frame(t(all.wide))

# For plotting purpose add a column for the RG we're comparing with. Seq8: RG1, Seq18: RG2 (actually RG2c, but here called RG2)
all.wide.t$clone2 <- 1
all.wide.t$clone5 <- 1
all.wide.t$clone6 <- 1

j <- 1
keep <- c()
for (j in 1:nrow(all.wide.t)){
  
  row <- abs(as.numeric(all.wide.t[j,]))
  
  if (any(row < 0.05)){
    
    keep <- c(keep, j)
  } else {
    
    next
  }
  
  
  
}

sub.all.wide.t <- all.wide.t[keep,]

### CONVERT TO MATRIX TO COMPLY WITH HEATMAP PROGRAM
allm <- data.matrix(sub.all.wide.t)

# subset and order allm
allm <- allm[,c("clone1", "clone2", "clone3", "clone4", "clone5", "clone6")]

### ADD IN BREAKS AND COLORS (FOR PVALUE)
colors = c(-1,-0.05,-0.001,-0.00001,0,0.00001,0.001,0.05,1)
my_palette <- c("white","#b3cde0","#005b96","#011f4b","black","#a70000","#ff0000","#ffbaba","white")

subclone.id <- c("clone1", "clone2", "clone3", "clone4", "clone5", "clone6")
marker.tables <- marker.files

### for merging marker files
for (i in 1:length(marker.tables)){
  
  # add the hallmark column to the file
  file <- marker.tables[[i]]
  file$genes <- rownames(file)
  
  # add a subclone column
  file$subclone <- subclone.id[i]
  
  # exchange the files
  marker.tables[[i]] <- file
}

# combine all the files together
marker.file.df <- bind_rows(marker.tables)
marker.file.df <- marker.file.df[marker.file.df$p_val < 0.05,]

# get the number of differential expressed genes per subclone
marker.df <- data.frame(table(marker.file.df$subclone))
rownames(marker.df) <- marker.df$Var1
marker.df$Var1 <- NULL

# manually add subclone 1
marker.df[nrow(marker.df)+1,] <- 0
rownames(marker.df)[nrow(marker.df)] <- "clone2"
marker.df[nrow(marker.df)+1,] <- 0
rownames(marker.df)[nrow(marker.df)] <- "clone5"
marker.df[nrow(marker.df)+1,] <- 0
rownames(marker.df)[nrow(marker.df)] <- "clone6"

# convert to matrix and reorder
marker.matrix <- as.matrix(marker.df)
marker.matrix <- as.matrix(marker.matrix[colnames(allm),])

# generate the first heatmap based on the score matrix
ht1 <- Heatmap(allm, 
               name = "score", 
               col = circlize::colorRamp2(colors, my_palette), 
               cluster_rows = F, 
               cluster_columns = FALSE,
               show_row_names = F, 
               row_names_side = "left",
               show_row_dend = FALSE,
               show_column_names = FALSE,
               row_names_gp = gpar(fontsize = 9), 
               rect_gp = gpar(col = "black", lwd = 1), 
               column_title_gp = gpar(fontsize = 15, fontface = "bold"),
               heatmap_legend_param = list(color_bar = "discrete",
                                           at = colors,
                                           title = ""),
               top_annotation = HeatmapAnnotation(column_barplot = anno_barplot(marker.matrix, 
                                                                                border = TRUE, 
                                                                                height = unit(1.5, "cm"),
                                                                                annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
                                                                                bar_width = 0.75, 
                                                                                gp = gpar(col = "darkgrey", fill = "blue", fontsize = 10, fontface = "bold"), 
                                                                                axis_param = list(at = c(0, 50, 100, 150, 200, 250),
                                                                                                  labels = c("0", "50", "100", "150",  "200", "250")),
                                                                                width = unit(2, "cm")),
                                                  show_annotation_name = F), 
               border = T)

pdf(paste0(sample.tmp, "/DEG/GSEA_heatmap_pval_ordered_hallmarks_scicone.pdf"),width=6,height = 5,pointsize=0.1)
print(ht1)
dev.off()


###########################################################################
#     CHECK THE POSITION OF ALL GENES WHICH ARE DIFFERENTIALLY EXPRESSED
###########################################################################

diff.exp.gene.df <- list()
for (i in 1:length(marker.files)){
  
  file <- marker.files[[i]]
  file$clone <- subclone.id[i]
  diff.exp.gene.df[[i]] <- file
  
}
diff.exp.gene.df <- bind_rows(diff.exp.gene.df)

# get mart object with certain attributes
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
ann <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"), mart  = mart)

# order cns_reads
gene_ordering_file <- ann
colnames(gene_ordering_file) <- c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id")
gene_ordering_file$chr <- paste0("chr", gene_ordering_file$chr)

# set accepted levels and remove other chromosomes
chrOrder<-c(paste("chr",1:22,sep=""))
gene_ordering_file$chr <- factor(gene_ordering_file$chr, levels = chrOrder)
gene_ordering_file <- gene_ordering_file[complete.cases(gene_ordering_file),]

# order file and remove duplicated values
gene_ordering_file <- gene_ordering_file[with(gene_ordering_file, order(chr, start)),]
gene_ordering_file <- gene_ordering_file[!duplicated(gene_ordering_file$hgnc_symbol),]

# merge both together
gene.exp.df <- merge(diff.exp.gene.df, gene_ordering_file, by = "hgnc_symbol", all.x = T)
gene.exp.df <- gene.exp.df[with(gene.exp.df, order(chr, start)),]

# count the number of genes on different chromosomes
chr.df <- prop.table(table(gene.exp.df$clone, gene.exp.df$chr), 2)
chr.df.melt <- melt(chr.df)
chr.df.melt$value <- round(chr.df.melt$value,4)

# calculate total number of samples
total.genes <- as.data.frame(colSums(table(gene.exp.df$clone, gene.exp.df$chr)))
total.genes$Var2 <- rownames(total.genes)
colnames(total.genes) <- c("total", "Var2")

# merge again
complete.melt <- merge(chr.df.melt, total.genes, by = "Var2")
colnames(complete.melt) <- c("chr", "Clones","proportion", "total")

# Small multiple
ggplot(complete.melt, aes(fill = Clones, y=proportion, x=chr)) + 
  geom_bar(position="stack", stat="identity",  color="black") +
  geom_text(aes(chr, 1.1, label = total, fill = NULL), data = complete.melt) +
  xlab("Chromosomes") + theme_classic() + ylab("Gene proportion") +
  scale_fill_manual(values = c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 14, face = "bold",),
        legend.text = element_text(colour="black", size=12, face="bold"), 
        legend.position = "top")
ggsave(paste0(sample.tmp, "/DEG/barplot_DEG_clone_chromosome.pdf"), width = 10, height = 5)


# count the number of genes on different chromosomes
chr.df <- prop.table(table(gene.exp.df$clone, gene.exp.df$chr), 1)
chr.df.melt <- melt(chr.df)
chr.df.melt$value <- round(chr.df.melt$value,4)

# calculate total number of samples
total.genes <- as.data.frame(rowSums(table(gene.exp.df$clone, gene.exp.df$chr)))
total.genes$Var2 <- rownames(total.genes)
colnames(total.genes) <- c("total", "Var1")

# merge again
complete.melt <- merge(chr.df.melt, total.genes, by = "Var1")
colnames(complete.melt) <- c("Clones", "chr", "proportion", "total")

nb.cols <- 22
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

# Small multiple
ggplot(complete.melt, aes(fill = chr, y=proportion, x=Clones)) + 
  geom_bar(position="stack", stat="identity",  color="black") +
  geom_text(aes(Clones, 1.1, label = total, fill = NULL), data = complete.melt) +
  xlab("Chromosomes") + theme_classic() + ylab("Gene proportion") +
  scale_fill_manual(values = mycolors) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 14, face = "bold",),
        legend.text = element_text(colour="black", size=12, face="bold"), 
        legend.position = "right")
ggsave(paste0(sample.tmp, "/DEG/barplot_DEG_chromosome_clone.pdf"), width = 10, height = 6)






