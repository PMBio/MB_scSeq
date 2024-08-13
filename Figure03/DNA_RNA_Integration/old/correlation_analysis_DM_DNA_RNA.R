#############################################################################################################################
##                                                                                                                      
##  MAKE VIOLINPLOTS FOR GENES IN DOUBLE MINUTES TARGETS
##                                                                                                                      
##  Date: 08 APRIL 2021                                                                                                                   
##  
##  Author: Moritz Przybilla
##
##                                                                                                                      
############################################################################################################################
# clear workspace
rm(list=ls())
set.seed(16011985) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "optparse", "BSgenome", "RColorBrewer", "ggplot2", "scales", "DescTools", "dendextend", "tidyverse", 
                      "Matrix", "devtools", "Matrix.utils", "matrixStats", "readr", "magrittr", "Signac", "BiocManager", 
                      "biomaRt", "httr", "ComplexHeatmap", "fields", "data.table", "Seurat", "ggrepel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

############################################################################
##                          READ IN THE DATA
############################################################################

# negate %in%
"%ni%" <- Negate("%in%")

# specify the output directory
out.dir <- "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA"

# read in the druggable targets
dm.genes <- read.table("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA/DM_cnv_values_per_gene.txt", header = T, stringsAsFactors = F)
samples <- unique(dm.genes$Sample)

# read in the expression files from scanpy
mtx.files <- list.files("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/scRNA_analysis/scanpy", pattern = "raw_matrix", full.names = T)
mtx.files <- mtx.files[grep("MB243|LFSMBP-Nuclei|RCMB18",mtx.files)]

# list the metadata files
metadata.files <- list.files('/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA', pattern = "_metadata.csv$", full.names = T)
metadata.files <- metadata.files[-grep("chromothripsisScore", metadata.files)]

# list to store
dm.gene.sample.list <- list()

i <- 1
for (i in 1:length(mtx.files)){
  
  # read in the matrix
  sample.tmp <- str_split_fixed(basename(mtx.files[i]), "_", 2)[,1]
  sample.mtx <- read.csv(mtx.files[i], header = T)
  rownames(sample.mtx) <- sample.mtx$X
  
  # read in the metadata
  sample.metadata <- fread(grep(sample.tmp, metadata.files, value = T), header = T)
  sample.metadata <- sample.metadata[,c("V1", "clone_id", "new_clusters")]
  colnames(sample.metadata) <- c("Cell_barcodes", "clone_id", "celltype")
  
  # get the dm genes of interest for this sample
  sample.dm.genes <- dm.genes[dm.genes$Sample == sample.tmp,]
  
  # # read in the raw 10x matrix
  # tenx.mtx <- Read10X(raw.mtx.files[grep(sample.tmp, raw.mtx.files)])
  # tenx.mtx <- tenx.mtx[,rownames(sample.mtx)]
  # dim(tenx.mtx)
  # 
  # # transpose
  # tenx.mtx <- t(tenx.mtx)
  # 
  # # # subset the matrix
  # dm.genes.mtx <- as.data.frame(as.matrix(tenx.mtx[,colnames(tenx.mtx) %in% sample.dm.genes$Gene]))
  # dim(dm.genes.mtx)
  # dm.genes.mtx$Cell_barcodes <- rownames(dm.genes.mtx)
  
  # subset the matrix
  dm.genes.mtx <- sample.mtx[,colnames(sample.mtx) %in% sample.dm.genes$Gene]
  dm.genes.mtx$Cell_barcodes <- rownames(dm.genes.mtx)
  
  # melt it 
  dm.genes.mtx.melt <- melt(dm.genes.mtx)
  
  # merge with clone id
  dm.genes.df <- merge(dm.genes.mtx.melt, sample.metadata, by = "Cell_barcodes")
  dm.genes.df$sample <- sample.tmp
  
  # only keep clones with more than 10 cells
  valid.clones <- names(table(sample.metadata$clone_id) > 10)[table(sample.metadata$clone_id) > 10]
  
  # exclude 0s
  dm.genes.df$clone_id <- factor(dm.genes.df$clone_id, levels = c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5", "Clone6"))
  dm.genes.df <- dm.genes.df[dm.genes.df$clone_id %in% valid.clones,]
  dm.genes.df <- dm.genes.df[dm.genes.df$value > 0,]
  dm.genes.df <- dm.genes.df[!is.na(dm.genes.df$clone_id),]
  
  # aggregate the expression by clone
  mean.dm.genes.df <- aggregate(dm.genes.df[,"value"], list(dm.genes.df[,"variable"], dm.genes.df[,"clone_id"]), FUN = mean)
  mean.dm.genes.df$sample <- sample.tmp
  mean.dm.genes.df <- mean.dm.genes.df[,c("sample", "Group.2", "Group.1", "x")]
  colnames(mean.dm.genes.df) <- c("Sample", "clone_id", "Gene", "Expression")
  
  # process the sample dm values to combine them
  sample.dm.genes <- separate(sample.dm.genes, "Clones_Mean", c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5", "Clone6"), sep = "-")
  sample.dm.genes <- sample.dm.genes[,c("Sample", "Gene", "Clone1", "Clone2", "Clone3", "Clone4", "Clone5", "Clone6")]
  sample.dm.melt <- melt(sample.dm.genes, id.vars = c("Sample", "Gene"))
  sample.dm.melt <- sample.dm.melt[sample.dm.melt$Gene %in% dm.genes.df$variable, ]
  sample.dm.melt <- sample.dm.melt[,c("Sample", "variable", "Gene", "value")]
  colnames(sample.dm.melt) <- c("Sample", "clone_id", "Gene", "cn_value")
  
  # merge both datasets
  complete.dm.genes.df <- merge(mean.dm.genes.df, sample.dm.melt, by = c("Sample", "clone_id", "Gene"))
  complete.dm.genes.df$cn_value <- as.numeric(complete.dm.genes.df$cn_value)
  
  # add to list
  dm.gene.sample.list[[i]] <- complete.dm.genes.df
}

complete.dm.gene.df <- bind_rows(dm.gene.sample.list)

complete.dm.gene.df$Sample <- factor(complete.dm.gene.df$Sample, levels = c("MB243-Nuclei", "LFSMBP-Nuclei", "RCMB18-PDX"))
complete.dm.gene.df <- complete.dm.gene.df[complete.dm.gene.df$Sample == "MB243-Nuclei",]

cols <- c(brewer.pal(n=6, "Set1"))
cols <- cols[-c(1,6)]

# plot the data
pdf(paste0(out.dir, "/MB243_scDNA_scRNA_DM_gene_cnv_correlation_wo_4M67.pdf"), width = 5, height = 4)
p <- ggplot(complete.dm.gene.df, aes(x=Expression, y = cn_value, label = Gene)) +
  geom_point(size = 0.75, aes(color = clone_id)) +
  theme_classic() +
  scale_color_manual(values = cols) +
  labs(x="Average normalized gene expression per integrated clone [scRNA]",
       y="Average copy number per clone [scDNA]") + 
  theme(legend.position="right") +
  
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.75),
        axis.text = element_text(colour = "black", size = 8, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 8, face = "bold" ),
        axis.title = element_text(colour = "black", size = 10, face = "bold" ),
        plot.title = element_text(colour = "black", size = 10, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"),
        legend.position = "none") +
  geom_text_repel(data          = complete.dm.gene.df %>% group_by(clone_id) %>% top_n(3, Expression),
                  size          = 3, fontface = "bold",
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x")
ggExtra::ggMarginal(p, type = "histogram", groupColour = TRUE, groupFill = TRUE, xparams = list(binwidth = 0.05), yparams = list(binwidth = 0.5))
dev.off()
ggsave(paste0(out.dir, "/MB243_scDNA_scRNA_DM_gene_cnv_correlation_wo_4M67.pdf"), width = 6, height = 5, dpi = 600)

