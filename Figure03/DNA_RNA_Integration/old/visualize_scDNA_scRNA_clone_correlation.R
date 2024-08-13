#############################################################################################################################
##                                                                                                                      
##  PLOT CONCORDANCE BETWEEN SCRNA AND SCDNA INTEGRATED CLONES
##                                                                                                                      
##  Date: 09 MARCH 2021                                                                                                                   
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
                      "biomaRt", "httr", "ComplexHeatmap", "fields", "data.table", "dbscan", "fpc", "Seurat", "ggpubr", "ggsci")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

############################################################################
##                          READ IN THE DATA
############################################################################

# negate %in%
"%ni%" <- Negate("%in%")

# specify the output directory
out.dir <- "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA/"

# list the file with clusters
scDNA.cluster.files <- list.files("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/DNA_processed_data_cellranger-dna_v.1.1.0/", pattern = "cell_ids_clusters_normalized", recursive = T, full.names =  T)
scDNA.cluster.files <- scDNA.cluster.files[grep("Final_Clustering_Trees", scDNA.cluster.files)]
scDNA.cluster.files <- scDNA.cluster.files[-c(1,7)]

# read in the files
scDNA.cluster.data <- lapply(scDNA.cluster.files, read.table, sep = "\t", header = T)

# list files from the scRNA scDNA integration
scRNA.cluster.files <- list.files(out.dir, pattern = "_scDNA_clones_filtered_cells.txt", recursive = T, full.names = T, all.files = T)
scRNA.cluster.files <- scRNA.cluster.files[-1]

# read in the files
scRNA.cluster.data <- lapply(scRNA.cluster.files, read.table, sep = "\t", header = T)

final.dataframe <- list()
i <- 1
for (i in 1:length(scDNA.cluster.files)){
  
  # get the data
  sample.tmp <- str_split_fixed(scDNA.cluster.files[i], "/", 12)[,10]
  scDNA.data <- scDNA.cluster.data[[i]]
  print(dim(scDNA.data))
  scRNA.data <- scRNA.cluster.data[[i]]
  print(dim(scRNA.data))
  
  # get the proportions of each cell type in the respective condition
  scDNA.prop <- prop.table(table(scDNA.data$cluster))
  scDNA.prop.melt <- melt(scDNA.prop)
  colnames(scDNA.prop.melt) <- c("clone_id", "scDNA_freq")
  scDNA.prop.melt$clone_id <- paste0("Clone", scDNA.prop.melt$clone_id)
  scDNA.prop.melt$sample <- sample.tmp
  
  # get the proportions of each cell type in the respective condition
  scRNA.data <- scRNA.data[scRNA.data$padj <= 0.05, ]
  scRNA.prop <- prop.table(table(scRNA.data$clone_id))
  scRNA.prop.melt <- melt(scRNA.prop)
  colnames(scRNA.prop.melt) <- c("clone_id", "scRNA_freq")
  scRNA.prop.melt$sample <- sample.tmp
  
  # combine both and store them
  complete.data <- merge(scDNA.prop.melt, scRNA.prop.melt, by = c("clone_id", "sample"))
  final.dataframe[[i]] <- complete.data
  
}

final.dataframe <- bind_rows(final.dataframe)
final.dataframe$scDNA_freq <- final.dataframe$scDNA_freq*100
final.dataframe$scRNA_freq <- final.dataframe$scRNA_freq*100

ggplot(final.dataframe, aes(x=scDNA_freq, y = scRNA_freq)) +
  geom_point(size = 2, aes( color = sample)) +
  theme_classic() +
  geom_smooth(method=lm, color="black", size = 1.5) +
  stat_cor(method = "pearson") +
  labs(x="Clone frequency in scDNA [%]",
       y="Clone frequency in scRNA [%]") + 
  theme(legend.position="right") +
  scale_color_lancet(palette = "lanonc") +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"))
ggsave(paste0(out.dir, "scDNA_scRNA_clone_correlation.pdf"), width = 10, height = 7, dpi = 600)
