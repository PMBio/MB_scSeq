#############################################################################################################################
##                                                                                                                      
##  PLOT CONTINGENCY TABLES ACROSS SAMPLES
##                                                                                                                      
##  Date: 25 MARCH 2021                                                                                                                   
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
                      "biomaRt", "httr", "ComplexHeatmap", "fields", "data.table", "dbscan", "fpc", "Seurat", "mclust")
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
out.dir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA"

# list files from the scRNA scDNA integration
scRNA.contingency.files <- list.files(out.dir, pattern = "contingency", recursive = T, full.names = T, all.files = T)
scRNA.contingency.files <- scRNA.contingency.files[-c(1, 4)]

# read in the files
scRNA.contingency.data <- lapply(scRNA.contingency.files, read.table, sep = ",", header = T)

i <- 1
for (i in 1:length(scRNA.contingency.files)){
  
  # get the data
  sample.tmp <- str_split_fixed(basename(scRNA.contingency.files[i]), "_", 2)[,1]
  scRNA.data <- scRNA.contingency.data[[i]]
  colnames(scRNA.data)[ncol(scRNA.data)] <- "Unassigned_Normal"
  rownames(scRNA.data) <- scRNA.data$new_clusters
  scRNA.data[,c("new_clusters", "Unassigned_Normal")] <- NULL
  print(dim(scRNA.data))
  
  # make it a frequency table
  freq.scRNA.data <- round(prop.table(as.matrix(scRNA.data), 1), 4)
  freq.scRNA.data <- freq.scRNA.data[complete.cases(freq.scRNA.data),]
  colnames(freq.scRNA.data) <- factor(colnames(freq.scRNA.data), levels = c(paste0("Clone", c(1:(ncol(scRNA.data))))))
  
  ### ADD IN BREAKS AND COLORS (FOR PVALUE)
  colors = c(0,0.5,1)
  my_palette <- c("white", "#005b96","#011f4b")
  
  # generate the first heatmap based on the score matrix
  ht1 <- Heatmap(freq.scRNA.data, 
                 name = "ContingencyTable", 
                 col = circlize::colorRamp2(colors, my_palette), 
                 cluster_rows = F, 
                 cluster_columns = FALSE,
                 show_row_names = TRUE, 
                 row_names_side = "left",
                 show_row_dend = FALSE,
                 show_column_names = TRUE,
                 row_names_gp = gpar(fontsize = 9), 
                 rect_gp = gpar(col = "black", lwd = 1), 
                 column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                 heatmap_legend_param = list(color_bar = "continous",
                                             at = c(0, 0.25, 0.5, 0.75, 1),
                                             title = ""),border = T)
  
  pdf(paste0(out.dir, "/", sample.tmp, "/Celltype_Clones_ContingencyTable.pdf"),width=10,height = 7,pointsize=0.1)
  print(ht1)
  dev.off()
  
  
}
