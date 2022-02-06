#############################################################################################################################
##                                                                                                                      
##  CREATE DATA FOR DRUGGABLE TARGET BOXPLOTS IN FIGURE 1
##                                                                                                                      
##  Date: 26 MARCH 2021                                                                                                                   
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

# read in the expression files from scanpy
mtx.files <- list.files("/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/scRNA_analysis/scanpy", pattern = "raw_matrix", full.names = T)

# list the metadata files
metadata.files <- list.files('/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA', pattern = "_metadata.csv$", full.names = T)
metadata.files <- metadata.files[-grep("chromothripsisScore", metadata.files)]

# list the raw expression matrix from 10x
raw.mtx.files <- list.files("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10XRNA5P/", pattern = "filtered_feature_bc_matrix", recursive = T, full.names = T, include.dirs = T)
raw.mtx.files <- raw.mtx.files[-grep("\\.h5|mm10|old", raw.mtx.files)]

# read in the druggable targets
druggable.targets <- read.table("/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/DruggableTargets/DruggableTargetsAurelie.txt", header = F, stringsAsFactors = F)

############################################################################
##                        PROCESS THE DATA
############################################################################

# list to store
druggable.sample.list <- list()

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
  
  # subset the matrix
  druggable.mtx <- sample.mtx[,colnames(sample.mtx) %in% druggable.targets$V1]
  druggable.mtx$Cell_barcodes <- rownames(druggable.mtx)

  # melt it 
  druggable.mtx.melt <- reshape2::melt(druggable.mtx)
  
  # merge with clone id
  druggable.df <- merge(druggable.mtx.melt, sample.metadata, by = "Cell_barcodes")
  druggable.df$sample <- sample.tmp
  
  # exclude 0s
  unassigned.tumor.cells.df <- druggable.df[druggable.df$clone_id == "nan",]
  unassigned.tumor.cells.df <- unassigned.tumor.cells.df[grep("malignant|neuronal development I|neuronal development II|neuronal development", unassigned.tumor.cells.df$celltype),]
  
  # exclude tumor cell barcodes which have not been assigned
  druggable.df <- druggable.df[druggable.df$Cell_barcodes %ni% unassigned.tumor.cells.df$Cell_barcodes, ]
  normal.celltypes <- unique(druggable.df[druggable.df$clone_id == "nan", "celltype"])
  druggable.df[druggable.df$clone_id == "nan", "clone_id"] <- druggable.df[druggable.df$clone_id == "nan", "celltype"]
  druggable.df$clone_id <- factor(druggable.df$clone_id, levels = c(normal.celltypes, "Clone1", "Clone2", "Clone3", "Clone4", "Clone5", "Clone6"))
  
  # exclude 0s
  # druggable.df <- druggable.df[druggable.df$value > 0,]
  druggable.df <- druggable.df[!is.na(druggable.df$clone_id),]
  
  # add to list
  druggable.sample.list[[i]] <- druggable.df
}

complete.druggable.targets.df <- bind_rows(druggable.sample.list)

# exclude non- and lowly expressed genes and write to file
druggable.non.zero.df <- complete.druggable.targets.df[complete.druggable.targets.df$value > 0,]
genes <- names(table(druggable.non.zero.df$variable))[table(druggable.non.zero.df$variable) > 50]
complete.druggable.targets.df <- complete.druggable.targets.df[complete.druggable.targets.df$variable %in% genes,]
colnames(complete.druggable.targets.df) <- c("Cell_barcodes", "Genes", "Expression", "Clone", "Celltype", "Sample")
write.table(complete.druggable.targets.df, paste0(out.dir, "/druggableGenes_allSamples_expression_w_zeros.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
