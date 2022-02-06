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
list.of.packages <- c("reshape2", "optparse", "BSgenome", "RColorBrewer", "ggplot2", "scales", "DescTools", 
                      "Matrix", "devtools", "Matrix.utils", "matrixStats", "readr", "magrittr", "BiocManager", 
                      "biomaRt", "httr", "data.table", "Seurat", "SeuratDisk", "fields", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))
httr::set_config(httr::config(ssl_verifypeer=0L))

############################################################################
##                          READ IN THE DATA
############################################################################

# read in the cell type information for AldingerEtal2020
aldinger.data <- read.csv('/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/AldingerEtal2020_Celltypes_pseudobulk.csv')
aldinger.data <- as.data.frame(t(aldinger.data))
colnames(aldinger.data) <- aldinger.data[1,]
aldinger.data <- aldinger.data[-1,]

# read in the DEGs per cell type
degs.aldinger <- read.csv("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/AldingerEtal2020_DEG_wilcoxon.csv")
colnames(degs.aldinger) <- gsub("X", "C", colnames(degs.aldinger)) 

# read in the scanpy anndata object
Convert("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/adata_nuclei.h5ad", dest = "h5seurat", overwrite = TRUE)
nuclei.data <- LoadH5Seurat("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/adata_nuclei.h5seurat")

# get the matrix and metadata
metadata <- read.csv('/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/Nuclei_metadata_final.csv')
cell.type <- metadata[, c("X", "Cell_types")]
rownames(cell.type) <- cell.type$X
cell.type$X <- NULL

# add metadata
nuclei.data <- AddMetaData(nuclei.data, cell.type, "Cell_types")

# calculate average expression per cell type
Idents(nuclei.data) <- "Cell_types"
average.nuclei.exp <- AverageExpression(nuclei.data)

# subset to malignant cell types
malignant.avg.exp <- as.data.frame(average.nuclei.exp$RNA)
malignant.avg.exp <- malignant.avg.exp[,grep("Malignant", colnames(malignant.avg.exp))]

# reduce to genes which are sufficiently expressed
gene.mean.exp <- rowMeans(malignant.avg.exp)
mean.gene.exp <- rowMeans(nuclei.data@assays$RNA)
q50.gene.exp <- quantile(mean.gene.exp, 0.50)

# remove genes below this threshold
final.genes <- names(gene.mean.exp)[gene.mean.exp >= q50.gene.exp] 
malignant.avg.exp <- malignant.avg.exp[rownames(malignant.avg.exp) %in% final.genes,]

# get cell types
malignant.celltypes <- colnames(malignant.avg.exp)

# initialize dataframe in which all information about the
# alignment will be stored
alignment.data <- data.frame("Celltype" = "Malignant", max_r2 = "pearson.correlation", "max_celltype" = NA, stringsAsFactors = F)
aldinger.celltypes <- colnames(aldinger.data) 
aldinger.r2.data <- as.data.frame(matrix(nrow = 1, ncol = length(aldinger.celltypes)))
colnames(aldinger.r2.data) <- aldinger.celltypes
alignment.r2.data <- cbind(alignment.data, aldinger.r2.data)

############################################################################
##            CORRELATE THE MALIGNANT CELLS TO THE CELLS TYPES
############################################################################

i <- 1
# iterate over each cell and compare its profile to the scDNA clones
for (i in 1:ncol(malignant.avg.exp)){
  
  # pick one cell at a time
  malignant.tmp <- colnames(malignant.avg.exp)[i]
  print(paste0("Cell type: ", malignant.tmp))
  
  malignant.data.tmp <- malignant.avg.exp
  malignant.data.tmp$hgnc_symbol <- rownames(malignant.data.tmp)
  malignant.data.tmp <- malignant.data.tmp[, c("hgnc_symbol", malignant.tmp)]
  alignment.r2.data[i, "Celltype"] <- malignant.tmp
  
  j <- 1
  for (j in 1:length(aldinger.celltypes)){
    
    aldinger.tmp <- aldinger.celltypes[j]
    print(aldinger.tmp)
    
    # get the aldinger data
    aldinger.data.tmp <- aldinger.data
    aldinger.data.tmp$hgnc_symbol <- rownames(aldinger.data.tmp)
    
    # 
    aldinger.data.tmp <- aldinger.data.tmp[, c("hgnc_symbol", aldinger.tmp)]
    
    # get the differentially expressed genes for this cell type
    degs.tmp <- degs.aldinger[, grep(str_split_fixed(aldinger.tmp, "-", 2)[,1], colnames(degs.aldinger))]
    
    # subset to the significantly expressed genes and positive log2FC
    degs.tmp <- degs.tmp[degs.tmp[,4] < 0.00001, ]
    common.genes <- degs.tmp[degs.tmp[,1] %in% malignant.data.tmp$hgnc_symbol, 1]
    
    # intersect genes
    intersected.genes <- aldinger.data.tmp[aldinger.data.tmp$hgnc_symbol %in% common.genes, "hgnc_symbol"]
    print(paste0("Number of Genes: ", length(intersected.genes)))
    # common.genes <- aldinger.data.tmp[aldinger.data.tmp$hgnc_symbol %in% malignant.data.tmp$hgnc_symbol, "hgnc_symbol"]
    
    # subset to common genes
    aldinger.vector <- aldinger.data.tmp[intersected.genes, ]
    malignant.vector <- malignant.data.tmp[intersected.genes,]
    aldinger.vector <- as.numeric(as.character(aldinger.vector[, aldinger.tmp]))
    malignant.vector <- as.numeric(as.character(malignant.vector[, malignant.tmp]))
    
    # calculate the correlation to each possible scDNA clone
    r2_correlation <- cor(malignant.vector, aldinger.vector, method = "pearson")
    plot(malignant.vector, aldinger.vector)
    
    # store correlation value
    alignment.r2.data[i,j+3] <- r2_correlation
    
  }
  
  # determine the maximum
  best.celltype.match <- colnames(alignment.r2.data)[grep(max(alignment.r2.data[i,4:ncol(alignment.r2.data)]), alignment.r2.data)]
  alignment.r2.data[i,"max_r2"] <- max(alignment.r2.data[i,4:ncol(alignment.r2.data)])
  alignment.r2.data[i,"max_celltype"] <- best.celltype.match
  
}
  
# wrangel into shape
alignment.mtx <- alignment.r2.data
alignment.mtx[,c("max_r2", "max_celltype")] <- NULL
rownames(alignment.mtx) <- alignment.mtx$Celltype
alignment.mtx$Celltype <- NULL
alignment.mtx <- as.matrix(alignment.mtx)

# check the minimum and maximum
max <- max(alignment.mtx)
min <- min(alignment.mtx)

### ADD IN BREAKS AND COLORS (FOR PVALUE)
colors = c(-0.3, -0.15, 0, 0.15, 0.3)
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(5)

# generate the first heatmap based on the score matrix
ht1 <- Heatmap(t(alignment.mtx), 
               name = "Celltype_Relation", 
               col = circlize::colorRamp2(colors, my_palette), 
               cluster_rows = T, 
               cluster_columns = T,
               show_row_names = TRUE, 
               row_names_side = "left",
               show_row_dend = FALSE,
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 12, fontface = "bold"), 
               column_names_gp = gpar(fontsize = 12, fontface = "bold"), 
               rect_gp = gpar(col = "black", lwd = 1), 
               column_title_gp = gpar(fontsize = 15, fontface = "bold"),
               heatmap_legend_param = list(color_bar = "continous",
                                           at = c(-0.3, -0.15, 0, 0.15, 0.3),
                                           title = "Correlation"),border = T)

pdf(paste0("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/Celltype_Clones_ContingencyTable.pdf"),width=5,height = 10,pointsize=0.1)
print(ht1)
dev.off()







