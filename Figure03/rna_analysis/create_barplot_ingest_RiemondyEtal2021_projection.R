#############################################################################################################################
##                                                                                                                      
##  VISUALIZE THE PROPORTION OF CELLS ASSIGNED TO EACH SUBGROUP PER SAMPLE
##                                                                                                                      
##  Date: 18 OCTOBER 2021                                                                                                                   
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
                      "Matrix", "devtools", "Matrix.utils", "matrixStats", "readr", "magrittr", "BiocManager", 
                      "biomaRt", "httr", "ComplexHeatmap", "fields", "data.table", "dbscan", "fpc", "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

############################################################################
##                          READ IN THE DATA
############################################################################

# read the csv file for the ingest analysis
ingest.projection.data <- read.csv("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/Nuclei_projection_Riemondy_metadata.csv")

# remove normal cell types
ingest.projection.data <- ingest.projection.data[grep("Malignant", ingest.projection.data$Cell_types),]

# get the proportions of each cell type in the respective condition
sample.projection <-  prop.table(table(ingest.projection.data$subgroup, ingest.projection.data$sample),2)
sample.projection.melt <- reshape2::melt(sample.projection)
sample.projection.melt$value <- round(sample.projection.melt$value,4)

# Small multiple
ggplot(sample.projection.melt, aes(fill = Var1, y=Var2, x=value)) + 
  geom_bar(position="stack", stat="identity",  color="black") +
  xlab("Cell proportion") + theme_classic() + ylab("") +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 14, face = "bold",),
        legend.text = element_text(colour="black", size=12, face="bold"), 
        legend.position = "top")
ggsave("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/ingest_barplot_Riemondy.pdf", width = 4, height = 3)


############################################################################
##                    VISUALIZE CELL TYPE MAPPING
############################################################################

ingest.projection.data <- read.csv("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/Nuclei_projection_Riemondy_metadata_Celltypes.csv")

celltype.table <- table(ingest.projection.data$Cell_types, ingest.projection.data$coarse_cell_type)
celltype.table <- round(prop.table(as.matrix(celltype.table[-1,]), 1), 4)
celltype.freq <- as.matrix.data.frame(celltype.table)
rownames(celltype.freq) <- c("Astro", "Endo", "Micro", "Cycling", "Neu Dev I", "Neu Dev II", "SHH I", "SHH II", "Men", "Purk")
colnames(celltype.freq) <- c("Lymphocytes", "Macrophages", "Malignant", "Oligodendrocytes")

### ADD IN BREAKS AND COLORS (FOR PVALUE)
colors = c(0,0.5,1)
my_palette <- c("white", "#005b96","#011f4b")

# generate the first heatmap based on the score matrix
ht1 <- Heatmap(as.matrix(celltype.freq), 
               name = "ContingencyTable", 
               col = circlize::colorRamp2(colors, my_palette), 
               cluster_rows = T, 
               cluster_columns = FALSE,
               show_row_names = TRUE, 
               row_names_side = "left",
               show_row_dend = FALSE,
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 12, fontface = "bold"), 
               column_names_gp = gpar(fontsize = 12, fontface = "bold", angle = 45), 
               rect_gp = gpar(col = "black", lwd = 2), 
               column_title_gp = gpar(fontsize = 15, fontface = "bold"),
               heatmap_legend_param = list(color_bar = "continous",
                                           at = c(0, 0.25, 0.5, 0.75, 1),
                                           title = ""),border = T)

pdf(paste0("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/Ingest_Riemondy_Nuclei_ContingencyTable.pdf"), width=5,height = 4,pointsize=0.1)
print(ht1)
dev.off()



