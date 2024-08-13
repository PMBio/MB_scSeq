#############################################################################################################################
##                                                                                                                      
##  GENERATE SUMMARY PLOTS FOR THE WHOLE GUT ATLAS DATASET
##                                                                                                                      
##  Date: 26 OCTOBER 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
############################################################################################################################
# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("BiocManager","tidyverse", "ggplot2", "ComplexHeatmap", "reshape2", 
                      "RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

############################################################################
##                            FIGURE 5E - NUCLEI
############################################################################
# read in the metadata of interest
metadata <- read.csv("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/Nuclei_metadata_final.csv")

# get the proportions of each cell type in the respective condition
celltype.df <- prop.table(table(metadata$Cell_types, metadata$sample),1)
celltype.df.melt <- reshape2::melt(celltype.df)
celltype.df.melt$value <- round(celltype.df.melt$value,4)
celltype.df.melt$Var1 <- factor(celltype.df.melt$Var1, levels = c("Astrocytes", "Endothelial Cells", "Macrophages", 
                                                                  "Malignant SHH I", "Malignant SHH II",
                                                                  'Malignant Cycling',
                                                                  'Malignant Neuronal Development I', 'Malignant Neuronal Development II',
                                                                  "Meninge Cells", "Purkinje Cells"))

celltype.df.melt$Var2 <- factor(celltype.df.melt$Var2, levels = c("LFSMBP", "LFSMB1R", "MB243"))

# Small multiple
ggplot(celltype.df.melt, aes(fill = Var2, y=value, x=Var1)) + 
  geom_bar(position="stack", stat="identity",  color="black") +
  xlab("") + theme_classic() + ylab("Cell proportion") +
  scale_fill_brewer(palette = "Blues") +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 45, hjust = 1),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 14, face = "bold",),
        legend.text = element_text(colour="black", size=12, face="bold"), 
        legend.position = "bottom")
ggsave("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/ParraEtal_Figure3_Nuclei_Barplot.pdf", width = 5, height = 6)

############################################################################
##                           FIGURE 5E - PDX
############################################################################

# read in the metadata of interest
metadata <- read.csv("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/PDX_metadata_final.csv")

# get the proportions of each cell type in the respective condition
celltype.df <- prop.table(table(metadata$Cell_types, metadata$sample),1)
celltype.df.melt <- melt(celltype.df)
celltype.df.melt$value <- round(celltype.df.melt$value,4)
celltype.df.melt$Var1 <- factor(celltype.df.melt$Var1, levels = c("Astrocytes", "Malignant Basal State", 'Malignant Cycling',
                                                                  'Malignant Granule-like Progenitor','Malignant Neuronal Development I',
                                                                  'Malignant Neuronal Development II', 'Malignant Neuronal Development III',
                                                                  "Granule Cells", "Microglia"))
celltype.df.melt$Var2 <- factor(celltype.df.melt$Var2, levels = c("LFSMBP", "LFSMB1R", "RCMB18", "BT084"))

# Small multiple
ggplot(celltype.df.melt, aes(fill = Var2, y=value, x=Var1)) + 
  geom_bar(position="stack", stat="identity",  color="black") +
  xlab("") + theme_classic() + ylab("Cell proportion") +
  scale_fill_brewer(palette = "Greens") +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 45, hjust = 1),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 14, face = "bold",),
        legend.text = element_text(colour="black", size=12, face="bold"), 
        legend.position = "bottom")
ggsave("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/ParraEtal_Figure3_PDX_Barplot.pdf", width = 5, height = 6)
