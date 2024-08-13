#############################################################################################################################
##                                                                                                                      
##  VISUALIZE CT SCORE IN NUCLEI AND PDX
##                                                                                                                      
##  Date: 05 FEBRUARY 2022                                                                                                                   
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

ct.score.data.nuclei <- read.csv("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/Nuclei_CT_scoring_bulk.csv")
ct.score.data.nuclei <- ct.score.data.nuclei[grep("Malignant", ct.score.data.nuclei$Cell_types),] 
ct.score.data.nuclei$sample <- factor(ct.score.data.nuclei$sample, levels = c("MB243", "STP", "ST1R"))
ct.score.data.nuclei <- ct.score.data.nuclei[,c("sample", "Upregulated_Chromothripsis_Score", "Downregulated_Chromothripsis_Score")]
ct.score.data.nuclei <- reshape2::melt(ct.score.data.nuclei)

my_comparisons <- list(c("MB243", "STP"), c("MB243", "ST1R"), c("STP", "ST1R"))

cols <- c(brewer.pal(n=6, "Set1"))
show_col(cols)
cols <- cols[c(2,3, 5)]

## COPY NUMBER
pdf("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/CT_boxplot_upregulated_samples.pdf", width = 5, height = 5)
ggplot(ct.score.data.nuclei, aes(sample, value, color = sample)) +
  # geom_violin(outlier.colour="black", outlier.shape= NA,
  #              outlier.size=2, notch=FALSE) +
  geom_boxplot(width = 0.4,outlier.colour="black", outlier.shape= NA,
               outlier.size=2, notch=FALSE) +
  geom_jitter(size = 0.5, width = 0.15, alpha = 0.05) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme_classic() + 
  scale_color_manual(values = cols) +
  facet_wrap(~ variable, scales='free_x', nrow = 2) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="Gene Expression Score") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.75),
        strip.text = element_text(face="bold", size=10, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1),
        # axis.ticks.x=element_blank(),
        axis.text.x = element_text(color = "black", size = 8, angle = 0, hjust = .5, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 8, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 10, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 8, hjust = .5, face = "bold"), 
        legend.position = "none")
dev.off()

## COPY NUMBER
pdf("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/CT_boxplot_downregulated_samples.pdf", width = 5, height = 2)
ggplot(ct.score.data.nuclei, aes(sample, Downregulated_Chromothripsis_Score, color = sample)) +
  # geom_violin(outlier.colour="black", outlier.shape= NA,
  #              outlier.size=2, notch=FALSE) +
  geom_boxplot(width = 0.4,outlier.colour="black", outlier.shape= NA,
               outlier.size=2, notch=FALSE) +
  geom_jitter(size = 0.5, width = 0.15, alpha = 0.1) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme_classic() + 
  # scale_color_lancet(palette = "lanonc") +
  scale_color_manual(values = cols) +
  # facet_wrap(~ sample, scales='free_x', nrow = 1) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="CT Score - Upregulated") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.75),
        strip.text = element_text(face="bold", size=6, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1),
        # axis.ticks.x=element_blank(),
        axis.text.x = element_text(color = "black", size = 8, angle = 0, hjust = .5, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 8, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 10, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 8, hjust = .5, face = "bold"), 
        legend.position = "none")
dev.off()

############################################################################
##                          SPLIT BY CELL TYPE
############################################################################

# define the comparisons
my_comparisons <- list(c("Malignant SHH II", "Malignant SHH I"), c("Malignant SHH II", "Malignant Neuronal Development II"), c("Malignant SHH II", "Malignant Cycling"), c("Malignant SHH II", "Malignant Neuronal Development I"),
                       c("Malignant SHH I", "Malignant Neuronal Development II"), c("Malignant SHH I", "Malignant Cycling"), c("Malignant SHH I", "Malignant Neuronal Development I"),
                       c("Malignant Neuronal Development II", "Malignant Cycling"), c("Malignant Neuronal Development II", "Malignant Neuronal Development I"),
                       c("Malignant Cycling", "Malignant Neuronal Development I"))

## COPY NUMBER
pdf("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/CT_boxplot_upregulated_celltypes.pdf", width = 9, height = 3)
ggplot(ct.score.data.nuclei, aes(Cell_types, Upregulated_Chromothripsis_Score, color = Cell_types)) +
  # geom_violin(outlier.colour="black", outlier.shape= NA,
  #              outlier.size=2, notch=FALSE) +
  geom_boxplot(outlier.colour="black", outlier.shape= NA,
               outlier.size=2, notch=FALSE) +
  geom_jitter(size = 0.5, width = 0.25) + 
  theme_classic() + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  scale_color_lancet(palette = "lanonc") +
  # facet_wrap(~ Cell_types, scales='free_x', nrow = 1) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="CT Score - Upregulated") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.75),
        strip.text = element_text(face="bold", size=6, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1),
        # axis.ticks.x=element_blank(),
        axis.text.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 6, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 8, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 8, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"), 
        legend.position = "none")
dev.off()


# define the comparisons
my_comparisons <- list(c("Malignant SHH II", "Malignant SHH I"), c("Malignant SHH II", "Malignant Neuronal Development II"), c("Malignant SHH II", "Malignant Cycling"), c("Malignant SHH II", "Malignant Neuronal Development I"),
                       c("Malignant SHH I", "Malignant Neuronal Development II"), c("Malignant SHH I", "Malignant Cycling"), c("Malignant SHH I", "Malignant Neuronal Development I"),
                       c("Malignant Neuronal Development II", "Malignant Cycling"), c("Malignant Neuronal Development II", "Malignant Neuronal Development I"),
                       c("Malignant Cycling", "Malignant Neuronal Development I"))

## COPY NUMBER
pdf("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/CT_boxplot_downregulated_celltypes.pdf", width = 9, height = 3)
ggplot(ct.score.data.nuclei, aes(Cell_types, Downregulated_Chromothripsis_Score, color = Cell_types)) +
  # geom_violin(outlier.colour="black", outlier.shape= NA,
  #              outlier.size=2, notch=FALSE) +
  geom_boxplot(outlier.colour="black", outlier.shape= NA,
               outlier.size=2, notch=FALSE) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  geom_jitter(size = 0.5, width = 0.25) + 
  theme_classic() + 
  scale_color_lancet(palette = "lanonc") +
  facet_wrap(~ Cell_types, scales='free_x', nrow = 1) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(x="",
       y="CT Score - Downregulated") + 
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.75),
        strip.text = element_text(face="bold", size=6, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1),
        # axis.ticks.x=element_blank(),
        axis.text.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 6, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 8, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 8, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"), 
        legend.position = "none")
dev.off()

