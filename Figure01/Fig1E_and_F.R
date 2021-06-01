# ---------------------------------------------------------
# This scripts is used to generate subpanels in Fig1A
# Autor: R. Gonzalo Parra
# May 16, 2021
#
# It reads previously calculated CT score matrix for all samples
# ---------------------------------------------------------


plot_sample_CT<-function(sample, min_consecutive)
{
  
  clusters <-read.table(paste0("/home/r511a/projects/Gonzalo/Meduloblastoma/ScriptsPaper/ClusterTables/",sample,".txt"), sep = "\t", stringsAsFactors = F, header = T)
  clusters[,"cell_id"]<-paste0("X",clusters[,"cell_id"])
  cluster_counts<-table(clusters$cluster)
  
  MetaCell_Matrix_chromatripsis_cells<-readRDS(paste0("/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/CT_Data/MetaCellsCT_",sample,"_",1,".RDS"))
  dim(MetaCell_Matrix_chromatripsis_cells)
  rownames(MetaCell_Matrix_chromatripsis_cells) <-paste0("c", 1:nrow(MetaCell_Matrix_chromatripsis_cells))
  
  # , " (n=",cluster_counts , ")"
  
  melted_ct <-reshape2::melt(t(as.matrix(MetaCell_Matrix_chromatripsis_cells)), varnames=c('Chromosome', 'Clone'), as.is=T)
  melted_ct$value<-as.numeric(melted_ct$value)
  melted_ct$Chromosome<-as.character(melted_ct$Chromosome)
  melted_ct$Chromosome[which(melted_ct$Chromosome=="23")]="X"
  melted_ct$Chromosome[which(melted_ct$Chromosome=="24")]="Y"
  melted_ct$Clone<-as.character(melted_ct$Clone)
  melted_ct$value[melted_ct$value==0]<-NA
  clone_order<-paste0("c", 1:length(cluster_counts))
  melted_ct$Clone <-factor(melted_ct$Clone,levels = rev(clone_order))
  melted_ct$Chromosome <- factor(melted_ct$Chromosome, levels = c(1:22,"X","Y"))
  colnames(melted_ct)[3]<-"CT Fraction"
  
  melted_ct$Clone_mod<-melted_ct$Clone
  melted_ct$Clone_mod[which(melted_ct$Clone=="c1")]="c1"
  melted_ct$Clone_mod[which(melted_ct$Clone=="c6")]="c2"
  melted_ct$Clone_mod[which(melted_ct$Clone=="c2")]="c3"
  melted_ct$Clone_mod[which(melted_ct$Clone=="c4")]="c4"
  melted_ct$Clone_mod[which(melted_ct$Clone=="c5")]="c5"
  melted_ct$Clone_mod[which(melted_ct$Clone=="c3")]="c6"
  melted_ct$Clone_mod <-factor(melted_ct$Clone_mod,levels = paste0("c", rev(c(1, 2, 3, 4, 5, 6))))
  my_cols <- c("yellow", "red", "darkred")
  
  # Visualization
  plot<-ggpubr::ggballoonplot(melted_ct, size = "CT Fraction", fill="CT Fraction", color = "black", x="Chromosome", y = "Clone_mod")+
    scale_fill_gradientn(colors = my_cols)+
    guides(size = FALSE)+ theme(axis.title = element_text(size = 18, family="Arial"), axis.text = element_text(size = 18, family="Arial"),
                                                                        legend.title = element_text(size = 14, family="Arial"),  
                                legend.text =element_text(size = 14, angle = 45, family="Arial"), legend.position = "top")
  
  return(plot)
}

# Panel E for LFSMBP-Nuclei
plot<- plot_sample_CT("STP-Nuclei", 1)
ggsave(file=paste0("/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1FData/Fig1E_STP-Nuclei.svg"), plot, width=140*2, height=50*2, dpi=300, units = "mm")

# -----------------------
# Panel F - All Samples Summary Figure
# -----------------------

Samples<- c("STP-Nuclei", "STP-PDX", "ST1R-Nuclei", "ST1R-PDX", "MB243-Nuclei", "RCMB18-PDX")

library(ggpubr)

meansCT_samples<-c()
wmeansCT_samples<-c()
sdsCT_samples<-c()
wsdsCT_samples<-c()
for( sample in Samples)
{
  if(sample %in% c("STP-Nuclei", "STP-PDX", "ST1R-Nuclei", "ST1R-PDX"))
  {
    MetaCell_Matrix_chromatripsis_cells<-readRDS(paste0("/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/CT_Data/MetaCellsCT_",sample,"_",1,".RDS"))
  }else if(sample %in% c("MB243-Nuclei", "RCMB18-PDX"))
  {
    MetaCell_Matrix_chromatripsis_cells<-readRDS(paste0("/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/CT_Data/MetaCellsCT_",sample,"_",30,".RDS"))
    
  }else{
    MetaCell_Matrix_chromatripsis_cells<-readRDS(paste0("/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/CT_Data/MetaCellsCT_",sample,"_",50,".RDS"))
    
  }
  dim(MetaCell_Matrix_chromatripsis_cells)
  
  colnames(MetaCell_Matrix_chromatripsis_cells) <- c(1:22, "X", "Y")
  rownames(MetaCell_Matrix_chromatripsis_cells)<- paste("comm", 1:nrow(MetaCell_Matrix_chromatripsis_cells))
  
  # read clusters
  Clusters<-read.table(paste0("/home/r511a/projects/Gonzalo/Meduloblastoma/ScriptsPaper/ClusterTables/",sample,".txt"), header = T, stringsAsFactors = F)
  nrow(MetaCell_Matrix_chromatripsis_cells)==length(unique(Clusters$cluster))
  cluster_counts<-table(Clusters$cluster)
  colnames(MetaCell_Matrix_chromatripsis_cells) <- c(1:22, "X", "Y")
  weighted_means <-c()
  weighted_sds<-c()
  for(i in colnames(MetaCell_Matrix_chromatripsis_cells))
  {
    # weighted_means<- c(weighted_means, sum(MetaCell_Matrix_chromatripsis_cells[,i] * cluster_counts)/sum(cluster_counts))
    weighted_means<- c(weighted_means,  weighted.mean(MetaCell_Matrix_chromatripsis_cells[,i], cluster_counts/sum(cluster_counts)))
    weighted_sds<- c(weighted_sds,  radiant.data::weighted.sd(MetaCell_Matrix_chromatripsis_cells[,i], cluster_counts/sum(cluster_counts)))
  }
  
  
  
  wmeansCT_samples<-rbind(wmeansCT_samples, weighted_means)
  meansCT_samples<-rbind(meansCT_samples, colMeans(MetaCell_Matrix_chromatripsis_cells))
  wsdsCT_samples<-rbind(wsdsCT_samples, weighted_sds)
  sdsCT_samples<-rbind(sdsCT_samples, apply(MetaCell_Matrix_chromatripsis_cells, 2, sd))
}

rownames(meansCT_samples)<-Samples
meansCT_samples<-as.data.frame.matrix(meansCT_samples)

rownames(wmeansCT_samples)<-Samples
colnames(wmeansCT_samples) <- c(1:22, "X", "Y")
wmeansCT_samples<-as.data.frame.matrix(wmeansCT_samples)

rownames(sdsCT_samples)<-Samples
sdsCT_samples<-as.data.frame.matrix(sdsCT_samples)

rownames(wsdsCT_samples)<-Samples
colnames(wsdsCT_samples) <- c(1:22, "X", "Y")
wsdsCT_samples<-as.data.frame.matrix(wsdsCT_samples)

SamplesDataFrame<-c()
for(chr in colnames(meansCT_samples))
{
  for(sample in rownames(meansCT_samples))
  {
    SamplesDataFrame<-rbind(SamplesDataFrame, c(sample, chr, meansCT_samples[sample, chr], sdsCT_samples[sample, chr], wmeansCT_samples[sample, chr], wsdsCT_samples[sample, chr]))
  }
}

colnames(SamplesDataFrame)<-c("sample", "chr", "mean CT Fraction", "var CT Fraction", "wmean CT Fraction", "wvar CT Fraction")
SamplesDataFrame<-as.data.frame(SamplesDataFrame)
SamplesDataFrame$`mean CT Fraction`<-as.numeric(SamplesDataFrame$`mean CT Fraction`)
SamplesDataFrame$`wmean CT Fraction`<-as.numeric(SamplesDataFrame$`wmean CT Fraction`)
SamplesDataFrame$`var CT Fraction`<-as.numeric(SamplesDataFrame$`var CT Fraction`)
SamplesDataFrame$`wvar CT Fraction`<-as.numeric(SamplesDataFrame$`wvar CT Fraction`)
SamplesDataFrame$order<-1:length(SamplesDataFrame$sample)
reorder(chr,order)

# Visualization
SamplesDataFrame$chr <- factor(SamplesDataFrame$chr, levels = colnames(meansCT_samples))
SamplesDataFrame$`wmean CT Fraction`[which(SamplesDataFrame$`wmean CT Fraction`==0)]=NA
SamplesDataFrame$sample<- factor(SamplesDataFrame$sample,levels = rev(c("STP-Nuclei", "STP-PDX", "ST1R-Nuclei", "ST1R-PDX", "MB243-Nuclei", "RCMB18-PDX")))

my_cols <- c("white", "lightgreen", "chartreuse3", "darkgreen", "black")

summary<-ggpubr::ggballoonplot(SamplesDataFrame, x="chr", y="sample", size = "wmean CT Fraction", fill = "wvar CT Fraction")+
  scale_fill_gradientn(colors=my_cols) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),
  legend.title = element_text(size = 14),  legend.text =element_text(size = 14, angle = 45), legend.position = "top")

ggsave(file=paste0("/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1FData/Fig1E_all.svg"), summary, width=140*2, height=50*2, dpi=300, units = "mm")
