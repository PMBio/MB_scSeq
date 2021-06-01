# ---------------------------------------------------------
# This scripts is used to generate subpanels in Fig1G
# Autor: R. Gonzalo Parra
# April 18, 2020
# ---------------------------------------------------------

BulkPaths<-c()
BulkPaths<-c(BulkPaths, "/icgc/dkfzlsdf/analysis/B060/share/B060_Stuttgart_case/cnv_CNmops/ST_LFS_1/XI046_ST_LFS_1_tumor/XI046_ST_LFS_1_tumor_control.bins.tsv")
BulkPaths<-c(BulkPaths, "/icgc/dkfzlsdf/analysis/B060/share/B060_Stuttgart_case/cnv_CNmops/ST_LFS_1/XI046_ST_LFS_1_tumor/XI046_ST_LFS_1_tumor_control.bins.tsv")
BulkPaths<-c(BulkPaths, "/icgc/dkfzlsdf/analysis/B060/share/B060_Stuttgart_case/cnv_CNmops/ST_LFS_1/XI046_ST_LFS_1_relapse01/XI046_ST_LFS_1_relapse01_control.bins.tsv")
BulkPaths<-c(BulkPaths, "/icgc/dkfzlsdf/analysis/B060/share/B060_Stuttgart_case/cnv_CNmops/ST_LFS_1/XI046_ST_LFS_1_relapse01/XI046_ST_LFS_1_relapse01_control.bins.tsv")

BulkPaths<-c(BulkPaths, "/icgc/dkfzlsdf/analysis/B060/share/B060_Stuttgart_case/cnv_CNmops/ICGC_MB/ICGC_MB243_tumor_control.bins.tsv")
BulkPaths<-c(BulkPaths, "/icgc/dkfzlsdf/analysis/B060/share/B060_Stuttgart_case/cnv_CNmops/ICGC_MB/RCMB18.bins.tsv")
BulkPaths<-c(BulkPaths, "/icgc/dkfzlsdf/analysis/B060/share/B060_Stuttgart_case/cnv_CNmops/ICGC_MB/4M67.bins.tsv")

names(BulkPaths)<-c("STP-Nuclei", "STP-PDX", "ST1R-Nuclei", "ST1R-PDX", "MB243-Nuclei", "RCMB18-PDX", "4M67-Nuclei")


getAnnotations<-function(gene_name)
{
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  listAttributes(ensembl, page="feature_page")
  
  annot <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                 filters = 'hgnc_symbol', 
                 values = gene_name, 
                 mart = ensembl)
  
  annot <- annot[which(annot[,"chromosome_name"] %in% c(1:22,"X","Y")),]
  
  unique_gene_ids<-which(!duplicated(annot[,"hgnc_symbol"]))
  annot <-annot[unique_gene_ids, ]
  rownames(annot) <-annot[,"hgnc_symbol"]
  
  return(annot)
}

# ----------------------------------
# CNV  Violin Plots - Panel G
# ----------------------------------
# Read Single Cell Data --------
# ----------------------------------
# # Samples: "STP-Nuclei", "STP-PDX", "ST1R-Nuclei", "ST1R-PDX", "MB243-Nuclei", "RCMB18-PDX", "4M67-Nuclei"

plots<-list()
p_i=1

for(Sample in c("STP-Nuclei", "ST1R-Nuclei"))
{
Segments_Matrix <- read.table(paste0("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/DNA_processed_data_cellranger-dna_v.1.1.0/", Sample, "/cnv.cells.mat.20kb-bin.corrected.txt"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
dim(Segments_Matrix)

ClusterAssignments<-read.table("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/DNA_processed_data_cellranger-dna_v.1.1.0/sample_color_mapping.txt", header=F, stringsAsFactors = F, skip = 1, sep="\t", comment.char = "")
fills<-ClusterAssignments[which(ClusterAssignments$V5==Sample), 11]

ClustersMapp<-ClusterAssignments[which(ClusterAssignments$V5==Sample), c(1,2)]
rownames(ClustersMapp)<-ClustersMapp[,1]
ClustersMapp[,2]<-gsub("cluster_", "", ClustersMapp[,2])
rownames(ClustersMapp)<-ClustersMapp[,1]

Clusters<-read.table(paste0("/home/r511a/projects/Gonzalo/Meduloblastoma/ScriptsPaper/ClusterTables/",Sample,".txt"), header = T, stringsAsFactors = F)
Clusters<-Clusters[, c("cell_id", "cluster")]

# Reassign clusters in needed
Clusters[,2]<-ClustersMapp[as.character(Clusters[,2]),2]

cluster_counts<-table(Clusters$cluster)
Clusters[,"cell_id"]<-paste0("X",Clusters[,"cell_id"])
rownames(Clusters)<-Clusters[,"cell_id"]

Segments_Matrix<-Segments_Matrix[, Clusters[,"cell_id"]]
dim(Segments_Matrix)

# we get the coordinates of each (row/segment)
splitted_ids <- strsplit(rownames(Segments_Matrix), split=":", fixed = FALSE, perl = FALSE, useBytes = FALSE)
splitted_ids <- matrix(unlist(splitted_ids), ncol=2, nrow=length(rownames(Segments_Matrix)), byrow = T)
segment_chromosomes <- splitted_ids[,1]
segment_coords <- splitted_ids[,2]

splitted_ids <- strsplit(segment_coords,  split="-", fixed = FALSE, perl = FALSE, useBytes = FALSE)
segment_coords <- matrix(unlist(splitted_ids), ncol=2, nrow=length(segment_coords), byrow = T)

chromosomes<- unique(segment_chromosomes)

# ------------

genes_test<-c("MYCN",  "SMO", "CD274")
library(biomaRt)
BulkTable_STP<-read.table(file = BulkPaths[Sample], stringsAsFactors = F, header = T)
head(BulkTable_STP)

BulkTable_STP_chr<-BulkTable_STP[,"chr"]

BulkTable_STP[which(BulkTable_STP[,1]=="X"),1]<-23
BulkTable_STP[which(BulkTable_STP[,1]=="Y"),1]<-24
unique(BulkTable_STP[,1])

# Non aggregated data
# order data per chromosome (it was alphabetically ordered in the stored data)
orderedSTP<-c()
for(i in 1:24)
{
  orderedSTP <- rbind(orderedSTP, BulkTable_STP[which(BulkTable_STP[,1]==i),])
}
# we transform the data to get approx cnv states
BulkData<- 2**(orderedSTP[,6])

BulkData <- 2**(BulkData)

violin_df_all<-c()
# violin_df_bulk<-c()
violin_df_bulk_all<-c()

for (gene_name in genes_test)
{
  violin_df<-c()
  violin_df_bulk<-c()
  # gene_name<-"SMO"
  annot<-getAnnotations(gene_name)
  
  if(!is.null(annot))
  {
    
  gene_chr<-annot$chromosome_name
  gene_start<-annot$start_position
  gene_end<-annot$end_position
  
  # take the segments within the gene is contained
  gene_segments<-c(which(as.numeric(segment_coords[,1]) <= gene_start & as.numeric(segment_coords[,2]) >= gene_start & segment_chromosomes==gene_chr), which(as.numeric(segment_coords[,1]) <= gene_end & as.numeric(segment_coords[,2]) >= gene_end & segment_chromosomes==gene_chr))
  gene_segments<-unique(gene_segments)
  
  if(length(gene_segments)>0){
    if(length(gene_segments)==1)
    {
      gene_segments<-c(gene_segments, gene_segments)
    }
    
  gene_profile<-colMeans(Segments_Matrix[gene_segments[1]:gene_segments[-1],])

  violin_df<-cbind(Clusters[names(gene_profile)[which(names(gene_profile) %in% Clusters$cell_id)],], gene_profile[names(gene_profile)[which(names(gene_profile) %in% Clusters$cell_id)]])
  violin_df<-cbind(rep(gene_name, nrow(violin_df)), violin_df)
  
  # BULK 
  gene_bulk<-which(BulkTable_STP$chr==gene_chr & BulkTable_STP$start >= gene_start & BulkTable_STP$start <= gene_end)
  
  # if the gene lies in the middle of two segments. 
  if(length(gene_bulk)==0)
  {
    gene_index<-which(BulkTable_STP$chr==gene_chr & BulkTable_STP$start >= gene_start)[1]
    gene_bulk<-c(gene_index-1, gene_index)
  }
  
  
  cnv_gene_bulk<-mean(BulkData[gene_bulk])
  # violin_df<- rbind(violin_df, c(gene_name, paste0("bulk_", gene_name), "7", cnv_gene_bulk))
  violin_df_bulk<- rbind(violin_df_bulk, c(gene_name, paste0("bulk_", gene_name), "bulk", cnv_gene_bulk))
  
  colnames(violin_df)<-c("gene", "cell_id", "cluster", "cnv")
  
  violin_df<-as.data.frame(violin_df)
  violin_df$cluster<-as.character(violin_df$cluster)
  violin_df$gene<-as.character(violin_df$gene)
  violin_df$cnv<-as.numeric(violin_df$cnv)
  
  # Add median normalization by chr median
  cell_median_ploidy<-apply(Segments_Matrix[which(segment_chromosomes==gene_chr),], 2, median)
  violin_df$median_chr<-cell_median_ploidy[violin_df$cell_id]
  violin_df$cnv_corrected <- as.numeric(violin_df$cnv/cell_median_ploidy[violin_df$cell_id])
  violin_df$log2_cnv_corrected<-log2(violin_df$cnv_corrected)

  
  colnames(violin_df_bulk)<-c("gene", "cell_id", "cluster", "cnv")
  violin_df_bulk<-as.data.frame(violin_df_bulk)
  violin_df_bulk$median_chr<-median(BulkData[which(BulkTable_STP_chr==gene_chr)])
  violin_df_bulk$cnv_corrected<-as.numeric(as.numeric(violin_df_bulk$cnv)/ median(BulkData[which(BulkTable_STP_chr==gene_chr)]) )
  violin_df_bulk$log2_cnv_corrected<-log2(violin_df_bulk$cnv_corrected)
  
  violin_df_all<-rbind(violin_df_all, violin_df)
  violin_df_bulk_all<-rbind(violin_df_bulk_all, violin_df_bulk)
  }
  }
}

library(ggplot2)
violin_df_bulk_all$cnv<-as.numeric(violin_df_bulk_all$cnv)
violin_df_bulk_all$cnv_corrected<-as.numeric(violin_df_bulk_all$cnv_corrected)

violin_df_all$cnv<- as.numeric(violin_df_all$cnv)
violin_df_all$cluster<- as.character(violin_df_all$cluster)
violin_df_all$cnv_corrected<-as.numeric(violin_df_all$cnv_corrected)


dodge <- position_dodge(width = 1)
dodge2<- position_dodge(width = 0.99)

theme_new <- theme_set(theme_bw())
theme_new <- theme_update(
panel.border = element_rect(color = "black", fill = NA, size = 0.8), strip.background=element_rect(colour="black"),  
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing=unit(.0001, "lines"), 
axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_text(size = 18, face = "italic", family = "Arial"), 
axis.title = element_text(size = 18), axis.text.y = element_text(size = 18), 
legend.title = element_text(size = 18),  legend.text =element_text(size = 18),  legend.position = "top")

g_ymax=max(c(1, max(violin_df_all$log2_cnv_corrected), max(violin_df_bulk$log2_cnv_corrected)))
g_ymin=max(c(-1, min(violin_df_all$log2_cnv_corrected), min(violin_df_bulk$log2_cnv_corrected)))

drug_genes<-ggplot(violin_df_all, aes(x=gene, y=log2(cnv_corrected), fill=cluster)) +
# drug_genes<-ggplot(violin_df_all, aes(x=gene, y=cnv, fill=cluster)) +
  # geom_dotplot(position = dodge2, binaxis='y', stackdir='center', stackratio=0.15, dotsize=1, color=NA)+
  geom_boxplot(position = dodge, width=0.9)+
  ylim(-1, g_ymax)+
  scale_fill_manual(values = scales::alpha(c(fills, "red"),1))+
  facet_grid(.~gene, scales="free", space="free") + facet_wrap(~gene, ncol=4, scales = "free")+
  geom_point(data = violin_df_bulk_all, col="red", size=4, pch=17)
  # theme(legend.position = "none")

plots[[p_i]]<-drug_genes
p_i=p_i+1
# plot(drug_genes)

}

targets<-gridExtra::grid.arrange(grobs = plots, nrow = 2)
plot(targets)
ggsave(file=paste0("/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1FData/Fig1F.svg"), plot=targets, width=160*2, height=30*4, dpi=300, units = "mm")

# ----- Expression 
# druggable.gene.expression.df <- read.table("/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA/druggableGenes_allSamples_expression.txt", header = T, sep = "\t")
druggable.gene.expression.df <- read.table("/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA/druggableGenes_allSamples_expression_w_zeros.txt", header = T, sep = "\t")

# druggable.gene.expression.df[which(druggable.gene.expression.df$Genes=="MYCN" & druggable.gene.expression.df$Clone=="Clone5" & druggable.gene.expression.df$Expression>3),]

dodge <- position_dodge(width = 1)
dodge2<- position_dodge(width = 0.99)

head(druggable.gene.expression.df)

plots_expression_clones<-list()
plots_expression_types<-list()

fills <- RColorBrewer::brewer.pal(n = 7, name = "Set1")
aux<-fills[1]
fills[1]<-fills[3]
fills[3]<-aux

p_i=1

STP_expression_drugabble_genes_all<-c()

for(Gene_plot in c("MYCN", "SMO"))
{
  STP_expression_drugabble_genes<-druggable.gene.expression.df[which(druggable.gene.expression.df$Sample=="STP-Nuclei" & druggable.gene.expression.df$Genes==Gene_plot),]
  unique(STP_expression_drugabble_genes$Clone)
  
  STP_expression_drugabble_genes<-STP_expression_drugabble_genes[which(STP_expression_drugabble_genes$Clone %in% c("purkinje cells","Clone2", "Clone3", "Clone4", "Clone5", "Clone6")),]
  
  STP_expression_drugabble_genes$Clone[which(STP_expression_drugabble_genes$Clone=="purkinje cells")]="Ref"
  STP_expression_drugabble_genes$Clone[which(STP_expression_drugabble_genes$Clone=="Clone2")]="C2"
  STP_expression_drugabble_genes$Clone[which(STP_expression_drugabble_genes$Clone=="Clone3")]="C3"
  STP_expression_drugabble_genes$Clone[which(STP_expression_drugabble_genes$Clone=="Clone4")]="C4"
  STP_expression_drugabble_genes$Clone[which(STP_expression_drugabble_genes$Clone=="Clone5")]="C5"
  STP_expression_drugabble_genes$Clone[which(STP_expression_drugabble_genes$Clone=="Clone6")]="C6"
  
  STP_expression_drugabble_genes<-STP_expression_drugabble_genes[which(!(STP_expression_drugabble_genes$Expression==0)),]
  
  
  CellTypes<-c("Ref", "C2", "C3", "C4", "C5", "C6")
  
  Missing<-CellTypes[which(!CellTypes %in% unique(STP_expression_drugabble_genes$Clone))]
  
  if(length(Missing==1))
  {
    Missing<-c(Missing, Missing)
  }
  
  for(aux in Missing)
  {
    STP_expression_drugabble_genes<-rbind(STP_expression_drugabble_genes, c(aux, Gene_plot, "0", aux, "-", "STP-Nuclei"))
  }
  
  STP_expression_drugabble_genes$Clone <-factor(STP_expression_drugabble_genes$Clone,levels = CellTypes)
  STP_expression_drugabble_genes$Expression<-as.numeric(STP_expression_drugabble_genes$Expression)
  
  p_clone<-ggplot(STP_expression_drugabble_genes, aes(x=Clone, y=Expression, fill=Clone)) + 
    geom_boxplot(position = dodge, width=0.9)+
    stat_summary(fun.y="mean", color="red", shape=15)+
    scale_fill_manual(values = scales::alpha(fills,1))+
    # geom_dotplot(position = dodge2, binaxis='y', stackdir='center', stackratio=0.15, dotsize=0.6, color=NA)+
    facet_grid(.~Genes, scales="free", space="free")+theme(legend.position = "none")
  
  plots_expression_clones[[p_i]]<-p_clone
  plots_expression_types[[p_i]]<-p_type
  
  p_i=p_i+1
  
  STP_expression_drugabble_genes_all<-rbind(STP_expression_drugabble_genes_all, STP_expression_drugabble_genes)
}

targets_clones<-gridExtra::grid.arrange(grobs = plots_expression_clones, nrow = 1)
# targets_types<-gridExtra::grid.arrange(grobs = plots_expression_types, nrow = 3)
# plot(targets_clones)

ggsave(file=paste0("/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1FData/Fig1F2.svg"), plot=targets_clones,  width=120*2, height=15*4, dpi=300, units = "mm")

