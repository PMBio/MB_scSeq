# ---------------------------------------------------------
# This scripts is used to generate subpanels in Fig1A
# Autor: R. Gonzalo Parra
# November 11, 2020
# ---------------------------------------------------------

getAnnotations<-function(gene_name)
{
  ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  biomaRt::listAttributes(ensembl, page="feature_page")
  
  annot <- biomaRt::getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                          filters = 'hgnc_symbol', 
                          values = gene_name, 
                          mart = ensembl)
  
  annot <- annot[which(annot[,"chromosome_name"] %in% c(1:22,"X","Y")),]
  
  unique_gene_ids<-which(!duplicated(annot[,"hgnc_symbol"]))
  annot <-annot[unique_gene_ids, ]
  rownames(annot) <-annot[,"hgnc_symbol"]
  
  return(annot)
}

modify_clone_colors<-function(plot_mod, chrs_num)
{
  library(grid)
  # Generate the ggplot2 plot grob
  g <- grid.force(ggplotGrob(plot_mod))
  # Get the names of grobs and their gPaths into a data.frame structure
  grobs_df <- do.call(cbind.data.frame, grid.ls(g, print = FALSE))
  # Build optimal gPaths that will be later used to identify grobs and edit them
  grobs_df$gPath_full <- paste(grobs_df$gPath, grobs_df$name, sep = "::")
  grobs_df$gPath_full <- gsub(pattern = "layout::", 
                              replacement = "", 
                              x = grobs_df$gPath_full, 
                              fixed = TRUE)
  
  # Get the gPaths of the strip background grobs
  strip_bg_gpath <- grobs_df$gPath_full[grepl(pattern = "strip\\.background.*", 
                                              x = grobs_df$gPath_full)]
  # Generate some color
  fills <- RColorBrewer::brewer.pal(n = 7, name = "Set1")
  aux<-fills[1]
  fills[1]<-fills[3]
  fills[3]<-aux
  
  for (i in (chrs_num+1):length(strip_bg_gpath)){
    g <- editGrob(grob = g, gPath = strip_bg_gpath[i], gp = gpar(fill = fills[i-chrs_num]))
  }
  # Draw the edited plot
  # grid.newpage(); grid.draw(g)
  return(g)
}

get_gene_segments<-function(gene_name, Segments_Matrix, segment_coords, segment_chromosomes)
{
  annot<-getAnnotations(gene_name)
  gene_chr<-annot$chromosome_name
  gene_start<-annot$start_position
  gene_end<-annot$end_position
  gene_segments<-c(which(as.numeric(segment_coords[,1]) <= gene_start & as.numeric(segment_coords[,2]) > gene_start & segment_chromosomes==gene_chr), which(as.numeric(segment_coords[,1]) <= gene_end & as.numeric(segment_coords[,2]) > gene_end & segment_chromosomes==gene_chr))
  return(rownames(Segments_Matrix)[gene_segments][1])
}


# Main script
Sample="STP-Nuclei"

# ----------------------
# We read the data
# ----------------------
Segments_Matrix <- read.table(paste0("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/DNA_processed_data_cellranger-dna_v.1.1.0/", Sample, "/cnv.cells.mat.20kb-bin.corrected.txt"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)

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

# ----------------------
# Plot de Full cells matrix
# ----------------------
# Segments_Matrix<-Segments_Matrix[,clusters[,"cell_id"]]
Segments_Matrix8 <-Segments_Matrix
Segments_colors <-Segments_Matrix

for(i in 1:dim(Segments_Matrix)[2])
{
  Segments_Matrix8[which(Segments_Matrix[,i] >8),i]<-8
}

# we plot the matrix
library(reshape2)
library(ggplot2)
library(ggrastr)
#
melted <-melt(as.matrix(Segments_Matrix8), varnames=c('Chromosome', 'Cells'), as.is=T)
head(melted)
#
melted$chrXY <- unlist(lapply(strsplit(melted$Chromosome, ':'), '[[',1))
melted$start <-unlist(lapply((strsplit(unlist(lapply(strsplit(melted$Chromosome, '-'), '[[',1)), ':')), '[[',2))
melted$end <- unlist(lapply(strsplit(melted$Chromosome, '-'), '[[',2))
melted$order<-1:length(melted$Chromosome)
melted$community<-Clusters[melted$Cells,"cluster"]
melted$community<-as.character(melted$community)
 
melted$clone<-paste0("c", melted$community)

MYCN<-get_gene_segments("MYCN", Segments_Matrix, segment_coords, segment_chromosomes)
TERT<-get_gene_segments("TERT", Segments_Matrix, segment_coords, segment_chromosomes)
HDAC9<-get_gene_segments("HDAC9", Segments_Matrix, segment_coords, segment_chromosomes)
CD274<-get_gene_segments("CD274", Segments_Matrix, segment_coords, segment_chromosomes)
SMO<-get_gene_segments("SMO", Segments_Matrix, segment_coords, segment_chromosomes)

midpoint=2
plotmatrix<- ggplot(melted) + rasterise(geom_tile(aes(x=reorder(Chromosome,order), y=reorder(Cells, order), fill=value)), dpi=300) +
  scale_fill_gradient2(low="darkblue", high="darkred",  midpoint = log1p(midpoint), mid = "lightgrey", trans = "log1p") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        panel.spacing=unit(.0001, "lines"), panel.border = element_rect(color = "black", fill = NA, size = 0.8), 
        strip.background=element_rect(colour="black"),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_blank(), strip.text = element_text(size = 18)) +
  facet_grid(factor(clone,levels = c("c1", "c2", "c3", "c4", "c5", "c6"))~reorder(chrXY,order), switch = "y", scales="free", space="free") + ggExtra::removeGrid() + labs(x="", y="") +
  geom_vline(xintercept = MYCN, col="red", linetype = "dotdash", lwd=1.1)+
  geom_vline(xintercept = TERT, col="blue", linetype = "dotdash", lwd=1.1)+
  geom_vline(xintercept = HDAC9, col="green", linetype = "dotdash", lwd=1.1)+
  geom_vline(xintercept = CD274, col="orange", linetype = "dotdash", lwd=1.1)+
  geom_vline(xintercept = SMO, col="violet", linetype = "dotdash", lwd=1.1)

# plotmatrix

# ggsave(file="/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1AData/middle.png", plotmatrix, width=160*4, height=90*2, dpi=300, units = "mm")
# ggsave(file="/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1AData/middle.svg", plotmatrix, width=160*4, height=90*2, dpi=300, units = "mm")

melted2<-melted[which(melted$chr==7),]
plotmatrix7<- ggplot(melted2) + rasterise(geom_tile(aes(x=reorder(Chromosome,order), y=Cells, fill=value)), dpi=300) +
  scale_fill_gradient2(low="darkblue", high="darkred",  midpoint = log1p(midpoint), mid = "lightgrey", trans = "log1p") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.spacing=unit(.0001, "lines"), panel.border = element_rect(color = "black", fill = NA, size = 0.8), strip.background=element_rect(colour="black"),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(),
  strip.text = element_text(size = 14)) +
  facet_grid(factor(clone,levels = c("c1", "c2", "c3", "c4", "c5", "c6"))~reorder(chrXY,order), switch = "y", scales="free", space="free") + ggExtra::removeGrid() + labs(x="", y="")+
  geom_vline(xintercept = HDAC9, col="green", linetype = "dotdash", lwd=1.1)+
  geom_vline(xintercept = SMO, col="violet", linetype = "dotdash", lwd=1.1)
plotmatrix7
ggsave(file="/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1AData/middle7.png", plotmatrix7, width=50*4, height=90*2, dpi=300, units = "mm")
ggsave(file="/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1AData/middle7.svg", plotmatrix7, width=50*4, height=90*2, dpi=300, units = "mm")

# ----------------------
# pseudobulk
# ---------------------- 
colors<-segment_chromosomes
colors[which(colors=="X")]=23
colors[which(colors=="Y")]=24

pseudobulk <- colMeans(t(Segments_Matrix))

colors_conf<-rep("darkblue", length(pseudobulk)) 

colors_scales<-rep("darkblue", length(pseudobulk)) 
colors_scales[which(pseudobulk>2.7)]<-"darkgreen"
colors_scales[which(pseudobulk<1.9)]<-"darkred"

pseudobulk.data <- data.frame(pseudobulk=pseudobulk, chr=segment_chromosomes, coordinates=1:length(pseudobulk))
pseudobulk.data$pseudobulk[which(pseudobulk.data$pseudobulk>8)]=8
pseudobulk.data$order<-1:length(pseudobulk.data$chr)

plotpseudobulk<- ggplot(pseudobulk.data, aes(coordinates, pseudobulk)) +
  geom_point(aes(colour = pseudobulk))+
  scale_colour_gradient2(low = "darkblue", high = "darkred", mid ="darkgray", midpoint = log1p(2), trans = "log1p" )+
  theme_bw()+theme(axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3), strip.background = element_blank(), strip.text = element_blank(), legend.position = "none")+
  facet_grid(.~reorder(chr,order), scales="free", space="free", switch = "both") + ggExtra::removeGrid() + labs(x="", y="")+
  scale_x_continuous(expand = c(0.01, 0.01))
# plotpseudobulk
ggsave(file="/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1AData/pbulk.png", plot=plotpseudobulk,width=160*4, height=25*2, dpi=300, units = "mm")

pseudobulk.data7<-pseudobulk.data[which(pseudobulk.data[,"chr"]==7),]
plotpseudobulk7<- ggplot(pseudobulk.data7, aes(coordinates, pseudobulk)) +
  geom_point(aes(colour = pseudobulk))+
  scale_colour_gradient2(low = "darkblue", high = "darkred", mid ="darkgray", midpoint = log1p(2), trans = "log1p" )+
  theme_bw()+theme(axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3),  strip.background = element_blank(), strip.text = element_blank(), legend.position = "none")+
  facet_grid(.~reorder(chr,order), scales="free", space="free", switch = "both") + ggExtra::removeGrid() + labs(x="", y="")+
  scale_x_continuous(expand = c(0.01, 0.01))
# plotpseudobulk7
# ggsave(file="/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1AData/pbulk7.png", plot=plotpseudobulk7, width=24, height=6, dpi=300, units = "cm")
ggsave(file="/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1AData/pbulk7.png", plot=plotpseudobulk7, width=50*4, height=25*2, dpi=300, units = "mm")

plotpseudobulk

rownames(Segments_Matrix)[which(bulk.data$chr==1 & bulk.data$bulk>4)]
pseudobulk.data[which(pseudobulk.data$chr==1 & pseudobulk.data$pseudobulk>6),]

rownames(Segments_Matrix)[which(pseudobulk.data$chr==6 & pseudobulk.data$pseudobulk>7)]

hist(rowMeans(Segments_Matrix[which(pseudobulk.data$chr==6 & pseudobulk.data$pseudobulk>6),]), xlab="Ploidy", main="", breaks = seq(1:100), xlim=c(0,30))

BulkTable_STP$CNV<-2**(log2(BulkData)) * 2

segment_start<-142600000
segment_end<-145220000
segment_chr<-1
  
BulkTable_STP[which(BulkTable_STP$start>=segment_start & BulkTable_STP$end<=segment_end &  BulkTable_STP$chr==segment_chr ),]

2**(log2(BulkData[which(BulkTable_STP$start>=segment_start & BulkTable_STP$end<=segment_end &  BulkTable_STP$chr==segment_chr )])) * 2

# ---------------------------
# Bulk 
# ---------------------------
# STP
BulkTable_STP<-read.table(file = "/icgc/dkfzlsdf/analysis/B060/share/B060_Stuttgart_case/cnv_CNmops/ST_LFS_1/XI046_ST_LFS_1_tumor/XI046_ST_LFS_1_tumor_control.bins.tsv", stringsAsFactors = F, header = T)
head(BulkTable_STP)


dim(BulkTable_STP)

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

# we calculate the tumor to control ratio
BulkData<-orderedSTP[,4]/orderedSTP[,5]

plot(2**(log2(BulkData)) * 2 )

xy<-which(BulkTable_STP[,1]=="23" | BulkTable_STP[,1]=="24")
notxy<-which(!(BulkTable_STP[,1]=="23" | BulkTable_STP[,1]=="24"))

BulkData[notxy] <- BulkData[notxy] * 2 
BulkData[xy]<- BulkData[xy] * 1 

bulk.data <- data.frame(bulk=BulkData, chr=BulkTable_STP[,1], coordinates=1:length(BulkData))
bulk.data$bulk[which(bulk.data$bulk>8)]=8
bulk.data$order<-1:length(bulk.data$chr)

bulk.data$chr[which(bulk.data$chr==23)]="X"
bulk.data$chr[which(bulk.data$chr==24)]="Y"

plotbulk<- ggplot(bulk.data, aes(coordinates, bulk)) +
  geom_point(aes(colour = bulk))+
  scale_colour_gradient2(low = "darkblue", high = "darkred", mid ="darkgray", midpoint = (2))+
  theme_bw()+theme(axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3), legend.position = "none", strip.background = element_blank(), strip.text = element_blank())+
  facet_grid(.~reorder(chr,order), scales="free", space="free") + ggExtra::removeGrid() + labs(x="", y="")+
  scale_x_continuous(expand = c(0.01, 0.01))
plotbulk
ggsave(file="/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1AData/bulk.png", plot=plotbulk, width=160*4, height=25*2, dpi=300, units = "mm")


bulk.data7<-bulk.data[which(bulk.data[,"chr"]==7),]
plotbulk7<- ggplot(bulk.data7, aes(coordinates, bulk)) +
  geom_point(aes(colour = bulk))+
  scale_colour_gradient2(low = "darkblue", high = "darkred", mid ="darkgray", midpoint = log1p(2), trans = "log1p" )+
  theme_bw()+theme(axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3), legend.position = "none", strip.background = element_blank(), strip.text = element_blank())+
  facet_grid(.~reorder(chr,order), scales="free", space="free") + ggExtra::removeGrid() + labs(x="", y="")+
  scale_x_continuous(expand = c(0.01, 0.01))
ggsave(file="/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1AData/bulk7.png", plot=plotbulk7, width=50*4, height=25*2, dpi=300, units = "mm")
