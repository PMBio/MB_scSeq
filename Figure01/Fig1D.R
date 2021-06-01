# ---------------------------------------------------------
# This scripts is used to generate subpanels in Fig1D
# Autor: R. Gonzalo Parra
# April 15, 2020
# ---------------------------------------------------------


pseudobulk_CT_scoring<-function(pseudo_cell, window_length, min_cnv_changes, min_consec_cnvs, segment_chromosomes)
{
  chromatriptic_cell<-FALSE
  chromosomes<- unique(segment_chromosomes)

  chromatriptic_chromosomes<-c()
  bins_windows<-rep(0, nrow(Segments_Matrix))
  chr_windows_n<-rep(0, length(chromosomes))
  names(chr_windows_n)<-chromosomes
  
  for (chromosome_i in chromosomes)
  {
    # start debug
    chromatriptic_chromosome <-0
    chr_windows<-0
    segments_in_chromosome_i=which(segment_chromosomes==chromosome_i)
    
    for (segment_i in segments_in_chromosome_i)
    {
      
      window_limit <- segments_in_chromosome_i[which(as.numeric(segment_coords[segments_in_chromosome_i,2]) - as.numeric(segment_coords[segment_i, 1]) > window_length)[1]]
      
      if(!is.na(window_limit)) # this check is for the case we reached the border condition
      {
        window_limit_length=as.numeric(segment_coords[window_limit,2]) - as.numeric(segment_coords[segment_i, 1])
        # This is the segment where the changes need to be quantified.
        if(segment_i>window_limit)
        {
          print(paste("shit ",segment_i, " ",window_limit))
        }
        # CNV_changes <- Segments_Matrix[segment_i:window_limit, cell_i]
        CNV_changes <-pseudo_cell[segment_i:window_limit]

        names(CNV_changes) <- NULL
        
        # We slide the window and count again until the end. 
        compressed_cnv <- rle(CNV_changes)
        window_changes<- length(compressed_cnv$values)
        cnv_states<-unique(compressed_cnv$values)
        
        # Clean vector
        # we ignore the windows that are too short
        shortwindows<-which(compressed_cnv$lengths<=min_consec_cnvs)
        if(length(shortwindows)>0)
        {
          window_changes <- window_changes - length(shortwindows)
        }
        #end clean vector
        chr_windows <- chr_windows +1 
        # if(window_changes>= (window_limit_length*min_cnv_changes)/window_length )
        if(window_changes>= min_cnv_changes)
        {
          chromatriptic_cell=TRUE
          chromatriptic_chromosome=chromatriptic_chromosome+1
          bins_windows[segment_i:window_limit]<-bins_windows[segment_i:window_limit] +1
        }
      }
    }
    # End debug
    if(chromatriptic_chromosome >0)
    {
      chromatriptic_chromosome <- chromatriptic_chromosome/chr_windows
    }
    chromatriptic_chromosomes=c(chromatriptic_chromosomes, chromatriptic_chromosome)
    chr_windows_n[chromosome_i]<-chr_windows
  }
  
  print(chromatriptic_chromosomes)
  return(list(chromatriptic_chromosomes, bins_windows, chr_windows_n))
}


# ---------------------------------------

Sample="STP-Nuclei"

# we use this matrix because it has less processing. 
Segments_Matrix <- read.table(paste0("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/DNA_processed_data_cellranger-dna_v.1.1.0/", Sample, "/cnv.cells.mat.20kb-bin.corrected.txt"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)

dim(Segments_Matrix)

clusters <-read.table(paste0("/home/r511a/projects/Gonzalo/Meduloblastoma/ScriptsPaper/ClusterTables/",Sample,".txt"), sep = "\t", stringsAsFactors = F, header = T, )
clusters[,"cell_id"]<-paste0("X",clusters[,"cell_id"])

Segments_Matrix<-Segments_Matrix[,clusters[,"cell_id"]]
dim(Segments_Matrix)

# we get the coordinates of each (row/segment)
splitted_ids <- strsplit(rownames(Segments_Matrix), split=":", fixed = FALSE, perl = FALSE, useBytes = FALSE)
splitted_ids <- matrix(unlist(splitted_ids), ncol=2, nrow=length(rownames(Segments_Matrix)), byrow = T)
segment_chromosomes <- splitted_ids[,1]
segment_coords <- splitted_ids[,2]

splitted_ids <- strsplit(segment_coords,  split="-", fixed = FALSE, perl = FALSE, useBytes = FALSE)
segment_coords <- matrix(unlist(splitted_ids), ncol=2, nrow=length(segment_coords), byrow = T)

chromosomes<- unique(segment_chromosomes)

# ---- We calculate the metaCells for each clone

Metacells<-c()

for(cluster_i in unique(clusters[,"cluster"]))
{
  print(paste0("creting metacell for clone ",  cluster_i, "..."))
  metacell<- apply(as.matrix(Segments_Matrix[,which(clusters[,"cluster"]==cluster_i)]), 1, function(x) as.numeric(names(which.max(table(as.numeric(x))))))
  Metacells<-rbind(Metacells, metacell)
}
rownames(Metacells)<-paste0("C", 1:6)
colnames(Metacells)<-rownames(Segments_Matrix)
# saveRDS(Metacells, "/home/r511a/projects/Gonzalo/Meduloblastoma/ScriptsPaper/MetaCellData/STP_Nuclei_Metacells.RDS")
# Metacells<-readRDS("/home/r511a/projects/Gonzalo/Meduloblastoma/ScriptsPaper/MetaCellData/STP_Nuclei_Metacells.RDS")

# chromosomes<-c(1:22, "X", "Y")

metacells_ct<-c()
windows_per_chromosomes<-c()
metacells_ct_bins<-c()

chromatriptic_chromosomes<-pseudobulk_CT_scoring(Metacells[metacell_i,], window_length, min_cnv_changes, min_consec_cnvs, segment_chromosomes)

window_length <-50000000
min_cnv_changes=10
min_consec_cnvs <- 1

for(metacell_i in 1:6)
{
  print(paste0("Calculating CT scores for metacell ",  cluster_i, "..."))
  chromatriptic_chromosomes<-pseudobulk_CT_scoring(Metacells[metacell_i,], window_length, min_cnv_changes, min_consec_cnvs, segment_chromosomes)
  metacells_ct<-rbind(metacells_ct, chromatriptic_chromosomes[[1]])
  metacells_ct_bins<-rbind(metacells_ct_bins, chromatriptic_chromosomes[[2]])
  windows_per_chromosomes<-rbind(windows_per_chromosomes, chromatriptic_chromosomes[[3]])
}

# Matrix_chromatripsis_cells_metacells <- pseudobulk_CT_scoring(pseudo_cell, window_length=50000000, min_cnv_changes=10, min_consec_cnvs=1, segment_chromosomes)
# saveRDS(Matrix_chromatripsis_cells_metacells, "/home/r511a/projects/Gonzalo/Meduloblastoma/ScriptsPaper/MetaCellData/STP_Nuclei_Metacells_CT.RDS")
# Matrix_chromatripsis_cells_metacells<-readRDS("/home/r511a/projects/Gonzalo/Meduloblastoma/ScriptsPaper/MetaCellData/STP_Nuclei_Metacells_CT.RDS")

CT_Metacells_new<-metacells_ct_bins
# Matrix_chromatripsis_cells_metacells[[3]]
rownames(CT_Metacells_new)<-rownames(Metacells)

# plot(CT_Metacells_new["C6",  segment_chromosomes==5]/Matrix_chromatripsis_cells_metacells[[3]][5], pch=16, ylab="CNV", xlab="Genomic Coordinates", main="", col=col1)

# -------------------------------------
# -------Plot All Clones --------------
# -------------------------------------
MetaCell_All_df<-c()

plots<-list()
p_i=1

# We have to reorder the cluster order to follow the order from the clone lineage tree
ClusterAssignments<-read.table("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/DNA_processed_data_cellranger-dna_v.1.1.0/sample_color_mapping.txt", header=F, stringsAsFactors = F, skip = 1, sep="\t", comment.char = "")
ClustersMapp<-ClusterAssignments[which(ClusterAssignments$V5==Sample), c(1,2)]
ClusterOrder<-ClustersMapp[,1]

library(ggplot2)
for(ClonePlot in paste0("C", ClusterOrder))
{
  MetaCell_df<-as.data.frame(cbind(rep(ClonePlot, length(segment_chromosomes)), colnames(Metacells)))
  colnames(MetaCell_df)<-c("Clone", "Chromosome")
  
  MetaCell_df$chrXY <- unlist(lapply(strsplit(MetaCell_df$Chromosome, ':'), '[[',1))
  MetaCell_df$start <-unlist(lapply((strsplit(unlist(lapply(strsplit(MetaCell_df$Chromosome, '-'), '[[',1)), ':')), '[[',2))
  MetaCell_df$end <- unlist(lapply(strsplit(MetaCell_df$Chromosome, '-'), '[[',2))
  MetaCell_df$CNV<-Metacells[ClonePlot,]
  MetaCell_df$CT_bin<-CT_Metacells_new[ClonePlot,]
  MetaCell_df$CT_bin_norm<-CT_Metacells_new[ClonePlot,]/windows_per_chromosomes[1,][MetaCell_df$chrXY]
  MetaCell_df$Coordinates<-1:length(MetaCell_df$CNV)
  
  MetaCell_df$CT_bin_norm[1] <-0
  MetaCell_df$CT_bin_norm[2]<-0.74
  
  MetaCell_df_filtered<-MetaCell_df
  # MetaCell_df_filtered<-MetaCell_df[MetaCell_df$chrXY %in% c(4,5,7,16,17,19,"X"),]
  MetaCell_df_filtered$CT_bin_norm[1] <-0
  MetaCell_df_filtered$CT_bin_norm[2]<-0.74
  MetaCell_df_filtered$CNV<-as.numeric(MetaCell_df_filtered$CNV)
  MetaCell_df_filtered$CNV[MetaCell_df_filtered$CNV>=10]<-10
  # MetaCell_df_filtered
  
  PlotClone<-ggplot(MetaCell_df_filtered, aes(Coordinates, CNV)) +
    ggrastr::rasterise(geom_point(aes(colour = (CT_bin_norm))), dpi=300)+
    scale_colour_gradient2(low = "yellow", high = "black", mid ="red", midpoint = (0.25))+  ylim(0, 10) + 
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))+
    facet_grid(.~reorder(chrXY,Coordinates), scales="free", space="free") + ggExtra::removeGrid() +
    scale_x_continuous(expand = c(0.01, 0.01))
  
  plots[[p_i]]<-PlotClone
  p_i=p_i+1
  
  # MetaCell_All_df<-rbind(MetaCell_All_df, MetaCell_df)
}

AllClones<-gridExtra::grid.arrange(grobs = plots, nrow = 6)
plot(AllClones)

ggsave(file="/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1AData/Fig1C_New.pdf", AllClones, width=210*4, height=120*2, dpi=300, units = "mm")
ggsave(file="/icgc/dkfzlsdf/analysis/B260/projects/Gonzalo/Meduloblastoma/ScriptsPaper/Fig1AData/Fig1C_New.png", AllClones, width=210*4, height=120*2, dpi=300, units = "mm")

