#############################################################################################################################
##                                                                                                                      
##  PERFORM DRUGGABLE TARGETS SUBCLONAL CNV/CT ANALYSIS
##                                                                                                                      
##  Date: 20 DECEMBER 2023                                                                                                                   
##  
##  Author: Moritz Przybilla
##
##                                                                                                                      
############################################################################################################################

# clear workspace
rm(list=ls()) 
set.seed(16011985) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "RColorBrewer", "ggplot2", "biomaRt", "httr", "data.table", 
                      "dplyr", "Biobase","tidyverse", "ggsci", "cowplot", "ggrepel", "memoise")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))
httr::set_config(httr::config(ssl_verifypeer=0L))

`%nin%` = Negate(`%in%`)

############################################################################
##                              FUNCTIONS
############################################################################

readDNATable1 <- function(Sample){
  if (Sample == "MB243-Nuclei" | Sample == "RCMB18-PDX"){
    
    Segments_Matrix <- read.table(paste0("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones/", str_split_fixed(Sample, "-", 2)[,1], "_hmmcopy_20kb_logR.csv"), sep = ",", stringsAsFactors = F, header = T, row.names = 1)
    Clusters<-read.table(paste0("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones/", str_split_fixed(Sample, "-", 2)[,1], "_cell_clone_assignments.csv"), sep = ",", header = T, stringsAsFactors = F)
    
  } else {
    
    Segments_Matrix <- read.table(paste0("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones/", gsub("-", "_", Sample), "_hmmcopy_20kb_logR.csv"), sep = ",", stringsAsFactors = F, header = T, row.names = 1)
    Clusters<-read.table(paste0("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones/", gsub("-", "_", Sample), "_cell_clone_assignments.csv"), sep = ",", header = T, stringsAsFactors = F)
  }
  
  rownames(Segments_Matrix) <- paste0(Segments_Matrix$chr, ":", Segments_Matrix$start, "-", Segments_Matrix$end)
  return(list(Segments_Matrix, Clusters))
}

cm <- cachem::cache_mem(max_size = 30 * 1024^3)

readDNATable <- memoise(readDNATable1, cache = cm)

# gene_name <- druggable.targets$V1
getAnnotations<-function(gene_name){
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

# Sample = "STP-PDX"
# gene_name = "PD-L1"
# FUNCTION FOR VISUALISING COPY NUMBER IN SCDNA
# Sample <- sample.tmp
# gene_name <- gene
plot_genes<-function(Sample, gene_name, gene_annot)
{
  
  # READ IN THE LOGR MATRIX
  # Segments_Matrix <- read.table(paste0("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones/MB243_hmmcopy_20kb_logR.csv"), sep = ",", stringsAsFactors = F, header = T, row.names = 1)
  
  tmp <- readDNATable(Sample)
  Segments_Matrix <- tmp[[1]]
  Clusters <- tmp[[2]]
  dim(Segments_Matrix)
  
  # READ IN THE CLUSTERING FILE
  # Clusters<-read.table(paste0("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones/MB243_cell_clone_assignments.csv"), sep = ",", header = T, stringsAsFactors = F)
  # Clusters<-read.table(paste0("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones/", str_split_fixed(Sample, "-", 2)[,1], "_cell_clone_assignments.csv"), sep = ",", header = T, stringsAsFactors = F)
  Clusters <- Clusters[, c("Clone", "cells")]
  Clusters <- Clusters[Clusters$Clone != "Not Assigned", ]
  colnames(Clusters)<- c("cluster", "barcode")
  cluster_counts<-table(Clusters$cluster)
  rownames(Clusters)<-Clusters[,"barcode"]
  
  # READ IN COLOURS
  ClusterAssignments<-read.table("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/DNA_processed_data_cellranger-dna_v.1.1.0/sample_color_mapping.txt", header=F, stringsAsFactors = F, skip = 1, sep="\t", comment.char = "")
  fills<-ClusterAssignments[which(ClusterAssignments$V5==Sample), 11]
  
  # SUBSET MATRIX TO ASSIGNED CELLS
  Segments_Matrix<-Segments_Matrix[, Clusters[,"barcode"]]
  dim(Segments_Matrix)
  
  # we get the coordinates of each (row/segment)
  splitted_ids <- strsplit(rownames(Segments_Matrix), split=":", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  splitted_ids <- matrix(unlist(splitted_ids), ncol=2, nrow=length(rownames(Segments_Matrix)), byrow = T)
  segment_chromosomes <- splitted_ids[,1]
  segment_coords <- splitted_ids[,2]
  
  splitted_ids <- strsplit(segment_coords,  split="-", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  segment_coords <- matrix(unlist(splitted_ids), ncol=2, nrow=length(segment_coords), byrow = T)
  chromosomes<- unique(segment_chromosomes)
  
  # READ IN THE BULK DATA AND WRANGLE INTO SHAPE
  # BulkTable_STP<-read.table(file = "/omics/groups/OE0014/internal/share/B060_Stuttgart_case/cnv_CNmops/ICGC_MB/ICGC_MB243_tumor_control.bins.tsv", stringsAsFactors = F, header = T)
  BulkTable_STP<-fread(file = BulkPaths[Sample], stringsAsFactors = F, header = T)
  head(BulkTable_STP)
  
  BulkTable_STP_chr<-BulkTable_STP[,"chr", with = FALSE]
  BulkTable_STP[which(BulkTable_STP[,1]=="X"),1]<-23
  BulkTable_STP[which(BulkTable_STP[,1]=="Y"),1]<-24
  # unique(BulkTable_STP[,1])
  
  BulkTable_STP[,chr := factor(chr, levels = 1:24, ordered=TRUE)]
  
  # Non aggregated data
  # order data per chromosome (it was alphabetically ordered in the stored data)
  # orderedSTP<-c()
  # for(i in 1:24)
  # {
  #   orderedSTP <- rbind(orderedSTP, BulkTable_STP[which(BulkTable_STP[,1]==i),])
  # }
  # 
  
  orderedSTP <- as.data.frame(BulkTable_STP[order(chr, start)])
  # we transform the data to get approx cnv states
  BulkData<- 2**(orderedSTP[,6])
  BulkData <- 2**(BulkData)
  
  # add correction for MB243
  if (Sample == "MB243-Nuclei"){
    
    for(clone in paste0("Clone_", 1:4)){
      cur.cells <- Clusters[Clusters$cluster == clone, "barcode"]
      if(clone=="Clone_1"){
        Segments_Matrix[, cur.cells] <- 2^(Segments_Matrix[, cur.cells]) * 3
      } else {
        Segments_Matrix[, cur.cells] <- 2^Segments_Matrix[, cur.cells] * 4
      }
    } 
    
  }
  
  # create list for storing information
  violin_df_all<-c()
  violin_df_bulk_all<-c()
  violin_df<-c()
  violin_df_bulk<-c()
  
  # get gene annotation
  annot<-gene_annot
  
  if(!is.null(annot))
  {
    
    gene_chr<-annot$chr
    gene_start<-annot$start
    gene_end<-annot$end
    
    # take the segments within the gene is contained
    gene_segments<-c(which(as.numeric(segment_coords[,1]) <= gene_start & as.numeric(segment_coords[,2]) >= gene_start & segment_chromosomes==gene_chr), which(as.numeric(segment_coords[,1]) <= gene_end & as.numeric(segment_coords[,2]) >= gene_end & segment_chromosomes==gene_chr))
    gene_segments<-unique(gene_segments)
    
    if(length(gene_segments)>0){
      if(length(gene_segments)==1)
      {
        gene_segments<-c(gene_segments, gene_segments)
      }
      
      gene_profile<-colMeans(Segments_Matrix[gene_segments[1]:gene_segments[-1],])
      
      violin_df<-cbind(Clusters[names(gene_profile)[which(names(gene_profile) %in% Clusters$barcode)],], gene_profile[names(gene_profile)[which(names(gene_profile) %in% Clusters$barcode)]])
      violin_df<-cbind(rep(gene_name, nrow(violin_df)), violin_df)
      colnames(violin_df) <- c("gene", "cluster", "cell_id", "cnv")
      
      # BULK 
      gene_bulk<-which(BulkTable_STP$chr==gene_chr & BulkTable_STP$start >= gene_start & BulkTable_STP$start <= gene_end)
      
      # if the gene lies in the middle of two segments. 
      if(length(gene_bulk)==0)
      {
        gene_index<-which(BulkTable_STP$chr==gene_chr & BulkTable_STP$start >= gene_start)[1]
        gene_bulk<-c(gene_index-1, gene_index)
      }
      
      
      cnv_gene_bulk<-mean(BulkData[gene_bulk])
      violin_df_bulk<- rbind(violin_df_bulk, c(gene_name, paste0("bulk_", gene_name), "bulk", cnv_gene_bulk))
      
      # change data types
      violin_df<-as.data.frame(violin_df)
      violin_df$cluster<-as.character(violin_df$cluster)
      violin_df$gene<-as.character(violin_df$gene)
      violin_df$cnv<-as.numeric(violin_df$cnv)
      
      # Add median normalization by cell median
      cell_median_ploidy<-apply(Segments_Matrix[which(!segment_chromosomes%in%c("X","Y")),], 2, median, na.rm = T)
      violin_df$median_chr<- cell_median_ploidy[violin_df$cell_id]
      violin_df$cnv_corrected <- as.numeric(violin_df$cnv/cell_median_ploidy[violin_df$cell_id])
      violin_df$log2_cnv_corrected <- log2(violin_df$cnv_corrected)
      
      # implement the same calculations for the bulk
      colnames(violin_df_bulk)<-c("gene", "cell_id", "cluster", "cnv")
      violin_df_bulk<-as.data.frame(violin_df_bulk)
      violin_df_bulk$median_chr<-median(BulkData[which(!BulkTable_STP_chr%in%c("X", "Y"))])
      violin_df_bulk$cnv_corrected<-as.numeric(as.numeric(violin_df_bulk$cnv)/ median(BulkData[which(BulkTable_STP_chr==gene_chr)]) )
      violin_df_bulk$log2_cnv_corrected<-log2(violin_df_bulk$cnv_corrected)
      
      violin_df_all<-rbind(violin_df_all, violin_df)
      violin_df_bulk_all<-rbind(violin_df_bulk_all, violin_df_bulk)
    }
  }
  
  # set the plotting parameters
  violin_df_bulk_all$cnv<-as.numeric(violin_df_bulk_all$cnv)
  violin_df_bulk_all$cnv_corrected<-as.numeric(violin_df_bulk_all$cnv_corrected)
  
  violin_df_all$cnv<- as.numeric(violin_df_all$cnv)
  violin_df_all$cluster<- as.character(violin_df_all$cluster)
  violin_df_all$cnv_corrected<-as.numeric(violin_df_all$cnv_corrected)
  violin_df_all[is.na(violin_df_all$log2_cnv_corrected), "log2_cnv_corrected"] <- 0
  
  dodge <- position_dodge(width = 1)
  dodge2<- position_dodge(width = 0.99)
  
  theme_new <- theme_set(theme_bw())
  theme_new <- theme_update(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8), strip.background=element_rect(colour="black"),  
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing=unit(.0001, "lines"), 
    axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_text(size = 18, face = "italic"), 
    axis.text.y = element_text(size = 18), axis.title = element_blank(), 
    legend.title = element_text(size = 18),  legend.text =element_text(size = 18),  legend.position = "top")
  
  g_ymax=max(c(1, max(violin_df_all$log2_cnv_corrected), max(violin_df_bulk$log2_cnv_corrected)))
  g_ymin=max(c(-1, min(c(min(violin_df_all$log2_cnv_corrected), min(violin_df_bulk$log2_cnv_corrected)))))
  
  colnames(violin_df_all)[2]<-"Clone"
  colnames(violin_df_bulk)[2]<-"Clone"
  
  if(Sample=="STP-Nuclei")
  {
    Sample="LFS-MBP_Nuclei"
    
  } else if (Sample=="STP-PDX") {
    
    Sample="LFS-MBP_PDX"
  }
  
  fills <- c(brewer.pal(n=length(c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5")),"Set1"))
  names(fills) <- c("Clone_1", "Clone_2", "Clone_3", "Clone_4", "Clone_5")
  
  # drug_genes <- ggplot(violin_df_all, aes(x=Clone, y=cnv_corrected, fill=Clone)) + 
  #   geom_dotplot(position = dodge2, binaxis='y', stackdir='center',
  #                stackratio=0.75, dotsize=0.5, color = "black") +
  #   geom_violin(position = dodge2, color = "black", width = 0.5, fill = NA) +
  #   scale_color_manual(values = scales::alpha(fills,1), "Clone")+
  #   scale_fill_manual(values = scales::alpha(c(fills, "red"),1))+
  #   facet_grid(.~gene, scales="free", space="free") + # ggtitle(Sample)+
  #   geom_point(data = violin_df_bulk, col="red", size=4, pch=17) + theme(legend.position = "right")
  
  drug_genes <- ggplot(violin_df_all, aes(x=gene, y=log2_cnv_corrected, fill=Clone)) + 
    geom_dotplot(position = dodge2, binaxis='y', stackdir='center', stackratio=0.01, dotsize=0.5, color=NA)+
    geom_boxplot(position = dodge, width=0.5, outlier.color = NA)+
    scale_x_discrete(drop = FALSE) +
    scale_color_manual(values = scales::alpha(fills,1), "Clone")+
    scale_fill_manual(values = scales::alpha(c(fills, "red"),1))+
    ylim(g_ymin, g_ymax) +
    facet_grid(.~gene, scales="free", space="free") + # ggtitle(Sample)+
    geom_point(data = violin_df_bulk, col="red", size=4, pch=17) + theme(legend.position = "right")
  
  return(drug_genes)
  
}

# FUNCTION FOR VISUALISING COPY NUMBER IN SCRNA
# Sample <- "MB243-Nuclei"
# gene_name <- "ABCE1"
expression_druggable_genes <- fread("infercnv/revision/scrna_analysis/infercnv_MB/scRNA_scDNA/druggableGenes_allSamples_expression_w_zeros.txt", header = T, data.table = F)
# expression_druggable_genes <- fread("infercnv/revision/scrna_analysis/infercnv_MB/scRNA_scDNA/allGenes_allSamples_expression_w_zeros.txt", header = T, data.table = F)


plot_gene_expression<-function(Sample, gene_name){
  
  # READ IN THE DATA
  # expression_druggable_genes <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_MB/scRNA_scDNA/druggableGenes_allSamples_expression_w_zeros.txt", header = T, data.table = F)
  # expression_druggable_genes <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_MB/scRNA_scDNA/selectedGenes_allSamples_expression_w_zeros.txt", header = T, data.table = F)
  # expression_druggable_genes <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_MB/scRNA_scDNA/druggableGenes_allSamples_inferCNV.txt", header = T, data.table = F)
  sample_druggable_genes <- expression_druggable_genes[which(expression_druggable_genes$Sample==Sample & expression_druggable_genes$Genes==gene_name),]
  
  if(nrow(sample_druggable_genes) > 0){
    
    sample_druggable_genes$Clone <-factor(sample_druggable_genes$Clone,levels = paste0("Clone", 1:length(unique(sample_druggable_genes$Clone))))
    
    withzeros <- table(sample_druggable_genes$Clone)
    withzeros<-withzeros[paste0("Clone", 1:length(unique(sample_druggable_genes$Clone)))]
    
    sample_druggable_genes.withoutzero <- sample_druggable_genes[which(!(sample_druggable_genes$Expression==0)),]
    
    withoutzeros <- table(sample_druggable_genes.withoutzero$Clone)
    withoutzeros <- withoutzeros[paste0("Clone", 1:length(unique(sample_druggable_genes$Clone)))]
    withoutzeros[which(is.na(withoutzeros))] <- 0
    
    counts_genes<-c()
    for(i in 1:length(withzeros))
    {
      counts_genes <- c(counts_genes, paste0(withoutzeros[i],"/",withzeros[i]))
    }
    
    sample_druggable_genes$Expression <- as.numeric(sample_druggable_genes$Expression)
    sample_druggable_genes <- sample_druggable_genes[sample_druggable_genes$Expression > 0,]
    sample_druggable_genes$Merged_Clone <- factor(sample_druggable_genes$Merged_Clone, levels = c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5"))
    
    if(Sample=="STP-Nuclei")
    {
      Sample="LFS-MBP_Nuclei"
    } else if (Sample=="STP-PDX") {
      
      Sample="LFS-MBP_PDX"
    }
    
    dodge <- position_dodge(width = 1)
    dodge2<- position_dodge(width = 0.99)
    
    theme_new <- theme_set(theme_bw())
    theme_new <- theme_update(
      panel.border = element_rect(color = "black", fill = NA, size = 0.8), strip.background=element_rect(colour="black"),  
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing=unit(.0001, "lines"), 
      axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text.x = element_text(size = 18, face = "italic"), 
      strip.text.y = element_text(size = 8),
      axis.text.y = element_text(size = 18), axis.title = element_blank(), 
      legend.title = element_text(size = 18),  legend.text =element_text(size = 18),  legend.position = "right")
    
    fills <- c(brewer.pal(n=length(c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5")),"Set1"))
    names(fills) <- c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5")
    
    g_ymax=max(c(1, max(sample_druggable_genes$Expression)))+0.1
    # g_ymin=0
    g_ymin = min(min(c(1, min(sample_druggable_genes$Expression)))-0.1)
    # browser()
    
    merged_clone_tables <- table(sample_druggable_genes$Merged_Clone, sample_druggable_genes$Clone)
    merged_clone_names <- names(which(rowSums(merged_clone_tables>0)>1))
    sample_druggable_genes$Merged_Clone <- as.character(sample_druggable_genes$Merged_Clone)
    for(merged_clone in merged_clone_names){
      collapsing_clones <- unique(sample_druggable_genes[sample_druggable_genes$Merged_Clone==merged_clone,"Clone"])
      names(fills)[names(fills)==merged_clone] <- paste(collapsing_clones, collapse=" + ")
      sample_druggable_genes[sample_druggable_genes$Merged_Clone==merged_clone,"Merged_Clone"] <- 
        paste(collapsing_clones, collapse=" + ")
    }
    
    drug_genes <- ggplot(sample_druggable_genes, aes(x=Merged_Clone, y=Expression, fill=Merged_Clone)) + 
      geom_dotplot(position = dodge2, binaxis='y', stackdir='center', stackratio=0.005, dotsize=0.5, color=NA)+
      geom_boxplot(position = dodge, width=0.35, outlier.color = NA) +
      ylim(g_ymin, g_ymax) +
      scale_x_discrete(drop = FALSE) +
      scale_color_manual(values = scales::alpha(fills,1), "Clone") +
      scale_fill_manual(values = scales::alpha(c(fills),1)) +
      facet_grid(. ~ Genes, scales="free", space="free")
    
    return(drug_genes)
    
  } else {
    
    print(paste0(gene, " is not expressed in Sample ", Sample))
  }
  
}

############################################################################
##                          READ IN THE DATA
############################################################################

# READ IN DRUGGABLE TARGET LIST
# druggable.targets <- read.table("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/druggable_targets_list.txt", header = F)

# READ IN DEG RESULTS FROM FF ANALYSIS
ff.deg.results <- read.table("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/MB_CT_vs_NCT_DESeq2_TCC_GROUP_wo_Immune.txt", header = T)
ff.deg.results <- ff.deg.results[,c("log2FoldChange", "pvalue", "padj", "gene_name")]
ff.deg.results$Type <- "FF"

# SUBSET TO GENES IN DRUGGABLE TARGET LIST
# ff.drug.vector <- ff.deg.results[ff.deg.results$gene_name %in% druggable.targets$V1,]
ff.drug.vector <- ff.deg.results
ff.drug.vector <- ff.drug.vector[ff.drug.vector$padj < 0.1, ]
ff.drug.vector <- ff.drug.vector[order(ff.drug.vector$log2FoldChange, decreasing = T),]

# READ IN DEG RESULTS FROM FFPE ANALYSIS
ffpe.deg.results <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/MB_TP53_vs_WT_DESeq2_wo_Immune.txt", header = T)
ffpe.deg.results <- ffpe.deg.results[,c("log2FoldChange", "pvalue", "padj", "gene_name")]
ffpe.deg.results$Type <- "FFPE"
ffpe.deg.results <- as.data.frame(ffpe.deg.results)

# SUBSET TO GENES IN DRUGGABLE TARGET LIST
# ffpe.drug.vector <- ffpe.deg.results[ffpe.deg.results$gene_name %in% druggable.targets$V1,]
ffpe.drug.vector <- ffpe.deg.results
ffpe.drug.vector <- ffpe.drug.vector[ffpe.drug.vector$padj < 0.1, ]
ffpe.drug.vector <- ffpe.drug.vector[order(ffpe.drug.vector$log2FoldChange, decreasing = T),]

# CREATE LIST OF BULK PATHS WITH SAMPLE NAMES
BulkPaths<-c()
BulkPaths<-c(BulkPaths, "/omics/groups/OE0014/internal/share/B060_Stuttgart_case/cnv_CNmops/ST_LFS_1/XI046_ST_LFS_1_tumor/XI046_ST_LFS_1_tumor_control.bins.tsv")
BulkPaths<-c(BulkPaths, "/omics/groups/OE0014/internal/share/B060_Stuttgart_case/cnv_CNmops/ST_LFS_1/XI046_ST_LFS_1_tumor/XI046_ST_LFS_1_tumor_control.bins.tsv")
BulkPaths<-c(BulkPaths, "/omics/groups/OE0014/internal/share/B060_Stuttgart_case/cnv_CNmops/ST_LFS_1/XI046_ST_LFS_1_relapse01/XI046_ST_LFS_1_relapse01_control.bins.tsv")
BulkPaths<-c(BulkPaths, "/omics/groups/OE0014/internal/share/B060_Stuttgart_case/cnv_CNmops/ST_LFS_1/XI046_ST_LFS_1_relapse01/XI046_ST_LFS_1_relapse01_control.bins.tsv")
BulkPaths<-c(BulkPaths, "/omics/groups/OE0014/internal/share/B060_Stuttgart_case/cnv_CNmops/ICGC_MB/ICGC_MB243_tumor_control.bins.tsv")
BulkPaths<-c(BulkPaths, "/omics/groups/OE0014/internal/share/B060_Stuttgart_case/cnv_CNmops/ICGC_MB/RCMB18.bins.tsv")
names(BulkPaths)<-c("STP-Nuclei", "STP-PDX", "ST1R-Nuclei", "ST1R-PDX", "MB243-Nuclei", "RCMB18-PDX")

# READ IN THE LIST OF ALTERED REGIONS PER CLONE
ct.altered.regions <- list.files("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones", pattern = "Altered_Regions.csv", all.files = T, full.names = T)

# READ IN DRUGGABLE TARGETS WHICH ARE EXPRESSED
# expression_druggable_genes <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_MB/scRNA_scDNA/druggableGenes_allSamples_expression_w_zeros.txt", header = T, data.table = F)
# expression_druggable_genes <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_MB/scRNA_scDNA/selectedGenes_allSamples_expression_w_zeros.txt", header = T, data.table = F)
# expression_druggable_genes <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_MB/scRNA_scDNA/druggableGenes_allSamples_inferCNV.txt", header = T, data.table = F)
# expression_druggable_genes <- fread("infercnv/revision/scrna_analysis/infercnv_MB/scRNA_scDNA/selectedGenes_allSamples_inferCNV.txt", header = T, data.table = F)
expressed.druggable.targets <- unique(expression_druggable_genes$Genes)

############################################################################
##                    VISUALISE THE TARGETS FROM BULK
############################################################################

# # STP-PDX - FF
# pdf("infercnv/Figures_PaperSTP-PDX_bulk_FF_druggable_targets_scDNA_scRNA.pdf", width = 10, height = 3)
# 
# i <- 2
# for(i in 1:nrow(ff.drug.vector)){
#   
#   gene <- ff.drug.vector[i,"gene_name"]
#   
#   gene.scDNA.plot <- plot_genes("STP-PDX", gene)
#   gene.scRNA.plot <- plot_gene_expression("STP-PDX", gene)
#   print(gene.scDNA.plot)
#   print(gene.scRNA.plot)
#   
# }
# 
# dev.off()
# 
# 
# # STP-PDX - FFPE
# pdf("infercnv/Figures_PaperSTP-PDX_bulk_FFPE_druggable_targets_scDNA_scRNA.pdf", width = 10, height = 3)
# 
# i <- 1
# for(i in 1:nrow(ffpe.drug.vector)){
#   
#   gene <- ffpe.drug.vector[i,"gene_name"]
#   
#   gene.scDNA.plot <- plot_genes("STP-PDX", gene)
#   gene.scRNA.plot <- plot_gene_expression("STP-PDX", gene)
#   print(gene.scDNA.plot)
#   print(gene.scRNA.plot)
#   
# }
# 
# dev.off()

############################################################################
##              GET THE TARGETS OVERLAPPING CT OR CNV REGIONS
############################################################################

# GET LOCATIONS FOR DRUGGABLE TARGETS
# drug.target.annotations <- getAnnotations(druggable.targets$V1)
drug.target.annotations <- getAnnotations(unique(expressed.druggable.targets))
colnames(drug.target.annotations) <- c("gene_name", "chr", "start", "end")
# drug.target.annotations <- drug.target.annotations[drug.target.annotations$gene_name %in% expressed.druggable.targets, ]
drug.target.GRange <- makeGRangesFromDataFrame(drug.target.annotations, keep.extra.columns = T)

# GET THE FILES AND WRANGLE INTO SHAPE
sample.ids <- str_split_fixed(basename(ct.altered.regions), "_Altered", 2)[,1]
ct.altered.region.data <- lapply(ct.altered.regions, read.table, sep = ",", header = T)
subclonal.altered.regions <- list()

j <- 1
for (j in 1:length(sample.ids)){
  
  sample.tmp <- sample.ids[j]
  print(sample.tmp)
  altered.region.data.tmp <- ct.altered.region.data[[j]]
  altered.region.data.tmp$sampleID <- sample.tmp
  
  # ONLY SELECT SUBCLONAL EVENTS
  region.data.tmp <- altered.region.data.tmp[grep("Subclonal|Clonal", altered.region.data.tmp$Status),]
  
  if(nrow(region.data.tmp) > 0){
    
    ct.region.GRange <- makeGRangesFromDataFrame(region.data.tmp[region.data.tmp$Status == "Subclonal CT Event" | region.data.tmp$Status == "Clonal CT Event",], keep.extra.columns = T)
    cnv.region.GRange <- makeGRangesFromDataFrame(region.data.tmp[region.data.tmp$Status == "Subclonal CNV" | region.data.tmp$Status == "Clonal CNV",], keep.extra.columns = T)
    
    # INTERSECT BY OVERLAPS FOR CNV AND CT RESPECTIVELY
    ct.overlap <- subsetByOverlaps(drug.target.GRange, ct.region.GRange)
    
    hits <- findOverlaps(drug.target.GRange, ct.region.GRange)
    status <- unique(CharacterList(split(ct.region.GRange$Status[subjectHits(hits)],
                                      queryHits(hits))))
    mcols(ct.overlap) <- DataFrame(mcols(ct.overlap), status)
    
    cnv.overlap <- subsetByOverlaps(drug.target.GRange, cnv.region.GRange)
    
    hits <- findOverlaps(drug.target.GRange, cnv.region.GRange)
    status <- unique(CharacterList(split(cnv.region.GRange$Status[subjectHits(hits)],
                                         queryHits(hits))))
    mcols(cnv.overlap) <- DataFrame(mcols(cnv.overlap), status)
    
    # MAKE DATAFRAME FROM IT KEEPING THE INFORMATION FOR THE SAMPLE
    ct.overlap.data <- as.data.frame(ct.overlap)
    ct.overlap.data$status <- as.character(ct.overlap.data$status)
    cnv.overlap.data <- as.data.frame(cnv.overlap)
    cnv.overlap.data$status <- as.character(cnv.overlap.data$status)
    altered.drug.targets <- rbind(ct.overlap.data, cnv.overlap.data)
    
    # STORE
    altered.drug.targets$sampleID <- sample.tmp

    
  } else {
    
    print(paste0(sample.tmp, " does not have any events."))
    next
  }
  
}

############################################################################
##            VISUALISE PIECHARTS PER SAMPLE ACROSS ALL TARGETS
############################################################################

# LIST TO STORE INFO
all.targets <- list()

j <- 1
for (j in 1:length(sample.ids)){
  
  sample.tmp <- sample.ids[j]
  print(sample.tmp)
  altered.region.data.tmp <- ct.altered.region.data[[j]]
  altered.region.data.tmp$sampleID <- sample.tmp
  altered.region.data.tmp <- altered.region.data.tmp[altered.region.data.tmp$Status != "Subclonal Low Confidence CT Event", ]
  
  region.GRange <- makeGRangesFromDataFrame(altered.region.data.tmp, keep.extra.columns = T)
  
  # INTERSECT BY OVERLAPS FOR CNV AND CT RESPECTIVELY
  overlap <- subsetByOverlaps(drug.target.GRange, region.GRange)
  
  hits <- findOverlaps(drug.target.GRange, region.GRange)
  status <- unique(CharacterList(split(region.GRange$Status[subjectHits(hits)],
                                       queryHits(hits))))
  mcols(overlap) <- DataFrame(mcols(overlap), status)
  
  # MAKE DATAFRAME FROM IT KEEPING THE INFORMATION FOR THE SAMPLE
  overlap.data <- as.data.frame(overlap)
  overlap.data$status <- as.character(as.vector(overlap$status))
  
  # STORE
  overlap.data$sampleID <- sample.tmp
  all.targets[[sample.tmp]] <- overlap.data

  # VISUALISE RESULTS AS A PIECHART
  target.count <- data.frame(table(overlap.data$status))
  print(sum(target.count$Freq))
  target.count$Var1 <- as.character(target.count$Var1)
  target.count <- target.count[order(target.count$Freq, decreasing = T),]
  target.count$Var1 <- factor(target.count$Var1, levels = c(target.count[order(target.count$Freq, decreasing = T), "Var1"]))
  
  # colours for the patients
  # nb.cols <- length(unique(target.count$Var1))
  mycolors <- c(
    "Not Altered" = "lightgrey",
    "Clonal CNV" = "#33741a",
    "Subclonal CNV" = "#6ba82f",
    "Subclonal CT Event" = "#de77af",
    "Clonal CT Event" = "#d01c8b",
    "Subclonal Low Confidence CT Event" = "gray50"
  )
  
  # Compute the position of labels
  pie.chart.data <- target.count %>% 
    arrange(desc(Var1)) %>%
    mutate(prop = Freq / sum(target.count$Freq)*100) %>%
    mutate(labels = paste0(round(prop,2), "%")) %>% 
    mutate(ypos = cumsum(prop)- 0.5*prop )
  print(pie.chart.data)
  # Basic piechart
  piechart <- ggplot(pie.chart.data, aes(x="", y=prop, fill=Var1)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="none",
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1)) +
    geom_label_repel(data = pie.chart.data %>% filter(prop > 0.5),
                     aes(label = paste0(labels), y = ypos), 
                     nudge_x = 0.6, nudge_y = 0.6, color = "white", fontface = "bold",
                     size = 7, show.legend = F) +
    scale_fill_manual(values = mycolors)
  
  pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_druggable_targets_scDNA_piechart.pdf"), width = 4, height = 4)
  # pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_all_genes_scDNA_piechart.pdf"), width = 4, height = 4)
  print(piechart)
  dev.off()
}

############################################################################
##              VISUALISE SUBCLONALLY ALTERED DRUG TARGETS
############################################################################

# GET THE FILES AND WRANGLE INTO SHAPE
sample.ids <- str_split_fixed(basename(ct.altered.regions), "_", 2)[,1]
ct.altered.region.data <- lapply(ct.altered.regions, read.table, sep = ",", header = T)
subclonal.altered.regions <- list()

j <- 1
for (j in 1:length(sample.ids)){
  
  sample.tmp <- sample.ids[j]
  print(sample.tmp)
  altered.region.data.tmp <- ct.altered.region.data[[j]]
  altered.region.data.tmp$sampleID <- sample.tmp
  
  # ONLY SELECT SUBCLONAL EVENTS
  subclonal.region.data.tmp <- altered.region.data.tmp[grep("Subclonal", altered.region.data.tmp$Status),]
  
  if(nrow(subclonal.region.data.tmp) > 0){
    
    subclonal.ct.region.GRange <- makeGRangesFromDataFrame(subclonal.region.data.tmp[subclonal.region.data.tmp$Status == "Subclonal CT Event",], keep.extra.columns = T)
    subclonal.cnv.region.GRange <- makeGRangesFromDataFrame(subclonal.region.data.tmp[subclonal.region.data.tmp$Status == "Subclonal CNV",], keep.extra.columns = T)
    
    # INTERSECT BY OVERLAPS FOR CNV AND CT RESPECTIVELY
    subclonal.ct.overlap <- subsetByOverlaps(drug.target.GRange, subclonal.ct.region.GRange)
    subclonal.cnv.overlap <- subsetByOverlaps(drug.target.GRange, subclonal.cnv.region.GRange)
    
    # MAKE DATAFRAME FROM IT KEEPING THE INFORMATION FOR THE SAMPLE
    subclonal.ct.overlap.data <- as.data.frame(subclonal.ct.overlap)
    if(nrow(subclonal.ct.overlap.data)) subclonal.ct.overlap.data$type <- "CT"
    subclonal.cnv.overlap.data <- as.data.frame(subclonal.cnv.overlap)
    if(nrow(subclonal.cnv.overlap.data)) subclonal.cnv.overlap.data$type <- "CNV"
    subclonal.altered.drug.targets <- rbind(subclonal.ct.overlap.data, subclonal.cnv.overlap.data)
    # subclonal.altered.drug.targets <- rbind(subclonal.cnv.overlap.data)
    
    # STORE
    subclonal.altered.drug.targets$sampleID <- sample.tmp
    subclonal.altered.regions[[sample.tmp]] <- subclonal.altered.drug.targets
    
  } else {
    
    print(paste0(sample.tmp, " does not have any subclonal events."))
    next
  }

}

# COMBINE LISTS
all.subclonal.drug.targets <- bind_rows(subclonal.altered.regions)
all.subclonal.drug.targets[all.subclonal.drug.targets$sampleID == "STPNuc", "sampleID"] <- "STP-Nuclei"
all.subclonal.drug.targets[all.subclonal.drug.targets$sampleID == "MB243", "sampleID"] <- "MB243-Nuclei"
all.subclonal.drug.targets[all.subclonal.drug.targets$sampleID == "ST1R", "sampleID"] <- "ST1R-PDX"
all.subclonal.drug.targets[all.subclonal.drug.targets$sampleID == "STP", "sampleID"] <- "STP-PDX"
sample.ids <- unique(all.subclonal.drug.targets$sampleID)

# ITERATE OVER EACH SAMPLE AND PRINT THE TARGETS AFFECTED BY WHICH TYPE
k <- 1
for (k in 1:length(sample.ids)){
  
  sample.tmp <- sample.ids[k]
  print(sample.tmp)
  dir.create(paste0("infercnv/Figures_Paper/", sample.tmp))
  
  # GET TARGETS
  cnv.targets.tmp <- all.subclonal.drug.targets[all.subclonal.drug.targets$sampleID == sample.tmp & all.subclonal.drug.targets$type == "CNV", "gene_name"]
  ct.targets.tmp <- all.subclonal.drug.targets[all.subclonal.drug.targets$sampleID == sample.tmp & all.subclonal.drug.targets$type == "CT", "gene_name"]
  
  
  sample_druggable_genes <- expression_druggable_genes[which(expression_druggable_genes$Sample==sample.tmp),]
  
  # sample_druggable_genes <- data.table(sample_druggable_genes)
  
  # sample_druggable_genes[,mean(Expression),.(Merged_Clone, Genes)][,sd(V1), Genes][,V1]
  
  # SUBCLONAL CNV TARGETS
  pdf(paste0("infercnv/Figures_Paper/", sample.tmp,"/", sample.tmp, "_CNV_druggable_targets_scDNA_scRNA.pdf"), width = 6, height = 8)

  # pdf(paste0("infercnv/Figures_Paper/", sample.tmp,"/", sample.tmp, "_CNV_druggable_targets_scDNA_scRNA_split_by_celltype.pdf"), width = 6, height = 8)
  # pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "/", sample.tmp, "_CNV_top_100_genes_scDNA_scRNA.pdf"), width = 6, height = 2.5)
  
  # cnv.targets.tmp <- sample_druggable_genes[Genes %in% cnv.targets.tmp,mean(Expression),.(Merged_Clone, Genes)][,sd(V1), Genes][order(V1, decreasing = TRUE)][1:100,Genes]
  
  i <- 1
  for(i in seq_along(cnv.targets.tmp)){
    gene <- cnv.targets.tmp[i]
    if(all(gene %in% sample_druggable_genes$Genes)){
      gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
      scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
      print(scDNA.gene.plot) 
      scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)
      print(scRNA.gene.plot)
    } else {
      print(print(paste0(gene, " is not expressed in Sample ", sample.tmp)))
    }

  }
  
  dev.off()
  
  pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "/", sample.tmp, "_CT_druggable_targets_scDNA_scRNA.pdf"), width = 6, height = 8)
  
  # pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "/", sample.tmp, "_CT_druggable_targets_scDNA_scRNA_split_by_celltype.pdf"), width = 6, height = 8)
  # pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "/", sample.tmp, "_CT_top_100_genes_scDNA_scRNA.pdf"), width = 6, height = 2.5)
  
  # ct.targets.tmp <- sample_druggable_genes[Genes %in% ct.targets.tmp,mean(Expression),.(Merged_Clone, Genes)][,sd(V1), Genes][order(V1, decreasing = TRUE)][1:100,Genes]
  
  
  i <- 1
  for(i in seq_along(ct.targets.tmp)){
    gene <- ct.targets.tmp[i]
    if(all(gene %in% sample_druggable_genes$Genes)){
      gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
      scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
      print(scDNA.gene.plot) 
      scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)
      print(scRNA.gene.plot)
    } else {
      print(print(paste0(gene, " is not expressed in Sample ", sample.tmp)))
    }
    
  }
  
  dev.off()
  
  # # SUBCLONAL CT TARGETS
  # # pdf(paste0("infercnv/Figures_Paper", sample.tmp, "_CT_druggable_targets_scDNA_scRNA.pdf"), width = 10, height = 3)
  # pdf(paste0("infercnv/Figures_Paper", sample.tmp, "/", sample.tmp, "_CT_selectedGenes_scDNA_scRNA.pdf"), width = 6, height = 2.5)
  # 
  # for(i in 1:length(ct.targets.tmp)){
  #   
  #   gene <- ct.targets.tmp[i]
  #   scDNA.gene.plot <- plot_genes(sample.tmp, gene)
  #   print(scDNA.gene.plot) 
  #   scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)
  #   print(scRNA.gene.plot)
  #   
  # }
  # 
  # dev.off()
  # 
  # # FOR FF
  # cnv.ff.targets.tmp <- cnv.targets.tmp[cnv.targets.tmp %in% ff.drug.vector$gene_name]
  # 
  # if(length(cnv.ff.targets.tmp) > 0){
  #   
  #   # SUBCLONAL CNV TARGETS
  #   # pdf(paste0("infercnv/Figures_Paper", sample.tmp, "_CNV_FF_druggable_targets_scDNA_scRNA.pdf"), width = 10, height = 3)
  #   pdf(paste0("infercnv/Figures_Paper", sample.tmp, "/", sample.tmp, "_CNV_FF_selectedGenes_scDNA_scRNA.pdf"), width = 6, height = 2.5)
  #   
  #   for(i in 1:length(cnv.ff.targets.tmp)){
  #     
  #     gene <- cnv.ff.targets.tmp[i]
  #     scDNA.gene.plot <- plot_genes(sample.tmp, gene)
  #     print(scDNA.gene.plot) 
  #     scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)
  #     print(scRNA.gene.plot)
  #     
  #   }
  #   
  #   dev.off()
  #   
  # }
  # 
  # ct.ff.targets.tmp <- ct.targets.tmp[ct.targets.tmp %in% ff.drug.vector$gene_name]
  # 
  # if(length(ct.ff.targets.tmp) > 0){
  #   # SUBCLONAL CT TARGETS
  #   # pdf(paste0("infercnv/Figures_Paper", sample.tmp, "_CT_FF_druggable_targets_scDNA_scRNA.pdf"), width = 10, height = 3)
  #   pdf(paste0("infercnv/Figures_Paper", sample.tmp, "/", sample.tmp, "_CT_FF_selectedGenes_scDNA_scRNA.pdf"), width = 6, height = 2.5)
  #   
  #   for(i in 1:length(ct.ff.targets.tmp)){
  #     
  #     gene <- ct.ff.targets.tmp[i]
  #     scDNA.gene.plot <- plot_genes(sample.tmp, gene)
  #     print(scDNA.gene.plot) 
  #     scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)
  #     print(scRNA.gene.plot)
  #     
  #   }
  #   
  #   dev.off()
  # }
  # 
  # # FOR FFPE
  # cnv.ffpe.targets.tmp <- cnv.targets.tmp[cnv.targets.tmp %in% ffpe.drug.vector$gene_name]
  # 
  # if(length(cnv.ffpe.targets.tmp) > 0){
  #   
  #   # SUBCLONAL CNV TARGETS
  #   # pdf(paste0("infercnv/Figures_Paper", sample.tmp, "_CNV_FFPE_druggable_targets_scDNA_scRNA.pdf"), width = 10, height = 3)
  #   pdf(paste0("infercnv/Figures_Paper", sample.tmp, "/", sample.tmp, "_CNV_FFPE_selectedGenes_scDNA_scRNA.pdf"), width = 6, height = 2.5)
  #   
  #   for(i in 1:length(cnv.ffpe.targets.tmp)){
  #     
  #     gene <- cnv.ffpe.targets.tmp[i]
  #     scDNA.gene.plot <- plot_genes(sample.tmp, gene)
  #     print(scDNA.gene.plot) 
  #     scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)
  #     print(scRNA.gene.plot)
  #     
  #   }
  #   
  #   dev.off()
  #   
  # }
  # 
  # ct.ffpe.targets.tmp <- ct.targets.tmp[ct.targets.tmp %in% ffpe.drug.vector$gene_name]
  # 
  # if(length(ct.ffpe.targets.tmp) > 0){
  #   
  #   # SUBCLONAL CT TARGETS
  #   # pdf(paste0("infercnv/Figures_Paper", sample.tmp, "_CT_FFPE_druggable_targets_scDNA_scRNA.pdf"), width = 10, height = 3)
  #   pdf(paste0("infercnv/Figures_Paper", sample.tmp, "/", sample.tmp, "_CT_FFPE_selectedGenes_scDNA_scRNA.pdf"), width = 6, height = 2.5)
  #   
  #   for(i in 1:length(ct.ffpe.targets.tmp)){
  #     
  #     gene <- ct.ffpe.targets.tmp[i]
  #     scDNA.gene.plot <- plot_genes(sample.tmp, gene)
  #     print(scDNA.gene.plot) 
  #     scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)
  #     print(scRNA.gene.plot)
  #     
  #   }
  #   
  #   dev.off()
  # }
  
}

## Make dotplot per sample
library(ggnewscale)

returnSampleCNV <- function(sample){
  
  sample.cnv.table <- switch(sample,
    "MB243-Nuclei" = readRDS("~/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones/MB243_clonal_cnv.rds"),
    "ST1R-PDX" = readRDS("~/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones/ST1RPDX_clonal_cnv.rds"),
    "STP-PDX" = readRDS("~/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones/STPPDX_clonal_cnv.rds"),
    "STP-Nuclei" = readRDS("~/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/scAbsolute_clones/STPNucres_1MB_clonal_cnv.rds")
  )
  return(sample.cnv.table)
}

k <- 1
for (k in 1:length(sample.ids)){
  
  sample.tmp <- sample.ids[k]
  print(sample.tmp)
  dir.create(paste0("infercnv/Figures_Paper/", sample.tmp))
  
  # GET TARGETS
  cnv.targets.tmp <- all.subclonal.drug.targets[all.subclonal.drug.targets$sampleID == sample.tmp & all.subclonal.drug.targets$type == "CNV", "gene_name"]
  ct.targets.tmp <- all.subclonal.drug.targets[all.subclonal.drug.targets$sampleID == sample.tmp & all.subclonal.drug.targets$type == "CT", "gene_name"]
  

  sample.cnv.table <- returnSampleCNV(sample.tmp)
  sample.cnv.table.m <- melt(sample.cnv.table, id.vars = c("chr", "start", "end"))
  sample.cnv.table.m[[4]] <- gsub(pat="_", rep="", x=sample.cnv.table.m[[4]])
  colnames(sample.cnv.table.m)[4:5] <- c("Clone", "TCN")
  sample.cnv.table.m <- data.table(sample.cnv.table.m)
  
  sample_druggable_genes <- data.table(expression_druggable_genes[which(expression_druggable_genes$Sample==sample.tmp),])
  sample_druggable_genes <- merge(sample_druggable_genes, drug.target.annotations, by.x="Genes", by.y="gene_name")
  
  merged_clone_tables <- table(sample_druggable_genes$Merged_Clone, sample_druggable_genes$Clone)
  merged_clone_names <- names(which(rowSums(merged_clone_tables>0)>1))
  sample_druggable_genes$Merged_Clone <- as.character(sample_druggable_genes$Merged_Clone)
  for(merged_clone in merged_clone_names){
    collapsing_clones <- unique(sample_druggable_genes[sample_druggable_genes$Merged_Clone==merged_clone,Clone])
    sample_druggable_genes[sample_druggable_genes$Merged_Clone==merged_clone,"Merged_Clone"] <- 
      paste(collapsing_clones, collapse=" + ")
  }
  
  
  # sample_druggable_genes <- data.table(sample_druggable_genes)
  
  # sample_druggable_genes[,mean(Expression),.(Merged_Clone, Genes)][,sd(V1), Genes][,V1]
  
  cnv_druggable_genes <- sample_druggable_genes[Genes %in% cnv.targets.tmp]
  
  # SUBCLONAL CNV TARGETS
  pdf(paste0("infercnv/Figures_Paper/", sample.tmp,"/", sample.tmp, "_CNV_druggable_targets_scDNA_scRNA_dotplot_with_legend.pdf"), width = length(unique(cnv_druggable_genes$Genes))/10*6, height = 8)
  
  # pdf(paste0("infercnv/Figures_Paper/", sample.tmp,"/", sample.tmp, "_CNV_druggable_targets_scDNA_scRNA_split_by_celltype.pdf"), width = 6, height = 8)
  # pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "/", sample.tmp, "_CNV_top_100_genes_scDNA_scRNA.pdf"), width = 6, height = 2.5)
  
  # cnv.targets.tmp <- sample_druggable_genes[Genes %in% cnv.targets.tmp,mean(Expression),.(Merged_Clone, Genes)][,sd(V1), Genes][order(V1, decreasing = TRUE)][1:100,Genes]
  
  cnv_druggable_genes[,Pct := mean(Expression > 0), .(Genes, Merged_Clone)]
  cnv_druggable_genes[,Ave := mean(Expression), .(Genes, Merged_Clone)]
  cnv_druggable_genes <- unique(cnv_druggable_genes[,.(Genes, Merged_Clone, Pct, Ave, chr, start, end, Clone)])
  
  
  cnv_druggable_genes$TCN <- apply(cnv_druggable_genes, 1, function(x){
    cur.cnv.table <- sample.cnv.table.m[Clone==x[["Clone"]]]
    myx <- which(cur.cnv.table$start<as.numeric(x[["start"]]) & cur.cnv.table$end>as.numeric(x[["start"]]) & cur.cnv.table$chr == x[["chr"]]):
      which(cur.cnv.table$start<as.numeric(x[["end"]]) & cur.cnv.table$end>as.numeric(x[["end"]]) & cur.cnv.table$chr == x[["chr"]])
    return(median(cur.cnv.table[myx][["TCN"]]))
    
    })
  
  cnv_druggable_genes[,Mean_TCN := mean(TCN), .(Merged_Clone, Genes)]
  
  cnv_druggable_genes <- unique(cnv_druggable_genes[,.(Genes, Merged_Clone, Pct, Ave, chr, start, end, Mean_TCN)])
  
  p <- ggplot(cnv_druggable_genes, aes(x=Genes, y=Merged_Clone)) + 
    geom_tile(aes(fill=Mean_TCN)) + 
    scale_fill_gradient2(guide = guide_colorbar(ticks.colour = "black",
                                                frame.colour = "black"),
                         name = "Average Copy Number", high="#b2182b", low="#4393c3", mid="#F6F6F6", midpoint = 2, limits=c(0,8)) +
    new_scale_fill() + 
    geom_point(aes(size=Pct, fill=Ave), color="black", shape=21) +
    scale_size("% detected", range = c(0,6),limits = c(0,1)) +
    scale_fill_gradientn(colours = viridisLite::mako(100),
                         guide = guide_colorbar(ticks.colour = "black",
                                                frame.colour = "black"),
                         name = "Average\nexpression", limits = c(0,4)) +
    ylab("Cluster") + xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title = element_text(size=14)) + coord_equal()
  print(p)
  dev.off()
  
  
  p <- p + theme(legend.position = "none")
  
  pdf(paste0("infercnv/Figures_Paper/", sample.tmp,"/", sample.tmp, "_CNV_druggable_targets_scDNA_scRNA_dotplot.pdf"), width = length(unique(cnv_druggable_genes$Genes))/10*6, height = 8)
  print(p)
  dev.off()
  
  png(paste0("infercnv/Figures_Paper/", sample.tmp,"/", sample.tmp, "_CNV_druggable_targets_scDNA_scRNA_dotplot.png"), 
      width = length(unique(cnv_druggable_genes$Genes))/10*6, height = 8, unit="in", res=600)
  print(p)
  dev.off()
  
  
  # pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "/", sample.tmp, "_CT_druggable_targets_scDNA_scRNA_split_by_celltype.pdf"), width = 6, height = 8)
  # pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "/", sample.tmp, "_CT_top_100_genes_scDNA_scRNA.pdf"), width = 6, height = 2.5)
  
  # ct.targets.tmp <- sample_druggable_genes[Genes %in% ct.targets.tmp,mean(Expression),.(Merged_Clone, Genes)][,sd(V1), Genes][order(V1, decreasing = TRUE)][1:100,Genes]
  
  ct_druggable_genes <- sample_druggable_genes[Genes %in% ct.targets.tmp]
  
  if(nrow(ct_druggable_genes)){
    pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "/", sample.tmp, "_CT_druggable_targets_scDNA_scRNA_dotplot_with_legend.pdf"), width = length(unique(ct_druggable_genes$Genes))/10*6, height = 8)
    
    ct_druggable_genes[,Pct := mean(Expression > 0), .(Genes, Merged_Clone)]
    ct_druggable_genes[,Ave := mean(Expression), .(Genes, Merged_Clone)]
    
    ct_druggable_genes <- unique(ct_druggable_genes[,.(Genes, Merged_Clone, Pct, Ave, chr, start, end, Clone)])
    
    
    ct_druggable_genes$TCN <- apply(ct_druggable_genes, 1, function(x){
      cur.ct.table <- sample.cnv.table.m[Clone==x[["Clone"]]]
      myx <- which(cur.ct.table$start<as.numeric(x[["start"]]) & cur.ct.table$end>as.numeric(x[["start"]]) & cur.ct.table$chr == x[["chr"]]):
        which(cur.ct.table$start<as.numeric(x[["end"]]) & cur.ct.table$end>as.numeric(x[["end"]]) & cur.ct.table$chr == x[["chr"]])
      return(median(cur.ct.table[myx][["TCN"]]))
      
    })
    
    ct_druggable_genes[,Mean_TCN := mean(TCN), .(Merged_Clone, Genes)]
    
    ct_druggable_genes <- unique(ct_druggable_genes[,.(Genes, Merged_Clone, Pct, Ave, chr, start, end, Mean_TCN)])
    
    p <- ggplot(ct_druggable_genes, aes(x=Genes, y=Merged_Clone)) + 
      geom_tile(aes(fill=Mean_TCN)) + 
      scale_fill_gradient2(guide = guide_colorbar(ticks.colour = "black",
                                                  frame.colour = "black"),
                           name = "Average Copy Number", high="#b2182b", low="#4393c3", mid="#F6F6F6", midpoint = 2, limits=c(0,8)) +
      new_scale_fill() + 
      geom_point(aes(size=Pct, fill=Ave), color="black", shape=21) +
      scale_size("% detected", range = c(0,6), limits = c(0,1)) +
      scale_fill_gradientn(colours = viridisLite::mako(100),
                           guide = guide_colorbar(ticks.colour = "black",
                                                  frame.colour = "black"),
                           name = "Average\nexpression", limits = c(0,4)) +
      ylab("Cluster") + xlab("") +
      theme_classic() +
      theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
            axis.text.y = element_text(size=12, color="black"),
            axis.title = element_text(size=14)) + coord_equal()
    print(p)
    
    dev.off()
    p <- p + theme(legend.position = "none")
    
    pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "/", sample.tmp, "_CT_druggable_targets_scDNA_scRNA_dotplot.pdf"), width = length(unique(ct_druggable_genes$Genes))/10*6, height = 8)
    print(p)
    dev.off()
    
    png(paste0("infercnv/Figures_Paper/", sample.tmp, "/", sample.tmp, "_CT_druggable_targets_scDNA_scRNA_dotplot.png"), 
        width = length(unique(ct_druggable_genes$Genes))/10*6, height = 8, units = "in", res=600)
    print(p)
    dev.off()
    
  }
    
 
  
  
}


############################################################################
##            VISUALISE INDIVIDUAL GENES FOR STP-PDX
############################################################################

sample.tmp <- "STP-PDX"
gene <- "MYCN"
gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)

# SUBCLONAL CNV TARGETS
pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_MYCN_druggable_targets_scDNA_scRNA_inferCNV.pdf"), width = 6, height = 4)
plot_grid(scDNA.gene.plot, scRNA.gene.plot, nrow = 2, align = "v")
dev.off()

gene <- "CDK9"
gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)

# SUBCLONAL CNV TARGETS
pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_CDK9_druggable_targets_scDNA_scRNA_inferCNV.pdf"), width = 6, height = 4)
plot_grid(scDNA.gene.plot, scRNA.gene.plot, nrow = 2, align = "v")
dev.off()


gene <- "EZH2"
gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)

# SUBCLONAL CNV TARGETS
pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_EZH2_druggable_targets_scDNA_scRNA_inferCNV.pdf"), width = 6, height = 4)
plot_grid(scDNA.gene.plot, scRNA.gene.plot, nrow = 2, align = "v")
dev.off()

gene <- "HDAC1"
gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)

# SUBCLONAL CNV TARGETS
pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_HDAC1_druggable_targets_scDNA_scRNA_inferCNV.pdf"), width = 6, height = 4)
plot_grid(scDNA.gene.plot, scRNA.gene.plot, nrow = 2, align = "v")
dev.off()

gene <- "HDAC3"
gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)

# SUBCLONAL CNV TARGETS
pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_HDAC3_druggable_targets_scDNA_scRNA_inferCNV.pdf"), width = 6, height = 4)
plot_grid(scDNA.gene.plot, scRNA.gene.plot, nrow = 2, align = "v")
dev.off()


gene <- "DLL3"
gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)

# SUBCLONAL CNV TARGETS
pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_DLL3_druggable_targets_scDNA_scRNA_inferCNV.pdf"), width = 6, height = 4)
plot_grid(scDNA.gene.plot, scRNA.gene.plot, nrow = 2, align = "v")
dev.off()

############################################################################
##            VISUALISE INDIVIDUAL GENES FOR MB243-NUCLEI
############################################################################

sample.tmp <- "MB243-Nuclei"
# 
# gene <- "CDK7"
# gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
# scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
# scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)

# # SUBCLONAL CNV TARGETS
# pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_CDK9_druggable_targets_scDNA_scRNA.pdf"), width = 6, height = 4)
# plot_grid(scDNA.gene.plot, scRNA.gene.plot, nrow = 2, align = "v")
# dev.off()

gene <- "RICTOR"
gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)

# SUBCLONAL CNV TARGETS
pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_HDAC1_druggable_targets_scDNA_scRNA.pdf"), width = 6, height = 4)
plot_grid(scDNA.gene.plot, scRNA.gene.plot, nrow = 2, align = "v")
dev.off()

gene <- "HDAC3"
gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)

# SUBCLONAL CNV TARGETS
pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_HDAC3_druggable_targets_scDNA_scRNA.pdf"), width = 6, height = 4)
plot_grid(scDNA.gene.plot, scRNA.gene.plot, nrow = 2, align = "v")
dev.off()

gene <- "HDAC7"
gene_annot <- drug.target.annotations[drug.target.annotations$gene_name==gene,]
scDNA.gene.plot <- plot_genes(sample.tmp, gene, gene_annot = gene_annot)
scRNA.gene.plot <- plot_gene_expression(sample.tmp, gene)

# SUBCLONAL CNV TARGETS
pdf(paste0("infercnv/Figures_Paper/", sample.tmp, "_HDAC7_druggable_targets_scDNA_scRNA.pdf"), width = 6, height = 4)
plot_grid(scDNA.gene.plot, scRNA.gene.plot, nrow = 2, align = "v")
dev.off()



