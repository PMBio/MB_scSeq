#############################################################################################################################
##                                                                                                                      
##  CREATE INPUT FILES FOR DRUGGABLE TARGET VISUALISATIONS
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
list.of.packages <- c("reshape2", "optparse", "RColorBrewer", "ggplot2", "scales", "dendextend", "tidyverse", 
                      "Matrix", "devtools", "readr","BiocManager", "biomaRt", "httr", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))


remove_uncertain_mapping = TRUE

############################################################################
##                          READ IN THE DATA
############################################################################

# negate %in%
"%ni%" <- Negate("%in%")

# specify the output directory
out.dir <- "infercnv/revision/scrna_analysis/infercnv_MB/scRNA_scDNA"

# read in the expression files from scanpy
mtx.files <- list.files("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/scRNA_analysis/scanpy", pattern = "raw_matrix", full.names = T)

# ## GET ALL THE SAMPLE INFORMATION REQUIRED FOR 10X RNA DATA FROM INFERCNV
# # list all scRNA matrices
# mtx.files <- list.files("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB", pattern = "infercnv.median_filtered.observations.txt", full.names = T, recursive = T)
# mtx.files <- mtx.files[grep("freeze", mtx.files)]
# mtx.files <- mtx.files[grep("s_tumor/|s_mouse", mtx.files)]
# mtx.files <- mtx.files[-c(1, 2, 5)]

# list the metadata files
metadata.files <- list.files('/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_MB/', pattern = "metadata_scAbsolute.csv", full.names = T, recursive = T, all.files = T)

# list the metadata files
clone.files <- list.files('/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_MB/scRNA_scDNA', pattern = "filtered", full.names = T, recursive = T, all.files = T)
clone.files <- clone.files[grep("_clones/", clone.files)]
clone.files <- clone.files[grep("DEG", clone.files)]

# list the raw expression matrix from 10x
# raw.mtx.files <- list.files("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/data/10XRNA5P/", pattern = "filtered_feature_bc_matrix", recursive = T, full.names = T, include.dirs = T)
# raw.mtx.files <- raw.mtx.files[-grep("\\.h5|mm10|old", raw.mtx.files)]

# READ IN DRUGGABLE TARGET LIST
druggable.targets <- read.table("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/scanpy/druggable_targets_list.txt", header = F)

# list to store
druggable.sample.list <- list()

# get gene ordering files
gene.order.files <- list.files("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/freeze", pattern = "gene_ordering_file.txt", recursive = T, full.names = T)
gene.order.files <- gene.order.files[grep("tumor", gene.order.files)]

# read in the hg19 chromosome cooordinates
x <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/cytoBand.txt.gz", col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
chr.arm.pos <- x[ , .(length = sum(chromEnd - chromStart)), by = .(chrom, arm = substring(name, 1, 1)) ]

# create start and end coordinates
chr.arm.pos$start <- 0
chr.arm.pos[chr.arm.pos$arm == "q", "start"] <- chr.arm.pos[chr.arm.pos$arm == "q", "start"] +chr.arm.pos[chr.arm.pos$arm == "p", "length"]
chr.arm.pos$end <- chr.arm.pos$length
chr.arm.pos[chr.arm.pos$arm == "q", "end"] <- chr.arm.pos[chr.arm.pos$arm == "q", "end"] +chr.arm.pos[chr.arm.pos$arm == "p", "end"]

# make a Grange object
chr.arm.GRange <- makeGRangesFromDataFrame(chr.arm.pos, keep.extra.columns = T)

i <- 1
for (i in 1:length(metadata.files)){
  
  # read in the matrix
  sample.tmp <- str_split_fixed(basename(metadata.files[i]), "_", 2)[,1]
  
  # # read in the matrix of interest
  sample.mtx <- fread(mtx.files[grep(sample.tmp, mtx.files)], header = T, sep = ",", data.table=FALSE)
  rownames(sample.mtx) <- sample.mtx$V1
  sample.mtx$X <- NULL
  
  # read in the matrix of interest
  # sample.mtx <- as.matrix(read.table(mtx.files[grep(sample.tmp, mtx.files)], header = T, sep = " "))
  # colnames(sample.mtx) <- paste0(str_split_fixed(colnames(sample.mtx), "_pos", 2)[,1], "-1")
  # sample.mtx <- as.data.frame(t(sample.mtx))
  
  # read in gene_ordering file
  gene.order.file <- gene.order.files[grep(sample.tmp, gene.order.files)]
  gene.ordering.file <- read.table(gene.order.file, header = F, col.names = c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id"),sep = "\t")
  gene.ordering.file$coordinates <- paste(gene.ordering.file$chr, gene.ordering.file$start, gene.ordering.file$end, sep = ":")
  gene.GRange <- makeGRangesFromDataFrame(gene.ordering.file, keep.extra.columns = T)
  
  # then merge each dataframe information by overlap
  # find the overlapping regions of the WGS segments with the genes
  gene.GR.merge <- mergeByOverlaps(chr.arm.GRange, gene.GRange)
  
  # get the essential columns and make a new dataframe
  gene.GR.df <- data.frame("coordinates" = gene.GR.merge$chr.arm.GRange,
                           "chr_arm" = gene.GR.merge$arm,
                           "hgnc_symbol" = gene.GR.merge$hgnc_symbol,
                           "coordinates" = gene.GR.merge$coordinates,
                           stringsAsFactors = F)
  
  # convert each row
  gene.GR.df$chr_arm <- paste0(gene.GR.df$coordinates.seqnames, gene.GR.df$chr_arm)
  
  # subset to chromosomes were there is a change expected
  if (sample.tmp == "STP-PDX"){
    
    # 
    variable.chroms <- c("chr2p", "chr2q", "chr15p", "chr15q", "chr19q", "chr22p", "chr22q")
    
  } else if (sample.tmp == "MB243-Nuclei"){
    
    variable.chroms <- c("chr1p", "chr4q", "chr5p", "chr7p", "chr19p", "chr19q")
    
  } else if (sample.tmp == "STP-Nuclei"){
    
    variable.chroms <- c("chr3q", "chr14q")
    
  } else if (sample.tmp == "ST1R-PDX"){
    
    variable.chroms <- c("chr8q", "chr12p", "chr12q")
    
  }
  
  # subset genes to the chromosomes of interest
  gene.GR.df <- gene.GR.df[gene.GR.df$chr_arm %in% variable.chroms,]
  genes <- as.character(gene.GR.df$hgnc_symbol)
  genes <- genes[order(genes)]
  
  # read in the metadata
  sample.metadata <- fread(grep(sample.tmp, metadata.files, value = T), header = T)
  sample.metadata <- sample.metadata[,c("V1", "clone_id", "new_clusters")]
  colnames(sample.metadata) <- c("Cell_barcodes", "clone_id", "celltype")
  
  # read in the clone data
  clone.metadata <- fread(grep(sample.tmp, clone.files, value = T), header = T)
  
  # Optionally, remove all cells that are still multi-mapping between clones
  
  if(remove_uncertain_mapping){
    clone.metadata$cell_3rd_best_match <- NA_real_
    clone.metadata$clone_3rd_best_match <- NA_character_
    clone.metadata$distance_3rd_best_match <- NA_real_
    for (ii in 1:nrow(clone.metadata)){
      
      cell.data.tmp <- clone.metadata[ii,]
      cell_best_matched <- cell.data.tmp$pearson.correlation
      cell_3rd_best_match <- suppressMessages(reshape2::melt(cell.data.tmp[,grepl("Cell_barcode|Clone[1-5]$", colnames(cell.data.tmp)), with=FALSE]))
      # cell_3rd_best_match <- reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3", "Clone4")])
      # cell_3rd_best_match <- reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3")])
      clone_3rd_best_match <- as.character(cell_3rd_best_match[order(cell_3rd_best_match$value, decreasing = T), "variable"])[3]
      cell_3rd_best_match_corr <- as.numeric(cell_3rd_best_match[order(cell_3rd_best_match$value, decreasing = T), "value"])[3]
      distance_3rd_best_match <- as.numeric(cell_3rd_best_match[order(cell_3rd_best_match$value, decreasing = T), "value"])[1] - as.numeric(cell_3rd_best_match[order(cell_3rd_best_match$value, decreasing = T), "value"])[3]
      
      clone.metadata[ii, cell_3rd_best_match:=cell_3rd_best_match_corr]
      clone.metadata[ii, clone_3rd_best_match := ..clone_3rd_best_match]
      clone.metadata[ii, distance_3rd_best_match := ..distance_3rd_best_match]
      
    }
    
    merged_clone_tables <- table(clone.metadata$merged_clone_id, clone.metadata$clone_id)
    
    merged_clone_names <- names(which(rowSums(merged_clone_tables>0)>1))
    clone.metadata[,dist_filter := distance]
    for(merged_clone in merged_clone_names){
      collapsing_clones <- unique(clone.metadata[clone.metadata$merged_clone_id==merged_clone,clone_id])
      clone.metadata[clone_id %in% collapsing_clones & clone_2nd_best_match %in% collapsing_clones, dist_filter := distance_3rd_best_match]
    }
    clone.metadata[,keep.cell := dist_filter >= 0.025]
    
    if(length(unique(clone.metadata[keep.cell==TRUE, merged_clone_id])) != length(unique(clone.metadata[, merged_clone_id]))){
      stop("Sample lost a full clone after filtering!")
    }
    clone.metadata <- clone.metadata[keep.cell == TRUE]
  }
  
  clone.metadata <- clone.metadata[,c("Cell_barcode", "clone_id", "merged_clone_id")]
  colnames(clone.metadata) <- c("Cell_barcodes", "clone_id", "merged_clone_id")
  
  sample.metadata <- as.data.frame(merge(sample.metadata, clone.metadata, by = c("Cell_barcodes", "clone_id")))
  sample.metadata <- sample.metadata[complete.cases(sample.metadata),]
  
  # subset the matrix
  druggable.mtx <- sample.mtx[, colnames(sample.mtx) %in% druggable.targets$V1]
  # druggable.mtx <- sample.mtx[, colnames(sample.mtx) %in% genes]
  # druggable.mtx <- sample.mtx
  druggable.mtx$Cell_barcodes <- rownames(druggable.mtx)
  
  # melt it 
  druggable.mtx.melt <- reshape2::melt(druggable.mtx)
  
  # merge with clone id
  druggable.df <- merge(druggable.mtx.melt, sample.metadata, by = "Cell_barcodes")
  druggable.df$sample <- sample.tmp
  druggable.df$clone_id <- factor(druggable.df$clone_id, levels = c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5"))
  druggable.df$merged_clone_id <- factor(druggable.df$merged_clone_id, levels = c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5"))
  
  # exclude NAs
  druggable.df <- druggable.df[!is.na(druggable.df$clone_id),]

  # check which genes are expressed in more than 5% of cells
  threshold <- round(nrow(sample.mtx)*0.05)
  druggable.non.zero.df <- druggable.df[druggable.df$value > 0,]
  genes <- names(table(druggable.non.zero.df$variable))[table(druggable.non.zero.df$variable) > threshold]
  
  druggable.df <- druggable.df[druggable.df$variable %in% genes,]
  
  sd.per.gene <- druggable.df %>%
    group_by(variable) %>%
    summarise_at(vars(value), list(name=sd))
  sd.per.gene <- sd.per.gene[order(sd.per.gene$name, decreasing = T),] 
  genes <- as.character(sd.per.gene[1:200, "variable"]$variable)
  
  # druggable.df <- druggable.df[druggable.df$variable %in% genes,]
  
  # add to list
  druggable.sample.list[[i]] <- druggable.df
  
}

complete.druggable.targets.df <- bind_rows(druggable.sample.list)
colnames(complete.druggable.targets.df) <- c("Cell_barcodes", "Genes", "Expression", "Clone", "Celltype", "Merged_Clone", "Sample")
write.table(complete.druggable.targets.df, paste0(out.dir, "/druggableGenes_allSamples_expression_w_zeros.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
# write.table(complete.druggable.targets.df, paste0(out.dir, "/allGenes_allSamples_expression_w_zeros.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
# write.table(complete.druggable.targets.df, paste0(out.dir, "/druggableGenes_allSamples_inferCNV.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
# write.table(complete.druggable.targets.df, paste0(out.dir, "/selectedGenes_allSamples_expression_w_zeros.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
# write.table(complete.druggable.targets.df, paste0(out.dir, "/selectedGenes_allSamples_raw.txt"), sep = "\t", quote = F, col.names = T, row.names = F)


