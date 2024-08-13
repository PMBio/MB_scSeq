#############################################################################################################################
##                                                                                                                      
##  Make scDNA genes x cells matrix                                                                                          
##                                                                                                                      
##  Date: 9 June 2020                                                                                                                   
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
                      "Matrix", "devtools", "Matrix.utils", "matrixStats", "readr", "magrittr", "Signac", "BiocManager", 
                      "biomaRt", "httr", "ComplexHeatmap")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

# negate %in%
"%ni%" <- Negate("%in%")

#####################################################################################
# READ IN DATA OF INTEREST
#####################################################################################
setwd("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/scDNA_gene_matrices")

# read gene ordering file containing gene locations
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org")
ann <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"), mart  = mart)

# order cns_reads
gene_ordering_file <- ann
colnames(gene_ordering_file) <- c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id")
gene_ordering_file$chr <- paste0("chr", gene_ordering_file$chr)

# set accepted levels and remove other chromosomes
chrOrder<-c(paste("chr",1:22,sep=""))
gene_ordering_file$chr <- factor(gene_ordering_file$chr, levels = chrOrder)
gene_ordering_file <- gene_ordering_file[complete.cases(gene_ordering_file),]

# order file and remove duplicated values
gene_ordering_file <- gene_ordering_file[with(gene_ordering_file, order(chr, start)),]

# remove empty gene values
gene_ordering_file <- gene_ordering_file[grep(".", gene_ordering_file$hgnc_symbol),]
gene_ordering_file <- gene_ordering_file[!duplicated(gene_ordering_file$hgnc_symbol),]

# make it a GRange object
gene.pos.GRange <- makeGRangesFromDataFrame(gene_ordering_file, keep.extra.columns = T)

# list all scDNA matrices
scDNA.files <- list.files("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/DNA_processed_data_cellranger-dna_v.1.1.0", pattern = "cnv.cells.mat.20kb-bin.corrected.txt", recursive = T, full.names = T)
scDNA.files <- scDNA.files[-grep("old", scDNA.files)]
scDNA.files <- file.path("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/HMMcopy_GandT/data/matrix.ploidy.bin_20kb.txt")

# and the sample ids
sample.ids <- str_split_fixed(scDNA.files, "/", 11)[,9]
sample.ids <- "GT_seq"

#####################################################################################
# ITERATE OVER ALL SAMPLES AND CREATE THE RESPECTIVE GENES X CELLS MATRICES
#####################################################################################

i <- 1
for (i in 1:length(sample.ids)){
  
  # which sample? 
  sample.tmp <- sample.ids[i]
  message(sample.tmp)
  
  # read in ploidy matrix from scDNA (construct_ploidyMatrix_scDNA.R)
  scDNA.matrix <- as.matrix(read.table(scDNA.files[i], sep = "\t", header = T, stringsAsFactors = F))
  
  # get dimensions
  print(dim(scDNA.matrix))
  
  # get the segments and make GRange from it 
  seg.coordinates <- data.frame(chr = paste0("chr", str_split_fixed(rownames(scDNA.matrix), ":",2)[,1]), 
                                start = str_split_fixed(str_split_fixed(rownames(scDNA.matrix), ":",2)[,2], "-", 2)[,1], 
                                end = str_split_fixed(str_split_fixed(rownames(scDNA.matrix), ":",2)[,2], "-", 2)[,2])
  seg.coordinates.GRange <- makeGRangesFromDataFrame(seg.coordinates, keep.extra.columns = T)
  
  # remove id column
  rownames(scDNA.matrix) <- paste(seg.coordinates$chr, seg.coordinates$start, seg.coordinates$end, sep = "-")
  # scDNA.matrix <- scDNA.matrix[,-1]
  
  # melt df/matrix with logR values
  melted.scDNA.matrix <- reshape2::melt(scDNA.matrix)
  colnames(melted.scDNA.matrix) <- c("seg_coordinates", "cell_id", "copy_number")
  melted.scDNA.matrix <- melted.scDNA.matrix %>% separate(seg_coordinates, c("chr", "start", "end"), "-")
  
  # make GRange object for scDNA
  scDNA.genes <- makeGRangesFromDataFrame(melted.scDNA.matrix, keep.extra.columns = T)
  
  # then merge each dataframe information by overlap
  # find the overlapping regions of the WGS segments with the genes
  scDNA.GR.merge <- mergeByOverlaps(gene.pos.GRange, scDNA.genes)
  
  # get the essential columns and make a new dataframe
  chr.gene.scDNA <- data.frame("hgnc_symbol" = scDNA.GR.merge$hgnc_symbol,
                              "cell_id" = scDNA.GR.merge$cell_id,
                              "copy_number" = scDNA.GR.merge$copy_number,
                              stringsAsFactors = F)
  
  # convert each row
  chr.gene.scDNA$chr_arm <- as.character(chr.gene.scDNA$hgnc_symbol)
  chr.gene.scDNA$cell_id <- as.character(chr.gene.scDNA$cell_id)
  
  # remove the weird formatting to have actual numbers instead
  new.cn.col <- gsub("  ", "", as.character(chr.gene.scDNA$copy_number))
  new.cn.col <- gsub(" ", "", new.cn.col)
  new.cn.col <- as.numeric(new.cn.col)
  chr.gene.scDNA$copy_number <- as.numeric(new.cn.col)
  
  # reshape the melted dataframe into an matrix 
  scDNA.gene.matrix <- reshape2::acast(chr.gene.scDNA, hgnc_symbol ~ cell_id , value.var = "copy_number", fun.aggregate = mean)
  
  # check the dims
  print(dim(scDNA.gene.matrix))
  
  # write scDNA gene-level matrix
  write.table(scDNA.gene.matrix, paste0("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/scDNA_gene_matrices/", sample.tmp, "_scDNA_gene_matrix_new.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
  
}

