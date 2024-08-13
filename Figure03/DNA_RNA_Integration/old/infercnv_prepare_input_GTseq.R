################################################################################################################################################
##                                                                                                                      
##  PREPARE INFERCNV INPUT ON THE STP-PDX DATA
##                                                                                                                      
##  Date: 03 NOVEMBER 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
##  Summary: This script shall produce the inferCNV inputs to use the inferCNV object afterwards.   
##           
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("Matrix", "tidyverse", "stringr", "biomaRt", "Seurat")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

library(httr)
set_config(config(ssl_verifypeer = 0L))

############################################################################
##                              Define Functions
############################################################################

# get count matrix with row and colnames
getCountMat <- function(cellranger_outs_folder){
  matrix_dir = file.path(cellranger_outs_folder,"filtered_feature_bc_matrix")
  barcode.path <- file.path(matrix_dir, "barcodes.tsv.gz")
  features.path <- file.path(matrix_dir, "features.tsv.gz")
  matrix.path <- file.path(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  return( mat )
}

# get the barcodes for the matrices
getBarcodeNames <- function(cellranger_outs_folder){
  matrix_dir = file.path(cellranger_outs_folder,"filtered_feature_bc_matrix")
  barcode.path <- file.path(matrix_dir, "barcodes.tsv.gz")
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  return( as.character(barcode.names[,1]) )
}

# merge the matrices together
check_and_merge <- function(sparse_matrix_list){
  if(length(sparse_matrix_list) < 2){
    return(sparse_matrix_list[[1]])
  }
  combtab <- combn(x = seq(length(sparse_matrix_list)),m = 2)
  
  genes_in_common <- list()
  for(i in seq(ncol(combtab))){
    genes_in_common[[i]] <- intersect(rownames(sparse_matrix_list[[combtab[1,i]]]), rownames(sparse_matrix_list[[combtab[2,i]]]))
  }
  
  cg <- genes_in_common[[1]]
  if(length(genes_in_common) > 1 ){
    for(i in 2:length(genes_in_common)){
      cg <- intersect(cg,genes_in_common[[i]])
    }
  }
  
  filter_genes <- function(x, genes){
    x <- x[genes,]
    return(x)
  }
  
  sparse_matrix_list <- lapply(sparse_matrix_list, filter_genes, cg)
  
  check.cols <- c()
  for(i in seq(ncol(combtab))){
    check.cols <- c(check.cols, identical(rownames(sparse_matrix_list[[combtab[1,i]]]), rownames(sparse_matrix_list[[combtab[2,i]]])) )
  }
  if(all(check.cols)){
    outmat <- do.call(cbind, sparse_matrix_list)
    return(outmat)
  } else {
    return(NULL)
  }
}

# negat &in&
`%!in%` = Negate(`%in%`)

is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}
############################################################################
##                              Read in data
############################################################################

# get input directory for tumor samples
input.dir <- "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis"

# get output directory
o.dir <- "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_MB"

# create output directory
dir.create(o.dir)

sample.id <- "GTseq_384_mouse_new"

# make list of working directories
w.dir <- file.path(o.dir, paste0(sample.id))

# create sample working/output directory
dir.create(w.dir)
message(w.dir)

# setwd 
setwd(w.dir)

# output cell annotation file
cell.annotation.name <- paste0("cellAnnotations_", basename(w.dir),".txt")

# output count file
counts.sparse.matrix.name <- paste0("sc.10x.counts_", basename(w.dir),".RData")

############################################################################
##   Get count matrices from which tumor barcodes will be extracted
############################################################################

# read in genes x cells matrix from 384 GT
mat.pos <- read.table('/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_384/gtseq_rna_counts.txt',header = T, stringsAsFactors = F)
dim(mat.pos)

# read in the fastq screen results
info.table <- read.table("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/20220201_fastq_screen_results.txt", header = T, sep = "\t")
info.table

ggplot(info.table, aes(x=DNA, y = RNA)) +
  geom_point() +
  theme_classic() +
  facet_grid( ~ Genome, scales = "free") + 
  geom_smooth(method=lm, color="darkblue") +
  stat_cor(method = "pearson") +
  labs(x="scDNA [%]", y="scRNA [%]") + 
  theme(legend.position="bottom") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        strip.background = element_rect(fill="white", colour="black", size=1),
        strip.text = element_text(face="bold", size=12, colour = "black",),
        axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"), 
        legend.position = "top")

# gate to the cells that have high mapping in RNA and DNA for mouse and human
info.table <- info.table[info.table$DNA >= 75 & info.table$DNA >= 75, ]

# split in mouse and human
human.table <- info.table[info.table$Genome == "Human", ]
mouse.table <- info.table[info.table$Genome == "Mouse", ]

# get the human cells from ht.mat.pos
ht.mat.pos <- mat.pos[, colnames(mat.pos) %in% human.table$Sample.y]
dim(ht.mat.pos)

# read in genes x cells matrix from 96 GT
seurat.obj <- readRDS("/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/analysis_Anna/scRNA/Seurat_objects/Seurat_obj_v3_71cells_cc-reg.rds")
lt.mat.pos <- seurat.obj@assays$RNA@counts
dim(lt.mat.pos)

# read in the reference PDX cells from Hovestadt et al.
# mat.neg <- read.table('/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/GSE119926_RAW/GSM3905435_Med2112FH.txt.gz',header = T, stringsAsFactors = F)
# dim(mat.neg)

# get the mouse cells from 384 GT
mat.neg <- mat.pos[, colnames(mat.pos) %in% mouse.table$Sample.y]
dim(mat.neg)

mat.pos.list <- list(ht.mat.pos, lt.mat.pos)
mat.pos.list <- list(ht.mat.pos)
mat.neg.list <- list(mat.neg)

# add label for cells from the same sample
for(i in seq(length(mat.pos.list))){
  colnames(mat.pos.list[[i]]) <- paste0(colnames(mat.pos.list[[i]]), "_pos")
}

for(i in seq(length(mat.neg.list))){
  colnames(mat.neg.list[[i]]) <- paste0(colnames(mat.neg.list[[i]]), "_neg")
}

############################################################################
##
############################################################################

# merge the matrix lists together
mat.pos <- check_and_merge(sparse_matrix_list = mat.pos.list)

# reference data
mat.neg <- check_and_merge(sparse_matrix_list = mat.neg.list)

# merge matrices
mat.matrix.sparse <- check_and_merge(sparse_matrix_list = list(mat.pos, mat.neg))
dim(mat.matrix.sparse)

# make a cell annotation table
cellAnnotations <- data.frame(Cell_Barcode = colnames(mat.matrix.sparse))
cellAnnotations$cluster_id <- NA

# add the malignant label
cellAnnotations[grep("_pos", cellAnnotations$Cell_Barcode), "cluster_id"] <- "malignant"
cellAnnotations$cluster_id[which(is.na(cellAnnotations$cluster_id))] <- "normal"
cellAnnotations <- distinct(cellAnnotations)

# write df
write.table(cellAnnotations, file= cell.annotation.name, col.names = F, row.names = F,quote=F, sep="\t")

# remove non-unique barcodes
mat.matrix.sparse <- mat.matrix.sparse[,unique(colnames(mat.matrix.sparse))]

# merge pos and neg
save(mat.matrix.sparse, file = counts.sparse.matrix.name, compress = T)

# make a gene ordering file
if(TRUE){
  
  # get mart object with certain attributes
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  ann <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"), "hgnc_symbol", rownames(mat.matrix.sparse), mart  = mart)
  
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
  gene_ordering_file <- gene_ordering_file[!duplicated(gene_ordering_file$hgnc_symbol),]
  
  # write gene ordering file
  write.table(gene_ordering_file, file= "gene_ordering_file.txt", col.names = F, row.names = F, quote=F, sep="\t")
}

