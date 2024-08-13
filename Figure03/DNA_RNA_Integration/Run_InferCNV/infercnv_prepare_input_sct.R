################################################################################################################################################
##                                                                                                                      
##  Prepare input for inferCNV                                                                                            
##                                                                                                                      
##  Date: 11 March 2020                                                                                                                    
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
list.of.packages <- c("Matrix", "tidyverse", "stringr", "biomaRt")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

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
tumor.dir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/scRNA_analysis"

# list files
tumor.dir <- list.files(tumor.dir, pattern = "_CellType_seurat_obj.rds", full.names = T, recursive = T)

# get output directory
o.dir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/freeze"

# create output directory
dir.create(o.dir)

# get sample names
sample.ids <- str_split_fixed(basename(tumor.dir), "_", 2)[,1]

# make list of working directories
w.dir_list <- file.path(o.dir, paste0(sample.ids, "_microglia"))

# w.dir <- w.dir_list[7]
for(w.dir in w.dir_list){
  
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
  
  # read in the output matrices from Gonzalo for with REF
  seurat.obj <- readRDS(grep(str_split_fixed(basename(w.dir), "_microglia",2)[,1], tumor.dir, value = T))
  
  # get matrix
  matrix.tmp <- seurat.obj@assays$RNA@counts
  colnames(matrix.tmp) <- str_split_fixed(colnames(matrix.tmp), "-1", 2)[,1]
  
  # read in the barcodes for the normal cells if using a reference
  cell.clusters <- data.frame("Cell_Barcode" = str_split_fixed(seurat.obj$Cell_Barcode, "-1", 2)[,1], "cluster" = seurat.obj@meta.data[,"adapted_cluster"], "CellType" = seurat.obj$CellType)
  
  # get microglia cells
  cell.clusters_WT <- as.character(cell.clusters[grep("Microglia|Immune|Purkinje|Endothelial", cell.clusters$CellType), "Cell_Barcode"])
  
  # make annnotations 
  new.annotation <- cell.clusters[cell.clusters$Cell_Barcode %!in% cell.clusters_WT,]
  new.annotation$Cell_Barcode <- paste0(new.annotation$Cell_Barcode, "_pos")
  
  # get microglia cells
  cell.clusters_glia <- as.character(cell.clusters[grep("(Glia Cells)", cell.clusters$CellType), "Cell_Barcode"])
  
  # get normal cells
  normal.cells <- cell.clusters_WT
  tumor.cells <-  colnames(matrix.tmp)[colnames(matrix.tmp) %!in% cell.clusters_WT & colnames(matrix.tmp) %!in% cell.clusters_glia]
  mat.pos.list <- list(matrix.tmp[,tumor.cells])
  mat.neg.list <- list(matrix.tmp[,normal.cells])
  
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
  
  # check on the dimensions of the matrix
  # they will be displayed as a whole string
  print(dim(mat.pos))
  print(dim(mat.neg))
  print(dim(mat.matrix.sparse))
  
  # make a cell annotation table
  cellAnnotations <- data.frame(Cell_Barcode = colnames(mat.matrix.sparse))
  colnames(new.annotation) <- c("Cell_Barcode", "cluster_id")
  cellAnnotations <- merge(cellAnnotations, new.annotation, by = "Cell_Barcode", all.x = T)
  
  cellAnnotations$cluster_id <- paste0("malignant_", cellAnnotations$cluster_id)
  cellAnnotations$cluster_id[which(cellAnnotations$cluster_id == "malignant_NA")] <- "normal"
  
  # cellAnnotations$cluster_id[which(cellAnnotations$cluster_id != "normal")] <- "malignant"
  colnames(cellAnnotations)[3] <- "CellType"
  
  # cellAnnotations[grep("_pos", cellAnnotations$Cell_Barcode),2] <- "malignant"
  cellAnnotations <- distinct(cellAnnotations)
  cellAnnotations$CellType <- NULL
  
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
  
}
