#############################################################################################################################
##                                                                                                                      
##  Make list of cells with amplification of chr8q in C1 R2 T7 and T12                                                                                       
##                                                                                                                      
##  Date: 20 April 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
############################################################################################################################

# clear workspace
rm(list=ls())
set.seed(14) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if nnot available
list.of.packages <- c("tidyverse", "Seurat", "biomaRt", "RColorBrewer", "hdf5r", "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg38",
                      "viridis", "fields")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

###########################################################################
#                             READ IN DATA
###########################################################################

# set the input to the DNA matrices
# w.dir <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/phylogenetics"
w.dir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/scDNA_gene_matrices"

# define the output directory for the heatmaps
o.dir <- "/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/scCNV_heatmaps/"

# list all the cnv.cells.mat.txt files
scDNA.matrices <- list.files(w.dir, pattern = "_scDNA_gene_matrix_new.txt", recursive = T, full.names = T)

# list the file with clusters
community.files <- list.files("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/DNA_processed_data_cellranger-dna_v.1.1.0/", pattern = "cell_ids_clusters", recursive = T, full.names =  T)
community.files <- community.files[grep("Final_Clustering_Trees", community.files)]
community.files <- community.files[grep("cell_ids_clusters_normalized", community.files)]

# get the sample ids
sample.ids <- str_split_fixed(basename(scDNA.matrices), "_", 2)[,1]

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
rownames(gene_ordering_file) <- gene_ordering_file$hgnc_symbol 

i <- 1
for (i in 1:length(scDNA.matrices)){
  
  # get sample name
  sample.tmp <- sample.ids[i]
  message(sample.tmp)

  ###########################################################################
  #                   READ IN DATA OF INTEREST FROM inferCNV
  ###########################################################################
  
  # read in ploidy matrix from scDNA
  scDNA.matrix <- read.table(scDNA.matrices[i], header = T, stringsAsFactors = F)
  colnames(scDNA.matrix) <- gsub("X", "", colnames(scDNA.matrix))
  
  # get dimensions
  print(dim(scDNA.matrix))
  
  # read in community files
  community.file <- read.table(community.files[i], header = T)
  community.file$cluster <- factor(community.file$cluster, levels = c(1:length(unique(community.file$cluster))))
  community.file <- community.file[order(community.file$cluster),]
  
  #####################################################################################
  # GENERATE A HEATMAP FOR THE LOG2FC
  #####################################################################################
  
  # only take genes which are present in the matrix
  # make it a GRange object
  genes <- makeGRangesFromDataFrame(gene_ordering_file, keep.extra.columns = T)
  genes <- genes[rownames(scDNA.matrix)]
  genes <- genes[order(genes@seqnames, genes@ranges)]
  
  # generate new object genes.interest with genes and respective coordinates
  genes.interest <- as.data.frame(genes)
  rownames(genes.interest) <- names(genes)
  
  # only take features which are present in the matrix
  # windowSummary is an grange object with windows names as rownames, the chromsome as seqnames and ranges
  # convert back to dataframe again with genes as rownames seqnames start end width
  mat <- as.matrix(scDNA.matrix)

  # check for each chromosome if there are genes present
  tl <- tapply(1:nrow(genes.interest),as.factor(genes.interest$seqnames),function(ii) {
    na.omit(mat[rownames(genes.interest)[ii[order((genes.interest[ii,]$start+genes.interest[ii,]$end)/2,decreasing=F)]],,drop=FALSE])
  })
  
  # subset to chromosome 3,4,9,13 and 20 (chr of interest in general)
  chrs <- paste0("chr", 1:22)
  tl <- tl[chrs]
  
  # we want to scale the output to the respective amount of genes and their widths
  # so here, the widths for each window are determined
  setWidths <- T
  chromSizes <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1:22]
  if(setWidths) {
    # widths <- sapply(tl, nrow); widths <- widths/max(widths)*100
    # Can also set by known chromosome size widths:
    widths <- end(chromSizes)[1:22]/1e7
  } else {
    widths <- rep(1, length(tl))
  }
  
  # generate the layout accordingly
  l <- layout(matrix(seq(1,length(tl)),1,length(tl),byrow=TRUE), widths=widths)
  
  # set parameters for the image
  tlsub <- tl
  # pcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(256) 
  pcol <- c("#053061", "#3C8ABE", "#ffffff", "#F4AA88", "#EB9273", "#D6604D", "#C7433F", "#B51F2E", "#67001F")
  zlim <- c(0, 8) # upper and lower limit
  
  # limit the gExp to maximums according to zlim
  tlsmooth <- lapply(names(tlsub),function(nam) {
    d <- tlsub[[nam]]
    d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
    return(d)
  })
  
  #####################################################################################
  #                   GENERATE A COLOUR CODE ACCORDING TO THE COMMUNITITES
  #####################################################################################
  
  # make an dataframe with an index for each ECB_RG
  order.file <- community.file
  order.file$index <- c(1:nrow(order.file))
  # vector <- order.file[order.file$cluster == 4,"cell_id"]
  
  # make a df with the index to color it accordingly
  community.col.bar <- data.frame(as.numeric(order.file$cluster),as.numeric(order.file$cluster))
  colnames(community.col.bar)<-c("rg1","rg2")
  # community.col.bar <- community.col.bar[community.col.bar$rg1 == 4,]
  
  # get the coloring according to RGs
  # cols <- c(as.character(unique(order.file$cluster)))
  cols <- c(brewer.pal(n=length(unique(community.file$cluster)),"Set1"))
  
  #####################################################################################
  #                                 Make the plot
  #####################################################################################

  png(paste0(o.dir, sample.tmp, "_10x_scCNV_gene_heatmap.png") , width = 26, height = 12, units = 'in', res = 600)
  par(mfrow=c(1,23), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
  image(t(community.col.bar),xlab="",ylab="", axes=FALSE, col=cols)
  ## plot chromosomes
  box()
  for (i in 1:length(tlsmooth)){
    message(paste0("chr", i))
    d <- tlsmooth[[i]]
    d <- d[, as.character(order.file$cell_id)]
    # image the respective chr 
    image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=F, main=paste0("chr", i), useRaster = T)
    box()
  }
  dev.off()
  
}

png(paste0(o.dir, "10x_scCNV_heatmap_legend.png") , width = 13, height = 4, units = 'in', res = 600, type = "cairo")
par(mfrow=c(1,23), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.5)
image.plot(legend.only=TRUE, zlim=zlim, col = pcol, legend.shrink = 0.4, legend.width = 10, smallplot=  c(.65, .9, .5, .8)) # first argument of small plot = width of legend, third is making the length, the higher the smaller
dev.off()
