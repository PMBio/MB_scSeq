#############################################################################################################################
##                                                                                                                      
##  READ IN FUNCTIONS FOR THE SCDNA/RNA ALIGNMENT
##                                                                                                                      
##  Date: 18 OCTOBER 2021                                                                                                                   
##  
##  Author: Moritz Przybilla
##
##                                                                                                                      
############################################################################################################################


############################################################################
##                   MERGE MATRICES TOGETHER FROM A LIST
############################################################################

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


############################################################################
##              TRANSFORM SCDNA MATRIX TO CLONE PROFILES 
############################################################################

## INPUT
# scDNA.gene.mtx - FILE PATH TO THE GENE X CELLS MATRIX GENERATED FROM THE 10X CNV CALLS
# scDNA.cluster.file - FILE PATH TO THE CLUSTER FILE WHICH INCLUDES THREE COLUMNS: <cell_id> <cluster_id>
# translation.file - FILE PATH TO THE PER CELL SUMMARY METRICS FROM 10X CNV TO ASSOCIATE BARCODES WITH THE CELL IDS
# normalize - SMALL NORMALIZATION FOR WHOLE GENOME DOUBLED POPULATIONS TO MAKE IT COMPATIBLE WITH SCRNA-SEQ DATA
# ploidy - SPECIFY TO WHICH PLOIDY YOU WANT TO NORMALIZE (DEFAULT = DIPLOID, 2)
## OUTPUT
# t.clone.gene.profile - DATA FRAME WITH GENES X CLONES 

create.scDNA.clone.profiles <- function(scDNA.gene.mtx, scDNA.cluster, seq_type = "10X", translation, normalize = FALSE, ploidy = 2){
  
  # read in scDNA matrix
  matrix <- read.table(scDNA.file, sep = "\t")
  colnames(matrix) <- gsub("X", "", colnames(matrix))
  
  # perform normalization to a ploidy of interest
  if (normalize == TRUE){
    
    # define ploidy to which you want to normalize to
    ploidy <- ploidy
    
    # normalize for ploidy before calculating distance
    matrix <- apply(matrix, 2, function(x) round( x*(if(ploidy/median(x) > 1) floor(ploidy/median(x)) else ploidy/median(x))) ) 
    
  } else {
    next
  }
  
  # read in cluster file
  # cluster.file <- read.table(scDNA.cluster, header = T, sep = "\t")
  cluster.file <- read.table(scDNA.cluster.file, header = T, sep = ",")
  colnames(cluster.file) <- c("cell_id", "cluster_id")
  
  if (seq_type == "10X"){
   
    # read in the file containing the barcodes corresponding to the cell ids
    translation.file <- read.table(translation, header = T)
    translation.file <- translation.file[,c("cell_id", "barcode")]  
    
    # merge cluster and translation file
    # creates a file with <cell_id> <cluster_id> <barcode>
    cluster.file <- merge(cluster.file, translation.file, by = "cell_id")
    
  }
  
  # transpose the scDNA matrix to cells x genes
  t.matrix <- as.data.frame(t(matrix))
  t.matrix$cell_id <- rownames(t.matrix)
  
  # merge the transposed single cell DNA matrix with the cluster id file
  cluster.scDNA.mtx <- merge(t.matrix, cluster.file, by = "cell_id", all.x = T)
  
  if (seq_type == "10X"){
    
    # clean up the matrix
    cluster.scDNA.mtx <- cluster.scDNA.mtx[!is.na(cluster.scDNA.mtx$barcode),]
    rownames(cluster.scDNA.mtx) <- cluster.scDNA.mtx$barcode
    
    cluster.scDNA.mtx[,c("cell_id", "barcode")] <- NULL
    cluster.scDNA.mtx <- cluster.scDNA.mtx[complete.cases(cluster.scDNA.mtx),]
    
  } else {
    
    rownames(cluster.scDNA.mtx) <- cluster.scDNA.mtx$cell_id
    cluster.scDNA.mtx[,c("cell_id")] <- NULL
    cluster.scDNA.mtx <- cluster.scDNA.mtx[complete.cases(cluster.scDNA.mtx),]
    
  }
  
  # calculate the average gene profile per clone
  clone.gene.profile <- aggregate(cluster.scDNA.mtx, list(cluster.scDNA.mtx$cluster_id), mean)
  colnames(clone.gene.profile)[1] <- "clone_id"
  clone.gene.profile$clone_id <- paste0("Clone", clone.gene.profile$clone_id)
  
  # transpose clone gene profile
  t.clone.gene.profile <- as.data.frame(t(clone.gene.profile))
  colnames(t.clone.gene.profile) <- c(paste0("Clone", c(1:ncol(t.clone.gene.profile))))
  t.clone.gene.profile <- t.clone.gene.profile[-c(1,nrow(t.clone.gene.profile)),]
  
  # make numeric for normalization
  t.clone.gene.profile <- as.data.frame(lapply(t.clone.gene.profile[,c(1:ncol(t.clone.gene.profile))], function(x) as.numeric(as.character(x))))
  rownames(t.clone.gene.profile) <- colnames(clone.gene.profile)[2:(ncol(clone.gene.profile)-1)]
  
  # calculate the standard deviation
  std.dev <- sd(as.matrix(t.clone.gene.profile))
  
  # transform the modified expression values to integer values
  # this assumes that we can only reliably call losses or gains
  # inferCNV modified expression values are centered around 1
  t.clone.gene.profile[which(t.clone.gene.profile > 2 + std.dev, arr.ind = T)] <- 3
  t.clone.gene.profile[which(t.clone.gene.profile <= 2 + std.dev & t.clone.gene.profile >= 2 - std.dev, arr.ind = T)] <- 2
  t.clone.gene.profile[which(t.clone.gene.profile < 2 - std.dev, arr.ind = T)] <- 1
  
  # return the final vector with genes x clones
  return(t.clone.gene.profile)
  
}

############################################################################
##              VISUALIZE SCDNA GENE HEATMAP
############################################################################

# # VISUALIZE THE DNA GENE MATRIX
# scDNA.gene.matrix = scDNA.file
# order.file = scDNA.cluster.file 
# sample_name = sample.tmp
# output_directory = out.dir

plot_scDNA_gene_heatmap <- function(scDNA.gene.matrix, gene.ordering.file, order.file, sample_name, output_directory = "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA/"){
  
  # read in the file of interest  
  norm.gExp.matrix <- read.table(scDNA.gene.matrix, header = T)
  
  # define ploidy to which you want to normalize to
  ploidy <- 2
  gexp.norm <- norm.gExp.matrix
  norm.gexp.norm <- apply(gexp.norm, 2, function(x) round( x*(if(ploidy/median(x) > 1) floor(ploidy/median(x)) else ploidy/median(x))) ) 
  
  # order the cells according to the clustering file
  order.data <- read.table(order.file, header = F, sep = ",")
  colnames(order.data) <- c("cell", "cluster")
  order.data$cluster <- order.data$cluster + 1 
  gexp.norm <- norm.gexp.norm[,order.data$cell]
  
  # get the clusters for the respective sample stored as a dataframe with a clone_id column
  clusters <- unique(order.data$cluster)
  valid.clones <- factor(clusters, levels = c(1:as.numeric(max(unique(clusters)))))
  valid.clones <- valid.clones[order(valid.clones)]
  
  # only take genes which are present in the matrix
  genes <- makeGRangesFromDataFrame(gene.ordering.file)
  genes <- genes[names(genes) %in% rownames(gexp.norm),]
  
  # generate new object genes.interest with genes and respective coordinates
  genes.interest <- as.data.frame(genes)
  rownames(genes.interest) <- names(genes)
  mat <- gexp.norm
  
  # check for each chromosome if there are genes present
  gExp.array <- tapply(1:nrow(genes.interest),as.factor(genes.interest$seqnames),function(ii) {
    na.omit(mat[rownames(genes.interest)[ii[order((genes.interest[ii,]$start+genes.interest[ii,]$end)/2,decreasing=F)]],,drop=FALSE])
  })
  
  # subset to chromosome 3,4,9,13 and 20 (chr of interest in general)
  chrs <- paste0("chr", 1:22)
  gExp.array <- gExp.array[chrs]
  
  # we want to scale the output to the respective amount of genes and their widths
  widths <- sapply(gExp.array, nrow); widths <- widths/max(widths)*100
  
  # generate the layout accordingly
  plot.layout <- layout(matrix(seq(1,length(gExp.array)),1,length(gExp.array),byrow=TRUE), widths=widths)
  
  # here we define the groups we want to cluster according to
  groups <- valid.clones
  ordered_names <- c()
  
  j <- 1
  # perform hierarchical clustering within each of these groups
  for (j in 1:length(groups)){
    
    # which group?
    print(groups[j])
    
    # get the barcodes from the respective group
    barcodes <- order.data[grep(paste(groups[j],"$", sep=""), order.data$cluster), "cell"]
    
    chr.list <- list()
    if (length(barcodes) > 1){
      # iterate over all chrs and subset the matrix
      
      for (i in 1:length(gExp.array)){

          # matrix 
          matrix <- gExp.array[[i]]
          matrix <- matrix[,colnames(matrix) %in% barcodes]
          
          chr.list[[paste0("chr",i)]] <- matrix

      }
      
      # calculate the col averages here for the subset
      avgd <- do.call(rbind, lapply(names(chr.list),function(nam) {
        d <- chr.list[[nam]]
        d <- colMeans(d)
        d
      }))
      
      # Order cells with hierarchical clustering
      dist.centered.matrix <- dist(as.matrix(t(avgd)), method = "euclidean")
      hc <- hclust(dist.centered.matrix, method = "ward.D2")
      
      # make a vector with the right order of barcodes per group
      ordered_names <- c(ordered_names, hc$labels[hc$order]) 
      
    } else {
      
      # add single barcode to the list
      ordered_names <- c(ordered_names, barcodes)
      cat(paste0("Clone", groups[j], " does only have 1 cell!\n"))
      next
    }
    
  }
  
  # set parameters for the image
  adapted.gExp.list <- gExp.array
  pcol <- c("#053061", "#3C8ABE", "#ffffff", "#F4AA88", "#EB9273", "#D6604D", "#C7433F", "#B51F2E", "#67001F")
  zlim <- c(0, 8) # upper and lower limit
  
  # limit the gExp to maximums according to zlim
  limit.gExp.list <- lapply(names(adapted.gExp.list),function(nam) {
    d <- adapted.gExp.list[[nam]]
    d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
    return(d)
  })
  
  #####################################################################################
  #                   GENERATE THE COLOUR BAR FOR THE HEATMAP
  #####################################################################################
  
  # the index is important for the order so we merge it with the seurat metadata
  heatmap.metadata <- order.data
  colnames(heatmap.metadata) <- c("Cell_barcode", "clone_number")
  rownames(heatmap.metadata) <- heatmap.metadata$Cell_barcode
  heatmap.metadata <- heatmap.metadata[ordered_names,]
  
  # make a df with the index to color it accordingly
  heatmap.col.bar <- data.frame(as.numeric(heatmap.metadata$clone_number),as.numeric(heatmap.metadata$clone_number))
  colnames(heatmap.col.bar)<-c("clone1","clone2")
  
  # get the coloring according to RGs
  cols <- c(brewer.pal(n=length(unique(heatmap.metadata$clone_number)), "Set1"))
  
  #####################################################################################
  #                 GENERATE THE HEATMAP ITSELF AND WRITE TO FILE
  #####################################################################################
  
  png(paste0(output_directory, sample_name, "/", sample_name, "_scDNA_gene_matrix_heatmap.png") , width = 26, height = 12, units = 'in', res = 600)
  par(mfrow=c(1,23), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
  image(t(heatmap.col.bar),xlab="",ylab="", axes=FALSE, col=cols)
  ## plot chromosomes
  box()
  ## plot chromosomes
  for (i in 1:length(limit.gExp.list)){
    message(paste0("chr", i))
    d <- limit.gExp.list[[i]]
    d <- d[, heatmap.metadata$Cell_barcode] 
    # image the respective chr 
    image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=F, main=paste0("chr", i), useRaster = T)
    box()
  }
  dev.off()
  
}


plot_scDNA_clone_heatmap <- function(scDNA.gene.matrix, gene.ordering.file, sample_name, output_directory = "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA/"){
  
  # read in the file of interest  
  norm.gExp.matrix <- scDNA.gene.matrix
  
  # define ploidy to which you want to normalize to
  ploidy <- 2
  gexp.norm <- norm.gExp.matrix
  norm.gexp.norm <- apply(gexp.norm, 2, function(x) round( x*(if(ploidy/median(x) > 1) floor(ploidy/median(x)) else ploidy/median(x))) ) 
  
  # # order the cells according to the clustering file
  # order.data <- read.table(order.file, header = F, sep = ",")
  # colnames(order.data) <- c("cell", "cluster")
  # order.data$cluster <- order.data$cluster + 1 
  # gexp.norm <- norm.gexp.norm[,order.data$cell]
  
  # get the clusters for the respective sample stored as a dataframe with a clone_id column
  # clusters <- unique(order.data$cluster)
  valid.clones <- colnames(norm.gexp.norm)
  
  # only take genes which are present in the matrix
  genes <- makeGRangesFromDataFrame(gene.ordering.file)
  genes <- genes[names(genes) %in% rownames(gexp.norm),]
  
  # generate new object genes.interest with genes and respective coordinates
  genes.interest <- as.data.frame(genes)
  rownames(genes.interest) <- names(genes)
  mat <- gexp.norm
  
  # check for each chromosome if there are genes present
  gExp.array <- tapply(1:nrow(genes.interest),as.factor(genes.interest$seqnames),function(ii) {
    na.omit(mat[rownames(genes.interest)[ii[order((genes.interest[ii,]$start+genes.interest[ii,]$end)/2,decreasing=F)]],,drop=FALSE])
  })
  
  # subset to chromosome 3,4,9,13 and 20 (chr of interest in general)
  chrs <- paste0("chr", 1:22)
  gExp.array <- gExp.array[chrs]
  
  # we want to scale the output to the respective amount of genes and their widths
  widths <- sapply(gExp.array, nrow); widths <- widths/max(widths)*100
  
  # generate the layout accordingly
  plot.layout <- layout(matrix(seq(1,length(gExp.array)),1,length(gExp.array),byrow=TRUE), widths=widths)
  
  # # here we define the groups we want to cluster according to
  # groups <- valid.clones
  # ordered_names <- c()
  # 
  # j <- 1
  # # perform hierarchical clustering within each of these groups
  # for (j in 1:length(groups)){
  #   
  #   # which group?
  #   print(groups[j])
  #   
  #   # get the barcodes from the respective group
  #   barcodes <- order.data[grep(paste(groups[j],"$", sep=""), order.data$cluster), "cell"]
  #   
  #   chr.list <- list()
  #   if (length(barcodes) > 1){
  #     # iterate over all chrs and subset the matrix
  #     
  #     for (i in 1:length(gExp.array)){
  #       
  #       # matrix 
  #       matrix <- gExp.array[[i]]
  #       matrix <- matrix[,colnames(matrix) %in% barcodes]
  #       
  #       chr.list[[paste0("chr",i)]] <- matrix
  #       
  #     }
  #     
  #     # calculate the col averages here for the subset
  #     avgd <- do.call(rbind, lapply(names(chr.list),function(nam) {
  #       d <- chr.list[[nam]]
  #       d <- colMeans(d)
  #       d
  #     }))
  #     
  #     # Order cells with hierarchical clustering
  #     dist.centered.matrix <- dist(as.matrix(t(avgd)), method = "euclidean")
  #     hc <- hclust(dist.centered.matrix, method = "ward.D2")
  #     
  #     # make a vector with the right order of barcodes per group
  #     ordered_names <- c(ordered_names, hc$labels[hc$order]) 
  #     
  #   } else {
  #     
  #     # add single barcode to the list
  #     ordered_names <- c(ordered_names, barcodes)
  #     cat(paste0("Clone", groups[j], " does only have 1 cell!\n"))
  #     next
  #   }
  #   
  # }
  
  # set parameters for the image
  adapted.gExp.list <- gExp.array
  pcol <- c("#053061", "#3C8ABE", "#ffffff", "#F4AA88", "#EB9273", "#D6604D", "#C7433F", "#B51F2E", "#67001F")
  zlim <- c(0, 8) # upper and lower limit
  
  # limit the gExp to maximums according to zlim
  limit.gExp.list <- lapply(names(adapted.gExp.list),function(nam) {
    d <- adapted.gExp.list[[nam]]
    d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
    return(d)
  })
  
  #####################################################################################
  #                   GENERATE THE COLOUR BAR FOR THE HEATMAP
  #####################################################################################
  
  # the index is important for the order so we merge it with the seurat metadata
  # heatmap.metadata <- order.data
  # colnames(heatmap.metadata) <- c("Cell_barcode", "clone_number")
  # rownames(heatmap.metadata) <- heatmap.metadata$Cell_barcode
  # heatmap.metadata <- heatmap.metadata[ordered_names,]
  clone.number <- gsub("Clone", "", valid.clones)
  
  # make a df with the index to color it accordingly
  heatmap.col.bar <- data.frame(as.numeric(clone.number),as.numeric(clone.number))
  colnames(heatmap.col.bar)<-c("clone1","clone2")
  
  # get the coloring according to RGs
  cols <- c(brewer.pal(n=length(unique(clone.number)), "Set1"))
  
  #####################################################################################
  #                 GENERATE THE HEATMAP ITSELF AND WRITE TO FILE
  #####################################################################################
  
  png(paste0(output_directory, sample_name, "/", sample_name, "_scDNA_gene_matrix_heatmap.png") , width = 26, height = 5, units = 'in', res = 600)
  par(mfrow=c(1,23), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
  image(t(heatmap.col.bar),xlab="",ylab="", axes=FALSE, col=cols)
  ## plot chromosomes
  box()
  ## plot chromosomes
  for (i in 1:length(limit.gExp.list)){
    message(paste0("chr", i))
    d <- limit.gExp.list[[i]]
    d <- d[,] 
    # image the respective chr 
    image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=F, main=paste0("chr", i), useRaster = T)
    box()
  }
  dev.off()
  
}


############################################################################
##            TRANSFORM INFERCNV RESULTS TO DISCRETE VALUES
############################################################################

## INPUT
# scRNA.mtx - PATH TO OUTPUT FILE FROM INFERCNV WITH EITHER MODIFIED EXPRESSION OR HMM COPY NUMBER VALUES
# HMM - SPECIFIES WHETHER THE HMM VALUES ARE USED (DEFAULT = FALSE)
## OUTPUT
# scRNA.mtx - MATRIX WITH INTEGER COPY NUMBER VALUES INSTEAD OF ARBITRAY EXPRESSION OR HMM VALUES


# # subset genes to the chromosomes of interest
# gene.GR.df <- gene.GR.df[gene.GR.df$chr_arm %in% variable.chroms,]
# genes <- as.character(gene.GR.df$hgnc_symbol)
# genes <- genes[order(genes)]
# genes <- genes[order(unique(subset.combined.profiles$hgnc_symbol))]

# # subset both the scRNA matrix and the scDNA profile
# subset.t.clone.gene.profile <- t.clone.gene.profile[rownames(t.clone.gene.profile) %in% genes, ]
# subset.t.clone.gene.profile <- t.clone.gene.profile
# subset.scRNA.mtx <- scRNA.mtx[rownames(scRNA.mtx) %in% genes,]
# subset.scRNA.mtx <- scRNA.mtx
# scRNA.mtx <- scRNA.file
transform.inferCNV.mtx <- function(scRNA.mtx, gene.order.file, HMM = FALSE){
  
  # read in the matrix of interest
  scRNA.mtx <- as.matrix(read.table(scRNA.mtx, header = T, sep = " "))
  # scRNA.mtx <- as.matrix(read.table(scRNA.mtx, header = T, sep = "\t"))
  colnames(scRNA.mtx) <- paste0(str_split_fixed(colnames(scRNA.mtx), "_pos", 2)[,1], "-1")
  
  # read in gene_ordering file
  gene.ordering.file <- read.table(gene.order.file, header = F, col.names = c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id"),sep = "\t")
  gene.ordering.file$coordinates <- paste(gene.ordering.file$chr, gene.ordering.file$start, gene.ordering.file$end, sep = ":")
  gene.ordering.file$gene_length <- gene.ordering.file$end - gene.ordering.file$start
  gene.GRange <- makeGRangesFromDataFrame(gene.ordering.file, keep.extra.columns = T)
  
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
  
  # then merge each dataframe information by overlap
  # find the overlapping regions of the WGS segments with the genes
  gene.GR.merge <- mergeByOverlaps(chr.arm.GRange, gene.GRange)
  
  # get the essential columns and make a new dataframe
  gene.GR.df <- data.frame("coordinates" = gene.GR.merge$chr.arm.GRange,
                           "chr_arm" = gene.GR.merge$arm,
                           "hgnc_symbol" = gene.GR.merge$hgnc_symbol,
                           "coordinates" = gene.GR.merge$coordinates,
                           "gene_length" = gene.GR.merge$gene_length,
                           stringsAsFactors = F)
  
  # convert each row
  gene.GR.df$chr_arm <- paste0(gene.GR.df$coordinates.seqnames, gene.GR.df$chr_arm)
  
  # transform to integer values according to standard deviation
  if (HMM == TRUE){
    
    # replace the HMM values from infercnv with integer copy number values
    # the HMM values are between 0.5 and 2.5
    scRNA.mtx[which(scRNA.mtx == 2.5, arr.ind = T)] <- 5
    scRNA.mtx[which(scRNA.mtx == 2.0, arr.ind = T)] <- 4
    scRNA.mtx[which(scRNA.mtx == 1.5, arr.ind = T)] <- 3
    scRNA.mtx[which(scRNA.mtx == 1.0, arr.ind = T)] <- 2
    scRNA.mtx[which(scRNA.mtx == 0.5, arr.ind = T)] <- 1
    
  } else {
    
    chr.arms <- unique(gene.GR.df$chr_arm)
    arm.df <- list()
    i <- 1
    
    for (i in 1:length(chr.arms)){
      
      arm <- chr.arms[i]
      print(arm)
      
      # get the information per arm level
      arm.gene.info <- gene.GR.df[gene.GR.df$chr_arm == arm, c("chr_arm", "coordinates.length", "hgnc_symbol", "gene_length")]
      
      # check how much of the chromosome arm is covered by genes
      covered.chr <- sum(arm.gene.info$gene_length)/unique(arm.gene.info$coordinates.length)
      
      # subset the matrix to the genes on this arm
      arm.scRNA.mtx <- scRNA.mtx[rownames(scRNA.mtx) %in% arm.gene.info$hgnc_symbol, ]
      dim(arm.scRNA.mtx)
    
      
      # calculate the standard deviation
      std.dev <- sd(arm.scRNA.mtx)
      
      # transform the modified expression values to integer values
      # this assumes that we can only reliably call losses or gains
      # inferCNV modified expression values are centered around 1
      arm.scRNA.mtx[which(arm.scRNA.mtx > 1 + std.dev, arr.ind = T)] <- 3
      arm.scRNA.mtx[which(arm.scRNA.mtx <= 1 + std.dev & arm.scRNA.mtx >= 1 - std.dev, arr.ind = T)] <- 2
      arm.scRNA.mtx[which(arm.scRNA.mtx < 1 - std.dev, arr.ind = T)] <- 1
      
      # replace the values in the scRNA.mtx
      scRNA.mtx[rownames(scRNA.mtx) %in% arm.gene.info$hgnc_symbol, ] <- arm.scRNA.mtx
    
    }
  
    majority.scRNA.mtx <- scRNA.mtx
    
    table(majority.scRNA.mtx[,1])
    
    
    
    
    # read in gene_ordering file
    gene.ordering.file <- read.table(gene.order.file, header = F, col.names = c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id"),sep = "\t")
    gene.ordering.file$coordinates <- paste(gene.ordering.file$chr, gene.ordering.file$start, gene.ordering.file$end, sep = ":")
    rownames(gene.ordering.file) <- gene.ordering.file$hgnc_symbol
    
    
    # create the heatmap for the new clustering
    plot_adapted_geneExp_heatmap(scRNA.mtx, gene.ordering.file, sample.tmp, output_directory = out.dir)
    
    
    
  }
  
  
  
  
  # return the addapted matrix
  return(scRNA.mtx)
  
}

############################################################################
##                     VISUALIZE A PSEUDOBULK PROFILE
############################################################################

## INPUT
# gene.profile - DATAFRAME OF CLONES X GENES WITH COPY NUMBER VALUES FROM SCDNA OR SCRNA
# gene.order - DATAFRAME WITH GENE LOCATIONS (SIX COLUMNS) WITH <hgnc_symbol> <chr> <start> <end> <ensembl_gene_id> <coordinates>
## OUTPUT
# pseudobulk_plot - GGPLOT OF CLONES ALONG COPY NUMBER VALUES PER CLONE ALONG THE GENOME
# plot_gene_profile - DATAFRAME USED TO CREATE THE PSEUDOBULK PLOT
plot.pseudobulk.profile <- function(gene.profile, gene.order){
  
  # get the gene profile to plot
  plot.gene.profile <- gene.profile
  
  # subset it to the intersect between the gene.profile and the ordering file
  genes.to.plot <- as.character(gene.order[gene.order$hgnc_symbol %in% colnames(plot.gene.profile), "hgnc_symbol"])
  plot.gene.profile <- plot.gene.profile[, genes.to.plot]
  
  if (ncol(t(gene.profile)) > 1){
    
    # melt the dataframe to make it compatible with ggplot
    plot.gene.profile <- reshape2::melt(plot.gene.profile)
    plot.gene.profile$value <- as.numeric(as.character(plot.gene.profile$value))
    colnames(plot.gene.profile) <- c("clone_id", "hgnc_symbol", "value")
    
  } else {
    
    # melt the dataframe to make it compatible with ggplot
    plot.gene.profile <- reshape2::melt(plot.gene.profile)
    plot.gene.profile$value <- as.numeric(as.character(plot.gene.profile$value))
    plot.gene.profile$clone_id <- "Clone1"
    plot.gene.profile$hgnc_symbol <- rownames(plot.gene.profile)
    
  }
  
  # combine it with the order information
  gene.coordinates <- gene.order[,c("hgnc_symbol", "chr", "start", "end")]
  plot.gene.profile <- merge(plot.gene.profile, gene.coordinates, by = "hgnc_symbol")
  
  # order the gene profile according to the chr and start position
  plot.gene.profile$chr <- factor(plot.gene.profile$chr, levels = c(paste0("chr", c(1:22))))
  plot.gene.profile <- plot.gene.profile[order(plot.gene.profile$chr, plot.gene.profile$start),]
  
  # limit maximum value to 8
  plot.gene.profile[plot.gene.profile$value > 8, "value"] <- 8
  
  # create a pseudobulk plot
  p <- ggplot(plot.gene.profile, aes(x=hgnc_symbol, y=value, group=clone_id)) +
    geom_point(aes(color=clone_id),  position = position_jitter(w = 0, h = 0.1)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme_bw() + 
    scale_color_brewer(palette="Set1") +
    facet_grid(clone_id ~ reorder(chr), scales="free_x", space="free") + 
    ggExtra::removeGrid() + labs(x="Chromosomes", y="Copy number values", color = "Clone") +
    theme(panel.grid.major = element_blank(),
          strip.text.x = element_text(face="bold", size=8, colour = "black",),
          strip.text.y = element_text(face="bold", size=12, colour = "black",),
          strip.background = element_rect(fill="white", colour="black", size=1), 
          axis.text = element_text(colour = "black", size = 8, face = "bold" ),
          axis.title = element_text(colour = "black", size = 12, face = "bold" ),
          plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
          legend.title = element_text(color = "black", size = 12, face = "bold",),
          legend.text = element_text(colour="black", size=8, face="bold"),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))
  
  # return the plot and the dataframe for plotting
  newList <- list("pseudobulk_plot" = p, "plot_gene_profile" = plot.gene.profile)
  return(newList)
}

############################################################################
##          ASSESS THE VARIABLE CHROMOSOMES IN THE SCDNA CLONES
############################################################################

## INPUT
# gene.profile - DATAFRAME OF CLONES X GENES WITH COPY NUMBER VALUES FROM SCDNA OR SCRNA
## OUTPUT
# keep - VECTOR OF CHROMOSOMES WHICH EXHIBIT MORE THAN 5% VARIABILITY
chromosome.variation <- function(gene.profile, cutoff = 0.05){
  
  # read in the hg19 chromosome cooordinates
  x <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/cytoBand.txt.gz", col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
  chr.arm.pos <- x[ , .(length = sum(chromEnd - chromStart)), by = .(chrom, arm = substring(name, 1, 1)) ]
  
  # create start and end coordinates
  chr.arm.pos$start <- 0
  chr.arm.pos[chr.arm.pos$arm == "q", "start"] <- chr.arm.pos[chr.arm.pos$arm == "q", "start"] +chr.arm.pos[chr.arm.pos$arm == "p", "length"]
  chr.arm.pos$end <- chr.arm.pos$length
  chr.arm.pos[chr.arm.pos$arm == "q", "end"] <- chr.arm.pos[chr.arm.pos$arm == "q", "end"] +chr.arm.pos[chr.arm.pos$arm == "p", "end"]
  
  # make a Grange object
  chr.arm.pos.GRange <- makeGRangesFromDataFrame(chr.arm.pos, keep.extra.columns = T)
  gene.profile.GRange <- makeGRangesFromDataFrame(gene.profile, keep.extra.columns = T)
  
  # then merge each dataframe information by overlap
  # find the overlapping regions of the WGS segments with the genes
  variation.GR.merge <- mergeByOverlaps(chr.arm.pos.GRange, gene.profile.GRange)
  
  # get the essential columns and make a new dataframe
  chr.arm.variation <- data.frame("chr_arm" = variation.GR.merge$arm,
                                  "hgnc_symbol" = variation.GR.merge$hgnc_symbol,
                                  "clone_id" = variation.GR.merge$clone_id,
                                  "copy_number_state" = variation.GR.merge$value,
                                  "gene_coordinates" = variation.GR.merge$gene.profile.GRange,
                                  stringsAsFactors = F)
  
  # convert each row
  chr.arm.variation$chr_arm <- as.character(chr.arm.variation$chr_arm)
  chr.arm.variation$chr_arm <- paste0(chr.arm.variation$gene_coordinates.seqnames, chr.arm.variation$chr_arm)
  
  # count the number of copy number alterations per clone
  variation.per.clone <- table(chr.arm.variation$chr_arm, chr.arm.variation$copy_number_state, chr.arm.variation$clone_id)
  
  # calculate the percentage of variation per chromsome and clone
  # and store it in a list
  perc.variation.per.clone <- list()
  
  for (i in 1:length(unique(gene.profile$clone_id))){
    
    # calculate the chr variation
    chr.variation.per.clone <- variation.per.clone[,,i]/rowSums(variation.per.clone[,,i])
    
    # melt it
    melted.var.per.clone <- reshape2::melt(chr.variation.per.clone)
    
    # add clone info
    melted.var.per.clone$cluster <- paste0("Clone", i)
    
    # store the information in a list
    perc.variation.per.clone[[i]] <- melted.var.per.clone
    
  }
  
  # combine it all into one dataframe
  complete.perc.var <- bind_rows(perc.variation.per.clone)
  colnames(complete.perc.var) <- c("chr", "copy_number_values", "value", "clone_id")
  
  # finally, check the variation in each of the chromosomes
  keep <- c()
  for(chr in unique(complete.perc.var$chr)){
    
    # print the chromosome
    print(chr)
    
    # grep that chromosome of interest
    chr.perc.var <- complete.perc.var[grep(paste(chr,"$", sep=""), complete.perc.var$chr),]
    
    # spread out the values and compare
    chr.perc.matrix <- acast(chr.perc.var, clone_id ~ copy_number_values, value.var="value")
    
    # check for the proportion
    if(max(colSds(chr.perc.matrix)) <= cutoff){
      
      next
      
    } else {
      
      keep <- c(keep, chr)  
    }
  }
  
  return(keep)
}

############################################################################
##    CHECK THE VARIABILITY OF THE CHROMOSOMES EVALUATED IN DNA IN RNA
############################################################################

## INPUT
# gene.profile - DATAFRAME OF CLONES X GENES WITH COPY NUMBER VALUES FROM SCDNA OR SCRNA
## OUTPUT
# keep - VECTOR OF CHROMOSOMES WHICH EXHIBIT MORE THAN 5% VARIABILITY

# variable.chroms
# gene.ordering.file
# scRNA.mtx

check.chr.variation <- function(scRNA.mtx, gene.ordering.file, variable.chroms, cutoff = 0.05){
  
  # read in the hg19 chromosome cooordinates
  x <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/cytoBand.txt.gz", col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
  chr.arm.pos <- x[ , .(length = sum(chromEnd - chromStart)), by = .(chrom, arm = substring(name, 1, 1)) ]
  
  # create start and end coordinates
  chr.arm.pos$start <- 0
  chr.arm.pos[chr.arm.pos$arm == "q", "start"] <- chr.arm.pos[chr.arm.pos$arm == "q", "start"] +chr.arm.pos[chr.arm.pos$arm == "p", "length"]
  chr.arm.pos$end <- chr.arm.pos$length
  chr.arm.pos[chr.arm.pos$arm == "q", "end"] <- chr.arm.pos[chr.arm.pos$arm == "q", "end"] +chr.arm.pos[chr.arm.pos$arm == "p", "end"]
  
  # make a Grange object for chr arms
  chr.arm.pos.GRange <- makeGRangesFromDataFrame(chr.arm.pos, keep.extra.columns = T)
  
  # get the gene positions
  gene.profile <- gene.ordering.file[rownames(scRNA.mtx),]
  gene.profile.GRange <- makeGRangesFromDataFrame(gene.profile, keep.extra.columns = T)
  
  # then merge each dataframe information by overlap
  # find the overlapping regions
  variation.GR.merge <- mergeByOverlaps(chr.arm.pos.GRange, gene.profile.GRange)
  
  # get the essential columns and make a new dataframe
  chr.arm.variation <- data.frame("chr_arm" = variation.GR.merge$arm,
                                  "hgnc_symbol" = variation.GR.merge$hgnc_symbol,
                                  "gene_coordinates" = variation.GR.merge$gene.profile.GRange,
                                  stringsAsFactors = F)
  
  # convert each row
  chr.arm.variation$chr_arm <- as.character(chr.arm.variation$chr_arm)
  chr.arm.variation$chr_arm <- paste0(chr.arm.variation$gene_coordinates.seqnames, chr.arm.variation$chr_arm)
  
  rna.pseudobulk <- t(data.frame("Clone1" = rowMeans(scRNA.mtx), stringsAsFactors = F))
  scRNA.pseudobulk$pseudobulk_plot <- plot.pseudobulk.profile(rna.pseudobulk, gene.ordering.file)
  
  # count the number of copy number alterations per clone
  variation.per.clone <- table(chr.arm.variation$chr_arm, chr.arm.variation$copy_number_state, chr.arm.variation$clone_id)
  
  # calculate the percentage of variation per chromsome and clone
  # and store it in a list
  perc.variation.per.clone <- list()
  
  for (i in 1:length(unique(gene.profile$clone_id))){
    
    # calculate the chr variation
    chr.variation.per.clone <- variation.per.clone[,,i]/rowSums(variation.per.clone[,,i])
    
    # melt it
    melted.var.per.clone <- reshape2::melt(chr.variation.per.clone)
    
    # add clone info
    melted.var.per.clone$cluster <- paste0("Clone", i)
    
    # store the information in a list
    perc.variation.per.clone[[i]] <- melted.var.per.clone
    
  }
  
  # combine it all into one dataframe
  complete.perc.var <- bind_rows(perc.variation.per.clone)
  colnames(complete.perc.var) <- c("chr", "copy_number_values", "value", "clone_id")
  
  # finally, check the variation in each of the chromosomes
  keep <- c()
  for(chr in unique(complete.perc.var$chr)){
    
    # print the chromosome
    print(chr)
    
    # grep that chromosome of interest
    chr.perc.var <- complete.perc.var[grep(paste(chr,"$", sep=""), complete.perc.var$chr),]
    
    # spread out the values and compare
    chr.perc.matrix <- acast(chr.perc.var, clone_id ~ copy_number_values, value.var="value")
    
    # check for the proportion
    if(max(colSds(chr.perc.matrix)) <= cutoff){
      
      next
      
    } else {
      
      keep <- c(keep, chr)  
    }
  }
  
  return(keep)
}
############################################################################
##           PERFORM SAMPLE ALIGNMENT BETWEEN SCDNA AND SCRNA
############################################################################

## INPUT
# scRNA.mtx - 
# scDNA.gene.profile - 
# correlation - 
## OUTPUT
# scRNA.mtx <- subset.scRNA.mtx
# scDNA.gene.profile <- subset.t.clone.gene.profile
# scDNA.gene.profile <- random.clone.gene.profile
integrate_scDNA_RNA <- function(scRNA.mtx, scDNA.gene.profile){
  
  # add clones to it
  clone.vector <- colnames(scDNA.gene.profile)
  names(clone.vector) <-colnames(scDNA.gene.profile)
  
  # make all columns numeric and integer values
  scDNA.gene.profile[,c(1:ncol(scDNA.gene.profile))] <- sapply(scDNA.gene.profile[,c(1:ncol(scDNA.gene.profile))], as.numeric)
  scDNA.gene.profile[,c(1:ncol(scDNA.gene.profile))] <- sapply(scDNA.gene.profile[,c(1:ncol(scDNA.gene.profile))], round)
  
  # correlate cells information with all possible clone profiles
  scrna.clone.correlation <- cor(scRNA.mtx, scDNA.gene.profile, method = "pearson")
  
  i <- 1
  scrna.max.correlation <- list()
  for (i in 1:nrow(scrna.clone.correlation)){
    
    if(i %% 100 == 0){
      message(sprintf("%s of %s", i, nrow(scrna.clone.correlation)))
    }
    
    if(any(is.na(scrna.clone.correlation[i,]) == TRUE)) {
      
      # create a cell vector with the barcode, clone_id and pearson.correlation
      cell.vector <- c(colnames(scRNA.mtx)[i], "NA", 0)
      names(cell.vector) <- c("Cell_barcode", "clone_id", "pearson.correlation")
      
    } else {
     
      # create a cell vector with the barcode, clone_id and pearson.correlation
      cell.vector <- c(colnames(scRNA.mtx)[i], colnames(scrna.clone.correlation)[which(scrna.clone.correlation[i,] == max(scrna.clone.correlation[i,]))], abs(max(scrna.clone.correlation[i,])))
      names(cell.vector) <- c("Cell_barcode", "clone_id", "pearson.correlation") 
      cell.vector <- cell.vector[c(1:3)]
    }
    
    # store as a row in alignment.data
    scrna.max.correlation[[i]] <- cell.vector
  }
  
  # combine the information together
  scrna.max.correlation <- bind_rows(scrna.max.correlation)
  alignment.data <- cbind(scrna.max.correlation, scrna.clone.correlation)
  alignment.data <- alignment.data[complete.cases(alignment.data),]
  
  # return the dataframe with the distance metrics and the maximum clone assignment
  return(alignment.data)
  
}

############################################################################
##        CALCULATE A P VALUE BETWEEN EMPIRICAL AND NULL DISTRIBUTION
############################################################################

# observed.clone.data <- subset.t.clone.gene.profile
generate.random.clone.profile <- function(observed.clone.data, sim_clones = 10000){
  
  # create random matrix
  random.clone.data <- matrix(NA, nrow = nrow(observed.clone.data), ncol = sim_clones)
  colnames(random.clone.data) <- c(paste0("Sim_Clone", c(1:ncol(random.clone.data))))
  rownames(random.clone.data) <- rownames(observed.clone.data)
  random.clone.data <- as.data.frame(random.clone.data)
  
  i <- 1
  for (i in 1:nrow(random.clone.data)){
    
    if(i %% 100 == 0){
      message(sprintf("%s of %s", i, nrow(random.clone.data), " Genes"))
    }
    
    random.cn.values <- observed.clone.data[i, ]
    gene.cn.value <- sample(as.vector(random.cn.values[c(1:ncol(observed.clone.data))]), 10000, replace = TRUE)
    names(gene.cn.value) <-  c(paste0("Sim_Clone", c(1:ncol(random.clone.data))))
    
    # add the value to the matrix
    random.clone.data[i,] <- gene.cn.value
  }
  
  return(random.clone.data)
}

############################################################################
##        CALCULATE A P VALUE BETWEEN EMPIRICAL AND NULL DISTRIBUTION
############################################################################

determine.clone.signficance <- function(observed.clone.data, random.clone.data){
  
  # create a matrix for storing the pvalues
  final.clone.correlation <- matrix(NA, nrow = nrow(observed.clone.data), ncol = length(unique(grep("Clone", colnames(observed.clone.data), value = T))))
  colnames(final.clone.correlation) <- c(paste0("Clone", c(1:length(unique(grep("Clone", colnames(observed.clone.data), value = T)))), "_pval"))
  rownames(final.clone.correlation) <- observed.clone.data$Cell_barcode
  
  # get the real and simulated data 
  observed.data.tmp <- observed.clone.data[, c(4:ncol(observed.clone.data))]
  random.data.tmp <- random.clone.data[, c(4:ncol(random.clone.data))]
  
  complete.correlation.dataframe <- list()
  
  j <- 1
  # iterate over each cell and clone and calculate a pvalue
  for(j in 1:nrow(observed.data.tmp)){
    
    if(j %% 100 == 0){
      message(sprintf("%s of %s", j, nrow(observed.data.tmp)))
    }
    
    # get the background correlation
    null.corr <- as.numeric(random.data.tmp[j,])
    
    k <- 1
    for (k in 1:ncol(observed.data.tmp)){
      
      # get the clone
      clone.tmp <- colnames(observed.data.tmp)[k]
      
      # get the correlations
      emp.corr <- as.numeric(observed.data.tmp[j,clone.tmp])
      
      # calculate the empirical significance
      pval <- (1+sum(null.corr > emp.corr))/length(null.corr)
      
      final.clone.correlation[j,k] <- pval
    }
    
    # get the clones with the min pval
    min.pval <- min(final.clone.correlation[j,])
    clones.of.interest <- names(final.clone.correlation[j,grep(min.pval, final.clone.correlation[j,])])
    
    # merge data with final correlation dataframe
    cell.obs.data <- observed.clone.data[j,]
    cell.pval.data <- as.data.frame(final.clone.correlation)[j,]
    cell.data <- cbind(cell.obs.data, cell.pval.data)
    cell.data$num_pval_clones <- length(clones.of.interest)
    cell.data$min_pval <- min.pval
    
    # add to the list
    complete.correlation.dataframe[[j]] <- cell.data
    
  }
  
  # combine them together
  complete.correlation.dataframe <- bind_rows(complete.correlation.dataframe)
  
  # adjust pvalue with Bonferroni correction
  complete.correlation.dataframe$padj <- complete.correlation.dataframe$min_pval * ncol(final.clone.correlation)
  
  # return the final dataframe
  return(complete.correlation.dataframe)
  
}


############################################################################
##                CREATE SCRNA CLONE SPECIFIC PROFILES
############################################################################


# INPUT
# scRNA.gene.mtx - FILE PATH TO THE GENE X CELLS MATRIX GENERATED FROM THE 10X CNV CALLS
# cluster.dataframe - FILE PATH TO THE CLUSTER FILE WHICH INCLUDES THREE COLUMNS: <cell_id> <cluster_id>
## OUTPUT
# t.clone.gene.profile - DATA FRAME WITH GENES X CLONES
create.scRNA.clone.profiles <- function(scRNA.gene.mtx, cluster.dataframe){
  
  # transpose the scRNA gene matrix
  t.scRNA.mtx <- as.data.frame(t(scRNA.gene.mtx))
  t.scRNA.mtx$Cell_barcode <- rownames(t.scRNA.mtx)
  
  # merge with the cluster dataframe information
  scRNA.plot.mtx <- merge(t.scRNA.mtx, cluster.dataframe[,c("Cell_barcode", "clone_id")], by = "Cell_barcode")
  rownames(scRNA.plot.mtx) <- scRNA.plot.mtx$Cell_barcode
  scRNA.plot.mtx$Cell_barcode <- NULL
  
  # create an average gene profile per clone in scRNA
  scRNA.plot.mtx[,c(1:(ncol(scRNA.plot.mtx)-1))] <- sapply(scRNA.plot.mtx[,c(1:(ncol(scRNA.plot.mtx)-1))], as.numeric)
  scRNA.clone.gene.profile <- aggregate(scRNA.plot.mtx[,c(1:(ncol(scRNA.plot.mtx)-1))], list(scRNA.plot.mtx$clone_id), FUN = mean, na.rm=TRUE, na.action=NULL)
  colnames(scRNA.clone.gene.profile)[1] <- "clone_id"
  
  # transpose clone gene profile
  t.scRNA.clone.gene.profile <- as.data.frame(t(scRNA.clone.gene.profile))
  clones <- unique(scRNA.clone.gene.profile$clone_id)
  colnames(t.scRNA.clone.gene.profile) <- clones
  t.scRNA.clone.gene.profile <- t.scRNA.clone.gene.profile[-c(1,nrow(t.scRNA.clone.gene.profile)),]
  
  # return the final vector with genes x clones
  return(t.scRNA.clone.gene.profile)
  
}

############################################################################
##                  PLOT DISCRETIZED INFERCNV HEATMAP
############################################################################
# norm.gExp.matrix <- scRNA.mtx

plot_adapted_geneExp_heatmap <- function(norm.gExp.matrix, gene.ordering.file, sample_name, output_directory = "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA/"){
  
  # center the infercnv output to 0
  gexp.norm <- norm.gExp.matrix
  
  # only take genes which are present in the matrix
  genes <- makeGRangesFromDataFrame(gene.ordering.file)
  genes <- genes[rownames(gexp.norm)]
  
  # generate new object genes.interest with genes and respective coordinates
  genes.interest <- as.data.frame(genes)
  rownames(genes.interest) <- names(genes)
  mat <- gexp.norm
  
  # check for each chromosome if there are genes present
  gExp.array <- tapply(1:nrow(genes.interest),as.factor(genes.interest$seqnames),function(ii) {
    na.omit(mat[rownames(genes.interest)[ii[order((genes.interest[ii,]$start+genes.interest[ii,]$end)/2,decreasing=F)]],,drop=FALSE])
  })
  
  # subset to chromosome 3,4,9,13 and 20 (chr of interest in general)
  chrs <- paste0("chr", 1:22)
  gExp.array <- gExp.array[chrs]
  
  # we want to scale the output to the respective amount of genes and their widths
  widths <- sapply(gExp.array, nrow); widths <- widths/max(widths)*100
  
  # generate the layout accordingly
  plot.layout <- layout(matrix(seq(1,length(gExp.array)),1,length(gExp.array),byrow=TRUE), widths=widths)
  
  # set parameters for the image
  adapted.gExp.list <- gExp.array
  pcol <- c("#053061", "#3C8ABE", "#ffffff", "#F4AA88", "#EB9273", "#D6604D", "#C7433F", "#B51F2E", "#67001F")
  zlim <- c(0, 8) # upper and lower limit
  
  # limit the gExp to maximums according to zlim
  limit.gExp.list <- lapply(names(adapted.gExp.list),function(nam) {
    d <- adapted.gExp.list[[nam]]
    d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
    return(d)
  })
  
  #####################################################################################
  #                 GENERATE THE HEATMAP ITSELF AND WRITE TO FILE
  #####################################################################################
  
  png(paste0(output_directory, sample_name, "/", sample_name, "_inferCNV_adapted_matrix_heatmap.png") , width = 26, height = 12, units = 'in', res = 600)
  par(mfrow=c(1,22), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
  ## plot chromosomes
  for (i in 1:length(limit.gExp.list)){
    message(paste0("chr", i))
    d <- limit.gExp.list[[i]]
    # image the respective chr 
    image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=F, main=paste0("chr", i), useRaster = T)
    box()
  }
  dev.off()
  
}

############################################################################
##                  PLOT FINAL ALIGNED INFERCNV HEATMAP
############################################################################

## INPUT
# norm.gExp.matrix - 
# gene.ordering.file 
# cluster_dataframe
# output_directory
# sample_name
## OUTPUT

# 
# cluster_dataframe <- complete.correlation.dataframe
# scDNA.clones <- unique(scDNA.plot.gene.profile$clone_id)
# chroms <- unique(str_split_fixed(variable.chroms, "p|q", 2)[,1])
# method = "correlation"
# cutoff = "selected_chromosomeArms_uncertainty_0.025_merged"
# sample_name <- sample.tmp
# distance = 0.025
# output_directory = out.dir

plot_geneExp_heatmap <- function(norm.gExp.matrix, gene.ordering.file, cluster_dataframe, scDNA.clones, chroms, distance = 0.05, method = "correlation", cutoff = "_10percent_woChr10_arm", sample_name, output_directory = "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/scRNA_scDNA/"){
  
  # get the clusters for the respective sample stored as a dataframe with a clone_id column
  clusters <- str_split_fixed(unique(cluster_dataframe$clone_id), "Clone", 2)[,2]
  valid.clones <- factor(clusters, levels = c(1:as.numeric(max(unique(clusters)))))
  valid.clones <- valid.clones[order(valid.clones)]
  
  # center the infercnv output to 0
  # gexp.norm <- norm.gExp.matrix -1
  gexp.norm <- norm.gExp.matrix
  
  # only take genes which are present in the matrix
  genes <- makeGRangesFromDataFrame(gene.ordering.file)
  genes <- genes[rownames(gexp.norm)]
  
  # generate new object genes.interest with genes and respective coordinates
  genes.interest <- as.data.frame(genes)
  rownames(genes.interest) <- names(genes)
  mat <- gexp.norm
  
  # check for each chromosome if there are genes present
  gExp.array <- tapply(1:nrow(genes.interest),as.factor(genes.interest$seqnames),function(ii) {
    na.omit(mat[rownames(genes.interest)[ii[order((genes.interest[ii,]$start+genes.interest[ii,]$end)/2,decreasing=F)]],,drop=FALSE])
  })
  
  # subset to chromosome 3,4,9,13 and 20 (chr of interest in general)
  chrs <- paste0("chr", 1:22)
  gExp.array <- gExp.array[chrs]
  
  # we want to scale the output to the respective amount of genes and their widths
  widths <- sapply(gExp.array, nrow); widths <- widths/max(widths)*100
  
  # generate the layout accordingly
  plot.layout <- layout(matrix(seq(1,length(gExp.array)),1,length(gExp.array),byrow=TRUE), widths=widths)
  
  # here we define the groups we want to cluster according to
  groups <- valid.clones
  chr.list <- list()
  ordered_names <- c()
  
  # perform hierarchical clustering within each of these groups
  for (j in 1:length(groups)){
    
    # which group?
    print(groups[j])
    
    # get the barcodes from the respective group
    barcodes <- cluster_dataframe[grep(paste("Clone",groups[j],"$", sep=""), cluster_dataframe$clone_id), "Cell_barcode"]
    
    if (length(barcodes) > 1){
      # iterate over all chrs and subset the matrix
      
      for (i in 1:length(gExp.array)){
        
        if (paste0("chr", i) %in% chroms){
          
          # matrix 
          matrix <- gExp.array[[i]]
          matrix <- matrix[,colnames(matrix) %in% barcodes]
          
          chr.list[[paste0("chr",i)]] <- matrix
        } else {
          next
        }
        
      }
      
      # calculate the col averages here for the subset
      avgd <- do.call(rbind, lapply(names(chr.list),function(nam) {
        d <- chr.list[[nam]]
        d <- colMeans(d)
        d
      }))
      
      # Order cells with hierarchical clustering
      dist.centered.matrix <- dist(as.matrix(t(avgd)), method = "euclidean")
      hc <- hclust(dist.centered.matrix, method = "ward.D2")
      
      # make a vector with the right order of barcodes per group
      ordered_names <- c(ordered_names, hc$labels[hc$order]) 
    } else {
      
      # add single barcode to the list
      ordered_names <- c(ordered_names, barcodes)
      cat(paste0("Clone", groups[j], " does only have 1 cell!\n"))
      next
    }
    
  }
  
  # set parameters for the image
  adapted.gExp.list <- gExp.array
  pcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(256) # color panel
  pcol_corr <- colorRampPalette(c("darkred", "lightgrey"))(256)
  zlim <- c(0.9, 1.1)
  zlim_corr <- c(0, 0.025)
  
  # limit the gExp to maximums according to zlim
  limit.gExp.list <- lapply(names(adapted.gExp.list),function(nam) {
    d <- adapted.gExp.list[[nam]]
    d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
    return(d)
  })
  
  show_col(pcol_corr)
  
  #####################################################################################
  #                   GENERATE THE COLOUR BAR FOR THE HEATMAP
  #####################################################################################
  
  # first, remove the clone id and create and index
  cluster_dataframe$clone_number <- as.numeric(str_split_fixed(cluster_dataframe$clone_id, "Clone", 2)[,2])
  cluster_dataframe <- cluster_dataframe[cluster_dataframe$padj <= 0.05, ]
  
  # the index is important for the order so we merge it with the seurat metadata
  heatmap.metadata <- cluster_dataframe[,c("Cell_barcode", "distance", "clone_number")]
  
  # order according to minimum and maximum correlation
  correlation.order <- c()
  for (clone in unique(heatmap.metadata$clone_number)){
    
    print(clone)
    data <- heatmap.metadata[heatmap.metadata$clone_number == clone, ]
    barcodes <- data[order(data$distance, decreasing = T), "Cell_barcode"]
    correlation.order <- c(correlation.order, barcodes)
    
  }
  
  rownames(heatmap.metadata) <- heatmap.metadata$Cell_barcode
  heatmap.metadata <- heatmap.metadata[correlation.order,]
  heatmap.metadata <- heatmap.metadata[order(heatmap.metadata$clone_number),]
  
  # make a df with the index to color it accordingly
  heatmap.col.bar <- data.frame(as.numeric(heatmap.metadata$clone_number),as.numeric(heatmap.metadata$clone_number))
  colnames(heatmap.col.bar)<-c("clone1","clone2")
  
  # get the coloring according to RGs
  cols <- c(brewer.pal(n=length(scDNA.clones), "Set1"))
  # cols <- cols[scDNA.clones %in% unique(cluster_dataframe$clone_id)]
  
  # get the correlation
  correlation.vector <- data.frame(as.numeric(heatmap.metadata$distance),as.numeric(heatmap.metadata$distance))
  
  #####################################################################################
  #                 GENERATE THE HEATMAP ITSELF AND WRITE TO FILE
  #####################################################################################
  
  png(paste0(output_directory, sample_name, "/", sample_name, "_inferCNV_scDNA_", cutoff,"_", method, "_heatmap.png") , width = 26, height = 12, units = 'in', res = 600)
  par(mfrow=c(1,24), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
  image(t(heatmap.col.bar),xlab="",ylab="", axes=FALSE, col=cols)
  ## plot chromosomes
  box()
  for (i in 1:length(limit.gExp.list)){
    message(chrs[i])
    d <- limit.gExp.list[[i]]
    d <- d[, heatmap.metadata$Cell_barcode] 
    # image the respective chr 
    image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=F, main=chrs[i], useRaster = T)
    box()
  }
  image(t(correlation.vector),xlab="",ylab="", col=pcol_corr, zlim=zlim_corr, axes=FALSE, main="Distance")
  box()
  dev.off()
  
  png(paste0(output_directory, sample_name, "/", sample_name, "_inferCNV_scDNA_clone_heatmap_legend.png") , width = 13, height = 4, units = 'in', res = 600, type = "cairo")
  par(mfrow=c(1,24), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.5)
  image.plot(legend.only=TRUE, zlim=zlim, col = pcol, legend.shrink = 0.4, legend.width = 10, smallplot=  c(.65, .9, .5, .8)) # first argument of small plot = width of legend, third is making the length, the higher the smaller
  dev.off()
  
  png(paste0(output_directory, sample_name, "/", sample_name, "_inferCNV_scDNA_clone_heatmap_uncertainty_legend.png") , width = 13, height = 4, units = 'in', res = 600, type = "cairo")
  par(mfrow=c(1,24), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.5)
  image.plot(legend.only=TRUE, zlim=zlim_corr, col = pcol_corr, legend.shrink = 0.4, legend.width = 10, smallplot=  c(.65, .9, .5, .8)) # first argument of small plot = width of legend, third is making the length, the higher the smaller
  dev.off()
  
}







