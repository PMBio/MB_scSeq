#############################################################################################################################
##                                                                                                                      
##  CORRELATE THE COPY NUMBER CLONES FROM SCDNA WITH THE COPY NUMBER FROM THE HMM FROM INFERCNV
##                                                                                                                      
##  Date: 12 NOVEMBER 2020                                                                                                                   
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
                      "Matrix", "devtools", "Matrix.utils", "matrixStats", "readr", "magrittr", "Signac", "BiocManager", "gridExtra",
                      "biomaRt", "httr", "ComplexHeatmap", "fields", "data.table", "dbscan", "fpc", "Seurat", "mclust")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
if(!"Matrix.utils" %in% installed.packages()[,"Package"]) remotes::install_github("cvarrichio/Matrix.utils")

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

############################################################################
##                  READ IN ALL THE FILES OF INTEREST
############################################################################

# negate %in%
"%ni%" <- Negate("%in%")

# get the functions to perform the alignment
source("~/scDNASeq_natgen_notebooks/infercnv/revision/scrna_analysis/scDNA_scRNA_alignment_functions.R")
debug(generate.random.clone.profile)
# set sample id
sample.tmp <- "MB243-Nuclei_clones"
sample.short <- "MB243-Nuclei"

# specify the input directory
input.dir <- "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB/freeze"

# specify the output directory
out.dir <- "~/scDNASeq_natgen_notebooks/infercnv/infercnv_MB/scRNA_scDNA/"
dir.create(out.dir)

# create a sample directory
dir.create(paste0(out.dir, sample.tmp))

# GET ALL THE SAMPLE INFORMATION REQUIRED FOR 10X CNV DATA
# list all scDNA gene x cell matrices
scDNA.files <- list.files("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/scDNA_gene_matrices", pattern = "_scDNA_gene_matrix_new.txt", recursive = T, full.names = T)
scDNA.file <- scDNA.files[grep(sample.tmp, scDNA.files)]

## GET ALL THE SAMPLE INFORMATION REQUIRED FOR 10X RNA DATA FROM INFERCNV
# list all scRNA matrices
scRNA.files <- list.files("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/infercnv_MB", pattern = "infercnv.median_filtered.observations.txt", full.names = T, recursive = T)
scRNA.files <- scRNA.files[grep("freeze", scRNA.files)]
scRNA.files <- scRNA.files[grep("s_tumor/|s_mouse", scRNA.files)]
scRNA.file <- scRNA.files[grep(sample.short, scRNA.files)]

# for the G&T-seq data
# scRNA.file <- file.path("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_MB/GTseq_normal/results/infercnv.median_filtered.observations.txt")

# and gene ordering files
gene.order.files <- list.files(input.dir, pattern = "gene_ordering_file.txt", recursive = T, full.names = T)
# gene.order.files <- gene.order.files[grep("Nuclei", gene.order.files)] 
gene.order.files <- gene.order.files[grep("tumor", gene.order.files)] 
gene.order.file <- gene.order.files[grep(sample.short, gene.order.files)]

############################################################################
##                  READ IN ALL THE FILES OF INTEREST
############################################################################

# read in gene_ordering file
gene.ordering.file <- read.table(gene.order.file, header = F, col.names = c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id"),sep = "\t")
gene.ordering.file$coordinates <- paste(gene.ordering.file$chr, gene.ordering.file$start, gene.ordering.file$end, sep = ":")
rownames(gene.ordering.file) <- gene.ordering.file$hgnc_symbol

# subset to the genes of interest
genes <- as.character(gene.ordering.file$hgnc_symbol)
genes <- genes[order(genes)]

# read in matrix form scDNA
# read in scDNA matrix
scdna.matrix <- read.table(scDNA.file, sep = "\t")
colnames(scdna.matrix) <- paste0("Clone", c(1:ncol(scdna.matrix)))

# define ploidy to which you want to normalize to
ploidy <- 2

# normalize for ploidy before calculating distance
scdna.matrix <- apply(scdna.matrix, 2, function(x) round( x*(if(ploidy/median(x) > 1) floor(ploidy/median(x)) else ploidy/median(x))) ) 

# visualise scDNA gene heatmap
plot_scDNA_clone_heatmap(scdna.matrix, gene.ordering.file, sample_name = sample.tmp, output_directory = out.dir)
  
# transform the inferCNV matrix into a compatible format to the 10x CNV data
scRNA.mtx <- transform.inferCNV.mtx(scRNA.file, gene.order.file = gene.order.file, HMM = FALSE)

# subset both the scRNA matrix and the scDNA profile
scdna.matrix <- scdna.matrix[rownames(scdna.matrix) %in% genes, ]
scRNA.mtx <- scRNA.mtx[rownames(scRNA.mtx) %in% genes,]

# subset clone profiles to genes which are present in the scRNA data
scdna.matrix <- scdna.matrix[rownames(scdna.matrix) %in% rownames(scRNA.mtx),]
scRNA.mtx <- scRNA.mtx[rownames(scRNA.mtx) %in% rownames(scdna.matrix),]

# order scRNA matrix according to scDNA.matrix
scRNA.mtx <- scRNA.mtx[order(rownames(scRNA.mtx)),]

# check the dimensions of both matrices
# scdna.matrix <- scdna.matrix[,-5] # STP-Nuclei
dim(scdna.matrix)
dim(scRNA.mtx)

############################################################################
##                  PERFORM AND VISUALISE UMAP
############################################################################

# create a pseudobulk visualization for each clone along the genome
scDNA.pseudobulk <- plot.pseudobulk.profile(t(scdna.matrix), gene.ordering.file)
scDNA.pseudobulk$pseudobulk_plot
ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_WG_scDNA_geneProfile.pdf"), height = 10, width = 26, dpi = 600)

# create the heatmap for the new clustering
plot_adapted_geneExp_heatmap(scRNA.mtx, gene.ordering.file, sample.tmp, output_directory = out.dir)

# determine gene-level CNV differences between clones
clone.gene.profile <- scDNA.pseudobulk$plot_gene_profile
clone.gene.profile$length <- clone.gene.profile$end - clone.gene.profile$start

clone.list <- unique(clone.gene.profile$clone_id)
clone <- "Clone1"
complete.gain.profiles <- list()
complete.loss.profiles <- list()

# iterate over each clone and determine the list of genes which are lost or gained
# compared to the other clones
for (clone in clone.list){

  # which clone are we starting with
  print(clone)
  clone.profile <- scdna.matrix

  # get the clones to compare to
  clones.to.compare <- colnames(clone.profile)
  clones.to.compare <- clones.to.compare[-grep(clone, clones.to.compare)]

  # initialize lists to store comparisons
  gain.profiles <- data.frame("hgnc_symbol" = NA, "clone_id" = NA, "value" = NA, "chr" = NA, "start" = NA, "end" = NA, "length" = NA, "type" = NA, "comparison" = NA)
  loss.profiles <- data.frame("hgnc_symbol" = NA, "clone_id" = NA, "value" = NA, "chr" = NA, "start" = NA, "end" = NA, "length" = NA, "type" = NA, "comparison" = NA)
  
  i <- 1
  # ITERATE OVER EACH CLONE AND CHECK THE UNIQUE DIFFERENCES PER CLONE
  for (i in 1:length(clones.to.compare)){

    # get the paired clone and determine the difference per gene
    paired.clone <- clones.to.compare[i]
    gene.cnv.diff <- clone.profile[,clone] - clone.profile[,paired.clone]

    # get the genes with unique losses and gains
    loss.genes <- names(gene.cnv.diff)[gene.cnv.diff == -1 & clone.profile[, paired.clone] != 3]
    gain.genes <- names(gene.cnv.diff)[gene.cnv.diff == 1 & clone.profile[, paired.clone] != 1]

    if (length(gain.genes) > 0){

      gain.clone.gene.profile <- clone.gene.profile[clone.gene.profile$clone_id == clone & clone.gene.profile$hgnc_symbol %in% gain.genes,]
      gain.clone.gene.profile$type <- "gain"
      gain.clone.gene.profile$comparison <- paired.clone
      gain.profiles <- rbind(gain.profiles, gain.clone.gene.profile)

    }

    if (length(loss.genes) > 0){

      loss.clone.gene.profile <- clone.gene.profile[clone.gene.profile$clone_id == clone & clone.gene.profile$hgnc_symbol %in% loss.genes,]
      loss.clone.gene.profile$type <- "loss"
      loss.clone.gene.profile$comparison <- paired.clone
      loss.clone.gene.profile$hgnc_symbol <- as.character(loss.clone.gene.profile$hgnc_symbol)
      loss.profiles <- rbind(loss.profiles, loss.clone.gene.profile)
    }

  }

  gain.profiles <- gain.profiles[-1, ]
  loss.profiles <- loss.profiles[-1, ]

  if (nrow(gain.profiles) > 0){
    complete.gain.profiles[[clone]] <- gain.profiles
  }

  if (nrow(loss.profiles) > 0){
    complete.loss.profiles[[clone]] <- loss.profiles
  }

}

complete.gain.profiles <- bind_rows(complete.gain.profiles)
complete.loss.profiles <- bind_rows(complete.loss.profiles)
combined.profiles <- rbind(complete.gain.profiles, complete.loss.profiles)

# calculate median gene size
combined.profiles <- combined.profiles %>% group_by(clone_id, comparison, type) %>% mutate(median_length = round(median(length)))

# check the density distribution of correlation coefficients
ggplot(combined.profiles, aes(as.numeric(length), group=type, fill=type)) +
  geom_histogram(bins = 100, color = "black", alpha = 0.25) +
  ylab("Density") + xlab("Gene size [bp]") +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  facet_grid(comparison ~ clone_id, space="free") +
  scale_x_continuous(labels = scales::comma) +
  geom_vline(data = combined.profiles[combined.profiles$type == "loss",], aes(xintercept = median_length), colour = 'darkblue', size = 1) +
  geom_text(data = combined.profiles[combined.profiles$type == "loss",], aes(1000000, 100, label = median_length, hjust = 1.5, vjust = -1), colour = 'darkblue', size = 4) +
  geom_vline(data = combined.profiles[combined.profiles$type == "gain",], aes(xintercept = median_length), colour = 'darkred', size = 1) +
  geom_text(data = combined.profiles[combined.profiles$type == "gain",], aes(1000000, 300, label = median_length, hjust = 1.5, vjust = -1), colour = 'darkred', size = 4) +
  scale_fill_brewer(palette="Set1") +
  theme(strip.text.x = element_text(face="bold", size=12, colour = "black",),
        strip.text.y = element_text(face="bold", size=12, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        legend.title = element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"),
        legend.position = "none")
ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_gene_level_CNV_difference.pdf"), height = 20, width = 20, dpi = 600)

############################################################################
##      FILTER THE GENES OF INTEREST TO GENES ON VARIABLE CHROMOSOMES
############################################################################

# # determine the variable chromosomes
# variable.chroms <- chromosome.variation(scDNA.pseudobulk$plot_gene_profile, cutoff = 0.15)

# this is the best possible solution

variable.chroms <- switch(sample.tmp, 
                          "STP-PDX_clones"= c("chr2p", "chr2q", "chr4q", "chr5p",  "chr5q", "chr12p", "chr10q", "chr15p", "chr16p", "chr16q", "chr17p", "chr17q", "chr19p", "chr19q", "chr22p", "chr22q"),
                          "STP-Nuclei_clones" = c("chr3p", "chr3q", "chr4q", "chr5p","chr9p", "chr13p", "chr13q",  "chr14p", "chr14q",  "chr17p", "chr19q", "chr20q"),
                          "MB243-Nuclei_clones" = c("chr1q", "chr2p", "chr2q", "chr4p", "chr4q", "chr5p", "chr5q", "chr7p", "chr7q", "chr18p", "chr8q",  "chr17p", "chr17q", "chr19p", "chr19q"),
                          "ST1R-PDX_clones" = c("chr8p", "chr8q",  "chr11p", "chr12p", "chr12q"))
  
  # c("chr2p", "chr2q", "chr4q", "chr5p",  "chr5q", "chr12p", "chr10q", "chr15p", "chr16p", "chr16q", "chr17p", "chr17q", "chr19p", "chr19q", "chr22p", "chr22q") # STP-PDX
# variable.chroms <- c("chr3p", "chr3q", "chr4q", "chr5p","chr9p", "chr13p", "chr13q",  "chr14p", "chr14q",  "chr17p", "chr19q", "chr20q") # STP-Nuclei 
# variable.chroms <- c("chr1q", "chr2p", "chr2q", "chr4p", "chr4q", "chr5p", "chr5q", "chr7p", "chr7q", "chr18p", "chr8q",  "chr17p", "chr17q", "chr19p", "chr19q") # MB243
# variable.chroms <- c("chr8p", "chr8q",  "chr11p", "chr12p", "chr12q") # ST1R-PDX

# read in gene_ordering file
gene.ordering.file <- read.table(gene.order.file, header = F, col.names = c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id"),sep = "\t")
gene.ordering.file$coordinates <- paste(gene.ordering.file$chr, gene.ordering.file$start, gene.ordering.file$end, sep = ":")
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
                         stringsAsFactors = F)

# convert each row
gene.GR.df$chr_arm <- paste0(gene.GR.df$coordinates.seqnames, gene.GR.df$chr_arm)

# subset genes to the chromosomes of interest
gene.GR.df <- gene.GR.df[gene.GR.df$chr_arm %in% variable.chroms,]
genes <- as.character(gene.GR.df$hgnc_symbol)
genes <- genes[order(genes)]
# genes <- genes[order(unique(subset.combined.profiles$hgnc_symbol))]

# subset both the scRNA matrix and the scDNA profile
subset.t.clone.gene.profile <- scdna.matrix[rownames(scdna.matrix) %in% genes, ]
# subset.t.clone.gene.profile <- t.clone.gene.profile
subset.scRNA.mtx <- scRNA.mtx[rownames(scRNA.mtx) %in% genes,]
# subset.scRNA.mtx <- scRNA.mtx

# check the dimensions of both matrices
dim(subset.t.clone.gene.profile)
dim(subset.scRNA.mtx)

# create a pseudobulk visualization for each clone along the genome
scDNA.pseudobulk.subset <- plot.pseudobulk.profile(t(subset.t.clone.gene.profile), gene.ordering.file)
scDNA.pseudobulk.subset$pseudobulk_plot
ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_WG_scDNA_geneProfile_selected_chromosomeArms.pdf"), height = 10, width = 18, dpi = 600)

############################################################################
##              CORRELATE CELLS FROM RNA TO CLONES FROM DNA 
############################################################################

# # remove clone 5 in STP-Nuclei as it is diploid
if(sample.tmp == "STP-Nuclei_clones"){
  # browser()
  
  subset.t.clone.gene.profile <- subset.t.clone.gene.profile[,-grep("Clone5", colnames(subset.t.clone.gene.profile))]
}

# create the correlation dataframe as well as the distance dataframe
correlation.data.frame <- integrate_scDNA_RNA(subset.scRNA.mtx, subset.t.clone.gene.profile)
correlation.data.frame$pearson.correlation <- (as.numeric(as.character(correlation.data.frame$pearson.correlation)))

# only keep cells in the RNA matrix which could get matched
subset.scRNA.mtx <- subset.scRNA.mtx[, correlation.data.frame$Cell_barcode %in% colnames(subset.scRNA.mtx)]

# write the final correlation to a file
write.table(correlation.data.frame, paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_clones_selected_chromosomeArms_pearson_correlation.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

# check the density distribution of correlation coefficients
ggplot(correlation.data.frame, aes(as.numeric(pearson.correlation), group=clone_id, fill=clone_id)) +
  geom_histogram(bins = 100) + 
  ylab("Density") + xlab("Pearson Correlation - R") +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  facet_grid( ~ clone_id, space="free") + 
  xlim(0,1) + 
  scale_fill_brewer(palette="Set1") +
  theme(strip.text.x = element_text(face="bold", size=12, colour = "black",),
        strip.text.y = element_text(face="bold", size=12, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        legend.title = element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"),
        legend.position = "none")
ggsave(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_selected_chromosomeArms_pearson_correlation_density.pdf"), height = 3, width = ncol(subset.t.clone.gene.profile)*3, dpi = 600)


# overall.correlation.data <- reshape2::melt(correlation.data.frame[, c("Cell_barcode", "Clone1", "Clone2", "Clone3", "Clone4", "Clone5")]) # STP-PDX
# overall.correlation.data <- reshape2::melt(correlation.data.frame[, c("Cell_barcode", "Clone1", "Clone2", "Clone3")]) # ST1R-PDX
# overall.correlation.data <- reshape2::melt(correlation.data.frame[, c("Cell_barcode", "Clone1", "Clone2", "Clone3",  "Clone4")])

overall.correlation.data <- switch(sample.tmp, 
                                   "STP-PDX_clones"= reshape2::melt(correlation.data.frame[, c("Cell_barcode", "Clone1", "Clone2", "Clone3", "Clone4", "Clone5")]),
                                   "STP-Nuclei_clones" = reshape2::melt(correlation.data.frame[, c("Cell_barcode", "Clone1", "Clone2", "Clone3",  "Clone4")]),
                                   "MB243-Nuclei_clones" = reshape2::melt(correlation.data.frame[, c("Cell_barcode", "Clone1", "Clone2", "Clone3",  "Clone4")]),
                                   "ST1R-PDX_clones" = reshape2::melt(correlation.data.frame[, c("Cell_barcode", "Clone1", "Clone2", "Clone3")]))


# check the density distribution of correlation coefficients
ggplot(overall.correlation.data, aes(as.numeric(value), group=variable, fill=variable)) +
  geom_histogram(bins = 100) + 
  ylab("Density") + xlab("Pearson Correlation - abs(R)") +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  facet_grid( ~ variable, space="free") + 
  xlim(0,1) + 
  scale_fill_brewer(palette="Set1") +
  theme(strip.text.x = element_text(face="bold", size=12, colour = "black",),
        strip.text.y = element_text(face="bold", size=12, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        legend.title = element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"),
        legend.position = "none")
ggsave(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_selected_chromosomeArms_pearson_correlation_all_clone_density.pdf"), height = 3, width = ncol(subset.t.clone.gene.profile)*3, dpi = 600)

# check how many cells are there per cluster
print(table(correlation.data.frame$clone_id))

#############################################################################
##              PERFORM PERMUTATION TO GENERATE NULL DISTRIBUTION
############################################################################

# create random clone profiles
random.clone.gene.profile <- generate.random.clone.profile(subset.t.clone.gene.profile)

# compare the random data to the clones of interest
random.correlation.data.frame <- integrate_scDNA_RNA(subset.scRNA.mtx, random.clone.gene.profile)

# check the density distribution of correlation coefficients
ggplot(random.correlation.data.frame, aes(as.numeric(pearson.correlation))) +
  geom_histogram(bins = 100, color = "black") + 
  ylab("Density") + xlab("Pearson Correlation - R") +
  theme_classic() +
  xlim(0,1) + 
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.title = element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"), 
        legend.position = "none")
ggsave(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_random_scDNA_clones_data.pdf"), height = 3, width = 3.5, dpi = 600)

random.correlation.data.frame$Clone <- correlation.data.frame[rownames(random.correlation.data.frame),"clone_id"]

# check the density distribution of correlation coefficients
ggplot(random.correlation.data.frame, aes(as.numeric(pearson.correlation), group=Clone, fill=Clone)) +
  geom_histogram(bins = 100) + 
  ylab("Density") + xlab("Pearson Correlation - R") +
  theme_classic() +
  xlim(0,1) + 
  scale_fill_brewer(palette="Set1") +
  facet_grid( ~ Clone, space="free") + 
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.title = element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"), 
        legend.position = "none")
ggsave(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_random_scDNA_clones_data_byclone.pdf"), height = 3, width = ncol(subset.t.clone.gene.profile)*3, dpi = 600)

random.correlation.data.frame$Clone <- NULL

# determine the empirical significance
complete.correlation.dataframe <- determine.clone.signficance(correlation.data.frame, random.correlation.data.frame)

# remove clones which do not have 5 cells
matched.clones <- unique(complete.correlation.dataframe$clone_id)
matched.clones <- matched.clones[order(matched.clones)]
valid.clones <- matched.clones[table(complete.correlation.dataframe$clone_id) >= 5]
complete.correlation.dataframe <- complete.correlation.dataframe[complete.correlation.dataframe$clone_id %in% valid.clones,]

# visualise p-values per clone and cell
overall.pval.correlation.data <- switch(sample.tmp, 
                                        "STP-PDX_clones"= reshape2::melt(complete.correlation.dataframe[, c("Cell_barcode", paste0(c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5"), "_pval"))]),
                                        "STP-Nuclei_clones" = reshape2::melt(complete.correlation.dataframe[, c("Cell_barcode", paste0(c("Clone1", "Clone2", "Clone3", "Clone4"), "_pval"))]),
                                        "MB243-Nuclei_clones" = reshape2::melt(complete.correlation.dataframe[, c("Cell_barcode", paste0(c("Clone1", "Clone2", "Clone3", "Clone4"), "_pval"))]),
                                        "ST1R-PDX_clones" = reshape2::melt(complete.correlation.dataframe[, c("Cell_barcode", paste0(c("Clone1", "Clone2", "Clone3"), "_pval"))]))

  # reshape2::melt(complete.correlation.dataframe[, c("Cell_barcode", paste0(c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5"), "_pval"))]) # STP-PDX
# overall.pval.correlation.data <- reshape2::melt(complete.correlation.dataframe[, c("Cell_barcode", paste0(c("Clone1", "Clone2", "Clone3"), "_pval"))]) # ST1R-PDX
# overall.pval.correlation.data <- reshape2::melt(complete.correlation.dataframe[, c("Cell_barcode", paste0(c("Clone1", "Clone2", "Clone3", "Clone4"), "_pval"))]) # STP-Nuclei

# check the density distribution of correlation coefficients
oneminussqrt <- trans_new(name = "oneminussqrt", transform = function(x) 1-sqrt(x), inverse = function(x) (1-x)^2, domain=c(0,1))


ggplot(complete.correlation.dataframe, aes(x=-log10(as.numeric(padj)), color=clone_id)) + 
  stat_ecdf(geom = "step") + scale_y_continuous(trans =oneminussqrt ) + 
  # geom_histogram(bins = 100) +
  ylab("Fraction") + xlab("Negative log10 Matching p value") +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  # facet_grid(.~clone_id) + 
  geom_vline(xintercept = -log10(0.05)) +
  theme(strip.text.x = element_text(face="bold", size=12, colour = "black",),
        strip.text.y = element_text(face="bold", size=12, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        legend.title = element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"),
        legend.position = "none")
ggsave(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_seleected_chromosomeArms_p_value_all_clone_ecdf.pdf"), height = 3, width = 4, dpi = 600)




############################################################################
##        CREATE A PSEUDOBULK PLOT ACROSS THE INTEGRATED CLONES
############################################################################

complete.correlation.dataframe <- read.table(paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_clones_filtered_cells.txt"), header = T, sep = "\t")

i <- 1
complete.correlation.dataframe$cell_2nd_best_match <- NA
complete.correlation.dataframe$clone_2nd_best_match <- NA
complete.correlation.dataframe$distance <- NA
for (i in 1:nrow(complete.correlation.dataframe)){
  
  cell.data.tmp <- complete.correlation.dataframe[i,]
  cell_best_matched <- cell.data.tmp$pearson.correlation
  cell_2nd_best_match <- switch(sample.tmp, 
         "STP-PDX_clones"= reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3", "Clone4", "Clone5")]),
         "STP-Nuclei_clones" = cell_2nd_best_match <- reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3", "Clone4")]),
         "MB243-Nuclei_clones" = cell_2nd_best_match <- reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3", "Clone4")]),
         "ST1R-PDX_clones" = reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3")]))
  
  
  # cell_2nd_best_match <- reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3", "Clone4")])
  # cell_2nd_best_match <- reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3")])
  clone_2nd_best_match <- as.character(cell_2nd_best_match[order(cell_2nd_best_match$value, decreasing = T), "variable"])[2]
  cell_2nd_best_match_corr <- as.numeric(cell_2nd_best_match[order(cell_2nd_best_match$value, decreasing = T), "value"])[2]
  distance <- as.numeric(cell_2nd_best_match[order(cell_2nd_best_match$value, decreasing = T), "value"])[1] - as.numeric(cell_2nd_best_match[order(cell_2nd_best_match$value, decreasing = T), "value"])[2]
  
  complete.correlation.dataframe[i, "cell_2nd_best_match"] <- cell_2nd_best_match_corr
  complete.correlation.dataframe[i, "clone_2nd_best_match"] <- clone_2nd_best_match
  complete.correlation.dataframe[i, "distance"] <- distance

}

write.table(complete.correlation.dataframe, paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_clones_filtered_cells.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

# check the density distribution of uncertainty per clone
ggplot(complete.correlation.dataframe, aes(as.numeric(distance), group=clone_id, fill=clone_id)) +
  geom_histogram(bins = 100) + 
  ylab("Density") + xlab("Distance R2_best - R2_2nd_best") +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  facet_grid( ~ clone_id, space="free") + 
  xlim(0,0.5) + 
  scale_fill_brewer(palette="Set1") +
  # geom_vline(xintercept = 0.05) +
  geom_vline(xintercept = 0.025) +
  # geom_vline(xintercept = 0.0125) +
  theme(strip.text.x = element_text(face="bold", size=12, colour = "black",),
        strip.text.y = element_text(face="bold", size=12, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        legend.title = element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"),
        legend.position = "none")
ggsave(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_selected_chromosomeArms_uncertainty_density.pdf"), height = 3, width = ncol(subset.t.clone.gene.profile)*3, dpi = 600)

# subset to cells that pass filter
complete.correlation.dataframe <- complete.correlation.dataframe[complete.correlation.dataframe$padj <= 0.05,]


# check the density distribution of uncertainty per clone
ggplot(complete.correlation.dataframe, aes(as.numeric(distance), group=clone_id, fill=clone_id)) +
  geom_histogram(bins = 100) + 
  ylab("Density") + xlab("Distance R2_best - R2_2nd_best") +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  facet_grid( ~ clone_id, space="free") + 
  xlim(0,0.5) + 
  scale_fill_brewer(palette="Set1") +
  # geom_vline(xintercept = 0.05) +
  geom_vline(xintercept = 0.025) +
  # geom_vline(xintercept = 0.0125) +
  theme(strip.text.x = element_text(face="bold", size=12, colour = "black",),
        strip.text.y = element_text(face="bold", size=12, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        legend.title = element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"),
        legend.position = "none")
ggsave(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_selected_chromosomeArms_uncertainty_density_post_filtering.pdf"), height = 3, width = ncol(subset.t.clone.gene.profile)*3, dpi = 600)


# create a gene profile vector for the scRNA data
t.scRNA.clone.gene.profile <- create.scRNA.clone.profiles(subset.scRNA.mtx, complete.correlation.dataframe)

# create pseudobulk plot for the scRNA-seq profile
scRNA.pseudobulk <- plot.pseudobulk.profile(t(t.scRNA.clone.gene.profile), gene.ordering.file)
scRNA.pseudobulk$pseudobulk_plot
ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_correlation_selected_chromosomeArms_scRNA_geneProfile.pdf"), height = 10, width = 18, dpi = 600)

scRNA.plot.gene.profile <- scRNA.pseudobulk$plot_gene_profile
scRNA.plot.gene.profile$type <- "scRNA"
scDNA.plot.gene.profile <- scDNA.pseudobulk.subset$plot_gene_profile
scDNA.plot.gene.profile$type <- "scDNA"

# bind dna and rna together
multi.plot.gene.profile <- rbind(scDNA.plot.gene.profile, scRNA.plot.gene.profile)
multi.plot.gene.profile <- multi.plot.gene.profile[multi.plot.gene.profile$clone_id %in% valid.clones,]
multi.plot.gene.profile[multi.plot.gene.profile$value > 3, "value"] <- 4

# 
multi.plot.gene.profile[multi.plot.gene.profile$type == "scRNA", "value"] <- round(multi.plot.gene.profile[multi.plot.gene.profile$type == "scRNA", "value"])

clone <- valid.clones[1]
plotlist <- list()
for (clone in valid.clones){
  
  multi.plot.gene.profile.tmp <- multi.plot.gene.profile[multi.plot.gene.profile$clone_id == clone,]
  
  # create a pseudobulk plot
  p <- ggplot(multi.plot.gene.profile.tmp, aes(x=hgnc_symbol, y=value, group=type)) +
    geom_point(aes(color=type), position = position_jitter(w = 0, h = 0.1)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme_bw() + 
    scale_color_brewer(palette="Dark2") +
    facet_grid(type ~ reorder(chr), scales="free_x", space="free", switch="y") + 
    scale_y_continuous("Copy number", position="right", breaks = c(0,2,4), limits = c(0,4)) +   # Put the y-axis labels on the right
    ggExtra::removeGrid() + labs(x="Chromosomes", y="Copy number", color = "Type") +
    theme(panel.grid.major = element_blank(),
          strip.text.x = element_text(face="bold", size=12, colour = "black",),
          strip.text.y = element_text(face="bold", size=14, colour = "black",),
          strip.background = element_rect(fill="white", colour="black", size=1), 
          axis.text = element_text(colour = "black", size = 14, face = "bold" ),
          axis.title = element_text(colour = "black", size = 0, face = "bold" ),
          plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
          legend.title = element_text(color = "black", size = 16, face = "bold",),
          legend.text = element_text(colour="black", size=14, face="bold"),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3),
          legend.position = "none")
  
  plotlist[[clone]] <- p
  
}


switch(sample.tmp, 
                  "STP-PDX_clones"= ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_scRNA_correlation_selected_chromosomeArms_geneProfile_wocells.pdf"), arrangeGrob(grobs = plotlist, ncol = 1), height = 10, width = 12, dpi = 600),
                  "STP-Nuclei_clones" = ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_scRNA_correlation_15percent_chromosomeArms_geneProfile_wocells.pdf"), arrangeGrob(grobs = plotlist, ncol = 1), height = 10, width = 18, dpi = 600),
                  "MB243-Nuclei_clones" = ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_scRNA_correlation_15percent_chromosomeArms_geneProfile_wocells.pdf"), arrangeGrob(grobs = plotlist, ncol = 1), height = 8, width = 20, dpi = 600),
                  "ST1R-PDX_clones" = ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_scRNA_correlation_15percent_chromosomeArms_geneProfile_wocells.pdf"), arrangeGrob(grobs = plotlist, ncol = 1), height = 4, width = 5, dpi = 600))

# MB243-Nuclei
# ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_scRNA_correlation_15percent_chromosomeArms_geneProfile_wocells.pdf"), arrangeGrob(grobs = plotlist, ncol = 1), height = 8, width = 20, dpi = 600)

# STP-Nuclei
# ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_scRNA_correlation_15percent_chromosomeArms_geneProfile_wocells.pdf"), arrangeGrob(grobs = plotlist, ncol = 1), height = 10, width = 18, dpi = 600)

# STP-PDX
# ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_scRNA_correlation_selected_chromosomeArms_geneProfile_wocells.pdf"), arrangeGrob(grobs = plotlist, ncol = 1), height = 10, width = 12, dpi = 600)

# ST1R-PDX
# ggsave(paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_scRNA_correlation_15percent_chromosomeArms_geneProfile_wocells.pdf"), arrangeGrob(grobs = plotlist, ncol = 1), height = 4, width = 5, dpi = 600)

############################################################################
##        USE THE NEWLY GENERATED CLUSTERING TO GENERATE A HEATMAP
############################################################################

# read in gene_ordering file
gene.ordering.file <- read.table(gene.order.file, header = F, col.names = c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id"),sep = "\t")
gene.ordering.file$coordinates <- paste(paste0("chr", gene.ordering.file$chr), gene.ordering.file$start, gene.ordering.file$end, sep = ":")
rownames(gene.ordering.file) <- gene.ordering.file$hgnc_symbol

# read in the raw matrix
norm.gExp.matrix = as.matrix(read.table(scRNA.file, header = T, sep = " "))
colnames(norm.gExp.matrix) <- paste0(str_split_fixed(colnames(norm.gExp.matrix), "_pos", 2)[,1], "-1")

# visualize
plot_geneExp_heatmap(norm.gExp.matrix, 
                     gene.ordering.file, 
                     complete.correlation.dataframe, 
                     scDNA.clones <- unique(scDNA.plot.gene.profile$clone_id),
                     chroms <- unique(str_split_fixed(variable.chroms, "p|q", 2)[,1]),
                     distance = 0.025,
                     method = "correlation",
                     cutoff = "selected_chromosomeArms_uncertainty_0.025",
                     sample.tmp, output_directory = out.dir)

# subset to cells that do not pass uncertainty filter
subset.correlation.dataframe <- complete.correlation.dataframe[complete.correlation.dataframe$distance <= 0.025,]

############################################################################
##          VISUALISE THE UNCERTAINTY CROSS OVER BETWEEN CLONES
############################################################################

# create a matrix
uncert.mtx <- as.data.frame.matrix(prop.table(table(subset.correlation.dataframe$clone_id, subset.correlation.dataframe$clone_2nd_best_match),1))

# VISUALISE PROPORTION OF MUTATED CELLS PER GENE PER SAMPLE
library(circlize)
col_fun = colorRamp2(c(0, 1), c("white", "#B63679FF"))

ht1 <- Heatmap(uncert.mtx, name = "p-value", col = col_fun,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(abs(abs(uncert.mtx[i, j]) > 0.5))
                   grid.text(sprintf("%.2f", uncert.mtx[i, j]), x, y, gp = gpar(fontsize = 6, fontface = "bold", col = "black"))
               },
               cluster_rows = F, 
               cluster_columns = F,
               show_row_names = TRUE, 
               row_names_side = "left",
               show_row_dend = FALSE,
               show_column_names = TRUE,
               column_title = "2nd best match",
               row_title = "Best match",
               column_names_rot = 45,
               row_names_gp = gpar(fontsize =6, fontface = "bold"), 
               column_names_gp = gpar(fontsize = 6, fontface = "bold"), 
               rect_gp = gpar(col = "white", lwd = 1),
               # border_gp = gpar(col = "black", lwd = 2),
               column_title_gp = gpar(fontsize = 15, fontface = "bold"),
               row_title_gp = gpar(fontsize = 15, fontface = "bold"),
               heatmap_legend_param = list(color_bar = "continous",
                                           at = c(0, 1),
                                           title = "Proportion"), 
               border = T)

pdf(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_selected_chromosomeArms_uncertainty_2nd_best_matching_clone.pdf"), width=4, height = 3, pointsize=0.1)
print(ht1)
dev.off()


## The threshold to merge clones we use is 50% of cells which do not pass this distance threshold. 

complete.correlation.dataframe$merged_clone_id <- complete.correlation.dataframe$clone_id


switch(sample.tmp, 
       "STP-PDX_clones"= complete.correlation.dataframe[complete.correlation.dataframe$clone_id == "Clone2", "merged_clone_id"] <- "Clone1",
       "STP-Nuclei_clones" = {
         complete.correlation.dataframe[complete.correlation.dataframe$clone_id == "Clone2", "merged_clone_id"] <- "Clone1"
         complete.correlation.dataframe[complete.correlation.dataframe$clone_id == "Clone4", "merged_clone_id"] <- "Clone3"
       },
       "MB243-Nuclei_clones" = complete.correlation.dataframe[complete.correlation.dataframe$clone_id == "Clone3", "merged_clone_id"] <- "Clone2",
       "ST1R-PDX_clones" = complete.correlation.dataframe[complete.correlation.dataframe$clone_id == "Clone2", "merged_clone_id"] <- "Clone1")



# MB243-Nuclei
# complete.correlation.dataframe[complete.correlation.dataframe$clone_id == "Clone3", "merged_clone_id"] <- "Clone2"

# STP-PDX
# complete.correlation.dataframe[complete.correlation.dataframe$clone_id == "Clone2", "merged_clone_id"] <- "Clone1"

# STP-Nuclei
# complete.correlation.dataframe[complete.correlation.dataframe$clone_id == "Clone2", "merged_clone_id"] <- "Clone1"
# complete.correlation.dataframe[complete.correlation.dataframe$clone_id == "Clone4", "merged_clone_id"] <- "Clone3"

# ST1R-PDX
# complete.correlation.dataframe[complete.correlation.dataframe$clone_id == "Clone2", "merged_clone_id"] <- "Clone1"

write.table(complete.correlation.dataframe, paste0(out.dir,"/", sample.tmp, "/", sample.tmp, "_scDNA_clones_filtered_cells_DEG.txt"), col.names = T, row.names = F, quote = F, sep = "\t")


complete.correlation.dataframe <- read.table(paste0(out.dir, sample.tmp, "/", sample.tmp, "_scDNA_clones_filtered_cells_DEG.txt"), header = T, sep = "\t")
# complete.correlation.dataframe$clone_id <- complete.correlation.dataframe$merged_clone_id

# visualize
plot_geneExp_heatmap(norm.gExp.matrix, 
                     gene.ordering.file, 
                     complete.correlation.dataframe, 
                     scDNA.clones <- unique(complete.correlation.dataframe$merged_clone_id),
                     chroms <- unique(str_split_fixed(variable.chroms, "p|q", 2)[,1]),
                     distance = 0.025,
                     method = "correlation",
                     cutoff = "selected_chromosomeArms_uncertainty_0.025_merged",
                     sample.tmp, output_directory = out.dir)


############################################################################
##                    VISUALISE BARPLOTS FOR CLONES PRE AND POST MERGE
############################################################################

# get the proportions of each cell type in the respective condition
prop.df <- prop.table(table(complete.correlation.dataframe$clone_id))
prop.df.melt <- reshape2::melt(prop.df)
prop.df.melt$value <- round(prop.df.melt$value,4)
prop.df.melt$type <- "pre"

prop.post.df <- prop.table(table(complete.correlation.dataframe$merged_clone_id))
prop.post.df.melt <- reshape2::melt(prop.post.df)
prop.post.df.melt$value <- round(prop.post.df.melt$value,4)
prop.post.df.melt$type <- "post"

prop.df.melt <- rbind(prop.df.melt, prop.post.df.melt)
prop.df.melt$Var1 <- factor(prop.df.melt$Var1, levels = c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5"))
prop.df.melt$type <- factor(prop.df.melt$type, levels = c("pre", "post"))

fills <- c(brewer.pal(n=length(c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5")),"Set1"))
names(fills) <- c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5")

# Small multiple
ggplot(prop.df.melt, aes(fill = Var1, y=value, x=type)) + 
  geom_bar(position="stack", stat="identity",  color="black") +
  xlab("") + theme_classic() + ylab("Cell proportion") +
  scale_fill_manual(values = scales::alpha(c(fills, "red"),1))+
  theme(axis.text = element_text(colour = "black", size = 8, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 8, face = "bold", angle = 45, hjust = 1),
        axis.title = element_text(colour = "black", size = 10, face = "bold" ),
        plot.title = element_text(colour = "black", size = 10, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 8, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"), 
        legend.position = "bottom")
ggsave(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_pre_post_merging.pdf"), width = 3, height = 4)







###############################################################################
##            VISUALISE BARPLOTS FOR CLONES PRE AND POST MERGE AND POST FILTER
###############################################################################




i <- 1
complete.correlation.dataframe$cell_3rd_best_match <- NA
complete.correlation.dataframe$clone_3rd_best_match <- NA
complete.correlation.dataframe$distance_3rd_best_match <- NA
for (i in 1:nrow(complete.correlation.dataframe)){
  
  cell.data.tmp <- complete.correlation.dataframe[i,]
  cell_best_matched <- cell.data.tmp$pearson.correlation
  cell_3rd_best_match <- switch(sample.tmp, 
                                "STP-PDX_clones"= reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3", "Clone4", "Clone5")]),
                                "STP-Nuclei_clones" = cell_2nd_best_match <- reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3", "Clone4")]),
                                "MB243-Nuclei_clones" = cell_2nd_best_match <- reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3", "Clone4")]),
                                "ST1R-PDX_clones" = reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3")]))
  # cell_3rd_best_match <- reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3", "Clone4")])
  # cell_3rd_best_match <- reshape2::melt(cell.data.tmp[,c("Cell_barcode", "Clone1", "Clone2", "Clone3")])
  clone_3rd_best_match <- as.character(cell_3rd_best_match[order(cell_3rd_best_match$value, decreasing = T), "variable"])[3]
  cell_3rd_best_match_corr <- as.numeric(cell_3rd_best_match[order(cell_3rd_best_match$value, decreasing = T), "value"])[3]
  distance_3rd_best_match <- as.numeric(cell_3rd_best_match[order(cell_3rd_best_match$value, decreasing = T), "value"])[1] - as.numeric(cell_3rd_best_match[order(cell_3rd_best_match$value, decreasing = T), "value"])[3]
  
  complete.correlation.dataframe[i, "cell_3rd_best_match"] <- cell_3rd_best_match_corr
  complete.correlation.dataframe[i, "clone_3rd_best_match"] <- clone_3rd_best_match
  complete.correlation.dataframe[i, "distance_3rd_best_match"] <- distance_3rd_best_match
  
}



merged_clone_tables <- table(complete.correlation.dataframe$merged_clone_id, complete.correlation.dataframe$clone_id)

merged_clone_names <- names(which(rowSums(merged_clone_tables>0)>1))
complete.correlation.dataframe <- data.table(complete.correlation.dataframe)
complete.correlation.dataframe[,dist_filter := distance]

for(merged_clone in merged_clone_names){
  collapsing_clones <- unique(complete.correlation.dataframe[complete.correlation.dataframe$merged_clone_id==merged_clone,clone_id])
  complete.correlation.dataframe[clone_id %in% collapsing_clones & clone_2nd_best_match %in% collapsing_clones, dist_filter := distance_3rd_best_match]
}
complete.correlation.dataframe[,keep.cell := dist_filter >= 0.025]

if(length(unique(complete.correlation.dataframe[keep.cell==TRUE, merged_clone_id])) != length(unique(complete.correlation.dataframe[, merged_clone_id]))){
  stop("Sample lost a full clone after filtering!")
}
# complete.correlation.dataframe <- complete.correlation.dataframe[keep.cell == TRUE]

# get the proportions of each cell type in the respective condition
prop.df <- prop.table(table(complete.correlation.dataframe$clone_id))
prop.df.melt <- reshape2::melt(prop.df)
prop.df.melt$value <- round(prop.df.melt$value,4)
prop.df.melt$type <- paste0("pre-filter (n=", length(complete.correlation.dataframe$clone_id), ")")


write.table(complete.correlation.dataframe, paste0(out.dir,"/", sample.tmp, "/", sample.tmp, "_scDNA_clones_filtered_cells_3rd_added.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

# 
# prop.post.df <- prop.table(table(complete.correlation.dataframe$merged_clone_id))
# prop.post.df.melt <- reshape2::melt(prop.post.df)
# prop.post.df.melt$value <- round(prop.post.df.melt$value,4)
# prop.post.df.melt$type <- paste0("post (n=", length(complete.correlation.dataframe$merged_clone_id), ")")
# 
# 

prop.filtered.df <- prop.table(table(complete.correlation.dataframe[keep.cell==TRUE,merged_clone_id]))
prop.filtered.df.melt <- reshape2::melt(prop.filtered.df)
prop.filtered.df.melt$value <- round(prop.filtered.df.melt$value,4)
prop.filtered.df.melt$type <- paste0("post-filter (n=", nrow(complete.correlation.dataframe[keep.cell==TRUE]), ")")


prop.df.melt <- rbind(prop.df.melt,  prop.filtered.df.melt)
prop.df.melt$Var1 <- factor(prop.df.melt$Var1, levels = c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5"))
prop.df.melt$type <- factor(prop.df.melt$type, levels = rev(sort(unique(prop.df.melt$type))))

fills <- c(brewer.pal(n=length(c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5")),"Set1"))
names(fills) <- c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5")
colnames(prop.df.melt) <- c("Clone", "value", "type")
# Small multiple
ggplot(prop.df.melt, aes(fill = Clone, y=value, x=type)) + 
  geom_bar(position="stack", stat="identity",  color="black") +
  xlab("") + theme_classic() + ylab("Cell proportion") +
  scale_fill_manual(values = scales::alpha(c(fills, "red"),1))+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 10, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"), 
        legend.position = "left")



ggsave(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_pre_post_merging_post_filter.pdf"), height = 3, width = 5)



#################################################################
## Visualize the clone merging filter
###############################################################
  


p <- ggplot(complete.correlation.dataframe, aes(x=dist_filter, color=merged_clone_id)) + 
  stat_ecdf(geom = "step")  + scale_y_reverse() + 
  # geom_histogram(bins = 100) +
  ylab("Fraction") + xlab("∆ correlation top two clones") +
  theme_classic() +
  scale_color_manual(values=fills) +
  scale_fill_manual(values=fills) +
  # facet_grid(.~clone_id) + 
  geom_vline(xintercept = 0.025, col="grey30") +
  theme(strip.text.x = element_text(face="bold", size=12, colour = "black",),
        strip.text.y = element_text(face="bold", size=12, colour = "black",),
        # strip.background = element_rect(fill="white", colour="black", size=1), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        legend.title = element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"),
        legend.position = "none")
# ggsave(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_selected_chromosomeArms_dist_ecdf.pdf"), 
#       height = 3, width = 4, dpi = 600, device = )
grDevices::cairo_pdf(paste0(out.dir, sample.tmp,  "/", sample.tmp, "_scDNA_clones_selected_chromosomeArms_dist_ecdf.pdf"), 
                     height = 3, width = 4)
print(p)
dev.off()

###############################################################################
##            VISUALISE HEATMAP FOR CLONES  POST FILTER
###############################################################################

scDNA.clones <- unique(complete.correlation.dataframe$clone_id)
complete.correlation.dataframe <- data.frame(complete.correlation.dataframe)


complete.correlation.dataframe$clone_id <- complete.correlation.dataframe$merged_clone_id
# visualize
plot_geneExp_heatmap(norm.gExp.matrix, 
                     gene.ordering.file, 
                     complete.correlation.dataframe, 
                     scDNA.clones,
                     chroms <- unique(str_split_fixed(variable.chroms, "p|q", 2)[,1]),
                     distance = 0.025,
                     method = "correlation",
                     cutoff = "selected_chromosomeArms_uncertainty_0.025_merged_3rd",
                     sample.tmp, output_directory = out.dir)


# myx <- complete.correlation.dataframe$distance > 0.025
# 
complete.correlation.dataframe <- complete.correlation.dataframe[complete.correlation.dataframe$keep.cell==TRUE,]


# visualize
plot_geneExp_heatmap(norm.gExp.matrix, 
                     gene.ordering.file, 
                     complete.correlation.dataframe, 
                     scDNA.clones,
                     chroms <- unique(str_split_fixed(variable.chroms, "p|q", 2)[,1]),
                     distance = 0.025,
                     method = "correlation",
                     cutoff = "selected_chromosomeArms_uncertainty_0.025_merged_3rd_filtered",
                     sample.tmp, output_directory = out.dir)




