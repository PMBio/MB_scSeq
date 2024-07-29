library(ComplexHeatmap)
library(data.table)
library(QDNAseq)
library(RColorBrewer)
library(HMMcopy)
library(doParallel)
library(gtools)
library(ggplot2)
library(GenomicRanges)


source("~/scAbsolute-main/R/visualization.R")
source("~/scAbsolute-main/R/mean-variance.R")
source("~/scAbsolute-main/R/scAbsolute.R")
source("~/scAbsolute-main/R/core.R")



source("~/MB_SCSEQ/Figure01/CloneInference/commonSegmentation.R")
source("~/MB_SCSEQ/Figure01/CloneInference/plotFunctionsCNV.R")
source("~/MB_SCSEQ/Figure01/CloneInference/inferClonalCNVProfiles.R")
# source("inferClonalCNV_continous.r")
source("~/MB_SCSEQ/Figure01/CloneInference/calculateBootstrappedStatistics.R")

MB243 <- readRDS("~/groups/OE0585/internal/p163v/scAbsoluteWorkflow/results/scale/MB243-Nuclei/500/predict/out.rds")
protocolData(MB243) <- AnnotatedDataFrame()
colnames(MB243) <- gsub(colnames(MB243), pat = "\\-[0-9]\\.filtered\\.sorted$", rep = "")


var.filter <- mad(MB243$observed.variance/MB243$expected.variance)*2+median(MB243$observed.variance/MB243$expected.variance) > MB243$observed.variance/MB243$expected.variance
gini.filter <- mad(MB243$gini_normalized)*2+median(MB243$gini_normalized) > MB243$gini_normalized
# error.filter <- mad(MB243$error)*2+median(MB243$error) > MB243$error
## Ploidy 0 doesn't work for ploidy normalization
ploidy.filter <- MB243$ploidy.mod > 0
MB243 <- MB243[,var.filter&gini.filter&ploidy.filter]


set.seed(420)# set the seed of random number generator 

list.of.packages <- c("HMMcopy", "doParallel", "gtools", "ggplot2")


# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))


window_length <- 5e7
min_cnv_changes <- 10

nBootSamples <- 50
min_consec_cnvs <- 5

nthread <- 8
registerDoParallel(nthread)

sample = 'MB243-Nuclei'



## bin size param
bin_kb = 20
bin = bin_kb*1000





# get the readcounts for the bin size of interest
path_cells <- paste0('/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/HMMcopy_10x_scDNA/all_samples/', 
                     sample,'/readCount_filtered_bam/bin_', bin_kb,'kb/')
path_to_readcounts <- paste0('/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/HMMcopy_10x_scDNA/all_samples/', 
                             sample,'/readCount_filtered_bam/bin_', bin_kb,'kb/')
files <- list.files(path_to_readcounts, pattern = '.seg', full.names = F)

# set the paths to the reference files
gfile = paste0('/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/HMMcopy_10x_scDNA/refdata/bin_', bin_kb, 'kb/hg19.', bin_kb, 'kb.gc.seg')
mfile =  paste0('/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/HMMcopy_10x_scDNA/refdata/bin_', bin_kb, 'kb/hg19.', bin_kb, 'kb.map.seg')



## Make list of cells
cellList <- list()

names(files) <- sapply(strsplit(files, split = "_|-"), `[[`, 3)

cellList <- foreach(cell = colnames(MB243)) %dopar% {
  print(paste0("Processing cell: ", cell))
  
  
  rfile <- paste(path_cells,files[[cell]],sep = "")
  cellData <- wigsToRangedData(rfile, gfile, mfile)
  
  (cellData)

}
names(cellList) <- colnames(MB243)

MB243$ploidy.mod[(MB243$ploidy > 3.4 & MB243$ploidy < 4)] <- 4
clonal.cnv.res.boot <- inferClonalCNVs(MB243, cellList, bootstrap = TRUE, nBoot = 101)
# saveRDS(clonal.cnv.res.boot, file="mb243_boot_res.rds")
saveRDS(clonal.cnv.res.boot,"mb243_boot_res_ploidyfixed.rds")
