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
source("inferClonalCNV_continous.r")
source("~/MB_SCSEQ/Figure01/CloneInference/calculateBootstrappedStatistics.R")


RCMB18PDX <- readRDS("~/groups/OE0585/internal/p163v/scAbsoluteWorkflow/results/scale/RCMB18-PDX/500/predict/out.rds")
protocolData(RCMB18PDX) <- AnnotatedDataFrame()
colnames(RCMB18PDX) <- gsub(colnames(RCMB18PDX), pat = "\\-[0-9]\\.filtered\\.sorted$", rep = "")

var.filter <- mad(RCMB18PDX$observed.variance/RCMB18PDX$expected.variance)*2+median(RCMB18PDX$observed.variance/RCMB18PDX$expected.variance) > RCMB18PDX$observed.variance/RCMB18PDX$expected.variance
gini.filter <- mad(RCMB18PDX$gini_normalized)*2+median(RCMB18PDX$gini_normalized) > RCMB18PDX$gini_normalized
# error.filter <- mad(RCMB18PDX$error)*2+median(RCMB18PDX$error) > RCMB18PDX$error
## Ploidy 0 doesn't work for ploidy normalization
ploidy.filter <- RCMB18PDX$ploidy.mod > 0
RCMB18PDX <- RCMB18PDX[,var.filter&gini.filter&ploidy.filter]

set.seed(42)# set the seed of random number generator 

list.of.packages <- c("HMMcopy", "doParallel", "gtools", "ggplot2")


# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))


window_length <- 5e7
min_cnv_changes <- 10

nBootSamples <- 50
min_consec_cnvs <- 5

nthread <- 8
registerDoParallel(nthread)

## set sample of interest
sample = 'RCMB18-PDX'



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

cellList <- foreach(cell = colnames(RCMB18PDX)) %dopar% {
  print(paste0("Processing cell: ", cell))
  
  
  rfile <- paste(path_cells,files[[cell]],sep = "")
  cellData <- wigsToRangedData(rfile, gfile, mfile)
  
  (cellData)

}
names(cellList) <- colnames(RCMB18PDX)




clonal.cnv.res.boot <- inferClonalCNVsCont(RCMB18PDX, cellList, bootstrap = TRUE, nBoot=101)
saveRDS(clonal.cnv.res.boot, file="RCMB18_PDX_boot_res_cont.rds")
