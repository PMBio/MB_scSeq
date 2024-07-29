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
source("~/MB_SCSEQ/Figure01/CloneInference/calculateBootstrappedStatistics.R")


list.of.folders <- list.files("~/groups/OE0585/internal/p163v/scAbsoluteWorkflow/results/scale/", pattern = "pell",full.names = T)

Clone2a <- readRDS("/home/p163v/groups/OE0585/internal/p163v/scAbsoluteWorkflow/results/scale//pellmann_rpe1_BR_CL_170825_S2_F9_CNV/500/predict/out.rds")

protocolData(Clone2a) <- AnnotatedDataFrame()

var.filter <- mad(Clone2a$observed.variance / Clone2a$expected.variance) * 2 + median(Clone2a$observed.variance / Clone2a$expected.variance) > Clone2a$observed.variance / Clone2a$expected.variance
gini.filter <- mad(Clone2a$gini_normalized) * 2 + median(Clone2a$gini_normalized) > Clone2a$gini_normalized
# error.filter <- mad(Clone2a$error)*2+median(Clone2a$error) > Clone2a$error
## Ploidy 0 doesn't work for ploidy normalization
ploidy.filter <- Clone2a$ploidy.mod > 0
Clone2a <- Clone2a[, var.filter & gini.filter & ploidy.filter]

set.seed(42)# set the seed of random number generator 

list.of.packages <- c("HMMcopy", "doParallel", "gtools", "ggplot2")


# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))


window_length <- 5e7
min_cnv_changes <- 10

nBootSamples <- 50
min_consec_cnvs <- 5

nthread <- 10
registerDoParallel(nthread)


## set sample of interest
sample = 'BR_CL_170825_S2_F9_CNV'



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

names(files) <- sapply(strsplit(files, split = "_|\\."), `[[`, 5)
# colnames(Clone2a) <- sapply(strsplit(colnames(Clone2a), split = "_|-"), `[[`, 2)


cellList <- foreach(cell = colnames(Clone2a)) %dopar% {
  print(paste0("Processing cell: ", cell))
  
  
  rfile <- paste(path_cells,files[[cell]],sep = "")
  cellData <- wigsToRangedData(rfile, gfile, mfile)
  
  # correctReadcount2(cellData)
  
}
names(cellList) <- colnames(Clone2a)



library(doParallel)
registerDoParallel(10)

clonal.cnv.res.boot <- inferClonalCNVs(Clone2a, cellList, bootstrap = TRUE, nBoot=101)
saveRDS(clonal.cnv.res.boot, file="BR_CL_170825_S2_F9_CNV_boot_res.rds")


