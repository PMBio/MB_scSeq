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

ST1R <- readRDS("~/groups/OE0585/internal/p163v/scAbsoluteWorkflow/results/scale/STPR-Nuclei/500/predict/out.rds")
protocolData(ST1R) <- AnnotatedDataFrame()
colnames(ST1R) <- gsub(colnames(ST1R), pat = "\\-[0-9]\\.filtered\\.sorted$", rep = "")



var.filter <- mad(ST1R$observed.variance/ST1R$expected.variance)*2+median(ST1R$observed.variance/ST1R$expected.variance) > ST1R$observed.variance/ST1R$expected.variance
gini.filter <- mad(ST1R$gini_normalized)*2+median(ST1R$gini_normalized) > ST1R$gini_normalized
# error.filter <- mad(ST1R$error)*2+median(ST1R$error) > ST1R$error
## Ploidy 0 doesn't work for ploidy normalization
ploidy.filter <- ST1R$ploidy.mod > 0
ST1R <- ST1R[,var.filter&gini.filter&ploidy.filter]


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
sample = 'ST1R-Nuclei'



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

cellList <- foreach(cell = colnames(ST1R)) %dopar% {
  print(paste0("Processing cell: ", cell))
  
  
  rfile <- paste(path_cells,files[[cell]],sep = "")
  cellData <- wigsToRangedData(rfile, gfile, mfile)
  
  (cellData)

}
names(cellList) <- colnames(ST1R)



clonal.cnv.res.boot <- inferClonalCNVsCont(ST1R, cellList, bootstrap = TRUE, nBoot=101)
saveRDS(clonal.cnv.res.boot, file="ST1R_Nuc_boot_res_cont.rds")
