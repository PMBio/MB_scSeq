################################################################################################################################################
##                                                                                                                      
##  CREATE SAMPLESHEET FOR FASTQ FILES FROM RNA-SEQ TO RUN KALLISTO
##                                                                                                                      
##  Date: 07 OCTOBER 2021                                                                                                                     
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
##           
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("BiocManager", "reshape2", "stringr", "readr", "tidyverse", "ggpubr", "ComplexHeatmap", "BuenColors", 
                      "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# CREATE XENOME SAMPLE SHEET FOR FASTQ FILES FROM BAM2FASTQ
#####################################################################################

# read sample metadata
metadata <- fread("/Users/mp34/dkfz/przybilm/medulloblastoma/revision/gt_seq/23806_meta.tsv", header = T)

# this is the path to the fastq files
file.path <- "/omics/gpcf/midterm/023806/data/211122_NB551291_0452_AHKWLCBGXK"

# and the folder names
samples <- unique(str_split_fixed(metadata$FASTQ_FILE, "_R1", 2)[,1])

scDNA.fastq.list <- data.frame("Sample" = NA, "R1" = NA, stringsAsFactors = F)
scRNA.fastq.list <- data.frame("Sample" = NA, "R1" = NA, stringsAsFactors = F)
for(i in 1:length(samples)){
  
  sample.tmp <- samples[i]
  true.sample.id <- metadata[grep(sample.tmp, metadata$FASTQ_FILE), "SAMPLE_NAME"]
  fastq.file <- metadata[grep(sample.tmp, metadata$FASTQ_FILE), FASTQ_FILE]
  
  
  if (length(grep("DNA", true.sample.id)) == 1){
    
    fastq.df <- data.frame("Sample" = true.sample.id$SAMPLE_NAME, 
                           "R1" = paste0(file.path, "/", sample.tmp, "/fastq/", fastq.file), 
                           stringsAsFactors = F)
    
    scDNA.fastq.list <- rbind(scDNA.fastq.list, fastq.df)
    
  } else {
    
    fastq.df <- data.frame("Sample" = true.sample.id$SAMPLE_NAME, 
                           "R1" = paste0(file.path, "/", sample.tmp, "/fastq/", fastq.file), 
                           stringsAsFactors = F)
    
    scRNA.fastq.list <- rbind(scRNA.fastq.list, fastq.df)
    
  }
  
}

scDNA.fastq.list <- scDNA.fastq.list[-1,]
scRNA.fastq.list <- scRNA.fastq.list[-1,]
write.table(scDNA.fastq.list, "/Users/mp34/dkfz/przybilm/medulloblastoma/revision/gt_seq/scDNA/scDNA_trimgalore_samplesheet.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(scRNA.fastq.list, "/Users/mp34/dkfz/przybilm/medulloblastoma/revision/gt_seq/scRNA/scRNA_trimgalore_samplesheet.txt", sep = "\t", col.names = F, row.names = F, quote = F)

#####################################################################################
# CREATE BWA SAMPLE SHEET
#####################################################################################

# this is the path to the fastq files
fastq.files <- list.files("/Users/mp34/dkfz/przybilm/medulloblastoma/revision/gt_seq", pattern = "\\.fq.gz$", recursive = T, all.files = T, full.names = T)

# and the folder names
samples <- unique(str_split_fixed(fastq.files, "/", 12)[,11])

scDNA.fastq.list <- data.frame("Sample" = NA, "R1" = NA, stringsAsFactors = F)
scRNA.fastq.list <- data.frame("Sample" = NA, "R1" = NA, stringsAsFactors = F)
i <- 1
for(i in 1:length(samples)){
  
  sample.tmp <- samples[i]
  fastq.file <- fastq.files[grep(paste0(sample.tmp, "/"), fastq.files)]
  
  
  if (length(grep("DNA", sample.tmp)) == 1){
    
    fastq.df <- data.frame("Sample" = sample.tmp,
                           "R1" = fastq.file,
                           stringsAsFactors = F)
    
    scDNA.fastq.list <- rbind(scDNA.fastq.list, fastq.df)
    
  } else {
    
    fastq.df <- data.frame("Sample" = sample.tmp,
                           "R1" = fastq.file,
                           stringsAsFactors = F)
    
    scRNA.fastq.list <- rbind(scRNA.fastq.list, fastq.df)
    
  }
  
}

scDNA.fastq.list <- scDNA.fastq.list[-1,]
scDNA.fastq.list$R1 <- gsub("/Users/mp34/dkfz/", "/omics/groups/OE0540/internal/projects/", scDNA.fastq.list$R1)
scRNA.fastq.list <- scRNA.fastq.list[-1,]
scRNA.fastq.list$R1 <- gsub("/Users/mp34/dkfz/", "/omics/groups/OE0540/internal/projects/", scRNA.fastq.list$R1)
write.table(scDNA.fastq.list, "/Users/mp34/dkfz/przybilm/medulloblastoma/revision/gt_seq/scDNA/scDNA_bwa_samplesheet.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(scRNA.fastq.list, "/Users/mp34/dkfz/przybilm/medulloblastoma/revision/gt_seq/scRNA/scRNA_star_samplesheet.txt", sep = "\t", col.names = F, row.names = F, quote = F)


