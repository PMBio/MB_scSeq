################################################################################################################################################
##                                                                                                                      
##  CREATE A COUNT MATRIX AND SEURAT OBJECT FROM THE SMARTSEQ2 DATA
##                                                                                                                      
##  Date: 05 DECEMBER 2021                                                                                                                     
##  
##  Author: Moritz Przybilla                                                                                                                    
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

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# READ IN THE DATA FROM HTSEQ
#####################################################################################

# read sample metadata
# detailed.metadata <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/23806_meta.tsv", header = T)
scdna.detailed.metadata <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/24197_meta.tsv", header = T)
scrna.detailed.metadata <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/24198_meta.tsv", header = T)
detailed.metadata <- as.data.frame(rbind(scdna.detailed.metadata, scrna.detailed.metadata))

smartseq.metadata <- fread("/omics/odcf/project/OE0585/chromothripsis/sequencing/sc_smartseq2_sequencing/view-by-pid/OE0014_chromothripsis_LFS_MB_P_PDX/tumor01-replicate02-x/0_all/OE0014_chromothripsis_LFS_MB_P_PDX_tumor01-replicate02-x_mapping.tsv", header = F)
smartseq.metadata$V3 <- basename(smartseq.metadata$V1)
sc.wgs.metadata <- fread("/omics/odcf/project/OE0585/chromothripsis/sequencing/sc_whole_genome_sequencing/view-by-pid/OE0014_chromothripsis_LFS_MB_P_PDX/tumor01-replicate02-x/0_all/OE0014_chromothripsis_LFS_MB_P_PDX_tumor01-replicate02-x_mapping.tsv", header = F)
sc.wgs.metadata$V3 <- basename(sc.wgs.metadata$V1)

# combine smartseq
smartseq.complete.metadata <- merge(smartseq.metadata, detailed.metadata, by.x = "V3", by.y = "FASTQ_FILE")
sc.wgs.complete.metadata <- merge(sc.wgs.metadata, detailed.metadata, by.x = "V3", by.y = "FASTQ_FILE")
complete.metadata <- as.data.frame(rbind(smartseq.complete.metadata, sc.wgs.complete.metadata))

# list files from scDNA and scRNA
htseq.files <- list.files("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_384/HTseq", pattern = "\\.htseq_counts.out", all.files = T, full.names = T)
htseq.data <- list()

# and the folder names
samples <- unique(str_split_fixed(basename(htseq.files), "\\.htseq", 2)[,1])

i <- 1
for (i in 1:length(samples)){
  
  sample.tmp <- samples[i]
  data.tmp <- read.table(htseq.files[grep(paste0(sample.tmp, ".htseq"), htseq.files)], header = F)
  colnames(data.tmp) <- c("ensembl_id", sample.tmp)  
  htseq.data[[i]] <- data.tmp
  
}

merged.htseq.data = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ensembl_id", all = TRUE), htseq.data)
merged.htseq.data <- merged.htseq.data[-c(1:5),]

# get mart object with certain attributes
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ann <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"), "ensembl_gene_id", merged.htseq.data$ensembl_id, mart  = mart)

# order cns_reads
gene_ordering_file <- ann
colnames(gene_ordering_file) <- c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id")
gene_ordering_file$chr <- paste0("chr", gene_ordering_file$chr)

# combine htseq 
complete.htseq.data <- merge(merged.htseq.data, gene_ordering_file, by.x = "ensembl_id", by.y = "ensembl_gene_id", all.x = T)

# set accepted levels and remove other chromosomes
chrOrder<-c(paste("chr",1:22,sep=""))
complete.htseq.data$chr <- factor(complete.htseq.data$chr, levels = chrOrder)
complete.htseq.data <- complete.htseq.data[complete.cases(complete.htseq.data),]

# order file and remove duplicated values
complete.htseq.data <- complete.htseq.data[with(complete.htseq.data, order(chr, start)),]
complete.htseq.data[, c("ensembl_id", "chr", "start", "end")] <- NULL
complete.htseq.data <- complete.htseq.data %>% group_by(hgnc_symbol) %>% summarise_all(.funs = sum,na.rm=T)
complete.htseq.data <- complete.htseq.data[-1,]

#
complete.htseq.data <- as.data.frame(complete.htseq.data)
rownames(complete.htseq.data) <- complete.htseq.data$hgnc_symbol
complete.htseq.data$hgnc_symbol <- NULL
write.table(complete.htseq.data, "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA/gtseq_rna_counts.txt", col.names = T, row.names = T,
            sep = "\t", quote = F)

complete.htseq.data <- read.table("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA/gtseq_rna_counts.txt", header = T)
counts <- as.data.frame(colSums(complete.htseq.data))
counts$Sample <- rownames(counts)
counts.complete <- merge(counts, complete.metadata, by.x = "Sample", by.y = "SAMPLE_NAME", all.x = T)
controls <- counts.complete[grep("A1$|A24$|P24$|E4$|I10$|L19$", counts.complete$V2), ]
controls$control <- "Negative"
controls[grep("A1|A24|P24", controls$V2), "control"] <- "Positive"

ggplot(controls, aes(as.numeric(controls[,2]), fill = control)) +
  geom_histogram(bins = 100, color = "black", size = 0.25) + 
  ylab("Frequency") + xlab("# Counts") +
  theme_classic() +
  scale_fill_manual(values = c("Darkred", "Darkgreen"))+
  theme(strip.text.x = element_text(face="bold", size=12, colour = "black",),
        strip.text.y = element_text(face="bold", size=12, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        legend.title = element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"),
        legend.position = "none") +
  scale_x_continuous(labels = scales::comma)
ggsave("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_library_size_control.pdf", width = 4, height = 3, dpi = 600)

mean_counts <- mean(counts$`colSums(complete.htseq.data)`)

# check the density distribution of correlation coefficients
ggplot(counts.complete, aes(as.numeric(colSums(complete.htseq.data)))) +
  geom_histogram(bins = 200, fill = "darkgrey", color = "black", size = 0.25) + 
  ylab("Frequency") + xlab("# Counts") +
  theme_classic() +
  geom_vline(aes(xintercept = mean_counts), colour = 'darkblue', size = 1) +
  theme(strip.text.x = element_text(face="bold", size=12, colour = "black",),
        strip.text.y = element_text(face="bold", size=12, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        legend.title = element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"),
        legend.position = "none") +
  scale_x_continuous(labels = scales::comma)
ggsave("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_matrix_librarySize.pdf", width = 5, height = 4)









