#############################################################################################################################
##                                                                                                                      
##  PERFORM BULK RNA-SEQ GSEA ANALYSIS FOR SUPPLEMENTARY FIGURE 6
##
##  Date: 11 MAY 2021                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                            
##                                                                                                                      
#############################################################################################################################

# clear workspace
rm(list = ls())
set.seed(14) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("Matrix", "biomaRt", "viridis", "GO.db", "HTSanalyzeR2", "org.Hs.eg.db", "KEGGREST", "igraph", 
                      "tidyr", "stats", "reshape2", "ggplot2", "forcats", "heatmap3", "gplots", "ggpubr", "dplyr",
                      "ComplexHeatmap", "ggplotify", "TxDb.Hsapiens.UCSC.hg38.knownGene", "regioneR", "Repitools")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

###########################################################################
#                       READ IN THE DATA OF INTEREST
###########################################################################

# define output directory
o.dir <- "/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision"

# create o.dir
dir.create(o.dir)
setwd(o.dir)

# get sample.ids
sample.tmp <- "Bulk_RNA"

# which sample? 
message(sample.tmp)
dir.create(paste0(o.dir, "/DEG"))

# read in files of interest
marker.files <- fread("/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/MB_CT_vs_NCT_DESeq2_TCC_GROUP_wo_Immune.txt", header = T)
marker.files <- marker.files[marker.files$baseMean > median(marker.files$baseMean), ]
marker.files <- marker.files[,c("gene_name", "log2FoldChange", "padj")]
marker.files <- marker.files[!duplicated(marker.files$gene_name),]

# filter the list of genes by logFC >1 and FDR == 5%
up.marker.files <- marker.files %>% top_n(1000, log2FoldChange)
up.marker.files <- up.marker.files[order(up.marker.files$log2FoldChange, decreasing = T),]
down.marker.files <- marker.files %>% top_n(-1000, log2FoldChange)
down.marker.files <- down.marker.files[order(down.marker.files$log2FoldChange, decreasing = F),]

# combine 
marker.files <- rbind(up.marker.files, down.marker.files)

# change colnames
colnames(marker.files) <- c("hgnc_symbol", "avg_logFC", "p_val")
marker.files$p_val <- as.numeric(marker.files$p_val)
rownames(marker.files) <- marker.files$hgnc_symbol

###########################################################################
#                       PERFORM THE GSEA ANALYSIS
###########################################################################

# GO THROUGH THE ANALYSIS STEPS
dfile_subs <- marker.files
phenotype <- as.vector(as.numeric(dfile_subs$avg_logFC)) ### USE LOG FOLD CHANGE AS PHENOTYPE VECTOR
names(phenotype) <- rownames(dfile_subs)

## specify the gene sets type you want to analyze

# HALLMARK
MSig_H <- MSigDBGeneSets(species = "Hs", collection = "H", subcategory = NULL) # Hallmarks!
ListGSC <- list(MSig_H=MSig_H)

## iniate a *GSCA* object
gsca <- GSCA(listOfGeneSetCollections=ListGSC, 
             geneList=phenotype)

## preprocess
gsca1 <- preprocess(gsca, species="Hs", initialIDs="SYMBOL",
                    keepMultipleMappings=TRUE, duplicateRemoverMethod="max",
                    orderAbsValue=FALSE)


## analysis
if (requireNamespace("doParallel", quietly=TRUE)) {
  doParallel::registerDoParallel(cores=4)
}  ## support parallel calculation using multiple cores


gsca2 <- analyze(gsca1, 
                 para=list(pValueCutoff=0.05, pAdjustMethod="BH",
                           nPermutations=100, minGeneSetSize=1,
                           exponent=1), 
                 doGSOA = FALSE)


## append gene sets terms
gsca3 <- appendGSTerms(gsca2, msigdbGSCs=c("MSig_H"))

## draw GSEA plot for a specific gene set
topGS <- getTopGeneSets(gsca3, resultName="GSEA.results",
                        gscs=c("MSig_H"),allSig=TRUE)#, allSig=TRUE)

# export the matrix with the differentially expressed pathways 
results <- getResult(gsca3)$GSEA.results$MSig_H
write.table(results, paste0(o.dir, "/scrna_analysis/Hallmark_GSEA_bulkRNA.txt"),sep="\t", quote=FALSE, row.names = TRUE, col.names=TRUE)


## HALLMARK_COAGULATION
pdf(paste0("DEG/GSEA_EnrichmentPlot_COAGULATION.pdf"),width=5,height=5,pointsize=0.1)
viewGSEA(gsca3, gscName="MSig_H", gsName="HALLMARK_COAGULATION", main.title = "COAGULATION")
dev.off()

## HALLMARK_COMPLEMENT
pdf(paste0("DEG/GSEA_EnrichmentPlot_COMPLEMENT.pdf"),width=5,height=5,pointsize=0.1)
viewGSEA(gsca3, gscName="MSig_H", gsName="HALLMARK_COMPLEMENT", main.title = "COMPLEMENT")
dev.off()

## HALLMARK_MTORC1_SIGNALING
pdf(paste0("DEG/GSEA_EnrichmentPlot_MTOR.pdf"),width=5,height=5,pointsize=0.1)
viewGSEA(gsca3, gscName="MSig_H", gsName="HALLMARK_MTORC1_SIGNALING", main.title = "MTOR SIGNALING")
dev.off()

## HALLMARK_HEDGEHOG_SIGNALING
pdf(paste0("DEG/GSEA_EnrichmentPlot_HEDGEHOG.pdf"),width=5,height=5,pointsize=0.1)
viewGSEA(gsca3, gscName="MSig_H", gsName="HALLMARK_HEDGEHOG_SIGNALING", main.title = "HEDGEHOG SIGNALING")
dev.off()

## HALLMARK_KRAS_SIGNALING_UP
pdf(paste0("DEG/GSEA_EnrichmentPlot_KRAS.pdf"),width=5,height=5,pointsize=0.1)
viewGSEA(gsca3, gscName="MSig_H", gsName="HALLMARK_KRAS_SIGNALING_UP", main.title = "KRAS SIGNALING UP")
dev.off()

## HALLMARK_UNFOLDED_PROTEIN_RESPONSE
pdf(paste0("DEG/GSEA_EnrichmentPlot_PROTEIN.pdf"),width=5,height=5,pointsize=0.1)
viewGSEA(gsca3, gscName="MSig_H", gsName="HALLMARK_UNFOLDED_PROTEIN_RESPONSE", main.title = "UNFOLDED PROTEIN RESPONSE")
dev.off()

## HALLMARK_GLYCOLYSIS
pdf(paste0("DEG/GSEA_EnrichmentPlot_GLYCOLYSIS.pdf"),width=5,height=5,pointsize=0.1)
viewGSEA(gsca3, gscName="MSig_H", gsName="HALLMARK_GLYCOLYSIS", main.title = "GLYCOLYSIS")
dev.off()




