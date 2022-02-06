# FIGURE 4

## G&T-seq 

**gtDNA_trimming.sh** & **gtRNA_trimming.sh** - Script to trim reads from G&T-seq data from both DNA and RNA.

**fastq_screen_scDNA.sh** & **fastq_screen_scRNA.sh** - Script to assess the percentage of mouse and human reads within G&T-seq data from both DNA and RNA.

**fastqc_scDNA.sh** & **fastqc_scRNA.sh** - Script to assess the quality of the reads within G&T-seq data from both DNA and RNA.

**picard_markduplicates_scDNA.sh** - Script to deduplicate the scDNA-seq reads. 

**bwa_alignment_hg19.sh** - Script to align the deduplicated scDNA reads to hg19. 

**STAR_align.sh** - Script to align the scRNA reads to hg19.

**HTseq.sh** - Script to count the scRNA reads from the STAR-aligned bam files.

## scDNA- & scRNA-seq copy number clone alignment

**infercnv_prepare_input_XX.R** - R script to prepare the inputs to run inferCNV on the **GTseq** and **nuclei** or single samples i.e. **MB243-Nuclei** data.

**infercnv_analysis_XX.R** - R script to analyse the **GTseq** and **nuclei** data with inferCNV.

**bsub_XX_infercnv.sh** - Lsf submission script to run infercnv on the **GTseq** and **nuclei** or single samples i.e. **MB243-Nuclei** data.

**make_GenexCell_matrix.R** - Transform the scDNA-seq data matrices from segment x cells to genes x cells matrices to align with the scRNA-seq data.

**scDNA_scRNA_alignment.R** - R script to align the scRNA-seq copy number data to the scDNA-seq copy number clones. 

**scDNA_scRNA_alignment_function.R** - Functions for the R script to align the scRNA-seq copy number data to the scDNA-seq copy number clones. 

**visualize_scDNA_scRNA_clone_correlation.R** - R script to visualize the results from the scDNA and scRNA copy number clone alignment.

**plot_scCNV_gene_heatmap.R** - R script use visualize the copy number heatmaps from the scDNA data in the same format as the scRNA-seq data.

## MB243-Nuclei analysis

**scanpy_MB243-Nuclei_analysis.ipynb** - Python script to analyze and process the single cell RNA sequencing data from MB243-Nuclei.

**MB243-Nuclei_GSEA_clone_analysis.R** - R script use the DEGs from a between clone analysis for Gene Set Enrichment Analysis.

**correlation_analysis_DM_DNA_RNA.R** - R script to assess the influence of double minute chromosomes on gene expression.
