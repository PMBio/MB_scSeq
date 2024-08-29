# FIGURE 1

## scDNAseq Data Processing

Prior to running the code in this folder, the 10x CNVKit sequencing data from all samples was processed by Cell Ranger v1.1.0, followed by `cellranger-dna bamslice` to create a BAM file for each detected cell (using per_cell_summary_metrics.csv to determine cell barcodes). Samtools was used to filter alignements based on the following flags:

read unmapped (0x4), not primary alignment (0x100), read fails platform/vendor quality checks (0x200), read is PCR or optical duplicate (0x400), as well as supplementary alignment (0x800) and a mapping quality ≤ 30. 

The filtered files were then used as input to `hmmcopy-utils readCounter` to create readCount files at 20kb resolution, as well as scAbsolute, which was run at 500kb resolution unless otherwise stated in the manuscript. 

The version of scAbsolute run was prior to public release of the tool, using v2.9.5 of the docker image and commit 896e8cf316d977533c9803ed5d29c0ddfbb1c02b of the scDNAseq-workflow (https://github.com/markowetzlab/scDNAseq-workflow). Copy states from 1 to 16 copies were allowed.  

## Clonal Inference

As a prerequisite, clone the scAbsolute github (https://github.com/markowetzlab/scAbsolute) version commit version 1b2cf53, and adjust the source file location in the notebooks.

The clonal inference was done using the notebooks in the folder CloneInference using input from scAbsolute and HMMUtils, as described above. 

## clonal CT Scoring

Scoring of Clonal CT events was done using the code in the CT_scoring subfolder

## Figure 1
Figure 1 was plotted using Figure1.Rmd

## Supplementary Figures

Supplementary Figures for Figure 1 were plotted using the notebooks in the Supplementary Figures folder.    



