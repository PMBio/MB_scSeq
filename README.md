# Title: "Multi-omic and single-cell profiling of chromothriptic medulloblastoma reveals genomic and transcriptomic consequences of genome instability"

## Authors: Petr Smirnov, Moritz J. Przybilla, Gonzalo Parra, & Hana Susak

## Date: 01/08/2024

This repository was created to store and provide all scripts used in the manuscript "Multi-omic and single-cell profiling of chromothriptic medulloblastoma reveals genomic and transcriptomic consequences of genome instability". This includes the scripts for generation of figures and for the analyses implemented. In line with the structure of the manuscript, the scripts are split according each individual Figures. 

For an example environment set up to run the related code, as well as an example of one analyzed sample for scDNAseq clone identification and CT calling, please see the asssociated Code Ocean capsule here: https://codeocean.com/capsule/5086193/tree 

To reproduce analyses for other samples, you can download the related processed data from this Zenodo repository: 10.5281/zenodo.13348419 

**Figure01** - Genetic heterogeneity in medulloblastomas with chromothripsis. 

**Figure02** - Chromothripsis is a major event for the formation of double-minute chromosomes.

**Figure03** - Distinct transcriptional programs dominate malignant cell types in medulloblastoma with chromothripsis & Integrating DNA and RNA on single cell level reveals copy number related pathway alterations.

**Figure04** - Combining single-nuclei DNA-seq and bulk whole-genome sequencing identifies early events potentially facilitating chromothripsis.

**Figure05** - Inactivation of SETD2 and TP53 in neural stem cells leads to genome instability.

**Supplementary_Material** - Additional Figures provided in the Supplementary Material of the manuscript.

Scripts for running the cellranger pipelines for [scRNA](https://support.10xgenomics.com/single-cell-gene-expression) and [scDNA](https://support.10xgenomics.com/single-cell-dna) were not included here, as they are simply mirroring single line bash commands from the 10X analysis pipeline.

## Reproducing locally

The code for reproducing the analysis of scDNAseq data was tested to run on Linux operating systems () under R 4.2.0. The recommended reproduction instruction is to clone the Code Ocean capsule above, or to download the docker image, and to download the processed HMMCopy and scAbsolute data from the Zenodo repository. 

The scRNAseq analysis should be reproducable using the provided scripts under python 3.12.5 using scanpy version 1.10.2. 
