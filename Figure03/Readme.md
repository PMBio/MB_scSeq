# FIGURE 3

## snRNA- & scRNA-seq analysis

**scanpy_Nuclei_analysis.ipynb** - Python script to analyze and process the single cell RNA sequencing data from the single-nuclei RNA-sequencing data from primary tumors.

**scanpy_PDX_analysis.ipynb** - Python script to analyze and process the single cell RNA sequencing data from the single-cell RNA-sequencing data from PDX models.

**create_barplot_Celltypes.R** - R script to create the barplot shown in Figure 3 for the cell type composition.

## Riemondy et al. projection

**scanpy_RiemondyEtal2021_analysis.ipynb** - Python script to analyze and process the single cell RNA sequencing data from published data by Riemondy et al., 2021. 

**create_barplot_ingest_RiemondyEtal2021_projection.R** - R script to create the barplot shown in Figure 3 for the MB subgroups.

## Bulk RNA-sequencing analysis

**bulkRNA_DESeq2_deconvolution** - R script for the analysis of bulk RNA-sequencing data from LFS and SHH medulloblastoma using DESeq2. A manual deconvolution was implemented using data from our snRNA-seq primary tumor data as well as Riemondy et al. data.

**bulkRNA_GSEA_analysis.R** - R script to run GSEA on the DEGs from the bulk RNA-sequencing data.

**XX.R** - R script to create violin plots from the upregulated and downregulated genes resulting from the bulk RNA DESeq2 analysis. The boxplots are shown in Figure3. 

