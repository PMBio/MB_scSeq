# SUPLLEMENTARY MATERIAL

## Aldinger et al. analysis

**scanpy_AldingerEtal2020_analysis.ipynb** - Python script to analyze and process the single cell RNA sequencing data from published data by Aldinger et al., 2020. 

**correlation_analysis_AldingerEtal2020.R** - R script to perform the correlation analysis of normal cell types and malignant primary tumor cells in order to assess the cell of origin.

## inferCNV analysis

**infercnv_prepare_input_XX.R** - R script to prepare the inputs to run inferCNV on each individual sample i.e. **LFSMBP-Nuclei** data.

**infercnv_analysis_ref.R** - R script to analyse the primary tumor and PDX data with inferCNV.

**bsub_XX_infercnv.sh** - Lsf submission script to run infercnv on each individual sample i.e. **LFSMBP-Nuclei** data.

## Individual sample snRNA- and scRNA-sequencing analysis

**scanpy_XX_analysis.ipynb** - Python script to analyze and process the single cell RNA sequencing data from each individual sample.

**XX_GSEA_clone_analysis.R** - R script use the DEGs from a between clone analysis for Gene Set Enrichment Analysis.

**visualize_celltype_clone_heatmap.R** - R script to assess the relationship between the cell types observed in each sample and the respective integrated clones from the procedure described in Figure 4.
