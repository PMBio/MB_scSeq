#!/bin/sh
# LSF directives 
#BSUB -J inferCNV_MB
#BSUB -e T-%J-%I.err                  # error file path, where % will be replaced with job ID
#BSUB -o T-%J-%I.out                  # outut file path, where % will be replaced with job ID
#BSUB -q verylong                     # queue 
#BSUB -N                              # end of the job enmail
#BSUB -u m.przybilla@dkfz-heidelberg.de   # user email
#BSUB -R rusage[mem=250GB]            # memory reservation limit
#BSUB -P Medulloblastoma

# load required module
module load R/4.0.0
module load jags/4.3.0

# Rscript /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/infercnv_analysis_ref.R /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/4M67_Nuclei_merged_WT
Rscript /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/infercnv_analysis_ref.R /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/STP_Nuclei_clusters_ungroup
# Rscript /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/infercnv_analysis_ref.R /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/STP_PDX_merged_WT
Rscript /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/infercnv_analysis_ref.R /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/ST1R_Nuclei_clusters_ungroup
# Rscript /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/infercnv_analysis_ref.R /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/ST1R_PDX_merged_WT
#Rscript /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/infercnv_analysis_ref.R /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/MB243_Nuclei_merged_WT
# Rscript /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/infercnv_analysis_ref.R /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/BT084_PDX_merged_WT
# Rscript /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/infercnv_analysis_ref.R /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/RCMB18_PDX_merged_WT