#!/bin/sh
# LSF directives 
#BSUB -J 4M67-Nuclei_infercnv
#BSUB -e T-%J-%I.err                  # error file path, where % will be replaced with job ID
#BSUB -o T-%J-%I.out                  # outut file path, where % will be replaced with job ID
#BSUB -q verylong                     # queue 
#BSUB -N                              # end of the job enmail
#BSUB -u m.przybilla@dkfz-heidelberg.de   # user email
#BSUB -R rusage[mem=75GB]            # memory reservation limit
#BSUB -P Medulloblastoma

# load required module
module load R/4.0.0
module load jags/4.3.0

Rscript /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/infercnv_analysis_nuclei.R /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/freeze/4M67-Nuclei_tumor
