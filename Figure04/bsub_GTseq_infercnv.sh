#!/bin/sh
# LSF directives 
#BSUB -J GTseq_complete_infercnv
#BSUB -e T-%J-%I.err                  # error file path, where % will be replaced with job ID
#BSUB -o T-%J-%I.out                  # outut file path, where % will be replaced with job ID
#BSUB -q verylong                     # queue 
#BSUB -N                              # end of the job enmail
#BSUB -u m.przybilla@dkfz-heidelberg.de   # user email
#BSUB -R rusage[mem=50GB]            # memory reservation limit
#BSUB -P Medulloblastoma

# load required module
module load R/4.0.0
module load jags/4.3.0

Rscript /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_analysis_GTseq.R /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/scrna_analysis/infercnv_MB/GTseq_384_mouse_new