#!/bin/sh
# LSF directives 
#BSUB -J Rtsne
#BSUB -e T-%J-%I.err                  # error file path, where % will be replaced with job ID
#BSUB -o T-%J-%I.out                  # outut file path, where % will be replaced with job ID
#BSUB -q verylong                     # queue 
#BSUB -N                              # end of the job enmail
#BSUB -u m.przybilla@dkfz-heidelberg.de   # user email
#BSUB -n 4                           # one CPU
#BSUB -M 200GB                        # memory limit
#BSUB -R rusage[mem=15GB]            # memory reservation limit
#BSUB -P Medulloblastoma


# load required module
module load R/3.6.2

Rscript /icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/infercnv_MB/hmm_tsne/hmm_clustering.R
