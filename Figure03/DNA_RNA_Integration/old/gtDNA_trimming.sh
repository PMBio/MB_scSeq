#!/bin/sh
# LSF directivesy
#BSUB -J scDNA_trim
#BSUB -e trim-%J.err                 		# error file path, where % will be replaced with job ID
#BSUB -o trim-%J.out                  		# outut file path, where % will be replaced with job ID
#BSUB -q verylong                     		# queue 
#BSUB -N                              		# end of the job enmail
#BSUB -u m.przybilla@dkfz-heidelberg.de   	# user email
#BSUB -n 4                            		# one CPU
#BSUB -M 10GB                        		# memory limit
#BSUB -R "rusage[mem=10GB]"
#BSUB -P chromotr_GT-seq_dna_trim-galore

module load fastqc
module load trim-galore/0.5.0
module unload python

# OUTPUT_PATH='/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA/trim-galore'
OUTPUT_PATH='/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/trim-galore'
mkdir $OUTPUT_PATH

SAMPLELIST=$(cat /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/scDNA_trimgalore_samplesheet.txt | cut -f 1 | sort -u)
for SAMPLE in $SAMPLELIST
do
    echo "starting with sample" $SAMPLE

    mkdir $OUTPUT_PATH/$SAMPLE
    READ1=$(grep -w $SAMPLE /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/scDNA_trimgalore_samplesheet.txt | cut -f 2 | sort -u)
    READ2=$(grep -w $SAMPLE /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/scDNA_trimgalore_samplesheet.txt | cut -f 3 | sort -u)
    
	  # trim_galore --nextera --fastqc $READ1 -o $OUTPUT_PATH/$SAMPLE
    trim_galore --paired --nextera --fastqc $READ1  $READ2 -o $OUTPUT_PATH/$SAMPLE
    
    echo "Done with sample $READ1"

done



