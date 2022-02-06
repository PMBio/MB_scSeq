#!/bin/sh
# LSF directivesy
#BSUB -J scRNA_fastqc
#BSUB -e scRNA_fastqc-%J.err                  # error file path, where % will be replaced with job ID
#BSUB -o scRNA_fastqc-%J.out                  # outut file path, where % will be replaced with job ID
#BSUB -q verylong                     # queue 
#BSUB -N                              # end of the job enmail
#BSUB -u m.przybilla@dkfz-heidelberg.de   # user email
#BSUB -n 4                            # one CPU
#BSUB -M 10GB                        # memory limit
#BSUB -R "rusage[mem=10GB]"

module load fastqc/0.11.5 

OUTPUT_PATH='/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_384/fastqc'
mkdir $OUTPUT_PATH

SAMPLELIST=$(cat /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_384/scRNA_trimgalore_samplesheet.txt | cut -f 1 | sort -u)
for SAMPLE in $SAMPLELIST
do
    echo "starting with sample" $SAMPLE

    mkdir $OUTPUT_PATH/$SAMPLE
    READ1=$(grep -w $SAMPLE /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_384/scRNA_trimgalore_samplesheet.txt | cut -f 2 | sort -u)
	
    fastqc -t 4 $READ1 -o $OUTPUT_PATH/$SAMPLE

    echo "Done with sample $READ1"

done



