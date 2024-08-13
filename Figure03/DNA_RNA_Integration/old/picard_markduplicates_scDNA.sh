#!/bin/sh
# LSF directivesy
#BSUB -J scDNA_picard_markduplicates
#BSUB -e picard_markduplicates-%J.err               		# error file path, where % will be replaced with job ID
#BSUB -o picard_markduplicates-%J.out                	# outut file path, where % will be replaced with job ID
#BSUB -q verylong                     		# queue 
#BSUB -N                              		# end of the job enmail
#BSUB -u m.przybilla@dkfz-heidelberg.de   	# user email
#BSUB -n 8                            		# one CPU
#BSUB -M 10GB                        		# memory limit
#BSUB -R "rusage[mem=10GB]"

echo + `date` $LSB_JOBNAME started on $HOSTNAME jobID=$LSB_JOBID

module unload python

DATA_PATH='/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/bwa'
OUTPUT_PATH='/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/MarkDuplicates'
mkdir $OUTPUT_PATH
cd $OUTPUT_PATH

SAMPLELIST=$(find $DATA_PATH -name '*.sorted.bam')
for SAMPLE in $SAMPLELIST
do echo $SAMPLE # std out PROJECT.ID
y=${SAMPLE%.sorted.bam} ; ## this gets the whole path before .bam
z=${y##*/} ;
picard MarkDuplicates I=$SAMPLE O=$OUTPUT_PATH/$z.sorted_dedup.bam M=$OUTPUT_PATH/$z.metrics.txt REMOVE_DUPLICATES=TRUE AS=TRUE
samtools index $OUTPUT_PATH/$z.sorted_dedup.bam
cd ..
done



