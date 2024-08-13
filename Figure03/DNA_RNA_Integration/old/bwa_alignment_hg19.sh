#!/bin/sh
# LSF directivesy
#BSUB -J bwa_align_hg19
#BSUB -e bwa_hg19-%J-%I.err                 # error file path, where % will be replaced with job ID
#BSUB -o bwa_hg19-%J-%I.out                 # outut file path, where % will be replaced with job ID
#BSUB -q verylong                     		# queue 
#BSUB -N                              		# end of the job enmail
#BSUB -u m.przybilla@dkfz-heidelberg.de   	# user email
#BSUB -n 8                            		# one CPU
#BSUB -M 20GB                        		# memory limit
#BSUB -R "rusage[mem=20GB]"

echo + `date` $LSB_JOBNAME started on $HOSTNAME jobID=$LSB_JOBID and taskID=$LSB_JOBINDEX

module load bwa
module load samtools

OUTPUT_PATH='/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/bwa'
mkdir $OUTPUT_PATH
REF='/omics/groups/OE0540/internal/resource/human/hg19/cellranger_reference/refdata-GRCh37-1.0.0/fasta/genome.fa'

## tab separated file: 1-sample name, 2-path to read 1, 3-path to read 2, 4- position on 96-well plate
SAMPLELIST=$(cat /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/scDNA_bwa_samplesheet.txt | cut -f 1 | sort -u)

for SAMPLE in $SAMPLELIST
do
    echo "starting with sample" $SAMPLE

    mkdir $OUTPUT_PATH/$SAMPLE
    READ1=$(grep -w $SAMPLE /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/scDNA_bwa_samplesheet.txt | cut -f 2 )
    READ2=$(grep -w $SAMPLE /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/scDNA_bwa_samplesheet.txt | cut -f 3 )

    cd $OUTPUT_PATH/$SAMPLE

    bwa mem -a -t 8 $REF $READ1 $READ2 > $OUTPUT_PATH/$SAMPLE/$SAMPLE.sam

    echo 'bwa sam file created'

    samtools import $REF.fai $OUTPUT_PATH//$SAMPLE/$SAMPLE.sam $OUTPUT_PATH/$SAMPLE/$SAMPLE.bam

    echo 'bam file created from sam'

    samtools sort -o $OUTPUT_PATH/$SAMPLE/$SAMPLE.sorted.bam $OUTPUT_PATH/$SAMPLE/$SAMPLE.bam

    echo 'bam file sorted'

    samtools index $OUTPUT_PATH/$SAMPLE/$SAMPLE.sorted.bam

    echo 'bam file indexed'

    echo 'aligned and done!'

    rm $OUTPUT_PATH/$SAMPLE/$SAMPLE.sam
    rm $OUTPUT_PATH/$SAMPLE/$SAMPLE.bam

    cd ..

done



