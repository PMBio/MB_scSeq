#!/bin/sh
# LSF directivesy
#BSUB -J STAR-align
#BSUB -e STAR-align-%J.err                  # error file path, where % will be replaced with job ID
#BSUB -o STAR-align-%J.out                  # outut file path, where % will be replaced with job ID
#BSUB -q verylong                     # queue 
#BSUB -N                              # end of the job enmail
#BSUB -u m.przybilla@dkfz-heidelberg.de   # user email
#BSUB -n 8                            # one CPU
#BSUB -M 20GB                        # memory limit
#BSUB -R "rusage[mem=20GB]"

echo + `date` $LSB_JOBNAME started on $HOSTNAME jobID=$LSB_JOBID and taskID=$LSB_JOBINDEX

module load STAR/2.5.3a
module load samtools/1.5

CPU=8
OUTPUT_PATH='/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_384/STAR'
mkdir $OUTPUT_PATH
GenomeDir="/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/analysis_Anna/scRNA/STAR2_align_2/GenomeForPass2"

## tab separated file: 1-sample name, 2-path to read 1, 3-path to read 2, 4- position on 96-well plate
SAMPLELIST=$(cat /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_384/scRNA_star_samplesheet.txt | cut -f 1 | sort -u)

for SAMPLE in $SAMPLELIST
do
    echo "starting with sample" $SAMPLE

    mkdir $OUTPUT_PATH/$SAMPLE
    READ1=$(grep -w $SAMPLE /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_384/scRNA_star_samplesheet.txt | cut -f 2 | sort -u)
    echo $READ1

    cd $OUTPUT_PATH/$SAMPLE

    STAR --runThreadN 8 --outSAMattributes All --genomeLoad LoadAndKeep --readFilesCommand zcat --genomeDir $GenomeDir --readFilesIn $READ1 --outFileNamePrefix $SAMPLE.

    echo 'STAR sam file created'

    samtools sort -m 10G -@ $CPU -O bam -o $SAMPLE'.Aligned.sortedByCoord.out.bam' $SAMPLE.Aligned.out.sam
    samtools index -@ $CPU  $SAMPLE'.Aligned.sortedByCoord.out.bam'

    echo 'bam file sorted and indexed'

    samtools stats $SAMPLE'.Aligned.sortedByCoord.out.bam' > $SAMPLE.samtools.stats.txt
    samtools flagstat $SAMPLE'.Aligned.sortedByCoord.out.bam' > $SAMPLE.samtools.flagstat.txt
    samtools idxstats $SAMPLE'.Aligned.sortedByCoord.out.bam' > $SAMPLE.samtools.idxstats.txt

    echo 'aligned and done!'

    cd ..
done
