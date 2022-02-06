#!/bin/sh
# LSF directivesy
#BSUB -J HTseq
#BSUB -e HTseq-%J.err                  # error file path, where % will be replaced with job ID
#BSUB -o HTseq-%J.out                  # outut file path, where % will be replaced with job ID
#BSUB -q verylong                     # queue 
#BSUB -N                              # end of the job enmail
#BSUB -u m.przybilla@dkfz-heidelberg.de   # user email
#BSUB -n 8                            # one CPU
#BSUB -M 20GB                        # memory limit
#BSUB -R "rusage[mem=20GB]"

echo + `date` $LSB_JOBNAME started on $HOSTNAME jobID=$LSB_JOBID and taskID=$LSB_JOBINDEX

module load samtools/1.5
module unload python

MAP_Q=10          # only unique reads (MAP_Q==255)
GTFfile='/omics/groups/OE0540/internal/resource/human/hg19/star-genome/hg19/genes/genes.gtf'
OUTPUT_PATH='/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_384/HTseq'
DATA_PATH='/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scRNA_384/STAR'
mkdir $OUTPUT_PATH
cd $OUTPUT_PATH

SAMPLELIST=$(find $DATA_PATH -name '*.Aligned.sortedByCoord.out.bam')
for SAMPLE in $SAMPLELIST
do echo $SAMPLE # std out PROJECT.ID
y=${SAMPLE%.Aligned.sortedByCoord.out.bam} ; ## this gets the whole path before .bam
z=${y##*/} ;
samtools view -q $MAP_Q $SAMPLE | python -m HTSeq.scripts.count --mode=intersection-strict --stranded=no --type=exon --idattr=gene_id -a $MAP_Q --order=pos - $GTFfile > $OUTPUT_PATH/${z}.htseq_counts.out
done
