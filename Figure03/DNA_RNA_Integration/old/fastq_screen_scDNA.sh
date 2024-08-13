#!/bin/sh
# LSF directivesy
#BSUB -J fastq_screen_gtDNA"
#BSUB -e fastq_screen_gtDNA-%J.err                 # error file path, where % will be replaced with job ID
#BSUB -o fastq_screen_gtDNA-%J.out                 # outut file path, where % will be replaced with job ID
#BSUB -q verylong                     		# queue 
#BSUB -N                              		# end of the job enmail
#BSUB -u m.przybilla@dkfz-heidelberg.de   	# user email
#BSUB -n 8                            		# one CPU
#BSUB -M 10GB                        		# memory limit
#BSUB -R "rusage[mem=10GB]"

echo + `date` $LSB_JOBNAME started on $HOSTNAME jobID=$LSB_JOBID 

module load bwa
module unload python

OUTPUT_PATH='/omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/fastq_screen'
mkdir $OUTPUT_PATH

CONFIG='/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/analysis_Anna/scDNA/fastq_screen/fastq_screen.conf'
## human ref: '/icgc/dkfzlsdf/analysis/B260/resource/human/hg19/cellranger_reference/refdata-GRCh37-1.0.0/fasta/genome.fa'
## mouse ref: /icgc/dkfzlsdf/analysis/B260/resource/mouse/mm10/index_ucsc_mm10/mm10.UCSC.genome.chr1-M.fa

SAMPLELIST=$(cat /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/scDNA_trimgalore_samplesheet.txt | cut -f 1 | sort -u)
for SAMPLE in $SAMPLELIST
do
    echo "starting with sample" $SAMPLE

    mkdir $OUTPUT_PATH/$SAMPLE
    READ1=$(grep -w $SAMPLE /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/scDNA_trimgalore_samplesheet.txt | cut -f 2 | sort -u)
    READ2=$(grep -w $SAMPLE /omics/groups/OE0540/internal/projects/przybilm/medulloblastoma/revision/gt_seq/scDNA_384/scDNA_trimgalore_samplesheet.txt | cut -f 3 | sort -u)
	
    fastq_screen --threads 8 --aligner 'bwa' --conf $CONFIG --outdir $OUTPUT_PATH $READ1 $READ2

    echo "Done with sample $READ1"

done