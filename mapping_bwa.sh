#!/bin/bash
#SBATCH --time=30
#SBATCH --mem=12000
#SBATCH --cpus-per-task=12
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=carolina.moitinhadealmeidapitabarros@wur.nl

module load bwa
module load samtools/gcc/64/1.5

CPUs=8
READ1=/lustre/backup/WUR/ABGC/shared/ABGC_Projects/Turkey_Assembly/SG_data/DTG-SG-188_R1_001.fastq.gz ## Fastq file with FW reads
READ2=/lustre/backup/WUR/ABGC/shared/ABGC_Projects/Turkey_Assembly/SG_data/DTG-SG-188_R2_001.fastq.gz ## Fastq file with RV reads
REF=/lustre/backup/WUR/ABGC/shared/ABGC_Projects/Turkey_Assembly/HiC_assembly/Mgal6.fa ## Reference genome fasta file
BAM=DTG-SG-188 ## Name prefix before .bam

 

time bwa mem -t $CPUs $REF $READ1 $READ2 | samblaster -r | samtools view -b - > $BAM.bam ##