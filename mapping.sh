module load bwa
module load samtools/gcc/64/1.5

# install samblaster and qualimap

## Read trimming (low quality)
#java -jar $TRIM/trimmomatic-0.38.jar PE $READ1 $READ2 trimmed_$READ1 trimmed_$READ2 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

 

unset DISPLAY

 

CPUs=8
READ1=$1 ## Fastq file with FW reads
READ2=$2 ## Fastq file with RV reads
REF=$3 ## Reference genome fasta file
BAM=$4 $$ ## Name prefix before .bam

 

time /bwa mem -t $CPUs $REF $READ1 $READ2 | /samblaster -r | samtools view -b - > $BAM.bam ## Process alignment and remove duplicates
samtools fixmate $BAM.bam $BAM.fixmate.bam ## Fill in mate coordinates, ISIZE and mate related flags from a name-sorted alignment.
samtools sort -m 2G -@ 6 -O bam $BAM.fixmate.bam > $BAM.fixmate.sort.bam ## Sort file 
samtools index -@ 4 $BAM.fixmate.sort.bam ## Index file
$QUALIMAP/qualimap bamqc -bam $BAM.fixmate.sort.bam --java-mem-size=5G -nt 1 -outformat PDF ## Quality report