configfile: "config.yaml"
import os

# wdir=os.getcwd()
# Sets the working directory to the directory where you want to save the results (set in the config file)
workdir: config["OUTDIR"]

## problem - if I remove the "subset" from {sample}.subset.fastq.gz, then both the subset and the others are going to be used ({sample}.subset.fastq.gz here in the next lines and also on the bwa_map rule) 
#- can instead have each read defined individually in the config file (FWD=PATH and REV=PATH)

# Get all reads that are in this dir
# reads, = glob_wildcards(os.path.join( config["READS_DIR"], "{sample}.subset.fastq.gz")) ### FOR SUBSET uncomment this line and comment out the next
reads, = glob_wildcards(os.path.join(config["READS_DIR"], "{sample}001.fastq.gz"))

ASSEMBLY=config["assembly"]

# Resources that can be set individually in each rule:
    # resources:
        # time_min=XXX (default 180),
        # cpus=XXX (default 16),
        # mem_mb=XXX (default 16000)


rule all:
    input:
        "variant_calling/var.vcf.gz", 
        "results/qualimap/genome_results.txt"

# localrules: move_qualimap_res


# run bwa index if the index files don't exist
if not os.path.isfile(os.path.join(ASSEMBLY+".amb")):
    rule bwa_index:
        input: 
            ASSEMBLY
        output:
            os.path.join(ASSEMBLY+".amb")
        group:
            "group_all"
        shell:
            "module load bwa && bwa index {input}"

rule bwa_map:
    # Index, align reads and remove duplicates
    input:
        assembly = ASSEMBLY,
        # reads=expand(os.path.join(config["READS_DIR"], "{sample}.subset.fastq.gz"), sample=reads) ### FOR SUBSET uncomment this line and comment out the next
        reads=expand(os.path.join(config["READS_DIR"], "{sample}001.fastq.gz"), sample=reads)
    output:
        temp(os.path.join("mapped_reads/", config["my_prefix"]+".bam"))
    resources: 
        cpus=16
    group:
        "group_all"
    message:
        "Rule {rule} processing"
    shell:
        "module load bwa samtools && bwa mem -t {resources.cpus} {input.assembly} {input.reads} | samblaster -r | samtools view -b - > {output}"

rule samtools_sort:
    input: 
        rules.bwa_map.output
    output: 
        os.path.join("sorted_reads/", config["my_prefix"] +".sort.bam")
    resources:
        cpus=7
    group:
        "group_all"
    message:
        "Rule {rule} processing"
    shell: 
        "module load samtools && samtools sort -m 2G -@ {resources.cpus} -O bam {input} > {output}"

rule samtools_index:
    input:
        rules.samtools_sort.output
    output:
        os.path.join("sorted_reads/", config["my_prefix"] +".sort.bam.bai")
    resources:
        cpus=16
    group:
        "group_all"
    message:
        "Rule {rule} processing"
    shell:
        "module load samtools && samtools index -@ {resources.cpus} {input}"

rule qualimap_report:
    input: 
        check=rules.samtools_index.output, # not used in the command, but it's here so snakemake knows to run the rule after the indexing
        bam=rules.samtools_sort.output
    output: 
        outdir=directory("results/qualimap/"),
        outfile="results/qualimap/genome_results.txt"
        # temp(os.path.join(wdir, "sorted_reads/", config["my_prefix"] +".sort_stats", "report.pdf"))
    # resources:
    #     time_min=10,
    #     cpus=1,
    #     mem_mb=2000
    group:
        "group_all"
    message:
        "Rule {rule} processing"
    shell: 
        "unset DISPLAY && qualimap bamqc -bam {input.bam} --java-mem-size=16G -nt 1 -outdir {output.outdir}"

rule freebayes_var:
    input: 
        reference= ASSEMBLY,
        bam = rules.samtools_sort.output, 
        bam_bai = rules.samtools_index.output # not used in the command, but it's here so snakemake knows to run the rule after the indexing
    output: 
        "variant_calling/var.vcf.gz"
    # resources:
    #     time_min=40
    group:
        "group_all"
    message:
        "Rule {rule} processing"
    shell:
        "module load freebayes samtools vcflib/gcc/64/0.00.2019.07.10 && freebayes -f {input.reference} --use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 --bam {input.bam} | vcffilter -f 'QUAL > 20' | bgzip -c > {output}"
