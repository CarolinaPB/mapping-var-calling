configfile: "config.yaml"
import os

wdir=os.getcwd()

# Get all reads that are in this dir
reads, = glob_wildcards(os.path.join(config["DATADIR"], config["READS_DIR"], "{sample}.subset.fastq.gz"))
ASSEMBLY=os.path.join(config["DATADIR"], config["assembly"])

# Resources that can be set individually in each rule:
    # resources:
        # time_min=XXX (default 180),
        # cpus=XXX (default 16),
        # mem_mb=XXX (default 16000)

rule all:
    input:
        "variant_calling/var.vcf.gz", 
        "results/qualimap/report.pdf"

localrules: move_qualimap_res


# run bwa index if the index files don't exist
if not os.path.isfile(os.path.join(ASSEMBLY+".amb")):
    rule bwa_index:
        input: 
            ASSEMBLY
        output:
            os.path.join(ASSEMBLY+".amb")
        shell:
            "module load bwa && bwa index {input}"

rule bwa_map:
    # Index, align reads and remove duplicates
    input:
        # check = "checks/bwa_index.txt",
        assembly = ASSEMBLY,
        reads=expand(os.path.join(config["DATADIR"], "SG_data/", "{sample}.subset.fastq.gz"), sample=reads)

    output:
        os.path.join("mapped_reads/", config["my_prefix"]+".bam")
    resources: 
        time_min=120,
        cpus=16,
        mem_mb=16000
    # conda:
    #     "envs/bwa.yaml" # for samblaster
    message:
        "Rule {rule} processing"
    shell:
        """
        module load bwa samtools && bwa mem -t {resources.cpus} {input.assembly} {input.reads} | samblaster -r | samtools view -b - > {output}
        """

rule samtools_fixmate:
    input: 
        rules.bwa_map.output
    output: 
        os.path.join("fixmate/", config["my_prefix"] +".fixmate.bam")
    resources:
        time_min=10
    message:
        "Rule {rule} processing"
    shell: 
        "module load samtools && samtools fixmate {input} {output}"


rule samtools_sort:
    input: 
        rules.samtools_fixmate.output
    output: 
        os.path.join("sorted_reads/", config["my_prefix"] +".fixmate.sort.bam")
    resources:
        time_min=5
    message:
        "Rule {rule} processing"
    shell: 
        "module load samtools && samtools sort -m 2G -@ 6 -O bam {input} > {output}"

rule samtools_index:
    input:
        rules.samtools_sort.output
    output:
        os.path.join("sorted_reads/", config["my_prefix"] +".fixmate.sort.bam.bai")
    resources:
        time_min=5
    message:
        "Rule {rule} processing"
    shell:
        "module load samtools && samtools index -@ 4 {input}"

rule qualimap_report:
    input: 
        check=rules.samtools_index.output, # not used in the command, but it's here so snakemake knows to run the rule after the indexing
        bam=rules.samtools_sort.output
    output: 
        temp(os.path.join(wdir, "sorted_reads/", config["my_prefix"] +".fixmate.sort_stats", "report.pdf"))
    # conda:
    #     "envs/qualimap.yaml"
    resources:
        time_min=10,
        cpus=1,
        mem_mb=2000
    message:
        "Rule {rule} processing"
    shell: 
        "unset DISPLAY && qualimap bamqc -bam {input.bam} --java-mem-size=5G -nt 1 -outformat PDF"

rule move_qualimap_res:
# Move the qualimap results to a dir easier to find
    input: 
        rules.qualimap_report.output
    output:
        file="results/qualimap/report.pdf",
        newdir= directory("results/qualimap/")
    resources:
        time_min=2,
        cpus=1
    shell:
        "mv sorted_reads/*_stats/* {output.newdir} && sleep 10"

rule freebayes_var:
    input: 
        reference= ASSEMBLY,
        bam = rules.samtools_sort.output, 
        bam_bai = rules.samtools_index.output # not used in the command, but it's here so snakemake knows to run the rule after the indexing
    output: 
        "variant_calling/var.vcf.gz"
    resources:
        time_min=40
    message:
        "Rule {rule} processing"
    shell:
        "module load freebayes samtools vcflib/gcc/64/0.00.2019.07.10 && freebayes -f {input.reference} --use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 --bam {input.bam} | vcffilter -f 'QUAL > 20' | bgzip -c > {output}"
