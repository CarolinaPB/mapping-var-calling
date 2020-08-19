configfile: "config.yaml"
import os

wdir=os.getcwd()
# workdir: /home/WUR/moiti001/snakemake-try

# Resources that can be set individually in each rule:
    # Resources:
        # time_min=XXX (default 180)
        # cpus=XXX (default 16)
        # mem_mb=XXX (default 16000)

rule all:
    input:
        os.path.join(wdir, "sorted_reads/", config["BAM_prefix"] +".fixmate.sort_stats", "qualimap_report.pdf"),
        "variant_calling/var.vcf.gz"

rule bwa_index:
    # check if the index files exist, if not, run bwa index
    input:
        os.path.join(config["DATADIR"], config["assembly"])
    output:
        "checks/bwa_index.txt"
    message:
        "Rule {rule} processing"
    shell:
        "scripts/bwa_index_check.sh {input} > {output}"

rule bwa_map:
    # Index, align reads and remove duplicates
    input:
        rules.bwa_index.output,
        assembly = os.path.join(config["DATADIR"], config["assembly"]),
        fwd = os.path.join(config["DATADIR"], config["fwd"]),
        rev = os.path.join(config["DATADIR"], config["rev"])
    output:
        os.path.join("mapped_reads/", config["BAM_prefix"]+".bam")
    resources: 
        cpus=16,
        mem_mb=16000
    message:
        "Rule {rule} processing"
    shell:
        """
        module load bwa samtools && bwa mem -t {resources.cpus} {input.assembly} {input.fwd} {input.rev}| samblaster -r | samtools view -b - > {output}
        """

rule samtools_fixmate:
    input: 
        os.path.join("mapped_reads/", config["BAM_prefix"]+".bam")
    output: 
        os.path.join("fixmate/", config["BAM_prefix"] +".fixmate.bam")
    envmodules: 
        "samtools/gcc/64/1.9"
    message:
        "Rule {rule} processing"
    shell: 
        "module load samtools && samtools fixmate {input} {output}"


rule samtools_sort:
    input: 
        os.path.join("fixmate/", config["BAM_prefix"] +".fixmate.bam")
    output: 
        os.path.join("sorted_reads/", config["BAM_prefix"] +".fixmate.sort.bam")
    envmodules:
        "samtools/gcc/64/1.9"
    message:
        "Rule {rule} processing"
    shell: 
        "module load samtools && samtools sort -m 2G -@ 6 -O bam {input} > {output}"

rule samtools_index:
    input:
        os.path.join("sorted_reads/", config["BAM_prefix"] +".fixmate.sort.bam")
    output:
        "checks/samtools_index.txt", 
        os.path.join("sorted_reads/", config["BAM_prefix"] +".fixmate.sort.bam.bai")
    message:
        "Rule {rule} processing"
    shell:
        "module load samtools && samtools index -@ 4 {input} && echo '{rule} done' > {output}"

rule qualimap_report:
    input: 
        check="checks/samtools_index.txt",
        bam=os.path.join("sorted_reads/", config["BAM_prefix"] +".fixmate.sort.bam")
    output: 
        os.path.join(wdir, "sorted_reads/", config["BAM_prefix"] +".fixmate.sort_stats", "qualimap_report.pdf")
    conda:
        "envs/qualimap.yaml"
    message:
        "Rule {rule} processing"
    shell: 
        #"unset DISPLAY && qualimap bamqc -bam {input.bam} --java-mem-size=5G -nt 1 -outdir qualimap_report -outfile {output}"
        "unset DISPLAY && qualimap bamqc -bam {input.bam} --java-mem-size=5G -nt 1 -outformat PDF"

rule freebayes_var:
    input: 
        reference= os.path.join(config["DATADIR"], config["assembly"]),
        bam = os.path.join("sorted_reads/", config["BAM_prefix"] +".fixmate.sort.bam"), 
        bam_bai = os.path.join("sorted_reads/", config["BAM_prefix"] +".fixmate.sort.bam.bai")
    output: 
        "variant_calling/var.vcf.gz"
    message:
        "Rule {rule} processing"
    shell:
        "module load freebayes samtools vcflib/gcc/64/0.00.2019.07.10 && freebayes -f {input.reference} --use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 --bam {input.bam} | vcffilter -f 'QUAL > 20' | bgzip -c > {output}"
