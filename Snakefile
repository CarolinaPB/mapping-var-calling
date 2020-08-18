configfile: "config.yaml"
import os

direction=["R1","R2"]
rule all:
    input:
        "qualimap_report/report.pdf"
        # "fixmate/DTG-SG-188_R1_001.fixmate.bam"
        # "sorted_reads/DTG-SG-188_R1_001.fixmate.sort.bam",
        # "sorted_reads/DTG-SG-188_R2_001.fixmate.sort.bam"
        # expand("sorted_reads/DTG-SG-188_{dir}_001.fixmate.sort.bam", dir=direction)
        # expand("mapped_reads/DTG-SG-188_{dir}_001.fastq.gz.bam", dir=direction)

localrules: qualimap_report



rule bwa_map:
    "Index, align reads and remove duplicates"
    input:
        assembly = os.path.join(config["DATADIR"], config["assembly"]),
        # reads = expand(os.path.join(config["DATADIR"],"SG_data/DTG-SG-188_{dir}_001.fastq.gz"), dir=direction)
        fwd = os.path.join(config["DATADIR"], config["fwd"]),
        rev = os.path.join(config["DATADIR"], config["rev"])
    output:
        # expand("mapped_reads/DTG-SG-188_{dir}_001.fastq.gz.bam", dir=direction)
        os.path.join("mapped_reads/", config["BAM_prefix"]+".bam")
    resources: 
        cpus=16,
        mem_mb=16000
    conda:
        "envs/bwa.yaml"
    # envmodules:
    #     "bwa/gcc/64/0.7.17", 
    #     "samtools/gcc/64/1.9"
    message:
        "Rule {rule} processing"
    shell:
        """
        bwa index {input.assembly} | 
        bwa mem -t {resources.cpus} {input.assembly} {input.fwd} {input.rev}| samblaster -r | samtools view -b - > {output}
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
    conda:
        "envs/bwa.yaml"  
    shell: 
        "samtools fixmate {input} {output}"


# rule samtools_sort_index:
#     input: 
#         "fixmate/{sample}.fixmate.bam"
#     output: 
#         "sorted_reads/{sample}.fixmate.sort.bam"
    # conda:
    #     "envs/bwa.yaml"
#     message:
#         "Rule {rule} processing"
#     shell: 
#         "samtools sort -m 2G -@ 6 -O bam {input} > {output} | samtools index -@ 4 {output}"

rule samtools_sort:
    input: 
        os.path.join("fixmate/", config["BAM_prefix"] +".fixmate.bam")
    output: 
        os.path.join("sorted_reads/", config["BAM_prefix"] +".fixmate.sort.bam")
    conda:
        "envs/bwa.yaml"
    envmodules:
        "samtools/gcc/64/1.9"
    message:
        "Rule {rule} processing"
    shell: 
        "samtools sort -m 2G -@ 6 -O bam {input} > {output}"

rule samtools_index:
    input:
        os.path.join("sorted_reads/", config["BAM_prefix"] +".fixmate.sort.bam")
    output:
        os.path.join("sorted_reads/", config["BAM_prefix"] +".fixmate.sort.bam")
    conda:
        "envs/bwa.yaml"
    envmodules:
        "samtools/gcc/64/1.9"
    message:
        "Rule {rule} processing"
    shell:
        "samtools index -@ 4 {input}"

rule qualimap_report:
    input: 
        os.path.join("sorted_reads/", config["BAM_prefix"] +".fixmate.sort.bam")
    output: 
        "qualimap_report/report.pdf"
    conda:
        "envs/qualimap.yaml"
    message:
        "Rule {rule} processing"
    shell: 
        "qualimap bamqc -bam {input} --java-mem-size=5G -nt 1 -outfile {output}"

    