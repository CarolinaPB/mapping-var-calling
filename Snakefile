configfile: "config.yaml"
import os

direction=["R1","R2"]

rule bwa_map:
    "Index, align reads and remove duplicates"
    input:
        assembly = os.path.join(config["DATADIR"], config["assembly"]),
        reads = expand(os.path.join(config["DATADIR"],"SG_data/DTG-SG-188_{dir}_001.fastq.gz"), dir=direction)
    output:
        expand("mapped_reads/DTG-SG-188_{dir}_001.fastq.gz.bam", dir=direction)
    threads: 8
    conda:
        "envs/bwa.yaml"
    message:
        "Rule {rule} processing"
    shell:
        """
        module load bwa
        module load samtools/gcc/64/1.5

        bwa index {input.assembly} | 
        bwa mem -t {threads} {input.assembly} {input.reads} | samblaster -r | samtools view -b - > {output}
        """

rule samtools_fixmate:
    input: 
        "mapped_reads/{sample}.bam"
    output: 
        temp("fixmate/{sample}.fixmate.bam")
    message:
        "Rule {rule} processing"
    shell: 
        "samtools fixmate {input} {output}"


# rule samtools_sort_index:
#     input: 
#         "fixmate/{sample}.fixmate.bam"
#     output: 
#         "sorted_reads/{sample}.fixmate.sort.bam"
    # message:
    #     "Rule {rule} processing"
#     shell: 
#         "samtools sort -m 2G -@ 6 -O bam {input} > {output} | samtools index -@ 4 {output}"


# rule qualimap_report:
#     input: 
#         "sorted_reads/{sample}.fixmate.sort.bam"
#     output: 
#         "qualimap_report/report.pdf"
    # message:
    #     "Rule {rule} processing"
#     shell: 
#         "qualimap bamqc -bam {input} --java-mem-size=5G -nt 1 -outfile {output}"

    