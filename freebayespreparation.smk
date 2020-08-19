
rule freebayes_var:
input: 
    reference= os.path.join(config["DATADIR"], config["assembly"]),
    bam = rules.samtools_index.output
output: 
    "var.vcf"
conda:
    "envs/freebayes.yaml"
message:
    "Rule {rule} processing"
shell:
    "freebayes -f {input.reference} min-base-quality 10 min-alternate-fraction 0.2 haplotype-length 0 ploidy 2 min-alternate-count 2 {input.bam} > {output}"