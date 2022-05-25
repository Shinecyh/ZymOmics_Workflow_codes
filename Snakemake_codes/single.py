configfile: "single.yaml"

rule all:
    input:
        expand("count/{sample}.count",sample=config["samples"]),
#rule all:
# input:
#  "count/merge.count"
#  expand("count/{sample}.count",sample=config["samples"])

# rule quality_check:
#     input:
#         lambda wildcards: config["samples"][wildcards.sample]
#     output:
#         "qc/{sample}"
#     shell:
#         "fastqc -o {output} {input}"
def get_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]
rule fastp_trim:
    input:
        get_input_fastqs
    output:
        trimmed = "trim/{sample}_trimmed.fq.gz",
        html = "trim/{sample}.html",
        json = "trim/{sample}.json"
    log:
        "log/fastp_{sample}.log"
    shell:
        "fastp --compression=6 --thread=4 -R {wildcards.sample}_report --html {output.html} --json {output.json} -i {input} -o {output.trimmed}"

rule bwa_index:
    input:
        INDEX = config["reference"]["bwaRef"]
    output:
        expand("{INDEX}.{IDX}",IDX=["amb", "ann", "bwt", "pac", "sa"],INDEX = config["reference"]["bwaRef"])
    shell:
        "bwa index {input}"

rule bwa_map:
    input:
        FILE = rules.fastp_trim.output.trimmed,
        INDEX= config["reference"]["bwaRef"],
        asd = expand("{INDEX}.{IDX}",IDX=["amb", "ann", "bwt", "pac", "sa"],INDEX = config["reference"]["bwaRef"])
    output:
        temp("bwa/{sample}.sam")
    log:
        "log/bwamap_{sample}.log"
    shell:
        "bwa mem -t 4 {input.INDEX}  {input.FILE} > {output}"

rule sam2bam:
    input:
        "bwa/{sample}.sam"
    output:
        "bam/{sample}.bam"
    log:
        "log/sam2bam_{sample}.log"
    shell:
        "samtools view -Sb {input} > {output}"

# rule sort:
#     input:
#         "bam/{sample}.bam"
#     output:
#         "sort/{sample}_sort.bam"
#     shell:
#         "samtools sort -l 9 -@ 4 -o {output} {input} "

# rule index:
#     input:
#         "sort/{sample}_sort.bam"
#     output:
#         "sort/{sample}.bam.bi"
#     shell:
#         "samtools index {input} {output}"

# rule count:
#     input:
#         sample="sort/{sample}_sort.bam",
#         annotation=config["annotation"]["gtf"]
#     output:
#         "count/{sample}.count"
#     params:
#         FEATURES=config["feature"]["features"]
#         GTF=config["annotation"]["gtf"]
#     shell:
#         "featureCounts -T 4 -O -t {params.FEATURES}  -a {params.GTF} -o {output} input.sample}"
# com=['featureCounts -T', str(threads),'-O','-t',feature_type,'-a',gtf_loc,'-o',i+".out",os.path.join(junc_dir,i+".bam")]
rule featurecounts_count:
    input:
        sample="bam/{sample}.bam"
    output:
        "count/{sample}.count"
    params:
        FEATURE=config["feature"]["features"],
        GTF=config["annotation"]["gtf"]
    log:
        "log/count_{sample}.log"
    shell:
        "featureCounts -T 4 -O -t {params.FEATURE}  -a {params.GTF} -o {output} {input.sample}"