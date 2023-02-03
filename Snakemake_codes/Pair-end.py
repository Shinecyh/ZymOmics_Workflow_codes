configfile: "paired.yaml"

rule all:
    input:
        expand("count/Result.count",sample=config["samples"])
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
        
# def get_input_fastqs(wildcards):
#     return config["samples"][wildcards.sample]
rule fastp_trim:
    input:
        # get_input_fastqs
        fq_1 = config["samples"]["paired_fq_1"],
        fq_2 = config["samples"]["paired_fq_2"]
    output:
        trimmed_1 = "trim/trimmed_R1.fq.gz",
        trimmed_2 = "trim/trimmed_R2.fq.gz",
        html = "trim/fastp.html",
        json = "trim/fastp.json"
    log:
        "log/fastp.log"
    shell:
        "fastp --compression=6 --thread=4 -R Trim_report --html {output.html} --json {output.json} -i {input.fq_1} -I {input.fq_2} -o {output.trimmed_1} -O {output.trimmed_2}"


rule subread_buildindex:
    input:
        INDEX = config["reference"]["bwaRef"]
    output:
        index_files = expand("ind.{IDX}",IDX=["00.b.array", "00.b.tab", "files", "log", "reads"])
    shell:
        "subread-buildindex -o ind {input}"

rule subread_map:
    input:
        trim_1 = rules.fastp_trim.output.trimmed_1,
        trim_2 = rules.fastp_trim.output.trimmed_2,
        INDEX= rules.subread_buildindex.output
    output:
        temp("map/mapped.sam")
    log:
        "log/bwamap.log"
    shell:
        "subjunc -T 4 -i ind -r {input.trim_1} -R {input.trim_2} -o {output}"

rule sam2bam:
    input:
        "map/mapped.sam"
    output:
        "bam/mapped.bam"
    log:
        "log/sam2bam.log"
    shell:
        "samtools view -@ 4 -Sb {input} > {output}"

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
        sample="bam/mapped.bam"
    output:
        "count/Result.count"
    params:
        FEATURE=config["feature"]["features"],
        GTF=config["annotation"]["gtf"]
    log:
        "log/count.log"
    shell:
        "featureCounts -T 4 -O -t {params.FEATURE}  -a {params.GTF} -o {output} {input.sample}"
