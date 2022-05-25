configfile: "Nanopore_SV.yaml"

rule all:
    input:
        "VCF/SVs.vcf",
        "QC/NanoQC_report",
        "QC/NanoPlot_before_Filt_report",
        "QC/NanoPlot_after_Filt_report",
        "Map/aln_sorted.bam.bai",
        expand("{REF}.{IDX}",IDX=["fai"],REF=config["reference"]["Ref"])
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

rule NanoQC_unfilted:
    input:
        fq = config["samples"]["single_fq"]
    output:
        report = directory("QC/NanoQC_report")
    log:
        "log/NanoQC.log"
    threads: 1
    shell:
        "nanoQC -o {output.report} {input}"
    
rule NanoPlot_before:
    input:
        fq = config["samples"]["single_fq"]
    output:
        report = directory("QC/NanoPlot_before_Filt_report")
    log:
        "log/NanoPlot_before.log"
    threads: 4
    shell:
        "NanoPlot -t {threads} --N50 --title NanoPlot_before_Filt --plots dot kde --fastq {input} -o {output.report}"
        
rule Filtlong:
    input:
        fq = config["samples"]["single_fq"]
    params:
        min_length = config["minlength"]["min_len"]
    output:
        "Filted/clean.filtlong.fq.gz"
    log:
        "log/Filtlong.log"
    threads: 4
    shell:
        "filtlong -p 90 --min_length {params.min_length} --min_mean_q 7 {input.fq} |  pigz -p {threads} > {output} "

rule NanoPlot_after:
    input:
        "Filted/clean.filtlong.fq.gz"
    output:
        report = directory("QC/NanoPlot_after_Filt_report")
    log:
        "log/NanoPlot_after.log"
    threads: 4
    shell:
        "NanoPlot -t {threads} --N50 --title NanoPlot_after_Filt --plots dot kde  --fastq {input} -o {output.report}"  

# rule NGMLR_map:
#     input:
#         filted_fq = "Filted/clean.filtlong.fq.gz",
#         ref = config["reference"]["Ref"]
#     output:
#         temp("Map/ngmlr.sam")
#     log:
#         "log/NGMLR.log"
#     shell:
#         "ngmlr -t 10 -r {input.ref} -q {input.filted_fq} -o {output}"
rule Minimap2_map:
    input:
        filted_fq = "Filted/clean.filtlong.fq.gz",
        ref = config["reference"]["Ref"]
    output:
        temp("Map/aln.sam")
    log:
        "log/aln.log"
    threads: 4
    shell:
        "minimap2 --MD -t {threads} -ax map-ont {input.ref} {input.filted_fq} > {output}"

#rule NGMLR_map:
#    input:
#        filted_fq = "Filted/clean.filtlong.fq.gz",
#        ref = config["reference"]["Ref"]
#    output:
#        temp("Map/aln.sam")
#    log:
#        "log/NGMLR.log"
#    shell:
#        "ngmlr -x ont --bam-fix -t 16 -r {input.ref} -q {input.filted_fq} -o {output}"
        
rule samtools_buildindex:
    input:
        ref = config["reference"]["Ref"]
    output:
        index_files = expand("{REF}.{IDX}",IDX=["fai"],REF=config["reference"]["Ref"])
    log:
        "log/build_index.log"
    threads: 1
    shell:
        "samtools faidx {input.ref}"

rule sam2bam:
    input:
        sam = "Map/aln.sam",
        ref = config["reference"]["Ref"]
    output:
        bam = temp("Map/aln.bam")
    log:
        "log/sam2bam.log"
    threads: 1
    shell:
        "samtools view -bh -t {input.ref}.fai -o {output.bam} {input.sam}"

        
rule samtools_sort:
    input:
        bam = "Map/aln.bam"
    output:
        sorted_bam = "Map/aln_sorted.bam"
    log:
        "log/samtools_sort.log"
    threads: 4
    shell:
        "samtools sort -@ {threads} -O bam -o {output.sorted_bam} {input.bam}"

rule samtools_index:
    input:
        sorted_bam = "Map/aln_sorted.bam"
    output:
        bam_ind = "Map/aln_sorted.bam.bai"
    log:
        "log/samtoos_index.log"
    threads: 1
    shell:
        "samtools index {input.sorted_bam} {output.bam_ind}"
        
rule sniffles:
    input:
        "Map/aln_sorted.bam"
    output:
        "VCF/SVs.vcf"
    log:
        "log/sniffles.log"
    threads: 4
    shell:
        "sniffles -t {threads} -m {input} -v {output}"
