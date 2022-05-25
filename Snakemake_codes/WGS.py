configfile: "WGS.yaml"

rule all:
    input:
        # index_files = expand("{INDEX}.{IDX}",IDX=["amb", "ann", "bwt", "pac", "sa"],INDEX = config["reference"]["bwaRef"]),
        expand("Call_variants/bwa_mem_snps_filtered.avinput.{IDX}",IDX=["exonic_variant_function", "variant_function"]),
        "Map/bwa_mem_sorted.bam.bai",
        "Map/bwa_mem_nopcr.bam.bai"

rule fastp_trim:
    input:
        fq_1 = config["samples"]["paired_fq_1"],
        fq_2 = config["samples"]["paired_fq_2"]
    output:
        trimmed_1 = "trim/trimmed_R1.fq.gz",
        trimmed_2 = "trim/trimmed_R2.fq.gz",
        html = "trim/fastp.html",
        json = "trim/fastp.json"
    log:
        "log/fastp.log"
    threads: 4
    shell:
        "fastp --compression=6 --thread={threads} -R Trim_report --html {output.html} --json {output.json} -i {input.fq_1} -I {input.fq_2} -f 15 -F 15-o {output.trimmed_1} -O {output.trimmed_2}"
        
rule BWA_buildindex:
    input:
        INDEX = config["reference"]["bwaRef"]
    output:
        index_files = expand("{INDEX}.{IDX}",IDX=["amb", "ann", "bwt", "pac", "sa"],INDEX = config["reference"]["bwaRef"])
    log:
        "log/bwa_index.log"
    threads: 1
    shell:
        "bwa index {input.INDEX}"
        
rule BWA_mem:
    input:
        trimmed_1 = "trim/trimmed_R1.fq.gz",
        trimmed_2 = "trim/trimmed_R2.fq.gz",
        INDEX = config["reference"]["bwaRef"],
        index_files = expand("{INDEX}.{IDX}",IDX=["amb", "ann", "bwt", "pac", "sa"],INDEX = config["reference"]["bwaRef"])
    output:
        temp("Map/bwa_mem.sam")
    log:
        "log/BWA_mem.log"
    threads: 4
    shell:
        "bwa mem -t {threads} {input.INDEX} {input.trimmed_1} {input.trimmed_2} > {output}"

rule samtools_buildindex:
    input:
        ref = config["reference"]["bwaRef"]
    output:
        index_files = expand("ref.{IDX}",IDX=["fai"])
    log:
        "log/build_index.log"
    threads: 1
    shell:
        "samtools faidx {input.ref}"
        
rule sam2bam:
    input:
        sam = "Map/bwa_mem.sam",
        ref = config["reference"]["bwaRef"]
    output:
        bam = temp("Map/bwa_mem.bam")
    log:
        "log/sam2bam.log"
    threads: 4
    shell:
        "samtools view -@ {threads} -bhS -t {input.ref}.fai -o {output.bam}  {input.sam}"

        
rule samtools_sort:
    input:
        bam = "Map/bwa_mem.bam"
    output:
        sorted_bam = temp("Map/bwa_mem_sorted.bam")
    log:
        "log/samtools_sort.log"
    threads: 4
    shell:
        "samtools sort -m 500M -@ {threads} -O bam -o {output.sorted_bam} {input.bam}"
        
rule samtools_index:
    input:
        sorted_bam = "Map/bwa_mem_sorted.bam"
    output:
        sorted_bam_index = "Map/bwa_mem_sorted.bam.bai"
    log:
        "log/samtools_index.log"
    threads: 4
    shell:
        "samtools index -@ {threads} {input.sorted_bam} {output.sorted_bam_index}"
        
rule samtools_stat:
    input:
        sorted_bam = "Map/bwa_mem_sorted.bam"
    output:
        depth = "Map/depth.txt",
        stat = "Map/stat.txt"
    log:
        "log/samtools_stat.log"
    threads: 1
    shell:
        "samtools depth {input.sorted_bam} >> {output.depth} ; samtools flagstat {input.sorted_bam} >> {output.stat}"

rule samtools_rmdup:
    input:
        sorted_bam = "Map/bwa_mem_sorted.bam"
    output:
        nopcr_bam = "Map/bwa_mem_nopcr.bam"
    log:
        "log/samtools_rmdup.log"
    threads: 1
    shell:
        "samtools rmdup {input.sorted_bam} {output.nopcr_bam}"
        # samtools markdup -r -@ 4 -O BAM --write-index Map/bwa_mem_sorted.bam Map/bwa_mem_nopcr.bam
        
rule samtools_index2:
    input:
        nopcr_bam = "Map/bwa_mem_nopcr.bam"
    output:
        nopcr_bam_index = "Map/bwa_mem_nopcr.bam.bai"
    log:
        "log/samtools_index2.log"
    threads: 4
    shell:
        "samtools index -@ {threads} {input.nopcr_bam} {output.nopcr_bam_index}"

rule samtools_callbcf:
    input:
        nopcr_bam = "Map/bwa_mem_nopcr.bam",
        ref = config["reference"]["bwaRef"]
    output:
        bcf = "Call_variants/bwa_mem.bcf"
    log:
        "log/samtools_callbcf.log"
    threads: 1
    shell:
        "bcftools mpileup -f {input.ref} {input.nopcr_bam} > {output.bcf}"

rule bcftools_call:
    input:
        bcf = "Call_variants/bwa_mem.bcf"
    output:
        variants_bcf = "Call_variants/bwa_mem_variants.bcf",
        snp_vcf = "Call_variants/bwa_mem_snps.vcf",
        filtered_vcf = "Call_variants/bwa_mem_snps_filtered.vcf"
    log:
        "log/bcftools_call.log"
    threads: 1
    shell:
        "bcftools  call -vm {input.bcf} -o {output.variants_bcf} && bcftools view -v snps,indels,mnps,ref,bnd {output.variants_bcf} > {output.snp_vcf} && bcftools filter -o {output.filtered_vcf} -i 'QUAL>20 &&DP>5' {output.snp_vcf}"
        
rule convert2annovar:
    input:
        filtered_vcf = "Call_variants/bwa_mem_snps_filtered.vcf"
    output:
        avinput = "Call_variants/bwa_mem_snps_filtered.avinput"
    log:
        "log/convert2annovar.log"
    threads: 1
    shell:
       " convert2annovar.pl -format vcf4 {input.filtered_vcf} > {output.avinput}"

rule make_refdb:
    input:
        gff = config["annotation"]["gff"],
        INDEX = config["reference"]["bwaRef"]
    output:
        gff_refGene = "annofile/refGene.txt",
        col1 = "annofile/column1.txt",
        col_else = "annofile/column_else.txt",
        ref_refGene =  "annofile/ref-genome_new_refGene.txt",
        ref_mrna = "annofile/ref-genome_new_refGeneMrna.fa"
    log:
        "log/make_refdb.log"        
    threads: 1
    shell:
        "gff3ToGenePred -useName {input.gff} {output.gff_refGene} && cut -f 12 {output.gff_refGene} > {output.col1} && cut -f 2-15 {output.gff_refGene} > {output.col_else} && paste {output.col1} {output.col_else} > {output.ref_refGene} && retrieve_seq_from_fasta.pl -format refGene -seqfile {input.INDEX} -outfile {output.ref_mrna} {output.ref_refGene} " 

rule annotation:
    input:
        avinput = "Call_variants/bwa_mem_snps_filtered.avinput",
        ref_refGene =  "annofile/ref-genome_new_refGene.txt",
        ref_mrna = "annofile/ref-genome_new_refGeneMrna.fa"
        
    output:
        files = expand("Call_variants/bwa_mem_snps_filtered.avinput.{IDX}",IDX=["exonic_variant_function", "variant_function"])
    log:
        "log/annotation.log"   
    threads: 1
    shell:
        "annotate_variation.pl --geneanno --dbtype refGene --buildver ref-genome_new {input.avinput} annofile/"
