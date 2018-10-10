# trimming reads
rule fastqc_init:
    input:
        get_fastq
    output:
        html="results/0_fastqc_init/{sample}-{unit}.html",
        zip="results/0_fastqc_init/{sample}-{unit}.zip"
    wrapper:
        "0.27.1/bio/fastqc"

rule trim_reads_se:
    input:
        get_fastq
    output:
        fastq="results/1_trimmed/{sample}-{unit}.fastq.gz",
        qc="results/1_trimmed/{sample}-{unit}.qc.txt"
    params:
        ampl_to_cutadapt_se, config["params"]["cutadapt"]["se"]
    log:
        "results/logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.27.1/bio/cutadapt/se"


rule trim_reads_pe:
    input:
        get_fastq
    output:
        fastq1="results/1_trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/1_trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/1_trimmed/{sample}-{unit}.qc.txt"
    params:
        ampl_to_cutadapt_pe, config["params"]["cutadapt"]["pe"]
    log:
        "results/logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.27.1/bio/cutadapt/pe"

rule fastqc_trim:
    input:
        get_trimmed_reads
    output:
        html="results/2_fastqc_trim/{sample}-{unit}.html",
        zip="results/2_fastqc_trim/{sample}-{unit}.zip"
    wrapper:
        "0.27.1/bio/fastqc"

# genome preparation and alignment
rule bwa_index:
    input:
        config["genome"]
    output:
        config["genome"] + ".amb",
        config["genome"] + ".ann",
        config["genome"] + ".bwt",
        config["genome"] + ".pac",
        config["genome"] + ".sa"
    log:
        expand("results/logs/bwa_index/{genome}.log", genome=config["genome"])
    wrapper:
        "0.27.1/bio/bwa/index"

rule samtools_faidx:
    input:
        config["genome"]
    output:
        config["genome"] + '.fai'
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input}"

rule map_reads_bwa_mem:
    input:
        reads=get_trimmed_reads
    output:
        "results/3_mapped/{sample}-{unit}.sorted.bam"
    log:
        "results/logs/bwa_mem/{sample}-{unit}.log"
    params:
        index=config["genome"],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: config["params"]["threads"]
    wrapper:
        "0.27.1/bio/bwa/mem"

# alignment filtering
rule mark_duplicates:
    input:
        "results/3_mapped/{sample}-{unit}.sorted.bam"
    output:
        bam=temp("results/4_dedup/{sample}-{unit}.bam"),
        metrics="results/4_dedup/{sample}-{unit}.dedup.txt"
    log:
        "results/logs/picard/dedup/{sample}-{unit}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.27.1/bio/picard/markduplicates"


rule samptools_filter:
    input:
        # get_dedup_bams
        "results/4_dedup/{sample}-{unit}.sorted.bam"
    output:
        "results/5_filtered/{sample}-{unit}.bam"
        #metrics="results/5_filtered/{sample}-{unit}.filter.txt"
    params:
        minq=config["params"]["filter"]["min_mapq"],
        minl=config["params"]["filter"]["min_len"]
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -h -q {params.minq} {input} | "
        "awk 'length($10) > {params.minl} || $1 ~ /^@/' | "
        "samtools view -bS - > {output}"

rule samtools_merge:
    input:
        # get_filtered_bams
        "results/5_filtered/{sample}-{unit}.bam"
    output:
        "results/6_merged/{sample}.bam"
    params:
        "" 
    threads: config["params"]["threads"]
    wrapper:
        "0.27.1/bio/samtools/merge"



# regions
rule regions:
    input:
        "results/6_merged/{sample}.bam",
        genome_fai=config["genome"] + '.fai'
    output:
        pos="results/7_positions/{sample}.bed",
        reg="results/7_regions/{sample}.tsv"
    params:
        sample="{sample}"
    conda:
        "../envs/regions.yaml"
    script:
        "../scripts/regions.py"

# statistics
# rule stats:
#     input:
#         "results/6_filtered/{sample}.bam"
#     output:
#         "results/stats.tsv"
#     conda:
#         "../envs/stats.yaml"
#     script:
#         "../scripts/stats.py"