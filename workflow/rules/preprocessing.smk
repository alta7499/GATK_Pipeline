# rules/preprocessing.smk
# Pre-processing rules: QC trimming, alignment, deduplication, and BQSR

rule fastp:
    """Quality control and adapter trimming of raw FASTQ reads."""
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["r2"]
    output:
        r1   = "results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2   = "results/trimmed/{sample}_R2.trimmed.fastq.gz",
        json = "results/trimmed/{sample}.fastp.json",
        html = "results/trimmed/{sample}.fastp.html"
    threads: 4
    resources:
        mem_mb = 4000,
        time   = "1:00:00"
    log:
        "logs/fastp/{sample}.log"
    benchmark:
        "benchmarks/fastp/{sample}.tsv"
    conda:
        "../envs/qc.yaml"
    shell:
        "(fastp "
        "   -i {input.r1} -I {input.r2} "
        "   -o {output.r1} -O {output.r2} "
        "   -j {output.json} -h {output.html} "
        "   --thread {threads}) 2> {log}"

rule align:
    """Align trimmed reads to the reference genome and sort the output BAM."""
    input:
        r1  = "results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2  = "results/trimmed/{sample}_R2.trimmed.fastq.gz",
        ref = config["reference"]["fasta"]
    output:
        bam = temp("results/aligned/{sample}.bam")
    threads: config["threads"]["bwa"]
    resources:
        mem_mb = 16000,
        time   = "8:00:00"
    log:
        "logs/align/{sample}.log"
    benchmark:
        "benchmarks/align/{sample}.tsv"
    conda:
        "../envs/align.yaml"
    shell:
        "(bwa-mem2 mem -t {threads} "
        "   -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA' "
        "   {input.ref} {input.r1} {input.r2} "
        "   | samtools sort -@ {threads} -o {output.bam} -) 2> {log}"

rule mark_duplicates:
    """Mark PCR and optical duplicates with GATK MarkDuplicates."""
    input:
        bam = "results/aligned/{sample}.bam"
    output:
        bam     = temp("results/preprocessing/{sample}_md.bam"),
        metrics = "results/preprocessing/{sample}_md_metrics.txt"
    threads: config["threads"]["gatk"]
    resources:
        mem_mb = 8000,
        time   = "4:00:00"
    log:
        "logs/mark_duplicates/{sample}.log"
    benchmark:
        "benchmarks/mark_duplicates/{sample}.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "(gatk --java-options '-Xmx{resources.mem_mb}m' MarkDuplicates "
        "   -I {input.bam} "
        "   -O {output.bam} "
        "   -M {output.metrics}) 2> {log}"

rule base_recalibrator:
    """Build a base quality score recalibration (BQSR) table."""
    input:
        bam          = "results/preprocessing/{sample}_md.bam",
        ref          = config["reference"]["fasta"],
        known_sites  = config["reference"]["known_sites"]
    output:
        recal_table = "results/preprocessing/{sample}_recal.table"
    resources:
        mem_mb = 8000,
        time   = "4:00:00"
    log:
        "logs/base_recalibrator/{sample}.log"
    benchmark:
        "benchmarks/base_recalibrator/{sample}.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "(gatk BaseRecalibrator "
        "   -I {input.bam} -R {input.ref} "
        "   --known-sites {input.known_sites} "
        "   -O {output.recal_table}) 2> {log}"

rule apply_bqsr:
    """Apply BQSR to produce the analysis-ready BAM."""
    input:
        bam         = "results/preprocessing/{sample}_md.bam",
        recal_table = "results/preprocessing/{sample}_recal.table",
        ref         = config["reference"]["fasta"]
    output:
        bam = protected("results/preprocessing/{sample}_ready.bam"),
        bai = "results/preprocessing/{sample}_ready.bai"
    resources:
        mem_mb = 8000,
        time   = "4:00:00"
    log:
        "logs/apply_bqsr/{sample}.log"
    benchmark:
        "benchmarks/apply_bqsr/{sample}.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "(gatk ApplyBQSR "
        "   -I {input.bam} -R {input.ref} "
        "   --bqsr-recal-file {input.recal_table} "
        "   -O {output.bam} && "
        "samtools index {output.bam}) 2> {log}"
