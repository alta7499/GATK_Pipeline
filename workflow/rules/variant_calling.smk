# rules/variant_calling.smk
# Variant discovery: HaplotypeCaller (per-sample GVCF), CombineGVCFs, GenotypeGVCFs

rule haplotype_caller:
    """Call germline variants per-sample in GVCF mode."""
    input:
        bam = "results/preprocessing/{sample}_ready.bam",
        bai = "results/preprocessing/{sample}_ready.bai",
        ref = config["reference"]["fasta"]
    output:
        gvcf = protected("results/variants/{sample}.g.vcf.gz")
    resources:
        mem_mb = 8000,
        time   = "8:00:00"
    log:
        "logs/haplotype_caller/{sample}.log"
    benchmark:
        "benchmarks/haplotype_caller/{sample}.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "(gatk HaplotypeCaller "
        "   -R {input.ref} -I {input.bam} "
        "   -O {output.gvcf} -ERC GVCF) 2> {log}"

rule combine_gvcfs:
    """Merge per-sample GVCFs into a single cohort GVCF."""
    input:
        gvcfs = expand("results/variants/{sample}.g.vcf.gz", sample=SAMPLES),
        ref   = config["reference"]["fasta"]
    output:
        combined_gvcf = "results/variants/cohort.combined.g.vcf.gz"
    resources:
        mem_mb = 8000,
        time   = "4:00:00"
    log:
        "logs/combine_gvcfs/combine.log"
    benchmark:
        "benchmarks/combine_gvcfs/combine.tsv"
    conda:
        "../envs/gatk.yaml"
    params:
        vcf_args = lambda wildcards, input: " ".join([f"-V {v}" for v in input.gvcfs])
    shell:
        "(gatk CombineGVCFs "
        "   -R {input.ref} "
        "   {params.vcf_args} "
        "   -O {output.combined_gvcf}) 2> {log}"

rule genotype_gvcfs:
    """Joint genotyping across cohort to produce raw variant calls."""
    input:
        gvcf = "results/variants/cohort.combined.g.vcf.gz",
        ref  = config["reference"]["fasta"]
    output:
        vcf = "results/variants/cohort_raw.vcf.gz"
    resources:
        mem_mb = 8000,
        time   = "4:00:00"
    log:
        "logs/genotype_gvcfs/genotype.log"
    benchmark:
        "benchmarks/genotype_gvcfs/genotype.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "(gatk GenotypeGVCFs "
        "   -R {input.ref} -V {input.gvcf} "
        "   -O {output.vcf}) 2> {log}"
