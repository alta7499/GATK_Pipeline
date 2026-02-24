# rules/variant_filtering.smk
# Hard filtering for SNPs and Indels following GATK Best Practices

rule select_snps:
    """Extract SNP variants from raw joint calls."""
    input:
        vcf = "results/variants/cohort_raw.vcf.gz",
        ref = config["reference"]["fasta"]
    output:
        vcf = temp("results/variants/cohort_snps_raw.vcf.gz")
    log:
        "logs/select_snps/select_snps.log"
    benchmark:
        "benchmarks/select_snps/select_snps.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "(gatk SelectVariants "
        "   -R {input.ref} -V {input.vcf} "
        "   -select-type SNP "
        "   -O {output.vcf}) 2> {log}"

rule filter_snps:
    """Apply GATK hard filters to SNPs (GATK Best Practices thresholds)."""
    input:
        vcf = "results/variants/cohort_snps_raw.vcf.gz",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/variants/cohort_snps_filtered.vcf.gz"
    log:
        "logs/filter_snps/filter_snps.log"
    benchmark:
        "benchmarks/filter_snps/filter_snps.tsv"
    conda:
        "../envs/gatk.yaml"
    params:
        filter_expr = "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0",
        filter_name = "snp_hard_filter"
    shell:
        "(gatk VariantFiltration "
        "   -R {input.ref} -V {input.vcf} "
        "   --filter-expression '{params.filter_expr}' "
        "   --filter-name '{params.filter_name}' "
        "   -O {output.vcf}) 2> {log}"

rule select_indels:
    """Extract Indel variants from raw joint calls."""
    input:
        vcf = "results/variants/cohort_raw.vcf.gz",
        ref = config["reference"]["fasta"]
    output:
        vcf = temp("results/variants/cohort_indels_raw.vcf.gz")
    log:
        "logs/select_indels/select_indels.log"
    benchmark:
        "benchmarks/select_indels/select_indels.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "(gatk SelectVariants "
        "   -R {input.ref} -V {input.vcf} "
        "   -select-type INDEL "
        "   -O {output.vcf}) 2> {log}"

rule filter_indels:
    """Apply GATK hard filters to Indels (GATK Best Practices thresholds)."""
    input:
        vcf = "results/variants/cohort_indels_raw.vcf.gz",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/variants/cohort_indels_filtered.vcf.gz"
    log:
        "logs/filter_indels/filter_indels.log"
    benchmark:
        "benchmarks/filter_indels/filter_indels.tsv"
    conda:
        "../envs/gatk.yaml"
    params:
        filter_expr = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0",
        filter_name = "indel_hard_filter"
    shell:
        "(gatk VariantFiltration "
        "   -R {input.ref} -V {input.vcf} "
        "   --filter-expression '{params.filter_expr}' "
        "   --filter-name '{params.filter_name}' "
        "   -O {output.vcf}) 2> {log}"

rule merge_filtered_vcfs:
    """Merge filtered SNPs and Indels into a final cohort VCF."""
    input:
        snps   = "results/variants/cohort_snps_filtered.vcf.gz",
        indels = "results/variants/cohort_indels_filtered.vcf.gz",
        ref    = config["reference"]["fasta"]
    output:
        vcf = protected("results/variants/cohort_filtered_final.vcf.gz")
    log:
        "logs/merge_filtered/merge.log"
    benchmark:
        "benchmarks/merge_filtered/merge.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "(gatk MergeVcfs "
        "   -I {input.snps} "
        "   -I {input.indels} "
        "   -O {output.vcf}) 2> {log}"
