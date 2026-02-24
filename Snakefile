# =============================================================================
# GATK Germline Short Variant Discovery Pipeline
# Best Practices: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932
# Snakemake: https://snakemake.readthedocs.io/
# =============================================================================

configfile: "config/config.yaml"

# Derive sample list directly from config
SAMPLES = list(config["samples"].keys())

# ---------------------------------------------------------------------------
# Include modular rule files
# ---------------------------------------------------------------------------
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/variant_calling.smk"
include: "workflow/rules/variant_filtering.smk"

# ---------------------------------------------------------------------------
# Target rule: defines final outputs of the full pipeline
# ---------------------------------------------------------------------------
rule all:
    input:
        # Per-sample QC reports
        expand("results/trimmed/{sample}.fastp.json", sample=SAMPLES),
        # Per-sample analysis-ready BAMs
        expand("results/preprocessing/{sample}_ready.bam", sample=SAMPLES),
        # Per-sample GVCFs
        expand("results/variants/{sample}.g.vcf.gz", sample=SAMPLES),
        # Final joint-filtered VCF
        "results/variants/cohort_filtered_final.vcf.gz"
