# GATK Snakemake Pipeline

A Snakemake pipeline following GATK Best Practices for Germline Short Variant Discovery.

## Installation

1. **Install Mamba/Conda** (if not already installed)
   We recommend using [Miniforge](https://github.com/conda-forge/miniforge) which includes `mamba`.

2. **Create the Snakemake environment**
   This environment contains Snakemake and the necessary plugins to execute the pipeline.
   ```bash
   mamba env create -f environment.yaml
   ```

3. **Activate the environment**
   ```bash
   conda activate gatk-snakemake
   ```

## Usage

1. Place your raw FASTQ files in `data/fastq/`.
2. Provide your reference genome in `resources/reference/`.
3. Update `config/config.yaml` with your sample names and file paths under the `samples:` key.
4. Run the pipeline (using conda to automatically manage tool dependencies):
   ```bash
   # Dry-run
   snakemake -n --use-conda

   # Run locally (e.g., 16 cores)
   snakemake --cores 16 --use-conda
   ```
