# AlleleFlux

AlleleFlux is a toolkit for analyzing allele frequencies in metagenomic data.

## Project Structure

```
AlleleFlux/
├── scripts/                 # All Python scripts
│   ├── utilities/           # Utility functions
│   ├── preprocessing/       # Data preparation scripts
│   ├── analysis/            # Analysis scripts
│   └── statistics/          # Statistical scripts
├── smk_workflow/            # Snakemake workflow
│   ├── config.yml           # Configuration file
│   ├── step1.smk            # Workflow step 1
│   ├── step2.smk            # Workflow step 2
│   └── rules/               # Snakemake rules
```

## Installation

To install the package in development mode:

```bash
# Clone the repository
git clone <repository-url>
cd AlleleFlux

# Install in development mode
pip install -e .
```

## Usage

To run the workflow:

```bash
# Step 1: Profile samples and generate eligibility table
cd smk_workflow
snakemake -s step1.smk --profile cornell_profile/

# Step 2: Analyze alleles and calculate scores
cd smk_workflow
snakemake -s step2.smk --profile cornell_profile/
```

### Command-line Scripts

After installation, you can use these commands directly from the command line:

**Analysis scripts:**
- `alleleflux-profile` - Profile MAGs using alignment files
- `alleleflux-allele-freq` - Analyze allele frequencies
- `alleleflux-scores` - Calculate scores based on allele frequencies
- `alleleflux-taxa-scores` - Calculate taxonomic group scores
- `alleleflux-gene-scores` - Calculate gene-level scores
- `alleleflux-outliers` - Detect outlier genes

**Preprocessing scripts:**
- `alleleflux-metadata` - Generate MAG metadata
- `alleleflux-qc` - Perform quality control
- `alleleflux-eligibility` - Generate eligibility tables
- `alleleflux-preprocess-two-sample` - Preprocess data for two-sample analysis

**Statistics scripts:**
- `alleleflux-lmm` - Run linear mixed models
- `alleleflux-single-sample` - Perform single sample statistical tests
- `alleleflux-two-sample-paired` - Perform paired two-sample tests
- `alleleflux-two-sample-unpaired` - Perform unpaired two-sample tests

For help with any command, use the `-h` or `--help` flag, e.g.:
```bash
alleleflux-profile --help
```

## Configuration

Edit `smk_workflow/config.yml` to customize:
- Input/output paths
- Analysis parameters
- Resource requirements

# 🖱️ Conda Installation 🖱️

1. Install [`mamba`](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install) using the following commdands:

    ```bash
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3-$(uname)-$(uname -m).sh
    ```
2. Glone the git repository and install the required packages

    ```bash
    git clone https://github.com/MoellerLabPU/AlleleFlux.git
    cd AlleleFlux
    mamba env create -f environment.yml
    ```

3. Activate the environment

    ```bash
    mamba activate AlleleFlux
    ```

Now you are ready to run the workflow !

# 🐍 Snakemake workflow 🐍 

## 📂 Relevant Files 📂 

- `/AlleleFlux/smk_workflow/snakefile`: This is the main file that defines the Snakemake workflow. It contains the rules and dependencies for the analysis.
- `/AlleleFlux/smk_workflow/config.yml`: This file contains configuration parameters for the workflow, such as input file paths and analysis settings.
- `/AlleleFlux/smk_workflow/cornell_profile/config.yaml`: Specific run parameters to be submitted to SLURM. No need to edit it it unless you know what you're doing.

## 🏃 Running the worklow 🏃

1. You only really need to edit `/AlleleFlux/smk_workflow/config.yml` to be able to run the workflow. Change the paths to the scripts, file and any input parameters.
2. Make sure that the correct environment is activate ie. `AlleleFlux` and you're in `/AlleleFlux/smk_workflow` directory.
3. To run the workflow do, `snakemake --profile cornell_profile` and let the magic happen 🪄 👨‍🔬 👩‍🔬 !

## 📁 Input files format 📁 ##

- **bamDir**: Directory with sorted and indexed `bam` files. The files should have the name in the format `<sampleID>.sorted.bam` and `<sampleID>.sorted.bam.bai`.
- **fasta**: A big combined FASTA file of all the contigs making up the representative MAGs that are assumed to be present in the samples. The FASTA header should have the format `<MAG_ID>.fa_<contig_ID>`.
- **prodigal**: nucleic acid ORFs predictions by Prodigal of the above **fasta** file. NOTE that it's important that prodigal is run on the above **fasta** file for the IDs to match properly.
- **metadata_file**: Metadata file with the following columns: `sample_id` (should match the sampleID used in **bamDir**), `replicate`, `subjectID`, `time`, `group`.
- **gtdb_file**: Path to the `gtdbtk.bac120.summary.tsv` output file produced by GTDB-Tk.