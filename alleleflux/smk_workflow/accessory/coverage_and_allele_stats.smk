"""
Coverage and Allele Statistics Calculation

Simplified Snakemake module to compute coverage/allele statistics per MAG with
low memory usage. Inputs are:
- A text file with MAG IDs (one per line), typically produced by alleleflux-list-mags
- rootDir: directory containing the per-MAG *_metadata.tsv files
- qc_dir (optional): directory with *_QC.tsv files

Each MAG is processed once as an independent job. Parallelism is handled by
Snakemake across MAGs; the Python script processes a single MAG using --mag_id.
"""

# Python imports used in rule definitions
import os

with open(config["mag_list"]) as f:
    MAGS = [line.strip() for line in f if line.strip()]

OUTDIR=config["output_dir"]

rule all:
    input:
        expand(os.path.join(OUTDIR, "{mag}_mean_coverage.tsv"), mag=MAGS)

rule compute_mean_coverage:
    input:
        metadata_dir=config["metadata_dir"],
    output:
        # Define the output file that this rule will create.
        mean_coverage_out=os.path.join(OUTDIR, "{mag}_mean_coverage.tsv"),
    params:
        # Conditionally add the --qc_dir argument. This checks if 'qc_dir' was provided in the config. 
        # If yes, it formats the argument string; otherwise, it returns an empty string.
        qc_option=(
            f"--qc_dir {config['qc_dir']}"
            if config.get("qc_dir")
            else ""
        )
    threads: 1
    shell:
        """
        alleleflux-coverage-allele-stats \
            --rootDir {input.metadata_dir} \
            --output_dir {OUTDIR} \
            --mag_id {wildcards.mag} \
            {params.qc_option}
        """