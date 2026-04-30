#!/usr/bin/env bash
# generate_bams.sh
# ----------------
# Regenerates the synthetic BAM files in bams/ from the bundled reference FASTA.
#
# What this script does
# ---------------------
# For each of the 8 samples it runs a three-step pipeline:
#
#   1. wgsim  — simulates paired-end Illumina reads from the reference FASTA
#   2. bwa mem — aligns those reads back to the reference
#   3. samtools sort + index — produces a sorted, indexed BAM
#
# Why separate seeds per sample
# ------------------------------
# All 8 samples are generated from the identical reference genome with
# identical parameters, but each uses a different random seed. This means
# the reads differ slightly between samples (different positions get covered,
# slightly different noise) but there is NO systematic treatment effect —
# treatment and control samples are statistically indistinguishable.
#
# To add a real treatment effect you would need to modify the reference FASTA
# itself before simulating (e.g. introduce SNPs at 8% of positions) and use
# the modified reference only for treatment post-timepoint samples.
# generate_synthetic_data.py takes the alternative approach: it injects the
# treatment effect directly into the allele count tables, bypassing reads
# entirely.
#
# Requirements
# ------------
# wgsim and bwa are NOT included in the alleleflux conda environment.
# Install them first:
#   conda install -c bioconda bwa wgsim
# samtools IS included in the alleleflux environment.
#
# Usage
# -----
#   cd docs/source/examples/example_data
#   bash generate_bams.sh
#
# Approximate runtime: < 5 seconds for the bundled 7 kb reference

set -euo pipefail   # exit on any error, treat unset vars as errors, propagate pipe failures

# Resolve the directory this script lives in regardless of where it is called from.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$SCRIPT_DIR/reference/combined_mags.fasta"
OUTDIR="$SCRIPT_DIR/bams"

echo "[generate_bams] Working in: $SCRIPT_DIR"

# ---- Dependency check --------------------------------------------------------
for cmd in wgsim bwa samtools; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: '$cmd' not found."
        echo "Install with: conda install -c bioconda bwa wgsim"
        exit 1
    fi
done

mkdir -p "$OUTDIR"

# ---- Index the reference (skip if already done) ------------------------------
# bwa needs 5 index files (.amb .ann .bwt .pac .sa) before it can align.
# We check for .bwt as a proxy for all five.
if [[ ! -f "$REF.bwt" ]]; then
    echo "[generate_bams] Indexing reference with bwa..."
    bwa index "$REF" 2>/dev/null
fi

# samtools faidx is needed to read contig lengths for the coverage calculation.
if [[ ! -f "$REF.fai" ]]; then
    samtools faidx "$REF"
fi

# ---- Compute number of read pairs for ~50× coverage -------------------------
# Formula: N_pairs = (target_coverage × genome_size) / (2 × read_length)
#          = (50 × TOTAL_BP) / (2 × 100)  = TOTAL_BP / 4
#
# Example: reference is ~7,173 bp → N_pairs = 7173 / 4 ≈ 1,793 pairs per sample
# At 100 bp per read: 1793 × 2 × 100 = 358,600 bp → 358,600 / 7,173 ≈ 50×
TOTAL_BP=$(awk '{s+=$2} END {print s}' "$REF.fai")
N_PAIRS=$(( (50 * TOTAL_BP) / 200 ))
echo "[generate_bams] Reference: ${TOTAL_BP} bp | Read pairs per sample: ${N_PAIRS} (~50×)"

# ---- Sample → seed mapping ---------------------------------------------------
# Using a fixed seed per sample makes the BAMs fully reproducible. Different
# seeds produce different random read sets (different coverage fluctuations,
# different sequencing errors) while keeping the same underlying genome.
declare -A SEEDS=(
    [control_subj1_pre]=1
    [control_subj1_post]=2
    [control_subj2_pre]=3
    [control_subj2_post]=4
    [treatment_subj3_pre]=5
    [treatment_subj3_post]=6
    [treatment_subj4_pre]=7
    [treatment_subj4_post]=8
)

# Use a temp directory for intermediate FASTQ files; clean up automatically on exit.
TMPDIR_LOCAL=$(mktemp -d)
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

# ---- Main loop: simulate → align → sort → index -----------------------------
for SAMPLE in "${!SEEDS[@]}"; do
    SEED="${SEEDS[$SAMPLE]}"
    BAM="$OUTDIR/${SAMPLE}.bam"
    echo "[generate_bams] $SAMPLE (seed=$SEED)"

    # Step 1 — wgsim: simulate paired-end reads from the reference
    #   -e 0.005  base error rate 0.5% (controls substitution frequency)
    #   -N        number of read pairs
    #   -1 100    read 1 length = 100 bp
    #   -2 100    read 2 length = 100 bp
    #   -r 0.001  SNP mutation rate (additional variants beyond the reference)
    #   -R 0.0    indel fraction = 0 (no indels, keeps alignment simple)
    #   -X 0.0    indel extension probability = 0
    #   -S        random seed (reproducibility)
    #
    # Output: r1.fq and r2.fq in the temp directory.
    # Note: wgsim assigns Phred quality ~23 to all bases regardless of error rate.
    # This is below the default alleleflux-profile threshold of 30, which is why
    # config_with_bams.yml sets min_base_quality: 0 for this synthetic dataset.
    wgsim \
        -e 0.005 \
        -N "$N_PAIRS" \
        -1 100 -2 100 \
        -r 0.001 -R 0.0 -X 0.0 \
        -S "$SEED" \
        "$REF" \
        "$TMPDIR_LOCAL/r1.fq" "$TMPDIR_LOCAL/r2.fq" \
        > /dev/null 2>&1

    # Step 2 — bwa mem: align reads to the reference
    #   -t 4          use 4 threads
    #   -R "@RG\t..."  read-group tag; samtools and many downstream tools expect
    #                  SM (sample name) to be set for per-sample identification
    #
    # Step 3 — samtools sort: convert SAM → coordinate-sorted BAM
    #   -@4  use 4 threads for sorting
    #   -o   output file path
    #   -    read SAM from stdin (piped directly from bwa mem)
    bwa mem -t 4 \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
        "$REF" "$TMPDIR_LOCAL/r1.fq" "$TMPDIR_LOCAL/r2.fq" 2>/dev/null \
    | samtools sort -@4 -o "$BAM" -

    # Step 4 — samtools index: create .bai index so alleleflux-profile can seek
    # to arbitrary genomic positions without reading the whole BAM.
    samtools index "$BAM"
done

echo "[generate_bams] Done. BAMs written to: $OUTDIR"
echo "[generate_bams] Run the full pipeline with:"
echo "  alleleflux run --config config_with_bams.yml --threads 4"
