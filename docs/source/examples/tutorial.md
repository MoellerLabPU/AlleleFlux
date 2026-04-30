# Tutorial

This tutorial uses the bundled synthetic dataset in
`docs/source/examples/example_data` to run AlleleFlux end-to-end on a normal
desktop computer. No HPC access or GPU is required.

---

## Dataset Overview

| Property | Value |
|---|---|
| MAGs | 2 (`MAG_001`, `MAG_002`) |
| Samples | 8 (4 control, 4 treatment; 2 subjects × 2 timepoints each) |
| Design | Longitudinal (pre → post) |
| Genome size | ~3,600–3,650 positions per MAG (2 contigs each) |

A simulated treatment effect is present in the pre-generated profiles: at ~8%
of positions in treatment post-timepoint samples, the dominant allele is shifted
to an alternative base. This signal is what AlleleFlux is designed to detect.

---

## Prerequisites

AlleleFlux installed (see [Installation](../getting_started/installation.md)):

```bash
conda activate alleleflux
```

Clone the repository and navigate to the example directory:

```bash
git clone https://github.com/MoellerLabPU/AlleleFlux.git
cd AlleleFlux/docs/source/examples/example_data
```

---

## Directory Structure

```
example_data/
├── config_example.yml              # Run from pre-generated profiles (Option A)
├── config_with_bams.yml            # Run from BAM files, includes profiling (Option B)
├── generate_bams.sh                # Regenerate BAMs from the reference FASTA
├── reference/
│   ├── combined_mags.fasta
│   ├── prodigal_genes.fna
│   ├── mag_mapping.tsv
│   └── gtdbtk_taxonomy.tsv
├── metadata/
│   ├── sample_metadata.tsv             # For Option A
│   └── sample_metadata_bams.tsv        # For Option B
├── profiles/                       # Pre-generated profiles with treatment effect
│   └── {sample}/{sample}_{MAG}_profiled.tsv.gz
├── bams/                           # Synthetic BAM files (~940 KB)
│   └── {sample}.bam[.bai]
└── significant_sites/
    └── significant_sites.tsv
```

---

## Option A — Pre-generated profiles (recommended)

The profiling step is skipped; pre-generated profiles are loaded directly. These
profiles were created with a treatment effect baked in (see
[Data generation details](#data-generation-details)), making the statistical
results meaningful even with a small dataset.

> **Caveat:** The profiles are synthetically generated, not derived from real
> BAM files. They are designed to demonstrate statistical detection, not to
> reflect real sequencing noise.

### Dry-run

```bash
alleleflux run --config config_example.yml --dry-run --threads 4
```

Expected: Snakemake prints 4 initial jobs (metadata, QC, eligibility, all) and
notes that checkpoint-based jobs will expand the DAG after running. No files
are created.

### Run end-to-end

```bash
alleleflux run --config config_example.yml --threads 4
```

Outputs land under `example_output/longitudinal/`.

### Inspect outputs

```bash
OUT=example_output/longitudinal

# QC eligibility table — both MAGs pass with 2 replicates/group
column -t -s $'\t' $OUT/eligibility_table_pre_post-treatment_control.tsv

# Allele frequency changes for MAG_001
zcat $OUT/allele_analysis/allele_analysis_pre_post-treatment_control/MAG_001_allele_frequency_changes_mean.tsv.gz | head -5

# Paired t-test results for MAG_001
zcat $OUT/significance_tests/two_sample_paired_pre_post-treatment_control/MAG_001_two_sample_paired.tsv.gz | head -5

# MAG-level parallel-evolution scores by phylum
column -t -s $'\t' \
  $OUT/scores/processed/combined/phylum/scores_two_sample_paired-pre_post-treatment_control-phylum.tsv

# dN/dS output directories (one per subject pair)
ls $OUT/dnds_analysis/pre_post-treatment_control/
```

---

## Option B — Full pipeline from BAM files

Runs the complete pipeline including the profiling step, starting from the
bundled synthetic BAM files.

> **Caveat:** The BAMs were generated with `wgsim` using the same parameters for
> all samples — there is **no treatment effect** at the read level. Any
> statistically significant sites in this run are artefacts of the small n rather
> than real biological signal. This option is provided to demonstrate that the
> profiling step runs correctly, not to demonstrate statistical detection.
>
> Additionally, with only 2 subjects per group, the paired longitudinal analysis
> yields n=2 paired differences per position. Non-parametric tests (Wilcoxon
> signed-rank, Mann-Whitney) require n≥5 to achieve p<0.05, so they will report
> zero significant sites. The paired t-test can produce sub-0.05 p-values at n=2
> but those results are statistically unreliable at this sample size.


### Run end-to-end

```bash
alleleflux run --config config_with_bams.yml --threads 4
```

Outputs land under `example_output_bams/longitudinal/`, including a `profiles/`
subdirectory with profiles generated from the BAMs.

---

## Expected Output Structure

```
longitudinal/
├── eligibility_table_pre_post-treatment_control.tsv   # Both MAGs: all tests eligible
├── QC/                                                 # Per-sample coverage stats
├── allele_analysis/                                    # Allele frequency changes per position
│   └── allele_analysis_pre_post-treatment_control/
├── significance_tests/                                 # Per-position p-values
│   ├── two_sample_paired_pre_post-treatment_control/
│   ├── two_sample_unpaired_pre_post-treatment_control/
│   └── single_sample_pre_post-treatment_control/
├── p_value_summary/                                    # FDR-corrected q-values
├── scores/                                             # Parallel-evolution scores by taxonomy
│   └── processed/combined/{phylum,genus,MAG}/
├── outlier_genes/                                      # Per-gene outlier scores
├── dnds_analysis/                                      # dN/dS per subject pair
│   └── pre_post-treatment_control/{subj1,subj2,subj3,subj4}/
└── preprocessing_eligibility/                          # Positions passing preprocessing filters
```

---

## Expected Run Time

Measured on a node pinned to 4 CPU cores:

| Option | Config | Wall time | Peak memory |
|---|---|---|---|
| A — Pre-generated profiles | `config_example.yml` | ~94 s | < 150 MB |
| B — Full pipeline from BAMs | `config_with_bams.yml` | ~74 s | < 150 MB |

Option B is faster in this demo because the wgsim-derived BAMs produce sparser
profiles (not every position achieves uniform coverage), so there is less data
in downstream steps. In a real experiment with whole-genome sequencing data,
profiling is the most computationally intensive step and will dominate runtime.


---

## Regenerating the BAM Files

To rebuild the BAMs from the reference FASTA:

```bash
bash generate_bams.sh
```
This uses the `reference/combined_mags.fasta` that is already committed in the
repository. `wgsim` and `bwa` must be installed (see above).

`wgsim` and `bwa` are **not** included in the `alleleflux` conda environment.
Install them if needed:

```bash
conda install -c bioconda bwa wgsim
```


## Data Generation Details

The pre-generated profiles were created with
[`generate_synthetic_data.py`](generate_synthetic_data.py) (seed = 42):

- Each MAG has 2 contigs (~800–2,000 bp, 3 genes per contig)
- Coverage drawn uniformly from 20–80× at every position
- The dominant allele receives 85–98% of reads; the remainder are split randomly
  among the other three bases

**Simulated treatment effect** (lines 364–366 of `generate_synthetic_data.py`):

```python
if is_treatment and is_post and random.random() < treatment_effect_rate:
    dominant_base = BASES[(pos % 4 + 1) % 4]
```

In treatment post-timepoint samples only, each position has an 8% probability of
shifting its dominant allele to the next base in the cycle (e.g., A → T). This
creates real but modest allele frequency differences between groups that AlleleFlux
is designed to detect.

The BAMs were generated with `wgsim` (100 bp paired-end, 0.5% error rate,
unique seed per sample) and aligned with `bwa mem`. All samples were generated
from the same reference with the same parameters, so no treatment effect is
present at the read level.


---

## File Format Reference

### Profile files (`*_profiled.tsv.gz`)

| Column | Type | Description |
|---|---|---|
| `contig` | string | Contig identifier |
| `position` | int | 0-based genomic position |
| `ref_base` | string | Reference base (A/C/G/T) |
| `total_coverage` | int | Total read depth |
| `A`, `C`, `G`, `T` | int | Per-base read counts |
| `N` | int | Ambiguous base count |
| `gene_id` | string | Overlapping gene ID (empty if intergenic) |

### Sample metadata

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier |
| `bam_path` | Path to BAM file or pre-generated profile directory |
| `subjectID` | Biological replicate / subject ID |
| `group` | Experimental group (`control` / `treatment`) |
| `replicate` | Replicate letter within group |
| `time` | Timepoint (`pre` / `post`) |
