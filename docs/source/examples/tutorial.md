# Tutorial

This walkthrough uses the bundled example dataset in `docs/source/examples/example_data` plus the ready-made `config_example.yml`. It shows how to dry-run, execute, and inspect the outputs for a single MAG/timepoint/group. Micro-examples at the end exercise individual tools on the tiny mock datasets in `tests/evolution/mock_data`.

## Prerequisites

- AlleleFlux installed and on `$PATH` (see [Installation](../getting_started/installation.md)).
- Working from the repo root (`AlleleFlux/`) so relative paths resolve.

## Step 1: Start from the template config

Option A: copy the example config (pre-populated for the bundled data):

```bash
cp docs/source/examples/example_data/config_example.yml ./config_example.yml
```

Option B: print the template then edit:

```bash
alleleflux init --template > config_example.yml
```

Open `config_example.yml` and confirm the paths point at `docs/source/examples/example_data`. The important bits (already set in the example):

```yaml
input:
  fasta_path: docs/source/examples/example_data/reference/combined_mags.fasta
  prodigal_path: docs/source/examples/example_data/reference/prodigal_genes.fna
  metadata_path: docs/source/examples/example_data/metadata/sample_metadata.tsv
  gtdb_path: docs/source/examples/example_data/reference/gtdbtk_taxonomy.tsv
  mag_mapping_path: docs/source/examples/example_data/reference/mag_mapping.tsv
output:
  root_dir: ./example_output
analysis:
  data_type: longitudinal
  timepoints_combinations:
    - timepoint: ["pre", "post"]
      focus: "post"
  groups_combinations:
    - ["treatment", "control"]
```

## Step 2: Dry-run the workflow

```bash
alleleflux run --config config_example.yml --dry-run
```

This builds the DAG and confirms inputs without running jobs.

## Step 3: Run the example end-to-end

```bash
# modest local resources to keep the run quick
alleleflux run --config config_example.yml --threads 4 --memory 8G
```

Outputs land under `example_output/longitudinal` (the pipeline appends `data_type`).

## Step 4: Inspect one MAG/timepoint/group

The example uses the label `pre_post` for timepoints and `treatment_control` for groups. Here is a quick inspection for `TEST_MAG_001`:

```bash
OUT=example_output/longitudinal

# QC eligibility table that drives downstream targets
column -t -s $'\t' $OUT/eligibility_table_pre_post-treatment_control.tsv | head

# Allele analysis outputs (mean changes and zero-diff filtered)
ls $OUT/allele_analysis/allele_analysis_pre_post-treatment_control \
   | grep TEST_MAG_001
zcat $OUT/allele_analysis/allele_analysis_pre_post-treatment_control/TEST_MAG_001_allele_frequency_changes_mean.tsv.gz \
     | head

# Significance tests for this MAG/timepoint/group
ls $OUT/significance_tests/two_sample_paired_pre_post-treatment_control \
   | grep TEST_MAG_001

# Combined MAG-level scores
column -t -s $'\t' \
  $OUT/scores/processed/combined/MAG/scores_two_sample_paired-pre_post-treatment_control-MAGs.tsv \
  | head

# dN/dS per subject (longitudinal only)
ls $OUT/dnds_analysis/pre_post-treatment_control
```

## Step 5: Micro-examples with mock data

These commands use the tiny bundles in `tests/evolution/mock_data` so you can exercise individual tools without running the full workflow.

### `alleleflux-preprocess-between-groups` (paired filter on mock mean changes)

Create a minimal mean-changes table from the bundled pre/post profiles for `MAG_001` and run the filter:

```bash
python - <<'PY'
import gzip, pandas as pd

pre = pd.read_csv("tests/evolution/mock_data/mock_dnds_1/PipelineMock/profiles/pre_sample/pre_sample_MAG_001_profiled.tsv.gz", sep="\t")
post = pd.read_csv("tests/evolution/mock_data/mock_dnds_1/PipelineMock/profiles/post_sample/post_sample_MAG_001_profiled.tsv.gz", sep="\t")

# compute per-position mean change for each nucleotide (post - pre) normalized by coverage
for df, label in [(pre, "pre"), (post, "post")]:
    for base in ["A", "T", "G", "C"]:
        df[f"{base}_frequency"] = df[base] / df["total_coverage"].clip(lower=1)
    df["timepoint"] = label

merged = pre.merge(
    post,
    on=["contig", "position", "gene_id", "ref_base"],
    suffixes=("_pre", "_post"),
)
out_cols = ["contig", "position", "gene_id"]
for base in ["A", "T", "G", "C"]:
    merged[f"{base}_frequency"] = merged[f"{base}_frequency_post"] - merged[f"{base}_frequency_pre"]
    out_cols.append(f"{base}_frequency")

merged[out_cols].to_csv("/tmp/mock_mean_changes.tsv", sep="\t", index=False)
print("Wrote /tmp/mock_mean_changes.tsv with", len(merged), "rows")
PY

alleleflux-preprocess-between-groups \
  --mean_changes_fPath /tmp/mock_mean_changes.tsv \
  --output_fPath /tmp/mock_mean_changes_preprocessed.tsv \
  --p_value_threshold 0.05 \
  --data_type longitudinal \
  --filter_type t-test \
  --mag_id MAG_001 \
  --min_positions 1 \
  --min_sample_num 1 \
  --status_dir /tmp
```

Outputs: `/tmp/mock_mean_changes_preprocessed.tsv` (filtered sites) and `/tmp/MAG_001_preprocessing_status.json` with eligibility counts.

### `alleleflux-dnds-from-timepoints` (mock significant sites → dN/dS)

Use the mock significant sites and profiles shipped with the tests:

```bash
alleleflux-dnds-from-timepoints \
  --input tests/evolution/mock_data/mock_dnds_1/PipelineMock/significant_sites.tsv \
  --output /tmp/mock_dnds \
  --mag_id MAG_001 \
  --profiles_dir tests/evolution/mock_data/mock_dnds_1/PipelineMock/profiles \
  --prodigal_fasta tests/evolution/mock_data/mock_dnds_1/PipelineMock/prodigal_genes.fasta \
  --fasta tests/evolution/mock_data/mock_dnds_1/PipelineMock/prodigal_genes.fasta \
  --p_value_column q_value \
  --p_value_threshold 0.05 \
  --test_type two_sample_paired_tTest
```

Outputs (under `/tmp/mock_dnds`): codon, gene, MAG, and global NG86 summaries ready to inspect or plot.
