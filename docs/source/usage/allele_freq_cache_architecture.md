# Allele Frequency Cache Architecture

This page explains how AlleleFlux computes allele frequencies in two stages to avoid
redundant work when multiple analysis combinations share the same timepoints. It covers
the group-independent cache design, how group filtering is applied downstream, and how
the pipeline selects the right QC file for each cache job.

---

## Overview

AlleleFlux splits allele frequency processing into two stages:

```
Stage 1: alleleflux-cache-allele-freq         one job per (MAG, timepoint)
         Reads: profile TSVs for all samples at this timepoint (all groups)
         Writes: allele_freq_cache/{timepoint}/{mag}_{timepoint}_allele_frequency.parquet

Stage 2: alleleflux-allele-freq               one job per (MAG, tp_combo, gr_combo)
         Reads: one or two Stage-1 Parquet cache files
         Filters: to the two groups in this combination via --groups
         Writes: allele_analysis/allele_analysis_{tp_combo}-{gr_combo}/{mag}_...parquet
```

The key design decision: **the cache is group-independent**. A single Parquet file per
`(MAG, timepoint)` holds every sample at that timepoint, regardless of group. Group
filtering happens in Stage 2.

---

## Why the Cache Is Group-Independent

### What "group-independent" means

Different group combinations (e.g. `1D_AL` vs `2D_AL`) compare different pairs of groups,
but they draw from the same underlying sample data. A mouse in group `1D` has identical
allele frequencies whether it will be compared to `AL` or to `2D`. The per-sample frequency
values are properties of the sample, not of which comparison it participates in.

**What varies per group combination is which rows you use**, not what the rows contain.
The cache therefore holds all rows for a timepoint, and the group split happens downstream
where it is cheap (a single `df[df["group"].isin(...)]` filter) rather than upstream where
it would require re-reading gigabytes of profile TSVs.

### DRIDO example

The active DRIDO config (`alleleflux_config_1.yml`) defines:

| Config dimension | Values |
|---|---|
| Timepoint combinations | `5mo_10mo`, `5mo_16mo`, `5mo_22mo` |
| Group combinations | `1D_AL`, `2D_AL`, `20_AL`, `40_AL`, `40_20`, `2D_1D` |
| Unique timepoints | `5mo`, `10mo`, `16mo`, `22mo` |
| Groups | `1D`, `2D`, `20`, `40`, `AL` |

**Cache-write jobs (per MAG):**

| Design | Formula | Jobs per MAG |
|---|---|---|
| Old (group-scoped) | 4 timepoints × 6 gr_combos | 24 |
| Current (group-independent) | 4 timepoints × 1 | **4** |

For ~274 DRIDO MAGs: 24 × 274 = 6,576 cache jobs → **4 × 274 = 1,096**.

The four cache files per MAG hold data for *all* groups:

```
allele_freq_cache/5mo/MRGM_1079_5mo_allele_frequency.parquet    ← groups: 1D, 2D, 20, 40, AL
allele_freq_cache/10mo/MRGM_1079_10mo_allele_frequency.parquet  ← same
allele_freq_cache/16mo/MRGM_1079_16mo_allele_frequency.parquet
allele_freq_cache/22mo/MRGM_1079_22mo_allele_frequency.parquet
```

All 6 group combinations read from these same 4 files.

---

## Stage 1: Writing the Cache

**Script:** `alleleflux-cache-allele-freq`
([`scripts/analysis/allele_frequency/allele_freq_cache.py`](../../../alleleflux/scripts/analysis/allele_frequency/allele_freq_cache.py))

**Invoked once per `(MAG, timepoint)`.** There is no groups argument.

### Example: `5mo` cache job for `MRGM_1079`

```bash
alleleflux-cache-allele-freq \
    --magID       MRGM_1079 \
    --qc_file     {OUTDIR}/QC/QC_{timepoint}/{mag}_QC.tsv \
    --timepoint   5mo \
    --data_type   longitudinal \
    --output_path {OUTDIR}/allele_freq_cache/5mo/MRGM_1079_5mo_allele_frequency.parquet \
    --cpus        8
```

### What the script does

1. Loads the MAG's QC table of the requested timepoint,
   `QC/QC_{timepoint}/{mag}_QC.tsv`. The table lists every sample of that one timepoint,
   across **all groups**.
2. Keeps the samples that passed QC (`coverage_threshold_passed`). Every row already
   belongs to the requested timepoint; the filter on `time` stays as a guard for a call by
   hand with a table that holds several timepoints.
3. Reads each sample's profile TSV in parallel (`--cpus` workers), computing per-base
   allele frequencies (`A_frequency`, `T_frequency`, `G_frequency`, `C_frequency`).
4. Concatenates all per-sample DataFrames and writes a single Parquet file.

The resulting Parquet has one row per `(contig, position, sample)` with the four base
counts, `total_coverage`, the derived `*_frequency` columns, `gene_id`, `ref_base` (the
reference genome's base at that position, kept for the dN/dS sanity check) and the
metadata columns `group`, `subjectID`, `replicate`, `time`.

### Why one QC file covers all groups

QC is built once per **single timepoint** and does not depend on groups or on comparisons.
Each QC job writes one table per MAG into `QC/QC_{timepoint}/`; the path has no comparison
and no groups component. The table lists the samples of every group at that timepoint, so
the Stage 1 cache job reads it directly to get its sample list.

---

## Stage 2: Filtering to a Group Combination

**Script:** `alleleflux-allele-freq`
([`scripts/analysis/allele_frequency/allele_freq.py`](../../../alleleflux/scripts/analysis/allele_frequency/allele_freq.py))

**Invoked once per `(MAG, tp_combo, gr_combo)`.** Receives `--groups` to filter the cache.

### Example: two runs on the same cache files

Both jobs below read the **same two Parquet cache files** but produce different outputs
because `--groups` selects different rows.

```bash
# Combination: 5mo_10mo vs 1D_AL
alleleflux-allele-freq \
    --magID       MRGM_1079 \
    --cache_files {OUTDIR}/allele_freq_cache/5mo/MRGM_1079_5mo_allele_frequency.parquet \
                  {OUTDIR}/allele_freq_cache/10mo/MRGM_1079_10mo_allele_frequency.parquet \
    --data_type   longitudinal \
    --groups      1D AL \
    --output_dir  {OUTDIR}/allele_analysis/allele_analysis_5mo_10mo-1D_AL/

# Combination: 5mo_10mo vs 2D_AL (same cache files, different --groups)
alleleflux-allele-freq \
    --magID       MRGM_1079 \
    --cache_files {OUTDIR}/allele_freq_cache/5mo/MRGM_1079_5mo_allele_frequency.parquet \
                  {OUTDIR}/allele_freq_cache/10mo/MRGM_1079_10mo_allele_frequency.parquet \
    --data_type   longitudinal \
    --groups      2D AL \
    --output_dir  {OUTDIR}/allele_analysis/allele_analysis_5mo_10mo-2D_AL/
```

### How `--groups` is applied

```python
# allele_freq.py (Stage 2) — immediately after loading the cache
if args.groups:
    allele_df["group"] = allele_df["group"].astype(str)
    allele_df = allele_df[allele_df["group"].isin(args.groups)].reset_index(drop=True)
    logger.info(f"Filtered cache to groups {args.groups}: {before:,} -> {len(allele_df):,} rows.")
```

The cache holds rows for `1D`, `2D`, `20`, `40`, `AL`. When `--groups 1D AL` is passed,
rows for `2D`, `20`, and `40` are dropped before any frequency diff computation. The
downstream diff, aggregation, and output Parquet files are all computed on the filtered
frame — there is no cross-group contamination.

The `groups_arg` in `allele_analysis.smk` is built directly from the Snakemake wildcard:

```python
# rules/allele_analysis.smk
params:
    groups_arg = lambda wildcards: "--groups " + " ".join(wildcards.groups.split("_"))
    # e.g. wildcards.groups = "1D_AL"  →  "--groups 1D AL"
    # e.g. wildcards.groups = "2D_AL"  →  "--groups 2D AL"
```

---

## Which QC Table a Cache Job Reads

There is one QC folder per single timepoint, so the answer needs no lookup: the cache job
of `(MAG, timepoint)` reads that MAG's QC table of that same timepoint,
`QC/QC_{timepoint}/{mag}_QC.tsv` (`get_timepoint_qc_file` in `shared/common.smk`).

Example: comparisons `pre_post` and `pre_end` name three timepoints.

| Cache job | QC table it reads |
|---|---|
| `(MAG_001, pre)` | `QC/QC_pre/MAG_001_QC.tsv` |
| `(MAG_001, post)` | `QC/QC_post/MAG_001_QC.tsv` |
| `(MAG_001, end)` | `QC/QC_end/MAG_001_QC.tsv` |

The `pre` samples are read by QC once and cached once, although two comparisons use them.

In the `allele_freq_cache` rule the declared input is the QC "done" marker of the timepoint
(`QC/QC_{timepoint}/.qc_done`), and the table's path is passed to the command as a
parameter. The QC step writes a table for every MAG of the run, header-only when the MAG
has no profiled sample at that timepoint, so the table is there once the marker is.

---

## Witness Caches and Zero-Row Caches

Normally a MAG gets a cache at the timepoints of the comparisons it is eligible in. One
case needs more. Under `baseline_scope: cohort`, with a collection date for every sample
and baseline presence or the gene-level de novo analysis switched on (longitudinal data),
a sample at **any** timepoint can show that an allele was already present. The workflow
then builds a cache for every MAG that is eligible in some comparison at **every witness
timepoint**: every compared timepoint, plus, under `cohort_witnesses: sheet` (the
default), every other timepoint of the sample sheet. See
[Which samples are the baseline](strain_turnover_analysis.md#which-samples-are-the-baseline-baseline_scope).

Example: comparisons `pre_post` and `pre_end`; the sheet also holds timepoint `mid`; MAG
`M1` is eligible in both comparisons, `M2` in `pre_post` only.

| | `M1` caches | `M2` caches |
|---|---|---|
| cohort condition not met | `pre`, `post`, `end` | `pre`, `post` |
| met, `cohort_witnesses: comparisons` | `pre`, `post`, `end` | `pre`, `post`, `end` |
| met, `cohort_witnesses: sheet` | `pre`, `post`, `end`, `mid` | `pre`, `post`, `end`, `mid` |

The list of cache files is fixed from the eligibility tables, before any QC table is read.
So a cache can be requested for a MAG that turns out to have no QC-passing sample at that
timepoint (`M2` at `mid`, say). The workflow runs the cache command with `--allow_empty`:
such a job writes a Parquet file with the normal columns and **zero rows**, plus a
header-only sidecar, logs a warning and exits 0. Baseline presence and the gene-level de
novo analysis read a zero-row cache as "no sample of this MAG at this timepoint"; a cache
file that is missing altogether stops them. Without `--allow_empty` (a call by hand) the
same case is an error, and a sample that passed QC but whose profile cannot be read is an
error either way.

Pairwise ANI and dN/dS never read the witness-only caches: they keep to the timepoints of
the comparisons a MAG is eligible in.

---

## End-to-End: One Combination

For `allele_analysis` of `(MRGM_1079, 5mo_22mo, 1D_AL)`:

```
allele_analysis (MRGM_1079, 5mo_22mo, 1D_AL)
    reads  → allele_freq_cache/5mo/MRGM_1079_5mo_allele_frequency.parquet   (all groups)
    reads  → allele_freq_cache/22mo/MRGM_1079_22mo_allele_frequency.parquet (all groups)
    filter → --groups 1D AL  (drops rows for 2D, 20, 40)
    writes → allele_analysis/allele_analysis_5mo_22mo-1D_AL/
               MRGM_1079_allele_frequency_changes_mean.parquet
               MRGM_1079_allele_frequency_changes_no_zero-diff.parquet

allele_freq_cache (MRGM_1079, 5mo)        ← shared by all 6 gr_combos + 3 tp_combos
    qc     → QC/QC_{timepoint}/{mag}_QC.tsv   (the QC table of this same timepoint)
    filter → time == "5mo"
    writes → allele_freq_cache/5mo/MRGM_1079_5mo_allele_frequency.parquet

allele_freq_cache (MRGM_1079, 22mo)
    qc     → QC/QC_{timepoint}/{mag}_QC.tsv   (the QC table of this same timepoint)
    filter → time == "22mo"
    writes → allele_freq_cache/22mo/MRGM_1079_22mo_allele_frequency.parquet
```

**Cache reuse in practice:** when `allele_analysis_5mo_10mo-1D_AL` runs first, it builds
the `5mo` and `10mo` cache files. When `allele_analysis_5mo_22mo-1D_AL` runs next,
Snakemake sees `5mo` already exists and skips that cache job — only `22mo` is built.
When `allele_analysis_5mo_22mo-2D_AL` runs after, both `5mo` and `22mo` already exist —
no additional cache writes, just reads and a different `--groups` filter.

---

## Output File Layout

```
{OUTDIR}/
├── QC/                                ← one folder per single timepoint, all groups
│   └── QC_{timepoint}/
│       └── {mag}_QC.tsv
│
├── allele_freq_cache/                 ← group-independent; shared across all gr_combos
│   ├── 5mo/
│   │   └── MRGM_1079_5mo_allele_frequency.parquet
│   ├── 10mo/
│   │   └── MRGM_1079_10mo_allele_frequency.parquet
│   ├── 16mo/
│   │   └── MRGM_1079_16mo_allele_frequency.parquet
│   └── 22mo/
│       └── MRGM_1079_22mo_allele_frequency.parquet
│
└── allele_analysis/                   ← one directory per (tp_combo, gr_combo)
    ├── allele_analysis_5mo_10mo-1D_AL/
    │   └── MRGM_1079_allele_frequency_changes_mean.parquet
    ├── allele_analysis_5mo_10mo-2D_AL/
    │   └── MRGM_1079_allele_frequency_changes_mean.parquet
    └── allele_analysis_5mo_10mo-20_AL/
        └── ...
```

---

## Downstream Consumers of the Cache

Statistical test rules that need per-sample data across both timepoints
(CMH across-time, LMM across-time) also read the Stage-1 Parquet files directly —
they pass two cache paths to `--input_df` and concatenate them internally via
`load_allele_freq_inputs()` in `utilities.py`, which handles `.parquet` vs `.tsv.gz`
by extension.

The `allele_analysis` rule (Stage 2) is different: it intentionally keeps the two
timepoints *separate* and merges on subject ID to compute per-subject diffs before
aggregating. It does not use `load_allele_freq_inputs` for the concatenation step.

---

## Source File Reference

| File | Role |
|---|---|
| [`shared/common.smk`](../../../alleleflux/smk_workflow/alleleflux_pipeline/shared/common.smk) | `get_timepoint_qc_file`, `get_allele_freq_cache_path`, `CACHE_TIMEPOINTS` (every timepoint that gets QC and caches) |
| [`shared/dynamic_targets.smk`](../../../alleleflux/smk_workflow/alleleflux_pipeline/shared/dynamic_targets.smk) | `_cache_universe`: which (MAG, timepoint) cache files a run builds and each consumer reads |
| [`rules/allele_analysis.smk`](../../../alleleflux/smk_workflow/alleleflux_pipeline/rules/allele_analysis.smk) | `allele_freq_cache` rule (Stage 1), `allele_analysis` rule (Stage 2) |
| [`scripts/analysis/allele_frequency/allele_freq_cache.py`](../../../alleleflux/scripts/analysis/allele_frequency/allele_freq_cache.py) | Stage 1: reads profiles for all groups at one timepoint, writes Parquet |
| [`scripts/analysis/allele_frequency/allele_freq.py`](../../../alleleflux/scripts/analysis/allele_frequency/allele_freq.py) | Stage 2: loads cache, applies `--groups` filter, computes diffs |
| [`scripts/utilities/utilities.py`](../../../alleleflux/scripts/utilities/utilities.py) | `load_allele_freq_inputs` — Parquet/TSV format detection for downstream consumers |
