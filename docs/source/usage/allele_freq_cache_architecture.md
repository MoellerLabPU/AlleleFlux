# Allele Frequency Cache Architecture

## What This Does (Plain Language)

For each group combination (e.g. `1D_AL`), AlleleFlux builds **one Parquet cache file per
timepoint** (e.g. `MAG1_1D_AL_5mo_allele_frequency.parquet`). That file is written once —
profile TSVs for `5mo` are read exactly once for `1D_AL` — and then every **timepoint
combination** that includes `5mo` for that **same group combination** reads from the same
cache instead of re-reading the raw profiles.

Cache files are **not shared across group combinations**. `MAG1_2D_AL_5mo_allele_frequency.parquet`
is a separate file built independently from `MAG1_1D_AL_5mo_allele_frequency.parquet`,
because the two group comparisons involve different sets of samples.

The rest of this page explains how the pipeline implements this and why a small bookkeeping
structure (`timepoint_gr_to_canonical_tp` in `shared/common.smk`) is needed to make it work
with Snakemake's checkpoint system.

---

This page also explains the two-stage allele frequency computation in detail.
Understanding it is useful if you are debugging the pipeline, tuning resource allocation,
or modifying the workflow rules.

---

## The Problem: Redundant Profile Reads

When a config defines multiple timepoint combinations that share a timepoint — e.g. `[5mo, 10mo]`
and `[5mo, 16mo]` both contain `5mo` — the naive approach is to read every sample's profile TSV
separately for each combination. With many combinations this becomes expensive.

**Example using the drido config** (15 tp_combos × 6 gr_combos):

```
timepoints_combinations:
  - [5mo, 10mo]    ← 5mo appears here
  - [5mo, 16mo]    ← 5mo appears here again
  - [5mo, 22mo]    ← and again
  - [5mo, 28mo]    ← ...
  - [5mo, 34mo]
  - [5mo, 40mo]
  - [10mo, 16mo]
  - [16mo, 22mo]
  - [22mo, 28mo]
  - [28mo, 34mo]
  - [8mo, 10mo]
  - [8mo, 16mo]
  - [8mo, 22mo]
  - [8mo, 28mo]
  - [8mo, 34mo]
```

Under the naive approach, `5mo` sample profiles would be read and frequency-computed **36 times**
(6 tp_combos containing `5mo` × 6 gr_combos). The resulting per-combination longitudinal TSV files
were also duplicated on disk (~600 MB per MAG per combination).

---

## The Solution: Two-Stage Cache

AlleleFlux splits allele frequency processing into two stages:

```
Stage 1 — compute_allele_freq_per_timepoint         (one job per MAG × gr_combo × single timepoint)
           alleleflux-cache-allele-freq
           Output: allele_freq_cache/{mag}_{groups}_{timepoint}_allele_frequency.parquet

Stage 2 — allele_analysis                           (one job per MAG × tp_combo × gr_combo)
           alleleflux-allele-freq
           Input:  two Stage-1 Parquet files (one per timepoint in the combo)
           Output: _allele_frequency_changes_mean.tsv.gz  (and optionally _no_zero-diff.tsv.gz)
```

For drido:
- **Before**: 15 tp_combos × 6 gr_combos × 2 timepoints = **180 profile-read operations** per MAG
- **After**: 8 unique timepoints × 6 gr_combos = **48 cache-write jobs** per MAG (**4× fewer reads**)

Statistical tests that previously read a `_allele_frequency_longitudinal.tsv.gz` — namely
`lmm_analysis_across_time`, `cmh_test_across_time`, and `cmh_test` — now read the Stage-1
Parquet files directly via `--input_df` (both scripts accept `nargs="+"` and handle `.parquet`
extension automatically through `load_allele_freq_inputs` in
[`alleleflux/scripts/utilities/utilities.py`](../../alleleflux/scripts/utilities/utilities.py)).

---

## How the Cache Rule Knows Which QC File to Use

The cache rule `compute_allele_freq_per_timepoint` in
[`rules/allele_analysis.smk`](../../alleleflux/smk_workflow/alleleflux_pipeline/rules/allele_analysis.smk)
needs a QC file to know which samples passed coverage/breadth thresholds at a given timepoint.

The complication is that AlleleFlux uses a **Snakemake checkpoint** (`eligibility_table`) that
resolves one `(tp_combo, gr_combo)` at a time. When the very first combination's jobs run,
only that combination's QC directory exists — `QC_5mo_10mo-1D_AL/` is present but
`QC_5mo_16mo-1D_AL/` may not be yet. The cache rule cannot safely list QC files from
combinations that haven't been resolved.

The fix: use **one canonical QC file** per `(gr_combo, timepoint)`, chosen as the first
tp_combo in config order that contains that timepoint. This file can always be built through
a normal rule dependency chain (not gated by another combination's checkpoint).

This canonical file selection is implemented by two data structures built in
[`shared/common.smk`](../../alleleflux/smk_workflow/alleleflux_pipeline/shared/common.smk):

---

## `unique_timepoints`: Extracting Individual Timepoints

```python
# shared/common.smk
unique_timepoints = []
_seen_tps = set()
for tp_label in timepoints_labels:
    if DATA_TYPE == "longitudinal":
        parts = tp_label.split("_")  # "5mo_10mo" → ["5mo", "10mo"]
    else:
        parts = [tp_label]           # "5mo" → ["5mo"]
    for tp in parts:
        if tp not in _seen_tps:
            _seen_tps.add(tp)
            unique_timepoints.append(tp)
```

`timepoints_labels` is the list of Snakemake wildcard labels — for drido these are the
15 combined labels like `"5mo_10mo"`, `"8mo_34mo"`, etc.

The loop splits each label on `"_"` to get the two constituent timepoints, then adds each
timepoint to `unique_timepoints` exactly once (using `_seen_tps` to skip duplicates).

**Concrete trace for drido** (first few iterations):

| `tp_label` | `parts` | `tp` | Already seen? | `unique_timepoints` after |
|---|---|---|---|---|
| `"5mo_10mo"` | `["5mo", "10mo"]` | `"5mo"` | No | `["5mo"]` |
| | | `"10mo"` | No | `["5mo", "10mo"]` |
| `"5mo_16mo"` | `["5mo", "16mo"]` | `"5mo"` | **Yes** | unchanged |
| | | `"16mo"` | No | `["5mo", "10mo", "16mo"]` |
| `"5mo_22mo"` | `["5mo", "22mo"]` | `"5mo"` | **Yes** | unchanged |
| | | `"22mo"` | No | `["5mo", "10mo", "16mo", "22mo"]` |
| `"10mo_16mo"` | `["10mo", "16mo"]` | `"10mo"` | **Yes** | unchanged |
| | | `"16mo"` | **Yes** | unchanged |
| `"8mo_10mo"` | `["8mo", "10mo"]` | `"8mo"` | No | `[..., "8mo"]` |
| | | `"10mo"` | **Yes** | unchanged |

**Final result** (in order of first appearance in the config):

```python
unique_timepoints = ["5mo", "10mo", "16mo", "22mo", "28mo", "34mo", "40mo", "8mo"]
#                                                                               ↑
#                            8mo appears last — it only occurs in Strategy C
```

This list becomes the Snakemake `{timepoint}` wildcard constraint (the single-timepoint
wildcard used in cache filenames), distinct from `{timepoints}` (the combo label like
`"5mo_10mo"`).

---

## `timepoint_gr_to_canonical_tp`: Selecting the Canonical QC Source

```python
# shared/common.smk
timepoint_gr_to_canonical_tp = {}
for tp_label in timepoints_labels:
    if DATA_TYPE == "longitudinal":
        parts = tp_label.split("_")  # "5mo_10mo" → ["5mo", "10mo"]
    else:
        parts = [tp_label]
    for tp in parts:
        for gr_label in groups_labels:
            key = (tp, gr_label)
            if key not in timepoint_gr_to_canonical_tp:
                timepoint_gr_to_canonical_tp[key] = tp_label  # first in config order wins
```

This builds a dict mapping `(individual_timepoint, gr_combo) → first tp_label that contains
that timepoint`. The `if key not in ...` guard means later tp_labels that also contain a
timepoint are silently ignored — first one wins.

**Concrete trace for drido** (showing `1D_AL` gr_combo only — all 6 behave identically):

```
tp_label = "5mo_10mo":   # First iteration
  tp="5mo",  gr="1D_AL"  → ("5mo",  "1D_AL"): "5mo_10mo"  ← recorded
  tp="10mo", gr="1D_AL"  → ("10mo", "1D_AL"): "5mo_10mo"  ← recorded

tp_label = "5mo_16mo":
  tp="5mo",  gr="1D_AL"  → key ALREADY SET, skip
  tp="16mo", gr="1D_AL"  → ("16mo", "1D_AL"): "5mo_16mo"  ← recorded

tp_label = "5mo_22mo":
  tp="5mo"   → skip (already set)
  tp="22mo", gr="1D_AL"  → ("22mo", "1D_AL"): "5mo_22mo"  ← recorded

tp_label = "10mo_16mo":  # Strategy B — both timepoints already seen
  tp="10mo"  → skip
  tp="16mo"  → skip

tp_label = "8mo_10mo":   # Strategy C — 8mo appears for the first time
  tp="8mo",  gr="1D_AL"  → ("8mo",  "1D_AL"): "8mo_10mo"  ← recorded
  tp="10mo"  → skip
```

**Final dict** (all 6 gr_combos shown for a few timepoints):

```python
timepoint_gr_to_canonical_tp = {
    ("5mo",  "1D_AL"): "5mo_10mo",   ("5mo",  "2D_AL"): "5mo_10mo",  # ... ×6
    ("10mo", "1D_AL"): "5mo_10mo",   ("10mo", "2D_AL"): "5mo_10mo",  # 10mo first in 5mo_10mo
    ("16mo", "1D_AL"): "5mo_16mo",   ("16mo", "2D_AL"): "5mo_16mo",  # 16mo first in 5mo_16mo
    ("22mo", "1D_AL"): "5mo_22mo",   ...
    ("28mo", "1D_AL"): "5mo_28mo",   ...
    ("34mo", "1D_AL"): "5mo_34mo",   ...
    ("40mo", "1D_AL"): "5mo_40mo",   ...
    ("8mo",  "1D_AL"): "8mo_10mo",   ("8mo",  "2D_AL"): "8mo_10mo",  # 8mo only in Strategy C
}
```

---

## `get_canonical_qc_file`: Turning the Map into a Path

```python
# shared/common.smk
def get_canonical_qc_file(mag_wildcard, timepoint, groups):
    canonical_tp = timepoint_gr_to_canonical_tp.get((timepoint, groups))
    return os.path.join(
        OUTDIR, "QC",
        f"QC_{canonical_tp}-{groups}",
        f"{mag_wildcard}_QC.tsv",
    )
```

**Examples:**

```python
get_canonical_qc_file("MAG1", "5mo", "1D_AL")
# canonical_tp = "5mo_10mo"
# → "out/QC/QC_5mo_10mo-1D_AL/MAG1_QC.tsv"

get_canonical_qc_file("MAG1", "16mo", "2D_AL")
# canonical_tp = "5mo_16mo"
# → "out/QC/QC_5mo_16mo-2D_AL/MAG1_QC.tsv"

get_canonical_qc_file("MAG1", "8mo", "40_AL")
# canonical_tp = "8mo_10mo"   ← 8mo is only in Strategy C
# → "out/QC/QC_8mo_10mo-40_AL/MAG1_QC.tsv"
```

This path points to a QC file from a specific, always-buildable combination — no checkpoint
racing, no missing directories.

---

## `get_allele_freq_cache_path`: Naming the Output Parquet

```python
# shared/common.smk
def get_allele_freq_cache_path(
    mag_wildcard="{mag}",
    groups_wildcard="{groups}",
    timepoint_wildcard="{timepoint}",
):
    return os.path.join(
        OUTDIR, "allele_freq_cache",
        f"{mag_wildcard}_{groups_wildcard}_{timepoint_wildcard}_allele_frequency.parquet",
    )
```

The gr_combo is part of the cache filename. This means the `5mo` timepoint produces **six
separate cache files** — one per gr_combo — rather than one shared file. This is necessary
because each gr_combo has a different set of samples (determined by its QC file), so the
per-sample allele frequency data differs between gr_combos.

```
allele_freq_cache/
  MAG1_1D_AL_5mo_allele_frequency.parquet   ← samples from 1D vs AL comparison
  MAG1_2D_AL_5mo_allele_frequency.parquet   ← samples from 2D vs AL comparison
  MAG1_20_AL_5mo_allele_frequency.parquet
  ...
```

**Cache files are NOT shared across gr_combos.** `MAG1_2D_AL_5mo_allele_frequency.parquet`
is a completely separate file from `MAG1_1D_AL_5mo_allele_frequency.parquet` and is built
by its own independent cache job. The two files may contain a different set of samples
because QC pass/fail is evaluated per gr_combo.

**Deduplication applies within a gr_combo, across tp_combos that share a timepoint.**
For `1D_AL`, `5mo` appears in six tp_combos (`5mo_10mo`, `5mo_16mo`, …, `5mo_40mo`).
All six read from the same `MAG1_1D_AL_5mo_allele_frequency.parquet` — it is written once
by the first tp_combo that triggers the cache job, and Snakemake skips it for the rest.

---

## End-to-End DAG: One Combination

When Snakemake needs the output of `allele_analysis_5mo_22mo-1D_AL` for `MAG1`:

```
allele_analysis (MAG1, 5mo_22mo, 1D_AL)
    needs → MAG1_1D_AL_5mo_allele_frequency.parquet
    needs → MAG1_1D_AL_22mo_allele_frequency.parquet

compute_allele_freq_per_timepoint (MAG1, 1D_AL, 5mo)    [Stage 1]
    QC input:  get_canonical_qc_file("MAG1", "5mo", "1D_AL")
             = "QC/QC_5mo_10mo-1D_AL/MAG1_QC.tsv"
    Script:    alleleflux-cache-allele-freq
                 --magID MAG1
                 --qc_files QC/QC_5mo_10mo-1D_AL/MAG1_QC.tsv
                 --timepoint 5mo              # longitudinal only
                 --output_path ...5mo...parquet

compute_allele_freq_per_timepoint (MAG1, 1D_AL, 22mo)   [Stage 1]
    QC input:  get_canonical_qc_file("MAG1", "22mo", "1D_AL")
             = "QC/QC_5mo_22mo-1D_AL/MAG1_QC.tsv"

allele_analysis (MAG1, 5mo_22mo, 1D_AL)                 [Stage 2]
    Script:    alleleflux-allele-freq
                 --magID MAG1
                 --cache_files ...5mo...parquet ...22mo...parquet
                 --groups 1D AL
                 --output_dir allele_analysis/allele_analysis_5mo_22mo-1D_AL/
```

When `allele_analysis_5mo_10mo-1D_AL` runs first and builds
`MAG1_1D_AL_5mo_allele_frequency.parquet`, then `allele_analysis_5mo_22mo-1D_AL` runs —
same gr_combo (`1D_AL`), different tp_combo — Snakemake sees that
`MAG1_1D_AL_5mo_allele_frequency.parquet` already exists and skips the `5mo` cache job.
Only `MAG1_1D_AL_22mo_allele_frequency.parquet` needs to be built.

`allele_analysis_5mo_22mo-2D_AL` is a **different gr_combo** (`2D_AL`) and will build
its own separate `MAG1_2D_AL_5mo_allele_frequency.parquet` independently — cache files
are never shared across gr_combos.

---

## Why CMH.py and LMM.py Now Receive Multiple Input Files

### Before the refactor

`allele_freq.py` wrote a single `_allele_frequency_longitudinal.tsv.gz` file that bundled
**all samples from both timepoints** in the combination into one table. Scripts that needed
data across both timepoints — `cmh_test`, `lmm_analysis_across_time`, and
`cmh_test_across_time` — each read this one file:

```
allele_analysis_5mo_10mo-1D_AL/
  MAG1_allele_frequency_longitudinal.tsv.gz   ← rows for 5mo AND 10mo samples, in one file
```

The corresponding Snakemake rule passed a single path:

```
--input_df allele_analysis_5mo_10mo-1D_AL/MAG1_allele_frequency_longitudinal.tsv.gz
```

### After the refactor

The longitudinal TSV is eliminated. The cache stores **one Parquet file per timepoint**, not
one file per combination. A combination such as `5mo_10mo` is represented by two separate
Parquet files:

```
allele_freq_cache/
  MAG1_1D_AL_5mo_allele_frequency.parquet    ← 5mo samples only
  MAG1_1D_AL_10mo_allele_frequency.parquet   ← 10mo samples only
```

Scripts that need data from both timepoints — exactly the `cmh_test`,
`lmm_analysis_across_time`, and `cmh_test_across_time` rules — must now receive **two paths**:

```
--input_df allele_freq_cache/MAG1_1D_AL_5mo_allele_frequency.parquet \
           allele_freq_cache/MAG1_1D_AL_10mo_allele_frequency.parquet
```

Inside each script, `load_allele_freq_inputs` (in
[`scripts/utilities/utilities.py`](../../alleleflux/scripts/utilities/utilities.py))
handles the concatenation before any statistical logic runs:

```python
# utilities.py — called at the start of CMH.py and LMM.py main()
def load_allele_freq_inputs(paths):
    frames = []
    for p in paths:
        if str(p).endswith(".parquet"):
            frames.append(pd.read_parquet(p))
        else:
            frames.append(pd.read_csv(p, sep="\t"))
    return pd.concat(frames, ignore_index=True)
```

The `--input_df` argument in both CMH.py and LMM.py uses `nargs="+"`, so it already accepted
multiple paths before the refactor — no changes to the statistical scripts were needed.

### Why Stage 2 (`allele_analysis`) is different

`allele_analysis` also reads two cache files, but it keeps the two timepoints **separate** on
purpose: it computes per-subject diffs (value at `post` minus value at `pre`) and only then
aggregates. So Stage 2 reads the two Parquet files, holds them in memory as distinct
DataFrames, merges on subject, computes the diff, and writes `_changes_mean.tsv.gz`. It does
**not** concatenate them blindly the way the statistical tests do.

The statistical tests (`cmh_test`, `lmm_analysis_across_time`, `cmh_test_across_time`)
receive the concatenated data and use the `time` column to distinguish timepoints internally,
which is why a simple `pd.concat` works for them.

---

## Relevant Source Files

| File | Role |
|---|---|
| [`shared/common.smk`](../../alleleflux/smk_workflow/alleleflux_pipeline/shared/common.smk) | Builds `unique_timepoints`, `timepoint_gr_to_canonical_tp`; defines `get_canonical_qc_file`, `get_allele_freq_cache_path` |
| [`rules/allele_analysis.smk`](../../alleleflux/smk_workflow/alleleflux_pipeline/rules/allele_analysis.smk) | `compute_allele_freq_per_timepoint` (Stage 1 rule), `allele_analysis` (Stage 2 rule) |
| [`scripts/analysis/allele_freq_cache.py`](../../alleleflux/scripts/analysis/allele_freq_cache.py) | Stage 1 script — loads QC, reads profiles in parallel, writes Parquet |
| [`scripts/analysis/_allele_freq_common.py`](../../alleleflux/scripts/analysis/_allele_freq_common.py) | Shared helpers: `process_mag_files`, `calculate_frequencies`, `CACHE_COLUMNS_*` |
| [`scripts/analysis/allele_freq.py`](../../alleleflux/scripts/analysis/allele_freq.py) | Stage 2 script — reads Parquet cache files, computes diffs and means |
| [`scripts/utilities/utilities.py`](../../alleleflux/scripts/utilities/utilities.py) | `load_allele_freq_inputs` — handles `.parquet` or `.tsv.gz`, single path or list |
