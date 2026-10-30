# Score Figures and Run Report

After the scores are written, the workflow draws one set of figures per
comparison and test family (`alleleflux-plot-scores`) and then folds every
figure, the sample-and-MAG funnel, the significant-sites counts and the
site-depth-filter summary into one self-contained HTML page
(`alleleflux-run-report`).  Both run in the workflow by default and can be run
by hand on any finished run.

---

## What is produced and where

```text
{root_dir}/{data_type}/figures/
├── {comparison}/                               # one directory per comparison, e.g. pre_post-treatment_control
│   ├── score_heatmap_{family}.png              # rows = MAGs, tiles = cells, colour = score, stars = q
│   ├── score_boxplot_{family}.png              # one box per cell, one dot per MAG
│   ├── score_boxplot_{family}_by_phylum.png    # the same, one panel per phylum
│   └── score_figures.done                      # the files written (the workflow's sentinel)
└── run_report.html                             # every figure embedded; no external files
```

Figures are drawn only for comparisons that have combined MAG score tables
(`scores/processed/combined/MAG/`), and only for the families whose tables
exist; a family with no table is skipped with one log line.  Nothing is drawn
when `analysis.allele_analysis_only` is on, because there are no scores.  Both
steps are local rules: they read small tables and take seconds, so they never
go to the scheduler.

---

## The three figures

A **family** is one statistic drawn consistently across the **cells** of a
comparison: parallelism in the first group, parallelism in the second group,
and divergence between them.  With groups `treatment` and `control` a family
has three cells: `treatment parallelism`, `control parallelism`, `divergence`.
The six families and the tests behind each cell:

| Family | Parallelism cells (per group) | Divergence cell |
|--------|-------------------------------|-----------------|
| `ttest` | `single_sample_tTest` | `two_sample_paired_tTest` |
| `wilcoxon` | `single_sample_Wilcoxon` | `two_sample_paired_Wilcoxon` |
| `ttest_unpaired` | `single_sample_tTest` (shared with `ttest`) | `two_sample_unpaired_tTest` |
| `mannwhitney` | `single_sample_Wilcoxon` (shared with `wilcoxon`) | `two_sample_unpaired_MannWhitney` |
| `lmm` | `LMM_across_time` | `LMM` |
| `cmh` | `CMH_across_time` | `CMH`, one cell per focus timepoint, labelled `divergence (post)` |

`ttest` and `wilcoxon` use the paired divergence test; `ttest_unpaired` and
`mannwhitney` use the unpaired one, which is the one to read when samples
cannot be paired across timepoints (different subjects at `pre` and `post`).
On `single` data there is no within-group test, so every family has the
divergence cell only.

A MAG is **tested** in a cell when it has a row in that cell's score table and
the number of sites behind the score is not zero.  The rows of every figure are
the union of the MAGs tested in any cell of the family; a cell the MAG was not
tested in is drawn **grey** (never 0, never blank).  When the pipeline wrote no
table at all for a cell (no MAG was eligible for that test), the column is
all grey and its label carries the note `(no eligible MAGs)` under the cell name.  The star on a tile
is the MAG's smallest BH-adjusted q over all its sites in that cell, taken from
the significant-sites summary
(`p_value_summary/significant_sites_summary/significant_sites_mag_cell_stats_long.tsv`):
`***` q < 0.001, `**` q < 0.01, `*` q below the run's threshold
(`statistics.p_value_threshold`), nothing otherwise.  The score table and the
summary are built from the same p-value file, so the tool stops with an error
if the tested MAGs of a cell differ between the two.

**Heatmap** (`score_heatmap_{family}`).  Three blocks from left to right: the
row labels (MAG id and species, species cut to 40 characters), a one-column
phylum strip with its own legend (unknown phylum grey), and the score tiles,
one column per cell, coloured 0-100 % on the `plasma` scale with the q stars in
white.  A colour bar labelled "score (% sites significant)" sits to the right.
Rows are in taxonomy order by default: phylum, class, order, family, genus,
species, then MAG id, with blanks last.  With a Newick tree (`tree_path` /
`--tree`, MAG ids as tip names) rows follow the tree's ladderized tip order
instead, and MAGs the tree does not contain are drawn below a thin divider.
Saved at 200 dpi so the 6-point labels stay readable.

**Boxplot** (`score_boxplot_{family}`).  One box per cell and one dot per MAG;
a MAG tested in both parallelism cells gets a grey line joining its two dots, so
a MAG that is parallel in one group and not the other is easy to spot.  The
number of MAGs in each box is written under it as `n=` (`n=0` for a cell with
no tested MAG).  Dots are jittered once per MAG, so a line lands exactly on both
of its dots.

**Boxplot by phylum** (`score_boxplot_{family}_by_phylum`).  The same boxplot,
one panel per phylum that has at least `min_mags_per_phylum` MAGs among the
rows, titled `<phylum> (n=K)`; every smaller phylum is pooled into one panel
titled `other (<N> phyla, n=M)`.

---

## The run report

`run_report.html` is one file with inline CSS, every image embedded and no
scripts, so it can be mailed or dropped into a shared folder as-is.  A section
whose input is missing says `not available in this run` in one line; the rest
of the page still renders.  A permuted (null) run gets a banner at the top
saying the group labels were shuffled before the tests.  Sections, in order:

**Run.**  Run name, run directory, data type, report date, AlleleFlux version,
code commit and the q threshold, then the `quality_control` block of the saved
runtime configuration, with the full configuration folded away below it.  The
configuration is the newest `config_runtime_*.yml` in the run root (the parent
of the data-type directory) unless `--config` names one.

**Funnel: samples and MAGs through the pipeline.**  One table per comparison,
top to bottom:

- `samples in the metadata sheet (treatment, pre)` — one row per group and
  timepoint in the sample sheet.
- `MAGs with profiles (pre_post)` — distinct MAGs in
  `inputMetadata_{timepoints}/summary_mag_index.tsv`.
- `MAGs with enough samples ...` — one row per yes/no flag of the comparison's
  eligibility table, worded by which tests the flag feeds: `for the
  within-group tests in treatment (parallelism)` (one row per group),
  `for the unpaired between-group tests and LMM (divergence)` and `for the
  paired between-group tests and CMH (divergence)`.
- `MAGs with enough sites after preprocessing ...` — the same wording, counted
  over the flags of the between-groups and within-groups preprocessing tables.

Flag columns are recognised by holding only `True`/`False`; the
`positions_eligible_*` site counts in the preprocessing tables are not flags
and are not reported.

**Significant sites per test cell.**  One table per comparison with a row per
test cell (`test type`, `group`): `MAGs tested` is the number of MAGs with at
least one site in that cell's p-value table, i.e. MAGs that passed eligibility
and preprocessing for the test; `MAGs with a significant site` is, of those,
the MAGs whose smallest BH-adjusted q over all their sites is below the q
threshold (one significant site is enough to count; this is exactly
`min_q_value < threshold` in the significant-sites summary); `% of tested` is
the previous column as a percentage, drawn again as a bar.  `group` is the group
tested for a within-group test, the focus timepoint for `CMH` divergence and
blank for the other between-group tests.  A definitions list under the tables
repeats these meanings and explains the `_abs` suffix (the test run on the
absolute size of each sample's frequency change).

**Site-depth filter.**  One row per timepoint summarising the per-(sample, MAG)
filter table `site_depth_filter/site_depth_filter_summary.tsv`:
`sample_mag_pairs_checked` (how many sample × MAG combinations the filter ran
on; each pair gets its own mean depth and therefore its own ceiling), `mags`
(distinct MAGs among those pairs), `sites_covered` (positions with any reads,
summed over all pairs, before the filter), `sites_removed` (positions dropped by
the ceiling plus the floor), `pct_removed` (removed as a percentage of covered:
the overall share of the data the filter took away), `median_row_pct` (for each
pair, the percentage of its own covered positions removed; the median over
pairs, i.e. what a typical pair lost) and `max_row_pct` (the same for the single
pair that lost the most; high with a low median points at one MAG with a large
repeat or cross-mapping region rather than a run-wide problem).  Below the
table, a histogram of the fraction removed per MAG and timepoint from the
`..._by_mag.tsv` rollup.

**Figures.**  An index table (comparison, family, links to its three figures)
followed by every `score_*.png` found under `figures/<comparison>/`, each with a
caption that restates the legend.

:::{note}
A PNG of every figure is always written, because the report embeds PNGs.
`format: pdf` adds a PDF copy of each figure next to it (`both` means the
same and is kept for older configurations).
:::

---

## Configuration

All five keys sit under `analysis.figures`; leaving the block out gives the
defaults.

```yaml
analysis:
  figures:
    enabled: True              # False = no figures, no report
    families: []               # [] = every family with a score table; else a subset of
                               # ttest, ttest_unpaired, wilcoxon, mannwhitney, lmm, cmh
    tree_path: ""              # optional Newick tree with MAG ids as tips: heatmap rows follow it
    format: "png"              # png | pdf (PNG always written; pdf adds a PDF copy)
    min_mags_per_phylum: 3     # by-phylum boxplot: smaller phyla are pooled into one "other" panel
```

The star threshold is not a figures key: it is the run's
`statistics.p_value_threshold`, so a star means the same thing as "significant"
everywhere else in the run.  A tree given as `tree_path` is an input of the
figure rule, so replacing the file redraws the figures.

---

## Running by hand on a finished run

Both commands read a finished run's data-type directory (the one holding
`scores/` and `p_value_summary/`); nothing is recomputed.

```bash
# Figures for one comparison; the two group labels in configuration order
alleleflux-plot-scores --run_dir run/longitudinal \
  --comparison pre_post-treatment_control --groups treatment control \
  --output_dir run/longitudinal/figures/pre_post-treatment_control \
  --families ttest lmm --tree mags.nwk --format both

# The report over every comparison's figures
alleleflux-run-report --run_dir run/longitudinal \
  --figures_dir run/longitudinal/figures \
  --output run/longitudinal/figures/run_report.html \
  --comparisons pre_post-treatment_control
```

The report finds the sample sheet and the `inputMetadata/` directory through
the newest `config_runtime_*.yml` in the run root (under `input.reuse_from`
when profiles were reused); pass `--metadata_path`, `--metadata_dir` or
`--config` when the run has been moved.  Add `--permuted` to a null run made
without the workflow's permutation settings.  Flag-by-flag details are in the
[CLI Reference](../reference/cli_reference.md).

---

## Reading the figures

The tile colour and the star answer two different questions.  The **score** is
the percentage of a MAG's tested sites whose raw p-value is below the
threshold; the **star** comes from the smallest BH-adjusted q over the same
sites.  A MAG with many sites can have a handful below p = 0.05 by chance
alone — a non-zero score — while no site survives the multiple-testing
correction, so it gets no star.  Read the score as "how much of the genome
moved" and the star as "is at least one site reliably significant"; the
report's `MAGs with a significant site` column counts the starred MAGs.
