# End-to-end regression test

Runs the bundled synthetic datasets through the whole pipeline and checks two
things: the outputs did not **change** (golden snapshot) and the outputs are
**right** (planted facts).  Both scenarios, longitudinal and single-timepoint,
run in profiles mode with every applicable feature switched on.  About four
minutes for both scenarios on four cores.

A second, independent workflow test lives in `visualization/`
(`test_visualization_workflow_e2e.py`); it runs the visualization Snakemake
workflow on a synthetic cache and is not described here.

## Commands

| Want to ... | Run |
|---|---|
| Fast unit loop, no pipeline | `pytest tests/ --ignore=tests/regression` |
| The regression, longitudinal only (default) | `pytest tests/regression/` |
| The regression, both scenarios | `pytest tests/regression/ --regression-scenario all` |
| One scenario | `pytest tests/regression/ --regression-scenario single` |
| Refresh the golden after an intentional output change | `pytest tests/regression/ --update-golden --regression-scenario all` |
| Compare two runs you made yourself | `python -m tests.regression.diff_runs <reference_tree> <candidate_tree>` |

Run from the repo root inside the `alleleflux` conda environment (the test
spawns the `alleleflux` command, so `conda run -n alleleflux python -m pytest
...` works too; the env's bare python binary does not).  Each run wipes and
rewrites `output/` and `.snakemake/` inside the dataset folders; neither is
committed.

## The two layers

**Golden snapshot** (`test_pipeline_e2e.py`).  The pipeline runs on the
scenario's dataset; every tracked output file is compared, row for row at the
comparator's strict tolerance (rtol 1e-9), against the committed copy under
`golden/<scenario>/`.  A golden file with no twin is `MISSING`; a file in the
scenario's `expected_absent` list that exists is `FORBIDDEN` (the single
scenario must not produce the longitudinal-only trees).  Extra files in the
run are tolerated.  Catches silent numerical drift and dropped outputs.

**Planted facts** (`test_planted_truth.py`).  The data generator decided every
significant site, the standing allele, the lone de novo allele, the block of
three, the strain-swapped sample and the depth spikes up front and wrote them
to `planted_truth.json` beside each dataset.  Six tests read that card and
look each fact up in the run's outputs: significant and rising in treatment;
standing with the planted mouse named at baseline; de novo with every carrier
absent before; one `block_event`; the swapped sample the only strain change;
the spikes removed for exactly the spiked samples.  A golden that was wrong on
the day it was captured would pass layer one forever; layer two is what catches
that.

Both layers read one pipeline run per scenario (`conftest.pipeline_run`), so
the second layer costs seconds.

## What is compared, what is not

The contract lives in `_scope.py` and is applied by the e2e test and
`diff_runs.py` alike.

Skipped (`FULL_EXCLUDE_SUBSTRINGS`):

| Pattern | Why |
|---|---|
| `_allele_frequency.parquet` | the allele-frequency cache; pinned through the site-filter sidecars beside it and every table read from it, not worth the golden's size |
| `profiles/` | copies of the input profiles |
| `inputMetadata/` | per-MAG sample sheets carrying absolute paths, so they differ between machines |
| `figures/` | PNG images and an HTML report, not comparable byte-wise |

Skipped on value only (`IGNORE_VALUE_SUBSTRINGS`): the LMM convergence
diagnostics (`*_warnings`, `*_converged_LMM`), whose values vary run to run
while the fitted statistics beside them do not.  The columns must still exist.

Everything else is compared: tables (`.tsv`, `.tsv.gz`, `.parquet`) row by
row after sorting both sides by their identity columns as plain strings, JSON
sidecars as parsed objects, anything else (the dN/dS `.ffn` sequences) byte
for byte.  The sign of a change column is compared like any other digit; a
former sign tolerance turned out to be a row-alignment bug in the comparator
(categorical keys sorted by code), fixed in `compare._sort_by_key_values`.

## Scenarios

`conftest.SCENARIOS` is the single source of truth: dataset folder, config,
output sub-folder, golden sub-folder, exclusions and forbidden trees.

| Scenario | Dataset | Config | Output under | Must be absent |
|---|---|---|---|---|
| `longitudinal` | `docs/source/examples/example_data_longitudinal/` | `config_generated_longitudinal.yml` | `output/longitudinal/` | nothing |
| `single` | `docs/source/examples/example_data_single/` | `config_generated_single.yml` | `output/single_timepoint/` | `dnds_analysis`, `baseline_presence`, `strain_turnover`, `regional_contrast` |

The configs switch on every feature the data type supports except LMM (about
58 CPU-minutes for two 5 kb genomes on this design, and a third of all
positions flagged at 4 mice per group; its own unit tests cover it).  Figures
are on for longitudinal and off for single.

The `config_with_bams_*.yml` configs and the `bams/` folders are the
tutorial's BAM-mode path.  Those BAMs are simulated from the unmodified
reference, so they carry none of the planted facts; the regression does not
use them.

## Comparing a branch with `main`

The golden under `golden/` is committed, so every branch carries main's
answers from the point it branched.  Running the regression on a branch
therefore **is** the main-versus-branch comparison:

1. On your branch, run `pytest tests/regression/ --regression-scenario all`.
   The pipeline runs with your code; the outputs are compared to the golden,
   which is main's.
2. Green: your branch produces exactly what main produced.  Red: the message
   lists every file and row that differs.
3. If the difference is intended, run the same command with
   `--update-golden`, read `git diff --stat tests/regression/golden`, and
   commit the new golden together with the change.  After the merge the
   golden on main is "main plus your change", ready for the next branch.

CI runs the regression on every push to main, so the committed golden is
always main's current output.

## Comparing any two runs

For two commits that have no golden between them, or two configs: run the
pipeline on each into its own folder, then

```bash
python -m tests.regression.diff_runs /path/to/run_A/longitudinal /path/to/run_B/longitudinal
```

The first tree is the reference.  Same comparator, same contract;
`--no-exclusions` compares the skipped files too.

## Adding a planted fact

1. In `docs/source/examples/generate_synthetic_data.py`, plant it in
   `plan_the_plants` (decide the ids, record them on the card) and write it
   in `generate_profile_data` (turn the ids into read counts).  Keep the plant
   clear of every other plant; the generator's `taken` set and its unit tests
   (`tests/examples/test_generate_synthetic_data.py`) enforce that.
2. Regenerate both datasets with the bundled seed, from
   `docs/source/examples/`:

   ```bash
   python generate_synthetic_data.py --output_dir example_data_longitudinal --data_type longitudinal --seed 42 --num_samples 16
   python generate_synthetic_data.py --output_dir example_data_single --data_type single --seed 42 --num_samples 8
   ```

   The generator does not delete sample folders it no longer writes; if the
   sample sheet shrank, remove the orphans by hand or the metadata step will
   refuse the run.
3. Add one test to `test_planted_truth.py` that reads the new card key and
   asserts where the fact surfaces.
4. Refresh the golden (`--update-golden`), review its diff, commit code, data,
   card and golden together.

## Refreshing the golden, and when not to

Refresh when an output change is intended: a new feature, a corrected
formula, a new column.  Review `git diff --stat tests/regression/golden`
before committing; the changed files should be the ones your change touches.

Do not refresh to make a red run green when you do not know why it is red.
The planted-fact tests are the tie-breaker: if they fail, the outputs are
wrong, not merely different.

The golden is about 22 MB on disk and about 6 MB in git; half of it is the
per-position `p_value_summary` tables.

## Files

| File | Role |
|---|---|
| `conftest.py` | `SCENARIOS`, the pytest options, and `pipeline_run`: one pipeline run per scenario, shared by every test here |
| `_scope.py` | the comparison contract: exclusions, value-ignored columns, forbidden trees |
| `compare.py` | the comparator: tables, JSON, bytes; `compare_trees` walks a golden tree |
| `test_pipeline_e2e.py` | golden snapshot test and `--update-golden` |
| `test_planted_truth.py` | the planted-fact assertions |
| `test_permutation_equivalence.py` | null (permuted) runs: reuse plus relabel equals a full rebuild; scoped to the statistics folders (`_scope.PERMUTATION_COMPARE_PREFIXES`) because the reuse path keeps label-stale QC summaries and the rebuild keeps the extras by design; three pipeline runs per scenario, LMM off, about 7 minutes for both |
| `test_compare.py` | unit tests of the comparator |
| `diff_runs.py` | compare any two output trees |
| `golden/` | the committed snapshots, one folder per scenario |

## Reading a failure

```
[scenario=longitudinal] 1 file(s) differ from the golden snapshot:

MAG_001_two_sample_paired.tsv.gz: 3 mismatched rows
  contig__actual  position__actual  p_value__actual  ...  p_value__golden
  ...
```

One block per file, the first ten mismatched rows side by side with `__actual`
and `__golden` suffixes, then `... and N more`.  `MISSING: <path>` means the
golden has a file the run did not produce; `FORBIDDEN (exists): <path>` means
the run produced a tree the scenario forbids.
