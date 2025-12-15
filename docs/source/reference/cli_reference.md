# CLI Reference

AlleleFlux ships a main `alleleflux` entrypoint plus console scripts that power the Snakemake workflow. Most users only need `alleleflux run`; the other tools are available for advanced or ad‑hoc use.

## Main commands

### `alleleflux run` — execute the workflow

```bash
alleleflux run --config config.yml [options] [-- <extra snakemake args>]
```

```{eval-rst}
.. list-table::
   :widths: 28 18 54
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``-c, --config``
     - (required)
     - Path to the AlleleFlux configuration YAML.
   * - ``-w, --working-dir``
     - ``.``
     - Working directory for Snakemake execution.
   * - ``-j, --jobs``
     - None
     - Max concurrent jobs (local only; ignored when ``--profile`` is set).
   * - ``-t, --threads``
     - None
     - Total threads available for local runs.
   * - ``-m, --memory``
     - None
     - Total memory for local runs (e.g., ``64G``).
   * - ``-p, --profile``
     - None
     - Snakemake profile directory for cluster/HPC execution.
   * - ``-n, --dry-run``
     - False
     - Plan the DAG without running jobs.
   * - ``--unlock``
     - False
     - Unlock a previously crashed working directory.
   * - ``--snakemake-args``
     - None
     - Quoted string of extra Snakemake flags (alternative to ``--``).
```

Pass additional Snakemake flags either after `--` or via `--snakemake-args` (e.g., `alleleflux run -c config.yml -- --forceall --reason`). See {doc}`../usage/running_workflow` for scheduling details.

### `alleleflux init` — create a config

```bash
alleleflux init [--template] [--output alleleflux_config.yml]
```

- `--template` prints the bundled template to stdout.
- Without `--template`, an interactive prompt writes the config to `--output` (default: `alleleflux_config.yml`).

### `alleleflux info` — show install paths

Print version, package location, and the packaged Snakefile. No options.

### `alleleflux tools` — list console scripts

```bash
alleleflux tools [--category {Analysis,Preprocessing,Statistics,Evolution,Accessory,Visualization}]
```

Lists every console script shipped with AlleleFlux, grouped by stage.

## Console scripts by stage

These are invoked automatically by the workflow but can be run manually for testing or custom tasks. Run any script with `--help` for full arguments.

### Analysis

- `alleleflux-profile` — profile BAMs into per-sample MAG tables.
- `alleleflux-allele-freq` — compute allele frequency tables per MAG.
- `alleleflux-scores` / `alleleflux-taxa-scores` / `alleleflux-gene-scores` — derive MAG, taxa, and gene-level scores.
- `alleleflux-outliers` — flag outlier genes.
- `alleleflux-cmh-scores` — CMH-specific score aggregation.

### Preprocessing

- `alleleflux-metadata` — build MAG metadata from profiles + sample sheet.
- `alleleflux-qc` — coverage/breadth QC on profiles.
- `alleleflux-eligibility` — QC-based MAG eligibility tables.
- `alleleflux-preprocess-between-groups` / `alleleflux-preprocess-within-group` — position-level filtering before tests.
- `alleleflux-preprocessing-eligibility` — aggregate preprocessing status into eligibility tables.
- `alleleflux-p-value-summary` — summarize p-values for downstream steps (e.g., dN/dS).

### Statistics

- `alleleflux-two-sample-unpaired` / `alleleflux-two-sample-paired` — two-sample tests.
- `alleleflux-single-sample` — within-group test.
- `alleleflux-lmm` — linear mixed models.
- `alleleflux-cmh` — Cochran–Mantel–Haenszel test.

### Evolution

- `alleleflux-dnds-from-timepoints` — dN/dS from significant sites (see {doc}`../usage/dnds_analysis`).

### Accessory

- `alleleflux-create-mag-mapping` — contig→MAG mapping and combined FASTA.
- `alleleflux-add-bam-path` — fill `bam_path` values in metadata.
- `alleleflux-coverage-allele-stats` — coverage/allele stats summary.
- `alleleflux-list-mags` — enumerate MAG IDs in a profiles directory.
- `alleleflux-positions-qc` — position-level QC filtering.
- `alleleflux-copy-profiles` — copy or symlink profile files.

### Visualization

- `alleleflux-prepare-metadata` — prep metadata for visualization inputs.
- `alleleflux-terminal-nucleotide` — terminal nucleotide analysis.
- `alleleflux-track-alleles` — track allele trajectories.
- `alleleflux-plot-trajectories` — plot tracked allele trajectories.

## Getting help

```bash
alleleflux-<tool> --help
```

For configuration details, see {doc}`configuration`. For how to run the workflow end to end, see {doc}`../usage/running_workflow`.
