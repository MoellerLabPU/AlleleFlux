# Generated AlleleFlux Test Data

Seed: 42 | MAGs: 2 | Samples: 16 (longitudinal)

All paths in the configs and metadata are RELATIVE to this directory.
`cd` into this directory before running the pipeline.

`planted_truth.json` lists every fact planted in the profiles
(significant sites, standing / de novo alleles, the strain-swapped
sample, the depth spikes); the regression test asserts them.

## Option A — pre-generated profiles (fastest; carries the planted facts)

```bash
cd /home/su2806/AlleleFlux-dev/docs/source/examples/example_data_longitudinal
alleleflux run --config config_generated_longitudinal.yml
```

## Option B — full pipeline from BAMs (requires wgsim + bwa)

The BAMs are simulated from the unmodified reference, so this path
exercises profiling but carries NONE of the planted facts.

```bash
# 1. Generate BAMs from the reference (one-time)
bash generate_bams.sh --data-dir /home/su2806/AlleleFlux-dev/docs/source/examples/example_data_longitudinal

# 2. Run the pipeline (includes profiling step)
cd /home/su2806/AlleleFlux-dev/docs/source/examples/example_data_longitudinal
alleleflux run --config config_with_bams_longitudinal.yml
```
