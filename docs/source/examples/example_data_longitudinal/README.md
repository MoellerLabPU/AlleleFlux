# Generated AlleleFlux Test Data

Seed: 42 | MAGs: 2 | Samples: 8 (longitudinal)

## Option A — pre-generated profiles (fastest)

```bash
alleleflux run --config /home/su2806/AlleleFlux-dev/docs/source/examples/example_data_longitudinal/config_generated_longitudinal.yml
```

## Option B — full pipeline from BAMs (requires wgsim + bwa)

```bash
# 1. Generate BAMs from the reference (one-time)
bash generate_bams.sh --data-dir /home/su2806/AlleleFlux-dev/docs/source/examples/example_data_longitudinal

# 2. Run the pipeline (includes profiling step)
alleleflux run --config /home/su2806/AlleleFlux-dev/docs/source/examples/example_data_longitudinal/config_with_bams_longitudinal.yml
```
