# Generated AlleleFlux Test Data

Seed: 42 | MAGs: 2 | Samples: 8 (single)

## Option A — pre-generated profiles (fastest)

```bash
alleleflux run --config /home/su2806/AlleleFlux-dev/docs/source/examples/example_data_single/config_generated_single.yml
```

## Option B — full pipeline from BAMs (requires wgsim + bwa)

```bash
# 1. Generate BAMs from the reference (one-time)
bash generate_bams.sh --data-dir /home/su2806/AlleleFlux-dev/docs/source/examples/example_data_single

# 2. Run the pipeline (includes profiling step)
alleleflux run --config /home/su2806/AlleleFlux-dev/docs/source/examples/example_data_single/config_with_bams_single.yml
```
