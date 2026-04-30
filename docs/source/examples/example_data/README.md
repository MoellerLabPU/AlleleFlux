# AlleleFlux Example Data

Synthetic dataset for testing AlleleFlux on a normal desktop computer.

For full instructions, expected output, run times, and file format details see
the [tutorial](../tutorial.md).

## Quick start

```bash
conda activate alleleflux
cd docs/source/examples/example_data

# Option A — pre-generated profiles (recommended, ~94 s on 4 cores)
alleleflux run --config config_example.yml --threads 4

# Option B — full pipeline from BAMs, includes profiling (~74 s on 4 cores)
alleleflux run --config config_with_bams.yml --threads 4
```
