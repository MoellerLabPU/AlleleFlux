# AlleleFlux — Claude Code Instructions

## Environment Setup

**Always load the module and activate the conda environment before running any Python command, test, or pipeline step.**

```bash
module load anaconda3/2025.12
conda activate alleleflux
```

For one-off commands without activating the shell:

```bash
module load anaconda3/2025.12
conda run -n alleleflux python -m pytest tests/ -v
```

**Python executable:** `/home/su2806/.conda/envs/alleleflux/bin/python`

The system Python (`/usr/bin/python`) does not have any project packages. Without loading the module first, `conda` itself is not on PATH.

---

## Running Tests

```bash
module load anaconda3/2025.12

# All tests
conda run -n alleleflux python -m pytest tests/ -v --tb=short

# Single module
conda run -n alleleflux python -m pytest tests/analysis/test_profile_mags.py -v

# With coverage
conda run -n alleleflux python -m pytest tests/ --cov=alleleflux --cov-report=html
```

See `.github/instructions/tests.instructions.md` for full test conventions (class naming, mocking patterns, data setup).

---

## Project Overview

AlleleFlux is a Snakemake workflow for analyzing allele frequencies in metagenomic data. It detects parallel evolution in MAG (Metagenome-Assembled Genome) populations by tracking allele frequency changes across samples/timepoints.

```
Input BAM files
    → Step 1: profile_mags.py    — pileup, base counts, gene mapping
    → QC & eligibility filtering
    → Step 2: statistical tests + scoring
    → Output: scores, p-values, dN/dS, outliers
```

---

## Repository Structure

```
alleleflux/
├── scripts/
│   ├── analysis/        # profile_mags.py, allele_freq.py, scores.py,
│   │                    # gene_scores.py, taxa_scores.py, outliers_genes.py
│   ├── preprocessing/   # mag_metadata.py, quality_control.py,
│   │                    # eligibility_table.py, preprocess_*.py
│   ├── statistics/      # LMM.py, CMH.py, two_sample_paired/unpaired,
│   │                    # single_sample
│   ├── accessory/       # create_mag_mapping.py, coverage_and_allele_stats.py
│   ├── utilities/       # utilities.py, logging_config.py
│   ├── visualization/   # plotting tools
│   └── evolution/       # dnds_from_timepoints.py
└── smk_workflow/
    └── alleleflux_pipeline/
        ├── Snakefile              # Unified entry point
        ├── shared/common.smk      # Config parsing, resource helpers
        ├── shared/dynamic_targets.smk  # Checkpoint-aware target generation
        └── rules/                 # 13 rule modules (one per pipeline stage)

tests/                   # Mirror of alleleflux/scripts/ structure
notebooks/               # Validation and benchmarking notebooks
config.template.yml      # Pipeline config template
environment.yml          # Conda environment spec
```

---

## Configuration

The pipeline is driven by a YAML config with these top-level sections:

- **`input`**: BAM files, FASTA, MAG mapping, Prodigal genes, GTDB taxonomy
- **`output.root_dir`**: Base output directory
- **`analysis`**:
  - `data_type`: `"single"` or `"longitudinal"` — **critical**, changes statistical tests, output structure, and timepoint handling
  - `timepoints_combinations`: List of `{timepoint: [...], focus: ...}` dicts (longitudinal only)
  - `groups_combinations`: Group pairs to compare, e.g. `[["fat", "control"]]`
  - `use_lmm`, `use_significance_tests`, `use_cmh`: Feature flags
- **`quality_control`**: `min_sample_num`, `breadth_threshold`
- **`statistics`**: p-value thresholds, FDR settings, preprocessing filters
- **`resources`**: Per-rule CPUs, memory, time limits

Runtime configs are auto-saved as `{output_dir}/config_runtime_{timestamp}.yml`.

**Data type implications:**
- `"single"`: Compares groups at a single timepoint, uses unpaired tests
- `"longitudinal"`: Tracks changes across timepoints, uses paired tests and CMH stratified by replicate

---

## Key Implementation Patterns

### Logging

```python
from alleleflux.scripts.utilities.logging_config import setup_logging
setup_logging()          # Call ONCE in main() — never reconfigure after this
logger = logging.getLogger(__name__)
```

### CLI Script Template

```python
#!/usr/bin/env python3
import argparse, logging, multiprocessing, functools
from pathlib import Path
from tqdm import tqdm
from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)

def process_single_mag(mag_metadata_path, **kwargs):
    mag_id = Path(mag_metadata_path).stem.split("_metadata")[0]
    # ... return result or None on error

def main():
    setup_logging()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--rootDir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--cpus", type=int, default=multiprocessing.cpu_count())
    args = parser.parse_args()

    metadata_files = list(Path(args.rootDir).glob("*_metadata.tsv"))
    worker = functools.partial(process_single_mag, param=args.param)
    with multiprocessing.Pool(processes=args.cpus) as pool:
        for result in tqdm(pool.imap_unordered(worker, metadata_files), total=len(metadata_files)):
            pass
```

### MAG-Based Processing

- Scripts discover `*_metadata.tsv` files (one per MAG) in an input directory
- MAG ID extraction: `mag_id = filename.split("_metadata")[0]`
- Metadata TSV columns: `sample_id, file_path, group, time` (`group`/`time` optional)
- Output naming: `{mag_id}_{output_type}.tsv` or `.tsv.gz`

### Profile File Format

TSV with columns: `contig, position, ref_base, total_coverage, A, C, G, T, N, gene_id`  
(`mapq_scores` was removed in the samtools mpileup optimization.)

### Parallelization

```python
from functools import partial
from multiprocessing import Pool
from tqdm import tqdm

worker = partial(run_test_function, param1=val1, param2=val2)
with Pool(processes=cpus) as pool:
    results = [r for r in tqdm(pool.imap_unordered(worker, items), total=len(items)) if r is not None]
```

### Memory Optimization

```python
dtype = {"contig": str, "position": int, "total_coverage": float,
         **{nuc: "int32" for nuc in ["A", "C", "G", "T"]}}
df = pd.read_csv(file, sep="\t", dtype=dtype)
df["group"] = df["group"].astype("category")  # categorical for grouping columns
```

### R Integration (CMH)

```python
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
ro.globalenv["r_df"] = pivoted_df
ro.r("""
    library(tidyr)
    tab3d <- xtabs(count ~ group + nucleotide + replicate, data = r_df)
    cmh_result <- mantelhaen.test(tab3d, alternative="two.sided", correct=FALSE)$p.value
""")
p_value = ro.globalenv["cmh_result"][0]
```

### Common Utilities

```python
from alleleflux.scripts.utilities.utilities import (
    load_mag_mapping,           # Load contig → MAG mapping TSV
    extract_mag_id,             # Get MAG ID from contig name
    load_and_filter_data,       # Load profile + apply preprocessing filters
    build_contig_length_index,  # For coverage-weighted calculations
)
```

---

## Snakemake Workflow

Rules call AlleleFlux CLI commands (not Python directly) and use centralized resource helpers:

```python
rule profile_mags:
    input: ...
    output: ...
    threads: get_threads("profile_mags")
    resources:
        mem_mb=get_mem_mb("profile_mags"),
        time=get_time("profile_mags"),
    shell: "alleleflux-profile --bam_path {input.bam} --cpus {threads} ..."
```

Two checkpoint stages control DAG resolution:
1. **`eligibility_table`** — QC after Step 1, determines which MAGs proceed
2. **`preprocessing_eligibility_*`** — Filter before Step 2 scoring

See `.github/instructions/snakemake.instructions.md` for full conventions.

---

## Data Paths

Reference data and BAM/FASTA files are on cluster scratch storage:

```
/scratch/gpfs/AMOELLER/sidd/
```

These paths are only accessible from Della HPC nodes, not local machines.

---

## Error Handling Conventions

- Worker functions return `None` on error (filtered out during results collection)
- Log errors with MAG ID context: `logger.error(f"MAG {mag_id}: {e}")`
- Always validate file existence before opening:
  ```python
  if not Path(file_path).exists():
      logger.warning(f"File not found: {file_path}")
      return None
  ```
- Check for empty DataFrames before writing output

---

## Common Pitfalls

- **Don't** use `/usr/bin/python` — always use the conda environment
- **Don't** assume file existence — check and handle gracefully
- **Don't** modify global config in worker functions — use `functools.partial`
- **Don't** mix `"single"` and `"longitudinal"` data in the same analysis
- **Don't** hardcode group/timepoint names — read from config or metadata
- **Don't** call `setup_logging()` more than once
- **Don't** use `print()` for output — use `logger`

---

## Reference

- Tests: `.github/instructions/tests.instructions.md`
- Snakemake: `.github/instructions/snakemake.instructions.md`
- Environment: `.github/instructions/environment.instructions.md`
- Development: `CONTRIBUTING.md`
- Config template: `config.template.yml`
