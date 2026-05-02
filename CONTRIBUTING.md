# Contributing to AlleleFlux

Thank you for your interest in contributing to AlleleFlux! This document provides guidelines for contributing to the project.

## Getting Started

### Development Environment Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/MoellerLabPU/AlleleFlux.git
   cd AlleleFlux
   ```

2. **Create a conda environment:**
   ```bash
   mamba env create -f environment.yml
   mamba activate AlleleFlux
   ```

3. **Install in development mode:**
   ```bash
   pip install -e .
   ```

4. **Verify installation:**
   ```bash
   alleleflux --help
   ```

## Project Structure

```
AlleleFlux/
├── alleleflux/
│   ├── __init__.py
│   ├── cli.py                  # Main CLI (alleleflux command)
│   ├── workflow.py             # Snakemake orchestration
│   ├── scripts/
│   │   ├── analysis/           # Core analysis tools
│   │   ├── preprocessing/      # Data preparation
│   │   ├── statistics/         # Statistical tests
│   │   ├── evolution/          # dN/dS analysis
│   │   ├── accessory/          # Helper tools
│   │   ├── visualization/      # Plotting tools
│   │   └── utilities/          # Shared functions
│   └── smk_workflow/           # Snakemake workflow
│       ├── alleleflux_pipeline/
│       │   ├── Snakefile
│       │   ├── rules/          # Rule definitions
│       │   └── shared/         # Common functions
│       ├── config.template.yml
│       └── slurm_profile/
├── docs/                       # Sphinx documentation
├── tests/                      # Test suite
├── pyproject.toml
├── environment.yml
└── CHANGELOG.md
```

## Contributing Guidelines

### Code Style

- Follow PEP 8 style guidelines
- Use type hints for function parameters and return values
- Add docstrings to all public functions and classes
- Use `logging` instead of `print()` for output — see [Code Style](#code-style-guide) below for details

### CLI Script Template

All CLI scripts should follow this structure:

```python
#!/usr/bin/env python3
"""Brief description of the script."""

import argparse
import logging
import multiprocessing
from pathlib import Path

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Script description",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Add arguments...
    args = parser.parse_args()
    
    # Implementation...


if __name__ == "__main__":
    main()
```

### Adding a New CLI Command

1. Create your script in the appropriate subdirectory under `alleleflux/scripts/`
2. Add the entry point to `pyproject.toml`:
   ```toml
   [project.scripts]
   alleleflux-your-command = "alleleflux.scripts.category.your_script:main"
   ```
3. Reinstall: `pip install -e .`
4. Add documentation in `docs/source/reference/cli_reference.md`

## Code Style Guide

### Logging

Always use the `logging` module — never use `print()` for status output. Call `setup_logging()` exactly once in each script's `main()` function:

```python
from alleleflux.scripts.utilities.logging_config import setup_logging

setup_logging()  # Call ONCE in main()
logger = logging.getLogger(__name__)
```

Do not reconfigure logging after `setup_logging()` — it handles all setup centrally. The configured format is:

```
[%(asctime)s %(levelname)s] %(name)s: %(message)s
```

### File Paths

Use `pathlib.Path` for all file operations instead of string concatenation or `os.path.join()`:

```python
from pathlib import Path

output_dir = Path(args.output_dir)
output_file = output_dir / f"{mag_id}_results.tsv"

if not output_file.exists():
    logger.warning(f"File not found: {output_file}")
```

### Argument Parsing

Use `argparse` with `ArgumentDefaultsHelpFormatter` so that default values are automatically shown in `--help` output:

```python
parser = argparse.ArgumentParser(
    description="Script description",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
```

### Memory Optimization

Use categorical dtypes for grouping columns and explicit dtypes when reading data:

```python
# Categorical dtypes for low-cardinality columns
df["group"] = df["group"].astype("category")
df["time"] = df["time"].astype("category")

# Explicit dtype specifications when reading
dtype = {
    "contig": str,
    "position": int,
    "total_coverage": float,
    **{nuc: "int32" for nuc in ["A", "C", "G", "T"]}
}
df = pd.read_csv(file, sep="\t", dtype=dtype)
```

## Testing

### Test Directory Structure

The test directory mirrors the source layout under `alleleflux/scripts/`:

```
tests/
├── __init__.py
├── accessory/              # Tests for alleleflux/scripts/accessory/
│   ├── __init__.py
│   ├── test_coverage_and_allele_stats.py
│   └── test_positions_qc.py
├── analysis/               # Tests for alleleflux/scripts/analysis/
│   ├── __init__.py
│   └── test_profile_mags.py
├── evolution/              # Tests for alleleflux/scripts/evolution/
│   ├── __init__.py
│   ├── test_dnds_from_timepoints.py
│   └── mock_data/          # Test fixtures for dN/dS tests
├── preprocessing/          # Tests for alleleflux/scripts/preprocessing/
│   ├── __init__.py
│   ├── test_p_value_summary.py
│   └── test_quality_control.py
├── statistics/             # Tests for alleleflux/scripts/statistics/
│   ├── __init__.py
│   ├── test_LMM.py
│   ├── test_two_sample_paired.py
│   └── test_two_sample_unpaired.py
├── utilities/              # Tests for alleleflux/scripts/utilities/
│   └── test_qc_metrics.py
└── visualization/          # Tests for alleleflux/scripts/visualization/
    └── test_terminal_nucleotide_analysis.py
```

### Running Tests

```bash
# Run all tests
pytest tests/

# Run all tests with verbose output and short tracebacks
pytest tests/ -v --tb=short

# Run tests for a specific module
pytest tests/statistics/

# Run a single test file
pytest tests/statistics/test_LMM.py

# Run a specific test function
pytest tests/statistics/test_LMM.py::test_function_name

# Run with coverage report
pytest tests/ --cov=alleleflux --cov-report=html
```

### Adding Tests for New Features

When adding a new script or feature, create a corresponding test file:

1. Create a test file in the matching subdirectory under `tests/`. For example, if you add `alleleflux/scripts/analysis/new_feature.py`, create `tests/analysis/test_new_feature.py`.
2. Ensure the subdirectory has an `__init__.py` file.
3. Name test functions with the `test_` prefix so pytest discovers them.
4. Place test fixtures (mock data files) in a `mock_data/` subdirectory alongside the tests.
5. Test edge cases: empty DataFrames, single-sample groups, missing files, zero coverage positions, and NaN handling.

Example test structure:

```python
import pytest
import pandas as pd
from pathlib import Path
from alleleflux.scripts.analysis.new_feature import process_data

def test_process_data_basic():
    """Test basic processing with valid input."""
    df = pd.DataFrame({"contig": ["c1"], "position": [100], "A": [10], "C": [5], "G": [3], "T": [2]})
    result = process_data(df)
    assert not result.empty
    assert "score" in result.columns

def test_process_data_empty():
    """Test that empty input is handled gracefully."""
    df = pd.DataFrame()
    result = process_data(df)
    assert result is None or result.empty
```

## Documentation

Documentation is built with [Sphinx](https://www.sphinx-doc.org/) using the [Furo](https://pradyunsg.me/furo/) theme and is hosted on [ReadTheDocs](https://alleleflux.readthedocs.io/).

### Building Docs Locally

```bash
cd docs
make html
# Open build/html/index.html in a browser
```

For live-reloading during development:

```bash
cd docs
make livehtml
```

### File Format

Documentation pages use **MyST Markdown** (`.md`), not reStructuredText (`.rst`). The MyST parser is configured in `docs/source/conf.py` with these extensions enabled:

- `colon_fence` — allows `::: directive` syntax
- `deflist` — enables definition lists

### Adding New Documentation Pages

1. Create a new `.md` file in the appropriate subdirectory under `docs/source/` (e.g., `docs/source/usage/new_guide.md`).
2. Add the file to the relevant `toctree` directive in `docs/source/index.md`:
   ````markdown
   ```{toctree}
   :caption: Usage Guide
   :maxdepth: 1

   usage/existing_page.md
   usage/new_guide.md
   ```
   ````
3. Use MyST Markdown syntax for directives and cross-references. For example:
   ````markdown
   ```{note}
   This is an admonition.
   ```
   ````
4. Add docstrings to new Python functions — Sphinx `autodoc` will pick them up automatically.

### ReadTheDocs Deployment

The documentation is automatically built and deployed on ReadTheDocs when changes are pushed. Configuration is in `.readthedocs.yaml` at the repository root. No manual deployment is needed.

## Snakemake Rules

The Snakemake workflow lives in `alleleflux/smk_workflow/alleleflux_pipeline/`. Rules are organized into individual `.smk` files under the `rules/` directory, with shared helper functions in `shared/common.smk`.

### Existing Rule Files

```
alleleflux/smk_workflow/alleleflux_pipeline/
├── Snakefile               # Main entry point, includes all rules
├── shared/
│   ├── common.smk          # Config parsing, wildcards, helper functions
│   └── dynamic_targets.smk # Dynamic target generation
└── rules/
    ├── profiling.smk
    ├── metadata.smk
    ├── quality_control.smk
    ├── eligibility.smk
    ├── preprocessing_eligibility.smk
    ├── allele_analysis.smk
    ├── scoring.smk
    ├── gene_analysis.smk
    ├── significance_between_groups.smk
    ├── significance_within_group.smk
    ├── p_value_summary.smk
    └── dnds_analysis.smk
```

### Adding a New Rule

1. **Create a rule file** in `alleleflux/smk_workflow/alleleflux_pipeline/rules/`:

   ```python
   # rules/your_analysis.smk

   rule your_analysis:
       input:
           metadata="{output_dir}/metadata/{mag_id}_metadata.tsv",
       output:
           results="{output_dir}/your_analysis/{timepoints}/{groups}/{mag_id}_results.tsv.gz",
       params:
           extra_param=config["analysis"]["your_param"],
       threads: config["resources"]["your_analysis"]["cpus"]
       resources:
           mem_mb=config["resources"]["your_analysis"]["memory"],
           runtime=config["resources"]["your_analysis"]["time"],
       log:
           "{output_dir}/logs/your_analysis/{timepoints}/{groups}/{mag_id}.log",
       shell:
           """
           alleleflux-your-command \
               --metadata {input.metadata} \
               --output {output.results} \
               --cpus {threads} \
               --param {params.extra_param} \
               2>&1 | tee {log}
           """
   ```

2. **Include the rule file** in the main `Snakefile`:

   ```python
   include: "rules/your_analysis.smk"
   ```

3. **Add resource configuration** to `config.template.yml`:

   ```yaml
   resources:
     your_analysis:
       cpus: 4
       memory: 8000
       time: 60
   ```

4. **Wire up targets** — if the rule should run as part of the standard pipeline, add its output to the appropriate target list in `shared/dynamic_targets.smk` or the `rule all` definition.

5. **Use wildcard constraints** — if your rule introduces new wildcards, add constraints in `shared/common.smk` to prevent ambiguous matches.

6. **Use helper functions** from `shared/common.smk` for dynamic input resolution:
   - `get_sample_info()` — retrieve sample metadata
   - `get_mags_by_eligibility(timepoints, groups, test_type)` — get eligible MAG IDs
   - `parse_metadata_for_timepoint_pairs()` — parse timepoint/group combinations from config

## Pull Request Process

1. **Create a feature branch:**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following the guidelines above

3. **Run tests** to ensure nothing is broken:
   ```bash
   pytest tests/
   ```

4. **Commit with clear messages:**
   ```bash
   git commit -m "Add feature: brief description"
   ```

5. **Push and create a pull request:**
   ```bash
   git push origin feature/your-feature-name
   ```

6. **In the PR description:**
   - Describe what the changes do
   - Reference any related issues
   - Note any breaking changes

## Reporting Issues

When reporting bugs, please include:

- AlleleFlux version (`alleleflux --version`)
- Python version (`python --version`)
- Operating system
- Steps to reproduce the issue
- Error messages (full traceback)
- Minimal example data if applicable

## Questions?

- Check the [documentation](https://alleleflux.readthedocs.io/)
- Open a [GitHub issue](https://github.com/MoellerLabPU/AlleleFlux/issues)
- Ask on [DeepWiki](https://deepwiki.com/MoellerLabPU/AlleleFlux)

Thank you for contributing!
