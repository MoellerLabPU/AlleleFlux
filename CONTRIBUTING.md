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

## Code Organization

```
alleleflux/
├── scripts/
│   ├── analysis/        # Core analysis (profiling, scoring, outliers)
│   ├── preprocessing/   # Data preparation (metadata, QC, eligibility)
│   ├── statistics/      # Statistical tests (LMM, CMH, two-sample)
│   ├── accessory/       # Helper utilities (MAG mapping, coverage stats)
│   ├── utilities/       # Shared functions and logging
│   ├── visualization/   # Plotting and visualization tools
│   └── evolution/       # Evolutionary analysis (dN/dS)
└── smk_workflow/        # Snakemake workflow files
```

## Contributing Guidelines

### Code Style

- Follow PEP 8 style guidelines
- Use type hints for function parameters and return values
- Add docstrings to all public functions and classes
- Use `logging` instead of `print()` for output

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
4. Add documentation in `docs/source/reference/cli_reference.rst`

### Testing

Run tests with pytest:

```bash
# Run all tests
pytest tests/

# Run tests for a specific module
pytest tests/statistics/

# Run with coverage
pytest tests/ --cov=alleleflux --cov-report=html
```

### Documentation

Documentation is built with Sphinx:

```bash
cd docs
make html
# Open build/html/index.html in a browser
```

When adding new features:
- Update relevant RST files in `docs/source/`
- Add docstrings to new functions
- Include usage examples in CLI reference

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
