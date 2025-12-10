"""
AlleleFlux: Identify genomic targets of natural selection in bacterial communities.

AlleleFlux is a Snakemake-based workflow for analyzing allele frequency changes
in metagenomic time-series data to detect parallel evolution.

Example usage:
    # CLI
    $ alleleflux run --config config.yml

    # Python
    >>> import alleleflux
    >>> print(alleleflux.__version__)
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("AlleleFlux")
except PackageNotFoundError:
    __version__ = "unknown"

# Public API exports
from alleleflux.workflow import execute_workflow, get_snakefile, validate_config

__all__ = [
    "__version__",
    "execute_workflow",
    "validate_config",
    "get_snakefile",
]
