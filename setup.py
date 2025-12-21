"""Backward-compatible setuptools entrypoint.

This repository is configured via `pyproject.toml` (PEP 621 + setuptools).
Keeping a tiny `setup.py` allows legacy tooling to keep working without
duplicating or accidentally overriding the pyproject configuration.
"""

from setuptools import setup

setup()
