# AlleleFlux Documentation

```{image} https://readthedocs.org/projects/alleleflux/badge/?version=latest
:alt: Documentation Status
:target: https://alleleflux.readthedocs.io/en/latest/?badge=latest
```

## Overview

Microbial communities are continually exposed to dynamic environmental conditions, which can drive rapid genetic adaptations. Monitoring changes in allele frequencies can provide critical insights into the selective pressures and evolutionary mechanisms shaping microbial communities, but quantifying and comparing allele-frequency shifts in the microbiota has remained challenging due to a lack of bioinformatics tools designed for longitudinal shotgun (meta)genomic datasets.

To address this, we developed **AlleleFlux**, a bioinformatic tool that analyzes site-specific allele frequency fluctuations and dN/dS ratios in microbial genomes across multiple timepoints and groups. AlleleFlux tests for deterministic allele frequency changes using two key statistics:

1. A **parallelism score** that tests for parallel changes in allele frequencies over time across replicates in each group.
2. A **divergence score** that quantifies the degree of allele-frequency divergence between groups.
3. A **dN/dS score** that calculates the ratio of non-synonymous to synonymous substitutions, providing insights into selective pressures acting on genes..

These scores enable direct comparisons of evolutionary dynamics across different taxa, genomes, and even genes. By comparing scores across the genome, AlleleFlux can identify genes exhibiting exceptionally high parallelism or divergence, indicative of strong selection. Hence, AlleleFlux provides a robust framework for understanding the genomics of adaptation in microbes and the selective forces shaping their genetic diversity.

```{toctree}
:caption: Getting Started
:maxdepth: 1

getting_started/overview
getting_started/installation
getting_started/quickstart
```

```{toctree}
:caption: Usage Guide
:maxdepth: 1

usage/input_preparation
usage/running_workflow
usage/visualization_guide
usage/interpreting_results
usage/dnds_analysis
```

```{toctree}
:caption: Examples
:maxdepth: 1

examples/tutorial
examples/use_cases
```

```{toctree}
:caption: Reference
:maxdepth: 1

reference/cli_reference
reference/configuration
reference/inputs
reference/outputs
```

# Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
