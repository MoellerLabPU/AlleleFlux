# Overview

## What is AlleleFlux?

AlleleFlux analyzes allele frequency changes in metagenomic data, quantifying evolutionary dynamics across timepoints and experimental groups. It provides statistical methods and scoring systems to identify genes under selection in microbial populations.

## Key Features

**Analysis Capabilities**

- Single-timepoint and longitudinal study designs
- Parallelism and divergence scoring
- Outlier gene detection (genes under strong selection)
- dN/dS ratio calculation for evolutionary rate estimation
- Taxonomic aggregation (phylum to species level)

**Statistical Methods**

- Two-sample tests (paired/unpaired)
- Single-sample analysis
- Linear Mixed Models (LMM) for complex designs
- Cochran-Mantel-Haenszel (CMH) tests

**Quality Control**

- Coverage-based sample filtering
- MAG eligibility criteria
- Configurable preprocessing thresholds

## Workflow

AlleleFlux implements a unified Snakemake pipeline:

```text
Input Files              Profile & QC           Statistical Analysis
━━━━━━━━━━━              ━━━━━━━━━━━━           ━━━━━━━━━━━━━━━━━━━━
• BAM files              • Extract alleles      • Two-sample tests
• Reference FASTA        • Quality control      • LMM / CMH tests
• Metadata TSV           • Eligibility checks   • dN/dS calculation
• Gene annotations                ↓
                         ┌─────────────────┐
                         │ Scoring & Viz   │
                         ├─────────────────┤
                         │ • Parallelism   │
                         │ • Divergence    │
                         │ • Outliers      │
                         │ • Trajectories  │
                         └─────────────────┘
```

**Key Steps:**

1. **Profiling**: Extract allele frequencies from BAM files for each MAG
2. **Quality Control**: Filter samples by coverage breadth; determine MAG eligibility
3. **Statistical Testing**: Apply appropriate tests based on experimental design
4. **Scoring**: Calculate parallelism/divergence scores and identify outlier genes
5. **Visualization**: Generate allele trajectory plots and summary statistics

For detailed workflow documentation, see [Running the Workflow](../usage/running_workflow.md).

## Use Cases

AlleleFlux is designed for:

- Experimental evolution studies
- Longitudinal microbiome studies (diet interventions, disease progression, environmental change)
- Comparative analysis of adaptation across treatments
- Identifying genes under selection in metagenomic data

For example workflows, see [Tutorial](../examples/tutorial.md) and [Interpreting Results](../usage/interpreting_results.md).
