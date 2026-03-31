# Overview

## What is AlleleFlux?

AlleleFlux is a bioinformatics toolkit for analyzing allele frequency changes in metagenomic data. It quantifies evolutionary dynamics across timepoints and experimental groups, providing statistical methods and scoring systems to identify genes under natural selection in microbial populations.

AlleleFlux was developed to fill a gap in the metagenomics toolkit: while many tools exist for taxonomic profiling and functional annotation, few provide the statistical framework needed to detect **parallel evolution** and **directional selection** from longitudinal shotgun metagenomic data.

## Key Features

### Analysis Capabilities

- **Single-timepoint and longitudinal** study designs
- **Parallelism scoring** -- detect consistent allele changes across replicates within a group
- **Divergence scoring** -- quantify allele frequency differences between experimental groups
- **Outlier gene detection** -- identify genes with exceptionally high scores (putative targets of selection)
- **dN/dS ratio calculation** -- measure selection pressure via non-synonymous to synonymous substitution rates (Nei-Gojobori 1986)
- **Taxonomic aggregation** -- summarize scores from phylum to species level using GTDB-Tk classifications

### Statistical Methods

| Test | Use Case | Data Type |
|------|----------|-----------|
| Two-sample unpaired (t-test, Mann-Whitney U) | Independent group comparison | Single or longitudinal |
| Two-sample paired (paired t-test, Wilcoxon) | Matched sample comparison | Longitudinal |
| Single-sample (one-sample t-test) | Within-group deviation from reference | Longitudinal |
| Linear Mixed Models (LMM) | Repeated measures with random effects | Longitudinal |
| Cochran-Mantel-Haenszel (CMH) | Stratified categorical analysis | Longitudinal |

### Quality Control

- Coverage breadth-based sample filtering
- MAG eligibility determination per statistical test
- Position-level preprocessing with configurable thresholds
- FDR correction (Benjamini-Hochberg) with optional per-MAG grouping

## How It Works

AlleleFlux implements a unified Snakemake pipeline that orchestrates the complete analysis in two major steps. The following diagram shows the full data flow from input files through final results:

```text
BAM files + FASTA + Metadata
    |
    v
Step 1: Profiling -> Metadata -> QC -> Eligibility (checkpoint)
    |
    v
Step 2: Allele Analysis -> Preprocessing -> Statistical Tests -> Scoring -> Gene Analysis -> dN/dS
                                                |
                                 +--------------+--------------+
                                 v              v              v
                           Between-groups  Within-group    CMH tests
                           (unpaired/LMM)  (single/paired) (stratified)
                                 |              |              |
                                 +--------------+--------------+
                                                v
                                      P-value Summary (FDR)
                                                v
                                       Scores + Outliers
                                                v
                                     Visualization (optional)
```

**Step 1: Profiling and Quality Control**

1. **Profiling** -- Extract per-position allele counts (A, C, G, T) from BAM files for each MAG using pysam. Each position records the total coverage and individual nucleotide counts.
2. **Metadata Generation** -- Build per-MAG sample metadata files linking profile paths to the experimental design (sample IDs, groups, timepoints).
3. **Quality Control** -- Filter samples by coverage breadth and depth thresholds. Samples that fail QC are excluded from downstream analysis.
4. **Eligibility** -- Determine which MAGs have sufficient passing samples for each statistical test type. This checkpoint gates entry into Step 2, ensuring only well-supported MAGs proceed.

**Step 2: Analysis and Scoring**

5. **Allele Frequency Analysis** -- Compute allele frequency changes across samples and timepoints for each eligible MAG.
6. **Preprocessing** -- Apply position-level filters (e.g., minimum coverage, minimum allele count) and prepare data structures for statistical testing.
7. **Statistical Testing** -- Apply the appropriate tests based on experimental design: between-group tests (unpaired t-test, Mann-Whitney U, LMM), within-group tests (paired t-test, Wilcoxon, single-sample), and CMH tests for stratified categorical analysis.
8. **P-value Summary** -- Aggregate test results across all positions and apply FDR correction (Benjamini-Hochberg) to control for multiple testing.
9. **Scoring** -- Calculate parallelism and divergence scores at MAG, taxon, and gene levels using the corrected p-values.
10. **Outlier Detection** -- Identify genes with significantly elevated scores using binomial and Poisson tests, flagging putative targets of selection.
11. **dN/dS Analysis** -- Calculate evolutionary rates for genes using the Nei-Gojobori method, distinguishing synonymous from non-synonymous substitutions.
12. **Visualization** -- Generate allele frequency trajectory plots (line, box, violin) for genes and positions of interest.

The workflow automatically parallelizes across samples and MAGs, handles checkpointing and restarts, and supports both local execution and HPC clusters (SLURM).

## Supported Study Designs

AlleleFlux supports two primary study designs, configured via the `data_type` parameter in the configuration file:

### Single-Timepoint (`data_type: "single"`)

Compare allele frequencies between experimental groups at a single timepoint. This design is appropriate when:

- You have cross-sectional data with two or more groups (e.g., treatment vs. control)
- No longitudinal sampling was performed
- You want to identify positions where allele frequencies differ significantly between groups

**Available tests:** Two-sample unpaired (t-test, Mann-Whitney U), LMM

**Scores produced:** Divergence scores quantifying between-group allele frequency differences

### Longitudinal (`data_type: "longitudinal"`)

Track allele frequency changes over multiple timepoints within and between groups. This design is appropriate when:

- You have repeated samples from the same subjects/replicates over time
- You want to detect parallel evolutionary trajectories across replicates
- You want to compare the direction and magnitude of allele frequency changes between groups

**Available tests:** Two-sample unpaired, two-sample paired, single-sample, LMM, CMH

**Scores produced:** Parallelism scores (within-group consistency), divergence scores (between-group differences), and dN/dS ratios

**Configuration:** Longitudinal analyses require specifying `timepoints_combinations` in the config -- a list of timepoint pairs (e.g., pre vs. post) with a `focus` timepoint indicating the direction of change.

## Key Concepts

### MAGs (Metagenome-Assembled Genomes)

AlleleFlux operates on MAGs -- draft genomes reconstructed from metagenomic assemblies. Each MAG represents a microbial population, and allele frequency analysis is performed independently per MAG. A MAG-to-contig mapping file links assembled contigs to their parent MAGs.

### Parallelism Scores

Parallelism scores measure the consistency of allele frequency changes across biological replicates within a group. A high parallelism score at a genomic position indicates that the same allele increased (or decreased) in frequency across multiple independent replicates, suggesting deterministic rather than stochastic change. Scores are computed at the position level and aggregated to gene, MAG, and taxonomic levels.

### Divergence Scores

Divergence scores quantify how much allele frequencies differ between experimental groups (e.g., treatment vs. control). High divergence at a position means the groups have shifted in opposite directions or to different magnitudes, indicating group-specific selective pressures. Like parallelism scores, these are aggregated across genomic levels.

### dN/dS Ratios

The ratio of non-synonymous (dN) to synonymous (dS) substitution rates provides insight into the type of selection acting on a gene:

- **dN/dS > 1** -- Positive (diversifying) selection: non-synonymous changes are favored
- **dN/dS = 1** -- Neutral evolution: no selective pressure distinguishing synonymous from non-synonymous changes
- **dN/dS < 1** -- Purifying (negative) selection: non-synonymous changes are removed

AlleleFlux calculates dN/dS using the Nei-Gojobori method, comparing allele frequencies between timepoints to estimate substitution rates.

### Eligibility

Not all MAGs have sufficient data for every statistical test. The eligibility system evaluates each MAG against test-specific criteria (e.g., minimum number of samples per group, paired samples across timepoints) and produces eligibility tables that gate which MAGs enter which analyses. This prevents unreliable results from underpowered tests.

## Architecture

```text
                            AlleleFlux Pipeline
 ┌─────────────────────────────────────────────────────────────────┐
 │                                                                 │
 │  STEP 1: Profiling & QC                                        │
 │  ━━━━━━━━━━━━━━━━━━━━━━                                        │
 │  BAM files ──► Profile MAGs ──► Generate Metadata ──► QC       │
 │                                                       │        │
 │                                              Eligibility Table  │
 │                                                       │        │
 │  STEP 2: Analysis & Scoring                           ▼        │
 │  ━━━━━━━━━━━━━━━━━━━━━━━━━━                                    │
 │  Allele Frequency Analysis ──► Preprocessing ──► Tests         │
 │                                                  │             │
 │                              ┌───────────────────┼──────┐      │
 │                              ▼                   ▼      ▼      │
 │                         Two-sample            LMM     CMH      │
 │                        (paired/unpaired)                       │
 │                              │                   │      │      │
 │                              └───────────────────┼──────┘      │
 │                                                  ▼             │
 │                              P-value Summary (FDR correction)  │
 │                                                  │             │
 │                              ┌───────────────────┼──────┐      │
 │                              ▼                   ▼      ▼      │
 │                           Scores            Gene     dN/dS     │
 │                         (MAG/Taxa)         Scores   Analysis   │
 │                                              │                 │
 │                                              ▼                 │
 │                                      Outlier Detection         │
 │                                              │                 │
 │                                              ▼                 │
 │                                      Visualization             │
 │                                   (Allele Trajectories)        │
 └─────────────────────────────────────────────────────────────────┘
```

## Use Cases

AlleleFlux is designed for:

- **Experimental evolution studies** -- Track microbial adaptation under controlled conditions
- **Longitudinal microbiome studies** -- Monitor allele frequency changes during diet interventions, antibiotic treatment, disease progression, or environmental change
- **Comparative group analysis** -- Identify differentially evolving genes between treatment and control groups
- **Selection target discovery** -- Find genes under strong positive or purifying selection in metagenomic data

For example workflows, see [Tutorial](../examples/tutorial.md), [Use Cases](../examples/use_cases.md), and [Interpreting Results](../usage/interpreting_results.md).

## Requirements

- Python >= 3.9
- R (for CMH tests via rpy2)
- Snakemake >= 8.0
- See [Installation](installation.md) for complete dependency information

## Next Steps

- [Installation](installation.md) -- Install AlleleFlux
- [Quickstart](quickstart.md) -- Run your first analysis
- [Input Preparation](../usage/input_preparation.md) -- Prepare your data
