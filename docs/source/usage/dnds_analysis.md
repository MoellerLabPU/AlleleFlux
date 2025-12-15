# dN/dS Analysis Guide

This guide covers dN/dS ratio analysis in AlleleFlux using the Nei-Gojobori (1986) method with path averaging.

## Overview

The **dN/dS ratio** (ω) measures evolutionary selection by comparing:
\- **dN**: Rate of non-synonymous substitutions (change amino acid)
\- **dS**: Rate of synonymous substitutions (silent mutations)

**Interpretation:**

```{eval-rst}
.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - dN/dS
     - Selection Mode
   * - < 1
     - **Purifying selection** – deleterious mutations removed
   * - ≈ 1
     - **Neutral evolution** – no selective pressure
   * - > 1
     - **Positive selection** – advantageous mutations favored
```

## Methodology

### Nei-Gojobori Method with Path Averaging

AlleleFlux implements the **NG86 method** to calculate dN/dS ratios from statistically significant sites identified by the workflow.

#### Key Features

1. **Potential Sites Calculation**
   \- For each codon, determine how many possible single-nucleotide changes are synonymous vs. non-synonymous
   \- Results in fractional S and N values (e.g., ATG: S=1.0, N=8.0)
2. **Path Averaging for Multi-Position Changes**
   \- When multiple positions change in a codon (k=2 or k=3), enumerate all possible mutational pathways
   \- Average S/NS classifications across paths
   \- Exclude pathways through intermediate stop codons
   \- Example: AAA→TTG has 6 pathways (3!), each potentially classifying mutations differently
3. **Codon-Level Analysis**
   \- Pre-computed cache for all 64 codons (potential sites)
   \- Pre-computed cache for all 4,096 codon pairs (path-averaged counts)
   \- O(1) lookups during analysis for efficiency

**Formula:**

```text
For each gene:
  Total_S = Σ(potential S sites for each codon)
  Total_N = Σ(potential N sites for each codon)

  Observed_S = Σ(path-averaged S count for changed codons)
  Observed_N = Σ(path-averaged NS count for changed codons)

  pS = Observed_S / Total_S
  pN = Observed_N / Total_N

  dN/dS = pN / pS
```

### Workflow Integration

dN/dS analysis follows statistical testing:

```text
Statistical tests → p_value_summary.py → significant sites
                                                   ↓
                  alleleflux-dnds-from-timepoints → dN/dS ratios
```

**Required inputs:**
1\. Significant sites file from `alleleflux-p-value-summary`
2\. Profile files for ancestral/derived timepoints
3\. Prodigal gene predictions
4\. Reference FASTA

## Usage

### Basic Command

```bash
alleleflux-dnds-from-timepoints \
    --input significant_sites_summary.tsv \
    --output output_dir/ \
    --mag_id Bacteroides_001 \
    --profiles_dir profiles/ \
    --prodigal_fasta genes.fna \
    --fasta reference.fasta \
    --ancestral_timepoint pre \
    --derived_timepoint post
```

### Key Parameters

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Parameter
     - Description
   * - ``--input``
     - Significant sites TSV from ``p_value_summary``
   * - ``--mag_id``
     - MAG identifier to analyze
   * - ``--p_value_column``
     - ``q_value`` (FDR-corrected) or ``min_p_value`` (default: ``q_value``)
   * - ``--p_value_threshold``
     - Significance cutoff (default: 0.05)
   * - ``--test_type``
     - Filter by test (e.g., ``two_sample_paired_tTest``)
   * - ``--group_analyzed``
     - Filter by group
   * - ``--ancestral_timepoint``
     - Earlier timepoint for sequence reconstruction
   * - ``--derived_timepoint``
     - Later timepoint for comparison
```

### Advanced Examples

**Filter by test type and group:**

```bash
alleleflux-dnds-from-timepoints \
    --input p_value_summary.tsv \
    --output dnds_results/ \
    --mag_id Lachnospira_002 \
    --test_type two_sample_paired_tTest \
    --group_analyzed treatment \
    --p_value_column q_value \
    --p_value_threshold 0.01 \
    --profiles_dir profiles/ \
    --prodigal_fasta genes.fna \
    --fasta ref.fasta \
    --ancestral_timepoint day0 \
    --derived_timepoint day30
```

**Process multiple MAGs:**

```bash
for mag in $(cat mag_list.txt); do
    alleleflux-dnds-from-timepoints \
        --input p_value_summary.tsv \
        --output dnds_results/ \
        --mag_id $mag \
        --profiles_dir profiles/ \
        --prodigal_fasta genes.fna \
        --fasta ref.fasta \
        --ancestral_timepoint pre \
        --derived_timepoint post
done
```

## Output Files

### Codon Events File

**Path:** `{output}/{mag}_codon_events_ng86.tsv.gz`

Codon-level substitution details with path-averaged counts.

```{eval-rst}
.. list-table::
   :widths: 20 15 65
   :header-rows: 1

   * - Column
     - Type
     - Description
   * - ``gene_id``
     - str
     - Gene identifier
   * - ``codon_number``
     - int
     - Position in gene (1-based)
   * - ``codon_before``
     - str
     - Ancestral codon
   * - ``codon_after``
     - str
     - Derived codon
   * - ``aa_before``, ``aa_after``
     - str
     - Amino acids
   * - ``k``
     - int
     - Number of changed positions (1, 2, or 3)
   * - ``fractional_S``
     - float
     - Path-averaged synonymous count
   * - ``fractional_NS``
     - float
     - Path-averaged non-synonymous count
   * - ``potential_S``, ``potential_N``
     - float
     - Expected S/N sites for this codon
```

### Gene Summary File

**Path:** `{output}/{mag}_gene_summary_ng86.tsv.gz`

Gene-level dN/dS ratios.

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Column
     - Description
   * - ``gene_id``
     - Gene identifier
   * - ``dN_dS``
     - dN/dS ratio
   * - ``pS``, ``pN``
     - Proportional synonymous/non-synonymous rates
   * - ``observed_S``, ``observed_N``
     - Total path-averaged counts
   * - ``potential_S``, ``potential_N``
     - Total expected S/N sites
   * - ``num_codons_changed``
     - Count of changed codons
   * - ``num_significant_sites``
     - Count of significant positions in gene
```

### MAG Summary File

**Path:** `{output}/{mag}_mag_summary_ng86.tsv.gz`

MAG-level dN/dS summary.

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Column
     - Description
   * - ``mag_id``
     - MAG identifier
   * - ``dN_dS``
     - Overall MAG dN/dS
   * - ``num_genes``
     - Total genes analyzed
   * - ``num_codons_changed``
     - Total changed codons across MAG
```

### Global Summary File

**Path:** `{output}/{mag}_global_summary_ng86.tsv`

Single-row summary of entire analysis.

Contains aggregate `dN_dS`, total sites analyzed, and overall statistics.

## Interpreting Results

### Understanding dN/dS Values

```{eval-rst}
.. list-table::
   :widths: 20 30 50
   :header-rows: 1

   * - dN/dS Range
     - Interpretation
     - Example Genes
   * - 0 - 0.3
     - Strong purifying selection
     - Essential housekeeping genes
   * - 0.3 - 0.7
     - Moderate constraint
     - Metabolic enzymes
   * - 0.7 - 1.3
     - Neutral/weak selection
     - Structural proteins
   * - 1.3 - 3.0
     - Positive selection
     - Antibiotic resistance genes
   * - > 3.0
     - Strong positive selection
     - Host-pathogen interaction genes
```

**NaN values:** Occur when no synonymous substitutions observed (pS = 0). Common in short genes or low-divergence comparisons.

### Understanding k-Values

The **k** value indicates how many positions changed in a codon:

- **k=1**: Single-position change (e.g., ATG→ATC) – straightforward S/NS classification
- **k=2**: Two positions changed (e.g., ATG→ACG) – 2! = 2 pathways averaged
- **k=3**: Three positions changed (e.g., ATG→CAA) – 3! = 6 pathways averaged

Higher k-values may indicate:
\- Longer time between samples
\- Stronger selection pressure
\- Less reliable dN/dS estimates (fewer intermediate observations)

### Biological Interpretation

**High dN/dS (>1) scenarios:**

1. **Antibiotic exposure**: Resistance gene mutations selected
2. **Nutrient shifts**: Metabolic adaptation to new resources
3. **Immune evasion**: Surface protein diversification
4. **Novel environment**: Rapid adaptation to colonization

**Low dN/dS (\<0.5) scenarios:**

1. **Core metabolism**: Essential enzymes under constraint
2. **Structural genes**: Protein folding requirements
3. **Information processing**: Ribosomal proteins, polymerases

## Common Use Cases

### Use Case 1: Genes Under Selection

**Goal:** Identify genes with dN/dS > 1 indicating positive selection.

**Workflow:**

```bash
# Run dN/dS analysis
alleleflux-dnds-from-timepoints --input p_val_summary.tsv --output dnds/ --mag_id MAG001 [...]

# Filter gene summary for high dN/dS
awk '$3 > 1' dnds/MAG001_gene_summary_ng86.tsv.gz | gzip -cd | sort -k3,3nr > genes_under_selection.tsv
```

**Interpretation:** Genes with dN/dS > 1 are candidates for adaptive evolution.

### Use Case 2: Comparing Selection Across MAGs

**Goal:** Compare selection patterns between microbial taxa.

**Workflow:**

```bash
# Process multiple MAGs
for mag in MAG001 MAG002 MAG003; do
    alleleflux-dnds-from-timepoints --input p_val_summary.tsv --mag_id $mag --output dnds/ [...]
done

# Extract MAG-level dN/dS
for mag in dnds/*_mag_summary_ng86.tsv.gz; do
    echo "$(basename $mag)"
    gzip -cd $mag | awk 'NR==2 {print $3}'
done
```

**Interpretation:** MAGs with higher dN/dS may be under stronger selective pressure.

### Use Case 3: Pathway-Level Analysis

**Goal:** Determine if functional pathways show coordinated selection.

**Workflow:**

```bash
# Annotate genes (e.g., with KEGG)
# Join gene summary with annotations

# Summarize by pathway
awk '{pathway=$5; dnds=$3; sum[pathway]+=dnds; count[pathway]++}
     END {for (p in sum) print p, sum[p]/count[p]}' annotated_genes.tsv
```

**Interpretation:** Pathways with elevated mean dN/dS suggest functional adaptation.

## Limitations and Considerations

### Low Coverage Regions

- Insufficient depth → unreliable allele frequency estimates
- Apply QC filters (`breadth_threshold` in config)

### Multiple Substitutions (k>1)

- Path averaging assumes all pathways equally likely
- High k-values may reduce accuracy
- Consider time between sampling points

### Stop Codons

- Intermediate stop codons in pathways are excluded (biologically unrealistic)
- Final stop codons allowed (represent valid termination)

### Sample Size

- Small sample sizes → less statistical power
- Ensure sufficient replicates for significance testing upstream

### Recombination

- dN/dS assumes vertical inheritance
- Horizontal gene transfer can inflate dN/dS estimates
- Consider phylogenetic context

### Time Between Samples

- Short intervals: fewer substitutions, many NaN values
- Long intervals: multiple substitutions (high k), path uncertainty
- Optimal: intermediate divergence (1-5% nucleotide difference)

## Troubleshooting

**Error: "Required columns not found"**

: Check significant sites file has columns: `mag_id`, `contig`, `position`, `gene_id`, `q_value` or `min_p_value`

**Warning: "Reference base mismatch"**

: Profile data reference doesn't match FASTA. Verify consistent reference used throughout workflow.

**Many NaN dN/dS values**
: - No synonymous substitutions observed (pS=0) – common in short genes or low divergence
  - Solution: Use longer time intervals or filter for genes with k≥2

**Script runs slowly**
: - Large number of significant sites → long runtime
  - Solution: Use stricter p-value threshold or parallelize across MAGs

### Validation

Sanity checks:

```bash
# Check potential sites make sense (should sum to ~3× gene length)
gzip -cd *_potential_sites_ng86.tsv.gz | awk '{s+=$2; n+=$3} END {print "S:", s, "N:", n, "Ratio:", n/s}'

# Expected N/S ratio ~2.5-3 for typical bacterial genes

# Check codon events have fractional counts summing correctly
gzip -cd *_codon_events_ng86.tsv.gz | awk '$9 != "" {print $9, $10}' | head
```

## Further Reading

**Primary References:**

- Nei & Gojobori (1986). *Mol Biol Evol* 3:418-426. [DOI:10.1093/oxfordjournals.molbev.a040410](https://doi.org/10.1093/oxfordjournals.molbev.a040410)
- Yang & Nielsen (2000). *Mol Biol Evol* 17:32-43. (Review of dN/dS methods)

**Related AlleleFlux Documentation:**

- {doc}`../usage/interpreting_results` – Allele frequency interpretation
- {doc}`../reference/outputs` – Output file specifications
- {doc}`../reference/cli_reference` – Command-line options

**External Tools:**

- [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) – Phylogeny-based dN/dS
- [HyPhy](http://hyphy.org/) – Hypothesis testing for positive selection
- [SNPGenie](https://github.com/chasewnelson/SNPGenie) – Within-population dN/dS
