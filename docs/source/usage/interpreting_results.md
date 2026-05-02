# Interpreting Results

This guide explains how to interpret the results produced by AlleleFlux.

## Output Structure

AlleleFlux organizes results by analysis type:

```text
output/
├── profiles/               # Per-sample allele counts
├── metadata/               # Per-MAG sample metadata
├── QC/                     # Quality control results
├── eligibility_table_*.tsv # MAG eligibility for tests
├── allele_analysis/        # Allele frequency analysis
├── significance_tests/     # Statistical test results
│   ├── two_sample_unpaired/
│   ├── two_sample_paired/
│   ├── single_sample/
│   ├── lmm/
│   └── cmh/
├── scores/
│   ├── intermediate/        # Per-MAG scores
│   └── processed/
│       ├── combined/        # Aggregated MAG/taxa scores
│       └── gene_scores/     # Gene-level scores
└── outliers/               # High-scoring outlier genes
```

## Quick Interpretation Guide

After running AlleleFlux, follow this decision tree to interpret your results:

1. **Check eligibility** → How many MAGs passed QC?
2. **Compare MAG scores** → Which MAGs show the strongest selection?
3. **Examine gene scores** → Which genes drive the MAG-level signal?
4. **Check outliers** → Are there genes under exceptionally strong selection?
5. **Validate with dN/dS** → Do outlier genes show positive selection (dN/dS > 1)?
6. **Cross-reference tests** → Are results consistent across statistical methods?

### Step-by-Step Checklist

Here's a practical checklist to work through after analysis:

- [ ] **Step 1: Check Eligibility** - Open `eligibility_table_*.tsv` and count how many MAGs passed QC filtering (have "true" in relevant test columns). If <10 MAGs eligible, lower QC thresholds or check coverage.
  
- [ ] **Step 2: Compare MAG Scores** - Load scores from `scores/processed/combined/MAG/` and identify top 10-20 MAGs by parallelism/divergence score. These are your primary candidates.
  
- [ ] **Step 3: Examine Gene Scores** - For top MAGs, inspect `scores/processed/gene_scores/{mag}_{test}_gene_scores_individual.tsv` to see which genes contribute most to the MAG signal.
  
- [ ] **Step 4: Check Outliers** - Review `outliers/{mag}_{test}_outlier_genes.tsv` files. Genes with multiple significant p-values (binomial, Poisson) are robust candidates.
  
- [ ] **Step 5: Validate with dN/dS** - If available, run `alleleflux-dnds --timepoint_data results/` on outlier genes. High scores with dN/dS > 1 indicate positive selection.
  
- [ ] **Step 6: Cross-Reference Tests** - Compare results across statistical methods (two_sample_paired, LMM, CMH). Genes significant in 2+ methods are more robust.

## Key Files

**1. Eligibility Table**

`eligibility_table_{timepoints}-{groups}.tsv` - Which MAGs qualify for each test based on coverage/samples

**2. Statistical Tests** (in `significance_tests/`)

Per-MAG files with p-values and test statistics:
\- `{mag}_two_sample_unpaired.tsv.gz` - Unpaired group comparisons
\- `{mag}_lmm.tsv.gz` - Linear mixed models
\- `{mag}_cmh.tsv.gz` - Cochran-Mantel-Haenszel tests

Key columns: `contig`, `position`, `gene_id`, `p_value_{test}`, `q_value_{test}`

**3. Scores** (in `scores/processed/combined/`)

- `scores_{test}-{tp}-{gr}-MAGs.tsv` - MAG-level parallelism/divergence scores
- `scores_{test}-{tp}-{gr}-{taxon}.tsv` - Taxonomic aggregations (phylum to species)

**4. Gene Scores** (in `scores/processed/gene_scores/`)

- `{mag}_{test}_gene_scores_individual.tsv` - Per-gene scores
- `{mag}_{test}_outlier_genes.tsv` - High-scoring genes under selection

## Score Interpretation

**Parallelism Score** (0-100%)
Measures consistent allele changes across replicates within a group. High scores → deterministic evolution (not random drift).

**Divergence Score** (0-100%)
Quantifies allele frequency differences between groups. High scores → differential selection between conditions.

**CMH Test**
Detects parallel allele changes across timepoints while controlling for individual variation. Particularly powerful for longitudinal studies.

## File Format Details

**Profile files** (`profiles/{sample}_{mag}_profiled.tsv.gz`):
`contig`, `position`, `ref_base`, `total_coverage`, `A`, `C`, `G`, `T`, `gene_id`

**Statistical test results** (`significance_tests/{test}/{mag}_{test}.tsv.gz`):
`contig`, `position`, `gene_id`, `p_value_{test}`, `q_value_{test}`

**Gene scores** (`scores/processed/gene_scores/{mag}_{test}_gene_scores_individual.tsv`):
`gene_id`, `total_sites`, `significant_sites`, `score_%`

**Outliers** (`outliers/{mag}_{test}_outlier_genes.tsv`):
`gene_id`, `gene_score_%`, `mag_score_%`, `p_value_binomial`, `p_value_poisson`

## Loading Results in Python

### Load and Display Top 10 MAGs by Parallelism Score

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load MAG-level scores
scores = pd.read_csv(
    "output/longitudinal/scores/processed/combined/MAG/"
    "scores_two_sample_paired-pre_post-treatment_control-MAGs.tsv",
    sep="\t"
)

# Top 10 MAGs by parallelism score
top_mags = scores.nlargest(10, "score_two_sample_paired_tTest (%)")
print("Top 10 MAGs by Parallelism Score:")
print(top_mags[["MAG_ID", "score_two_sample_paired_tTest (%)", 
                "total_sites_per_group_two_sample_paired_tTest", "phylum", "family"]])
```

### Plot Score Distribution

```python
# Visualize score distribution across all MAGs
fig, ax = plt.subplots(figsize=(10, 6))
scores["score_two_sample_paired_tTest (%)"].hist(bins=50, ax=ax, edgecolor="black")
ax.set_xlabel("Parallelism Score (%)")
ax.set_ylabel("Number of MAGs")
ax.set_title("Distribution of Parallelism Scores Across MAGs")
ax.axvline(scores["score_two_sample_paired_tTest (%)"].median(), 
           color="red", linestyle="--", label="Median")
plt.legend()
plt.tight_layout()
plt.savefig("score_distribution.png", dpi=150)
plt.show()

# Summary statistics
print(f"Mean parallelism score: {scores['score_two_sample_paired_tTest (%)'].mean():.2f}%")
print(f"Median parallelism score: {scores['score_two_sample_paired_tTest (%)'].median():.2f}%")
print(f"95th percentile: {scores['score_two_sample_paired_tTest (%)'].quantile(0.95):.2f}%")
```

### Compare Scores Across Different Tests

```python
# Load scores from different statistical tests
output_base = "output/longitudinal/scores/processed/combined/MAG/"
paired = pd.read_csv(
    f"{output_base}scores_two_sample_paired-pre_post-treatment_control-MAGs.tsv", 
    sep="\t"
)
lmm = pd.read_csv(
    f"{output_base}scores_lmm-pre_post-treatment_control-MAGs.tsv", 
    sep="\t"
)
cmh = pd.read_csv(
    f"{output_base}scores_cmh-pre_post-treatment_control-MAGs.tsv", 
    sep="\t"
)

# Merge on MAG_ID
merged = paired[["MAG_ID", "score_two_sample_paired_tTest (%)"]].merge(
    lmm[["MAG_ID", "score_lmm (%)"]],
    on="MAG_ID",
    how="outer"
).merge(
    cmh[["MAG_ID", "score_CMH (%)"]],
    on="MAG_ID",
    how="outer"
)

# Fill NaN for MAGs missing from certain tests
merged.fillna(0, inplace=True)

# Correlation between test results
print("Correlation between statistical tests:")
print(merged[["score_two_sample_paired_tTest (%)", "score_lmm (%)", "score_CMH (%)"]].corr())

# Visualize consistency
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].scatter(merged["score_two_sample_paired_tTest (%)"], merged["score_lmm (%)"], alpha=0.6)
axes[0].set_xlabel("Two-Sample Paired (%)"); axes[0].set_ylabel("LMM (%)")
axes[1].scatter(merged["score_two_sample_paired_tTest (%)"], merged["score_CMH (%)"], alpha=0.6)
axes[1].set_xlabel("Two-Sample Paired (%)"); axes[1].set_ylabel("CMH (%)")
axes[2].scatter(merged["score_lmm (%)"], merged["score_CMH (%)"], alpha=0.6)
axes[2].set_xlabel("LMM (%)"); axes[2].set_ylabel("CMH (%)")
plt.tight_layout()
plt.savefig("test_consistency.png", dpi=150)
plt.show()

# Identify robust MAGs (high scores in multiple tests)
threshold = 50  # 50th percentile
robust_mags = merged[
    (merged["score_two_sample_paired_tTest (%)"] > threshold) |
    (merged["score_lmm (%)"] > threshold) |
    (merged["score_CMH (%)"] > threshold)
]
print(f"Found {len(robust_mags)} MAGs with score > {threshold}% in at least one test")
```

### Load and Investigate Outlier Genes for a Specific MAG

```python
# Load outlier genes for a specific MAG
import gzip

mag_id = "Bacteroides_001"
outlier_file = (
    f"output/longitudinal/outliers/pre_post-treatment_control/"
    f"{mag_id}_two_sample_paired_outlier_genes.tsv.gz"
)

# Handle gzipped files
if outlier_file.endswith('.gz'):
    outliers = pd.read_csv(outlier_file, sep="\t", compression='gzip')
else:
    outliers = pd.read_csv(outlier_file, sep="\t")

print(f"Found {len(outliers)} outlier genes for {mag_id}")

# Display top outliers sorted by gene score
top_outliers = outliers.nlargest(20, "gene_score_%")[
    ["gene_id", "gene_score_%", "mag_score_%", "p_value_binomial", "p_value_poisson"]
]
print("\nTop 20 Outlier Genes:")
print(top_outliers)

# Filter by significance (both p-values < 0.05)
significant = outliers[
    (outliers["p_value_binomial"] < 0.05) & 
    (outliers["p_value_poisson"] < 0.05)
]
print(f"\nGenes significant in both tests: {len(significant)}")
print(significant[["gene_id", "gene_score_%", "p_value_binomial", "p_value_poisson"]])

# Visualize outlier gene scores
fig, ax = plt.subplots(figsize=(12, 6))
outliers_sorted = outliers.sort_values("gene_score_%", ascending=False).head(30)
ax.barh(range(len(outliers_sorted)), outliers_sorted["gene_score_%"])
ax.set_yticks(range(len(outliers_sorted)))
ax.set_yticklabels(outliers_sorted["gene_id"], fontsize=9)
ax.set_xlabel("Gene Score (%)")
ax.set_title(f"Top 30 Outlier Genes in {mag_id}")
plt.tight_layout()
plt.savefig(f"{mag_id}_outliers.png", dpi=150)
plt.show()
```

### Merge Taxonomic Information with Scores

```python
# Load MAG scores with taxonomic metadata
scores = pd.read_csv(
    "output/longitudinal/scores/processed/combined/MAG/"
    "scores_two_sample_paired-pre_post-treatment_control-MAGs.tsv",
    sep="\t"
)

# Group by phylum and calculate aggregated statistics
phylum_stats = scores.groupby("phylum").agg({
    "score_two_sample_paired_tTest (%)": ["mean", "max", "count"],
    "total_sites_per_group_two_sample_paired_tTest": "mean"
}).round(2)
phylum_stats.columns = ["mean_score", "max_score", "num_mags", "avg_sites"]
phylum_stats = phylum_stats.sort_values("mean_score", ascending=False)

print("Parallelism Scores by Phylum:")
print(phylum_stats)

# Similarly for family level
family_stats = scores.groupby("family").agg({
    "score_two_sample_paired_tTest (%)": ["mean", "max", "count"]
}).round(2)
family_stats.columns = ["mean_score", "max_score", "num_mags"]
family_stats = family_stats[family_stats["num_mags"] >= 3]  # Filter families with 3+ MAGs
family_stats = family_stats.sort_values("mean_score", ascending=False)

print("\nTop Families by Mean Parallelism Score:")
print(family_stats.head(10))

# Visualization
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Phylum-level scores
phylum_stats.head(10)["mean_score"].plot(kind="barh", ax=ax1, color="steelblue")
ax1.set_xlabel("Mean Parallelism Score (%)")
ax1.set_title("Top 10 Phyla by Mean Score")

# Family-level scores
family_stats.head(10)["mean_score"].plot(kind="barh", ax=ax2, color="coral")
ax2.set_xlabel("Mean Parallelism Score (%)")
ax2.set_title("Top 10 Families by Mean Score")

plt.tight_layout()
plt.savefig("taxonomic_scores.png", dpi=150)
plt.show()
```

## Loading Results in R

### Load MAG Scores with tidyverse

```r
library(tidyverse)
library(ggplot2)

# Load MAG-level scores
scores <- read_tsv("output/longitudinal/scores/processed/combined/MAG/scores_two_sample_paired-pre_post-treatment_control-MAGs.tsv")

# Display top MAGs by score
scores %>%
  arrange(desc(`score_two_sample_paired_tTest (%)`)) %>%
  select(MAG_ID, `score_two_sample_paired_tTest (%)`, total_sites_per_group_two_sample_paired_tTest, phylum, family) %>%
  head(10)
```

### Arrange by Score and Select Columns

```r
# More focused view with key columns
scores %>%
  select(MAG_ID, `score_two_sample_paired_tTest (%)`, phylum, family, class) %>%
  arrange(desc(`score_two_sample_paired_tTest (%)`)) %>%
  filter(`score_two_sample_paired_tTest (%)` > 50) %>%
  head(20)

# Save top MAGs to CSV
top_mags <- scores %>%
  arrange(desc(`score_two_sample_paired_tTest (%)`)) %>%
  slice_head(n = 20)
write_csv(top_mags, "top_20_mags.csv")
```

### Group by Taxonomic Family and Summarize

```r
# Aggregate statistics by family
family_summary <- scores %>%
  group_by(family) %>%
  summarise(
    mean_score = mean(`score_two_sample_paired_tTest (%)`, na.rm = TRUE),
    max_score = max(`score_two_sample_paired_tTest (%)`, na.rm = TRUE),
    median_score = median(`score_two_sample_paired_tTest (%)`, na.rm = TRUE),
    n_mags = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))

print(family_summary)

# Visualize top families
family_summary %>%
  filter(n_mags >= 3) %>%  # Filter families with 3+ MAGs
  slice_head(n = 15) %>%
  ggplot(aes(x = reorder(family, mean_score), y = mean_score)) +
  geom_col(fill = "steelblue", color = "black") +
  geom_text(aes(label = paste0("n=", n_mags)), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(x = "Family", y = "Mean Parallelism Score (%)", 
       title = "Top Families by Mean Parallelism Score") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("family_scores.png", width = 10, height = 6, dpi = 150)
```

### Aggregate by Phylum and Compare

```r
# Phylum-level aggregation
phylum_summary <- scores %>%
  group_by(phylum) %>%
  summarise(
    mean_score = mean(`score_two_sample_paired_tTest (%)`, na.rm = TRUE),
    sd_score = sd(`score_two_sample_paired_tTest (%)`, na.rm = TRUE),
    max_score = max(`score_two_sample_paired_tTest (%)`, na.rm = TRUE),
    n_mags = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))

print(phylum_summary)

# Boxplot by phylum
scores %>%
  ggplot(aes(x = reorder(phylum, `score_two_sample_paired_tTest (%)`, FUN = median), 
             y = `score_two_sample_paired_tTest (%)`)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
  coord_flip() +
  labs(x = "Phylum", y = "Parallelism Score (%)",
       title = "Score Distribution by Phylum") +
  theme_minimal()

ggsave("phylum_boxplot.png", width = 10, height = 6, dpi = 150)
```

### Join Multiple Test Results

```r
# Load multiple test results
paired <- read_tsv("output/longitudinal/scores/processed/combined/MAG/scores_two_sample_paired-pre_post-treatment_control-MAGs.tsv") %>%
  select(MAG_ID, score_paired = `score_two_sample_paired_tTest (%)`)

lmm <- read_tsv("output/longitudinal/scores/processed/combined/MAG/scores_lmm-pre_post-treatment_control-MAGs.tsv") %>%
  select(MAG_ID, score_lmm = `score_lmm (%)`)

cmh <- read_tsv("output/longitudinal/scores/processed/combined/MAG/scores_cmh-pre_post-treatment_control-MAGs.tsv") %>%
  select(MAG_ID, score_cmh = `score_CMH (%)`)

# Join all results
merged_scores <- paired %>%
  full_join(lmm, by = "MAG_ID") %>%
  full_join(cmh, by = "MAG_ID") %>%
  replace_na(list(score_paired = 0, score_lmm = 0, score_cmh = 0))

# Identify robust MAGs (high scores in multiple tests)
robust_mags <- merged_scores %>%
  mutate(num_high_scores = rowSums(across(starts_with("score_"), ~. > 50))) %>%
  filter(num_high_scores >= 2) %>%
  arrange(desc(num_high_scores), desc(score_paired))

print(robust_mags)

# Correlation between tests
merged_scores %>%
  select(starts_with("score_")) %>%
  cor()

# Visualize correlation
library(corrplot)
scores_matrix <- merged_scores %>% select(starts_with("score_")) %>% as.matrix()
corrplot(cor(scores_matrix), type = "lower", diag = FALSE)
ggsave("test_correlation.png", width = 8, height = 8, dpi = 150)
```

### Load and Analyze Gene-Level Data

```r
# Load gene scores for a specific MAG
mag_id <- "Bacteroides_001"
gene_scores <- read_tsv(
  glue::glue("output/longitudinal/scores/processed/gene_scores/{mag_id}_two_sample_paired_gene_scores_individual.tsv")
)

# Top genes by score
gene_scores %>%
  arrange(desc(score_%)) %>%
  select(gene_id, score_%, total_sites, significant_sites) %>%
  head(20)

# Load outlier genes
outliers <- read_tsv(
  glue::glue("output/longitudinal/outliers/pre_post-treatment_control/{mag_id}_two_sample_paired_outlier_genes.tsv")
)

# Filter for highly significant outliers
outliers_sig <- outliers %>%
  filter(p_value_binomial < 0.05, p_value_poisson < 0.05) %>%
  arrange(desc(gene_score_%))

print(outliers_sig)

# Visualization
outliers %>%
  arrange(desc(gene_score_%)) %>%
  slice_head(n = 25) %>%
  ggplot(aes(x = reorder(gene_id, gene_score_%), y = gene_score_%)) +
  geom_col(fill = "coral", color = "black") +
  geom_point(aes(color = p_value_binomial < 0.05), size = 3, alpha = 0.7) +
  coord_flip() +
  labs(x = "Gene ID", y = "Gene Score (%)", 
       title = glue::glue("Top 25 Outlier Genes in {mag_id}"),
       color = "Binomial p < 0.05") +
  theme_minimal()

ggsave(glue::glue("{mag_id}_outliers.png"), width = 12, height = 8, dpi = 150)
```

## Biological Validation

After identifying candidate MAGs and genes, validate findings with these approaches:

### Functional Annotation

- **KEGG Annotation** - Assign enzymes (EC numbers) and metabolic pathways to outlier genes. Look for enrichment in pathways relevant to your experimental conditions:
  - Antibiotic exposure → Check for antibiotic resistance genes (aminoglycoside modifying enzymes, β-lactamases, efflux pumps)
  - Nutrient limitation → Look for transport proteins, biosynthetic pathways
  - High temperature → Molecular chaperones, thermostability proteins
  
- **COG (Cluster of Orthologous Groups)** - Functional classification by COG family. Count significant genes per COG category (e.g., "Energy production" vs. "Defense mechanisms") and test for over-representation using Fisher's exact test.

- **Pfam Domains** - Protein family annotation identifies functional domains within outlier genes. Multiple domains may indicate multifunctional or regulatory proteins.

- **Tools for annotation**:
  ```bash
  # Using eggNOG-mapper (recommended for MAGs)
  emapper.py -i outlier_genes.fasta --output annotation_results -m diamond --data_dir /path/to/eggnog_data
  
  # Using InterProScan for detailed domain analysis
  interproscan.sh -i outlier_genes.fasta -f tsv -o annotation_results.tsv
  ```

### Cross-referencing with Literature

- **Curated Databases**:
  - **CARD (Comprehensive Antibiotic Resistance Database)** - If studying antibiotic resistance
  - **VFDB (Virulence Factor Database)** - For pathogenicity studies
  - **BioCyc** - Metabolic pathways and reactions
  - **TCDB** - Membrane transport proteins

- **Literature Mining**:
  - Search PubMed/Google Scholar for MAG taxa (genus/family) + experimental condition (e.g., "Bacteroides antibiotic resistance")
  - Compare high-scoring genes with previously published selection signatures in same organisms
  - Check if same functional categories repeatedly appear in independent studies

- **Example validation script**:
  ```python
  import pandas as pd
  
  # Load outliers and CARD annotations
  outliers = pd.read_csv("outlier_genes.tsv", sep="\t")
  card = pd.read_csv("card_annotations.tsv", sep="\t")
  
  # Find overlap
  resistance_genes = outliers.merge(card, on="gene_id", how="inner")
  print(f"Outlier genes matching CARD database: {len(resistance_genes)}")
  print(resistance_genes[["gene_id", "ARO", "phenotype"]])
  ```

### Checking Genomic Context

- **Operon Structure** - Genes functioning in the same pathway often cluster in operons. Check if multiple outlier genes are co-localized:
  ```bash
  # Extract genomic regions around outlier genes
  samtools faidx reference.fa contig:start-end > region.fa
  # Use BLAST or alignment tools to compare regions across samples
  ```

- **Mobile Genetic Elements** - Outlier genes near transposases, integrases, or IS elements suggest recent horizontal gene transfer and potential adaptive introgression:
  ```python
  # Check proximity to identified mobile elements
  outliers_nearby_mobile = outliers[
      outliers["distance_to_nearest_mobile_element"] < 5000  # 5 kb window
  ]
  print(f"Outliers near mobile elements: {len(outliers_nearby_mobile)}")
  ```

- **Genomic Islands** - Large segments with unusual GC content or codon usage may indicate acquired genes. Tools like SIGI-HMM or IslandViewer can identify these regions.

- **Conservation Analysis** - Highly conserved genes under selection are more likely functionally important:
  ```python
  # Calculate sequence identity with close relatives
  # High identity + positive selection → functionally constrained
  # Low identity + positive selection → recent adaptation
  ```

### Validating with dN/dS Analysis

- **Positive Selection Signature** - dN/dS > 1 indicates positive (diversifying) selection; dN/dS < 1 indicates purifying selection:
  ```bash
  # AlleleFlux includes dN/dS calculation from timepoint data
  alleleflux-dnds --timepoint_data step2_results/ --output_dir dnds_results/ --cpus 16
  ```

- **Interpreting Results**:
  - High score + dN/dS > 1 → Gene evolving adaptively under strong selection
  - High score + dN/dS < 1 → Adaptive structural changes without amino acid replacement (e.g., regulatory mutations)
  - High score + dN/dS ≈ 1 → Neutral evolution or recent selective sweep

- **Combined Analysis**:
  ```python
  import pandas as pd
  
  # Load outliers and dN/dS results
  outliers = pd.read_csv("outlier_genes.tsv", sep="\t")
  dnds = pd.read_csv("dnds_results.tsv", sep="\t")
  
  # Merge and filter
  combined = outliers.merge(dnds, on="gene_id", how="inner")
  positive_selection = combined[combined["dnds_ratio"] > 1]
  
  print(f"Genes with adaptive evolution signature: {len(positive_selection)}")
  print(positive_selection[["gene_id", "gene_score_%", "dnds_ratio"]])
  ```

### Experimental Validation

- **Allele Frequency Confirmation** - Use targeted sequencing or qPCR to independently validate predicted allele frequency changes:
  ```bash
  # Design primers flanking polymorphic sites
  # qPCR quantify allele-specific transcripts or genomic DNA
  # Expected: High frequency allele in enriched samples, low in control
  ```

- **Phenotypic Testing** - Culture representative strains carrying different alleles and measure phenotype under experimental conditions:
  - **Antibiotic resistance**: Growth inhibition assays with/without antibiotics
  - **Nutrient metabolism**: Growth rates on different substrates
  - **Stress tolerance**: Growth at extreme pH, temperature, or osmolarity

- **Longitudinal Tracking** - For time-series data, verify that predicted allele trajectories match experimental observations:
  ```python
  # Extract predicted vs. observed allele frequencies
  predicted = pd.read_csv("allele_predictions.tsv", sep="\t")
  observed = pd.read_csv("observed_frequencies.tsv", sep="\t")
  
  # Calculate R² between predicted and observed
  correlation = predicted.merge(observed, on=["sample", "gene", "position"]).corr()
  print(f"Prediction accuracy (R²): {correlation.iloc[0,1]**2:.3f}")
  ```

- **Heterologous Expression** - Express outlier genes in model organisms (e.g., *E. coli*) to assess function:
  - Clone gene into expression vector
  - Measure phenotypic effect (growth, metabolite production, resistance)
  - Compare wildtype vs. allelic variants if polymorphisms encode different proteins

## Analysis Workflow

**Step 1: Check Eligibility**

```bash
cat eligibility_table_pre_post-treatment_control.tsv
```

Identify MAGs with sufficient coverage for statistical tests.

**Step 2: Examine Scores**

```bash
# MAG-level scores
head scores_two_sample_unpaired-pre_post-treatment_control-MAGs.tsv

# Taxonomic aggregation (family level)
head scores_two_sample_unpaired-pre_post-treatment_control-family.tsv
```

Focus on MAGs/taxa with high parallelism or divergence scores.

**Step 3: Investigate Genes**

```bash
# Gene scores for a high-scoring MAG
head MAG123_two_sample_unpaired_gene_scores_individual.tsv

# Outlier genes
head MAG123_two_sample_unpaired_outlier_genes.tsv
```

Identify candidate genes under strong selection.

**Step 4: Compare Tests**

Check consistency across statistical approaches (two-sample, LMM, CMH). Genes significant in multiple tests are most robust.

**Step 5: Functional Analysis**

- Annotate outlier genes (KEGG, COG, Pfam)
- Check biological relevance to experimental conditions
- Consider genomic context (operons, mobile elements)

## Troubleshooting

**No results / empty files**
\- Check eligibility table: MAGs may not meet `min_sample_num` or `breadth_threshold`
\- Verify input file paths in configuration
\- Check log files in `logs/` directory

**Low scores across all MAGs**
\- Insufficient selective pressure or inappropriate timepoints
\- Try lowering `p_value_threshold` (e.g., 0.1 instead of 0.05)
\- Check if experimental conditions are strong enough

**Inconsistent results between tests**
\- LMM is sensitive to experimental design complexity
\- Two-sample tests affected by unbalanced groups
\- CMH best for detecting consistent directional changes
\- Use multiple tests for robust conclusions

**Missing gene IDs**
\- Ensure Prodigal predictions match reference FASTA contig names
\- Verify `prodigal_path` in configuration
\- Check gene FASTA headers match contig naming

For visualization of results, see [Visualization Guide](visualization_guide.md).
