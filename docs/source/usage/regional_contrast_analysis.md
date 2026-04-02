# Regional Contrast Analysis

## Overview

Regional contrast analysis detects **genes or genomic windows with differential allele-frequency evolution** between treatment and control groups across paired hosts (e.g., biological replicates). It is designed to identify regions under divergent evolutionary pressures or parallel adaptation between conditions.

### Key Features

- **Longitudinal only**: Requires multiple timepoints to measure allele-frequency changes
- **Region-level analysis**: Operates on genes or sliding windows (not individual sites)
- **Paired-host comparison**: Contrasts evolution in treatment vs. control within each host, then tests consistency across hosts
- **No preprocessing required**: Works directly on raw allele-frequency changes, avoiding selection bias
- **Multiple comparison correction**: Applies Benjamini–Hochberg FDR control per region type and test direction
- **Multiple group support**: Fully supports analyzing multiple treatment-control pairs simultaneously

## When to Use

Use regional contrast analysis when:

- Your data is **longitudinal** (multiple timepoints per host)
- You have **paired hosts** (e.g., matched biological replicates, crossover designs)
- You want to detect regions with **consistent directional evolution differences** between conditions
- You're looking for **genes under adaptive pressure** or **genomic regions diverging between groups**

**Do NOT use** if:

- Data is single-timepoint only
- Unpaired hosts (single observation per condition)
- Looking for site-level significance (use statistical tests instead)

## Configuration

Regional contrast is configured in the `regional_contrast` section of `config.yml`:

```yaml
analysis:
  use_regional_contrast: true  # Enable/disable the entire module
  
  data_type: longitudinal      # MUST be longitudinal
  
  groups_combinations:
    - treatment: "high_fat"
      control: "control"
    - treatment: "high_fat"
      control: "standard"

regional_contrast:
  mode: "both"                      # "gene", "window", or "both"
  window_size: 1000                 # bp; used if mode includes "window"
  agg_method: "median"              # "median" (default), "mean", or "trimmed_mean"
  min_informative_sites: 5          # Min. sites per region (0 = no filter)
  min_informative_fraction: 0.0     # Min. fraction of region coverage (0.0 = no filter)
  use_fisher: true                  # Include Fisher combined p-values (exploratory)
  use_regional_contrast: true       # Overall enable/disable flag
```

### Configuration Parameters Explained

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `mode` | `both` | Analyze genes (`gene`), sliding windows (`window`), or both (`both`) |
| `window_size` | `1000` | Non-overlapping tile width in base pairs. Smaller (e.g., 500) = fine-grained; larger (e.g., 5000) = coarse-grained |
| `agg_method` | `median` | How to aggregate site scores per region. `median` = robust; `mean` = simple; `trimmed_mean` = robust mean |
| `min_informative_sites` | `5` | Minimum variable sites per region. Regions with fewer are excluded. Set to `0` to disable. **Note:** A region with all `site_score == 0` (stasis) still has full informative sites. |
| `min_informative_fraction` | `0.0` | Minimum coverage: `(sites_observed) / (region_length)`. Set to `0.0` to disable. |
| `use_fisher` | `true` | Compute Fisher combined p-values from percentiles (secondary/exploratory). Set to `false` to skip. |
| `use_regional_contrast` | `true` | Global enable/disable flag. Set to `false` to skip this analysis entirely. |

## Workflow Overview

Regional contrast follows a 10-step pipeline:

1. **Load & Validate**: Read allele-frequency-change file; filter to treatment/control groups
2. **Site Scores**: Compute Total-Variation distance at each site
3. **Region Definitions**: Map sites to genes and/or sliding windows
4. **Aggregation**: Summarize site scores per (host, group, region)
5. **Percentiles**: Rank region scores within each (host, group) to normalize host-specific variation
6. **Reshape**: Pivot so each (host, region) has treatment and control scores side-by-side
7. **Test Contrasts**: Run one-sided Wilcoxon and t-tests across hosts
8. **Fisher P-values** (optional): Combine percentile-derived empirical p-values
9. **FDR Correction**: Apply Benjamini–Hochberg correction per region type and test direction
10. **Write Outputs**: Produce per-host and summary result files

## Understanding the Analysis

### Total-Variation Distance (Site Score)

For each genomic site, the allele-frequency change is captured as a **Total-Variation (TV) distance**:

$$\text{site\_score} = \frac{1}{2} \times (|\Delta A| + |\Delta C| + |\Delta G| + |\Delta T|)$$

where $\Delta X$ is the change in frequency of nucleotide X from timepoint 1 to timepoint 2.

- **Range**: 0 (no change) to 1 (complete reversal)
- **Intuition**: Captures magnitude of allele-frequency change in a single number

### Region Aggregation & Normalization

Site scores within each region are aggregated (median/mean/trimmed-mean), then normalized to **genome-wide percentiles** within each (host, group):

- **Why percentiles?** Hosts differ in overall evolutionary rates. Percentiles account for host-specific variation, focusing on *relative* significance.
- **Example**: Gene A at 75th percentile in treatment but 30th in control → consistent positive contrast

### Treatment vs. Control Contrast

For each region and host:

$$\text{contrast} = \text{percentile\_treatment} - \text{percentile\_control}$$

**Question**: Do hosts consistently show higher (or lower) evolution in the treatment group?

### Across-Host Statistical Test

For each region, contrasts across hosts are tested with:

- **Wilcoxon signed-rank test** (primary, non-parametric): H₁ = median contrast > 0
- **One-sample t-test** (secondary, parametric): H₁ = mean contrast > 0

Results are:
- **One-sided**: Both test whether treatment evolves *more* than control
- **Independent**: Wilcoxon and t-test use separate α budgets; each has its own FDR correction
- **Directional**: Output columns explicitly label direction (`_greater` = treatment evolved more, `_less` = control evolved more)

## Multiple Group Combinations

AlleleFlux fully supports analyzing **multiple treatment-control pairs simultaneously**.

### Example Configuration

```yaml
groups_combinations:
  - treatment: "high_fat"
    control: "control"
  - treatment: "high_fat"
    control: "standard"
```

### What Happens

1. **Separate eligibility**: Each (treatment, control) pair gets its own eligibility table based on sample QC metrics
2. **Independent analysis**: Regional contrast runs separately for each pair
3. **Isolated outputs**: Results are written to separate directories:
   - `regional_contrast/regional_contrast_{timepoints}-high_fat_control/`
   - `regional_contrast/regional_contrast_{timepoints}-high_fat_standard/`
4. **Wildcard constraints**: Snakemake restricts `{treatment}` and `{control}` to valid values from your config, preventing cross-contamination

### Interpreting Multiple-Pair Results

- Compare significance across pairs: Which genes are differentially evolving in *both* treatment comparisons?
- Identify group-specific signals: Genes significant only in one pair but not the other suggest condition-specific adaptation
- Robustness check: High-confidence findings often appear in multiple pair comparisons

## Output Files

### Per-Host Region Table (`*_per_host_region.tsv.gz`)

One row per (host, region_id, region_type) tuple. Gzip-compressed due to size.

**Key columns:**
- `replicate`: Host/paired-unit ID
- `region_id`: Gene ID or window ID
- `region_type`: `"gene"` or `"window"`
- `region_score_treatment`, `region_score_control`: Aggregated scores
- `percentile_treatment`, `percentile_control`: Genome-wide percentile ranks
- `n_informative_sites_treatment`, `n_informative_sites_control`: Site counts
- `contrast`: `percentile_treatment - percentile_control`

### Region Summary Table (`*_region_summary.tsv`)

One row per region. Plain-text (uncompressed) for easy viewing.

**Key columns:**
- `region_id`, `region_type`, `contig`, `region_start`, `region_end`: Region metadata
- `n_hosts`: Number of hosts with data for this region
- `mean_contrast`, `median_contrast`: Effect sizes across hosts
- `p_value_wilcoxon_greater`, `p_value_wilcoxon_less`: Wilcoxon p-values (directional)
- `p_value_ttest_greater`, `p_value_ttest_less`: t-test p-values (directional)
- `q_value_wilcoxon_greater`, `q_value_wilcoxon_less`: FDR-adjusted Wilcoxon (primary)
- `q_value_ttest_greater`, `q_value_ttest_less`: FDR-adjusted t-test (secondary)
- `fisher_p_value_treatment`, `fisher_p_value_control` (if `use_fisher: true`): Fisher combined p-values
- `fisher_q_value_treatment`, `fisher_q_value_control` (if `use_fisher: true`): FDR-adjusted Fisher q-values

## Interpreting Results

### High-Significance Regions (q < 0.05)

- **Positive mean/median contrast** + low `q_value_wilcoxon_greater`: Treatment group evolves consistently *more* than control
  - Interpretation: Gene/window may be under stronger selective pressure, undergoing faster adaptation, or diverging more in treatment
- **Negative mean/median contrast** + low `q_value_wilcoxon_less`: Control group evolves more than treatment
  - Interpretation: Opposite pattern; control adaptation is stronger

### Effect Sizes

- **`median_contrast`**: Robust measure; use for interpretability ("on average, treatment regions are X percentile points higher")
- **`mean_contrast`**: Sensitive to outlier hosts; use cautiously but informative for overall magnitude
- **`n_hosts`**: Number of paired hosts contributing to the test. Fewer hosts → less statistical power

### Robustness Indicators

- **Agreement between Wilcoxon and t-test**: Both test types significant → robust finding
- **Wilcoxon only**: Suggests skewed or outlier-affected contrasts (non-parametric is more robust)
- **t-test only**: May indicate extreme host values; validate manually
- **Both insignificant**: No consistent signal; region likely under similar pressures in both groups

### Region Coverage

- **`informative_fraction` near 1.0**: Region well-sampled; results reliable
- **`informative_fraction` near 0**: Region sparsely sampled; treat result with caution
- Regions below `min_informative_fraction` threshold are excluded entirely

## Common Workflows

### Workflow 1: Default Gene-Level Analysis

```yaml
regional_contrast:
  mode: gene
  agg_method: median
  min_informative_sites: 5
  min_informative_fraction: 0.0
  use_fisher: false
```

Best for: Quick screening of gene-level signals without extra computation.

### Workflow 2: Gene + Window, Stringent Filtering

```yaml
regional_contrast:
  mode: both
  window_size: 1000
  agg_method: trimmed_mean
  trim_fraction: 0.15
  min_informative_sites: 10
  min_informative_fraction: 0.5
  use_fisher: false
```

Best for: High-confidence discoveries; requires dense coverage and multiple hosts.

### Workflow 3: Fine-Grained Windows with Exploratory Fisher

```yaml
regional_contrast:
  mode: window
  window_size: 500
  agg_method: median
  min_informative_sites: 3
  min_informative_fraction: 0.0
  use_fisher: true
```

Best for: Exploratory analysis at fine resolution; includes Fisher p-values for secondary validation.

## Troubleshooting

### Issue: No regions found after filtering

**Causes:** `min_informative_sites` or `min_informative_fraction` too strict.

**Solutions:**
- Lower thresholds (e.g., `min_informative_sites: 1`, `min_informative_fraction: 0.0`)
- Check input data quality
- Verify that your timepoint pair has multiple hosts with coverage

### Issue: Many NaN p-values

**Causes:** Fewer than 3 hosts, or insufficient site-level variation within regions.

**Solutions:**
- Increase number of biological replicates if possible
- Lower filtering thresholds
- Check that treatment and control groups are properly labeled

### Issue: All q-values = 1.0

**Causes:** Raw p-values uniformly distributed; no significant signal detected.

**Solutions:**
- Verify treatment/control assignment is correct
- Check input allele-frequency data quality
- Consider whether the groups truly differ evolutionarily
- Lower filtering thresholds to allow more regions into testing

### Issue: Output files missing for some group combinations

**Causes:** Configuration mismatch; group values not matching metadata.

**Solutions:**
- Verify `groups_combinations` values match exactly those in metadata
- Check that metadata has samples in both treatment and control groups
- Review logs for eligibility warnings

## References

- **Total-Variation Distance**: Classical metric in optimal transport; widely used in population genetics
- **Wilcoxon Signed-Rank Test**: Non-parametric paired test; robust to non-normal distributions
- **Benjamini–Hochberg FDR**: Controls expected proportion of false discoveries; more powerful than Bonferroni
- **Fisher's Method**: Combines independent p-values via $-2 \sum \ln(p_i) \sim \chi^2(2k)$; used for exploratory multi-perspective testing

## See Also

- [Configuration Reference](../reference/configuration.md) — Detailed parameter documentation
- [Statistical Tests Reference](../reference/statistical_tests.md) — Overview of all statistical methods
- [Regional Contrast Guide](https://github.com/MoellerLabPU/AlleleFlux/blob/main/docs/regional_contrast_guide.md) — In-depth technical guide
- [Regional Contrast Code Walkthrough](https://github.com/MoellerLabPU/AlleleFlux/blob/main/docs/regional_contrast_code_explanation.md) — Function-by-function code explanation
