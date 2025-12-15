# Use Cases

Real-world applications of AlleleFlux for metagenomic evolution studies.

## Antibiotic Resistance Evolution

**Question**: How do bacterial communities evolve under antibiotic treatment?

**Design**: Longitudinal fecal samples from antibiotic-treated vs. control mice (pre, during, post-treatment)

**AlleleFlux Analysis**:
\- High parallelism scores in treatment group → deterministic resistance evolution
\- Divergence scores identify differentially selected genes between groups
\- Outlier genes highlight novel resistance mechanisms

**Expected Findings**: Antibiotic resistance genes with high scores; correlation with treatment outcomes

## Diet-Microbiome Adaptation

**Question**: How do microbiomes adapt to dietary changes?

**Design**: High-fat vs. standard diet groups sampled at baseline, week 1, week 4

**AlleleFlux Analysis**:
\- Identify MAGs with diet-specific selection signatures
\- Focus on metabolic gene evolution
\- Track temporal adaptation patterns

**Expected Findings**: Nutrient acquisition genes under selection; taxonomic differences in adaptation rate

## Host-Microbe Co-evolution

**Question**: How do host genotypes shape microbial evolution?

**Design**: Multiple timepoints from different host genotypes (e.g., WT vs. knockout)

**AlleleFlux Analysis**:
\- Compare allele changes between host types
\- Identify host-specific selective pressures
\- Focus on host-interaction genes

**Expected Findings**: Host genotype-specific microbial adaptations; immune evasion genes

## Environmental Adaptation

**Question**: How do communities adapt to pollution?

**Design**: Contaminated vs. pristine sites sampled over time

**AlleleFlux Analysis**:
\- Parallel evolution of degradation pathways
\- Contaminant-specific gene selection
\- Biomarker discovery

**Expected Findings**: Degradation genes with high scores; potential bioremediation candidates

## Analysis Tips

**Data Quality**
\- Aim for ≥10x average coverage per MAG
\- ≥4 biological replicates per group
\- Match sequencing depth across samples

**Parameter Tuning**
\- Start with default QC thresholds (`breadth_threshold: 0.1`)
\- Use LMM for complex experimental designs
\- Enable CMH for detecting parallel changes

**Result Interpretation**
\- Compare multiple statistical approaches for robustness
\- Focus on genes with consistent signals across tests
\- Validate high-scoring genes with functional data

## Configuration Examples

See [Input Preparation](../usage/input_preparation.md) and [Configuration Reference](../reference/configuration.md) for detailed configuration guides.

**Longitudinal study** (antibiotic treatment):

```yaml
data_type: longitudinal
timepoints_combinations:
  - timepoint: [pre, during, post]
    focus: post
groups_combinations:
  - [antibiotic, control]
quality_control:
  breadth_threshold: 0.2
  min_sample_num: 6
```

**Single timepoint** (disease vs. healthy):

```yaml
data_type: single
timepoints_combinations:
  - timepoint: [baseline]
groups_combinations:
  - [disease, healthy]
quality_control:
  breadth_threshold: 0.1
  min_sample_num: 4
```

For complete worked examples, see [Tutorial](tutorial.md) and [Interpreting Results](../usage/interpreting_results.md).
