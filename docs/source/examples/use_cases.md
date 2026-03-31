# Use Cases

Real-world applications of AlleleFlux for metagenomic evolution studies.

## Antibiotic Resistance Evolution

**Question**: How do bacterial communities evolve under antibiotic treatment?

**Design**: Longitudinal fecal samples from antibiotic-treated vs. control mice (pre, during, post-treatment). 10 treated, 8 control mice with 3 samples each = 54 total samples. To analyze with AlleleFlux, we compare two timepoints at a time using separate configurations: pre→during and during→post.

**Complete Configuration (pre→during comparison)**:

```yaml
input:
   bam_dir: data/bam/
   fasta: data/reference/combined_mags.fa
   metadata: data/metadata/sample_metadata.tsv
   prodigal_path: data/prodigal/combined_genes.gff
   gtdb_taxonomy: data/taxonomy/gtdb_taxonomy.tsv
   mag_mapping: data/mapping/mag_contig_mapping.tsv

 output:
   root_dir: output/antibiotic_resistance/pre_during/

 quality_control:
   breadth_threshold: 0.2           # Require 20% genome coverage
   min_sample_num: 6                # At least 6 samples per group-timepoint

 analysis:
   data_type: longitudinal
   timepoints_combinations:
     - timepoint: [pre, during]       # Compare pre → during treatment
       focus: during                  # Focus on during-treatment state
   groups_combinations:
     - [antibiotic, control]
   use_lmm: true
   use_significance_tests: true
   use_cmh: true

 statistics:
   p_value_threshold: 0.05
   fdr_method: fdr_bh
   fdr_threshold: 0.1

 resources:
   allele_freq:
     cpus: 8
     mem_mb: 8000
     time_min: 60
   statistical_tests:
     cpus: 16
     mem_mb: 16000
     time_min: 120
```

**Additional Configuration (during→post comparison)**:

For the post-treatment phase, create a second config file with:
```yaml
timepoints_combinations:
  - timepoint: [during, post]      # Compare during → post treatment
    focus: post                     # Focus on post-treatment state
```

**Expected command**:

```bash
# Phase 1: Early response to antibiotic
alleleflux run --config config_antibiotic_pre_during.yml --threads 16

# Phase 2: Late response and stabilization
alleleflux run --config config_antibiotic_during_post.yml --threads 16
```

**Output files to examine**:
- `pre_during/step1_output/eligibility_tables/eligibility_table_pre_during-antibiotic_control.tsv` - Which MAGs pass QC
- `pre_during/step2_output/scores/mag_level_scores.tsv` - Parallelism scores during initial treatment response
- `during_post/step2_output/scores/mag_level_scores.tsv` - Parallelism scores during stabilization phase
- Combined gene-level scores identifying resistance candidates across both phases
- Combined statistical results controlling for mouse identity

**Result interpretation**:
1. **Early response phase (pre→during)**: MAGs with parallelism scores >2-3% show rapid adaptation to antibiotic stress. Expect lower scores than late phase
2. **Late stabilization phase (during→post)**: Scores may increase further if selection continues, or plateau if community has stabilized
3. **Divergence between phases**: Compare scores to identify which genes are selected early (rapid response) vs. late (fine-tuning). Genes present in both phases indicate core resistance mechanisms
4. **Cross-group comparison**: Compare antibiotic group evolution to control group (expect <2% parallelism in controls for both phases)
5. **Gene-level outliers**: Look for genes in CARD database (antibiotic resistance genes). Unknown genes are candidates for novel resistance mechanisms
6. **CMH significance**: p-values < 0.05 indicate consistent allele changes across individual mice, ruling out within-mouse noise
7. **Example interpretation**: If *Bacteroides* MAG shows parallelism score 3% in pre→during (treated group) with early outliers including `tetR`, and then 6% in during→post with additional outliers in a porin gene, this indicates: (1) early selection for direct resistance, (2) late selection for drug efflux optimization

## Diet-Microbiome Adaptation

**Question**: How do microbiomes adapt to dietary changes?

**Design**: Two diet groups (high-fat vs. standard diet), 12 mice per group, sampled at baseline, week 1, and week 4 = 72 total samples

**Complete Configuration**:

```yaml
input:
  bam_dir: data/bam/
  fasta: data/reference/combined_mags.fa
  metadata: data/metadata/diet_metadata.tsv
  prodigal_path: data/prodigal/combined_genes.gff
  gtdb_taxonomy: data/taxonomy/gtdb_taxonomy.tsv
  mag_mapping: data/mapping/mag_contig_mapping.tsv

output:
  root_dir: output/diet_adaptation/baseline_week1/

quality_control:
  breadth_threshold: 0.1           # Lower threshold (may have diet-dependent coverage)
  min_sample_num: 4                # 12 mice per group, require 4 across samples

analysis:
  data_type: longitudinal
  timepoints_combinations:
    - timepoint: [baseline, week1]
      focus: week1                 # Focus on initial dietary response
  groups_combinations:
    - [highfat, standard]
  use_lmm: true
  use_significance_tests: true
  use_cmh: true

statistics:
  p_value_threshold: 0.05
  fdr_method: fdr_bh
  fdr_threshold: 0.1

resources:
  allele_freq:
    cpus: 12
    mem_mb: 12000
    time_min: 90
  statistical_tests:
    cpus: 12
    mem_mb: 12000
    time_min: 90
```

**Expected command**:

```bash
alleleflux run --config config_diet.yml --threads 12
```

**Additional Configuration (week1→week4 comparison)**:

For the longer-term adaptation trajectory, create a second config file with:
```yaml
output:
  root_dir: output/diet_adaptation/week1_week4/

timepoints_combinations:
  - timepoint: [week1, week4]      # Compare week1 → week4 adaptation
    focus: week4                     # Focus on steady-state adaptations
```

**Output files to examine**:
- Phase 1: `baseline_week1/step1_output/eligibility_tables/eligibility_table_baseline_week1-highfat_standard.tsv`
- Phase 2: `week1_week4/step1_output/eligibility_tables/eligibility_table_week1_week4-highfat_standard.tsv`
- Phase 1 scores: `baseline_week1/step2_output/scores/mag_level_scores.tsv` - Early dietary response
- Phase 2 scores: `week1_week4/step2_output/scores/mag_level_scores.tsv` - Stabilization phase
- `step2_output/scores/gene_scores.tsv` - Metabolic genes under selection
- `step2_output/scores/taxa_scores.tsv` - Phylum/Family level aggregation for dietary responders
- `step2_output/outliers/outlier_genes.tsv` - Genes under selection, mapped to KEGG pathways

**Result interpretation**:
1. **Two-phase adaptation**: Phase 1 (baseline→week1) shows rapid initial response; Phase 2 (week1→week4) shows slower stabilization
2. **Metabolic focus**: Filter `gene_scores.tsv` for genes annotated with: carbohydrate transport, lipid metabolism, short-chain fatty acid synthesis
3. **Taxa-level patterns**: In high-fat group, expect high scores for Bacteroidetes (fiber fermenters) and Faecalibacterium (SCFA producers). Standard diet shows stable low scores in both phases
4. **Phase comparison**: Scores may increase, plateau, or decrease from phase 1 to phase 2 depending on adaptation kinetics
5. **Example interpretation**: A Roseburia MAG may show parallelism 3% in phase 1 (baseline→week1) with early outliers in carbohydrate transport, then 5% in phase 2 (week1→week4) with additional outliers in propionate synthesis genes. This indicates initial colonization of gut with subsequent metabolic fine-tuning

## Host-Microbe Co-evolution

**Question**: How do host genotypes shape microbial evolution?

**Design**: Longitudinal samples from WT vs. knockout mice, 10 mice per group, 3 timepoints = 60 total samples

**Complete Configuration**:

```yaml
input:
  bam_dir: data/bam/
  fasta: data/reference/combined_mags.fa
  metadata: data/metadata/host_genotype_metadata.tsv
  prodigal_path: data/prodigal/combined_genes.gff
  gtdb_taxonomy: data/taxonomy/gtdb_taxonomy.tsv
  mag_mapping: data/mapping/mag_contig_mapping.tsv

output:
  root_dir: output/host_microbe/week0_week4/

quality_control:
  breadth_threshold: 0.15
  min_sample_num: 5                # 10 mice per group

analysis:
  data_type: longitudinal
  timepoints_combinations:
    - timepoint: [week0, week4]
      focus: week4                 # Compare early → mid-stage genotype effects
  groups_combinations:
    - [wildtype, knockout]
  use_lmm: true                    # Critical: captures mouse-level random effects
  use_significance_tests: true
  use_cmh: true

statistics:
  p_value_threshold: 0.05
  fdr_method: fdr_bh
  fdr_threshold: 0.1

resources:
  allele_freq:
    cpus: 8
    mem_mb: 8000
    time_min: 60
  statistical_tests:
    cpus: 16
    mem_mb: 16000
    time_min: 120
```

**Expected command**:

```bash
alleleflux run --config config_host.yml --threads 16
```

**Additional Configuration (week4→week8 comparison)**:

For the longer-term co-evolution analysis, create a second config file with:
```yaml
output:
  root_dir: output/host_microbe/week4_week8/

timepoints_combinations:
  - timepoint: [week4, week8]
    focus: week8                   # Compare mid-stage → late-stage genotype effects
```

**Output files to examine**:
- Phase 1: `week0_week4/step2_output/scores/mag_level_scores.tsv` - Early genotype-dependent selection
- Phase 2: `week4_week8/step2_output/scores/mag_level_scores.tsv` - Late genotype-dependent selection
- Phase 1 LMM: `week0_week4/step2_output/statistical_tests/lmm_results.tsv` - Test genotype effect early
- Phase 2 LMM: `week4_week8/step2_output/statistical_tests/lmm_results.tsv` - Test genotype effect late
- `step2_output/outliers/outlier_genes.tsv` - Genotype-dependent adaptive genes
- `step2_output/scores/gene_scores.tsv` - Filtered for surface proteins, secretion systems

**Result interpretation**:
1. **Genotype-specific MAGs**: MAGs with high scores in KO but low/zero in WT indicate host-genotype-dependent selection
2. **LMM significance**: p-values < 0.05 in LMM show genotype effect while controlling for individual mouse variation
3. **Gene annotation**: Focus outliers on: outer membrane proteins, secretion systems (T6SS, Sec), immune-related factors
4. **Convergent evolution**: If multiple MAGs show similar outliers (e.g., same flagellar genes), this suggests host-driven convergence
5. **Temporal dynamics**: Compare phase 1 (week0→week4) vs. phase 2 (week4→week8) to reveal whether selection is rapid early or gradual throughout
6. **Example interpretation**: If Bacteroides MAG shows parallelism 0% in WT but 4% in IL-10KO during phase 1, then 0% in WT but 7% in IL-10KO during phase 2, with outliers in flagellar genes and mucus-degrading glycosidases in both phases, this indicates sustained and strengthening IL-10-dependent selection for motile mucus-colonizers

## Environmental Adaptation

**Question**: How do communities adapt to pollution?

**Design**: Contaminated vs. pristine soil sites, sampled at month 0, 3, 6. Multiple sites per treatment = 30+ samples

**Complete Configuration**:

```yaml
input:
  bam_dir: data/bam/
  fasta: data/reference/combined_mags.fa
  metadata: data/metadata/environmental_metadata.tsv
  prodigal_path: data/prodigal/combined_genes.gff
  gtdb_taxonomy: data/taxonomy/gtdb_taxonomy.tsv
  mag_mapping: data/mapping/mag_contig_mapping.tsv

output:
  root_dir: output/environmental/month0_month3/

quality_control:
  breadth_threshold: 0.1           # Environmental samples often have lower coverage
  min_sample_num: 3                # 4-5 replicates per group-timepoint

analysis:
  data_type: longitudinal
  timepoints_combinations:
    - timepoint: [month0, month3]
      focus: month3
  groups_combinations:
    - [contaminated, pristine]
  use_lmm: false                   # Site effects complex; pair-wise comparisons instead
  use_significance_tests: true
  use_cmh: true

statistics:
  p_value_threshold: 0.05
  fdr_method: fdr_bh
  fdr_threshold: 0.1

resources:
  allele_freq:
    cpus: 8
    mem_mb: 8000
    time_min: 60
  statistical_tests:
    cpus: 12
    mem_mb: 12000
    time_min: 90
```

**Expected command**:

```bash
alleleflux run --config config_env.yml --threads 8
```

**Additional Configuration (month3→month6 comparison)**:

For the longer-term environmental adaptation analysis, create a second config file with:
```yaml
output:
  root_dir: output/environmental/month3_month6/

timepoints_combinations:
  - timepoint: [month3, month6]
    focus: month6                  # Track sustained adaptation to contamination
```

**Output files to examine**:
- Phase 1: `month0_month3/step2_output/scores/mag_level_scores.tsv` - Initial contamination response
- Phase 2: `month3_month6/step2_output/scores/mag_level_scores.tsv` - Sustained adaptation or stabilization
- Phase 1 CMH: `month0_month3/step2_output/statistical_tests/cmh_results.tsv` - Parallel response across sites
- Phase 2 CMH: `month3_month6/step2_output/statistical_tests/cmh_results.tsv` - Sustained parallel effects
- `step2_output/outliers/outlier_genes.tsv` - Pollutant-specific genes (e.g., heavy metal resistance, hydrocarbon degradation)
- `step2_output/scores/gene_scores.tsv` - Pathway annotation for functional interpretation

**Result interpretation**:
1. **Biphasic adaptation**: Phase 1 (month0→month3) shows rapid initial response; Phase 2 (month3→month6) tracks sustained or enhanced adaptation
2. **CMH significance**: p < 0.05 confirms parallel evolution across contaminated replicates in each phase, not chance
3. **Temporal progression**: Contaminated site scores may increase across both phases; pristine site scores remain low and flat in both
4. **Gene-level biomarkers**: Outliers encoding heavy metal efflux (CzcA, CopA), hydrocarbon degradation (alkane hydroxylase, cytochrome P450), or xenobiotic pathways are credible pollutant-responsive genes
5. **Functional categories**: Use KEGG pathway annotation to identify complete degradation operons under selection
6. **Example interpretation**: Arthrobacter MAG shows parallelism 5% at contaminated site during phase 1 (month0→month3), with early outliers for mercury resistance (merA). During phase 2 (month3→month6), parallelism increases to 9.5%, with additional outliers for arsenic efflux (arsB) and PCB degradation. This indicates multi-stage adaptation: initial mercury response followed by broader multi-contaminant tolerance. Genes absent from pristine site samples in both phases confirm pollution-driven selection

## Fecal Microbiota Transplant (FMT) Study

**Question**: How does the donor microbiome adapt and stabilize after transplant into recipients?

**Design**: 15 FMT recipients sampled pre-FMT, day 1, week 1, month 1, month 3 post-FMT. Include 5 donor samples for baseline = 80 total samples. Track whether donor-derived taxa establish and whether they evolve to match recipient genetics.

**Complete Configuration**:

```yaml
input:
  bam_dir: data/bam/
  fasta: data/reference/combined_mags.fa
  metadata: data/metadata/fmt_metadata.tsv
  prodigal_path: data/prodigal/combined_genes.gff
  gtdb_taxonomy: data/taxonomy/gtdb_taxonomy.tsv
  mag_mapping: data/mapping/mag_contig_mapping.tsv

output:
  root_dir: output/fmt_adaptation/pre_fmt_month3_post/

quality_control:
  breadth_threshold: 0.15          # Clinical samples vary; balance coverage vs. MAG count
  min_sample_num: 5                # ~15 recipients give ample replicates

analysis:
  data_type: longitudinal
  timepoints_combinations:
    - timepoint: [pre_fmt, month3_post]
      focus: month3_post           # Final steady state
  groups_combinations:
    - [recipient, donor]           # Compare recipient+donor evolution
  use_lmm: true                    # Critical: each recipient is unique host environment
  use_significance_tests: true
  use_cmh: true

statistics:
  p_value_threshold: 0.05
  fdr_method: fdr_bh
  fdr_threshold: 0.1

resources:
  allele_freq:
    cpus: 12
    mem_mb: 16000
    time_min: 120
  statistical_tests:
    cpus: 16
    mem_mb: 20000
    time_min: 180
```

**Expected command**:

```bash
alleleflux run --config config_fmt.yml --threads 16
```

**Fine-grained temporal tracking**:

For detailed temporal resolution of early establishment and late adaptation phases, create additional config files:
```yaml
# Early phase: initial colonization shock
output:
  root_dir: output/fmt_adaptation/pre_fmt_day1_post/
timepoints_combinations:
  - timepoint: [pre_fmt, day1_post]
    focus: day1_post

# Mid phase: establishment
output:
  root_dir: output/fmt_adaptation/day1_post_week1_post/
timepoints_combinations:
  - timepoint: [day1_post, week1_post]
    focus: week1_post

# Late phase: stabilization
output:
  root_dir: output/fmt_adaptation/week1_post_month3_post/
timepoints_combinations:
  - timepoint: [week1_post, month3_post]
    focus: month3_post
```

**Output files to examine**:
- Long-term: `pre_fmt_month3_post/step2_output/scores/mag_level_scores.tsv` - Overall recipient vs donor trajectory
- Early phase: `pre_fmt_day1_post/step2_output/scores/mag_level_scores.tsv` - Colonization shock signature
- Establishment: `day1_post_week1_post/step2_output/scores/mag_level_scores.tsv` - Early stabilization
- Late phase: `week1_post_month3_post/step2_output/scores/mag_level_scores.tsv` - Final adaptation
- All phases: `step2_output/statistical_tests/lmm_results.tsv` - Recipient-level variation
- All phases: `step2_output/outliers/outlier_genes.tsv` - Recipient-specific adaptations (phase-dependent)

**Result interpretation**:
1. **Bifurcation in scores**: Long-term analysis (pre_fmt→month3_post) should show donor group scores stable ~0-2% across all timepoints; recipient group scores increase from ~1% (pre-FMT) to 3-8% by month3_post
2. **Temporal kinetics across phases**:
   - **Early phase (pre→day1)**: Colonization shock; large allele frequency swings as stressed donor cells meet new environment; weak but detectable parallelism (~1-2%)
   - **Establishment phase (day1→week1)**: Intermediate dynamics; scores increase toward stabilization plateau
   - **Late phase (week1→month3)**: Plateau phase; scores stabilize as donor microbiota equilibrates; CMH p-values strengthen (p < 0.01) indicating coordinated adaptation across recipients
3. **Gene-level biomarkers vary by phase**:
   - **Early (day1_post)**: High allele frequency variance; genes involved in stress response, adhesion initiation
   - **Late (month3_post)**: Stabilized genes in colonization (adhesins, mucin-binding), nutrient scavenging (vitamin synthesis, carbohydrate transport), immune evasion (flagellar reduction, LPS modification)
4. **Cross-recipient parallelism**: High CMH significance in late phase indicates similar donor lineages undergo similar selective pressures in different recipients, revealing universal recipient-environment constraints
5. **Example interpretation**: Donor-derived Faecalibacterium MAG shows parallelism 0.5% (pre-FMT, long-term view) or complex signal in early phases (pre→day1 shows 1%, day1→week1 shows 2%, week1→month3 shows 5% incremental). Outlier genes at month3 include butyrate production pathway genes and flagellar genes (selected out). Compare to donor samples (score ~0% across all phases), confirming sustained recipient-driven selection

## Experimental Evolution in Bioreactors

**Question**: How do microbial communities evolve during serial passage in controlled bioreactor conditions?

**Design**: Replicate bioreactors (e.g., n=4) under two conditions: high temperature (40°C) vs. standard (37°C). Samples taken every 50 generations for 200 generations = 5 timepoints × 2 conditions × 4 replicates = 40 samples.

**Complete Configuration**:

```yaml
input:
  bam_dir: data/bam/
  fasta: data/reference/combined_mags.fa
  metadata: data/metadata/bioreactor_metadata.tsv
  prodigal_path: data/prodigal/combined_genes.gff
  gtdb_taxonomy: data/taxonomy/gtdb_taxonomy.tsv
  mag_mapping: data/mapping/mag_contig_mapping.tsv

output:
  root_dir: output/bioreactor_evolution/gen0_gen200/

quality_control:
  breadth_threshold: 0.2           # Experimental design = uniform high coverage
  min_sample_num: 6                # 4 replicates per condition per timepoint

analysis:
  data_type: longitudinal
  timepoints_combinations:
    - timepoint: [gen0, gen200]
      focus: gen200                # Compare ancestral vs. final evolved state
  groups_combinations:
    - [high_temp, standard]
  use_lmm: true                    # Captures bioreactor-specific effects
  use_significance_tests: true
  use_cmh: true

statistics:
  p_value_threshold: 0.05
  fdr_method: fdr_bh
  fdr_threshold: 0.1

resources:
  allele_freq:
    cpus: 8
    mem_mb: 8000
    time_min: 60
  statistical_tests:
    cpus: 16
    mem_mb: 16000
    time_min: 120
```

**Expected command**:

```bash
alleleflux run --config config_bioreactor.yml --threads 16
```

**Intermediate timepoint analysis**:

For resolution of evolutionary kinetics across generational transitions, create two additional configs:
```yaml
# Early-to-mid evolution
output:
  root_dir: output/bioreactor_evolution/gen0_gen100/
timepoints_combinations:
  - timepoint: [gen0, gen100]
    focus: gen100

# Mid-to-late evolution
output:
  root_dir: output/bioreactor_evolution/gen100_gen200/
timepoints_combinations:
  - timepoint: [gen100, gen200]
    focus: gen200
```

**Output files to examine**:
- Long-term: `gen0_gen200/step2_output/scores/mag_level_scores.tsv` - Overall 200-generation trajectory per condition
- Early-mid: `gen0_gen100/step2_output/scores/mag_level_scores.tsv` - Initial adaptation phase
- Mid-late: `gen100_gen200/step2_output/scores/mag_level_scores.tsv` - Secondary adaptation or plateau
- All phases: `step2_output/statistical_tests/cmh_results.tsv` - Parallel evolution across bioreactor replicates
- All phases: `step2_output/outliers/outlier_genes.tsv` - Adaptive genes, timestamped by phase
- All phases: `step2_output/evolution/dnds_results.tsv` - dN/dS ratio for adaptive validation (if enabled)

**Result interpretation**:
1. **Strong parallelism signature** (long-term, gen0→gen200): Expect parallelism scores 10-20% in high_temp group (cleaner signal than in vivo due to controlled conditions). Standard group should stay <2%
2. **CMH statistical strength**: High-replicate experimental design (n=4 bioreactors) yields very strong CMH signals (p < 0.001), confirming reproducible evolution, not sampling noise
3. **Temporal kinetics** (multi-phase view):
   - **Phase 1 (gen0→gen100)**: Weak to moderate signal (parallelism ~1-5%), stochastic drift→early favorable mutations sweep
   - **Phase 2 (gen100→gen200)**: Scores plateau or increase further (parallelism 10-15% additional), late-stage stabilization or secondary mutations
   - **Long-term (gen0→gen200)**: Combined signal showing cumulative evolution (10-20% total)
4. **Gene appearance timeline**: Parse `outlier_genes.tsv` by phase:
   - **Early outliers (gen0-gen100)**: First-hit beneficial genes (e.g., heat shock protein upregulators)
   - **Late outliers (gen100-gen200)**: Second-site compensatory mutations or fine-tuning
5. **Predicted temperature adaptation genes**: Heat shock proteins (GroEL, DnaK), membrane lipid remodeling, oxidative stress resistance (catalase, superoxide dismutase)
6. **dN/dS validation**: If outlier genes show dN/dS > 1, confirms positive selection. dN/dS near 0-0.2 may indicate relaxed purifying selection (hitchhiking) or structural evolution
7. **Example interpretation**: A core metabolic MAG (e.g., Achromobacter) shows: Phase 1 (gen0→gen100) parallelism 2% with early outlier flagellar gene (likely thermotaxis); Phase 2 (gen100→gen200) parallelism 13% additional with outliers in GroEL and transporter gene. CMH p-value 1e-6 confirms all 4 high-temp bioreactors independently selected this gene set in both phases. Compare to standard condition (parallelism <1% in both phases), confirming temperature-driven selection. This pattern suggests: (1) initial thermotaxis response (phase 1), (2) proteostasis and nutrient transport optimization at elevated temperature (phase 2)

## Configuration Strategy Guide

Choosing the right configuration for your study design is critical. This section summarizes the key decisions and tradeoffs.

### Data Type: `single` vs. `longitudinal`

| Aspect | `data_type: "single"` | `data_type: "longitudinal"` |
|--------|----------------------|---------------------------|
| **Use when** | One timepoint only; comparing cross-sectional groups | Multiple timepoints; tracking evolution over time |
| **Sample design** | Disease vs. healthy; treatment A vs. treatment B | Pre/during/post; day0 → day7 → day30 |
| **Statistical tests** | Unpaired two-sample test, single-sample tests | Paired two-sample test, CMH, LMM |
| **Power requirements** | More samples needed (n=10-15 per group) | Fewer samples (n=4-8 per group) - paired design has higher power |
| **Output structure** | Flat: one score per MAG per group | Hierarchical: one score per MAG per timepoint per group |
| **Example** | Gut microbiota in IBD vs. control | Antibiotic resistance during treatment: pre/during/post |

**Decision rule**: If you have multiple timepoints, use `longitudinal` for better statistical power and ability to detect temporal dynamics.

### Statistical Tests: When to Enable Each

| Test | Enable when | Key output | Notes |
|------|-------------|-----------|-------|
| **Unpaired two-sample** (`use_significance_tests: true`) | Any design with 2+ groups | `p_value` per MAG | Always compute; foundational test |
| **LMM** (`use_lmm: true`) | Unbalanced designs, repeated measures, covariates | `lmm_p_value`, `effect_size` | Use for mouse/host/individual variation; essential for clinical studies |
| **CMH** (`use_cmh: true`) | High replicates (n≥4), stratified designs | `cmh_p_value`, `stratified_odds_ratio` | Detects consistent allele changes across replicates; excellent for experimental replicates |

**Decision rules**:
- **Longitudinal + LMM**: Clinical/in vivo studies (e.g., FMT, antibiotic treatment). LMM handles repeated measures per subject
- **Longitudinal + CMH**: Experimental designs with multiple independent replicates (e.g., bioreactors, replicate mice)
- **All three**: High-power designs (n≥10 samples per group, n≥4 temporal points) - maximize signal detection
- **Single + unpaired only**: Low-budget studies, cross-sectional designs

### Quality Control Thresholds: `breadth_threshold` and `min_sample_num`

| Scenario | `breadth_threshold` | `min_sample_num` | Rationale |
|----------|------------------|-----------------|-----------|
| **High-quality isolated culture** | 0.5-1.0 | 3 | Uniform coverage; can be strict |
| **Human gut microbiota (high biomass)** | 0.2-0.3 | 5-6 | Deep sequencing + abundant organisms |
| **Environmental samples** | 0.05-0.1 | 3-4 | Sparse coverage; retain rare MAGs |
| **Ultra-low biomass (lung, blood)** | 0.01-0.05 | 2-3 | Contamination risk; threshold critical |
| **Mixed/stressed microbiota** | 0.15 | 4 | Intermediate; account for uneven sampling |

**Interpretation of thresholds**:
- **`breadth_threshold`**: Fraction of genomic positions with ≥1 coverage. A threshold of 0.2 means "MAG must be present at 20%+ of its genome in a sample"
- **`min_sample_num`**: Minimum sample count required for a MAG to pass QC per group-timepoint combination

**Tuning strategy**:
1. Start conservative (breadth 0.2, min_sample 6) to get clean signal
2. If too few MAGs pass QC, relax to 0.15/5 or 0.1/4
3. Never go below 0.05 breadth (high false positive risk from mapping errors)
4. Never set min_sample_num below the number of replicates in your smallest group

### Resource Allocation Recommendations

```yaml
resources:
  allele_freq:          # Profiling phase (step 1)
    cpus: 8-16          # Scales linearly with # MAGs and samples
    mem_mb: 8000-16000  # Mostly BAM loading
    time_min: 60-120    # Depends on sequencing depth

  statistical_tests:    # Analysis phase (step 2)
    cpus: 12-32         # Parallelizes well across genes/MAGs
    mem_mb: 16000-32000 # Keeps large matrices in memory
    time_min: 120-240   # Scales with # positions, # tests, # replicates
```

**Scaling rules**:
- **CPUs**: Use `--threads` equal to available cores. AlleleFlux parallelizes MAG and gene processing
- **Memory**: 
  - Minimum 8 GB (single MAG, small dataset)
  - Standard 16 GB (typical microbiome study, 50-200 MAGs)
  - High-demand 32+ GB (>500 MAGs, 1000+ samples, many-replicate CMH tests)
- **Time**:
  - Step 1 profiling: ~5-10 min per sample per core
  - Step 2 analysis: ~1-10 min per MAG (depends on test complexity, # positions)
  - CMH tests dominate runtime on high-replicate designs

**Cluster submission (SLURM)**:
```bash
# Typical microbiome study
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=240

alleleflux run --config config.yml --threads 16
```

### Matching Config to Study Characteristics

**Large animal study** (e.g., cattle microbiome, n=30 cattle, 3 timepoints):
```yaml
data_type: longitudinal
breadth_threshold: 0.2
min_sample_num: 8
use_lmm: true         # Accounts for individual animal variation
use_cmh: true         # 30 animals = many replicates
```

**Clinical trial** (e.g., probiotic intervention, n=20 subjects, pre/post):
```yaml
data_type: longitudinal
breadth_threshold: 0.15
min_sample_num: 5
use_lmm: true         # Essential: subject-level random effects
use_significance_tests: true
use_cmh: false        # Only 20 subjects; LMM better handles imbalance
```

**Bioreactor experiment** (n=4 replicates, 5 timepoints, controlled):
```yaml
data_type: longitudinal
breadth_threshold: 0.2
min_sample_num: 6     # Exceed replicate count for robustness
use_lmm: true
use_cmh: true         # 4 replicates perfect for CMH
```

**Environmental survey** (e.g., soil sites, low coverage, n=10 sites):
```yaml
data_type: single     # One-time sampling
breadth_threshold: 0.1
min_sample_num: 3     # Conservative with coverage
use_lmm: false
use_significance_tests: true
```

---

For complete worked examples, see [Tutorial](tutorial.md) and [Interpreting Results](../usage/interpreting_results.md).
