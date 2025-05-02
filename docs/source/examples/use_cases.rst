Use Cases
=========

This page describes real-world use cases for AlleleFlux.

Case Study 1: Antibiotic Resistance Evolution
----------------------------------------------

**Research Question**: How do bacterial communities evolve in response to antibiotic treatment?

**Experimental Design**:
- Treatment group: Mice treated with antibiotics
- Control group: Mice without antibiotic treatment
- Timepoints: Pre-treatment, during treatment, post-treatment
- Samples: Fecal samples collected at each timepoint
- Sequencing: Shotgun metagenomic sequencing

**AlleleFlux Analysis**:
1. Profile allele frequencies in bacterial genomes across all samples
2. Identify MAGs with high parallelism scores in the treatment group
3. Compare parallelism and divergence scores between treatment and control groups
4. Identify genes with exceptionally high scores, focusing on known antibiotic resistance genes

**Key Findings**:
- Higher parallelism scores in antibiotic resistance genes in the treatment group
- Identification of previously unknown genes involved in antibiotic resistance
- Correlation between parallelism scores and clinical outcomes

Case Study 2: Adaptation to Dietary Changes
-------------------------------------------

**Research Question**: How do gut microbiomes adapt to changes in diet?

**Experimental Design**:
- Groups: High-fat diet vs. standard diet
- Timepoints: Before diet change, 1 week after, 4 weeks after
- Samples: Fecal samples
- Sequencing: Shotgun metagenomic sequencing

**AlleleFlux Analysis**:
1. Compare allele frequency changes between diet groups
2. Identify MAGs with high divergence scores between diet groups
3. Focus on genes involved in metabolism and nutrient acquisition
4. Correlate gene-level scores with diet-specific nutrients

**Key Findings**:
- Identification of specific bacterial genes undergoing selection in response to dietary changes
- Temporal patterns of adaptation in different taxonomic groups
- Metabolic pathways under strongest selection

Case Study 3: Host-Microbiome Co-evolution
-------------------------------------------

**Research Question**: How do host-specific factors influence microbial evolution?

**Experimental Design**:
- Groups: Different host genotypes (e.g., wild-type vs. immune-deficient)
- Timepoints: Multiple timepoints over development
- Samples: Gut microbiome samples
- Sequencing: Shotgun metagenomic sequencing

**AlleleFlux Analysis**:
1. Compare allele frequency changes between host genotypes
2. Identify host-specific patterns of microbial evolution
3. Focus on microbial genes involved in host-microbe interactions
4. Correlate microbial evolution with host phenotypes

**Key Findings**:
- Host genotype-specific selection on microbial genes
- Temporal dynamics of host-microbe co-evolution
- Identification of key genes mediating host-microbe interactions

Case Study 4: Environmental Adaptation
---------------------------------------

**Research Question**: How do environmental microbial communities adapt to pollution exposure?

**Experimental Design**:
- Groups: Contaminated vs. pristine environments
- Timepoints: Before contamination, multiple timepoints after
- Samples: Soil or water samples
- Sequencing: Shotgun metagenomic sequencing

**AlleleFlux Analysis**:
1. Profile allele frequencies in microbial genomes across samples
2. Compare parallelism scores between contaminated and pristine sites
3. Identify genes under strongest selection in response to contaminants
4. Correlate gene-level scores with contaminant concentrations

**Key Findings**:
- Identification of microbial genes involved in contaminant degradation
- Temporal patterns of community adaptation
- Potential biomarkers for environmental monitoring

Implementation Tips
----------------

For all use cases:

1. **Data Quality**:
   - Ensure sufficient sequencing depth for accurate allele frequency estimation
   - Carefully design experimental groups and timepoints
   - Collect detailed metadata for correlation analyses

2. **Analysis Parameters**:
   - Adjust breadth threshold based on sequencing depth
   - Consider different statistical tests for different experimental designs
   - Balance sensitivity and specificity in p-value thresholds

3. **Interpretation**:
   - Validate key findings with independent methods
   - Consider functional annotations when interpreting gene-level results
   - Integrate results with other omics data when available