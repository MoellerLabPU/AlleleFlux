Overview
========

What is AlleleFlux?
------------------

AlleleFlux is a bioinformatic tool designed to analyze site-specific allele frequency fluctuations in microbial genomes across multiple timepoints and experimental groups. It addresses the challenge of quantifying and comparing allele-frequency shifts in microbiota by providing robust statistical methods specifically developed for longitudinal shotgun (meta)genomic datasets.

Key Features
-----------

1. **Flexible Data Analysis**
   
   * Supports both single-timepoint and longitudinal studies
   * Handles different experimental groups and timepoints
   * Works with metagenome-assembled genomes (MAGs)

2. **Statistical Testing Framework**
   
   * **Parallelism Score**: Measures consistent allele frequency changes over time across replicates in each group
   * **Divergence Score**: Quantifies the degree of allele-frequency divergence between experimental groups
   * Multiple statistical testing approaches:
     - Two-sample tests (paired and unpaired)
     - Single-sample analysis
     - Linear Mixed Models (LMM)
     - Cochran-Mantel-Haenszel (CMH) test

3. **Outlier Detection**
   
   * Identifies genes exhibiting exceptionally high parallelism or divergence scores
   * Pinpoints genomic regions under strong selection

4. **Taxonomic Analysis**
   
   * Aggregates scores across different taxonomic levels (phylum, class, order, family, genus, species)
   * Enables comparison of evolutionary dynamics across different taxa

5. **Quality Control**
   
   * Filtering options for allele frequency data
   * Eligibility tests for MAGs to ensure meaningful statistical analysis

Basic Workflow
-------------

AlleleFlux's workflow is implemented as a Snakemake pipeline and consists of several key steps:

.. figure:: /_static/images/alleleflux_workflow.svg
   :width: 700px
   :align: center
   :alt: AlleleFlux workflow diagram
   
   **Figure 1:** Overview of the AlleleFlux workflow showing the key steps from input data to results interpretation.

.. note::
   If the image above is not displaying correctly, here is a text-based representation of the workflow:
   
   .. code-block:: none
      :caption: AlleleFlux Workflow
   
      Input Data        →        Preprocessing        →        Allele Analysis
      ▔▔▔▔▔▔▔▔▔▔▔               ▔▔▔▔▔▔▔▔▔▔▔▔▔                ▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔
      MAG sequences              Quality control               Compute frequencies
      Sample metadata            Eligibility tables            Analyze changes
      BAM files                  Metadata preparation          Filter (optional)
                                                                      ↓
      ▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁           ▁▁▁▁▁▁▁▁▁▁▁▁▁▁          ▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁
      Results Analysis  ←        Outlier Detection    ←      Statistical Testing
      ▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔           ▔▔▔▔▔▔▔▔▔▔▔▔▔▔          ▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔  
      Summary statistics         High score genes             Two-sample tests
      Visualizations             Statistical                  Single-sample
      Functional analysis        assessment                   LMM / CMH tests
                                      ↑↓
                               ▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁
                               Score Generation
                               ▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔
                               Parallelism scores
                               Divergence scores
                               Taxonomic aggregation

The workflow consists of the following steps:

1. **Input Preparation**
   
   * MAG sequences and annotations
   * Sample metadata with timepoints and experimental groups
   * Read mapping data (BAM files)

2. **Preprocessing**
   
   * Quality control of input data
   * Generation of eligibility tables for MAGs
   * Metadata preparation

3. **Allele Frequency Analysis**
   
   * Computation of allele frequencies at each genomic position
   * Analysis of allele frequency changes over time
   * Filtering constant alleles (optional)

4. **Statistical Testing**
   
   * Application of appropriate statistical tests based on study design:
     - Two-sample tests for comparing groups
     - Single-sample analysis for within-group changes
     - Linear Mixed Models for complex experimental designs
     - CMH test for position-by-position assessment of allele frequency changes

5. **Score Generation**
   
   * Calculation of parallelism and divergence scores for each MAG
   * Aggregation of scores across genes and taxonomic levels

6. **Outlier Detection**
   
   * Identification of genes with significantly high scores
   * Statistical assessment of outliers

7. **Results Interpretation**
   
   * Generation of summary statistics and visualizations
   * Integration with functional annotations for biological interpretation

Use Cases
--------

AlleleFlux is particularly useful for studies involving:

* Experimental evolution of microbial communities
* Host-associated microbiomes under different conditions (e.g., diet, medication)
* Environmental microbiomes responding to changing conditions
* Comparative analysis of microbial adaptation across different experimental treatments
* Identification of genes under selection in longitudinal studies

By providing a comprehensive framework for analyzing allele frequency dynamics, AlleleFlux enables researchers to gain deeper insights into the selective pressures and evolutionary mechanisms shaping microbial communities over time.
