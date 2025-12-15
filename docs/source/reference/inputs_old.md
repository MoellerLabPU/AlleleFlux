# Input Files Reference

This page provides detailed specifications for all input files required by AlleleFlux.

\---

## BAM Files

Sorted and indexed BAM files from metagenomic samples aligned to reference MAGs.

**Requirements:**

- Format: `.bam` with accompanying `.bam.bai` index file
- Must be sorted by coordinate
- Contains alignments of reads to reference MAG contigs

**Naming:**

BAM files can have any naming convention. The sample identifier is obtained from:

1. The `bam_path` column in the metadata file
2. Or automatically extracted from the BAM filename

**Quality Considerations:**

- Higher mapping quality reads produce more reliable allele frequencies
- Duplicate reads should be marked (AlleleFlux does not remove duplicates)

\---

## Reference FASTA

A combined FASTA file containing all MAG contigs.

**Format Requirements:**

- File extension: `.fa` or `.fasta`
- Header format: `<MAG_ID>.fa_<contig_ID>`
- Must contain all contigs referenced in the BAM files

**Example Headers:**

```text
>Bacteroides_001.fa_k141_1234
ATCGATCGATCGATCG...
>Bacteroides_001.fa_k141_5678
GCTAGCTAGCTAGCTA...
>Lachnospira_002.fa_k141_9012
TTAAGCCTTAGGCCTT...
```

**Creating Combined FASTA:**

Use the `alleleflux-create-mag-mapping` utility:

```bash
alleleflux-create-mag-mapping \
    --dir /path/to/individual_mag_fastas \
    --extension fa \
    --output-fasta combined_mags.fasta \
    --output-mapping mag_mapping.tsv
```

\---

## Prodigal Genes FASTA

Gene predictions from Prodigal in nucleotide FASTA format.

**Format Requirements:**

- File extension: `.fna`
- Generated using Prodigal with nucleotide output (`-d` option)
- Gene IDs must match contig IDs in the reference FASTA

**Prodigal Header Format:**

```text
>Bacteroides_001.fa_k141_1234_1 # 100 # 450 # 1 # ID=1_1;partial=00;...
>Bacteroides_001.fa_k141_1234_2 # 460 # 900 # -1 # ID=1_2;partial=00;...
```

The header contains:
\- Gene ID (e.g., `Bacteroides_001.fa_k141_1234_1`)
\- Start position (1-based)
\- End position (1-based)
\- Strand (1 or -1)
\- Additional metadata

**Generating Prodigal Predictions:**

```bash
# For a single MAG
prodigal -i mag.fasta -d genes.fna -a genes.faa -o genes.gff -f gff

# For combined FASTA (recommended)
prodigal -i combined_mags.fasta -d prodigal_genes.fna -a prodigal_genes.faa -p meta
```

\---

## Metadata File

A tab-separated file containing sample information.

**Required Columns:**

```{eval-rst}
.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Description
   * - ``sample_id``
     - Unique identifier for each sample. Must be unique across all samples.
   * - ``bam_path``
     - **Full absolute path** to the BAM file for this sample.
   * - ``subjectID``
     - Identifier for the biological replicate or subject (e.g., mouse ID, patient ID).
   * - ``group``
     - Experimental group (e.g., "treatment", "control", "high_fat", "standard").
   * - ``replicate``
     - Replicate identifier within each group (e.g., "A", "B", "C").
```

**Additional Columns for Longitudinal Data:**

```{eval-rst}
.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Description
   * - ``time``
     - Timepoint identifier (e.g., "pre", "post", "day1", "week4").
```

**Example Metadata File:**

```text
sample_id   bam_path        subjectID       group   replicate       time
S1  /data/bams/sample1.bam  mouse1  control A       pre
S2  /data/bams/sample2.bam  mouse2  control B       pre
S3  /data/bams/sample3.bam  mouse3  treatment       A       pre
S4  /data/bams/sample4.bam  mouse4  treatment       B       pre
S5  /data/bams/sample5.bam  mouse1  control A       post
S6  /data/bams/sample6.bam  mouse2  control B       post
S7  /data/bams/sample7.bam  mouse3  treatment       A       post
S8  /data/bams/sample8.bam  mouse4  treatment       B       post
```

**Important Notes:**

- For longitudinal analysis, each `subjectID` must have entries at multiple timepoints
- The `replicate` column is used for stratification in CMH tests
- Group names should match those specified in the config file

**Adding BAM Paths:**

Use the utility script to add BAM paths to existing metadata:

```bash
alleleflux-add-bam-path \
    --metadata metadata_without_bam.tsv \
    --output metadata_with_bam.tsv \
    --bam-dir /path/to/bam_files \
    --bam-extension .bam
```

\---

## GTDB Taxonomy File

Taxonomy classification from GTDB-Tk.

**File:** `gtdbtk.bac120.summary.tsv` (bacterial) or `gtdbtk.ar53.summary.tsv` (archaeal)

**Required Columns:**

```{eval-rst}
.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Description
   * - ``user_genome``
     - MAG identifier (must match MAG IDs in reference FASTA)
   * - ``classification``
     - Full GTDB taxonomy string (e.g., ``d__Bacteria;p__Firmicutes;c__Clostridia;...``)
```

**Generating GTDB Classification:**

```bash
gtdbtk classify_wf --genome_dir /path/to/mags --out_dir gtdbtk_output --cpus 16
```

\---

## MAG Mapping File

A tab-separated file mapping contig names to MAG identifiers.

**Required Columns:**

```{eval-rst}
.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Description
   * - ``contig_name``
     - Full contig name as it appears in the reference FASTA
   * - ``mag_id``
     - MAG identifier
```

**Example:**

```text
contig_name mag_id
Bacteroides_001.fa_k141_1234        Bacteroides_001
Bacteroides_001.fa_k141_5678        Bacteroides_001
Lachnospira_002.fa_k141_9012        Lachnospira_002
```

**Creating Mapping File:**

```bash
alleleflux-create-mag-mapping \
    --dir /path/to/mag_fastas \
    --extension fa \
    --output-fasta combined_mags.fasta \
    --output-mapping mag_mapping.tsv
```

\---

## Input File Checklist

Before running AlleleFlux, verify:

```text
☐ BAM files are sorted and indexed (.bai files exist)
☐ Reference FASTA has correct header format (<MAG_ID>.fa_<contig_ID>)
☐ Prodigal genes match reference FASTA contigs
☐ Metadata file has all required columns
☐ BAM paths in metadata are absolute paths
☐ Each subjectID has samples at all required timepoints (longitudinal)
☐ At least min_sample_num samples per group
☐ GTDB taxonomy file contains all MAG IDs
☐ MAG mapping file covers all contigs
```

\---

## Profile Files (Step 2 Input)

Profile files are generated by the profiling step and serve as input for subsequent analysis. They contain base-level allele information extracted from BAM files.

**Location:** `profiles/{sample}/{sample}_{mag}_profiled.tsv.gz`

**Format:** Gzip-compressed TSV

**Columns:**

```{eval-rst}
.. list-table::
   :widths: 20 15 65
   :header-rows: 1

   * - Column
     - Type
     - Description
   * - ``contig``
     - string
     - Contig identifier (matches reference FASTA)
   * - ``position``
     - int
     - 0-based genomic position
   * - ``gene_id``
     - string
     - Overlapping gene identifier (empty if intergenic)
   * - ``ref_base``
     - string
     - Reference base at this position (A, C, G, T)
   * - ``A``
     - int
     - Count of adenine bases from aligned reads
   * - ``T``
     - int
     - Count of thymine bases from aligned reads
   * - ``G``
     - int
     - Count of guanine bases from aligned reads
   * - ``C``
     - int
     - Count of cytosine bases from aligned reads
```

**Example:**

```text
contig              position    gene_id                     ref_base    A    T    G    C
MAG_001.fa_contig1  120         MAG_001.fa_contig1_gene1    A           45   2    1    2
MAG_001.fa_contig1  121         MAG_001.fa_contig1_gene1    T           1    46   1    2
MAG_001.fa_contig1  122         MAG_001.fa_contig1_gene1    G           0    1    48   1
```

**Notes:**

- Total coverage at a position = A + T + G + C
- Allele frequency for base X = X / (A + T + G + C)
- Positions with zero coverage may be omitted

\---

## MAG-Specific Metadata Files

During the workflow, per-MAG metadata files are generated to track sample information for each MAG.

**Location:** `inputMetadata/inputMetadata_{timepoints}-{groups}/{mag}_metadata.tsv`

**Columns:**

```{eval-rst}
.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Description
   * - ``sample_id``
     - Sample identifier
   * - ``file_path``
     - Path to the sample's profile file for this MAG
   * - ``group``
     - Experimental group
   * - ``time``
     - Timepoint identifier (longitudinal only)
   * - ``subjectID``
     - Subject/replicate identifier
   * - ``replicate``
     - Replicate letter
```

\---

## Visualization Workflow Inputs

The visualization workflow has specific input requirements. See {doc}`../usage/visualization_guide` for the complete workflow.

### Significant Sites File

Input for terminal nucleotide analysis, typically from the p-value summary.

**Required Columns:**

```{eval-rst}
.. list-table::
   :widths: 20 15 65
   :header-rows: 1

   * - Column
     - Type
     - Description
   * - ``mag_id``
     - string
     - MAG identifier
   * - ``contig``
     - string
     - Contig identifier
   * - ``position``
     - int
     - 0-based genomic position
   * - ``gene_id``
     - string
     - Gene identifier
   * - ``test_type``
     - string
     - Statistical test used (e.g., ``two_sample_paired_tTest``)
   * - ``min_p_value``
     - float
     - Minimum p-value across tests
   * - ``q_value``
     - float
     - FDR-adjusted p-value
```

**Example:**

```text
mag_id      contig                  position    gene_id                     test_type               min_p_value    q_value
MAG_001     MAG_001.fa_contig1      120         MAG_001.fa_contig1_gene1    two_sample_paired_tTest 1.5e-06        2.3e-05
MAG_001     MAG_001.fa_contig1      145         MAG_001.fa_contig1_gene1    two_sample_paired_tTest 3.2e-05        1.8e-04
```

### Visualization Metadata File

Enhanced metadata for the visualization workflow with profile directory paths.

**Required Columns:**

```{eval-rst}
.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Description
   * - ``sample_id``
     - Unique sample identifier
   * - ``sample_profile_dir``
     - Path to the sample's profile directory (containing ``*_profiled.tsv.gz`` files)
   * - ``group``
     - Experimental group
   * - ``time``
     - Timepoint identifier
   * - ``subjectID``
     - Subject/replicate identifier for pairing
```

**Optional Columns:**

```{eval-rst}
.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Description
   * - ``day``
     - Numeric day for continuous time axis plotting
   * - ``replicate``
     - Replicate identifier within group
```

**Example:**

```text
sample_id           sample_profile_dir              group       time    subjectID
control_subj1_pre   profiles/control_subj1_pre      control     pre     subj1
control_subj1_post  profiles/control_subj1_post     control     post    subj1
treatment_subj2_pre profiles/treatment_subj2_pre    treatment   pre     subj2
```

\---

## Example Data

AlleleFlux provides example data for testing and learning the workflow:

**Location:** `docs/source/examples/example_data/`

**Contents:**

- `reference/combined_mags.fasta` - Reference FASTA with 2 test MAGs
- `reference/prodigal_genes.fna` - Gene predictions
- `reference/mag_mapping.tsv` - Contig-to-MAG mapping
- `reference/gtdbtk_taxonomy.tsv` - Mock GTDB taxonomy
- `metadata/sample_metadata.tsv` - Sample metadata (8 samples)
- `profiles/` - Pre-generated profile files
- `significant_sites/significant_sites.tsv` - Example significant sites
- `config_example.yml` - Working configuration file

**Usage:**

```bash
cd docs/source/examples/example_data
alleleflux run --config config_example.yml
```

For generating larger synthetic datasets:

```bash
python docs/source/examples/generate_synthetic_data.py \
    --output_dir my_test_data \
    --num_mags 10 \
    --num_samples 20
```
