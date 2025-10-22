
# dN/dS Analysis from Timepoint Comparisons (NG86 Path Averaging)

## 1. Overview

The [`dnds_from_timepoints.py`](../../../alleleflux/scripts/evolution/dnds_from_timepoints.py) script performs **dN/dS ratio analysis** to detect evolutionary selection patterns by comparing ancestral and derived states at statistically significant sites identified by AlleleFlux.

**Important Update:** As of the latest version, this script implements the **Nei-Gojobori (1986) path averaging method** for codons where multiple positions change. This provides more accurate dN/dS estimates than the previous per-site approach.

### Purpose

The **dN/dS ratio** (also written as ω) is the ratio of the rate of **non-synonymous substitutions** (dN) to the rate of **synonymous substitutions** (dS) in protein-coding sequences. This ratio is a fundamental measure in molecular evolution that helps distinguish between different modes of natural selection:

- **dN/dS < 1**: **Purifying (negative) selection** - Most mutations are deleterious and are removed by natural selection. This indicates functional constraint on the protein.
- **dN/dS ≈ 1**: **Neutral evolution** - Mutations are neither beneficial nor harmful, accumulating at random.
- **dN/dS > 1**: **Positive (adaptive) selection** - Non-synonymous mutations confer a fitness advantage and are favored by natural selection.

### Role in AlleleFlux Workflow

This script integrates into the AlleleFlux workflow by:
1. Accepting significant sites identified by [`p_value_summary.py`](../../../alleleflux/scripts/preprocessing/p_value_summary.py)
2. Reconstructing ancestral sequences from allele frequency profiles
3. Identifying substitutions between timepoints
4. Calculating dN/dS ratios at multiple levels (codon, gene, MAG, global)

### Nei-Gojobori Method with Path Averaging

The script implements the **Nei-Gojobori (1986)** method for calculating dN/dS ratios with proper **path averaging** for multi-position changes. This method:

1. **Calculates potential synonymous (S) and non-synonymous (N) sites** for each codon by considering all possible single-nucleotide changes
2. **Groups sites by codon** to identify codons where multiple positions have changed (k=1, 2, or 3)
3. **For single-position changes (k=1)**: Directly classifies as synonymous or non-synonymous
4. **For multi-position changes (k>1)**: Enumerates all k! mutational pathways and averages S/NS classifications across paths, accounting for path-dependent amino acid changes
5. **Computes fractional counts** for each codon and sums across codons:
   - pN = (sum of fractional NS counts) / (potential non-synonymous sites)
   - pS = (sum of fractional S counts) / (potential synonymous sites)
   - dN/dS = pN / pS

**Why Path Averaging Matters:** When multiple positions change in a codon, the classification of each change (S vs NS) depends on which mutations occurred first. The NG86 method accounts for this by considering all possible orders and averaging, providing more accurate estimates than treating each position independently.

The Nei-Gojobori method is computationally efficient and appropriate for analyzing within-species polymorphism data from metagenomic samples.

---

## 2. Detailed Methodology

### 2.1 Nei-Gojobori Method Implementation

#### Potential Sites Calculation

For each codon, the method determines the **potential** number of synonymous and non-synonymous sites by examining all possible single-nucleotide changes.

**Algorithm:**

For a given codon (e.g., "ATG"):

1. Iterate through each of the 3 positions in the codon
2. For each position, consider all 3 possible alternative nucleotides (excluding the original)
3. Create the mutant codon and translate it to an amino acid
4. Compare the mutant amino acid to the original:
   - If they match → count as **synonymous** (S)
   - If they differ → count as **non-synonymous** (NS)
5. Since there are 9 total changes (3 positions × 3 alternatives), each position contributes 1/3 to the final count

**Fractional Site Counts:**

The final potential S and N sites for a codon are fractional values calculated as:
- **S = (count of synonymous changes) / 3**
- **N = (count of non-synonymous changes) / 3**

This averaging approach accounts for the fact that each nucleotide position can have different proportions of synonymous vs. non-synonymous changes.

#### Path Averaging for Multi-Nucleotide Changes (k > 1)

**Key Innovation:** When a codon has multiple positions changed (k=2 or k=3), the script enumerates all k! mutational pathways and averages the S/NS classifications.

**Why This Matters:**

Consider codon "AAA" (Lysine) changing to "TTG" (Leucine) - 3 positions changed:
- Path 1: AAA → TAA (stop) → TTA (Leu) → TTG (Leu)
- Path 2: AAA → ATA (Ile) → TTA (Leu) → TTG (Leu)
- Path 3: AAA → AAT (Asn) → TAT (Tyr) → TTG (Leu)
- ... (3! = 6 total pathways)

Each pathway may classify the mutations differently. The NG86 method:
1. Enumerates all pathways
2. Classifies each mutation in each pathway
3. **Excludes pathways through intermediate stop codons**
4. Averages the S/NS counts across valid pathways
5. Reports fractional counts (e.g., 1.5 NS, 1.5 S if half the paths yield different classifications)

**Stop Codon Handling:**
- Pathways that create intermediate stop codons are excluded (biologically unrealistic)
- Final stop codons are allowed (represent valid termination mutations)

**Formulas:**

For a gene with multiple codons:
```
Total_potential_S = Σ(S_codon_i)
Total_potential_N = Σ(N_codon_i)
```

For observed substitutions (with fractional counts from path averaging):
```
pN = Σ(fractional_NS_codon_i) / Total_potential_N
pS = Σ(fractional_S_codon_i) / Total_potential_S
dN/dS = pN / pS
```

#### Codon-by-Codon Analysis with NG86 Cache

The implementation uses **two pre-computed caches** to maximize efficiency:

1. **Potential sites cache** (`_CODON_SITE_CACHE`): All 64 codons with their S/N potential sites
2. **NG86 path averaging cache** (`_NG86_CODON_PAIR_CACHE`): All 4,096 codon pairs (64 × 64) with pre-computed fractional S/NS counts

**Cache building at startup:**
- The `_precompute_codon_sites_cache()` function calculates S and N values for all 64 codons
- The `_build_ng86_codon_pair_cache()` function pre-computes path-averaged S/NS counts for all codon pairs
- During analysis, lookups are O(1) for instant retrieval

**Benefits:**
- Avoids redundant pathway enumerations during analysis
- Significant performance improvement for large datasets
- Guarantees consistent results across all analyses

### 2.2 Data Processing Pipeline

#### Step 1: Loading and Filtering Significant Sites

The `setup_and_load_data()` function:

1. Reads the output TSV file from [`p_value_summary.py`](../../../alleleflux/scripts/preprocessing/p_value_summary.py)
2. Filters for the specified MAG ID
3. Optionally filters by test type (e.g., "two_sample_paired_tTest")
4. Optionally filters by group analyzed
5. Applies significance threshold using either:
   - `min_p_value`: Raw p-values from statistical tests
   - `q_value`: FDR-corrected q-values (recommended)
6. Removes sites without valid gene IDs or with multiple overlapping genes

**Required columns in input file:**
- `mag_id`: Metagenome-assembled genome identifier
- `contig`: Contig name
- `position`: 0-based position on the contig
- `gene_id`: Gene identifier from Prodigal
- Selected p-value column (`min_p_value` or `q_value`)

#### Step 2: Loading Profile Data

Profile data from [`profile_mags.py`](../../../alleleflux/scripts/analysis/profile_mags.py) provides allele frequency information:

**Profile format (TSV):**
```
contig  position  ref_base  gene_id  A  T  G  C  coverage
contig1 1000      G         gene_1   5  2  150 3  160
contig1 1001      A         gene_1   148 4  2   6  160
```

Each row represents a variable site with counts for each nucleotide.

#### Step 3: Ancestral Sequence Reconstruction

The `reconstruct_ancestral_sequences()` function:

1. Starts with the **reference sequence** from the Prodigal FASTA file
2. For each variable site in the ancestral profile:
   - Determines the **major allele** using the `get_major_allele()` function
   - Handles ties by preferring the reference base or choosing randomly
3. Substitutes major alleles into the reference sequence
4. Handles strand orientation properly (see Section 2.3)

**Major Allele Determination:**
- Counts for each nucleotide (A, T, G, C) are compared
- The allele with the highest count is selected
- In case of ties, the reference base is preferred
- If the reference isn't part of the tie, a random selection breaks the tie

#### Step 4: Identifying Substitutions with NG86 Path Averaging

The `find_substitutions()` function (now a wrapper) calls `analyze_codon_substitutions_with_ng86_paths()`:

1. **Groups sites by codon:**
   - For each significant site, determine which codon it belongs to
   - Group all significant sites within the same codon together
   - Identify how many positions changed (k=1, 2, or 3)

2. **For each codon with substitutions:**
   - Retrieve ancestral and derived sequences
   - Extract the codon before and after
   - Look up pre-computed fractional S/NS counts from NG86 cache
   - If not in cache (rare), compute via path enumeration

3. **Path Enumeration (when needed):**
   - The `_enumerate_ng86_paths()` function generates all k! mutation pathways
   - Each pathway is checked for intermediate stop codons (excluded if found)
   - Each valid pathway's mutations are classified as S or NS
   - Fractional counts are averaged across all valid pathways

4. **Output:**
   - Returns codon-level events with fractional S/NS counts
   - Each row represents one codon, even if multiple positions changed
   - Preserves all positional information for tracking

**Example with k=2 (two positions changed):**

Codon "ATG" (Met) → "CTT" (Leu):
- Positions 0 and 2 changed
- 2! = 2 pathways:
  1. ATG → CTG → CTT: NS (Met→Leu), S (Leu→Leu) = 1 NS, 1 S
  2. ATG → ATT → CTT: NS (Met→Ile), NS (Ile→Leu) = 2 NS, 0 S
- Average: **1.5 NS, 0.5 S** (fractional counts)

#### Step 5: Global dN/dS Calculation with Fractional Counts

The `calculate_global_dnds_for_sites()` function implements **codon-level summation** of fractional counts:

**Algorithm:**
1. Load codon-level substitution results (from NG86 analysis)
2. For each codon event:
   - Extract potential S/N sites for the codon
   - Add `frac_S` and `frac_N` (fractional observed counts) to totals
   - Track unique codons using `(gene_id, codon_start_index)` tuples
3. Sum fractional counts across all codons:
   - `Total_observed_S = Σ(frac_S_codon_i)`
   - `Total_observed_NS = Σ(frac_N_codon_i)`
4. Calculate rates and ratio:
   - `pS = Total_observed_S / Total_potential_S`
   - `pN = Total_observed_NS / Total_potential_N`
   - `dN/dS = pN / pS`

**Key Difference from Previous Version:**
- Now sums **fractional** S/NS counts (not integer counts)
- Properly handles multi-position changes via path averaging
- Still maintains anti-double-counting for potential sites

**Example:**

| Codon | k | frac_S | frac_N | Potential S | Potential N |
|-------|---|--------|--------|-------------|-------------|
| ATG→CTT | 2 | 0.5 | 1.5 | 0.0 | 3.0 |
| GTA→GTG | 1 | 1.0 | 0.0 | 0.67 | 2.33 |
| **Totals** | | **1.5** | **1.5** | **0.67** | **5.33** |

- pS = 1.5 / 0.67 = 2.24
- pN = 1.5 / 5.33 = 0.28
- dN/dS = 0.28 / 2.24 = 0.125

### 2.3 Coordinate Systems and Strand Handling

#### Coordinate Systems

The script handles two coordinate systems:

1. **Contig coordinates (0-based):** Positions from the allele frequency profiles
2. **Gene coordinates (0-based):** Positions within individual gene sequences

**Prodigal coordinates** are 1-based in the FASTA headers but are converted to 0-based for internal calculations.

#### Forward Strand Genes

For forward strand genes (strand = +1):

```
Contig:  0    1    2    3    4    5    6    7    8    9   ...
Gene:              0    1    2    3    4    5
         ←gene start (1-based = 3)
```

Position in gene = contig_position - (gene_start - 1)

**Example:** If a SNV is at contig position 5, and the gene starts at position 3 (1-based):
- Position in gene = 5 - (3 - 1) = 5 - 2 = 3

#### Reverse Strand Genes

For reverse strand genes (strand = -1):

The gene sequence in the FASTA file is the **reverse complement** of the genomic sequence. The first base of this sequence corresponds to the gene's **end** position on the contig.

```
Contig:  0    1    2    3    4    5    6    7    8    9   ...
Gene:                            5    4    3    2    1    0
                                              ←gene end (1-based = 8)
```

Position in gene = (gene_end - 1) - contig_position

**Example:** If a SNV is at contig position 5, and the gene ends at position 8 (1-based):
- Position in gene = (8 - 1) - 5 = 7 - 5 = 2

#### Strand-Specific Allele Handling

**Major Alleles vs. Effective Alleles:**

- **Major alleles** from profiles are always on the **forward strand**
- **Effective alleles** for gene sequences must match the strand orientation

For **reverse strand genes**, alleles must be **complemented**:
- A ↔ T
- G ↔ C

The `get_codon_from_site()` function handles this coordinate conversion automatically.

---

## 3. Worked Examples

### Example 1: Calculating Potential Sites for Codon "ATG" (Methionine/Start)

**Original codon:** ATG → Amino acid: **M** (Methionine)

**Step-by-step analysis of all 9 possible changes:**

**Position 0 (A):**
- A→**T**: TTG → **L** (Leucine) - **NS**
- A→**G**: GTG → **V** (Valine) - **NS**
- A→**C**: CTG → **L** (Leucine) - **NS**

**Position 1 (T):**
- T→**A**: AAG → **K** (Lysine) - **NS**
- T→**G**: AGG → **R** (Arginine) - **NS**
- T→**C**: ACG → **T** (Threonine) - **NS**

**Position 2 (G):**
- G→**A**: ATA → **I** (Isoleucine) - **NS**
- G→**T**: ATT → **I** (Isoleucine) - **NS**
- G→**C**: ATC → **I** (Isoleucine) - **NS**

**Summary:**
- Synonymous changes: **0** out of 9
- Non-synonymous changes: **9** out of 9

**Final calculation:**
- S = 0 / 3 = **0.0** potential synonymous sites
- N = 9 / 3 = **3.0** potential non-synonymous sites

**Interpretation:** ATG is the start codon and codes uniquely for methionine. Every single-nucleotide change results in a different amino acid, making it highly constrained.

### Example 2: Calculating Potential Sites for Codon "TTT" (Phenylalanine)

**Original codon:** TTT → Amino acid: **F** (Phenylalanine)

**Step-by-step analysis:**

**Position 0 (T):**
- T→**A**: ATT → **I** (Isoleucine) - **NS**
- T→**G**: GTT → **V** (Valine) - **NS**
- T→**C**: CTT → **L** (Leucine) - **NS**

**Position 1 (T):**
- T→**A**: TAT → **Y** (Tyrosine) - **NS**
- T→**G**: TGT → **C** (Cysteine) - **NS**
- T→**C**: TCT → **S** (Serine) - **NS**

**Position 2 (T):**
- T→**A**: TTA → **L** (Leucine) - **NS**
- T→**G**: TTG → **L** (Leucine) - **NS**
- T→**C**: TTC → **F** (Phenylalanine) - **S** ✓

**Summary:**
- Synonymous changes: **1** out of 9 (position 2: T→C)
- Non-synonymous changes: **8** out of 9

**Final calculation:**
- S = 1 / 3 = **0.33** potential synonymous sites
- N = 8 / 3 = **2.67** potential non-synonymous sites

**Interpretation:** TTT has partial synonymy at the third position (wobble position), where one change (T→C) maintains phenylalanine.

### Example 3: Complete Substitution Analysis with Path Averaging (k=2)

Let's work through a complete workflow with a multi-position change.

**Scenario:**
- Gene: `gene_12345` on forward strand
- Two adjacent positions both have significant changes
- Position 1020: Ancestral G (150 reads) → Derived A (145 reads)
- Position 1021: Ancestral T (148 reads) → Derived C (150 reads)
- These positions fall in the same codon

**Step 1: Extract ancestral codon context**

From the ancestral sequence at positions 1020-1022:
- Position in gene for 1020 = 1020 - (gene_start - 1) = 1020 - 999 = 21
- Codon start = (21 // 3) × 3 = 21 (positions 21, 22, 23)
- Ancestral codon = "**GT**A" (positions 0,1,2 of codon map to contig 1020,1021,1022)

**Step 2: Create derived codon**

Replace positions 0 and 1 with new alleles A and C:
- Derived codon = "**AC**A"

**Step 3: Identify k value**

- Positions 0 and 1 changed
- k = 2, so we need path averaging (2! = 2 pathways)

**Step 4: Enumerate pathways**

Using the `_enumerate_ng86_paths()` function:

**Pathway 1:** Change position 0 first, then position 1
- GTA (Val) → **A**TA (Ile) - **NS** (Val ≠ Ile)
- ATA (Ile) → A**C**A (Thr) - **NS** (Ile ≠ Thr)
- Result: **2 NS, 0 S**

**Pathway 2:** Change position 1 first, then position 0
- GTA (Val) → G**C**A (Ala) - **NS** (Val ≠ Ala)
- GCA (Ala) → **A**CA (Thr) - **NS** (Ala ≠ Thr)
- Result: **2 NS, 0 S**

**Step 5: Average across pathways**

- Average NS = (2 + 2) / 2 = **2.0 NS**
- Average S = (0 + 0) / 2 = **0.0 S**

**Step 6: Record the codon-level event**

```
mag_id: MAG_001
contig: contig_1
position: 1020,1021  (both positions listed)
gene_id: gene_12345
codon_before: GTA
codon_after: ACA
aa_before: V
aa_after: T
k: 2
frac_S: 0.0
frac_N: 2.0
potential_S_sites_codon: 0.67
potential_N_sites_codon: 2.33
```

**Interpretation:** Even though both pathways yielded the same classification (all NS), the NG86 approach correctly accounts for the path-dependent nature of multi-position changes.

### Example 4: Global dN/dS Calculation with Fractional Counts

Let's calculate a global dN/dS ratio with NG86 path-averaged data.

**Data:**

| Codon Event | Gene | Codon Start | k | Codon Before→After | frac_S | frac_N | S_codon | N_codon |
|-------------|------|-------------|---|-------------------|--------|--------|---------|---------|
| Event 1 | gene_1 | 21 | 2 | GTA→ACA | 0.0 | 2.0 | 0.67 | 2.33 |
| Event 2 | gene_1 | 24 | 1 | CCT→CCA | 1.0 | 0.0 | 0.33 | 2.67 |
| Event 3 | gene_2 | 15 | 3 | ATG→TTC | 0.83 | 2.17 | 0.00 | 3.00 |

**Step 1: Identify unique codons** (anti-double-counting for potential sites)

All three events are in different codons, so all contribute to potential sites.

**Step 2: Sum potential sites for all codons**

- Total potential S = 0.67 + 0.33 + 0.00 = **1.00**
- Total potential N = 2.33 + 2.67 + 3.00 = **8.00**

**Step 3: Sum fractional observed counts**

- Total observed S = 0.0 + 1.0 + 0.83 = **1.83**
- Total observed NS = 2.0 + 0.0 + 2.17 = **4.17**

**Step 4: Calculate rates**

- pS = 1.83 / 1.00 = **1.830**
- pN = 4.17 / 8.00 = **0.521**

**Step 5: Calculate dN/dS**

- dN/dS = 0.521 / 1.830 = **0.285**

**Interpretation:** The dN/dS ratio of 0.285 indicates **purifying selection**. Note how fractional counts from path averaging (Event 3 with k=3) are properly incorporated into the global calculation.

---

## 4. Functional Breakdown

### 4.1 Core Calculation Functions

#### `_enumerate_ng86_paths(codon1, codon2, table)`

**Purpose:** Enumerate all mutational pathways between two codons and compute path-averaged S/NS counts.

**Input:**
- `codon1` (str): Starting codon (3 bases)
- `codon2` (str): Ending codon (3 bases)
- `table` (CodonTable): Biopython genetic code table

**Output:**
- Tuple of (frac_S, frac_N, k_value) where:
  - `frac_S`: Average synonymous count across valid pathways
  - `frac_N`: Average non-synonymous count across valid pathways
  - `k_value`: Number of positions that differ

**Algorithm:**
1. Identify which positions differ (k = 1, 2, or 3)
2. Generate all k! permutations of mutation orders
3. For each pathway:
   - Apply mutations in order
   - Check for intermediate stop codons (exclude pathway if found)
   - Classify each mutation as S or NS
   - Track cumulative S/NS counts for the pathway
4. Average counts across all valid (non-stop) pathways
5. Return fractional results

**Example:**
```python
from Bio.Data import CodonTable
table = CodonTable.unambiguous_dna_by_id[11]
frac_s, frac_n, k = _enumerate_ng86_paths("ATG", "CTT", table)
# Returns: (0.5, 1.5, 2) for a k=2 change with path averaging
```

#### `_build_ng86_codon_pair_cache(table_id=11)`

**Purpose:** Pre-compute fractional S/NS counts for all 4,096 codon pairs (64×64).

**Input:**
- `table_id` (int): NCBI genetic code table ID (default: 11)

**Output:**
- Dictionary mapping `(codon1, codon2)` tuples to `(frac_S, frac_N, k)` tuples

**Key points:**
- Called once at script startup
- Stores results in global `_NG86_CODON_PAIR_CACHE`
- Enables O(1) lookup during analysis
- Includes all possible codon transitions including multi-position changes

#### `_get_ng86_cache()`

**Purpose:** Lazy initialization wrapper for the NG86 cache.

**Output:**
- Reference to the global `_NG86_CODON_PAIR_CACHE` (builds if not yet initialized)

**Usage:** Called internally by analysis functions to ensure cache exists before lookups.

#### `_calculate_codon_sites(codon, table)`

**Purpose:** Calculate potential S and N sites for a single codon.

**Input:**
- `codon` (str): 3-base DNA string (e.g., "ATG")
- `table` (CodonTable): Biopython genetic code table (default: table 11 for bacteria)

**Output:**
- Tuple of (S_sites, N_sites) as floats

**Algorithm:**
1. Determine the amino acid for the original codon
2. For each of the 3 positions:
   - Test all 3 alternative bases
   - Translate the mutant codon
   - Compare amino acids to classify as S or NS
3. Return (total_S / 3, total_N / 3)

**Example usage:**
```python
from Bio.Data import CodonTable
table = CodonTable.unambiguous_dna_by_id[11]
s_sites, n_sites = _calculate_codon_sites("ATG", table)
# Returns: (0.0, 3.0)
```

#### `_precompute_codon_sites_cache(table_id=11)`

**Purpose:** Pre-compute S/N potential sites for all 64 codons to maximize performance.

**Input:**
- `table_id` (int): NCBI genetic code table ID (default: 11)

**Output:**
- Dictionary mapping codon strings to (S, N) tuples

**Key points:**
- Called once at script startup
- Stores results in the global `_CODON_SITE_CACHE` variable
- Eliminates redundant calculations during analysis

#### `calculate_potential_sites_for_gene(gene_seq_str)`

**Purpose:** Sum potential S and N sites across an entire gene sequence.

**Input:**
- `gene_seq_str` (str): DNA sequence of the gene

**Output:**
- Dictionary with keys "S" and "N" containing totals

**Algorithm:**
1. Iterate through sequence in 3-base steps
2. Look up pre-computed values from cache for each codon
3. Skip partial codons (sequence length not divisible by 3)
4. Sum values across all complete codons

**Example:**
```python
gene_seq = "ATGAAACGTTAA"  # 12 bases = 4 codons
potential_sites = calculate_potential_sites_for_gene(gene_seq)
# Returns: {"S": sum of S values, "N": sum of N values}
```

#### `_group_sites_into_codon_events(sites_df, ancestral_seqs, derived_profile_df, prodigal_records, ancestral_major_alleles)`

**Purpose:** Transform per-site data into per-codon structure for NG86 analysis.

**Input:**
- `sites_df` (DataFrame): Significant sites
- `ancestral_seqs` (dict): Reconstructed ancestral sequences
- `derived_profile_df` (DataFrame): Derived timepoint profiles
- `prodigal_records` (dict): Gene metadata
- `ancestral_major_alleles` (dict): Pre-computed ancestral alleles

**Output:**
- DataFrame with one row per codon, columns include:
  - Codon identification: `gene_id`, `codon_start_index`
  - Position lists: `positions_in_codon` (which of 0,1,2 changed)
  - Sequences: `codon_before`, `codon_after`
  - All positional information for tracking

**Algorithm:**
1. For each significant site, determine its codon
2. Group sites by `(gene_id, codon_start_index)`
3. For each group:
   - Extract ancestral and derived codons
   - Record which positions changed
   - Preserve all relevant metadata
4. Return codon-level DataFrame

#### `analyze_codon_substitutions_with_ng86_paths(codon_events_df, prodigal_records)`

**Purpose:** Apply NG86 path averaging to codon-level events and compute fractional S/NS counts.

**Input:**
- `codon_events_df` (DataFrame): Codon events from grouping function
- `prodigal_records` (dict): Gene metadata

**Output:**
- DataFrame with codon-level results including:
  - `k`: Number of positions changed
  - `frac_S`, `frac_N`: Fractional S/NS counts from path averaging
  - `potential_S_sites_codon`, `potential_N_sites_codon`
  - Translations: `aa_before`, `aa_after`

**Algorithm:**
1. For each codon event:
   - Look up fractional counts from NG86 cache
   - If not cached, compute via `_enumerate_ng86_paths()`
   - Get potential sites from codon sites cache
   - Translate codons to amino acids
2. Compile comprehensive results per codon
3. Return enriched DataFrame

**Key Feature:** Handles k=1, 2, and 3 uniformly through cache lookup or dynamic computation.

#### `calculate_global_dnds_for_sites(substitution_results_df, ancestral_seqs, prodigal_records)`

**Purpose:** Calculate a single global dN/dS ratio using fractional counts with proper codon deduplication.

**Input:**
- `substitution_results_df` (DataFrame): Codon-level substitution results with fractional counts
- `ancestral_seqs` (dict): Reconstructed ancestral sequences
- `prodigal_records` (dict): Gene metadata

**Output:**
- Tuple of (summary_df, codon_specific_df):
  - `summary_df`: Global dN/dS metrics including k-value distribution
  - `codon_specific_df`: Per-codon potential S/N values

**Key Algorithm:**
1. Use a set to track unique `(gene_id, codon_start_index)` tuples
2. For each codon event:
   - Extract potential S/N sites
   - Add to codon-specific table
   - Only add to global potential totals if codon not previously seen
3. Sum fractional `frac_S` and `frac_N` across all codons (not deduplicated - each codon counted once)
4. Calculate pN, pS, and dN/dS
5. Include statistics on k-value distribution

**Handles edge cases:**
- Multiple codon events (already deduplicated at codon level)
- Invalid codons (logs warning and skips)
- Division by zero (returns NaN)
- Fractional counts (proper summation)

### 4.2 Data Loading and Reconstruction Functions

#### `setup_and_load_data(significant_sites_path, mag_id, p_value_column, p_value_threshold, test_type, group_analyzed)`

**Purpose:** Load and filter significant sites data.

**Input:**
- `significant_sites_path` (Path): Path to p_value_summary output TSV
- `mag_id` (str): MAG identifier to filter for
- `p_value_column` (str): "min_p_value" or "q_value"
- `p_value_threshold` (float): Significance threshold
- `test_type` (str, optional): Statistical test to filter for
- `group_analyzed` (str, optional): Group name to filter for

**Output:**
- DataFrame of filtered significant sites, or None if no sites pass filters

**Filtering steps:**
1. Validate required columns exist
2. Filter for specified MAG ID
3. Optionally filter by test type
4. Optionally filter by group analyzed
5. Apply significance threshold
6. Remove sites without valid single gene IDs

#### `get_major_allele(row, ref_base)`

**Purpose:** Determine the most frequent allele at a position.

**Input:**
- `row` (Series): DataFrame row with A, T, G, C count columns
- `ref_base` (str): Reference allele for tie-breaking

**Output:**
- Single character string representing the major allele

**Tie-breaking logic:**
1. Find maximum count across A, T, G, C
2. If only one allele has max count, return it
3. If multiple alleles tie:
   - Prefer reference base if it's in the tie
   - Otherwise, choose randomly to avoid systematic bias

**Example:**
```python
row = pd.Series({"A": 50, "T": 50, "G": 5, "C": 5})
major = get_major_allele(row, "A")  # Returns "A" (ref base breaks tie)

row = pd.Series({"A": 45, "T": 50, "G": 5, "C": 5})
major = get_major_allele(row, "G")  # Returns "T" (clear winner)
```

#### `reconstruct_ancestral_sequences(unique_genes, prodigal_records, profile_by_gene)`

**Purpose:** Reconstruct ancestral sequences by substituting major alleles into reference sequences.

**Input:**
- `unique_genes` (list): Gene IDs to reconstruct
- `prodigal_records` (dict): Gene metadata and sequences
- `profile_by_gene` (DataFrameGroupBy): Ancestral profile data grouped by gene

**Output:**
- Tuple of:
  - `ancestral_orfs` (dict): Maps gene_id to reconstructed sequence string
  - `ancestral_major_alleles` (dict): Maps (contig, position) to forward-strand major allele

**Algorithm:**
1. Start with reference sequence from Prodigal
2. For each variable site in the gene's profile:
   - Get major allele from profile
   - Convert coordinates to gene-relative position
   - Handle strand orientation (complement for reverse strand)
   - Substitute into sequence
3. Store both the sequence and the forward-strand allele

**Coordinate handling:**
- Uses the `get_codon_from_site()` function for coordinate conversion
- Properly handles forward and reverse strand genes

### 4.3 Substitution Analysis Functions

#### `find_substitutions(sites_df, ancestral_seqs, prodigal_records, derived_profile_df, ancestral_major_alleles)`

**Purpose:** Wrapper function that calls the NG86 path averaging analysis pipeline.

**Input:**
- `sites_df` (DataFrame): Significant sites to check
- `ancestral_seqs` (dict): Reconstructed ancestral sequences
- `prodigal_records` (dict): Gene metadata
- `derived_profile_df` (DataFrame): Derived timepoint profile data
- `ancestral_major_alleles` (dict): Pre-computed ancestral major alleles

**Output:**
- DataFrame with codon-level substitution events including fractional S/NS counts

**Algorithm:**
1. Call `_group_sites_into_codon_events()` to transform per-site data to per-codon
2. Call `analyze_codon_substitutions_with_ng86_paths()` to apply NG86 method
3. Return enriched DataFrame with fractional counts and k-values

**Note:** This is now a wrapper that orchestrates the NG86 pipeline. The previous per-site analysis has been replaced with codon-centric NG86 path averaging.

#### `analyze_mutation_effect(gene_info, ancestral_seq, position, allele_after)`

**Purpose:** Classify a specific mutation as synonymous or non-synonymous.

**Input:**
- `gene_info` (dict): Gene metadata (start, end, strand)
- `ancestral_seq` (str): Ancestral gene sequence
- `position` (int): Contig position of mutation
- `allele_after` (str): Derived allele (forward strand)

**Output:**
- Dictionary with mutation details, or empty dict if analysis fails

**Algorithm:**
1. Extract ancestral codon using the `get_codon_from_site()` function
2. Handle strand orientation (complement allele for reverse strand)
3. Create derived codon by substituting the new allele
4. Translate both codons using Biopython
5. Compare amino acids to determine S vs NS

**Example output:**
```python
{
    "codon_before": "GTA",
    "codon_after": "ATA",
    "aa_before": "V",
    "aa_after": "I",
    "mutation_type": "NS"
}
```

#### `get_codon_from_site(position, gene_info, sequence)`

**Purpose:** Convert contig position to gene position and extract the containing codon.

**Input:**
- `position` (int): 0-based contig position
- `gene_info` (dict): Gene metadata (start, end, strand)
- `sequence` (list or str): Gene sequence

**Output:**
- Tuple of (codon_string, pos_in_gene) or (None, None) if invalid

**Coordinate conversion:**
- **Forward strand:** `pos_in_gene = position - (gene_start - 1)`
- **Reverse strand:** `pos_in_gene = (gene_end - 1) - position`

**Codon extraction:**
1. Validate position is within gene bounds
2. Find codon start: `(pos_in_gene // 3) × 3`
3. Extract 3-base codon from sequence
4. Validate codon is complete (length = 3)

---

## 5. Usage Guidelines

### 5.1 Required Input Files

#### Significant Sites File (from `p_value_summary.py`)

**Format:** Tab-separated values (TSV)

**Required columns:**
- `mag_id`: Metagenome-assembled genome identifier
- `contig`: Contig name
- `position`: 0-based position on contig
- `gene_id`: Gene identifier from Prodigal annotation
- `min_p_value`: Minimum raw p-value across tests
- `q_value`: FDR-corrected q-value (recommended for filtering)

**Optional columns:**
- `test_type`: Statistical test used (e.g., "two_sample_paired_tTest")
- `group_analyzed`: Group name if multiple groups were analyzed

**Example:**
```tsv
mag_id	contig	position	gene_id	test_type	min_p_value	q_value
MAG_001	contig_1	1020	gene_12345	two_sample_paired_tTest	0.001	0.02
MAG_001	contig_1	1021	gene_12345	two_sample_paired_tTest	0.003	0.04
MAG_001	contig_2	5500	gene_67890	two_sample_paired_tTest	0.0001	0.005
```

**Source:** Generated by [`p_value_summary.py`](../../../alleleflux/scripts/preprocessing/p_value_summary.py) in the AlleleFlux workflow.

#### Profile Data Files

**Format:** Tab-separated values (TSV), gzipped

**Location pattern:**
```
{profile_dir}/{sample_id}/{sample_id}_{mag_id}_profiled.tsv.gz
```

**Required columns:**
- `contig`: Contig name
- `position`: 0-based position on contig
- `ref_base`: Reference allele (forward strand)
- `gene_id`: Gene identifier
- `A`, `T`, `G`, `C`: Read counts for each nucleotide
- `coverage`: Total coverage at the position

**Example:**
```tsv
contig	position	ref_base	gene_id	A	T	G	C	coverage
contig_1	1020	G	gene_12345	5	2	150	3	160
contig_1	1021	T	gene_12345	3	148	4	5	160
```

**Source:** Generated by [`profile_mags.py`](../../../alleleflux/scripts/analysis/profile_mags.py) in the AlleleFlux workflow.

#### Prodigal FASTA File

**Format:** FASTA with Prodigal-format headers

**Header format:**
```
>gene_id # start # end # strand # partial # start_type # stop_codon
```

**Example:**
```fasta
>gene_12345 # 1000 # 1500 # 1 # 00 # ATG # TAA
ATGAAACGTAAAGCTTAG...
>gene_67890 # 2000 # 2400 # -1 # 00 # ATG # TAG
ATGGCATTATCCGGTAAG...
```

**Important notes:**
- Start and end are **1-based** positions
- Strand is 1 (forward) or -1 (reverse)
- Sequences are always 5' to 3' in the coding direction
- For reverse genes, the sequence is the reverse complement

**Source:** Generated by Prodigal during gene prediction step in AlleleFlux.

### 5.2 Command-Line Usage

#### Basic Usage

```bash
python dnds_from_timepoints.py \
    --significant_sites results/p_value_summary_results.tsv \
    --mag_ids MAG_001 MAG_002 MAG_003 \
    --p_value_column q_value \
    --p_value_threshold 0.05 \
    --ancestral_sample_id sample_T0 \
    --derived_sample_id sample_T1 \
    --profile_dir results/profiles/ \
    --prodigal_fasta annotations/genes.fasta \
    --outdir results/dnds/ \
    --prefix experiment1 \
    --cpus 4
```

#### Parameter Descriptions

**Required Parameters:**

- `--significant_sites`: Path to the TSV file from p_value_summary.py
- `--mag_ids`: One or more MAG identifiers to analyze (space-separated)
- `--ancestral_sample_id`: Sample ID representing the earlier timepoint
- `--derived_sample_id`: Sample ID representing the later timepoint
- `--profile_dir`: Directory containing profile subdirectories for each sample
- `--prodigal_fasta`: Path to Prodigal gene predictions FASTA
- `--outdir`: Output directory for results

**Optional Parameters:**

- `--p_value_column`: Column for filtering ["min_p_value" or "q_value", default: "q_value"]
  - Use "q_value" for FDR-corrected significance (recommended)
  - Use "min_p_value" for raw p-values
- `--p_value_threshold`: Significance threshold [default: 0.05]
- `--test-type`: Filter for specific statistical test
  - Options: "two_sample_unpaired_tTest", "two_sample_paired_tTest", "cmh", "lmm", etc.
- `--group-analyzed`: Filter for specific group name
- `--prefix`: Prefix for output filenames [default: "ancestral_dnds"]
- `--cpus`: Number of parallel processes [default: all available CPUs]

#### Advanced Usage Examples

**Example 1: Analyze single MAG with strict threshold**
```bash
python dnds_from_timepoints.py \
    --significant_sites results/significant_sites.tsv \
    --mag_ids MAG_001 \
    --p_value_column q_value \
    --p_value_threshold 0.01 \
    --ancestral_sample_id baseline \
    --derived_sample_id followup \
    --profile_dir profiles/ \
    --prodigal_fasta genes.fasta \
    --outdir results/
```

**Example 2: Filter by test type and group**
```bash
python dnds_from_timepoints.py \
    --significant_sites results/significant_sites.tsv \
    --mag_ids MAG_001 MAG_002 \
    --p_value_column q_value \
    --p_value_threshold 0.05 \
    --test-type two_sample_paired_tTest \
    --group-analyzed treatment_group \
    --ancestral_sample_id pre_treatment \
    --derived_sample_id post_treatment \
    --profile_dir profiles/ \
    --prodigal_fasta genes.fasta \
    --outdir results/treatment_analysis/
```

**Example 3: Parallel processing of many MAGs**
```bash
python dnds_from_timepoints.py \
    --significant_sites all_mags_significant.tsv \
    --mag_ids $(cat mag_list.txt) \
    --ancestral_sample_id T0 \
    --derived_sample_id T6 \
    --profile_dir profiles/ \
    --prodigal_fasta all_genes.fasta \
    --outdir dnds_results/ \
    --cpus 16
```

### 5.3 Output Files

For each MAG analyzed, the script creates a subdirectory and generates multiple output files with NG86-specific naming patterns.

#### Codon Events File (NG86)

**Filename:** `{prefix}_codon_events_ng86.tsv`

**Content:** Codon-level substitution events with NG86 path-averaged fractional counts

**Columns:**
- Position information: `mag_id`, `contig`, `position` (comma-separated if k>1), `gene_id`, `pos_in_gene`, `codon_start_index`, `position_in_codon` (comma-separated), `strand`
- Alleles: `major_allele_ancestral`, `major_allele_derived` (comma-separated if k>1)
- Codons and amino acids: `codon_before`, `codon_after`, `aa_before`, `aa_after`
- NG86 metrics: `k` (number of positions changed), `frac_S` (fractional synonymous count), `frac_N` (fractional non-synonymous count)
- Potential sites: `potential_S_sites_codon`, `potential_N_sites_codon`

**Example:**
```tsv
mag_id	contig	position	gene_id	codon_before	codon_after	aa_before	aa_after	k	frac_S	frac_N	potential_S_sites_codon	potential_N_sites_codon
MAG_001	contig_1	1020	gene_12345	GTA	ATA	V	I	1	0.0	1.0	0.67	2.33
MAG_001	contig_1	1021,1023	gene_12345	GTA	ACA	V	T	2	0.5	1.5	0.67	2.33
MAG_001	contig_2	5500,5501,5502	gene_67890	ATG	TTC	M	F	3	0.83	2.17	0.00	3.00
```

**Use cases:**
- Identify codons with multiple position changes (k>1)
- Verify NG86 path averaging calculations
- Understand fractional S/NS counts
- Map mutations to protein structure

**Key difference from previous version:** Reports codon-level events with fractional counts instead of per-site integer counts.

#### Gene Summary File (NG86)

**Filename:** `{prefix}_gene_summary_ng86.tsv`

**Content:** Per-gene dN/dS ratios calculated from fractional counts

**Columns:**
- `mag_id`: MAG identifier
- `gene_id`: Gene identifier
- `s_count`: Sum of fractional synonymous counts (`Σ frac_S`)
- `ns_count`: Sum of fractional non-synonymous counts (`Σ frac_N`)
- `potential_S_sites`: Total potential synonymous sites
- `potential_N_sites`: Total potential non-synonymous sites
- `dN`: Non-synonymous substitution rate
- `dS`: Synonymous substitution rate
- `dN_dS_ratio`: Gene-level dN/dS ratio

**Example:**
```tsv
mag_id	gene_id	s_count	ns_count	potential_S_sites	potential_N_sites	dN	dS	dN_dS_ratio
MAG_001	gene_12345	5.5	2.3	45.67	123.33	0.0186	0.1204	0.1546
MAG_001	gene_67890	0.83	8.17	23.00	89.00	0.0918	0.0361	2.5432
```

**Note:** `s_count` and `ns_count` may now be non-integer values due to path averaging.

**Use cases:**
- Identify genes under positive selection (dN/dS > 1)
- Identify highly constrained genes (dN/dS << 1)
- Compare selection patterns across genes

#### MAG Summary File (NG86)

**Filename:** `{prefix}_mag_summary_ng86.tsv`

**Content:** Aggregate statistics at the MAG level with fractional counts

**Columns:**
- `mag_id`: MAG identifier
- `total_s_count`: Sum of fractional synonymous counts across all genes
- `total_ns_count`: Sum of fractional non-synonymous counts across all genes
- `potential_S_sites`: Sum of potential S sites across all genes
- `potential_N_sites`: Sum of potential N sites across all genes
- `dN`: Overall non-synonymous rate
- `dS`: Overall synonymous rate
- `dN_dS_ratio`: MAG-level dN/dS ratio

**Example:**
```tsv
mag_id	total_s_count	total_ns_count	potential_S_sites	potential_N_sites	dN	dS	dN_dS_ratio
MAG_001	45.67	67.83	1234.50	3890.67	0.0174	0.0370	0.4705
```

**Use cases:**
- Compare evolutionary patterns across different MAGs
- Identify MAGs experiencing strong selection
- Aggregate-level interpretation of selection

#### Global dN/dS Summary File (NG86)

**Filename:** `{prefix}_global_dnds_ng86_summary.tsv`

**Content:** Single global dN/dS calculation with k-value distribution statistics

**Format:**
```tsv
Metric	Value
Total Potential Non-Synonymous Sites (n)	125.67
Total Potential Synonymous Sites (s)	34.33
Observed Non-Synonymous Substitutions (fractional sum)	15.83
Observed Synonymous Substitutions (fractional sum)	8.17
pN (NS / n)	0.1260
pS (S / s)	0.2380
Global dN/dS (pN/pS)	0.5294
Number of codon events (k=1)	8
Number of codon events (k=2)	3
Number of codon events (k=3)	1
Total codon events	12
```

**Important:** 
- Fractional counts properly reflect NG86 path averaging
- K-value distribution shows breakdown of single vs. multi-position changes
- Proper codon-level deduplication applied

**Use cases:**
- Overall assessment of selection at significant sites
- Understanding the distribution of k-values (single vs. multi-position changes)
- Comparison with per-gene ratios
- Publication-quality summary statistic

#### Codon-Specific Potential Sites File (NG86)

**Filename:** `{prefix}_codon_specific_sites_ng86.tsv`

**Content:** Detailed table showing potential S/N sites for each codon with a substitution

**Columns:**
- Position info: `mag_id`, `contig`, `position`, `gene_id`, `pos_in_gene`, `codon_start_index`, `position_in_codon`, `strand`
- Alleles and codons: `major_allele_ancestral`, `major_allele_derived`, `codon_before`, `codon_after`
- Translation: `aa_before`, `aa_after`
- NG86 metrics: `k`, `frac_S`, `frac_N`
- Codon potential: `potential_S_sites_codon`, `potential_N_sites_codon`

**Example:**
```tsv
contig	position	gene_id	codon_before	codon_after	k	frac_S	frac_N	potential_S_sites_codon	potential_N_sites_codon
contig_1	1020	gene_12345	GTA	ATA	1	0.0	1.0	0.67	2.33
contig_1	1021,1023	gene_12345	GTA	ACA	2	0.5	1.5	0.67	2.33
contig_2	5500,5501,5502	gene_67890	ATG	TTC	3	0.83	2.17	0.00	3.00
```

**Use cases:**
- Verify codon-level calculations
- Understand path averaging for multi-position changes
- Debug global dN/dS calculations
- Validate NG86 implementation

#### Reconstructed Ancestral Sequences File

**Filename:** `{prefix}_ancestral_orfs.ffn`

**Format:** FASTA

**Content:** Reconstructed ancestral sequences for all analyzed genes

**Example:**
```fasta
>ancestral_gene_12345 # 1000 # 1500 # 1 # 00 # ATG # TAA
ATGAAACGTAAAGCTTAGGCATTACCGGTAAGTCATAA
>ancestral_gene_67890 # 2000 # 2400 # -1 # 00 # ATG # TAG
ATGGCATTATCCGGTAAGTTAACGGTAAG
```

**Use cases:**
- Further downstream analysis (e.g., protein modeling)
- Validation of reconstruction accuracy
- Input for other evolutionary analysis tools

### 5.4 Interpreting Results

#### Understanding dN/dS Values

The dN/dS ratio provides insight into the selective pressures acting on protein-coding genes:

**dN/dS < 1: Purifying (Negative) Selection**
- Most non-synonymous mutations are deleterious
- Natural selection removes them from the population
- Indicates functional constraint
- Example: dN/dS = 0.2 means NS mutations occur at 20% the rate of S mutations

**dN/dS ≈ 1: Neutral Evolution**
- Non-synonymous and synonymous mutations occur at similar rates
- Mutations are neither beneficial nor harmful
- May indicate non-functional sequence regions
- Range typically 0.8 - 1.2 considered neutral

**dN/dS > 1: Positive (Adaptive) Selection**
- Non-synonymous mutations are favored
- Natural selection promotes amino acid changes
- Indicates active adaptation
- Example: dN/dS = 2.0 means NS mutations occur at twice the rate of S mutations

#### Understanding Fractional Counts (NG86)

**What are fractional counts?**

Fractional S/NS counts arise from NG86 path averaging when k > 1 (multiple positions changed in a codon):
- **k=1**: Counts are always integers (0 or 1 per position)
- **k=2**: Counts can be 0.0, 0.5, 1.0, 1.5, or 2.0 (averaging over 2 pathways)
- **k=3**: Counts can range from 0.0 to 3.0 in increments of 1/6 (averaging over 6 pathways)

**Example interpretation:**

If a codon event shows `frac_S = 1.5` and `frac_N = 0.5`:
- This is a k=2 change (two positions)
- Pathway 1: 2 synonymous, 0 non-synonymous
- Pathway 2: 1 synonymous, 1 non-synonymous
- Average: (2+1)/2 = 1.5 S, (0+1)/2 = 0.5 NS

**Why this matters:**

- Traditional per-site methods would count each position separately, ignoring path dependency
- NG86 recognizes that the order of mutations affects their classification
- Fractional counts represent the expected classification averaged over all possible orders
- More accurate for codons where multiple positions change

**In practice:**

- Most codon events (k=1) will have integer counts
- Fractional counts indicate multi-position changes requiring careful interpretation
- Global dN/dS properly sums these fractional values
- Gene-level summaries aggregate fractional counts across codons

#### Understanding k-Values

**k = number of positions that differ between ancestral and derived codons**

**k=1 (Single-position changes):**
- Most common scenario
- Straightforward interpretation
- No path dependency
- Direct S vs NS classification

**k=2 (Two-position changes):**
- Less common but biologically relevant
- Could indicate:
  - Two independent mutations close in time
  - Compensatory mutations
  - Rapid evolution
- NG86 averages over 2! = 2 pathways
- Classification depends on mutation order

**k=3 (Three-position changes):**
- Rare but possible scenarios:
  - Very rapid evolution
  - Recombination events
  - Long time intervals between samples
  - Multiple selective sweeps
- NG86 averages over 3! = 6 pathways
- Complex path dependencies

**Biological interpretation:**

- High proportion of k>1 events may indicate:
  - Strong selective pressure (multiple compensatory changes)
  - Long time intervals between timepoints
  - Possible recombination
  - Multiple sweeps at nearby positions

- Check k-value distribution in global summary:
  ```tsv
  Number of codon events (k=1)	150
  Number of codon events (k=2)	12
  Number of codon events (k=3)	2
  ```

- If k>1 events are common (>20%), consider:
  - Shortening time intervals between samples
  - Checking for recombination
  - Examining spatial distribution of k>1 events

#### Statistical Significance

**Important considerations:**

1. **Sample size matters:** Low numbers of substitutions can produce unreliable ratios
   - Look for genes with at least 10-20 total substitutions for robust estimates
   - Small sample sizes can produce extreme ratios due to stochastic effects

2. **Zero denominators:** If pS = 0, dN/dS is undefined (reported as NaN)
   - This can occur by chance in small datasets
   - May indicate very strong purifying selection on synonymous sites (rare)

3. **Multiple testing:** When examining many genes, apply multiple testing correction
   - Use FDR correction for declaring genes under selection
   - Consider the biological context

#### Biological Interpretation

**Purifying selection (dN/dS < 1):**
- Expected for most genes, especially housekeeping genes
- Lower values indicate stronger constraint
- Core metabolic genes typically show dN/dS < 0.3

**Genes under positive selection (dN/dS > 1):**
- Often involved in host-pathogen interactions
- Immune system genes
- Surface proteins
- Environmental adaptation genes

**Context-dependent interpretation:**
- Consider the biological role of the gene
- Compare to neutral expectation for the organism
- Look at spatial distribution of substitutions (clustered vs. distributed)

### 5.5 Common Use Cases

#### Use Case 1: Detecting Genes Under Selection

**Goal:** Identify genes experiencing adaptive evolution

**Approach:**
1. Run analysis with appropriate significance threshold (e.g., q_value ≤ 0.05)
2. Examine `{prefix}_gene_summary.tsv`
3. Filter for genes with dN_dS_ratio > 1.5 and sufficient substitutions
4. Investigate biological function using gene annotations

**Example R code for post-processing:**
```r
library(tidyverse)

# Load gene summary
gene_stats <- read_tsv("experiment1_MAG_001_gene_summary.tsv")

# Filter for positive selection
adaptive_genes <- gene_stats %>%
  filter(dN_dS_ratio > 1.5,
         (s_count + ns_count) >= 10)  # Minimum 10 substitutions

# Order by dN/dS ratio
adaptive_genes %>%
  arrange(desc(dN_dS_ratio)) %>%
  print(n = 20)
```

#### Use Case 2: Comparing Selection Patterns Across MAGs

**Goal:** Compare evolutionary pressures on orthologous genes across different MAGs

**Approach:**
1. Run analysis for multiple MAGs
2. Merge gene summaries using orthology information
3. Compare dN/dS ratios for the same gene across MAGs

**Example:**
```bash
# Run for all MAGs
python dnds_from_timepoints.py \
    --mag_ids MAG_001 MAG_002 MAG_003 \
    --significant_sites all_sites.tsv \
    --ancestral_sample_id T0 \
    --derived_sample_id T6 \
    --profile_dir profiles/ \
    --prodigal_fasta genes.fasta \
    --outdir comparative_analysis/

# Combine results in R
library(tidyverse)
gene_files <- list.files("comparative_analysis/", 
                         pattern = "_gene_summary.tsv$", 
                         recursive = TRUE, 
                         full.names = TRUE)
all_genes <- map_df(gene_files, read_tsv)

# Compare across MAGs
all_genes %>%
  group_by(gene_id) %>%
  summarize(mean_dnds = mean(dN_dS_ratio, na.rm = TRUE),
            sd_dnds = sd(dN_dS_ratio, na.rm = TRUE),
            n_mags = n())
```

#### Use Case 3: Identifying Adaptive Evolution in Specific Lineages

**Goal:** Track selection patterns over multiple timepoints

**Approach:**
1. Run pairwise comparisons between consecutive timepoints
2. Track how dN/dS changes over time for specific genes
3. Identify when positive selection begins or ends

**Example:**
```bash
# T0 vs T1
python dnds_from_timepoints.py \
    --ancestral_sample_id T0 --derived_sample_id T1 \
    --prefix T0_T1 --outdir results/

# T1 vs T2
python dnds_from_timepoints.py \
    --ancestral_sample_id T1 --derived_sample_id T2 \
    --prefix T1_T2 --outdir results/

# T2 vs T3
python dnds_from_timepoints.py \
    --ancestral_sample_id T2 --derived_sample_id T3 \
    --prefix T2_T3 --outdir results/
```

### 5.6 Potential Pitfalls and Limitations

#### Low Coverage Regions

**Problem:** Low sequencing coverage can lead to:
- Inaccurate allele frequency estimates
- Incorrect identification of major alleles
- False positive or false negative substitutions

**Solution:**
- Use AlleleFlux's quality filters before dN/dS analysis
- Ensure adequate coverage (typically ≥10× minimum)
- Consider coverage information when interpreting results

#### Multiple Substitutions in Same Codon

**Current behavior (NG86):** 
- Multiple positions changing in the same codon are handled as a single codon-level event
- NG86 path averaging properly accounts for all k! mutational pathways
- Fractional S/NS counts reflect averaged classifications across pathways
- Stop codon pathways are automatically excluded

**What this means:**
- If positions 0 and 1 of a codon both change, reported as one k=2 event
- Fractional counts (e.g., 1.5 NS, 0.5 S) represent pathway averaging
- No need for manual interpretation of multiple adjacent sites

**Check for such cases:**
- Look at `k` column in codon events file
- k=1: Single position changed
- k=2: Two positions changed (2 pathways averaged)
- k=3: Three positions changed (6 pathways averaged)

**Interpretation:**
- High k-values may indicate:
  - Strong selection pressure (compensatory mutations)
  - Longer time intervals
  - Rapid evolution
  - Possible recombination

**Note:** This represents a **breaking change** from the previous per-site implementation. The old version reported each position separately with integer counts; the new version reports codon-level events with fractional counts and proper NG86 path averaging.

#### Stop Codons

**Problem:** Mutations creating or removing stop codons are special cases

**Handled by the script:**
- Stop codons are represented as "*" in translations
- Changes
 to/from stop codons are classified as non-synonymous
- May indicate pseudogenization or reading frame errors

**Recommendation:**
- Review genes with stop codon mutations manually
- Consider filtering out premature stop codons if they indicate sequencing/assembly errors

#### Partial Codons at Gene Boundaries

**Problem:** Genes whose length is not a multiple of 3 have incomplete codons

**Handled by the script:**
- Partial codons are automatically skipped with a warning
- Only complete 3-base codons are analyzed
- Logged warnings indicate affected genes

**Impact:**
- Minimal effect on most analyses
- May slightly reduce potential site counts for very short genes

#### Recombination and Horizontal Gene Transfer

**Problem:** The Nei-Gojobori method assumes vertical inheritance

**Limitation:**
- Recombination can violate assumptions about independent sites
- Horizontal gene transfer creates complex evolutionary histories
- May lead to misleading dN/dS interpretations

**Recommendation:**
- Consider recombination detection tools for genes with unusual patterns
- Be cautious with genes known to undergo frequent HGT
- Interpret results in the context of known biology

#### Sample Size Considerations

**Guideline for minimum data requirements:**

- **Per gene:** At least 5-10 substitutions for reliable estimates
- **For global dN/dS:** At least 50-100 total substitutions across genes
- **Coverage:** Minimum 10× coverage recommended, 20× preferred

**Warning signs of insufficient data:**
- Extreme dN/dS values (>10 or <0.01)
- Many genes with dN/dS = NaN (due to pS = 0)
- High variance between replicate analyses

#### Time Between Samples

**Consideration:** The timespan between ancestral and derived samples affects:

- Number of substitutions observed (longer time = more substitutions)
- Multiple hits at the same site (saturation effects)
- Reliability of ancestral state reconstruction

**Recommendations:**
- Ideal timespan: Sufficient for mutations to accumulate but not so long that multiple substitutions occur at the same site
- For bacterial evolution experiments: Several hundred to thousand generations
- Very short timespans may yield too few substitutions for robust analysis
- Very long timespans may violate single-substitution assumptions

---

## 6. Integration with AlleleFlux Workflow

### Workflow Position

The dN/dS analysis script is designed to run **after** the core AlleleFlux statistical analysis:

```
AlleleFlux Workflow:
1. Quality control and MAG profiling
2. Statistical significance testing
3. p_value_summary.py ← Generates input for dN/dS
4. dnds_from_timepoints.py ← This script
5. Interpretation and visualization
```

### Required Upstream Steps

Before running dN/dS analysis, ensure these steps are complete:

1. **Quality Control:** [`quality_control.py`](../../../alleleflux/scripts/preprocessing/quality_control.py)
2. **MAG Profiling:** [`profile_mags.py`](../../../alleleflux/scripts/analysis/profile_mags.py)
3. **Statistical Testing:** Various test scripts (e.g., [`two_sample_paired.py`](../../../alleleflux/scripts/statistics/two_sample_paired.py))
4. **P-value Summary:** [`p_value_summary.py`](../../../alleleflux/scripts/preprocessing/p_value_summary.py)

### Typical Workflow Integration

**Example Snakemake rule integration:**

```python
rule dnds_analysis:
    input:
        significant_sites="results/p_value_summary_results.tsv",
        profile_dir="results/profiles/",
        prodigal_fasta="annotations/genes.fasta"
    output:
        "results/dnds/{mag_id}/{prefix}_all_substitutions.tsv",
        "results/dnds/{mag_id}/{prefix}_gene_summary.tsv",
        "results/dnds/{mag_id}/{prefix}_global_dnds_summary.tsv"
    params:
        mag_ids=config["mag_ids"],
        ancestral=config["ancestral_sample"],
        derived=config["derived_sample"]
    shell:
        """
        python dnds_from_timepoints.py \
            --significant_sites {input.significant_sites} \
            --mag_ids {params.mag_ids} \
            --ancestral_sample_id {params.ancestral} \
            --derived_sample_id {params.derived} \
            --profile_dir {input.profile_dir} \
            --prodigal_fasta {input.prodigal_fasta} \
            --outdir results/dnds/ \
            --cpus {threads}
        """
```

---

## 7. Performance and Optimization

### Computational Complexity

**Time complexity:**
- Pre-computation of codon cache: O(1) - done once at startup
- Per-gene potential sites: O(L) where L = gene length
- Substitution analysis: O(S) where S = number of significant sites
- Overall: Linear in the number of sites and genes

**Memory usage:**
- Profile data for two timepoints per MAG
- Reconstructed sequences for all genes
- Typically requires 1-4 GB RAM for most datasets

### Parallelization Strategy

The script uses **multiprocessing at the MAG level:**

```python
# Each MAG is processed independently
with mp.Pool(processes=num_cpus) as pool:
    pool.map(process_mag, mag_ids)
```

**Benefits:**
- Near-linear speedup with number of CPUs
- No communication overhead between processes
- Scales well to large numbers of MAGs

**Recommendations:**
- Use `--cpus` equal to the number of MAGs (up to available cores)
- For very large datasets, process MAGs in batches

### Optimization Tips

1. **Pre-filter significant sites:** More stringent thresholds reduce computation time
2. **Use compressed profiles:** Gzipped TSV files save disk I/O
3. **Process similar MAGs together:** Similar coverage patterns benefit from caching
4. **Monitor memory:** Large genomes may require more RAM per process

---

## 8. Troubleshooting

### Common Errors and Solutions

#### Error: "Required columns not found in significant_sites file"

**Cause:** Input file from p_value_summary.py is missing expected columns

**Solution:**
- Verify the input file contains: `mag_id`, `contig`, `position`, `gene_id`, and your chosen p-value column
- Check that p_value_summary.py ran successfully
- Ensure you're specifying the correct `--p_value_column`

#### Error: "Profile file(s) for MAG not found"

**Cause:** Missing or incorrectly named profile files

**Solution:**
- Verify profile directory structure: `{profile_dir}/{sample_id}/{sample_id}_{mag_id}_profiled.tsv.gz`
- Check that sample IDs match exactly (case-sensitive)
- Ensure profile_mags.py completed successfully for both samples

#### Warning: "Reference base mismatch"

**Cause:** Discrepancy between Prodigal reference and profile reference

**Possible reasons:**
- Different reference assemblies used
- Strand orientation issues
- Data processing errors

**Solution:**
- Usually safe to ignore if infrequent (script uses Prodigal base)
- If frequent, verify that profiles and Prodigal annotations are from the same assembly
- Check for strand annotation errors

#### Output: Many genes with dN/dS = NaN

**Cause:** No synonymous substitutions observed (pS = 0)

**Solutions:**
- This is expected for genes with very few substitutions
- Use more relaxed significance thresholds to include more sites
- Combine data across multiple timepoint comparisons
- Focus interpretation on genes with sufficient substitutions

#### Performance: Script runs slowly

**Solutions:**
- Reduce number of significant sites using stricter thresholds
- Increase `--cpus` parameter
- Process MAGs in smaller batches
- Ensure adequate RAM is available (avoid swapping)

### Validation Steps

**To verify correct execution:**

1. **Check log output:**
   - Number of sites loaded and filtered
   - Number of genes reconstructed
   - Number of substitutions identified

2. **Inspect output files:**
   - All expected files are created
   - Files contain data (not empty)
   - Gene and MAG summaries have reasonable values

3. **Sanity checks:**
   - Total substitutions = S_count + NS_count
   - Potential S + Potential N ≈ 3 × number of codons
   - dN/dS values are within reasonable range (0.01 - 10)

4. **Spot check specific genes:**
   - Manually verify a few substitutions
   - Check codon translations
   - Confirm S vs NS classifications

---

## 9. References and Further Reading

### Primary References

1. **Nei, M., & Gojobori, T. (1986).** Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. *Molecular Biology and Evolution*, 3(5), 418-426.
   - Original description of the method implemented in this script

2. **Yang, Z., & Nielsen, R. (2000).** Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. *Molecular Biology and Evolution*, 17(1), 32-43.
   - Improvements and extensions to dN/dS methods

3. **Kryazhimskiy, S., & Plotkin, J. B. (2008).** The population genetics of dN/dS. *PLoS Genetics*, 4(12), e1000304.
   - Interpretation of dN/dS in population genetics contexts

### Related Tools

- **PAML (Phylogenetic Analysis by Maximum Likelihood):** More sophisticated dN/dS analysis with phylogenetic models
- **HyPhy:** Flexible framework for evolutionary analysis
- **KaKs_Calculator:** Alternative implementation of various dN/dS methods

### AlleleFlux Documentation

- [Getting Started](../getting_started/overview.rst)
- [Input Preparation](../usage/input_preparation.rst)
- [Interpreting Results](../usage/interpreting_results.rst)

---

## 10. Example Analysis Workflow

### Complete Step-by-Step Example

Here's a complete workflow from raw data to dN/dS results:

```bash
# Step 1: Run AlleleFlux core workflow (quality control, profiling, statistics)
snakemake --configfile config.yml --cores 16

# Step 2: Generate p-value summary
python alleleflux/scripts/preprocessing/p_value_summary.py \
    --score_dir results/scores/ \
    --output results/p_value_summary_results.tsv

# Step 3: Run dN/dS analysis for multiple MAGs
python alleleflux/scripts/evolution/dnds_from_timepoints.py \
    --significant_sites results/p_value_summary_results.tsv \
    --mag_ids MAG_001 MAG_002 MAG_003 MAG_004 MAG_005 \
    --p_value_column q_value \
    --p_value_threshold 0.05 \
    --test-type two_sample_paired_tTest \
    --ancestral_sample_id T0_baseline \
    --derived_sample_id T6_followup \
    --profile_dir results/profiles/ \
    --prodigal_fasta annotations/all_genes.fasta \
    --outdir results/dnds_analysis/ \
    --prefix experiment1 \
    --cpus 8

# Step 4: Examine results
ls -lh results/dnds_analysis/MAG_001/

# Step 5: Summarize across all MAGs (example R script)
Rscript analyze_dnds_results.R results/dnds_analysis/
```

### Example R Script for Post-Processing

```r
#!/usr/bin/env Rscript
library(tidyverse)

# Load all gene summaries
gene_files <- list.files("results/dnds_analysis/", 
                         pattern = "_gene_summary.tsv$",
                         recursive = TRUE,
                         full.names = TRUE)

all_genes <- map_df(gene_files, read_tsv)

# Summary statistics
cat("=== Overall dN/dS Summary ===\n")
all_genes %>%
  summarize(
    mean_dnds = mean(dN_dS_ratio, na.rm = TRUE),
    median_dnds = median(dN_dS_ratio, na.rm = TRUE),
    n_genes = n(),
    n_positive_selection = sum(dN_dS_ratio > 1.5, na.rm = TRUE)
  ) %>%
  print()

# Genes under positive selection
cat("\n=== Top 20 Genes Under Positive Selection ===\n")
all_genes %>%
  filter(!is.na(dN_dS_ratio), dN_dS_ratio > 1) %>%
  arrange(desc(dN_dS_ratio)) %>%
  select(mag_id, gene_id, dN_dS_ratio, ns_count, s_count) %>%
  head(20) %>%
  print()

# Distribution plot
ggplot(all_genes, aes(x = dN_dS_ratio)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "Distribution of dN/dS Ratios",
       x = "dN/dS (log scale)",
       y = "Number of Genes") +
  theme_minimal()

ggsave("dnds_distribution.pdf", width = 8, height = 6)
```

---

## Summary

The [`dnds_from_timepoints.py`](../../../alleleflux/scripts/evolution/dnds_from_timepoints.py) script provides a robust implementation of the Nei-Gojobori method for calculating dN/dS ratios from AlleleFlux results. Key features include:

- **Accurate ancestral reconstruction** from allele frequency data
- **Proper handling of coordinate systems** and strand orientations
- **Anti-double-counting mechanisms** for global dN/dS calculation
- **Comprehensive output** at multiple levels (site, gene, MAG, global)
- **Efficient parallel processing** for multiple MAGs

The script enables researchers to:
- Identify genes under positive selection
- Detect purifying selection on essential genes  
- Compare evolutionary patterns across lineages
- Quantify selection pressures in microbial populations

For questions or issues, please refer to the AlleleFlux repository or contact the development team.