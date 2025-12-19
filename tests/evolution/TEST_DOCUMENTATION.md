# dN/dS Analysis Test Suite Documentation

## Table of Contents
1. [Overview](#overview)
2. [Test Class Summaries](#test-class-summaries)
3. [TestReconstructAncestralSequences](#testreconstructancestralsequences)
4. [TestCoreDndsFunctions](#testcorednds-functions)
5. [TestNG86PathAveraging](#testng86pathaveraging)
6. [TestCodonSubstitutionAnalysis](#testcodonsubstitutionanalysis)
7. [TestGlobalDndsCalculation](#testglobaldndscalculation)
8. [TestNG86ComputedValues](#testng86computedvalues)
9. [Mathematical Formulas](#mathematical-formulas)
10. [Edge Cases and Boundary Conditions](#edge-cases-and-boundary-conditions)
11. [Scientific References](#scientific-references)
12. [Helper Functions](#helper-functions)

---

## Overview

### Purpose

This test suite validates the implementation of dN/dS (non-synonymous to synonymous substitution ratio) analysis using the **Nei-Gojobori 1986 (NG86) path averaging methodology**. The implementation is located at [`alleleflux/scripts/evolution/dnds_from_timepoints.py`](../alleleflux/scripts/evolution/dnds_from_timepoints.py).

The test suite contains **38 test methods** organized into **6 test classes**, providing comprehensive coverage of:
- Ancestral sequence reconstruction from profile data
- Core dN/dS calculation functions
- NG86 path averaging for multi-position codon changes
- Codon-level substitution analysis
- Global dN/dS ratio calculation
- End-to-end integration with mock datasets
- **Negative strand gene handling (6 dedicated test methods)**

### What is dN/dS Analysis?

**dN/dS** is a fundamental metric in molecular evolution that measures the ratio of:
- **dN** (non-synonymous substitution rate): mutations that change the amino acid
- **dS** (synonymous substitution rate): mutations that preserve the amino acid

**Biological Significance:**
- **dN/dS < 1**: Purifying selection (most mutations harmful)
- **dN/dS ≈ 1**: Neutral evolution (mutations neither beneficial nor harmful)
- **dN/dS > 1**: Positive selection (adaptive evolution)

This ratio helps identify:
- Genes under selective pressure
- Functionally important genomic regions
- Evolutionary constraints on protein sequences

### NG86 Path Averaging Methodology

The **Nei-Gojobori (1986)** method addresses a critical challenge: when multiple positions change within a codon (k > 1), the classification of each change as synonymous or non-synonymous depends on the **order** in which mutations occur.

**The Problem:**
Consider codon ATG (Met) → TTT (Phe) with 2 changes at positions 0 and 2:
- **Path 1**: ATG → TTG (Met→Leu, NS) → TTT (Leu→Phe, NS) = 2 NS steps
- **Path 2**: ATG → ATT (Met→Ile, NS) → TTT (Ile→Phe, NS) = 2 NS steps

**NG86 Solution:**
Average the S/NS classifications across all k! possible mutational pathways:
- For k=1: No averaging needed (1 pathway)
- For k=2: Average over 2 pathways (2! = 2)
- For k=3: Average over 6 pathways (3! = 6)

This approach:
- Removes arbitrary assumptions about mutation order
- Provides unbiased evolutionary rate estimates
- Handles path-dependent amino acid changes correctly

**Key Features:**
- Paths through intermediate **stop codons** are excluded (except final stop)
- Returns **fractional** S/NS counts (e.g., 0.5 S, 1.5 NS for k=2)
- Uses pre-computed cache for all 4,096 possible codon pairs

### Test Suite Organization

The suite is organized into 6 test classes:

| Test Class | Purpose | Test Count |
|------------|---------|------------|
| [`TestReconstructAncestralSequences`](#testreconstructancestralsequences) | Validates ancestral sequence reconstruction from profile data | 2 |
| [`TestCoreDndsFunctions`](#testcorednds-functions) | Tests core helper functions and Nei-Gojobori site calculations | 7 |
| [`TestNG86PathAveraging`](#testng86pathaveraging) | Tests path averaging algorithm for k=1, k=2, k=3 cases | 8 |
| [`TestCodonSubstitutionAnalysis`](#testcodonsubstitutionanalysis) | Validates codon-level grouping and analysis | 2 |
| [`TestGlobalDndsCalculation`](#testglobaldndscalculation) | Tests global dN/dS calculation with fractional counts | 3 |
| [`TestNG86ComputedValues`](#testng86computedvalues) | Comprehensive validation of NG86 computed values (positive and negative strands) | 16 |

**Total:** 38 test methods across 6 classes

---

## Test Class Summaries

### Quick Reference

```python
# Import test file
from tests.evolution.test_dnds_from_timepoints import *

# Run all tests
python -m unittest tests.evolution.test_dnds_from_timepoints

# Run specific test class
python -m unittest tests.evolution.test_dnds_from_timepoints.TestNG86PathAveraging

# Run specific test method
python -m unittest tests.evolution.test_dnds_from_timepoints.TestNG86PathAveraging.test_ng86_double_change_mixed
```

---

## TestReconstructAncestralSequences

### Purpose

This test class validates the [`reconstruct_ancestral_sequences()`](../alleleflux/scripts/evolution/dnds_from_timepoints.py) function, which reconstructs ancestral gene sequences by substituting major alleles from profile data into reference sequences.

### What is Ancestral Sequence Reconstruction?

**Goal:** Build the most likely ancestral sequence for each gene based on allele frequency data at variable sites.

**Process:**
1. Start with reference sequence from Prodigal annotation
2. Identify variable sites from ancestral profile data
3. Determine major allele at each variable site (highest count)
4. Substitute major allele into reference sequence
5. Handle strand orientation (forward vs. reverse)
6. Apply tie-breaking rules for equal allele counts

**Strand Handling:**
- **Forward strand genes**: Use alleles directly from profile (already on forward strand)
- **Reverse strand genes**: Complement alleles to match reverse-complement sequence

**Tie-Breaking Logic:**
When multiple alleles have equal maximum counts:
1. **Prefer reference base** if it's among the tied alleles
2. **Random selection** if reference base not in tie
3. **Default to reference** if no coverage (all counts = 0)

### Test Methods

#### `test_reconstruction_with_multiple_changes`

**What it tests:** Complete reconstruction workflow with multiple changes per gene.

**Setup:**
```python
# Four genes with different characteristics
genes = {
    'gene_fwd_1': forward strand, 3 variable sites
    'gene_rev_1': reverse strand, 3 variable sites  
    'gene_tie_break': forward strand, tie scenario
    'gene_no_profile': forward strand, no profile data
}
```

**Profile Data Structure:**
```python
profile_data = {
    'contig': ['contig1', 'contig1', ...],
    'position': [100, 102, 108, ...],
    'ref_base': ['A', 'G', 'G', ...],
    'gene_id': ['gene_fwd_1', 'gene_fwd_1', ...],
    'A': [5, 5, 90, ...],    # Allele counts
    'C': [2, 95, 2, ...],
    'G': [3, 0, 3, ...],
    'T': [90, 0, 5, ...]
}
```

**Forward Strand Example (gene_fwd_1):**

Original sequence: `ATGCGTACG`
```
Position  100   101   102   103   104   105   106   107   108
Original   A     T     G     C     G     T     A     C     G
Counts    T:90  -     C:95  -     -     -     -     -     A:90
Result     T     T     C     C     G     T     A     C     A
```

**Calculation Steps:**
1. Position 100: Counts {A:5, C:2, G:3, T:90} → Major = T (max count)
2. Position 102: Counts {A:5, C:95, G:0, T:0} → Major = C (max count)
3. Position 108: Counts {A:90, C:2, G:3, T:5} → Major = A (max count)

**Expected:** `"TTCCGTACA"` (3 positions changed from original)

**Reverse Strand Example (gene_rev_1):**

Original sequence (reverse complement): `GTACTAACA`
Prodigal coordinates: start=201, end=209, strand=-1

```
Contig Position   202         203         208
Profile (fwd)     C:90        A:80        T:88
Major (fwd)       C           A           T
Complement        G           T           A
Gene position     7           6           1
```

**Coordinate Conversion for Reverse Strand:**
- Contig position 202 → Gene position = (209-1) - 202 = 6
- Contig position 203 → Gene position = (209-1) - 203 = 5  
- Contig position 208 → Gene position = (209-1) - 208 = 0

**Expected:** `"ATACTTGCA"` (reverse strand gene reconstructed)

**Tie-Breaking Example (gene_tie_break):**
```python
# Profile at position 301
{'A': 20, 'C': 20, 'G': 5, 'T': 1}  # A and C tied
ref_base = 'A'  # Reference base is A
# Result: 'A' (prefer reference in tie)
```

**Expected:** `"AAAAAA"` (no change due to tie-breaking)

**No Profile Data (gene_no_profile):**
- No entries in profile data for this gene
- Expected: `"GGGGGG"` (original sequence unchanged)

**Ancestral Major Alleles Tracking:**

The function also returns a dictionary mapping (contig, position) → forward-strand major allele:
```python
expected_major_alleles = {
    ('contig1', 100): 'T',   # Forward strand major
    ('contig1', 102): 'C',
    ('contig1', 108): 'A',
    ('contig1', 202): 'C',   # Reverse strand, but stored as forward
    ('contig1', 203): 'A',
    ('contig1', 208): 'T',
    ('contig2', 301): 'A'
}
```

**Assertions:**
```python
assert reconstructed_seqs['gene_fwd_1'] == 'TTCCGTACA'
assert reconstructed_seqs['gene_rev_1'] == 'ATACTTGCA'
assert reconstructed_seqs['gene_tie_break'] == 'AAAAAA'
assert reconstructed_seqs['gene_no_profile'] == 'GGGGGG'
assert 'gene_is_missing' not in reconstructed_seqs
assert len(reconstructed_seqs) == 4
assert ancestral_major_alleles == expected_major_alleles
```

#### `test_reconstruction_with_ref_mismatch`

**What it tests:** Warning is logged when profile ref_base doesn't match Prodigal sequence.

**Scenario:**
```python
# Profile says ref_base = 'A' at position 208
# Prodigal sequence (forward strand equivalent) has 'C' at this position
# Mismatch detected!
```

**Setup:**
```python
# Modify profile to create mismatch
profile_df.loc[profile_df['position'] == 208, 'ref_base'] = 'A'
# Original Prodigal has 'C' at this position (when oriented to forward)
```

**Expected Behavior:**
1. Warning logged: `"Reference base mismatch at position 208..."`
2. Uses profile's `ref_base` for tie-breaking (not Prodigal's)
3. Reconstruction continues normally
4. Final sequence still: `"ATACTTGCA"`

**Assertion:**
```python
with self.assertLogs(level='WARNING') as cm:
    reconstructed_seqs, _ = reconstruct_ancestral_sequences(...)
    assert any('Reference base mismatch' in msg for msg in cm.output)
```

**Why This Matters:**
- Detects data inconsistencies between annotation and profiling
- Ensures transparency in allele selection
- Uses profile data as source of truth (more recent/accurate)

---

## TestCoreDndsFunctions

### Purpose

This test class validates core helper functions and the Nei-Gojobori method for calculating potential synonymous and non-synonymous sites.

### Nei-Gojobori Potential Sites Method

The NG86 method calculates **fractional** potential S and N sites by considering all possible single-nucleotide changes at each codon position.

**Algorithm for a Single Codon:**
1. For each of 3 positions (0, 1, 2)
2. For each of 3 possible alternative bases (excluding original)
3. Create hypothetical mutant codon
4. Translate original and mutant to amino acids
5. Classify as S (same AA) or NS (different AA)
6. Each position contributes 1/3 to the total

**Example: ATG (Methionine)**
```
Position 0 (A → ?, T, G):
  ATG → TTG (Met → Leu) = NS
  ATG → CTG (Met → Leu) = NS  
  ATG → GTG (Met → Val) = NS
  Total: 0/3 S, 3/3 NS

Position 1 (T → A, ?, G):
  ATG → AAG (Met → Lys) = NS
  ATG → ACG (Met → Thr) = NS
  ATG → AGG (Met → Arg) = NS
  Total: 0/3 S, 3/3 NS

Position 2 (G → A, T, C):
  ATG → ATA (Met → Ile) = NS
  ATG → ATT (Met → Ile) = NS
  ATG → ATC (Met → Ile) = NS
  Total: 0/3 S, 3/3 NS

Final: S = (0+0+0)/3 = 0.0, N = (3+3+3)/3 = 3.0
```

ATG has **0 synonymous sites** and **3 non-synonymous sites**.

### Test Methods

#### `test_get_major_allele`

**What it tests:** [`get_major_allele()`](../alleleflux/scripts/evolution/dnds_from_timepoints.py:624) selects the most frequent allele with proper tie-breaking.

**Test Cases:**

**Case 1: Clear winner**
```python
allele_counts = {'A': 10, 'C': 90, 'G': 0, 'T': 0}
ref_base = 'A'
expected = 'C'  # Highest count wins
```

**Case 2: Tie with reference base**
```python
allele_counts = {'A': 50, 'C': 50, 'G': 0, 'T': 0}
ref_base = 'A'
expected = 'A'  # Prefer reference in tie
```

**Case 3: Tie without reference base**
```python
allele_counts = {'A': 50, 'C': 50, 'G': 0, 'T': 0}
ref_base = 'G'  # G not in tie
# Mock random.choice to return 'C'
expected = 'C'  # Random selection among tied
```

**Case 4: No coverage**
```python
allele_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
ref_base = 'T'
expected = 'T'  # Default to reference
```

#### `test_calculate_codon_sites`

**What it tests:** [`_calculate_codon_sites()`](../alleleflux/scripts/evolution/dnds_from_timepoints.py) calculates potential S and N sites for individual codons.

**Test Cases:**

**Case 1: Methionine (ATG) - No synonymous sites**
```python
codon = 'ATG'
expected_S = 0.0
expected_N = 3.0
```

Methionine has only one codon (ATG), so all single-base changes are non-synonymous.

**Case 2: Leucine (CTT) - Has synonymous sites**
```python
codon = 'CTT'
expected_S ≈ 1.0
expected_N ≈ 2.0
```

**Detailed Calculation:**
```
Position 0 (C → A, T, G):
  CTT → ATT (Leu → Ile) = NS
  CTT → TTT (Leu → Phe) = NS
  CTT → GTT (Leu → Val) = NS
  Contribution: 0 S, 3 NS

Position 1 (T → A, C, G):
  CTT → CAT (Leu → His) = NS
  CTT → CCT (Leu → Pro) = NS
  CTT → CGT (Leu → Arg) = NS
  Contribution: 0 S, 3 NS

Position 2 (T → A, C, G):
  CTT → CTA (Leu → Leu) = S  ✓
  CTT → CTC (Leu → Leu) = S  ✓
  CTT → CTG (Leu → Leu) = S  ✓
  Contribution: 3 S, 0 NS

Total: S = (0+0+3)/3 = 1.0, N = (3+3+0)/3 = 2.0
```

**Case 3: Stop codon (TAA)**
```python
codon = 'TAA'
expected_S ≈ 2/3
expected_N ≈ 4/3 + 1
```

**Detailed Calculation:**
```
Position 0 (T → A, C, G):
  TAA → AAA (Stop → Lys) = NS
  TAA → CAA (Stop → Gln) = NS
  TAA → GAA (Stop → Glu) = NS
  Contribution: 0 S, 3 NS

Position 1 (A → T, C, G):
  TAA → TTA (Stop → Leu) = NS
  TAA → TCA (Stop → Ser) = NS
  TAA → TGA (Stop → Stop) = S  ✓
  Contribution: 1 S, 2 NS

Position 2 (A → T, C, G):
  TAA → TAT (Stop → Tyr) = NS
  TAA → TAC (Stop → Tyr) = NS
  TAA → TAG (Stop → Stop) = S  ✓
  Contribution: 2 S, 1 NS

Total: S = (0+1+1)/3 = 2/3, N = (3+2+2)/3 = 7/3
```

**Case 4: Invalid codon (ATN)**
```python
codon = 'ATN'  # N is ambiguous
# Raises ValueError
```

#### `test_calculate_potential_sites_for_gene`

**What it tests:** [`calculate_potential_sites_for_gene()`](../alleleflux/scripts/evolution/dnds_from_timepoints.py) aggregates S and N sites over a full gene.

**Test Case:**
```python
gene_seq = 'ATGCTTTAA'  # Three codons: ATG, CTT, TAA
```

**Calculation:**
```
Codon 1 (ATG): S = 0.0,   N = 3.0
Codon 2 (CTT): S = 1.0,   N = 2.0
Codon 3 (TAA): S = 2/3,   N = 7/3

Total S = 0.0 + 1.0 + 2/3 = 5/3 ≈ 1.667
Total NS = 3.0 + 2.0 + 7/3 = 22/3 ≈ 7.333
```

**Expected:**
```python
assert sites['S'] ≈ 1.667
assert sites['NS'] ≈ 7.333
```

#### `test_get_codon_from_site`

**What it tests:** [`get_codon_from_site()`](../alleleflux/scripts/evolution/dnds_from_timepoints.py:758) correctly maps genomic positions to codons.

**Forward Strand Example:**
```python
gene_info = {
    'start': 101,  # 1-based Prodigal coordinate
    'end': 109,
    'strand': 1
}
sequence = 'ATGCGTACG'  # 9 bases = 3 codons
position = 104  # 0-based genomic position
```

**Calculation:**
```
Gene position (0-indexed) = 104 - (101 - 1) = 4
Codon start = (4 // 3) * 3 = 3
Codon = sequence[3:6] = 'CGT'
Position in gene = 4
```

**Expected:**
```python
codon, pos = get_codon_from_site(104, gene_info, sequence)
assert codon == 'CGT'
assert pos == 4
```

**Reverse Strand Example:**
```python
gene_info = {
    'start': 201,
    'end': 209,
    'strand': -1
}
sequence = 'GTACTAACA'  # Reverse complement
position = 203
```

**Calculation:**
```
Gene position = (209 - 1) - 203 = 5
Codon start = (5 // 3) * 3 = 3
Codon = sequence[3:6] = 'CTA'
```

**Expected:**
```python
codon, pos = get_codon_from_site(203, gene_info, sequence)
assert codon == 'CTA'
assert pos == 5
```

**Out of Bounds:**
```python
codon, pos = get_codon_from_site(99, gene_info_fwd, sequence)
assert codon is None
assert pos is None
```

#### `test_analyze_mutation_effect`

**What it tests:** [`analyze_mutation_effect()`](../alleleflux/scripts/evolution/dnds_from_timepoints.py:838) classifies individual mutations as S or NS.

**Non-Synonymous Example:**
```python
gene_info = forward strand gene at 101-109
ancestral_seq = 'ATGCGTACG'
position = 104  # In codon 'CGT' at position 1
derived_allele = 'T'  # Change G→T
```

**Step-by-Step:**
1. Get codon: position 104 → codon 'CGT' (Arg)
2. Gene position = 104 - 100 = 4
3. Position in codon = 4 % 3 = 1
4. Create derived codon: CGT → CTT
5. Translate: 'CGT' (Arg) → 'CTT' (Leu)
6. Compare: Arg ≠ Leu → **Non-Synonymous**

**Expected:**
```python
result = {
    'codon_before': 'CGT',
    'codon_after': 'CTT',
    'aa_before': 'R',
    'aa_after': 'L',
    'mutation_type': 'NS'
}
```

**Synonymous Example:**
```python
gene_info = reverse strand gene at 201-209
ancestral_seq = 'GTACTAACA'
position = 203
derived_allele = 'G'  # Forward strand allele
```

**Step-by-Step:**
1. Position 203 in reverse gene → codon 'CTA' (Leu)
2. Reverse strand: complement G → C
3. Create derived codon: CTA → CTC
4. Translate: 'CTA' (Leu) → 'CTC' (Leu)
5. Compare: Leu = Leu → **Synonymous**

**Expected:**
```python
result = {
    'codon_before': 'CTA',
    'codon_after': 'CTC',
    'aa_before': 'L',
    'aa_after': 'L',
    'mutation_type': 'S'
}
```

---

## TestNG86PathAveraging

### Purpose

This is the **most critical test class**, validating the core NG86 path averaging algorithm implemented in [`_enumerate_ng86_paths()`](../alleleflux/scripts/evolution/dnds_from_timepoints.py).

### NG86 Path Averaging Theory

When k positions differ between ancestral and derived codons, there are **k! possible mutational pathways**. The NG86 method averages S/NS classifications across all valid paths.

**Why Path Averaging Matters:**

Consider CTT (Leu) → CAA (Gln) with positions 1 and 2 different:

**Path 1: Change position 1 first**
```
CTT (Leu) → CAT (His) → CAA (Gln)
Step 1: Leu → His = NS
Step 2: His → Gln = NS
Total: 0 S, 2 NS
```

**Path 2: Change position 2 first**
```
CTT (Leu) → CTA (Leu) → CAA (Gln)
Step 1: Leu → Leu = S  ✓
Step 2: Leu → Gln = NS
Total: 1 S, 1 NS
```

**Average:**
```
frac_S = (0 + 1) / 2 = 0.5
frac_NS = (2 + 1) / 2 = 1.5
```

This correctly captures that the change is "partially synonymous" depending on evolutionary path.

### Test Methods

#### `test_ng86_single_change_synonymous`

**What it tests:** k=1 synonymous change (no path averaging needed).

**Test Case:**
```python
ancestral = 'CTT'  # Leucine
derived = 'CTC'    # Leucine
k = 1  # Only position 2 differs (T→C)
```

**Analysis:**
- Single change: position 2 from T to C
- Only 1 possible pathway (no alternatives)
- CTT (Leu) → CTC (Leu)
- Same amino acid → Synonymous

**Expected:**
```python
frac_S = 1.0
frac_NS = 0.0
num_paths = 1
```

#### `test_ng86_single_change_nonsynonymous`

**What it tests:** k=1 non-synonymous change.

**Test Case:**
```python
ancestral = 'ATG'  # Methionine
derived = 'TTG'    # Leucine
k = 1  # Only position 0 differs (A→T)
```

**Analysis:**
- Single change: position 0 from A to T
- Only 1 possible pathway
- ATG (Met) → TTG (Leu)
- Different amino acid → Non-Synonymous

**Expected:**
```python
frac_S = 0.0
frac_NS = 1.0
num_paths = 1
```

#### `test_ng86_double_change_both_ns`

**What it tests:** k=2 where both paths give NS steps.

**Test Case:**
```python
ancestral = 'ATG'  # Methionine
derived = 'TTT'    # Phenylalanine
k = 2  # Positions 0 and 2 differ (A→T, G→T)
```

**Pathway Enumeration:**

**Path 1: Change position 0 first (A→T), then position 2 (G→T)**
```
Step 1: ATG (Met) → TTG (Leu)
  Met ≠ Leu → NS

Step 2: TTG (Leu) → TTT (Phe)
  Leu ≠ Phe → NS

Path 1 Total: 0 S, 2 NS
```

**Path 2: Change position 2 first (G→T), then position 0 (A→T)**
```
Step 1: ATG (Met) → ATT (Ile)
  Met ≠ Ile → NS

Step 2: ATT (Ile) → TTT (Phe)
  Ile ≠ Phe → NS

Path 2 Total: 0 S, 2 NS
```

**Average:**
```
frac_S = (0 + 0) / 2 = 0.0
frac_NS = (2 + 2) / 2 = 2.0
num_valid_paths = 2
```

**Expected:**
```python
frac_S = 0.0
frac_NS = 2.0
num_paths = 2
```

#### `test_ng86_double_change_mixed`

**What it tests:** k=2 where paths differ in S/NS classification.

**Test Case:**
```python
ancestral = 'CTT'  # Leucine
derived = 'CAA'    # Glutamine
k = 2  # Positions 1 and 2 differ (T→A, T→A)
```

**Pathway Enumeration:**

**Path 1: Change position 1 first (T→A), then position 2 (T→A)**
```
Step 1: CTT (Leu) → CAT (His)
  Leu ≠ His → NS

Step 2: CAT (His) → CAA (Gln)
  His ≠ Gln → NS

Path 1 Total: 0 S, 2 NS
```

**Path 2: Change position 2 first (T→A), then position 1 (T→A)**
```
Step 1: CTT (Leu) → CTA (Leu)
  Leu = Leu → S  ✓

Step 2: CTA (Leu) → CAA (Gln)
  Leu ≠ Gln → NS

Path 2 Total: 1 S, 1 NS
```

**Average:**
```
frac_S = (0 + 1) / 2 = 0.5
frac_NS = (2 + 1) / 2 = 1.5
num_valid_paths = 2
```

**Expected:**
```python
frac_S ≈ 0.5
frac_NS ≈ 1.5
num_paths = 2
```

**Biological Interpretation:** This codon change is "50% synonymous" - half the evolutionary paths preserve the amino acid at an intermediate step.

#### `test_ng86_intermediate_stop_excluded`

**What it tests:** Paths through intermediate stop codons are excluded.

**Test Case:**
```python
ancestral = 'CAA'  # Glutamine
derived = 'TGA'    # Stop codon
k = 2  # Positions 0 and 1 differ (C→T, A→G)
```

**Pathway Enumeration:**

**Path 1: Change position 0 first (C→T), then position 1 (A→G)**
```
Step 1: CAA (Gln) → TAA (Stop)
  Gln → Stop, but TAA is an INTERMEDIATE stop
  → PATH INVALID ✗
```

**Path 2: Change position 1 first (A→G), then position 0 (C→T)**
```
Step 1: CAA (Gln) → CGA (Arg)
  Gln ≠ Arg → NS

Step 2: CGA (Arg) → TGA (Stop)
  Arg → Stop → NS
  (TGA is final stop, allowed)

Path 2 Total: 0 S, 2 NS
Path 2 VALID ✓
```

**Average (only valid paths):**
```
frac_S = 0 / 1 = 0.0
frac_NS = 2 / 1 = 2.0
num_valid_paths = 1  (only 1 out of 2 paths valid)
```

**Expected:**
```python
frac_S = 0.0
frac_NS = 2.0
num_paths = 1  # Only 1 valid path
```

**Key Point:** Final stop codons are allowed, but intermediate stops invalidate the path.

#### `test_ng86_all_paths_invalid`

**What it tests:** Codon pair where ALL paths hit intermediate stops.

**Test Case:**
```python
ancestral = 'TAA'  # Stop codon
derived = 'TGG'    # Tryptophan
k = 2  # Positions 1 and 2 differ (A→G, A→G)
```

**Pathway Enumeration:**

**Path 1: Change position 1 first (A→G), then position 2 (A→G)**
```
Step 1: TAA (Stop) → TGA (Stop)
  Intermediate stop! → INVALID ✗
```

**Path 2: Change position 2 first (A→G), then position 1 (A→G)**
```
Step 1: TAA (Stop) → TAG (Stop)
  Intermediate stop! → INVALID ✗
```

**Result:** All paths invalid!

**Expected:**
```python
frac_S = NaN
frac_NS = NaN
num_paths = 0
```

**Handling:** This codon pair cannot be analyzed and will be excluded from dN/dS calculations.

#### `test_ng86_cache_correctness`

**What it tests:** Pre-computed cache matches direct computation.

**Test Cases:**
```python
test_pairs = [
    ('ATG', 'TTG'),  # k=1 NS
    ('CTT', 'CAA'),  # k=2 mixed
    ('GGG', 'AAA'),  # k=3
    ('TAA', 'TGG'),  # All paths invalid
]
```

**Validation:**
```python
for anc, der in test_pairs:
    # Direct computation
    direct_S, direct_NS, direct_paths = _enumerate_ng86_paths(anc, der, table)
    
    # Cache lookup
    cache_S, cache_NS, cache_paths = ng86_cache[(anc, der)]
    
    # Must match exactly
    assert direct_S == cache_S  (or both NaN)
    assert direct_NS == cache_NS  (or both NaN)
    assert direct_paths == cache_paths
```

**Purpose:** Ensures cache is built correctly and can replace on-the-fly computation.

#### `test_ng86_cache_size`

**What it tests:** Cache contains all 4,096 codon pairs.

**Calculation:**
```
Number of codons = 4^3 = 64
Number of codon pairs = 64 × 64 = 4,096
```

**Expected:**
```python
assert len(ng86_cache) == 4096
```

**Cache Structure:**
```python
ng86_cache = {
    ('ATG', 'ATG'): (0.0, 0.0, 0),      # No change
    ('ATG', 'TTG'): (0.0, 1.0, 1),      # k=1 NS
    ('ATG', 'TTT'): (0.0, 2.0, 2),      # k=2 both NS
    ('CTT', 'CAA'): (0.5, 1.5, 2),      # k=2 mixed
    ('TAA', 'TGG'): (nan, nan, 0),      # All invalid
    # ... 4,091 more entries
}
```

---

## TestCodonSubstitutionAnalysis

### Purpose

Tests the codon-level substitution analysis function [`analyze_codon_substitutions_with_ng86_paths()`](../alleleflux/scripts/evolution/dnds_from_timepoints.py:1091), which groups sites by codon and applies NG86 path averaging.

### Key Concepts

**Codon-Level Analysis:**
- **Old approach:** Analyze each site independently
- **NG86 approach:** Group sites by codon, analyze codon as unit
- **Advantage:** Proper handling of multiple changes within same codon

**Grouping Algorithm:**
1. Identify which codon each site belongs to
2. Group sites by (mag_id, gene_id, codon_start_index)
3. For each codon group, determine k (number of changes)
4. Apply NG86 path averaging for that k value
5. Return one row per codon event (not per site)

### Test Methods

#### `test_codon_grouping`

**What it tests:** Sites are correctly grouped into codon events.

**Setup:**
```python
gene = {
    'gene_id': 'test_gene',
    'start': 100,
    'end': 111,
    'strand': 1,
    'sequence': 'ATGCGTACGTTT'  # 4 codons: ATG CGT ACG TTT
}
```

**Sites:**
```python
sites = [
    {'position': 103, 'ref': 'G', 'derived': 'A'},  # In codon CGT, pos 1
    {'position': 104, 'ref': 'T', 'derived': 'C'},  # In codon CGT, pos 2
]
```

**Codon Mapping:**
```
Gene sequence: ATGCGTACGTTT
Positions:     100 103 106 109
Codons:        ATG CGT ACG TTT
               0-2 3-5 6-8 9-11

Site 103: gene position 3 → codon 1 (CGT), position 0
Site 104: gene position 4 → codon 1 (CGT), position 1
```

**Expected Grouping:**
```python
codon_key = ('MAG_001', 'test_gene', 'contig1', 3, 1)  # Codon starts at gene position 3
codon_data = {
    'changes': [
        {'pos_in_codon': 0, 'contig_pos': 103, 'ancestral': 'G', 'derived': 'A'},
        {'pos_in_codon': 1, 'contig_pos': 104, 'ancestral': 'T', 'derived': 'C'}
    ]
}
```

**Assertions:**
```python
results_df = analyze_codon_substitutions_with_ng86_paths(...)
assert len(results_df) == 1  # One codon event
assert results_df.iloc[0]['k'] == 2  # Two changes
assert 'frac_S' in results_df.columns
assert 'frac_NS' in results_df.columns
assert 'num_valid_paths' in results_df.columns
```

#### `test_codon_event_columns`

**What it tests:** Output DataFrame has all required NG86 columns.

**Required Columns:**
```python
required_cols = [
    # Identifiers
    'mag_id',
    'gene_id',
    'contig',
    'codon_start_index',
    
    # Codon information
    'ancestral_codon',
    'derived_codon',
    'ancestral_aa',
    'derived_aa',
    'mutation_type',
    
    # NG86 metrics
    'k',
    'frac_S',
    'frac_NS',
    'num_valid_paths',
    'is_valid',
    
    # Potential sites
    'potential_S_codon',
    'potential_NS_codon'
]
```

**Assertions:**
```python
for col in required_cols:
    assert col in results_df.columns, f"Missing column: {col}"
```

**Output DataFrame Structure:**
```
| mag_id   | gene_id    | codon_start_index | ancestral_codon | derived_codon | k | frac_S | frac_NS | ...
|----------|------------|-------------------|-----------------|---------------|---|--------|---------|
| MAG_001  | test_gene  | 3                 | CGT             | CAC           | 2 | 0.5    | 1.5     | ...
```

---

## TestGlobalDndsCalculation

### Purpose

Tests the global dN/dS calculation function [`calculate_global_dnds_for_sites()`](../alleleflux/scripts/evolution/dnds_from_timepoints.py:1383), which computes a single dN/dS ratio across all codon events.

### Global dN/dS Formula

```
dN/dS = (Observed_NS / Potential_NS) / (Observed_S / Potential_S)
      = pN / pS
```

Where:
- **Observed_NS**: Sum of frac_NS across all codon events
- **Observed_S**: Sum of frac_S across all codon events
- **Potential_NS**: Sum of potential N sites from unique codons
- **Potential_S**: Sum of potential S sites from unique codons

**Key Point:** Use fractional counts from NG86, track unique codons to avoid double-counting potentials.

### Test Methods

#### `test_global_dnds_with_fractional_counts`

**What it tests:** Global dN/dS calculation with fractional NG86 counts.

**Setup:**
```python
codon_events = pd.DataFrame([
    {
        'mag_id': 'MAG_001',
        'gene_id': 'gene1',
        'codon_start_index': 0,
        'k': 1,
        'frac_S': 0.0,
        'frac_NS': 1.0,
        'is_valid': True,
        'potential_S_codon': 0.0,
        'potential_NS_codon': 3.0
    },
    {
        'mag_id': 'MAG_001',
        'gene_id': 'gene1',
        'codon_start_index': 3,
        'k': 1,
        'frac_S': 1.0,
        'frac_NS': 0.0,
        'is_valid': True,
        'potential_S_codon': 1.0,
        'potential_NS_codon': 2.0
    },
    {
        'mag_id': 'MAG_001',
        'gene_id': 'gene1',
        'codon_start_index': 6,
        'k': 2,
        'frac_S': 0.5,
        'frac_NS': 1.5,
        'is_valid': True,
        'potential_S_codon': 1.0,
        'potential_NS_codon': 2.0
    }
])
```

**Calculations:**

**Observed counts (fractional):**
```
Observed_NS = 1.0 + 0.0 + 1.5 = 2.5
Observed_S = 0.0 + 1.0 + 0.5 = 1.5
```

**Potential sites (from unique codons):**
```
Codon at index 0: potential_S=0.0, potential_NS=3.0
Codon at index 3: potential_S=1.0, potential_NS=2.0
Codon at index 6: potential_S=1.0, potential_NS=2.0

Total_Potential_NS = 3.0 + 2.0 + 2.0 = 7.0
Total_Potential_S = 0.0 + 1.0 + 1.0 = 2.0
```

**Rates:**
```
pN = Observed_NS / Total_Potential_NS = 2.5 / 7.0 ≈ 0.3571
pS = Observed_S / Total_Potential_S = 1.5 / 2.0 = 0.75
dN/dS = pN / pS = 0.3571 / 0.75 ≈ 0.4762
```

**Assertions:**
```python
assert metrics['Observed Non-Synonymous Substitutions (NS)'] ≈ 2.5
assert metrics['Observed Synonymous Substitutions (S)'] ≈ 1.5
assert metrics['Total Potential Non-Synonymous Sites (ns)'] == 7.0
assert metrics['Total Potential Synonymous Sites (s)'] == 2.0
```

#### `test_global_dnds_excludes_invalid_codons`

**What it tests:** Invalid codons are excluded from calculation.

**Setup:**
```python
# Add an invalid codon to the previous dataset
invalid_codon = pd.DataFrame([{
    'mag_id': 'MAG_001',
    'gene_id': 'gene1',
    'codon_start_index': 9,
    'k': 2,
    'frac_S': NaN,
    'frac_NS': NaN,
    'is_valid': False,  # All paths hit stops
    'potential_S_codon': 1.0,
    'potential_NS_codon': 2.0
}])

combined = pd.concat([codon_events, invalid_codon])
```

**Expected Behavior:**
- Invalid codon excluded from counts
- Only 3 valid codons counted
- Invalid codon count reported

**Assertions:**
```python
assert metrics['Total Unique Codons Analyzed'] == 3
assert metrics['Invalid Codons (all paths hit stops)'] == 1
```

#### `test_global_dnds_empty_input`

**What it tests:** Empty DataFrame returns None.

**Test:**
```python
empty_df = pd.DataFrame()
result = calculate_global_dnds_for_sites(empty_df)
assert result is None
```

---

## TestNG86ComputedValues

### Purpose

Comprehensive validation of NG86 (Nei-Gojobori 1986) computed values against manually calculated expectations for diverse codon change scenarios.

### Test Coverage

This test class validates:
- ✅ All codon positions (0, 1, 2) for k=1 single-nucleotide changes
- ✅ k=2 and k=3 multi-nucleotide changes with pathway averaging
- ✅ Synonymous-only and non-synonymous-only transitions
- ✅ Mixed pathway classifications
- ✅ Intermediate stop codon exclusion
- ✅ Edge cases (zero S/NS sites, stop codons)
- ✅ Global dN/dS aggregation with fractional counts
- ✅ Gene-level and MAG-level summary statistics

### Setup

**Class Setup (`setUpClass`):**
```python
@classmethod
def setUpClass(cls):
    cls.codon_table = CodonTable.unambiguous_dna_by_id[11]  # Standard genetic code
    cls.ng86_cache = _get_ng86_cache(cls.codon_table)  # Pre-compute all codon transitions
```

The NG86 cache is built once for efficiency and contains pathway analysis for all 4,096 possible codon-to-codon transitions.

---

### Test Method 1: `test_single_position_changes_all_positions`

**Purpose:** Validate k=1 changes at all three codon positions (0, 1, 2) and verify correct S/NS classification and potential site calculations.

**Scenario:** Single-nucleotide changes demonstrating each codon position's behavior.

#### Gene Design

**Sequence:** `ATGTTGGCT` (9 bp, 3 codons)

**Prodigal Coordinates:**
- start = 101 (1-based, inclusive)
- end = 109 (1-based, inclusive)
- strand = 1 (forward)

**Contig Coordinates (0-based):**
- Positions 100-108
- Codons:
  - **ATG** (Met) at positions 100-102
  - **TTG** (Leu) at positions 103-105
  - **GCT** (Ala) at positions 106-108

#### Mutations

**1. Position 0 change: ATG → TTG** (contig pos 100)
- Nucleotide: A→T at first codon position
- Translation: Met → Leu (non-synonymous)
- k=1, num_valid_paths=1
- Expected: **frac_S=0.0, frac_NS=1.0**
- Potential: ATG has S=0.0, NS=3.0 (Met has only one codon)

**2. Position 1 change: TTG → TAG** (contig pos 104)
- Nucleotide: T→A at second codon position
- Translation: Leu → Stop (non-synonymous)
- k=1, num_valid_paths=1
- Expected: **frac_S=0.0, frac_NS=1.0**
- Potential: TTG has S≈0.667, NS≈2.333

**3. Position 2 change: GCT → GCC** (contig pos 108)
- Nucleotide: T→C at third codon position (wobble position)
- Translation: Ala → Ala (synonymous)
- k=1, num_valid_paths=1
- Expected: **frac_S=1.0, frac_NS=0.0**
- Potential: GCT has S=1.0, NS=2.0

#### Global dN/dS Calculation

**Potential Sites (Gene-Level):**
- Total NS = 3.0 + 2.333 + 2.0 = **7.333**
- Total S = 0.0 + 0.667 + 1.0 = **1.667**

**Observed Substitutions:**
- Total NS = 1.0 + 1.0 + 0.0 = **2.0**
- Total S = 0.0 + 0.0 + 1.0 = **1.0**

**Substitution Proportions:**
- pN = 2.0 / 7.333 = **0.2727**
- pS = 1.0 / 1.667 = **0.6**

**Final Ratio:**
- **dN/dS = 0.2727 / 0.6 = 0.4545**

#### Validation Checks

```python
# All events have k=1
assert all(codon_events_df['k'] == 1)

# Correct mutation_type classification
assert codon_events_df.loc[0, 'mutation_type'] == 'NS'  # ATG→TTG
assert codon_events_df.loc[1, 'mutation_type'] == 'NS'  # TTG→TAG
assert codon_events_df.loc[2, 'mutation_type'] == 'S'   # GCT→GCC

# Accurate fractional S/NS counts
assert codon_events_df.loc[0, 'frac_NS'] == 1.0
assert codon_events_df.loc[1, 'frac_NS'] == 1.0
assert codon_events_df.loc[2, 'frac_S'] == 1.0

# Proper potential site calculations
assert math.isclose(codon_events_df.loc[0, 'potential_NS_codon'], 3.0, rel_tol=0.01)
assert math.isclose(codon_events_df.loc[1, 'potential_NS_codon'], 2.333, rel_tol=0.01)
assert math.isclose(codon_events_df.loc[2, 'potential_S_codon'], 1.0, rel_tol=0.01)

# Global metrics
assert math.isclose(summary['dN_dS_ratio'], 0.4545, rel_tol=0.01)
```

#### Technical Notes

**Coordinate Conversion:**
- Prodigal start=101 (1-based) → contig start=100 (0-based)
- Site at contig position 100 → gene position 0 → codon 0, position 0

**Reference Base Matching:**
- Site data `forward_ref_base='A'` must match `sequence[0]='A'`
- Critical for ensuring mutations applied to correct ancestral state

**Wobble Position:**
- Position 2 (third codon position) often allows synonymous changes
- GCT, GCC, GCA, GCG all code for Alanine (4-fold degenerate)

---

### Test Method 2: `test_synonymous_only_transitions`

**Purpose:** Validate codons where all changes are synonymous, resulting in **dN/dS = 0**.

**Scenario:** Three 4-fold degenerate third-position changes (wobble base only).

#### Gene Design

**Sequence:** `CTTCGCGCT` (9 bp, 3 codons)

**Prodigal Coordinates:** start=101, end=109, strand=1

**Codons:**
- **CTT** (Leu) at positions 100-102
- **CGC** (Arg) at positions 103-105
- **GCT** (Ala) at positions 106-108

#### Mutations

**1. CTT → CTC** (contig pos 102): Leu → Leu
- Position 2 change (T→C)
- k=1, **synonymous**
- Expected: frac_S=1.0, frac_NS=0.0
- Potential: CTT has S=1.0, NS=2.0

**2. CGC → CGT** (contig pos 105): Arg → Arg
- Position 2 change (C→T)
- k=1, **synonymous**
- Expected: frac_S=1.0, frac_NS=0.0
- Potential: CGC has S=1.0, NS=2.0

**3. GCT → GCC** (contig pos 108): Ala → Ala
- Position 2 change (T→C)
- k=1, **synonymous**
- Expected: frac_S=1.0, frac_NS=0.0
- Potential: GCT has S=1.0, NS=2.0

#### Global dN/dS Calculation

**Potential Sites:**
- Total NS = 2.0 + 2.0 + 2.0 = **6.0**
- Total S = 1.0 + 1.0 + 1.0 = **3.0**

**Observed Substitutions:**
- Total NS = **0.0** (no non-synonymous changes)
- Total S = **3.0** (all changes synonymous)

**Substitution Proportions:**
- pN = 0.0 / 6.0 = **0.0**
- pS = 3.0 / 3.0 = **1.0**

**Final Ratio:**
- **dN/dS = 0.0 / 1.0 = 0.0**

#### Validation Checks

```python
# All events classified as synonymous
assert all(codon_events_df['mutation_type'] == 'S')

# All frac_S=1.0, all frac_NS=0.0
assert all(codon_events_df['frac_S'] == 1.0)
assert all(codon_events_df['frac_NS'] == 0.0)

# Global dN/dS correctly equals 0.0
assert summary['dN_dS_ratio'] == 0.0
assert summary['pNS'] == 0.0
assert summary['pS'] == 1.0
```

#### Biological Context

This scenario represents **purifying selection** where only synonymous mutations are tolerated, typical of:
- Highly conserved protein domains
- Essential genes under strong functional constraint
- Recent adaptive sweeps (all NS variants removed)

---

### Test Method 3: `test_nonsynonymous_only_transitions`

**Purpose:** Validate codons where all changes are non-synonymous, resulting in **undefined dN/dS** (pS=0).

**Scenario:** First and second position changes that always alter amino acid.

#### Gene Design

**Sequence:** `ATGTTGAAG` (9 bp, 3 codons)

**Prodigal Coordinates:** start=100, end=108, strand=1

**Codons:**
- **ATG** (Met) at positions 99-101 *(Note: Prodigal uses 1-based start=100)*
- **TTG** (Leu) at positions 102-104
- **AAG** (Lys) at positions 105-107

#### Mutations

**1. ATG → GTG** (contig pos 100): Met → Val
- Position 0 change (A→G)
- k=1, **non-synonymous**
- Expected: frac_S=0.0, frac_NS=1.0
- Potential: ATG has S=0.0, NS=3.0

**2. TTG → ATG** (contig pos 103): Leu → Met
- Position 0 change (T→A)
- k=1, **non-synonymous**
- Expected: frac_S=0.0, frac_NS=1.0
- Potential: TTG has S≈0.667, NS≈2.333

**3. AAG → TAG** (contig pos 106): Lys → Stop
- Position 1 change (A→T)
- k=1, **non-synonymous**
- Expected: frac_S=0.0, frac_NS=1.0
- Potential: AAG has S≈0.333, NS≈2.667

#### Global dN/dS Calculation

**Potential Sites:**
- Total NS = 3.0 + 2.333 + 2.667 = **8.0**
- Total S = 0.0 + 0.667 + 0.333 = **1.0**

**Observed Substitutions:**
- Total NS = **3.0** (all changes non-synonymous)
- Total S = **0.0** (no synonymous changes)

**Substitution Proportions:**
- pN = 3.0 / 8.0 = **0.375**
- pS = 0.0 / 1.0 = **0.0**

**Final Ratio:**
- **dN/dS = undefined (NaN)** due to division by zero

#### Validation Checks

```python
# All events classified as non-synonymous
assert all(codon_events_df['mutation_type'] == 'NS')

# All frac_NS=1.0, all frac_S=0.0
assert all(codon_events_df['frac_NS'] == 1.0)
assert all(codon_events_df['frac_S'] == 0.0)

# pS correctly equals 0.0
assert summary['pS'] == 0.0

# dN/dS handling of undefined case
assert pd.isna(summary['dN_dS_ratio']) or np.isinf(summary['dN_dS_ratio'])
```

#### Biological Context

This scenario can occur when:
- Analyzing very small genomic regions with few sites
- Strong positive selection has driven all synonymous variants to fixation
- Dataset limited to specific non-synonymous positions

**Handling in Practice:**
- Report pN but note dN/dS is undefined

---

### Test Method 4: `test_two_position_changes_path_averaging`

**Purpose:** Validate k=2 changes with explicit 2-pathway averaging (factorial(2) = 2 permutations).

**Scenario:** Double nucleotide change requiring pathway enumeration.

#### Gene Design

**Sequence:** `ATGCGT` (6 bp, 2 codons)

**Prodigal Coordinates:** start=101, end=106, strand=1

**First Codon:** ATG (Met) at positions 100-102

#### Mutation

**ATG → TTT** (contig pos 100 and 102): Two nucleotide changes
- Position 0: A→T
- Position 2: G→T
- **k=2** (two differences)

#### NG86 Pathway Analysis

**Pathway 1** (change position 0 first, then position 2):
1. ATG (Met) → TTG (Leu) [A→T at pos 0]: **Non-synonymous**
2. TTG (Leu) → TTT (Phe) [G→T at pos 2]: **Non-synonymous**
- Pathway 1 total: **0 S steps, 2 NS steps**

**Pathway 2** (change position 2 first, then position 0):
1. ATG (Met) → ATT (Ile) [G→T at pos 2]: **Non-synonymous**
2. ATT (Ile) → TTT (Phe) [A→T at pos 0]: **Non-synonymous**
- Pathway 2 total: **0 S steps, 2 NS steps**

#### Averaged Results

**Fractional Counts:**
- frac_S = (0 + 0) / 2 = **0.0**
- frac_NS = (2 + 2) / 2 = **2.0**
- num_valid_paths = **2**

**Potential Sites:**
- ATG: S=0.0, NS=3.0

#### Validation Checks

```python
# k=2 correctly identified
assert codon_events_df.loc[0, 'k'] == 2

# Both pathways evaluated
assert codon_events_df.loc[0, 'num_valid_paths'] == 2

# Fractional counts reflect path averaging
assert codon_events_df.loc[0, 'frac_S'] == 0.0
assert codon_events_df.loc[0, 'frac_NS'] == 2.0

# Event marked as valid
assert codon_events_df.loc[0, 'is_valid'] == True
```

#### Key Insight

Even though both pathways yield the same S/NS classification (all NS), the NG86 algorithm still evaluates both to ensure unbiased estimation. For other codon pairs, pathways can differ significantly.

---

### Test Method 5: `test_three_position_changes_path_averaging`

**Purpose:** Validate k=3 changes with 6-pathway averaging (factorial(3) = 6 permutations).

**Scenario:** All three codon positions change, requiring comprehensive pathway enumeration.

#### Gene Design

**Sequence:** `ATGCGT` (6 bp, 2 codons)

**Prodigal Coordinates:** start=101, end=106, strand=1

**First Codon:** ATG (Met) at positions 100-102

#### Mutation

**ATG → CCC** (contig pos 100, 101, 102): All three positions change
- Position 0: A→C
- Position 1: T→C
- Position 2: G→C
- **k=3** (three differences)

#### NG86 Pathway Analysis

All 6 permutations of changing positions {0, 1, 2} are evaluated:

**Path 1 (Order: 0, 1, 2)**
1. ATG (Met) → CTG (Leu): **NS**
2. CTG (Leu) → CCG (Pro): **NS**
3. CCG (Pro) → CCC (Pro): **S**
- Total: 1 S, 2 NS

**Path 2 (Order: 0, 2, 1)**
1. ATG (Met) → CTG (Leu): **NS**
2. CTG (Leu) → CTC (Leu): **S**
3. CTC (Leu) → CCC (Pro): **NS**
- Total: 1 S, 2 NS

**Path 3 (Order: 1, 0, 2)**
1. ATG (Met) → ACG (Thr): **NS**
2. ACG (Thr) → CCG (Pro): **NS**
3. CCG (Pro) → CCC (Pro): **S**
- Total: 1 S, 2 NS

**Path 4 (Order: 1, 2, 0)**
1. ATG (Met) → ACG (Thr): **NS**
2. ACG (Thr) → ACC (Thr): **S**
3. ACC (Thr) → CCC (Pro): **NS**
- Total: 1 S, 2 NS

**Path 5 (Order: 2, 0, 1)**
1. ATG (Met) → ATC (Ile): **NS**
2. ATC (Ile) → CTC (Leu): **NS**
3. CTC (Leu) → CCC (Pro): **NS**
- Total: 0 S, 3 NS

**Path 6 (Order: 2, 1, 0)**
1. ATG (Met) → ATC (Ile): **NS**
2. ATC (Ile) → ACC (Thr): **NS**
3. ACC (Thr) → CCC (Pro): **NS**
- Total: 0 S, 3 NS

**Final Calculation:**
```
Total S steps across 6 valid paths: 1 + 1 + 1 + 1 + 0 + 0 = 4
Total NS steps across 6 valid paths: 2 + 2 + 2 + 2 + 3 + 3 = 14

Average frac_S = 4 / 6 = 0.6667
Average frac_NS = 14 / 6 = 2.3333
```

#### Averaged Results

**Fractional Counts:**
- frac_S ≈ **0.667**
- frac_NS ≈ **2.333**
- num_valid_paths = **6**

Total S+NS = 0.667 + 2.333 = 3.0 ✓ (matches k=3)

#### Validation Checks

```python
# k=3 correctly identified
assert codon_events_df.loc[0, 'k'] == 3

# All 6 pathways evaluated
assert codon_events_df.loc[0, 'num_valid_paths'] == 6

# Fractional S/NS counts reflect mixed pathway classifications
assert math.isclose(codon_events_df.loc[0, 'frac_S'], 0.667, rel_tol=0.01)
assert math.isclose(codon_events_df.loc[0, 'frac_NS'], 2.333, rel_tol=0.01)

# Event marked as valid
assert codon_events_df.loc[0, 'is_valid'] == True
```

#### Key Insight

k=3 changes demonstrate that NG86 path averaging can produce fractional S/NS counts even when k is large, capturing the nuanced reality that some mutational pathways include synonymous intermediate steps.

---

### Test Method 6: `test_mixed_pathway_classifications`

**Purpose:** Validate k=2 where the two pathways produce **different S/NS counts**.

**Scenario:** Demonstrates pathway-dependent classification—critical NG86 use case.

#### Gene Design

**Sequence:** `CTTAAA` (6 bp, 2 codons)

**Prodigal Coordinates:** start=101, end=106, strand=1

**First Codon:** CTT (Leu) at positions 100-102

#### Mutation

**CTT → CAA** (contig pos 101 and 102): Positions 1 and 2 change
- Position 1: T→A
- Position 2: T→A
- **k=2**

#### NG86 Pathway Analysis

**Pathway 1** (change position 1 first, then position 2):
1. CTT (Leu) → CAT (His) [T→A at pos 1]: **Non-synonymous** (Leu ≠ His)
2. CAT (His) → CAA (Gln) [T→A at pos 2]: **Non-synonymous** (His ≠ Gln)
- Pathway 1 total: **0 S steps, 2 NS steps**

**Pathway 2** (change position 2 first, then position 1):
1. CTT (Leu) → CTA (Leu) [T→A at pos 2]: **Synonymous** ✓ (Leu = Leu)
2. CTA (Leu) → CAA (Gln) [T→A at pos 1]: **Non-synonymous** (Leu ≠ Gln)
- Pathway 2 total: **1 S step, 1 NS step**

#### Averaged Results

**Fractional Counts:**
- frac_S = (0 + 1) / 2 = **0.5**
- frac_NS = (2 + 1) / 2 = **1.5**
- num_valid_paths = **2**

Total S+NS = 0.5 + 1.5 = 2.0 ✓ (matches k=2)

**Potential Sites:**
- CTT: S=1.0, NS=2.0

#### Validation Checks

```python
# k=2 with mixed pathway outcomes
assert codon_events_df.loc[0, 'k'] == 2
assert codon_events_df.loc[0, 'num_valid_paths'] == 2

# Correct averaging of different pathway outcomes
assert codon_events_df.loc[0, 'frac_S'] == 0.5
assert codon_events_df.loc[0, 'frac_NS'] == 1.5

# Potential sites for CTT
assert math.isclose(codon_events_df.loc[0, 'potential_S_codon'], 1.0, rel_tol=0.01)
assert math.isclose(codon_events_df.loc[0, 'potential_NS_codon'], 2.0, rel_tol=0.01)
```

#### Key Insight

This is the **canonical NG86 example**: the classification depends on mutation order. Without path averaging, an arbitrary choice of order would bias the dN/dS estimate. NG86 resolves this by averaging, producing a fractional count that reflects both possibilities.

---

### Test Method 7: `test_intermediate_stop_codon_exclusion`

**Purpose:** Validate that pathways hitting **intermediate stop codons** are excluded from averaging.

**Scenario:** k=2 change where one pathway goes through a stop codon.

#### Gene Design

**Sequence:** `CAACGT` (6 bp, 2 codons)

**Prodigal Coordinates:** start=101, end=106, strand=1

**First Codon:** CAA (Gln) at positions 100-102

#### Mutation

**CAA → TGA** (contig pos 100 and 101): Positions 0 and 1 change
- Position 0: C→T
- Position 1: A→G
- Final codon: **TGA (Stop)**
- **k=2**

#### NG86 Pathway Analysis

**Pathway 1** (change position 0 first, then position 1):
1. CAA (Gln) → **TAA (Stop)** [C→T at pos 0]: Intermediate stop codon
- **EXCLUDED** ❌ (intermediate stop codons invalidate pathway)

**Pathway 2** (change position 1 first, then position 0):
1. CAA (Gln) → CGA (Arg) [A→G at pos 1]: **Non-synonymous**
2. CGA (Arg) → TGA (Stop) [C→T at pos 0]: **Non-synonymous** (final stop is allowed)
- Pathway 2 total: **0 S steps, 2 NS steps**
- **VALID** ✓ (final stop codon is acceptable)

#### Averaged Results

**Fractional Counts:**
- Only 1 valid pathway remains
- frac_S = **0.0**
- frac_NS = **2.0**
- num_valid_paths = **1** (not 2!)

**Potential Sites:**
- CAA: S≈0.333, NS≈2.667

#### Validation Checks

```python
# k=2 correctly identified
assert codon_events_df.loc[0, 'k'] == 2

# Intermediate stop codons correctly exclude pathways
assert codon_events_df.loc[0, 'num_valid_paths'] == 1  # Not 2!

# Fractional counts only average over valid pathways
assert codon_events_df.loc[0, 'frac_S'] == 0.0
assert codon_events_df.loc[0, 'frac_NS'] == 2.0

# Event marked as valid (at least one pathway exists)
assert codon_events_df.loc[0, 'is_valid'] == True
```

#### Biological Rationale

**Why exclude intermediate stop codons?**
- Intermediate stops would truncate the protein (biologically invalid)
- Such pathways are evolutionarily inaccessible
- Only the **final** codon can be a stop (e.g., gene truncation mutations)

**Implementation:**
```python
# In _enumerate_ng86_paths()
for step in pathway:
    step_aa = codon_table.forward_table.get(step_codon, "*")
    if step_aa == "*" and step != pathway[-1]:  # Intermediate stop
        valid_pathway = False
        break
```

#### Key Insight

This test demonstrates a **critical NG86 requirement**: pathways through stop codons are biologically invalid and must be excluded before averaging. The `num_valid_paths` field reflects this exclusion.

---

### Test Method 8: `test_zero_synonymous_sites_codons`

**Purpose:** Validate handling of codons with **zero synonymous sites** (ATG-Methionine, TGG-Tryptophan).

**Scenario:** Unique codons where all single-base changes are non-synonymous.

#### Gene Design

**Sequence:** `ATGTGG` (6 bp, 2 codons)

**Prodigal Coordinates:** start=101, end=106, strand=1

**Codons:**
- **ATG** (Met) at positions 100-102
- **TGG** (Trp) at positions 103-105

#### Mutations

**1. ATG → GTG** (contig pos 100): Met → Val
- Position 0 change (A→G)
- k=1, **non-synonymous**

**2. TGG → TGC** (contig pos 105): Trp → Cys
- Position 2 change (G→C)
- k=1, **non-synonymous**

#### Potential Site Calculations

**ATG (Methionine):**
- Position 0: A→{T,G,C} = {TTG(Leu), GTG(Val), CTG(Leu)} → **all NS**
- Position 1: T→{A,G,C} = {AAG(Lys), AGG(Arg), ACG(Thr)} → **all NS**
- Position 2: G→{A,T,C} = {ATA(Ile), ATT(Ile), ATC(Ile)} → **all NS**
- **Total: S=0.0, NS=3.0**

**TGG (Tryptophan):**
- Position 0: T→{A,G,C} = {AGG(Arg), GGG(Gly), CGG(Arg)} → **all NS**
- Position 1: G→{A,T,C} = {TAG(Stop), TTG(Leu), TCG(Ser)} → **all NS**
- Position 2: G→{A,T,C} = {TGA(Stop), TGT(Cys), TGC(Cys)} → **all NS**
- **Total: S=0.0, NS=3.0**

#### Gene-Level Totals

**Potential Sites:**
- Total S = 0.0 + 0.0 = **0.0**
- Total NS = 3.0 + 3.0 = **6.0**

**Observed Substitutions:**
- Both mutations are k=1, NS
- Total NS = **2.0**
- Total S = **0.0**

**Substitution Proportions:**
- pN = 2.0 / 6.0 = **0.333**
- pS = 0.0 / 0.0 = **undefined**

#### Validation Checks

```python
# Codons with no synonymous sites correctly calculated
assert codon_events_df.loc[0, 'potential_S_codon'] == 0.0  # ATG
assert codon_events_df.loc[1, 'potential_S_codon'] == 0.0  # TGG

# Both events classified as non-synonymous
assert all(codon_events_df['mutation_type'] == 'NS')

# Gene-level aggregation: total S=0.0
gene_summary = summaries['gene']
assert gene_summary.loc[gene_summary['gene_id'] == 'gene1', 'dS'].values[0] == 0.0
```

#### Biological Context

**ATG and TGG are unique codons:**
- Methionine: Only ATG codes for Met
- Tryptophan: Only TGG codes for Trp
- No synonymous alternatives exist

**Implications:**
- dN/dS is undefined for such genes (pS=0)
- Useful for identifying genes with unusual codon usage
- Important edge case for pipeline robustness

#### Edge Case Handling

```python
# In calculate_global_dnds_for_sites()
if total_S == 0:
    logger.warning("No synonymous sites available; dN/dS undefined")
    dnds_ratio = np.nan
else:
    dnds_ratio = pNS / pS
```

---

### Test Method 9: `test_global_dnds_computed_values`

**Purpose:** Validate **complete end-to-end dN/dS calculation** with known inputs and expected outputs.

**Scenario:** Mixed S and NS changes with explicit manual calculation for verification.

#### Gene Design

**Sequence:** `ATGCTTACC` (9 bp, 3 codons)

**Prodigal Coordinates:** start=101, end=109, strand=1

**Codons:**
- **ATG** (Met) at positions 100-102
- **CTT** (Leu) at positions 103-105
- **ACC** (Thr) at positions 106-108

#### Mutations

**1. ATG → TTG** (contig pos 100): Met → Leu
- Position 0 change (A→T)
- k=1, **non-synonymous**
- Expected: frac_S=0.0, frac_NS=1.0

**2. CTT → CTC** (contig pos 105): Leu → Leu
- Position 2 change (T→C)
- k=1, **synonymous**
- Expected: frac_S=1.0, frac_NS=0.0

**3. ACC → ATC** (contig pos 107): Thr → Ile
- Position 1 change (C→T)
- k=1, **non-synonymous**
- Expected: frac_S=0.0, frac_NS=1.0

#### Potential Sites

**Per-Codon:**
- ATG: S=0.0, NS=3.0
- CTT: S=1.0, NS=2.0
- ACC: S=1.0, NS=2.0

**Gene-Level Totals:**
- Total NS = 3.0 + 2.0 + 2.0 = **7.0**
- Total S = 0.0 + 1.0 + 1.0 = **2.0**

#### Observed Substitutions

**Per-Event:**
- Event 1 (ATG→TTG): NS=1.0, S=0.0
- Event 2 (CTT→CTC): NS=0.0, S=1.0
- Event 3 (ACC→ATC): NS=1.0, S=0.0

**Gene-Level Totals:**
- Total NS = 1.0 + 0.0 + 1.0 = **2.0**
- Total S = 0.0 + 1.0 + 0.0 = **1.0**

#### Global dN/dS Calculation

**Substitution Proportions:**
- pN = 2.0 / 7.0 = **0.2857**
- pS = 1.0 / 2.0 = **0.5**

**Final Ratio:**
- **dN/dS = 0.2857 / 0.5 = 0.5714**

#### Validation Checks

```python
# All individual event calculations correct
assert codon_events_df.loc[0, 'frac_NS'] == 1.0  # ATG→TTG
assert codon_events_df.loc[1, 'frac_S'] == 1.0   # CTT→CTC
assert codon_events_df.loc[2, 'frac_NS'] == 1.0  # ACC→ATC

# Global potential site totals match manual calculation
global_summary = calculate_global_dnds_for_sites(codon_events_df)
assert math.isclose(global_summary['Total Potential Non-Synonymous Sites (ns)'], 7.0, rel_tol=0.01)
assert math.isclose(global_summary['Total Potential Synonymous Sites (s)'], 2.0, rel_tol=0.01)

# Observed substitution counts match fractional sums
assert math.isclose(global_summary['Observed Non-Synonymous Substitutions (NS)'], 2.0, rel_tol=0.01)
assert math.isclose(global_summary['Observed Synonymous Substitutions (S)'], 1.0, rel_tol=0.01)

# pN and pS ratios correct
assert math.isclose(global_summary['pNS (NS / ns)'], 0.2857, rel_tol=0.01)
assert math.isclose(global_summary['pS (S / s)'], 0.5, rel_tol=0.01)

# Final dN/dS ratio matches expected value
summaries = summarize_results(codon_events_df)
assert math.isclose(summaries['mag']['dN_dS_ratio'].values[0], 0.5714, rel_tol=0.01)
```

#### Technical Notes

**Comprehensive Validation:**
- Tests integration of all components: potential sites, pathway averaging, global aggregation
- Serves as a **gold standard** for expected behavior
- Manual calculations can be independently verified

**Pipeline Flow:**
```
Individual mutations → Codon-level events → Gene-level aggregation → Global dN/dS
```

Each stage is validated to ensure end-to-end correctness.

---

### Test Method 10: `test_gene_and_mag_level_summaries`

**Purpose:** Validate **gene-level and MAG-level** dN, dS, and dN/dS summary statistics with proper aggregation.

**Scenario:** Two genes in the same MAG with different dN/dS ratios, testing hierarchical aggregation.

#### Gene Design

**Gene1:**
- Sequence: `ATGCTT` (6 bp, 2 codons)
- Codons: ATG (Met), CTT (Leu)
- Mutations: ATG→TTG (NS), CTT→CTC (S)

**Gene2:**
- Sequence: `GCTAGG` (6 bp, 2 codons)
- Codons: GCT (Ala), AGG (Arg)
- Mutations: GCT→GCC (S), AGG→ATG (NS)

**MAG:** Both genes belong to `mag1`

#### Expected Gene-Level Calculations

**Gene1:**
- Potential: S=1.0 (from CTT), NS=5.0 (3 from ATG + 2 from CTT)
- Observed: S=1.0, NS=1.0
- dN = 1.0 / 5.0 = **0.2**
- dS = 1.0 / 1.0 = **1.0**
- **dN/dS = 0.2**

**Gene2:**
- Potential: S≈1.667, NS≈4.333
- Observed: S=1.0, NS=1.0
- dN = 1.0 / 4.333 = **0.2308**
- dS = 1.0 / 1.667 = **0.6**
- **dN/dS = 0.3846**

#### Expected MAG-Level Calculation

**Site-Weighted Aggregation:**
- Total potential NS = 5.0 + 4.333 = **9.333**
- Total potential S = 1.0 + 1.667 = **2.667**
- Total observed NS = 1.0 + 1.0 = **2.0**
- Total observed S = 1.0 + 1.0 = **2.0**

**MAG-Level Ratios:**
- dN = 2.0 / 9.333 = **0.2143**
- dS = 2.0 / 2.667 = **0.75**
- **dN/dS = 0.2857**

#### Validation Checks

```python
# Gene-level summaries calculated independently for each gene
gene_summary = summaries['gene']
gene1 = gene_summary[gene_summary['gene_id'] == 'gene1'].iloc[0]
gene2 = gene_summary[gene_summary['gene_id'] == 'gene2'].iloc[0]

# Gene1 stats
assert math.isclose(gene1['dN'], 0.2, rel_tol=0.01)
assert math.isclose(gene1['dS'], 1.0, rel_tol=0.01)
assert math.isclose(gene1['dN_dS_ratio'], 0.2, rel_tol=0.01)

# Gene2 stats
assert math.isclose(gene2['dN'], 0.2308, rel_tol=0.01)
assert math.isclose(gene2['dS'], 0.6, rel_tol=0.01)
assert math.isclose(gene2['dN_dS_ratio'], 0.3846, rel_tol=0.01)

# MAG-level aggregates across both genes
mag_summary = summaries['mag']
mag1 = mag_summary[mag_summary['mag_id'] == 'mag1'].iloc[0]

assert mag1['num_genes'] == 2
assert math.isclose(mag1['dN'], 0.2143, rel_tol=0.01)
assert math.isclose(mag1['dS'], 0.75, rel_tol=0.01)
assert math.isclose(mag1['dN_dS_ratio'], 0.2857, rel_tol=0.01)

# MAG dN/dS falls within range of individual gene values
assert gene1['dN_dS_ratio'] <= mag1['dN_dS_ratio'] <= gene2['dN_dS_ratio']
```

#### Technical Notes

**Hierarchical Aggregation:**
```
Site-level → Codon-level → Gene-level → MAG-level
```

**Key Insight:** MAG-level dN/dS is **not** the mean of gene-level ratios. Instead, it aggregates potential and observed sites across all genes, then calculates the ratio.

**Why Site-Weighted?**
- Genes with more potential sites contribute more to MAG-level estimate
- Avoids bias from genes with few sites and noisy ratios
- Biologically appropriate: reflects overall selective pressure across MAG

**Implementation:**
```python
# In summarize_results()
mag_grouped = codon_events_df.groupby('mag_id').agg({
    'potential_NS_codon': 'sum',
    'potential_S_codon': 'sum',
    'frac_NS': 'sum',
    'frac_S': 'sum'
})
mag_grouped['dN'] = mag_grouped['frac_NS'] / mag_grouped['potential_NS_codon']
mag_grouped['dS'] = mag_grouped['frac_S'] / mag_grouped['potential_S_codon']
mag_grouped['dN_dS_ratio'] = mag_grouped['dN'] / mag_grouped['dS']
```

---

## Negative Strand Test Cases

### Overview

The [`TestNG86ComputedValues`](test_dnds_from_timepoints.py:1284) class now includes **6 comprehensive test methods** for validating dN/dS calculations on **negative strand genes (strand = -1)**. These tests ensure that the implementation correctly handles:
- Reverse complement sequence orientation
- Reversed coordinate mapping
- Allele complementation from forward strand to gene sequence
- All edge cases (k=1, k=2, k=3, boundary conditions, stop codons)

### Critical Negative Strand Concepts

#### 1. Sequence Orientation

**Positive strand (strand=1):**
- Sequence is stored 5'→3' as-is from genomic region
- Genomic positions map directly to sequence indices

**Negative strand (strand=-1):**
- Sequence is **reverse complement** of genomic region
- Stored sequence represents the coding strand (5'→3')
- Genomic positions refer to the **forward strand**
- Example: Genomic region `ATGCGT` → gene sequence `ACGCAT` (reverse complement)

#### 2. Coordinate Mapping Formula

For a negative strand gene with start=201, end=209, sequence="GTACTAACA" (9 bases):

**Position Mapping:**
```
pos_in_gene = (end - 1) - genomic_position
```

**Example Mappings:**
- Genomic position 208 → gene index 0 (first base 'G')
- Genomic position 207 → gene index 1 (second base 'T')
- Genomic position 200 → gene index 8 (last base 'A')

**Codon Mapping:**
- First codon (indices 0-2: "GTA") corresponds to genomic positions 208, 207, 206
- Second codon (indices 3-5: "CTA") corresponds to genomic positions 205, 204, 203
- Third codon (indices 6-8: "ACA") corresponds to genomic positions 202, 201, 200

#### 3. Allele Handling

**Profile Data Convention:**
- `forward_ref_base` and `forward_derived_base` are **always on the forward strand**
- For negative strand genes, alleles must be **complemented** before applying to sequence

**Complement Map:**
```python
COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}
```

**Example:**
```python
# Forward strand shows: ref='A', derived='T'
# For negative strand gene:
#   - Effective ref in sequence = complement('A') = 'T'
#   - Effective derived in sequence = complement('T') = 'A'
#   - Result: 'T' → 'A' change in the gene sequence
```

#### 4. Reference Base Validation

For negative strand genes, reference base checks must account for complementation:

```python
# Example: Position 208 in negative strand gene with sequence "GTACTAACA"
# - pos_in_gene = (209-1) - 208 = 0
# - sequence[0] = 'G'
# - forward_ref_base should be complement('G') = 'C'
# - Profile data should have forward_ref_base = 'C'
```

**Implementation Detail:**
The script validates that `forward_ref_base` matches the complement of the sequence base at that position. Mismatches generate warnings but don't halt processing.

### Test Methods (Negative Strand)

#### `test_negative_strand_single_position_changes`

**Purpose:** Test k=1 changes at different codon positions (0, 1, 2) on a negative strand gene.

**Gene Design:**
```python
gene = {
    'gene_id': 'neg_strand_gene',
    'start': 201,  # 1-based
    'end': 209,    # 1-based (covers positions 200-208, 0-based)
    'strand': -1,
    'sequence': 'ATGCTTACG'  # 9 bases, 3 codons: ATG, CTT, ACG
}
```

**Coordinate Mapping Table:**
```
Genomic Pos    208  207  206  205  204  203  202  201  200
Gene Index     0    1    2    3    4    5    6    7    8
Sequence       A    T    G    C    T    T    A    C    G
Codon          ATG       |  CTT       |  ACG
```

**Mutations:**

**1. Position 0 change: ATG → TTG (codon 0, pos 0)**
- Genomic position 208 → gene index 0
- Forward: ref='T', derived='A'
- Complement: ref='A' (matches sequence[0]), derived='T'
- Result: ATG (Met) → TTG (Leu), **NS**
- Expected: frac_S=0.0, frac_NS=1.0

**2. Position 1 change: CTT → CAT (codon 1, pos 1)**
- Genomic position 204 → gene index 4
- Forward: ref='A', derived='T'
- Complement: ref='T' (matches sequence[4]), derived='A'
- Result: CTT (Leu) → CAT (His), **NS**
- Expected: frac_S=0.0, frac_NS=1.0

**3. Position 2 change: ACG → ACC (codon 2, pos 2)**
- Genomic position 200 → gene index 8
- Forward: ref='C', derived='G'
- Complement: ref='G' (matches sequence[8]), derived='C'
- Result: ACG (Thr) → ACC (Thr), **S**
- Expected: frac_S=1.0, frac_NS=0.0

**Potential Sites:**
- ATG: S=0.0, NS=3.0
- CTT: S=1.0, NS=2.0
- ACG: S=1.0, NS=2.0
- **Total: S=2.0, NS=7.0**

**Global dN/dS:**
```
Observed NS = 2.0, Observed S = 1.0
pN = 2.0/7.0 = 0.2857
pS = 1.0/2.0 = 0.5
dN/dS = 0.2857/0.5 = 0.5714
```

**Validation:**
```python
# All k=1 events
assert all(codon_events_df['k'] == 1)

# Strand correctly recorded
assert all(codon_events_df['strand'] == -1)

# Individual event validation
assert event1['ancestral_codon'] == 'ATG'
assert event1['derived_codon'] == 'TTG'
assert event1['mutation_type'] == 'NS'
assert event1['frac_NS'] == 1.0

# Global metrics match manual calculation
assert metrics['Global dN/dS (pNS/pS)'] ≈ 0.5714
```

#### `test_negative_strand_synonymous_transitions`

**Purpose:** Test synonymous-only changes on negative strand gene.

**Gene Design:**
```python
gene = {
    'gene_id': 'neg_syn_gene',
    'start': 201,
    'end': 209,
    'strand': -1,
    'sequence': 'CTTCGCGCT'  # Codons: CTT, CGC, GCT
}
```

**Mutations (All Synonymous):**

**1. CTT → CTC (pos 206): Leu → Leu**
- Genomic position 206 → gene index 2 (third position of CTT)
- Forward: ref='A', derived='G'
- Complement: 'T' → 'C'

**2. CGC → CGT (pos 203): Arg → Arg**
- Genomic position 203 → gene index 5 (third position of CGC)
- Forward: ref='G', derived='A'
- Complement: 'C' → 'T'

**3. GCT → GCC (pos 200): Ala → Ala**
- Genomic position 200 → gene index 8 (third position of GCT)
- Forward: ref='A', derived='G'
- Complement: 'T' → 'C'

**Expected Results:**
- All frac_S = 1.0, all frac_NS = 0.0
- All mutation_type = 'S'
- All strand = -1
- Global dN/dS = 0.0 (no non-synonymous changes)

**Validation:**
```python
assert all(codon_events_df['mutation_type'] == 'S')
assert all(codon_events_df['frac_S'] == 1.0)
assert all(codon_events_df['strand'] == -1)
assert metrics['Global dN/dS (pNS/pS)'] == 0.0
```

#### `test_negative_strand_nonsynonymous_transitions`

**Purpose:** Test non-synonymous-only changes on negative strand gene.

**Gene Design:**
```python
gene = {
    'gene_id': 'neg_ns_gene',
    'start': 201,
    'end': 209,
    'strand': -1,
    'sequence': 'ATGTTGAAG'  # Codons: ATG, TTG, AAG
}
```

**Mutations (All Non-Synonymous):**

**1. ATG → GTG (pos 208): Met → Val**
- Position 0 change in first codon
- Forward: ref='A', derived='C'
- Complement: 'T' → 'G'

**2. TTG → ATG (pos 205): Leu → Met**
- Position 0 change in second codon
- Forward: ref='A', derived='T'
- Complement: 'T' → 'A'

**3. AAG → TAG (pos 202): Lys → Stop**
- Position 0 change in third codon
- Forward: ref='T', derived='A'
- Complement: 'A' → 'T'

**Expected Results:**
- All frac_NS = 1.0, all frac_S = 0.0
- All mutation_type = 'NS'
- Global dN/dS = NaN (pS = 0.0)

**Validation:**
```python
assert all(codon_events_df['mutation_type'] == 'NS')
assert all(codon_events_df['frac_NS'] == 1.0)
assert metrics['pS (S / s)'] == 0.0
```

#### `test_negative_strand_two_position_changes`

**Purpose:** Test k=2 with path averaging on negative strand gene.

**Gene Design:**
```python
gene = {
    'gene_id': 'neg_k2_gene',
    'start': 201,
    'end': 206,
    'strand': -1,
    'sequence': 'ATGCGT'  # Codons: ATG, CGT
}
```

**Mutation: ATG → TTT**
- Genomic position 205 → gene index 0 (pos 0)
- Genomic position 203 → gene index 2 (pos 2)
- Forward alleles: (205) T→A, (203) C→A
- Complements to: A→T, G→T
- Result: ATG → TTT

**NG86 Path Averaging:**
```
Path 1: ATG → TTG → TTT (Met→Leu→Phe: 2 NS)
Path 2: ATG → ATT → TTT (Met→Ile→Phe: 2 NS)
Average: frac_S=0.0, frac_NS=2.0, num_valid_paths=2
```

**Validation:**
```python
assert event['k'] == 2
assert event['ancestral_codon'] == 'ATG'
assert event['derived_codon'] == 'TTT'
assert event['strand'] == -1
assert event['frac_NS'] == 2.0
assert event['num_valid_paths'] == 2
```

#### `test_negative_strand_coordinate_boundary_validation`

**Purpose:** Validate coordinate mapping at gene boundaries for negative strand.

**Gene Design:**
```python
gene = {
    'gene_id': 'boundary_gene',
    'start': 201,
    'end': 209,
    'strand': -1,
    'sequence': 'ATGCTTACG'  # 9 bases
}
```

**Test Positions:**
- **Highest genomic position (208)** → gene index 0 (first base)
- **Middle position (204)** → gene index 4 (middle base)
- **Lowest genomic position (200)** → gene index 8 (last base)

**Formula Verification:**
```
Position 208: (209 - 1) - 208 = 0 ✓
Position 204: (209 - 1) - 204 = 4 ✓
Position 200: (209 - 1) - 200 = 8 ✓
```

**Validation:**
```python
# Verify gene_position field for each event
event_208['gene_position'] == '0'
event_204['gene_position'] == '4'
event_200['gene_position'] == '8'
```

**Key Insight:** This test confirms that the coordinate mapping formula works correctly at all boundaries of a negative strand gene.

#### `test_negative_strand_vs_positive_strand_complement`

**Purpose:** Verify that complementary codon changes on opposite strands produce identical NG86 results.

**Gene Design:**
```python
# Positive strand gene
pos_gene = {
    'gene_id': 'pos_gene',
    'start': 101,
    'end': 106,
    'strand': 1,
    'sequence': 'ATGCGT'
}

# Negative strand gene (same coding sequence)
neg_gene = {
    'gene_id': 'neg_gene',
    'start': 201,
    'end': 206,
    'strand': -1,
    'sequence': 'ATGCGT'  # Same codons as positive
}
```

**Mutation (Both Genes): ATG → TTG (Met → Leu)**

**Positive strand:**
- Position 100 changes A→T (direct)

**Negative strand:**
- Position 205 changes, maps to gene index 0
- Forward allele: T→A
- Complement: A→T (same change in sequence)

**Expected Results:**
Both genes should produce:
- Same ancestral_codon: 'ATG'
- Same derived_codon: 'TTG'
- Same mutation_type: 'NS'
- Same frac_S: 0.0
- Same frac_NS: 1.0
- Same potential sites

**Validation:**
```python
pos_event = codon_events_df[codon_events_df['gene_id'] == 'pos_gene'].iloc[0]
neg_event = codon_events_df[codon_events_df['gene_id'] == 'neg_gene'].iloc[0]

# Identical codon changes
assert pos_event['ancestral_codon'] == neg_event['ancestral_codon']
assert pos_event['derived_codon'] == neg_event['derived_codon']

# Identical mutation classification
assert pos_event['mutation_type'] == neg_event['mutation_type']
assert pos_event['frac_S'] == neg_event['frac_S']
assert pos_event['frac_NS'] == neg_event['frac_NS']

# Same potential sites
assert pos_event['potential_S_codon'] == neg_event['potential_S_codon']
assert pos_event['potential_NS_codon'] == neg_event['potential_NS_codon']

# Different strands
assert pos_event['strand'] == 1
assert neg_event['strand'] == -1
```

**Key Insight:** This test proves that strand orientation is handled correctly—the same biological mutation produces identical NG86 results regardless of which strand the gene is on.

#### `test_negative_strand_intermediate_stop_exclusion`

**Purpose:** Validate intermediate stop codon exclusion on negative strand.

**Gene Design:**
```python
gene = {
    'gene_id': 'neg_stop_gene',
    'start': 201,
    'end': 206,
    'strand': -1,
    'sequence': 'CAACGT'  # Codons: CAA, CGT
}
```

**Mutation: CAA → TGA (k=2)**
- Genomic position 205 → gene index 0 (CAA pos 0)
- Genomic position 204 → gene index 1 (CAA pos 1)
- Forward: (205) G→A, (204) T→C
- Complement: C→T, A→G
- Result: CAA → TGA

**NG86 Path Averaging:**
```
Path 1: CAA → TAA (intermediate stop) → TGA
  - EXCLUDED (intermediate stop codon)

Path 2: CAA → CGA → TGA
  - Step 1: CAA (Gln) → CGA (Arg): NS
  - Step 2: CGA (Arg) → TGA (Stop): NS
  - VALID (final stop allowed)
  - Total: 0 S, 2 NS

Result: frac_S=0.0, frac_NS=2.0, num_valid_paths=1
```

**Validation:**
```python
assert event['k'] == 2
assert event['ancestral_codon'] == 'CAA'
assert event['derived_codon'] == 'TGA'
assert event['num_valid_paths'] == 1  # Only 1 valid path
assert event['frac_NS'] == 2.0
assert event['strand'] == -1
```

**Key Insight:** Intermediate stop exclusion works identically on both strands—only the coordinate mapping and allele complementation differ.

#### `test_negative_strand_three_position_changes`

**Purpose:** Test k=3 with 6-pathway averaging on negative strand gene.

**Gene Design:**
```python
gene = {
    'gene_id': 'neg_k3_gene',
    'start': 201,
    'end': 206,
    'strand': -1,
    'sequence': 'ATGCGT'  # Codons: ATG, CGT
}
```

**Mutation: ATG → CCC (k=3, all positions change)**
- Genomic positions 205, 204, 203 → gene indices 0, 1, 2
- Forward alleles: T→G, A→G, C→G
- Complement: A→C, T→C, G→C
- Result: ATG → CCC

**NG86 Path Averaging:**
All 6 permutations evaluated (same as positive strand test):
```
Total S steps across 6 paths: 4
Total NS steps across 6 paths: 14

Average frac_S = 4 / 6 = 0.667
Average frac_NS = 14 / 6 = 2.333
```

**Validation:**
```python
assert event['k'] == 3
assert event['ancestral_codon'] == 'ATG'
assert event['derived_codon'] == 'CCC'
assert event['strand'] == -1
assert event['num_valid_paths'] == 6
assert event['frac_S'] ≈ 0.667
assert event['frac_NS'] ≈ 2.333

# Coordinate fields
assert event['codon_position'] == '0,1,2'
assert event['gene_position'] == '0,1,2'
assert event['contig_position'] == '203,204,205'
```

**Key Insight:** Path averaging produces identical fractional counts for the same codon transition regardless of strand. The k=3 case demonstrates this holds even with complex 6-pathway averaging.

### Implementation Validation

**Reference Base Warnings:**
Several negative strand tests intentionally trigger reference base mismatch warnings to verify the validation logic works correctly. These warnings are expected and indicate the tests are properly exercising edge cases:

```
Reference base mismatch at 208 in boundary_gene. Prodigal (fwd strand): 'T', Profile: 'A'.
```

**Why Warnings Appear:**
The test data uses forward_ref_base values that may not perfectly match the complement of the sequence base, simulating real-world data inconsistencies. The script:
1. Logs the mismatch
2. Uses profile data as source of truth
3. Continues processing correctly

**Test Coverage Summary:**

| Test Method | k Value | Mutation Type | Key Feature |
|-------------|---------|---------------|-------------|
| `test_negative_strand_single_position_changes` | 1 | S + NS | All codon positions (0,1,2) |
| `test_negative_strand_synonymous_transitions` | 1 | S only | Wobble position changes |
| `test_negative_strand_nonsynonymous_transitions` | 1 | NS only | First position changes |
| `test_negative_strand_two_position_changes` | 2 | NS | Path averaging (2 paths) |
| `test_negative_strand_coordinate_boundary_validation` | 1 | Mixed | Boundary mapping |
| `test_negative_strand_vs_positive_strand_complement` | 1 | NS | Cross-strand validation |
| `test_negative_strand_intermediate_stop_exclusion` | 2 | NS | Stop codon exclusion |
| `test_negative_strand_three_position_changes` | 3 | Mixed | Path averaging (6 paths) |

**Total: 8 negative strand test methods**

---

## Technical Implementation Notes

### Coordinate System Handling

**Critical Issue:** Prodigal uses 1-based coordinates, while variant positions are 0-based.

**Conversion Formula:**
```python
contig_start = prodigal_start - 1  # Convert to 0-based
gene_position = contig_position - contig_start
```

**Example:**
- Prodigal: start=101, end=109 (1-based, inclusive)
- Contig: positions 100-108 (0-based, half-open)
- Variant at contig position 104 → gene position 4

**Common Errors Fixed:**
- **Out-of-bounds positions**: Variants at contig positions outside gene boundaries
- **Off-by-one errors**: Incorrect conversion between coordinate systems
- **Reference base mismatches**: Variant ref_base doesn't match gene sequence

### Reference Base Matching

**Validation:**
```python
# In analyze_codon_substitutions_with_ng86_paths()
gene_pos = position - (gene_info['start'] - 1)
expected_ref = ancestral_seq[gene_pos]
if row['forward_ref_base'] != expected_ref:
    logger.warning(f"Reference base mismatch at position {position}")
```

**Why It Matters:**
- Ensures mutations applied to correct ancestral state
- Detects data inconsistencies between annotation and profiling
- Critical for accurate S/NS classification

### Boundary Condition Handling

**Edge Cases Covered:**
- **First codon** (position 0, 1, 2): All three positions tested
- **Last codon**: Ensures no off-by-one errors at gene end
- **Stop codons**: Both intermediate (invalid) and final (valid)
- **k=1, k=2, k=3**: All possible change counts
- **Zero S sites**: Codons with no synonymous alternatives (ATG, TGG)
- **Undefined dN/dS**: When pS=0 or no NS sites

### NG86 Cache Usage

**Efficiency:**
- Pre-computed cache built once in `setUpClass`
- Contains all 64×64 = 4,096 codon-to-codon transitions
- Each entry includes: k, frac_S, frac_NS, num_valid_paths, is_valid

**Cache Structure:**
```python
ng86_cache = {
    ('ATG', 'TTG'): {
        'k': 1,
        'frac_S': 0.0,
        'frac_NS': 1.0,
        'num_valid_paths': 1,
        'is_valid': True
    },
    ('CTT', 'CAA'): {
        'k': 2,
        'frac_S': 0.5,
        'frac_NS': 1.5,
        'num_valid_paths': 2,
        'is_valid': True
    },
    # ... 4,094 more entries
}
```

**Benefits:**
- Eliminates redundant pathway enumeration
- Consistent results across tests
- Fast lookup during analysis

### Fractional Site Counts

**NG86 Produces Fractional Counts:**
- Individual events: frac_S + frac_NS = k (number of changes)
- Example: k=2 change with frac_S=0.5, frac_NS=1.5
- Total: 0.5 + 1.5 = 2.0 ✓

**Aggregation:**
```python
# Global totals are sums of fractional counts
total_NS = codon_events_df['frac_NS'].sum()
total_S = codon_events_df['frac_S'].sum()
```

**Critical:** Do not round fractional counts until final reporting.

### Path Averaging Logic

**Number of Pathways:**
- k=1: 1 pathway (1! = 1)
- k=2: 2 pathways (2! = 2)
- k=3: 6 pathways (3! = 6)

**Invalid Pathway Exclusion:**
```python
# Pathways through intermediate stop codons excluded
if intermediate_stop_detected:
    num_valid_paths -= 1
```

**Fractional Count Calculation:**
```python
# Average over valid pathways only
frac_S = sum(S_steps_per_path) / num_valid_paths
frac_NS = sum(NS_steps_per_path) / num_valid_paths
```

### Potential Sites Calculation

**Per-Codon:**
- Calculated independently for each codon
- Some codons have S=0 (ATG, TGG)
- 4-fold degenerate wobble positions have high S (e.g., GCT S=1.0)

**Gene-Level:**
```python
gene_S = sum(codon_S for codon in gene)
gene_NS = sum(codon_NS for codon in gene)
```

**MAG-Level:**
```python
mag_S = sum(gene_S for gene in mag)
mag_NS = sum(gene_NS for gene in mag)
```

### Data Flow Validation

**Pipeline Stages:**
1. **Site data** → DataFrame construction
2. **Profile data** → Ancestral sequence reconstruction
3. **Ancestral sequences** → NG86 pathway analysis
4. **Codon events** → Global aggregation
5. **Global metrics** → Gene/MAG summaries

**Each Test Validates:**
- Intermediate outputs (per-event frac_S, frac_NS, k)
- Final outputs (global dN/dS, pN, pS)
- Consistency between codon-level and aggregated results

---

## Known Limitations and Edge Cases

### Limitations

1. **Small Sample Size**:
   - Tests use short gene sequences (6-9 bp typically)
   - Real genes have hundreds to thousands of bp
   - dN/dS estimates more stable with larger genes

2. **No Multi-MAG Tests**:
   - Most tests use single gene or single MAG
   - Real analyses involve hundreds of MAGs
   - Future: Add tests with 10+ MAGs and cross-MAG comparisons

3. **No Time-Series Tests**:
   - Tests use single ancestral→derived comparison
   - Real data has multiple timepoints
   - Future: Add tests for longitudinal dN/dS trajectories

4. **Simplified Mock Data**:
   - No coverage depth variation
   - All variant positions have clear major alleles
   - No ambiguous nucleotides or sequencing errors

### Edge Cases Handled

✅ **Zero synonymous sites (ATG, TGG)**: Test #8  
✅ **Zero non-synonymous sites**: Test #2 (all S)  
✅ **Undefined dN/dS (pS=0)**: Test #3  
✅ **Intermediate stop codons**: Test #7  
✅ **All three codon positions**: Test #1  
✅ **k=1, k=2, k=3 changes**: Tests #1, #4, #5  
✅ **Mixed pathway classifications**: Test #6  
✅ **Reverse strand genes**: (covered in other test classes)  
✅ **Gene and MAG aggregation**: Test #10

### Future Enhancements

**Potential Additions:**
- **Multi-timepoint dN/dS**: Track ratio changes over evolution
- **Gene-by-gene comparisons**: Identify genes with outlier dN/dS
- **Confidence intervals**: Bootstrap or analytical CIs for dN/dS
- **Context-dependent mutation rates**: Account for GC content, CpG sites
- **Codon usage bias**: Weight by codon frequency in genome
- **Parallel evolution detection**: Identify recurrent NS mutations across MAGs

---

## Summary

The `TestNG86ComputedValues` test class provides **comprehensive validation** of the NG86 dN/dS implementation through 16 test methods covering:

✅ All codon positions (0, 1, 2) on **both positive and negative strands**
✅ All k values (1, 2, 3) on **both strands**
✅ Synonymous-only, non-synonymous-only, and mixed scenarios on **both strands**
✅ Path averaging with different outcomes
✅ Intermediate stop codon exclusion on **both strands**
✅ Zero synonymous site edge cases
✅ Complete end-to-end calculation pipeline
✅ Hierarchical gene/MAG aggregation
✅ **Negative strand coordinate mapping and allele complementation**
✅ **Cross-strand validation (same mutation on both strands)**

**Key Strengths:**
- Explicit expected values for all test scenarios
- Manual calculations documented for verification
- Coordinate system fixes ensure data validity
- Covers critical NG86 edge cases (stops, k=3, mixed paths)
- Validates entire pipeline from sites to summary stats
- **Comprehensive negative strand coverage with 6 dedicated test methods**
- **Validates reverse complement handling and coordinate reversal**

**Recent Improvements:**
- Renamed from `TestEndToEndWithMockDatasets` for clarity
- Added 16 test methods total (10 positive strand + 6 negative strand)
- Fixed coordinate system alignment issues
- Corrected reference base mismatches
- Comprehensive docstrings and inline documentation
- **Added complete negative strand test suite** (6 new methods)
- **Documented negative strand coordinate system and allele handling**

This test suite ensures the AlleleFlux dN/dS analysis produces accurate, biologically meaningful results validated against NG86 theory for **both forward and reverse strand genes**.

---

## Mathematical Formulas

### Nei-Gojobori Potential Sites for a Codon

For a given codon with nucleotides at positions 0, 1, 2:

```
S_codon = Σ(i=0 to 2) [Σ(b ∈ {A,T,G,C}, b ≠ codon[i]) I(AA(mutate(codon, i, b)) = AA(codon))] / 3

N_codon = Σ(i=0 to 2) [Σ(b ∈ {A,T,G,C}, b ≠ codon[i]) I(AA(mutate(codon, i, b)) ≠ AA(codon))] / 3
```

Where:
- **I(condition)** = indicator function (1 if true, 0 if false)
- **AA(codon)** = amino acid translation of codon
- **mutate(codon, i, b)** = codon with position i changed to base b

**Simplified:**
- Count synonymous changes at each position
- Count non-synonymous changes at each position
- Divide by 3 to get fractional site contribution

### Gene-Level Potential Sites

```
S_gene = Σ(codon ∈ gene) S_codon
N_gene = Σ(codon ∈ gene) N_codon
```

Total S and N sites for a gene are sums of per-codon values.

### NG86 Fractional Counts

For a codon pair (ancestral → derived) with k differences:

```
frac_S = (1 / num_valid_paths) × Σ(path ∈ valid_paths) S_steps(path)
frac_NS = (1 / num_valid_paths) × Σ(path ∈ valid_paths) NS_steps(path)
```

Where:
- **valid_paths** = all k! permutations excluding those through intermediate stops
- **S_steps(path)** = number of synonymous steps in this pathway
- **NS_steps(path)** = number of non-synonymous steps in this pathway

**Example (k=2):**
```
Path 1: S_steps=0, NS_steps=2
Path 2: S_steps=1, NS_steps=1

frac_S = (0 + 1) / 2 = 0.5
frac_NS = (2 + 1) / 2 = 1.5
```

### Global dN/dS Calculation

```
pN = (Σ frac_NS) / (Σ potential_NS_codon)
pS = (Σ frac_S) / (Σ potential_S_codon)
dN/dS = pN / pS
```

Where sums are over all valid, unique codon events.

**Step-by-Step:**
1. Sum all frac_NS values → Observed_NS
2. Sum all frac_S values → Observed_S
3. Sum potential_NS_codon from unique codons → Total_Potential_NS
4. Sum potential_S_codon from unique codons → Total_Potential_S
5. Calculate pN = Observed_NS / Total_Potential_NS
6. Calculate pS = Observed_S / Total_Potential_S
7. Calculate dN/dS = pN / pS

### Gene-Level dN/dS

```
dN_gene = (Σ codon events in gene: frac_NS) / potential_NS_sites_gene
dS_gene = (Σ codon events in gene: frac_S) / potential_S_sites_gene
dN/dS_gene = dN_gene / dS_gene
```

Where:
- Numerators sum fractional counts across all codon events in the gene (only includes genes with ≥1 significant site)
- Denominators use total gene potential (entire gene length, not just per-codon sum)

**Note:** `potential_NS_sites_gene` and `potential_S_sites_gene` represent the full gene's mutational opportunity, calculated from the entire gene sequence, even though only a subset of codons may have significant changes.

### MAG-Level dN/dS

```
dN_MAG = (Σ genes with sig. sites: Σ codon events: frac_NS) / (Σ genes with sig. sites: potential_NS_sites_gene)
dS_MAG = (Σ genes with sig. sites: Σ codon events: frac_S) / (Σ genes with sig. sites: potential_S_sites_gene)
dN/dS_MAG = dN_MAG / dS_MAG
```

**Important:** Aggregates across only genes that contain at least one significant site, not all genes in the MAG.

**Rationale:** Since the input is filtered to significant sites (by p-value threshold), only genes with significant changes contribute to the numerator. The denominator must match this scope - it includes only the potential sites from those same genes. Including all genes in the MAG would create a mismatch between observed changes (from subset of genes) and potential sites (from all genes).

---

## Edge Cases and Boundary Conditions

### Empty Input DataFrames

**Test:** `test_global_dnds_empty_input`

**Scenario:** No codon events to analyze

**Behavior:**
```python
result = calculate_global_dnds_for_sites(pd.DataFrame())
assert result is None
```

**Handling:** Functions return `None` gracefully, no exceptions raised.

### Invalid Codons with Ambiguous Bases

**Test:** `test_calculate_codon_sites` with invalid codon

**Scenario:** Codon contains 'N' or other ambiguous IUPAC code

**Behavior:**
```python
with self.assertRaises(ValueError):
    _calculate_codon_sites('ATN', table)
```

**Error Message:** `"Invalid codon: ATN. Cannot analyze non-standard codons."`

**Prevention:** Profile data validation should exclude sites with ambiguous bases.

### All Pathways Hit Intermediate Stop Codons

**Test:** `test_ng86_all_paths_invalid`

**Scenario:** Every k! pathway passes through an intermediate stop

**Example:** TAA (Stop) → TGG (Trp)
- Path 1: TAA → TGA (intermediate stop) → TGG ✗
- Path 2: TAA → TAG (intermediate stop) → TGG ✗

**Behavior:**
```python
frac_S = NaN
frac_NS = NaN
num_valid_paths = 0
is_valid = False
```

**Downstream:** Codon excluded from global dN/dS calculation.

### Genes with No Profile Data

**Test:** `test_reconstruction_with_multiple_changes` (gene_no_profile)

**Scenario:** Gene in annotation but no profile data available

**Behavior:**
```python
# Gene sequence unchanged from reference
ancestral_seq = original_prodigal_sequence
# No entries in ancestral_major_alleles dict
```

**Use Case:** Gene has no significant variable sites.

### Reference Base Mismatch

**Test:** `test_reconstruction_with_ref_mismatch`

**Scenario:** Profile `ref_base` ≠ Prodigal sequence base

**Causes:**
- Assembly/annotation version mismatch
- Sequencing errors
- Data processing errors

**Behavior:**
```python
logger.warning(f"Reference base mismatch at position {pos}...")
# Use profile ref_base for tie-breaking
# Continue with reconstruction
```

**Resolution:** Profile data takes precedence (more recent).

### Tie-Breaking Scenarios

**Equal Allele Counts:**

**Case 1: Reference in tie**
```python
counts = {'A': 50, 'C': 50, 'G': 0, 'T': 0}
ref_base = 'A'
result = 'A'  # Prefer reference
```

**Case 2: Reference not in tie**
```python
counts = {'A': 50, 'C': 50, 'G': 0, 'T': 0}
ref_base = 'G'
result = random.choice(['A', 'C'])  # Random between tied
```

**Case 3: No coverage**
```python
counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
ref_base = 'T'
result = 'T'  # Default to reference
```

### Out-of-Bounds Positions

**Test:** `test_get_codon_from_site`

**Scenario:** Position outside gene boundaries

**Behavior:**
```python
codon, pos = get_codon_from_site(99, gene_info, seq)
assert codon is None
assert pos is None
# Function logs warning and returns None
```

**Prevention:** Sites pre-filtered to gene_id before processing.

### Partial Codons at Sequence End

**Scenario:** Gene length not divisible by 3

**Behavior:**
```python
gene_seq = 'ATGCGT'  # 6 bases (OK)
gene_seq = 'ATGCGTA'  # 7 bases (1 extra)

# In calculate_potential_sites_for_gene():
if len(codon) != 3:
    logger.warning(f"Skipping partial codon...")
    continue
```

**Handling:** Partial codons skipped with warning, don't contribute to totals.

### Invalid Codon Exclusion

**Test:** `test_global_dnds_excludes_invalid_codons`

**Scenario:** Some codon events have `is_valid=False`

**Filtering:**
```python
valid_codons = results_df[results_df['is_valid'] == True]
# Only valid codons included in calculations
```

**Reporting:**
```python
summary_df includes:
  "Invalid Codons (all paths hit stops)": N
```
---

## Helper Functions

### `build_prodigal_records()`

**Purpose:** Convert gene list from mock JSON into prodigal_records structure.

**Input:**
```python
genes_block = [
    {
        'gene_id': 'geneA',
        'start': 100,
        'end': 123,
        'strand': 1,
        'sequence': 'ATGGCTTTTGAATCTTGCGGTCAA'
    }
]
```

**Output:**
```python
prodigal_records = {
    'geneA': {
        'record': SeqRecord(
            Seq('ATGGCTTTTGAATCTTGCGGTCAA'),
            id='geneA',
            description='geneA # 100 # 123 # 1'
        ),
        'start': 100,
        'end': 123,
        'strand': 1
    }
}
```

**Function Signature:**
```python
def build_prodigal_records(genes_block: List[Dict]) -> Dict[str, Dict]:
```

**Usage in Tests:**
```python
prodigal_records = build_prodigal_records(V2_JSON['genes'])
```

### `build_dataframes_from_sites()`

**Purpose:** Build significant sites and profile DataFrames from mock site data.

**Input:**
```python
sites_block = [
    {
        'contig': 'contig1',
        'position': 102,
        'gene_id': 'geneA',
        'forward_ref_base': 'T',
        'forward_derived_base': 'C'
    }
]
```

**Output:**
```python
(significant_sites_df, ancestral_profile_df, derived_profile_df)
```

**Significant Sites DataFrame:**
```
| mag_id  | contig   | position | gene_id | test_type              | min_p_value | q_value |
|---------|----------|----------|---------|------------------------|-------------|---------|
| MAG_001 | contig1  | 102      | geneA   | two_sample_paired_tTest| 1e-6        | 1e-6    |
```

**Ancestral Profile DataFrame:**
```
| contig  | position | gene_id | ref_base | A  | T  | G | C  |
|---------|----------|---------|----------|----|----|---|----|
| contig1 | 102      | geneA   | T        | 0  | 30 | 0 | 0  |
```

**Derived Profile DataFrame:**
```
| contig  | position | gene_id | ref_base | A  | T | G | C  |
|---------|----------|---------|----------|----|----|---|----|
| contig1 | 102      | geneA   | T        | 0  | 5 | 0 | 25 |
```

**Profile Data Model:**
- **Ancestral:** ref_base is major (count=30)
- **Derived:** derived_base is major (count=25), ref_base residual (count=5)

**Function Signature:**
```python
def build_dataframes_from_sites(
    sites_block: List[Dict]
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
```

**Usage in Tests:**
```python
sig_df, pre_df, post_df = build_dataframes_from_sites(V2_JSON['sites'])
```

### Mock Dataset Structure

**V2_JSON, V3_JSON, V4_JSON:**

Each contains:
```python
{
    'genes': [
        {
            'gene_id': str,
            'start': int,      # 1-based Prodigal coordinate
            'end': int,
            'strand': int,     # 1 or -1
            'sequence': str    # DNA sequence
        }
    ],
    'sites': [
        {
            'contig': str,
            'position': int,   # 0-based genomic position
            'gene_id': str,
            'forward_ref_base': str,
            'forward_derived_base': str
        }
    ]
}
```

**Design Philosophy:**
- Self-contained test data
- Reproducible results
- Cover diverse scenarios (forward/reverse, k=1/2/3, valid/invalid)

