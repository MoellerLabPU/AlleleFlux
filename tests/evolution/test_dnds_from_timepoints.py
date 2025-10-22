#!/usr/bin/env python3
"""
Consolidated test suite for dN/dS analysis with NG86 path averaging.

Updated for NG86 methodology:
- Codon-level analysis (not per-site)
- Fractional S/NS counts from path averaging
- New columns: frac_S, frac_NS, k, is_valid, num_valid_paths
- Helper functions updated to call analyze_codon_substitutions_with_ng86_paths()
"""

import json
import logging
import random
import unittest
from typing import Dict, List, Tuple
from unittest import mock

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Import functions from the updated script
from alleleflux.scripts.evolution.dnds_from_timepoints import (
    _calculate_codon_sites,
    _enumerate_ng86_paths,
    _get_ng86_cache,
    analyze_codon_substitutions_with_ng86_paths,
    analyze_mutation_effect,
    calculate_global_dnds_for_sites,
    calculate_potential_sites_for_gene,
    get_codon_from_site,
    get_major_allele,
    reconstruct_ancestral_sequences,
    summarize_results,
)

# =============================================================================
# HELPER FUNCTIONS (Updated for NG86)
# =============================================================================


def build_prodigal_records(genes_block: List[Dict]) -> Dict[str, Dict]:
    """
    Convert gene list into prodigal_records structure.

    Args:
        genes_block: List of gene dictionaries with keys: gene_id, start, end, strand, sequence

    Returns:
        Dict mapping gene_id to gene metadata including record, start, end, strand
    """
    records = {}
    for g in genes_block:
        gene_id = g["gene_id"]
        seq = g["sequence"].upper()
        records[gene_id] = {
            "record": SeqRecord(
                Seq(seq),
                id=gene_id,
                description=f"{gene_id} # {g['start']} # {g['end']} # {g['strand']}",
            ),
            "start": g["start"],
            "end": g["end"],
            "strand": g["strand"],
        }
    return records


def build_dataframes_from_sites(
    sites_block: List[Dict],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build significant sites and profile DataFrames from mock data.

    Args:
        sites_block: List of site dictionaries with keys: contig, position, gene_id,
                     forward_ref_base, forward_derived_base

    Returns:
        Tuple of (significant_sites_df, ancestral_profile_df, derived_profile_df)
    """
    sig_rows = []
    pre_rows = []
    post_rows = []

    for s in sites_block:
        contig = s["contig"]
        pos = s["position"]
        gene_id = s["gene_id"]
        ref = s["forward_ref_base"].upper()
        der = s["forward_derived_base"].upper()

        # Significant site entry
        sig_rows.append(
            {
                "mag_id": "MAG_001",
                "contig": contig,
                "position": pos,
                "gene_id": gene_id,
                "test_type": "two_sample_paired_tTest",
                "min_p_value": 1e-6,
                "q_value": 1e-6,
            }
        )

        # Ancestral profile: ref is major allele (count=30)
        pre_count = {"A": 0, "T": 0, "G": 0, "C": 0}
        pre_count[ref] = 30
        pre_rows.append(
            {
                "contig": contig,
                "position": pos,
                "gene_id": gene_id,
                "ref_base": ref,
                **pre_count,
            }
        )

        # Derived profile: derived is major (count=25), ref residual (count=5)
        post_count = {"A": 0, "T": 0, "G": 0, "C": 0}
        post_count[der] = 25
        post_count[ref] += 5
        post_rows.append(
            {
                "contig": contig,
                "position": pos,
                "gene_id": gene_id,
                "ref_base": ref,
                **post_count,
            }
        )

    sig_df = pd.DataFrame(sig_rows)
    pre_df = pd.DataFrame(pre_rows)
    post_df = pd.DataFrame(post_rows)
    return sig_df, pre_df, post_df


# =============================================================================
# TEST CLASSES
# =============================================================================


class TestReconstructAncestralSequences(unittest.TestCase):
    """
    Unit tests for the reconstruct_ancestral_sequences function.
    Tests ancestral sequence reconstruction with multiple changes per gene.
    """

    def setUp(self):
        """
        Set up mock data and disable logging for clean test output.
        This method is run before each test in this class.
        """
        # --- 1. Mock Prodigal Records (Gene Templates) ---
        # A dictionary mapping gene IDs to their sequence, coordinates, and strand.
        # The description string is formatted to match real Prodigal output.
        self.prodigal_records = {
            "gene_fwd_1": {
                "record": SeqRecord(
                    Seq("ATGCGTACG"),
                    id="gene_fwd_1",
                    description="gene_fwd_1 # 101 # 109 # 1",
                ),
                "start": 101,
                "end": 109,
                "strand": 1,
            },
            "gene_rev_1": {
                "record": SeqRecord(
                    Seq("GTACTAACA"),
                    id="gene_rev_1",
                    description="gene_rev_1 # 201 # 209 # -1",
                ),
                "start": 201,
                "end": 209,
                "strand": -1,
            },
            "gene_tie_break": {
                "record": SeqRecord(
                    Seq("AAAAAA"),
                    id="gene_tie_break",
                    description="gene_tie_break # 301 # 306 # 1",
                ),
                "start": 301,
                "end": 306,
                "strand": 1,
            },
            "gene_no_profile": {
                "record": SeqRecord(
                    Seq("GGGGGG"),
                    id="gene_no_profile",
                    description="gene_no_profile # 401 # 406 # 1",
                ),
                "start": 401,
                "end": 406,
                "strand": 1,
            },
        }

        # --- 2. Mock Ancestral Profile Data (Corrected) ---
        # The 'position' values now correctly target the bases needed to
        # satisfy the assertions in the test method.
        profile_data = {
            "contig": [
                "contig1",
                "contig1",
                "contig1",
                "contig1",
                "contig1",
                "contig1",
                "contig2",
            ],
            "position": [
                100,  # Changed positions for gene_fwd_1
                102,  # Changed positions for gene_fwd_1
                108,  # Changed positions for gene_fwd_1
                202,  # Changed positions for gene_rev_1
                203,  # Changed positions for gene_rev_1
                208,  # Changed positions for gene_rev_1
                301,  # Changed position for gene_tie_break
            ],
            # ref_base from profile. Corresponds to the forward-strand base.
            "ref_base": ["A", "G", "G", "T", "T", "C", "A"],
            "gene_id": [
                "gene_fwd_1",
                "gene_fwd_1",
                "gene_fwd_1",
                "gene_rev_1",
                "gene_rev_1",
                "gene_rev_1",
                "gene_tie_break",
            ],
            # Allele counts are set up to cause the specific changes we want to test.
            "A": [5, 5, 90, 10, 80, 2, 20],
            "C": [2, 95, 2, 90, 2, 3, 20],
            "G": [3, 0, 3, 8, 8, 4, 5],
            "T": [90, 0, 5, 0, 0, 88, 1],
        }
        ancestral_profile_df = pd.DataFrame(profile_data)

        # The function expects a grouped DataFrame, so we create it here.
        self.profile_by_gene = ancestral_profile_df.groupby("gene_id")

        # --- 3. List of Genes to Process ---
        self.unique_genes = [
            "gene_fwd_1",
            "gene_rev_1",
            "gene_tie_break",
            "gene_no_profile",
            "gene_is_missing",
        ]

    def test_reconstruction_with_multiple_changes(self):
        """
        Main test method to run the reconstruction and assert outcomes.
        """
        # --- Run the function with our mock data ---
        reconstructed_seqs, ancestral_major_alleles = reconstruct_ancestral_sequences(
            self.unique_genes, self.prodigal_records, self.profile_by_gene
        )

        # --- Assertions: Check if the output is as expected ---

        # Case 1: Forward strand gene with multiple changes
        # Original: ATGCGTACG
        # Change 1: pos 100 (index 0, 'A') -> major allele 'T' -> TTGCGTACG
        # Change 2: pos 102 (index 2, 'G') -> major allele 'C' -> TTCCGTACG
        # Change 3: pos 108 (index 8, 'G') -> major allele 'A' -> TTCCGTACA
        self.assertIn("gene_fwd_1", reconstructed_seqs)
        self.assertEqual(reconstructed_seqs["gene_fwd_1"], "TTCCGTACA")

        # Case 2: Reverse strand gene with multiple changes
        # Original: GTACTAACA
        # Change 1: pos 208 (index 0, 'G') -> major allele 'T' (complement 'A') -> ATACTAACA
        # Change 2: pos 203 (index 5, 'A') -> major allele 'A' (complement 'T') -> ATACTTACA
        # Change 3: pos 202 (index 6, 'A') -> major allele 'C' (complement 'G') -> ATACTTGCA
        self.assertIn("gene_rev_1", reconstructed_seqs)
        self.assertEqual(reconstructed_seqs["gene_rev_1"], "ATACTTGCA")

        # Case 3: Tie-breaking rule (reference base is in the tie)
        # Original: AAAAAA. Position 301 (index 0). Ref base is 'A'.
        # Profile has a tie between 'A' (20) and 'C' (20).
        # The script should prefer the reference base 'A'. No change is expected.
        self.assertIn("gene_tie_break", reconstructed_seqs)
        self.assertEqual(reconstructed_seqs["gene_tie_break"], "AAAAAA")

        # Case 4: Gene with no entry in the profile data
        # Expected: GGGGGG (the original sequence)
        self.assertIn("gene_no_profile", reconstructed_seqs)
        self.assertEqual(reconstructed_seqs["gene_no_profile"], "GGGGGG")

        # Case 5: Gene missing from prodigal records
        # The function should skip this gene, so it shouldn't be in the output.
        self.assertNotIn("gene_is_missing", reconstructed_seqs)

        # Check total number of reconstructed genes
        self.assertEqual(len(reconstructed_seqs), 4)

        # --- Assertions for Ancestral Major Alleles Dictionary ---
        expected_major_alleles = {
            ("contig1", 100): "T",
            ("contig1", 102): "C",
            ("contig1", 108): "A",
            ("contig1", 202): "C",
            ("contig1", 203): "A",
            ("contig1", 208): "T",
            ("contig2", 301): "A",  # Tie broken by profile 'ref_base'
        }
        self.assertDictEqual(ancestral_major_alleles, expected_major_alleles)

    def test_reconstruction_with_ref_mismatch(self):
        """Test that warning is logged when ref_base mismatches prodigal."""
        # --- Modify a single ref_base to create a mismatch ---
        # For gene_rev_1 at pos 208, prodigal reverse base is 'G', so fwd is 'C'.
        # We set the profile's ref_base to 'A' to force a mismatch.
        profile_df = self.profile_by_gene.obj.copy()
        profile_df.loc[profile_df["position"] == 208, "ref_base"] = "A"
        profile_by_gene_mismatch = profile_df.groupby("gene_id")

        with self.assertLogs(
            "alleleflux.scripts.evolution.dnds_from_timepoints", level="WARNING"
        ) as cm:
            reconstructed_seqs, _ = reconstruct_ancestral_sequences(
                self.unique_genes, self.prodigal_records, profile_by_gene_mismatch
            )
            self.assertTrue(any("Reference base mismatch" in msg for msg in cm.output))

        # Even with the warning, the reconstruction should still complete correctly
        # using the major allele from the profile data.
        self.assertEqual(reconstructed_seqs["gene_rev_1"], "ATACTTGCA")


class TestCoreDndsFunctions(unittest.TestCase):
    """Unit tests for core dN/dS logic and helper functions."""

    def setUp(self):
        self.gene_info_fwd = {
            "record": SeqRecord(Seq("ATGCGTACG"), id="gene_fwd_1"),
            "start": 101,
            "end": 109,
            "strand": 1,
        }
        self.gene_info_rev = {
            "record": SeqRecord(Seq("GTACTAACA"), id="gene_rev_1"),
            "start": 201,
            "end": 209,
            "strand": -1,
        }
        self.table = CodonTable.unambiguous_dna_by_id[11]

    def test_get_major_allele(self):
        """Tests the logic for selecting the major allele from coverage data."""
        self.assertEqual(get_major_allele({"A": 10, "C": 90}, "A"), "C")
        self.assertEqual(get_major_allele({"A": 50, "C": 50}, "A"), "A")
        with mock.patch("random.choice", return_value="C"):
            self.assertEqual(get_major_allele({"A": 50, "C": 50}, "G"), "C")
        self.assertEqual(get_major_allele({"A": 0, "C": 0}, "T"), "T")

    def test_calculate_codon_sites(self):
        """Test calculation of potential S and N sites for single codons."""
        # Methionine (ATG) has 0 synonymous sites
        s, n = _calculate_codon_sites("ATG", self.table)
        self.assertEqual(s, 0.0)
        self.assertEqual(n, 3.0)

        # Leucine (CTT) is 1/3 synonymous, 2/3 non-synonymous
        s, n = _calculate_codon_sites("CTT", self.table)
        self.assertAlmostEqual(s, 1.0)
        self.assertAlmostEqual(n, 2.0)

        # Stop codon (TAA)
        s, n = _calculate_codon_sites("TAA", self.table)
        self.assertAlmostEqual(s, 2 / 3)
        self.assertAlmostEqual(n, 4 / 3 + 1)

        # An invalid codon should raise an error
        with self.assertRaises(ValueError):
            _calculate_codon_sites("ATN", self.table)

    def test_calculate_potential_sites_for_gene(self):
        """Tests the aggregation of S and N sites over a full gene sequence."""
        # Gene = ATG CTT TAA (Met, Leu, Stop)
        # S sites = 0 (ATG) + 1 (CTT) + 2/3 (TAA) = 1.666...
        # N sites = 3 (ATG) + 2 (CTT) + 1 + 4/3 (TAA) = 7.333...
        gene_seq = "ATGCTTTAA"
        sites = calculate_potential_sites_for_gene(gene_seq)
        self.assertAlmostEqual(sites["S"], 1 + (2 / 3))
        self.assertAlmostEqual(sites["NS"], 6.0 + 4 / 3)
        # Gene = GAC ACA GCG GTT (Asp, Thr, Ala, Val)
        # S sites = 0.333... (GAC) + 1.0 (ACA) + 1.0 (GCG) + 1.0 (GTT) = 3 + 1/3
        # N sites = 2.667... (GAC) + 2.0 (ACA) + 2.0 (GCG) + 2.0 (GTT) = 8 + 2/3
        gene_seq = "GACACAGCGGTT"
        sites = calculate_potential_sites_for_gene(gene_seq)
        self.assertAlmostEqual(sites["S"], 3 + (1 / 3))
        self.assertAlmostEqual(sites["NS"], 8 + (2 / 3))

    def test_get_codon_from_site(self):
        """Test mapping of genomic position to correct codon."""
        seq_fwd = "ATGCGTACG"  # Codons: ATG, CGT, ACG
        codon, pos = get_codon_from_site(104, self.gene_info_fwd, seq_fwd)
        self.assertEqual(codon, "CGT")
        self.assertEqual(pos, 4)

        # Reverse strand
        seq_rev = "GTACTAACA"  # Codons: GTA, CTA, ACA
        codon, pos = get_codon_from_site(203, self.gene_info_rev, seq_rev)
        self.assertEqual(codon, "CTA")
        self.assertEqual(pos, 5)

        # Out of bounds
        codon, pos = get_codon_from_site(99, self.gene_info_fwd, seq_fwd)
        self.assertIsNone(codon)

    def test_analyze_mutation_effect(self):
        """Test classification of mutations as S or NS."""
        # Forward strand, non-synonymous (Arg -> Leu)
        ancestral_seq = "ATGCGTACG"
        result = analyze_mutation_effect(self.gene_info_fwd, ancestral_seq, 104, "T")
        self.assertEqual(result["mutation_type"], "NS")
        self.assertEqual(result["codon_before"], "CGT")
        self.assertEqual(result["codon_after"], "CTT")
        self.assertEqual(result["aa_before"], "R")
        self.assertEqual(result["aa_after"], "L")

        # Reverse strand, synonymous (Leu -> Leu)
        # This tests that the derived allele is correctly complemented.
        ancestral_seq_rev = "GTACTAACA"
        result = analyze_mutation_effect(
            self.gene_info_rev, ancestral_seq_rev, 203, "G"
        )
        self.assertEqual(result["mutation_type"], "S")
        self.assertEqual(result["codon_before"], "CTA")
        self.assertEqual(result["codon_after"], "CTC")
        self.assertEqual(result["aa_before"], "L")
        self.assertEqual(result["aa_after"], "L")


class TestNG86PathAveraging(unittest.TestCase):
    """
    Unit tests for NG86 path averaging functions.

    Tests the core NG86 functionality including:
    - Single nucleotide changes (k=1)
    - Double changes (k=2) with path averaging
    - Triple changes (k=3) with path averaging
    - Invalid paths (intermediate stop codons)
    """

    def setUp(self):
        self.table = CodonTable.unambiguous_dna_by_id[11]
        # Build NG86 cache once for all tests
        self.ng86_cache = _get_ng86_cache()

    def test_ng86_single_change_synonymous(self):
        """Test k=1 synonymous change (no path averaging needed)."""
        # CTT (Leu) -> CTC (Leu) - synonymous
        frac_S, frac_NS, num_paths = _enumerate_ng86_paths("CTT", "CTC", self.table)
        self.assertEqual(frac_S, 1.0)
        self.assertEqual(frac_NS, 0.0)
        self.assertEqual(num_paths, 1)

    def test_ng86_single_change_nonsynonymous(self):
        """Test k=1 non-synonymous change."""
        # ATG (Met) -> TTG (Leu) - non-synonymous
        frac_S, frac_NS, num_paths = _enumerate_ng86_paths("ATG", "TTG", self.table)
        self.assertEqual(frac_S, 0.0)
        self.assertEqual(frac_NS, 1.0)
        self.assertEqual(num_paths, 1)

    def test_ng86_double_change_both_ns(self):
        """Test k=2 where both paths give NS steps."""
        # ATG (Met) -> TTT (Phe), positions 0 and 2 differ
        # Path 1: ATG -> TTG -> TTT (Met->Leu->Phe: 2 NS)
        # Path 2: ATG -> ATT -> TTT (Met->Ile->Phe: 2 NS)
        # Average: (0+0)/2 = 0.0 S, (2+2)/2 = 2.0 NS
        frac_S, frac_NS, num_paths = _enumerate_ng86_paths("ATG", "TTT", self.table)
        self.assertEqual(frac_S, 0.0)
        self.assertEqual(frac_NS, 2.0)
        self.assertEqual(num_paths, 2)

    def test_ng86_double_change_mixed(self):
        """Test k=2 where paths differ in S/NS classification."""
        # CTT (Leu) -> CAA (Gln), positions 1 and 2 differ
        # Path 1: CTT -> CAT (H) -> CAA (Q): NS, NS = 2 NS
        # Path 2: CTT -> CTA (L) -> CAA (Q): S, NS = 1 S, 1 NS
        # Average: 0.5 S, 1.5 NS
        frac_S, frac_NS, num_paths = _enumerate_ng86_paths("CTT", "CAA", self.table)
        self.assertAlmostEqual(frac_S, 0.5)
        self.assertAlmostEqual(frac_NS, 1.5)
        self.assertEqual(num_paths, 2)

    def test_ng86_intermediate_stop_excluded(self):
        """Test that paths through intermediate stops are excluded."""
        # CAA (Gln) -> TGA (Stop), positions 0 and 1 differ
        # Path 1: CAA -> TAA (Stop) -> TGA (invalid - intermediate stop)
        # Path 2: CAA -> CGA (Arg) -> TGA (valid)
        # Only path 2 is counted
        frac_S, frac_NS, num_paths = _enumerate_ng86_paths("CAA", "TGA", self.table)
        # Should only have 1 valid path (the one avoiding TAA)
        self.assertEqual(num_paths, 1)
        # That path should be: CAA->CGA->TGA (Gln->Arg->Stop: 2 NS steps)
        self.assertEqual(frac_S, 0.0)
        self.assertEqual(frac_NS, 2.0)

    def test_ng86_all_paths_invalid(self):
        """Test codon pair where all paths hit intermediate stops."""
        # Find a codon pair where all paths are invalid
        # This is rare but can happen with certain stop codon transitions
        # Let's manually check: TAA (Stop) -> TGG (Trp), positions 1 and 2
        # Path 1: TAA -> TGA (Stop) -> TGG (intermediate stop)
        # Path 2: TAA -> TAG (Stop) -> TGG (intermediate stop)
        # Both invalid!
        frac_S, frac_NS, num_paths = _enumerate_ng86_paths("TAA", "TGG", self.table)
        self.assertEqual(num_paths, 0)
        self.assertTrue(np.isnan(frac_S))
        self.assertTrue(np.isnan(frac_NS))

    def test_ng86_cache_correctness(self):
        """Test that cache lookups match direct computation."""
        # Test a few random codon pairs
        test_pairs = [
            ("ATG", "TTG"),
            ("CTT", "CAA"),
            ("GGG", "AAA"),
            ("TAA", "TGG"),
        ]

        for anc, der in test_pairs:
            # Direct computation
            direct_S, direct_NS, direct_paths = _enumerate_ng86_paths(
                anc, der, self.table
            )
            # Cache lookup
            cache_S, cache_NS, cache_paths = self.ng86_cache[(anc, der)]

            # Should match exactly
            if np.isnan(direct_S):
                self.assertTrue(np.isnan(cache_S))
            else:
                self.assertEqual(direct_S, cache_S)

            if np.isnan(direct_NS):
                self.assertTrue(np.isnan(cache_NS))
            else:
                self.assertEqual(direct_NS, cache_NS)

            self.assertEqual(direct_paths, cache_paths)

    def test_ng86_cache_size(self):
        """Test that cache contains all 4096 codon pairs."""
        self.assertEqual(len(self.ng86_cache), 64 * 64)


class TestCodonSubstitutionAnalysis(unittest.TestCase):
    """
    Tests for analyze_codon_substitutions_with_ng86_paths().

    Tests the updated analysis function that:
    - Groups sites by codon
    - Returns codon-level events (not per-site)
    - Includes fractional S/N counts
    - Handles k=1, k=2, k=3 cases
    """

    def setUp(self):
        """Set up test data for codon-level analysis."""
        self.ng86_cache = _get_ng86_cache()

        # Simple test case: one gene with a k=2 codon change
        self.genes = [
            {
                "gene_id": "test_gene",
                "start": 100,
                "end": 111,
                "strand": 1,
                "sequence": "ATGCGTACGTTT",  # ATG CGT ACG TTT
            }
        ]

        # Sites at positions 103 and 104 (within second codon CGT)
        self.sites = [
            {
                "contig": "contig1",
                "position": 103,  # Position 1 in second codon
                "gene_id": "test_gene",
                "forward_ref_base": "G",
                "forward_derived_base": "A",
            },
            {
                "contig": "contig1",
                "position": 104,  # Position 2 in second codon (same codon as 103)
                "gene_id": "test_gene",
                "forward_ref_base": "T",
                "forward_derived_base": "C",
            },
        ]

        self.prodigal_records = build_prodigal_records(self.genes)
        self.sig_df, self.pre_df, self.post_df = build_dataframes_from_sites(self.sites)

    def test_codon_grouping(self):
        """Test that sites are correctly grouped into codon events."""
        # Reconstruct ancestral sequences
        profile_by_gene = self.pre_df.groupby("gene_id")
        unique_genes = ["test_gene"]
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            unique_genes, self.prodigal_records, profile_by_gene
        )

        # Calculate potential sites
        for gene_id, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            self.prodigal_records[gene_id]["potential_S_sites"] = pot["S"]
            self.prodigal_records[gene_id]["potential_NS_sites"] = pot["NS"]

        # Analyze codons
        results_df = analyze_codon_substitutions_with_ng86_paths(
            self.sig_df,
            ancestral_seqs,
            self.prodigal_records,
            self.post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have exactly 1 codon event (k=2)
        self.assertEqual(len(results_df), 1)

        # Check k value
        self.assertEqual(results_df.iloc[0]["k"], 2)

        # Check that it has fractional counts
        self.assertIn("frac_S", results_df.columns)
        self.assertIn("frac_NS", results_df.columns)
        self.assertIn("num_valid_paths", results_df.columns)

    def test_codon_event_columns(self):
        """Test that output DataFrame has all required NG86 columns."""
        profile_by_gene = self.pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["test_gene"], self.prodigal_records, profile_by_gene
        )

        for gene_id, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            self.prodigal_records[gene_id]["potential_S_sites"] = pot["S"]
            self.prodigal_records[gene_id]["potential_NS_sites"] = pot["NS"]

        results_df = analyze_codon_substitutions_with_ng86_paths(
            self.sig_df,
            ancestral_seqs,
            self.prodigal_records,
            self.post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Check for new NG86-specific columns
        required_cols = [
            "mag_id",
            "gene_id",
            "contig",
            "codon_start_index",
            "ancestral_codon",
            "derived_codon",
            "ancestral_aa",
            "derived_aa",
            "mutation_type",
            "k",
            "frac_S",
            "frac_NS",
            "num_valid_paths",
            "is_valid",
            "potential_S_codon",
            "potential_NS_codon",
        ]

        for col in required_cols:
            self.assertIn(col, results_df.columns, f"Missing column: {col}")


class TestGlobalDndsCalculation(unittest.TestCase):
    """
    Unit tests for calculate_global_dnds_for_sites() with NG86.

    Tests the updated global dN/dS calculation that:
    - Works with codon-level events (not per-site)
    - Uses fractional S/NS counts
    - Avoids double-counting codons
    - Reports k distribution
    """

    def setUp(self):
        """Set up mock codon events for global dN/dS testing."""
        # Create mock codon events DataFrame
        self.codon_events = pd.DataFrame(
            [
                {
                    "mag_id": "MAG_001",
                    "gene_id": "gene1",
                    "codon_start_index": 0,
                    "k": 1,
                    "frac_S": 0.0,
                    "frac_NS": 1.0,
                    "is_valid": True,
                    "potential_S_codon": 0.0,
                    "potential_NS_codon": 3.0,
                },
                {
                    "mag_id": "MAG_001",
                    "gene_id": "gene1",
                    "codon_start_index": 3,
                    "k": 1,
                    "frac_S": 1.0,
                    "frac_NS": 0.0,
                    "is_valid": True,
                    "potential_S_codon": 1.0,
                    "potential_NS_codon": 2.0,
                },
                {
                    "mag_id": "MAG_001",
                    "gene_id": "gene1",
                    "codon_start_index": 6,
                    "k": 2,
                    "frac_S": 0.5,
                    "frac_NS": 1.5,
                    "is_valid": True,
                    "potential_S_codon": 1.0,
                    "potential_NS_codon": 2.0,
                },
            ]
        )

    def test_global_dnds_with_fractional_counts(self):
        """Test global dN/dS calculation with fractional counts."""
        summary_df = calculate_global_dnds_for_sites(self.codon_events)

        self.assertIsNotNone(summary_df)
        metrics = summary_df.set_index("Metric")["Value"].to_dict()

        # Check metrics exist
        self.assertIn("Total Unique Codons Analyzed", metrics)
        self.assertIn("Observed Non-Synonymous Substitutions (NS)", metrics)
        self.assertIn("Observed Synonymous Substitutions (S)", metrics)

        # Observed counts should be sums of fractional values
        # NS: 1.0 + 0.0 + 1.5 = 2.5
        # S: 0.0 + 1.0 + 0.5 = 1.5
        self.assertAlmostEqual(
            float(metrics["Observed Non-Synonymous Substitutions (NS)"]), 2.5, places=3
        )
        self.assertAlmostEqual(
            float(metrics["Observed Synonymous Substitutions (S)"]), 1.5, places=3
        )

        # Potential sites: sum of unique codons
        # NS: 3.0 + 2.0 + 2.0 = 7.0
        # S: 0.0 + 1.0 + 1.0 = 2.0
        self.assertEqual(
            float(metrics["Total Potential Non-Synonymous Sites (ns)"]), 7.0
        )
        self.assertEqual(float(metrics["Total Potential Synonymous Sites (s)"]), 2.0)

    def test_global_dnds_excludes_invalid_codons(self):
        """Test that invalid codons are excluded from calculation."""
        # Add an invalid codon
        invalid_codon = pd.DataFrame(
            [
                {
                    "mag_id": "MAG_001",
                    "gene_id": "gene1",
                    "codon_start_index": 9,
                    "k": 2,
                    "frac_S": np.nan,
                    "frac_NS": np.nan,
                    "is_valid": False,
                    "potential_S_codon": 1.0,
                    "potential_NS_codon": 2.0,
                }
            ]
        )

        combined = pd.concat([self.codon_events, invalid_codon], ignore_index=True)
        summary_df = calculate_global_dnds_for_sites(combined)

        metrics = summary_df.set_index("Metric")["Value"].to_dict()

        # Should only count 3 valid codons
        self.assertEqual(int(metrics["Total Unique Codons Analyzed"]), 3)
        self.assertEqual(int(metrics["Invalid Codons (all paths hit stops)"]), 1)

    def test_global_dnds_empty_input(self):
        """Test that empty DataFrame returns None."""
        empty_df = pd.DataFrame()
        result = calculate_global_dnds_for_sites(empty_df)
        self.assertIsNone(result)


class TestNG86ComputedValues(unittest.TestCase):
    """
    Comprehensive validation of NG86 computed values against manually
    calculated expectations for diverse codon change scenarios.

    This test class replaces basic structural validation with explicit
    value assertions that verify the NG86 implementation produces correct
    results for all types of codon transitions.

    Tests cover:
    - All codon positions (0, 1, 2) for k=1 changes
    - k=2 and k=3 changes with path averaging
    - Synonymous-only and non-synonymous-only changes
    - Mixed pathway classifications
    - Intermediate stop codon exclusion
    - Edge cases (zero S/NS sites, stop codons)
    - Global dN/dS aggregation with fractional counts
    """

    @classmethod
    def setUpClass(cls):
        """Build NG86 cache once for all tests."""
        cls.ng86_cache = _get_ng86_cache()
        cls.table = CodonTable.unambiguous_dna_by_id[11]

    def test_single_position_changes_all_positions(self):
        """
        Test k=1 changes at codon positions 0, 1, and 2.
        Validates S/NS classification and potential sites for each position.

        MOCK DATA DESIGN:
        - Gene with 3 codons: ATG (Met), TTG (Leu), GCT (Ala)
        - Sites at positions 0, 1, 2 of different codons to test all three positions

        EXPECTED CALCULATIONS:

        Position 0 change: ATG → TTG (Met → Leu)
          - k=1, position 0 changes A→T
          - ATG potential: S=0.0, NS=3.0
          - Translation: Met → Leu (non-synonymous)
          - Expected: frac_S=0.0, frac_NS=1.0, num_valid_paths=1

        Position 1 change: TTG → TAG (Leu → Stop)
          - k=1, position 1 changes T→A
          - TTG potential: S=0.667, NS=2.333
          - Translation: Leu → Stop (non-synonymous)
          - Expected: frac_S=0.0, frac_NS=1.0, num_valid_paths=1

        Position 2 change: GCT → GCC (Ala → Ala)
          - k=1, position 2 changes T→C
          - GCT potential: S=1.0, NS=2.0
          - Translation: Ala → Ala (synonymous)
          - Expected: frac_S=1.0, frac_NS=0.0, num_valid_paths=1

        GLOBAL dN/dS:
          Total potential NS = 3.0 + 2.333 + 2.0 = 7.333
          Total potential S = 0.0 + 0.667 + 1.0 = 1.667
          Observed NS = 1.0 + 1.0 + 0.0 = 2.0
          Observed S = 0.0 + 0.0 + 1.0 = 1.0
          pN = 2.0 / 7.333 = 0.2727
          pS = 1.0 / 1.667 = 0.6
          dN/dS = 0.2727 / 0.6 = 0.4545
        """
        # Create gene with ATGTTGGCT (3 codons)
        # Prodigal coordinates are 1-based: start=101, end=109
        # Contig positions are 0-based: 100-108
        genes = [
            {
                "gene_id": "test_gene",
                "start": 101,
                "end": 109,
                "strand": 1,
                "sequence": "ATGTTGGCT",  # Codons: ATG (100-102), TTG (103-105), GCT (106-108)
            }
        ]

        # Three sites, one per codon position
        sites = [
            {
                "contig": "contig1",
                "position": 100,  # Pos 0 in first codon ATG
                "gene_id": "test_gene",
                "forward_ref_base": "A",  # First base of ATG
                "forward_derived_base": "T",
            },
            {
                "contig": "contig1",
                "position": 104,  # Pos 1 in second codon TTG
                "gene_id": "test_gene",
                "forward_ref_base": "T",  # Second base of TTG
                "forward_derived_base": "A",
            },
            {
                "contig": "contig1",
                "position": 108,  # Pos 2 in third codon GCT
                "gene_id": "test_gene",
                "forward_ref_base": "T",  # Third base of GCT
                "forward_derived_base": "C",
            },
        ]

        # Build data structures
        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        # Reconstruct and analyze
        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["test_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have 3 codon events (one per site, all k=1)
        self.assertEqual(len(codon_events_df), 3)

        # All should be k=1
        self.assertTrue((codon_events_df["k"] == 1).all())

        # Find each specific event and validate
        # Event 1: ATG → TTG (position 0)
        event1 = codon_events_df[codon_events_df["codon_start_index"] == 0].iloc[0]
        self.assertEqual(event1["ancestral_codon"], "ATG")
        self.assertEqual(event1["derived_codon"], "TTG")
        self.assertEqual(event1["mutation_type"], "NS")
        self.assertAlmostEqual(event1["frac_S"], 0.0, places=4)
        self.assertAlmostEqual(event1["frac_NS"], 1.0, places=4)
        self.assertEqual(event1["num_valid_paths"], 1)
        self.assertAlmostEqual(event1["potential_S_codon"], 0.0, places=4)
        self.assertAlmostEqual(event1["potential_NS_codon"], 3.0, places=4)

        # Event 2: TTG → TAG (position 1)
        event2 = codon_events_df[codon_events_df["codon_start_index"] == 3].iloc[0]
        self.assertEqual(event2["ancestral_codon"], "TTG")
        self.assertEqual(event2["derived_codon"], "TAG")
        self.assertEqual(event2["mutation_type"], "NS")
        self.assertAlmostEqual(event2["frac_S"], 0.0, places=4)
        self.assertAlmostEqual(event2["frac_NS"], 1.0, places=4)
        self.assertAlmostEqual(event2["potential_S_codon"], 0.667, places=2)
        self.assertAlmostEqual(event2["potential_NS_codon"], 2.333, places=2)

        # Event 3: GCT → GCC (position 2)
        event3 = codon_events_df[codon_events_df["codon_start_index"] == 6].iloc[0]
        self.assertEqual(event3["ancestral_codon"], "GCT")
        self.assertEqual(event3["derived_codon"], "GCC")
        self.assertEqual(event3["mutation_type"], "S")
        self.assertAlmostEqual(event3["frac_S"], 1.0, places=4)
        self.assertAlmostEqual(event3["frac_NS"], 0.0, places=4)
        self.assertAlmostEqual(event3["potential_S_codon"], 1.0, places=4)
        self.assertAlmostEqual(event3["potential_NS_codon"], 2.0, places=4)

        # Global dN/dS calculation
        summary_df = calculate_global_dnds_for_sites(codon_events_df)
        metrics = summary_df.set_index("Metric")["Value"].to_dict()

        # Total potential sites = sum of unique codons
        # NS: 3.0 (ATG) + 2.333 (TTG) + 2.0 (GCT) = 7.333
        # S: 0.0 (ATG) + 0.667 (TTG) + 1.0 (GCT) = 1.667
        self.assertAlmostEqual(
            float(metrics["Total Potential Non-Synonymous Sites (ns)"]), 7.333, places=2
        )
        self.assertAlmostEqual(
            float(metrics["Total Potential Synonymous Sites (s)"]), 1.667, places=2
        )

        # Observed counts = sum of fractional values
        self.assertAlmostEqual(
            float(metrics["Observed Non-Synonymous Substitutions (NS)"]), 2.0, places=3
        )
        self.assertAlmostEqual(
            float(metrics["Observed Synonymous Substitutions (S)"]), 1.0, places=3
        )

        # pN = 2.0/7.333 = 0.2727, pS = 1.0/1.667 = 0.6
        self.assertAlmostEqual(float(metrics["pNS (NS / ns)"]), 0.2727, places=3)
        self.assertAlmostEqual(float(metrics["pS (S / s)"]), 0.6, places=3)

        # dN/dS = 0.2727/0.6 = 0.4545
        self.assertAlmostEqual(
            float(metrics["Global dN/dS (pNS/pS)"]), 0.4545, places=3
        )

    def test_synonymous_only_transitions(self):
        """
        Test codons where all changes are synonymous (frac_S=1, frac_NS=0).

        MOCK DATA DESIGN:
        - Gene with synonymous-only changes
        - CTT → CTC (Leu → Leu), CGC → CGT (Arg → Arg), GCT → GCC (Ala → Ala)

        EXPECTED CALCULATIONS:
        Each is k=1 synonymous:
          - CTT → CTC: S=1.0, NS=0.0
          - CGC → CGT: S=1.0, NS=0.0
          - GCT → GCC: S=1.0, NS=0.0

        CTT potential: S=1.0, NS=2.0
        CGC potential: S=1.0, NS=2.0
        GCT potential: S=1.0, NS=2.0

        GLOBAL dN/dS:
          Total potential NS = 2.0 + 2.0 + 2.0 = 6.0
          Total potential S = 1.0 + 1.0 + 1.0 = 3.0
          Observed NS = 0.0
          Observed S = 3.0
          pN = 0.0 / 6.0 = 0.0
          pS = 3.0 / 3.0 = 1.0
          dN/dS = 0.0 / 1.0 = 0.0
        """
        genes = [
            {
                "gene_id": "syn_gene",
                "start": 101,  # 1-based
                "end": 109,  # 1-based
                "strand": 1,
                "sequence": "CTTCGCGCT",  # Codons: CTT (100-102), CGC (103-105), GCT (106-108)
            }
        ]

        sites = [
            {
                "contig": "contig1",
                "position": 102,  # Third base of CTT
                "gene_id": "syn_gene",
                "forward_ref_base": "T",  # Matches sequence[2]
                "forward_derived_base": "C",
            },
            {
                "contig": "contig1",
                "position": 105,  # Third base of CGC
                "gene_id": "syn_gene",
                "forward_ref_base": "C",  # Matches sequence[5]
                "forward_derived_base": "T",
            },
            {
                "contig": "contig1",
                "position": 108,  # Third base of GCT
                "gene_id": "syn_gene",
                "forward_ref_base": "T",  # Matches sequence[8]
                "forward_derived_base": "C",
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["syn_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # All events should be synonymous
        self.assertTrue((codon_events_df["mutation_type"] == "S").all())

        # All frac_S should be 1.0, frac_NS should be 0.0
        self.assertTrue((codon_events_df["frac_S"] == 1.0).all())
        self.assertTrue((codon_events_df["frac_NS"] == 0.0).all())

        # Global dN/dS
        summary_df = calculate_global_dnds_for_sites(codon_events_df)
        metrics = summary_df.set_index("Metric")["Value"].to_dict()

        self.assertAlmostEqual(
            float(metrics["Observed Synonymous Substitutions (S)"]), 3.0, places=3
        )
        self.assertAlmostEqual(
            float(metrics["Observed Non-Synonymous Substitutions (NS)"]), 0.0, places=3
        )
        self.assertAlmostEqual(float(metrics["pNS (NS / ns)"]), 0.0, places=4)
        self.assertAlmostEqual(float(metrics["Global dN/dS (pNS/pS)"]), 0.0, places=4)

    def test_nonsynonymous_only_transitions(self):
        """
        Test codons where all changes are non-synonymous.

        MOCK DATA DESIGN:
        - Gene: ATGTTGAAG (Met, Leu, Lys)
        - Changes: ATG→GTG (Met→Val), TTG→ATG (Leu→Met), AAG→TAG (Lys→Stop)

        EXPECTED CALCULATIONS:
        All k=1 non-synonymous:
          - ATG → GTG: Met → Val (NS), frac_S=0.0, frac_NS=1.0
          - TTG → ATG: Leu → Met (NS), frac_S=0.0, frac_NS=1.0
          - AAG → TAG: Lys → Stop (NS), frac_S=0.0, frac_NS=1.0

        Potentials (NCBI table 11):
          - ATG: S=0.0, NS=3.0
          - TTG: S≈0.667, NS≈2.333
          - AAG: S≈0.333, NS≈2.667
          - Totals: S=1.0, NS=8.0

        Global ratios:
          - Observed NS = 3.0, Observed S = 0.0
          - pN = 3.0 / 8.0 = 0.375
          - pS = 0.0 / 1.0 = 0.0
          - dN/dS = undefined (pS=0, expect NaN)
        """
        genes = [
            {
                "gene_id": "ns_gene",
                "start": 101,  # 1-based
                "end": 109,  # 1-based
                "strand": 1,
                "sequence": "ATGTTGAAG",
            }
        ]

        sites = [
            {
                "contig": "contig1",
                "position": 100,
                "gene_id": "ns_gene",
                "forward_ref_base": "A",
                "forward_derived_base": "G",
            },
            {
                "contig": "contig1",
                "position": 103,
                "gene_id": "ns_gene",
                "forward_ref_base": "T",
                "forward_derived_base": "A",
            },
            {
                "contig": "contig1",
                "position": 106,
                "gene_id": "ns_gene",
                "forward_ref_base": "A",
                "forward_derived_base": "T",
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["ns_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # All events should be non-synonymous
        self.assertTrue((codon_events_df["mutation_type"] == "NS").all())

        # All frac_NS should be 1.0, frac_S should be 0.0
        self.assertTrue((codon_events_df["frac_NS"] == 1.0).all())
        self.assertTrue((codon_events_df["frac_S"] == 0.0).all())

        # Global dN/dS
        summary_df = calculate_global_dnds_for_sites(codon_events_df)
        metrics = summary_df.set_index("Metric")["Value"].to_dict()

        self.assertAlmostEqual(
            float(metrics["Observed Non-Synonymous Substitutions (NS)"]), 3.0, places=3
        )
        self.assertAlmostEqual(
            float(metrics["Observed Synonymous Substitutions (S)"]), 0.0, places=3
        )

        # pS should be 0, making dN/dS undefined (NaN)
        self.assertAlmostEqual(float(metrics["pS (S / s)"]), 0.0, places=4)

    def test_two_position_changes_path_averaging(self):
        """
        Test k=2 with explicit 2-path averaging.

        MOCK DATA DESIGN:
        - Gene: ATGCGT (ATG, CGT)
        - Change ATG → TTT (positions 0 and 2 change: A→T, G→T)

        EXPECTED NG86 PATH AVERAGING:

        Changed positions: 0 and 2
        k = 2, so 2 factorial = 2 pathways

        Pathway 1 (change position 0 first, then position 2):
          Step 1: ATG (Met) → TTG (Leu) [position 0: A→T]
                  Met ≠ Leu → NS
          Step 2: TTG (Leu) → TTT (Phe) [position 2: G→T]
                  Leu ≠ Phe → NS
          Pathway 1 total: 0 S steps, 2 NS steps

        Pathway 2 (change position 2 first, then position 0):
          Step 1: ATG (Met) → ATT (Ile) [position 2: G→T]
                  Met ≠ Ile → NS
          Step 2: ATT (Ile) → TTT (Phe) [position 0: A→T]
                  Ile ≠ Phe → NS
          Pathway 2 total: 0 S steps, 2 NS steps

        Average across pathways:
          frac_S = (0 + 0) / 2 = 0.0
          frac_NS = (2 + 2) / 2 = 2.0
          num_valid_paths = 2

        ATG potential sites: S=0.0, NS=3.0

        GLOBAL dN/dS:
          Total potential NS = 3.0
          Total potential S = 0.0
          Observed NS = 2.0
          Observed S = 0.0
          pN = 2.0 / 3.0 = 0.6667
          pS = 0.0 / 0.0 = NaN
          dN/dS = NaN
        """
        genes = [
            {
                "gene_id": "k2_gene",
                "start": 101,  # 1-based: covers contig positions 100-105
                "end": 106,
                "strand": 1,
                "sequence": "ATGCGT",
            }
        ]

        # Two sites in the same codon (ATG)
        sites = [
            {
                "contig": "contig1",
                "position": 100,  # Pos 0 in ATG
                "gene_id": "k2_gene",
                "forward_ref_base": "A",  # sequence[0]
                "forward_derived_base": "T",
            },
            {
                "contig": "contig1",
                "position": 102,  # Pos 2 in ATG
                "gene_id": "k2_gene",
                "forward_ref_base": "G",  # sequence[2]
                "forward_derived_base": "T",
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["k2_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have exactly 1 codon event (k=2)
        self.assertEqual(len(codon_events_df), 1)

        event = codon_events_df.iloc[0]
        self.assertEqual(event["k"], 2)
        self.assertEqual(event["ancestral_codon"], "ATG")
        self.assertEqual(event["derived_codon"], "TTT")
        self.assertEqual(event["mutation_type"], "NS")

        # Validate path averaging results
        self.assertAlmostEqual(event["frac_S"], 0.0, places=4)
        self.assertAlmostEqual(event["frac_NS"], 2.0, places=4)
        self.assertEqual(event["num_valid_paths"], 2)

    def test_three_position_changes_path_averaging(self):
        """
        Test k=3 with 6-pathway averaging (3 factorial = 6 permutations).

        MOCK DATA DESIGN:
        - Gene: ATGCGT (ATG, CGT)
        - Change ATG → CCC (all three positions change: A→C, T→C, G→C)

        EXPECTED NG86 PATH AVERAGING:

        Changed positions: 0, 1, 2
        k = 3, so 3 factorial = 6 pathways

        All 6 permutations:
        --- Paths Analysis ---

        Path 1 (Order: 0, 1, 2)
        1. ATG (Met) → CTG (Leu) : 1 NS
        2. CTG (Leu) → CCG (Pro) : 1 NS
        3. CCG (Pro) → CCC (Pro) : 1 S
        Total: 1 S, 2 NS

        Path 2 (Order: 0, 2, 1)
        1. ATG (Met) → CTG (Leu) : 1 NS
        2. CTG (Leu) → CTC (Leu) : 1 S
        3. CTC (Leu) → CCC (Pro) : 1 NS
        Total: 1 S, 2 NS

        Path 3 (Order: 1, 0, 2)
        1. ATG (Met) → ACG (Thr) : 1 NS
        2. ACG (Thr) → CCG (Pro) : 1 NS
        3. CCG (Pro) → CCC (Pro) : 1 S
        Total: 1 S, 2 NS

        Path 4 (Order: 1, 2, 0)
        1. ATG (Met) → ACG (Thr) : 1 NS
        2. ACG (Thr) → ACC (Thr) : 1 S
        3. ACC (Thr) → CCC (Pro) : 1 NS
        Total: 1 S, 2 NS

        Path 5 (Order: 2, 0, 1)
        1. ATG (Met) → ATC (Ile) : 1 NS
        2. ATC (Ile) → CTC (Leu) : 1 NS
        3. CTC (Leu) → CCC (Pro) : 1 NS
        Total: 0 S, 3 NS

        Path 6 (Order: 2, 1, 0)
        1. ATG (Met) → ATC (Ile) : 1 NS
        2. ATC (Ile) → ACC (Thr) : 1 NS
        3. ACC (Thr) → CCC (Pro) : 1 NS
        Total: 0 S, 3 NS

        --- Final Calculation ---

        Total S steps across 6 valid paths: 1 + 1 + 1 + 1 + 0 + 0 = 4
        Total NS steps across 6 valid paths: 2 + 2 + 2 + 2 + 3 + 3 = 14

        Average frac_S = 4 / 6 = 0.6667
        Average frac_N = 14 / 6 = 2.3333

        ATG potential sites: S=0.0, NS=3.0
        """
        genes = [
            {
                "gene_id": "k3_gene",
                "start": 101,  # 1-based: covers contig positions 100-105
                "end": 106,
                "strand": 1,
                "sequence": "ATGCGT",
            }
        ]

        # Three sites in the same codon (ATG)
        sites = [
            {
                "contig": "contig1",
                "position": 100,  # Pos 0
                "gene_id": "k3_gene",
                "forward_ref_base": "A",  # sequence[0]
                "forward_derived_base": "C",
            },
            {
                "contig": "contig1",
                "position": 101,  # Pos 1
                "gene_id": "k3_gene",
                "forward_ref_base": "T",  # sequence[1]
                "forward_derived_base": "C",
            },
            {
                "contig": "contig1",
                "position": 102,  # Pos 2
                "gene_id": "k3_gene",
                "forward_ref_base": "G",  # sequence[2]
                "forward_derived_base": "C",
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["k3_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have exactly 1 codon event (k=3)
        self.assertEqual(len(codon_events_df), 1)

        event = codon_events_df.iloc[0]
        self.assertEqual(event["k"], 3)
        self.assertEqual(event["ancestral_codon"], "ATG")
        self.assertEqual(event["derived_codon"], "CCC")

        # Validate codon_position and gene_position
        # Gene starts at Prodigal position 101 (1-based) = contig position 100 (0-based)
        # Codon 0 (ATG) spans gene positions 0-2 (contig positions 100-102)
        # Since all three positions in codon changed, we expect:
        # - codon_position: "0,1,2" (all positions within the codon changed)
        # - gene_position: "0,1,2" (gene positions 0, 1, 2 changed)
        self.assertEqual(
            event["codon_position"], "0,1,2"
        )  # All positions in codon changed
        self.assertEqual(event["gene_position"], "0,1,2")  # Gene positions 0, 1, 2

        # Also validate the contig_position field
        self.assertEqual(event["contig_position"], "100,101,102")  # Contig positions

        # Validate ancestral and derived alleles (forward strand)
        self.assertEqual(event["ancestral_allele"], "A,T,G")  # Original bases
        self.assertEqual(event["derived_allele"], "C,C,C")  # Changed bases

        # Validate path averaging results (Met → Pro with some synonymous intermediate steps)
        self.assertAlmostEqual(event["frac_S"], 0.667, places=2)
        self.assertAlmostEqual(event["frac_NS"], 2.333, places=2)
        self.assertEqual(event["num_valid_paths"], 6)

    def test_mixed_pathway_classifications(self):
        """
        Test k=2 where pathway 1 and pathway 2 have different S/NS counts.

        MOCK DATA DESIGN:
        - Gene: CTTAAA (CTT, AAA)
        - Change CTT → CAA (positions 1 and 2: T→A, T→A)

        EXPECTED NG86 PATH AVERAGING:

        Changed positions: 1 and 2
        k = 2, so 2 pathways

        Pathway 1 (change position 1 first, then position 2):
          Step 1: CTT (Leu) → CAT (His) [position 1: T→A]
                  Leu ≠ His → NS
          Step 2: CAT (His) → CAA (Gln) [position 2: T→A]
                  His ≠ Gln → NS
          Pathway 1 total: 0 S steps, 2 NS steps

        Pathway 2 (change position 2 first, then position 1):
          Step 1: CTT (Leu) → CTA (Leu) [position 2: T→A]
                  Leu = Leu → S
          Step 2: CTA (Leu) → CAA (Gln) [position 1: T→A]
                  Leu ≠ Gln → NS
          Pathway 2 total: 1 S step, 1 NS step

        Average across pathways:
          frac_S = (0 + 1) / 2 = 0.5
          frac_NS = (2 + 1) / 2 = 1.5
          num_valid_paths = 2

        CTT potential sites: S=1.0, NS=2.0
        """
        genes = [
            {
                "gene_id": "mixed_gene",
                "start": 101,  # 1-based: covers contig positions 100-105
                "end": 106,
                "strand": 1,
                "sequence": "CTTAAA",
            }
        ]

        # Two sites in the same codon (CTT)
        sites = [
            {
                "contig": "contig1",
                "position": 101,  # Pos 1 in CTT
                "gene_id": "mixed_gene",
                "forward_ref_base": "T",  # sequence[1]
                "forward_derived_base": "A",
            },
            {
                "contig": "contig1",
                "position": 102,  # Pos 2 in CTT
                "gene_id": "mixed_gene",
                "forward_ref_base": "T",  # sequence[2]
                "forward_derived_base": "A",
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["mixed_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have exactly 1 codon event (k=2)
        self.assertEqual(len(codon_events_df), 1)

        event = codon_events_df.iloc[0]
        self.assertEqual(event["k"], 2)
        self.assertEqual(event["ancestral_codon"], "CTT")
        self.assertEqual(event["derived_codon"], "CAA")

        # Validate mixed pathway averaging
        self.assertAlmostEqual(event["frac_S"], 0.5, places=4)
        self.assertAlmostEqual(event["frac_NS"], 1.5, places=4)
        self.assertEqual(event["num_valid_paths"], 2)

    def test_intermediate_stop_codon_exclusion(self):
        """
        Test that pathways hitting intermediate stop codons are excluded.

        MOCK DATA DESIGN:
        - Gene: CAACGT (CAA, CGT)
        - Change CAA → TGA (positions 0 and 1: C→T, A→G)

        EXPECTED NG86 PATH AVERAGING:

        Changed positions: 0 and 1
        k = 2, so 2 possible pathways

        Pathway 1 (change position 0 first, then position 1):
          Step 1: CAA (Gln) → TAA (Stop) [position 0: C→T]
                  Intermediate stop codon → EXCLUDED
          This pathway is INVALID

        Pathway 2 (change position 1 first, then position 0):
          Step 1: CAA (Gln) → CGA (Arg) [position 1: A→G]
                  Gln ≠ Arg → NS
          Step 2: CGA (Arg) → TGA (Stop) [position 0: C→T]
                  Arg ≠ Stop → NS
          Pathway 2 total: 0 S steps, 2 NS steps
          This pathway is VALID (final stop is allowed)

        Only 1 valid pathway remains:
          frac_S = 0.0
          frac_NS = 2.0
          num_valid_paths = 1

        CAA potential sites: S≈0.333, NS≈2.667
        """
        genes = [
            {
                "gene_id": "stop_gene",
                "start": 101,  # 1-based: covers contig positions 100-105
                "end": 106,
                "strand": 1,
                "sequence": "CAACGT",
            }
        ]

        # Two sites in CAA codon
        sites = [
            {
                "contig": "contig1",
                "position": 100,  # Pos 0 in CAA
                "gene_id": "stop_gene",
                "forward_ref_base": "C",  # sequence[0]
                "forward_derived_base": "T",
            },
            {
                "contig": "contig1",
                "position": 101,  # Pos 1 in CAA
                "gene_id": "stop_gene",
                "forward_ref_base": "A",  # sequence[1]
                "forward_derived_base": "G",
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["stop_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have 1 codon event
        self.assertEqual(len(codon_events_df), 1)

        event = codon_events_df.iloc[0]
        self.assertEqual(event["k"], 2)
        self.assertEqual(event["ancestral_codon"], "CAA")
        self.assertEqual(event["derived_codon"], "TGA")

        # Only 1 valid pathway (the other hit intermediate stop TAA)
        self.assertEqual(event["num_valid_paths"], 1)
        self.assertAlmostEqual(event["frac_S"], 0.0, places=4)
        self.assertAlmostEqual(event["frac_NS"], 2.0, places=4)

    def test_zero_synonymous_sites_codons(self):
        """
        Test codons with zero synonymous sites (ATG, TGG).

        MOCK DATA DESIGN:
        - Gene: ATGTGG (ATG-Met, TGG-Trp)
        - Changes: ATG→GTG, TGG→TGC

        EXPECTED CALCULATIONS:

        ATG (Methionine) potential sites calculation:
          Position 0: A→{T,G,C} = {TTG(Leu), GTG(Val), CTG(Leu)} = all NS
          Position 1: T→{A,G,C} = {AAG(Lys), AGG(Arg), ACG(Thr)} = all NS
          Position 2: G→{A,T,C} = {ATA(Ile), ATT(Ile), ATC(Ile)} = all NS
          Total: S = 0/3 = 0.0, NS = 9/3 = 3.0

        TGG (Tryptophan) potential sites:
          Position 0: T→{A,G,C} = {AGG(Arg), GGG(Gly), CGG(Arg)} = all NS
          Position 1: G→{A,T,C} = {TAG(Stop), TTG(Leu), TCG(Ser)} = all NS
          Position 2: G→{A,T,C} = {TGA(Stop), TGT(Cys), TGC(Cys)} = all NS
          Total: S = 0/3 = 0.0, NS = 9/3 = 3.0

        Both codons have S=0.0, NS=3.0 as expected
        """
        genes = [
            {
                "gene_id": "zero_s_gene",
                "start": 101,  # 1-based: covers contig positions 100-105
                "end": 106,
                "strand": 1,
                "sequence": "ATGTGG",
            }
        ]

        sites = [
            {
                "contig": "contig1",
                "position": 100,  # ATG → GTG (Met → Val)
                "gene_id": "zero_s_gene",
                "forward_ref_base": "A",  # sequence[0]
                "forward_derived_base": "G",
            },
            {
                "contig": "contig1",
                "position": 105,  # TGG → TGC (Trp → Cys)
                "gene_id": "zero_s_gene",
                "forward_ref_base": "G",  # sequence[5]
                "forward_derived_base": "C",
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["zero_s_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        # Verify both codons have zero S sites
        self.assertAlmostEqual(
            prodigal_records["zero_s_gene"]["potential_S_sites"], 0.0, places=4
        )
        self.assertAlmostEqual(
            prodigal_records["zero_s_gene"]["potential_NS_sites"], 6.0, places=4
        )

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have 2 events (both k=1, both NS)
        self.assertEqual(len(codon_events_df), 2)
        self.assertTrue((codon_events_df["k"] == 1).all())
        self.assertTrue((codon_events_df["mutation_type"] == "NS").all())

        # Both should have zero S potential
        for _, event in codon_events_df.iterrows():
            self.assertAlmostEqual(event["frac_S"], 0.0, places=4)
            self.assertAlmostEqual(event["frac_NS"], 1.0, places=4)
            self.assertAlmostEqual(event["potential_S_codon"], 0.0, places=4)
            self.assertAlmostEqual(event["potential_NS_codon"], 3.0, places=4)
            self.assertAlmostEqual(event["potential_S_sites_gene"], 0.0, places=4)
            self.assertAlmostEqual(event["potential_NS_sites_gene"], 6.0, places=4)

        # Global dN/dS
        summary_df = calculate_global_dnds_for_sites(codon_events_df)
        metrics = summary_df.set_index("Metric")["Value"].to_dict()

        # Validate global dN/dS metrics
        # Total potential NS = 6.0 (from both codons)
        # Total potential S = 0.0 (both codons have no synonymous sites)
        # Observed NS = 2.0 (both events are NS)
        # Observed S = 0.0 (no synonymous changes)
        # pNS = 2.0 / 6.0 = 0.333...
        # pS = 0.0 / 0.0 = NaN (division by zero)
        # dN/dS = pNS / pS = NaN (pS is NaN)

        self.assertAlmostEqual(
            float(metrics["Observed Non-Synonymous Substitutions (NS)"]), 2.0, places=3
        )
        self.assertAlmostEqual(
            float(metrics["Observed Synonymous Substitutions (S)"]), 0.0, places=3
        )
        self.assertAlmostEqual(float(metrics["pNS (NS / ns)"]), 0.333, places=3)

        # pS should be NaN (0/0)
        self.assertTrue(pd.isna(float(metrics["pS (S / s)"])))

        # Global dN/dS should be NaN (division by NaN)
        self.assertTrue(pd.isna(float(metrics["Global dN/dS (pNS/pS)"])))

    def test_global_dnds_computed_values(self):
        """
        Test complete dN/dS calculation with known inputs.

        MOCK DATA DESIGN:
        - Gene: ATGCTTACC (ATG-Met, CTT-Leu, ACC-Thr)
        - Change 1: ATG→TTG (position 0: A→T, k=1, NS)
        - Change 2: CTT→CTC (position 2: T→C, k=1, S)
        - Change 3: ACC→ATC (position 1: C→T, k=1, NS)

        EXPECTED POTENTIAL SITES:
        ATG (Met): S=0.0, NS=3.0
        CTT (Leu): S=1.0, NS=2.0
        ACC (Thr): S=1.0, NS=2.0

        EXPECTED OBSERVED COUNTS:
        Event 1 (ATG→TTG): frac_S=0.0, frac_NS=1.0
        Event 2 (CTT→CTC): frac_S=1.0, frac_NS=0.0
        Event 3 (ACC→ATC): frac_S=0.0, frac_NS=1.0

        GLOBAL dN/dS CALCULATION:
        Total potential NS = 3.0 + 2.0 + 2.0 = 7.0
        Total potential S = 0.0 + 1.0 + 1.0 = 2.0
        Observed NS = 1.0 + 0.0 + 1.0 = 2.0
        Observed S = 0.0 + 1.0 + 0.0 = 1.0
        pN = 2.0 / 7.0 = 0.2857
        pS = 1.0 / 2.0 = 0.5
        dN/dS = 0.2857 / 0.5 = 0.5714
        """
        genes = [
            {
                "gene_id": "calc_gene",
                "start": 101,  # 1-based: covers contig positions 100-108
                "end": 109,
                "strand": 1,
                "sequence": "ATGCTTACC",
            }
        ]

        sites = [
            {
                "contig": "contig1",
                "position": 100,  # ATG pos 0
                "gene_id": "calc_gene",
                "forward_ref_base": "A",  # sequence[0]
                "forward_derived_base": "T",
            },
            {
                "contig": "contig1",
                "position": 105,  # CTT pos 2
                "gene_id": "calc_gene",
                "forward_ref_base": "T",  # sequence[5]
                "forward_derived_base": "C",
            },
            {
                "contig": "contig1",
                "position": 107,  # ACC pos 1
                "gene_id": "calc_gene",
                "forward_ref_base": "C",  # sequence[7]
                "forward_derived_base": "T",
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["calc_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Validate individual events
        self.assertEqual(len(codon_events_df), 3)

        # Global dN/dS
        summary_df = calculate_global_dnds_for_sites(codon_events_df)
        metrics = summary_df.set_index("Metric")["Value"].to_dict()

        # Validate computed values against manual calculations
        self.assertAlmostEqual(
            float(metrics["Total Potential Non-Synonymous Sites (ns)"]), 7.0, places=2
        )
        self.assertAlmostEqual(
            float(metrics["Total Potential Synonymous Sites (s)"]), 2.0, places=2
        )
        self.assertAlmostEqual(
            float(metrics["Observed Non-Synonymous Substitutions (NS)"]), 2.0, places=3
        )
        self.assertAlmostEqual(
            float(metrics["Observed Synonymous Substitutions (S)"]), 1.0, places=3
        )
        self.assertAlmostEqual(float(metrics["pNS (NS / ns)"]), 0.2857, places=4)
        self.assertAlmostEqual(float(metrics["pS (S / s)"]), 0.5, places=4)
        self.assertAlmostEqual(
            float(metrics["Global dN/dS (pNS/pS)"]), 0.5714, places=4
        )

    def test_gene_and_mag_level_summaries(self):
        """
        Test gene-level and MAG-level dN, dS, dN/dS calculations using exact NG86 values.

        MOCK DATA DESIGN (NCBI table 11):
        - Gene1: ATGCTT (ATG=Met, CTT=Leu)
        Changes:
            - ATG→TTG (pos 0)  -> NS (Met→Leu)
            - CTT→CTC (pos 5)  -> S  (Leu→Leu)
        NG86 per-codon potentials:
            - ATG: S=0.000, N=3.000  (Met is unique; all 9 single changes are NS → N=9/3=3)
            - CTT: S=1.000, N=2.000  (3 syn, 6 nonsyn → S=3/3=1, N=6/3=2)
        Gene1 potentials: S=1.000, N=5.000
        Gene1 observed:  S=1.000, N=1.000
        => dN = 1.0/5.0   = 0.2
            dS = 1.0/1.0   = 1.0
            dN/dS          = 0.2

        - Gene2: GCTAGG (GCT=Ala, AGG=Arg)
        Changes:
            - GCT→GCC (pos 2)  -> S  (Ala→Ala)
            - AGG→ATG (pos 3)  -> NS (Arg→Met)
        NG86 per-codon potentials:
            - GCT: S=1.000, N=2.000  (Ala has 4 codons; 3 syn changes)
            - AGG: S=0.667, N=2.333  (2 syn, 7 nonsyn → S=2/3, N=7/3)
        Gene2 potentials: S=1.667, N=4.333
        Gene2 observed:  S=1.000, N=1.000
        => dN = 1.0/4.333... = 3/13 ≈ 0.230769
            dS = 1.0/1.667... = 3/5   = 0.6
            dN/dS = (3/13)/(3/5) = 5/13 ≈ 0.384615

        MAG-LEVEL (sum over genes):
        Potentials: S = 1.000 + 1.667 = 2.667 = 8/3
                    N = 5.000 + 4.333 = 9.333 = 28/3
        Observed:   S = 2.000
                    N = 2.000
        => dN = 2 / (28/3) = 3/14 ≈ 0.2142857
            dS = 2 / ( 8/3) = 3/4  = 0.75
            dN/dS = (3/14)/(3/4) = 2/7 ≈ 0.285714
        """
        genes = [
            {
                "gene_id": "gene1",
                "start": 101,  # 1-based: covers contig positions 100-105
                "end": 106,
                "strand": 1,
                "sequence": "ATGCTT",
            },
            {
                "gene_id": "gene2",
                "start": 201,  # 1-based: covers contig positions 200-205
                "end": 206,
                "strand": 1,
                "sequence": "GCTAGG",
            },
        ]

        sites = [
            # Gene1 changes
            {
                "contig": "contig1",
                "position": 100,
                "gene_id": "gene1",
                "forward_ref_base": "A",  # sequence[0] in ATG
                "forward_derived_base": "T",  # ATG->TTG (NS)
            },
            {
                "contig": "contig1",
                "position": 105,
                "gene_id": "gene1",
                "forward_ref_base": "T",  # sequence[5] in CTT
                "forward_derived_base": "C",  # CTT->CTC (S)
            },
            # Gene2 changes
            {
                "contig": "contig1",
                "position": 202,
                "gene_id": "gene2",
                "forward_ref_base": "T",  # sequence[2] in GCT
                "forward_derived_base": "C",  # GCT->GCC (S)
            },
            {
                "contig": "contig1",
                "position": 203,
                "gene_id": "gene2",
                "forward_ref_base": "A",  # sequence[3] in AGG
                "forward_derived_base": "T",  # AGG->ATG (NS)
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["gene1", "gene2"], prodigal_records, profile_by_gene
        )

        # Not strictly required by summarize_results (it uses per-codon potentials),
        # but keep parity with the pipeline setup:
        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Quick sanity checks on the codon-level observed fractions
        # We expect 4 single-step (k=1) events: two S, two NS
        self.assertAlmostEqual(codon_events_df["frac_S"].sum(), 2.0, places=6)
        self.assertAlmostEqual(codon_events_df["frac_NS"].sum(), 2.0, places=6)

        # Generate summaries
        summaries = summarize_results(codon_events_df)
        gene_stats = summaries["gene"]
        mag_stats = summaries["mag"]

        # Should have 2 genes
        self.assertEqual(len(gene_stats), 2)

        # ---- Verify gene1 exact values ----
        gene1 = gene_stats[gene_stats["gene_id"] == "gene1"].iloc[0]
        # Potentials per gene
        self.assertAlmostEqual(gene1["potential_S_sites_gene"], 1.0, places=6)
        self.assertAlmostEqual(gene1["potential_NS_sites_gene"], 5.0, places=6)
        # Observed (fractional) per gene
        self.assertAlmostEqual(gene1["total_frac_s"], 1.0, places=6)
        self.assertAlmostEqual(gene1["total_frac_ns"], 1.0, places=6)
        # dN, dS, dN/dS
        self.assertAlmostEqual(gene1["dN"], 0.2, places=6)
        self.assertAlmostEqual(gene1["dS"], 1.0, places=6)
        self.assertAlmostEqual(gene1["dN_dS_ratio"], 0.2, places=6)

        # ---- Verify gene2 exact values ----
        gene2 = gene_stats[gene_stats["gene_id"] == "gene2"].iloc[0]
        self.assertAlmostEqual(gene2["potential_S_sites_gene"], 1.667, places=3)  # 5/3
        self.assertAlmostEqual(
            gene2["potential_NS_sites_gene"], 4.333, places=3
        )  # 13/3
        self.assertAlmostEqual(gene2["total_frac_s"], 1.0, places=6)
        self.assertAlmostEqual(gene2["total_frac_ns"], 1.0, places=6)
        self.assertAlmostEqual(gene2["dN"], 3.0 / 13.0, places=6)  # ≈ 0.230769
        self.assertAlmostEqual(gene2["dS"], 3.0 / 5.0, places=6)  # 0.6
        self.assertAlmostEqual(gene2["dN_dS_ratio"], 5.0 / 13.0, places=6)  # ≈ 0.384615

        # ---- MAG-level should aggregate both genes ----
        self.assertEqual(len(mag_stats), 1)
        mag = mag_stats.iloc[0]
        self.assertEqual(mag["num_genes"], 2)

        # Aggregate potentials
        self.assertAlmostEqual(
            mag["potential_S_sites"], 8.0 / 3.0, places=6
        )  # 2.666666...
        self.assertAlmostEqual(
            mag["potential_NS_sites"], 28.0 / 3.0, places=6
        )  # 9.333333...

        # Aggregate observed (fractional)
        self.assertAlmostEqual(mag["total_frac_s"], 2.0, places=6)
        self.assertAlmostEqual(mag["total_frac_ns"], 2.0, places=6)

        # dN, dS, dN/dS
        self.assertAlmostEqual(
            mag["dN"], (2.0) / (28.0 / 3.0), places=6
        )  # 3/14 ≈ 0.2142857
        self.assertAlmostEqual(mag["dS"], (2.0) / (8.0 / 3.0), places=6)  # 3/4  = 0.75
        self.assertAlmostEqual(mag["dN_dS_ratio"], (2.0 / 7.0), places=6)  # ≈ 0.285714

    def test_negative_strand_single_position_changes(self):
        """
        Test k=1 changes at different codon positions (0, 1, 2) on a negative strand gene.

        CRITICAL NEGATIVE STRAND CONCEPTS:
        - Prodigal sequence is stored as reverse complement (5'→3' coding strand)
        - Genomic positions refer to the forward strand
        - Position mapping: pos_in_gene = (end - 1) - position
        - Alleles from profile (forward strand) must be complemented before applying

        MOCK DATA DESIGN:
        - Gene with strand=-1, sequence="ATGCTTACG" (9 bases, 3 codons)
        - Prodigal: start=201, end=209 (covers contig positions 200-208)
        - Sequence covers genomic positions 200-208 (0-based)

        COORDINATE MAPPING:
        - Genomic position 208 → gene index 0 (first base 'A')
        - Genomic position 207 → gene index 1 (second base 'T')
        - Genomic position 200 → gene index 8 (last base 'G')

        MUTATIONS:
        1. Position 0 change: ATG → TTG (codon 0, pos 0)
           - Genomic position 208 → gene index 0
           - Forward ref='T', derived='A'
           - Complement: ref='A' (in sequence), derived='T' (applied)
           - Result: ATG (Met) → TTG (Leu), NS
           - Expected: frac_S=0.0, frac_NS=1.0

        2. Position 1 change: CTT → CAT (codon 1, pos 1)
           - Genomic position 204 → gene index 4
           - Forward ref='A', derived='T'
           - Complement: ref='T' (in sequence), derived='A' (applied)
           - Result: CTT (Leu) → CAT (His), NS
           - Expected: frac_S=0.0, frac_NS=1.0

        3. Position 2 change: ACG → ACC (codon 2, pos 2)
           - Genomic position 200 → gene index 8
           - Forward ref='C', derived='G'
           - Complement: ref='G' (in sequence), derived='C' (applied)
           - Result: ACG (Thr) → ACC (Thr), S
           - Expected: frac_S=1.0, frac_NS=0.0

        POTENTIAL SITES:
        - ATG: S=0.0, NS=3.0
        - CTT: S=1.0, NS=2.0
        - ACG: S=1.0, NS=2.0
        - Total: S=2.0, NS=7.0

        GLOBAL dN/dS:
        - Observed NS = 2.0, Observed S = 1.0
        - pN = 2.0/7.0 = 0.2857
        - pS = 1.0/2.0 = 0.5
        - dN/dS = 0.2857/0.5 = 0.5714
        """
        genes = [
            {
                "gene_id": "neg_strand_gene",
                "start": 201,  # 1-based
                "end": 209,  # 1-based (covers positions 200-208, 0-based)
                "strand": -1,
                "sequence": "ATGCTTACG",  # Reverse complement of genomic region
            }
        ]

        # Three sites at different codon positions
        sites = [
            {
                "contig": "contig1",
                "position": 208,  # Maps to gene index 0 (codon 0, pos 0)
                "gene_id": "neg_strand_gene",
                "forward_ref_base": "T",  # Complement of sequence[0]='A'
                "forward_derived_base": "A",  # Will complement to 'T', changing A→T
            },
            {
                "contig": "contig1",
                "position": 204,  # Maps to gene index 4 (codon 1, pos 1)
                "gene_id": "neg_strand_gene",
                "forward_ref_base": "A",  # Complement of sequence[4]='T'
                "forward_derived_base": "T",  # Will complement to 'A', changing T→A
            },
            {
                "contig": "contig1",
                "position": 200,  # Maps to gene index 8 (codon 2, pos 2)
                "gene_id": "neg_strand_gene",
                "forward_ref_base": "C",  # Complement of sequence[8]='G'
                "forward_derived_base": "G",  # Will complement to 'C', changing G→C
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["neg_strand_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have 3 codon events (one per site, all k=1)
        self.assertEqual(len(codon_events_df), 3)

        # All should be k=1
        self.assertTrue((codon_events_df["k"] == 1).all())

        # Find each specific event and validate
        # Event 1: ATG → TTG (position 0)
        event1 = codon_events_df[codon_events_df["codon_start_index"] == 0].iloc[0]
        self.assertEqual(event1["ancestral_codon"], "ATG")
        self.assertEqual(event1["derived_codon"], "TTG")
        self.assertEqual(event1["mutation_type"], "NS")
        self.assertAlmostEqual(event1["frac_S"], 0.0, places=4)
        self.assertAlmostEqual(event1["frac_NS"], 1.0, places=4)
        self.assertEqual(event1["num_valid_paths"], 1)
        self.assertEqual(event1["strand"], -1)

        # Event 2: CTT → CAT (position 1)
        event2 = codon_events_df[codon_events_df["codon_start_index"] == 3].iloc[0]
        self.assertEqual(event2["ancestral_codon"], "CTT")
        self.assertEqual(event2["derived_codon"], "CAT")
        self.assertEqual(event2["mutation_type"], "NS")
        self.assertAlmostEqual(event2["frac_S"], 0.0, places=4)
        self.assertAlmostEqual(event2["frac_NS"], 1.0, places=4)
        self.assertEqual(event2["strand"], -1)

        # Event 3: ACG → ACC (position 2)
        event3 = codon_events_df[codon_events_df["codon_start_index"] == 6].iloc[0]
        self.assertEqual(event3["ancestral_codon"], "ACG")
        self.assertEqual(event3["derived_codon"], "ACC")
        self.assertEqual(event3["mutation_type"], "S")
        self.assertAlmostEqual(event3["frac_S"], 1.0, places=4)
        self.assertAlmostEqual(event3["frac_NS"], 0.0, places=4)
        self.assertEqual(event3["strand"], -1)

        # Global dN/dS calculation
        summary_df = calculate_global_dnds_for_sites(codon_events_df)
        metrics = summary_df.set_index("Metric")["Value"].to_dict()

        # Validate global metrics
        self.assertAlmostEqual(
            float(metrics["Observed Non-Synonymous Substitutions (NS)"]), 2.0, places=3
        )
        self.assertAlmostEqual(
            float(metrics["Observed Synonymous Substitutions (S)"]), 1.0, places=3
        )
        self.assertAlmostEqual(float(metrics["pNS (NS / ns)"]), 0.2857, places=3)
        self.assertAlmostEqual(float(metrics["pS (S / s)"]), 0.5, places=3)
        self.assertAlmostEqual(
            float(metrics["Global dN/dS (pNS/pS)"]), 0.5714, places=3
        )

    def test_negative_strand_synonymous_transitions(self):
        """
        Test synonymous-only changes on negative strand gene.

        MOCK DATA DESIGN:
        - Gene: CTTCGCGCT (9 bases, 3 codons: CTT, CGC, GCT)
        - Prodigal: start=201, end=209, strand=-1
        - All changes are synonymous (wobble position)

        COORDINATE MAPPING:
        - Position 208 → gene index 0 (first codon, third position)
        - Position 205 → gene index 3 (second codon, first position) - ERROR, should be pos 2!
        - Position 200 → gene index 8 (third codon, third position)

        Actually, let me recalculate:
        - Genomic positions 200-208 (9 positions)
        - Gene sequence "CTTCGCGCT"
        - Position 208 → index (209-1) - 208 = 0 → 'C'
        - Position 207 → index (209-1) - 207 = 1 → 'T'
        - Position 206 → index (209-1) - 206 = 2 → 'T' (third position of CTT)
        - Position 205 → index (209-1) - 205 = 3 → 'C'
        - Position 204 → index (209-1) - 204 = 4 → 'G'
        - Position 203 → index (209-1) - 203 = 5 → 'C' (third position of CGC)
        - Position 202 → index (209-1) - 202 = 6 → 'G'
        - Position 201 → index (209-1) - 201 = 7 → 'C'
        - Position 200 → index (209-1) - 200 = 8 → 'T' (third position of GCT)

        MUTATIONS:
        1. CTT → CTC (pos 206): Leu → Leu (synonymous)
           - Forward: C→T, Complement: G→A applied
        2. CGC → CGT (pos 203): Arg → Arg (synonymous)
           - Forward: C→T, Complement: G→A applied
        3. GCT → GCC (pos 200): Ala → Ala (synonymous)
           - Forward: T→C, Complement: A→G applied

        Wait, I need to think about this more carefully. The forward_ref_base should match the complement of what's in the sequence.

        For neg strand gene with sequence "CTTCGCGCT":
        - Position 206 maps to gene index 2 ('T', third position of CTT)
        - Forward strand at pos 206 has complement of 'T' = 'A'
        - So forward_ref_base = 'A'
        - To change CTT→CTC, we need T→C in the sequence
        - Forward derived = 'G' (complement of 'C')

        Let me recalculate properly:
        """
        genes = [
            {
                "gene_id": "neg_syn_gene",
                "start": 201,  # 1-based
                "end": 209,  # 1-based
                "strand": -1,
                "sequence": "CTTCGCGCT",  # Codons: CTT, CGC, GCT
            }
        ]

        # Three synonymous changes at third codon positions
        sites = [
            {
                "contig": "contig1",
                "position": 206,  # Maps to gene index 2 (CTT pos 2)
                "gene_id": "neg_syn_gene",
                "forward_ref_base": "A",  # Complement of sequence[2]='T'
                "forward_derived_base": "G",  # Complements to 'C', making CTT→CTC
            },
            {
                "contig": "contig1",
                "position": 203,  # Maps to gene index 5 (CGC pos 2)
                "gene_id": "neg_syn_gene",
                "forward_ref_base": "G",  # Complement of sequence[5]='C'
                "forward_derived_base": "A",  # Complements to 'T', making CGC→CGT
            },
            {
                "contig": "contig1",
                "position": 200,  # Maps to gene index 8 (GCT pos 2)
                "gene_id": "neg_syn_gene",
                "forward_ref_base": "A",  # Complement of sequence[8]='T'
                "forward_derived_base": "G",  # Complements to 'C', making GCT→GCC
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["neg_syn_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # All events should be synonymous
        self.assertTrue((codon_events_df["mutation_type"] == "S").all())
        self.assertTrue((codon_events_df["frac_S"] == 1.0).all())
        self.assertTrue((codon_events_df["frac_NS"] == 0.0).all())

        # Verify strand is -1 for all events
        self.assertTrue((codon_events_df["strand"] == -1).all())

        # Global dN/dS should be 0.0
        summary_df = calculate_global_dnds_for_sites(codon_events_df)
        metrics = summary_df.set_index("Metric")["Value"].to_dict()
        self.assertAlmostEqual(float(metrics["Global dN/dS (pNS/pS)"]), 0.0, places=4)

    def test_negative_strand_nonsynonymous_transitions(self):
        """
        Test non-synonymous-only changes on negative strand gene.

        MOCK DATA DESIGN:
        - Gene: ATGTTGAAG (9 bases, 3 codons: ATG, TTG, AAG)
        - Prodigal: start=201, end=209, strand=-1
        - All changes are non-synonymous

        COORDINATE MAPPING:
        - Position 208 → gene index 0 (ATG pos 0)
        - Position 205 → gene index 3 (TTG pos 0)
        - Position 202 → gene index 6 (AAG pos 0)

        MUTATIONS:
        1. ATG → GTG (pos 208): Met → Val
           - Forward: A→C, Complement: T→G
        2. TTG → ATG (pos 205): Leu → Met
           - Forward: A→T, Complement: T→A
        3. AAG → TAG (pos 202): Lys → Stop
           - Forward: T→A, Complement: A→T (changes codon first base A→T)

        POTENTIAL SITES (NCBI table 11):
          - ATG: S=0.0, NS=3.0
          - TTG: S≈0.667, NS≈2.333
          - AAG: S≈0.333, NS≈2.667
          - Totals: S=1.0, NS=8.0

        GLOBAL dN/dS EXPECTATIONS:
          - Observed NS = 3.0, Observed S = 0.0
          - pN = 3.0 / 8.0 = 0.375
          - pS = 0.0 / 1.0 = 0.0
          - dN/dS = undefined (pS=0, expect NaN)
        """
        genes = [
            {
                "gene_id": "neg_ns_gene",
                "start": 201,
                "end": 209,
                "strand": -1,
                "sequence": "ATGTTGAAG",
            }
        ]

        sites = [
            {
                "contig": "contig1",
                "position": 208,  # ATG pos 0
                "gene_id": "neg_ns_gene",
                "forward_ref_base": "A",  # Complement of 'A' is 'T'
                "forward_derived_base": "C",  # Complements to 'G', making ATG→GTG
            },
            {
                "contig": "contig1",
                "position": 205,  # TTG pos 0
                "gene_id": "neg_ns_gene",
                "forward_ref_base": "A",  # Complement of 'T' is 'A'
                "forward_derived_base": "T",  # Complements to 'A', making TTG→ATG
            },
            {
                "contig": "contig1",
                "position": 202,  # AAG pos 0
                "gene_id": "neg_ns_gene",
                "forward_ref_base": "T",  # Complement of 'A' is 'T'
                "forward_derived_base": "A",  # Complements to 'T', making AAG→TAG
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["neg_ns_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # All events should be non-synonymous
        self.assertTrue((codon_events_df["mutation_type"] == "NS").all())
        self.assertTrue((codon_events_df["frac_NS"] == 1.0).all())
        self.assertTrue((codon_events_df["frac_S"] == 0.0).all())

        # Verify strand is -1
        self.assertTrue((codon_events_df["strand"] == -1).all())

        # Global dN/dS should be NaN (pS=0)
        summary_df = calculate_global_dnds_for_sites(codon_events_df)
        metrics = summary_df.set_index("Metric")["Value"].to_dict()
        self.assertAlmostEqual(float(metrics["pS (S / s)"]), 0.0, places=4)

    def test_negative_strand_two_position_changes(self):
        """
        Test k=2 with path averaging on negative strand gene.

        MOCK DATA DESIGN:
        - Gene: ATGCGT (6 bases, 2 codons: ATG, CGT)
        - Prodigal: start=201, end=206, strand=-1
        - Change ATG → TTT (positions 0 and 2)

        COORDINATE MAPPING:
        - Position 205 → gene index 0 (ATG pos 0)
        - Position 203 → gene index 2 (ATG pos 2)

        MUTATIONS:
        - Genomic pos 205: Forward T→A, Complement: A→T (A→T in sequence)
        - Genomic pos 203: Forward C→A, Complement: G→T (G→T in sequence)
        - Result: ATG → TTT

        NG86 PATH AVERAGING:
        Path 1: ATG → TTG → TTT (Met→Leu→Phe: 2 NS)
        Path 2: ATG → ATT → TTT (Met→Ile→Phe: 2 NS)
        Average: frac_S=0.0, frac_NS=2.0, num_valid_paths=2

        POTENTIAL SITES:
        - ATG: S=0.0, NS=3.0
        """
        genes = [
            {
                "gene_id": "neg_k2_gene",
                "start": 201,
                "end": 206,
                "strand": -1,
                "sequence": "ATGCGT",
            }
        ]

        sites = [
            {
                "contig": "contig1",
                "position": 205,  # Maps to gene index 0 (ATG pos 0)
                "gene_id": "neg_k2_gene",
                "forward_ref_base": "T",  # Complement of sequence[0]='A'
                "forward_derived_base": "A",  # Complements to 'T'
            },
            {
                "contig": "contig1",
                "position": 203,  # Maps to gene index 2 (ATG pos 2)
                "gene_id": "neg_k2_gene",
                "forward_ref_base": "C",  # Complement of sequence[2]='G'
                "forward_derived_base": "A",  # Complements to 'T'
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["neg_k2_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have exactly 1 codon event (k=2)
        self.assertEqual(len(codon_events_df), 1)

        event = codon_events_df.iloc[0]
        self.assertEqual(event["k"], 2)
        self.assertEqual(event["ancestral_codon"], "ATG")
        self.assertEqual(event["derived_codon"], "TTT")
        self.assertEqual(event["mutation_type"], "NS")
        self.assertEqual(event["strand"], -1)

        # Validate path averaging results
        self.assertAlmostEqual(event["frac_S"], 0.0, places=4)
        self.assertAlmostEqual(event["frac_NS"], 2.0, places=4)
        self.assertEqual(event["num_valid_paths"], 2)

    def test_negative_strand_coordinate_boundary_validation(self):
        """
        Test that coordinate mapping is correct at gene boundaries for negative strand.

        MOCK DATA DESIGN:
        - Gene: ATGCTTACG (9 bases)
        - Prodigal: start=201, end=209, strand=-1
        - Test positions at start, middle, and end

        COORDINATE VALIDATION:
        - Highest genomic position (208) → gene index 0 (first base)
        - Middle position (204) → gene index 4 (middle base)
        - Lowest genomic position (200) → gene index 8 (last base)

        FORMULA: pos_in_gene = (end - 1) - genomic_position
        - Position 208: (209 - 1) - 208 = 0 ✓
        - Position 204: (209 - 1) - 204 = 4 ✓
        - Position 200: (209 - 1) - 200 = 8 ✓
        """
        genes = [
            {
                "gene_id": "boundary_gene",
                "start": 201,
                "end": 209,
                "strand": -1,
                "sequence": "ATGCTTACG",
            }
        ]

        # Three sites at boundaries
        sites = [
            {
                "contig": "contig1",
                "position": 208,  # First position
                "gene_id": "boundary_gene",
                "forward_ref_base": "A",
                "forward_derived_base": "T",
            },
            {
                "contig": "contig1",
                "position": 204,  # Middle position
                "gene_id": "boundary_gene",
                "forward_ref_base": "T",
                "forward_derived_base": "A",
            },
            {
                "contig": "contig1",
                "position": 200,  # Last position
                "gene_id": "boundary_gene",
                "forward_ref_base": "C",
                "forward_derived_base": "G",
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["boundary_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have 3 events
        self.assertEqual(len(codon_events_df), 3)

        # Verify coordinate mapping by checking gene_position field
        # Position 208 → gene position 0
        event_208 = codon_events_df[codon_events_df["contig_position"] == "208"].iloc[0]
        self.assertEqual(event_208["gene_position"], "0")

        # Position 204 → gene position 4
        event_204 = codon_events_df[codon_events_df["contig_position"] == "204"].iloc[0]
        self.assertEqual(event_204["gene_position"], "4")

        # Position 200 → gene position 8
        event_200 = codon_events_df[codon_events_df["contig_position"] == "200"].iloc[0]
        self.assertEqual(event_200["gene_position"], "8")

    def test_negative_strand_vs_positive_strand_complement(self):
        """
        Test that complementary sequences on opposite strands give consistent results.

        MOCK DATA DESIGN:
        - Two genes with complementary sequences
        - Positive strand: ATGCGT (positions 100-105)
        - Negative strand: ACGCAT (positions 200-205, reverse complement of ATGCGT)
        - Same mutation in both: ATG → TTG (Met → Leu)

        The negative strand gene's sequence is ACGCAT, which is the reverse complement
        of ATGCGT. When we read the negative gene 5'→3', we get ACGCAT. The first
        codon ACG codes for Threonine.

        For the mutation ATG→TTG:
        - Positive strand: position 100 changes A→T
        - Negative strand: to get same codon change...
          - Sequence is ATGCGT, first codon ATG
          - Position 205 maps to index 0
          - To change A→T in sequence, forward strand has complement T→A
          - So forward_ref='T', forward_derived='A'
        """
        genes = [
            {
                "gene_id": "pos_gene",
                "start": 101,
                "end": 106,
                "strand": 1,
                "sequence": "ATGCGT",
            },
            {
                "gene_id": "neg_gene",
                "start": 201,
                "end": 206,
                "strand": -1,
                "sequence": "ATGCGT",  # Same coding sequence
            },
        ]

        sites = [
            # Positive strand: ATG → TTG
            {
                "contig": "contig1",
                "position": 100,
                "gene_id": "pos_gene",
                "forward_ref_base": "A",
                "forward_derived_base": "T",
            },
            # Negative strand: ATG → TTG (same codon change)
            {
                "contig": "contig1",
                "position": 205,  # Maps to gene index 0
                "gene_id": "neg_gene",
                "forward_ref_base": "T",  # Complement of 'A'
                "forward_derived_base": "A",  # Complements to 'T'
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["pos_gene", "neg_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have 2 events (one per gene)
        self.assertEqual(len(codon_events_df), 2)

        # Both should have identical codon changes
        pos_event = codon_events_df[codon_events_df["gene_id"] == "pos_gene"].iloc[0]
        neg_event = codon_events_df[codon_events_df["gene_id"] == "neg_gene"].iloc[0]

        # Same codon transition
        self.assertEqual(pos_event["ancestral_codon"], neg_event["ancestral_codon"])
        self.assertEqual(pos_event["derived_codon"], neg_event["derived_codon"])

        # Same mutation type
        self.assertEqual(pos_event["mutation_type"], neg_event["mutation_type"])

        # Same fractional counts
        self.assertAlmostEqual(pos_event["frac_S"], neg_event["frac_S"], places=4)
        self.assertAlmostEqual(pos_event["frac_NS"], neg_event["frac_NS"], places=4)

        # Same potential sites
        self.assertAlmostEqual(
            pos_event["potential_S_codon"], neg_event["potential_S_codon"], places=4
        )
        self.assertAlmostEqual(
            pos_event["potential_NS_codon"], neg_event["potential_NS_codon"], places=4
        )

        # Verify strands are different
        self.assertEqual(pos_event["strand"], 1)
        self.assertEqual(neg_event["strand"], -1)

    def test_negative_strand_intermediate_stop_exclusion(self):
        """
        Test intermediate stop exclusion on negative strand gene.

        MOCK DATA DESIGN:
        - Gene: CAACGT (6 bases, strand=-1)
        - Prodigal: start=201, end=206
        - Change CAA → TGA (positions 0 and 1)

        COORDINATE MAPPING:
        - Position 205 → gene index 0 (CAA pos 0)
        - Position 204 → gene index 1 (CAA pos 1)

        MUTATIONS:
        - Pos 205: Forward C→A, Complement: G→T (C→T in sequence)
        - Pos 204: Forward T→C, Complement: A→G (A→G in sequence)
        - Result: CAA → TGA

        NG86 PATH AVERAGING:
        Path 1: CAA → TAA (intermediate stop) → EXCLUDED
        Path 2: CAA → CGA → TGA (valid, 2 NS steps)
        Result: frac_S=0.0, frac_NS=2.0, num_valid_paths=1
        """
        genes = [
            {
                "gene_id": "neg_stop_gene",
                "start": 201,
                "end": 206,
                "strand": -1,
                "sequence": "CAACGT",
            }
        ]

        sites = [
            {
                "contig": "contig1",
                "position": 205,  # CAA pos 0
                "gene_id": "neg_stop_gene",
                "forward_ref_base": "G",  # Complement of 'C'
                "forward_derived_base": "A",  # Complements to 'T'
            },
            {
                "contig": "contig1",
                "position": 204,  # CAA pos 1
                "gene_id": "neg_stop_gene",
                "forward_ref_base": "T",  # Complement of 'A'
                "forward_derived_base": "C",  # Complements to 'G'
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["neg_stop_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have 1 codon event
        self.assertEqual(len(codon_events_df), 1)

        event = codon_events_df.iloc[0]
        self.assertEqual(event["k"], 2)
        self.assertEqual(event["ancestral_codon"], "CAA")
        self.assertEqual(event["derived_codon"], "TGA")
        self.assertEqual(event["strand"], -1)

        # Only 1 valid pathway
        self.assertEqual(event["num_valid_paths"], 1)
        self.assertAlmostEqual(event["frac_S"], 0.0, places=4)
        self.assertAlmostEqual(event["frac_NS"], 2.0, places=4)

    def test_negative_strand_three_position_changes(self):
        """
        Test k=3 with 6-pathway averaging on negative strand gene.

        MOCK DATA DESIGN:
        - Gene: ATGCGT (6 bases, strand=-1)
        - Prodigal: start=201, end=206
        - Change ATG → CCC (all three positions)

        COORDINATE MAPPING:
        - Position 205 → gene index 0 (pos 0)
        - Position 204 → gene index 1 (pos 1)
        - Position 203 → gene index 2 (pos 2)

        MUTATIONS:
        - Pos 205: Forward A→G, Complement: T→C (A→C in sequence)
        - Pos 204: Forward T→G, Complement: A→C (T→C in sequence)
        - Pos 203: Forward C→G, Complement: G→C (G→C in sequence)
        - Result: ATG → CCC

        NG86 PATH AVERAGING:
        6 pathways, same as positive strand test
        Average: frac_S≈0.667, frac_NS≈2.333, num_valid_paths=6
        """
        genes = [
            {
                "gene_id": "neg_k3_gene",
                "start": 201,
                "end": 206,
                "strand": -1,
                "sequence": "ATGCGT",
            }
        ]

        sites = [
            {
                "contig": "contig1",
                "position": 205,  # ATG pos 0
                "gene_id": "neg_k3_gene",
                "forward_ref_base": "T",  # Complement of sequence[0]='A'
                "forward_derived_base": "G",  # Complements to 'C'
            },
            {
                "contig": "contig1",
                "position": 204,  # ATG pos 1
                "gene_id": "neg_k3_gene",
                "forward_ref_base": "A",  # Complement of sequence[1]='T'
                "forward_derived_base": "G",  # Complements to 'C'
            },
            {
                "contig": "contig1",
                "position": 203,  # ATG pos 2
                "gene_id": "neg_k3_gene",
                "forward_ref_base": "C",  # Complement of sequence[2]='G'
                "forward_derived_base": "G",  # Complements to 'C'
            },
        ]

        prodigal_records = build_prodigal_records(genes)
        sig_df, pre_df, post_df = build_dataframes_from_sites(sites)

        profile_by_gene = pre_df.groupby("gene_id")
        ancestral_seqs, ancestral_major = reconstruct_ancestral_sequences(
            ["neg_k3_gene"], prodigal_records, profile_by_gene
        )

        for gid, seq in ancestral_seqs.items():
            pot = calculate_potential_sites_for_gene(seq)
            prodigal_records[gid]["potential_S_sites"] = pot["S"]
            prodigal_records[gid]["potential_NS_sites"] = pot["NS"]

        codon_events_df = analyze_codon_substitutions_with_ng86_paths(
            sig_df,
            ancestral_seqs,
            prodigal_records,
            post_df,
            ancestral_major,
            ng86_cache=self.ng86_cache,
        )

        # Should have exactly 1 codon event (k=3)
        self.assertEqual(len(codon_events_df), 1)

        event = codon_events_df.iloc[0]
        self.assertEqual(event["k"], 3)
        self.assertEqual(event["ancestral_codon"], "ATG")
        self.assertEqual(event["derived_codon"], "CCC")
        self.assertEqual(event["strand"], -1)

        # Validate path averaging (same as positive strand)
        self.assertAlmostEqual(event["frac_S"], 0.667, places=2)
        self.assertAlmostEqual(event["frac_NS"], 2.333, places=2)
        self.assertEqual(event["num_valid_paths"], 6)

        # Verify coordinate fields
        self.assertEqual(event["codon_position"], "0,1,2")
        self.assertEqual(event["gene_position"], "0,1,2")
        self.assertEqual(event["contig_position"], "203,204,205")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    unittest.main(verbosity=2)
