# AlleleFlux/tests/test_dnds_analysis.py

import logging
import random
import unittest
from unittest import mock

import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --- Import the function to be tested
# This assumes you run the tests from the root 'AlleleFlux' directory.
# To run: python -m unittest discover tests
from alleleflux.scripts.evolution.dnds_from_timepoints import (
    _calculate_codon_sites,
    analyze_mutation_effect,
    calculate_global_dnds_for_sites,
    calculate_potential_sites_for_gene,
    get_codon_from_site,
    get_major_allele,
    reconstruct_ancestral_sequences,
)


class TestReconstructAncestralSequences(unittest.TestCase):
    """
    Unit tests for the reconstruct_ancestral_sequences function.
    This version includes multiple changes per gene to increase test coverage.
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
            "gene_is_missing",  # This gene is not in self.prodigal_records
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
        """
        Tests that a warning is logged when ref_base in profile mismatches prodigal.
        """
        # --- Modify a single ref_base to create a mismatch ---
        # For gene_rev_1 at pos 208, prodigal reverse base is 'G', so fwd is 'C'.
        # We set the profile's ref_base to 'A' to force a mismatch.
        profile_df = self.profile_by_gene.obj.copy()
        profile_df.loc[profile_df["position"] == 208, "ref_base"] = "A"
        profile_by_gene_mismatch = profile_df.groupby("gene_id")

        # Use assertLogs to capture logging output
        with self.assertLogs("root", level="WARNING") as cm:
            reconstructed_seqs, _ = reconstruct_ancestral_sequences(
                self.unique_genes, self.prodigal_records, profile_by_gene_mismatch
            )
            # Check that the correct warning message was logged
            self.assertIn(
                "WARNING:alleleflux.scripts.evolution.dnds_from_timepoints"
                ":Reference base mismatch at 208 in gene_rev_1. "
                "Prodigal (fwd strand): 'C', Profile: 'A'. "
                "Using Prodigal base for reconstruction.",
                cm.output,
            )

        # Even with the warning, the reconstruction should still complete correctly
        # using the major allele from the profile data.
        self.assertEqual(reconstructed_seqs["gene_rev_1"], "ATACTTGCA")


class TestCoreDndsFunctions(unittest.TestCase):
    """Unit tests for the core dN/dS logic and helper functions."""

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
        """Tests the calculation of potential S and N sites for single codons."""
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
        self.assertAlmostEqual(sites["N"], 6.0 + 4 / 3)
        # Gene = GAC ACG CGG TTT (Asp, Thr, Ala, Val)
        # S sites = 3 (for the three 4-fold sites)+ 0.33 (for the single 2-fold position) + 0 (for the 8 nondegenerate positions)
        # N sites =  8 nondegenerate positions + 2/3 for the 2-fold position
        gene_seq = "GACACAGCGGTT"
        sites = calculate_potential_sites_for_gene(gene_seq)
        self.assertAlmostEqual(sites["S"], 3 + (1 / 3))
        self.assertAlmostEqual(sites["N"], 8 + (2 / 3))

    def test_get_codon_from_site(self):
        """Tests the mapping of a genomic position to the correct codon."""
        # Forward strand
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
        """Tests the classification of mutations as synonymous or non-synonymous."""
        # Forward strand, non-synonymous (Arg -> Leu)
        ancestral_seq = "ATGCGTACG"
        result = analyze_mutation_effect(self.gene_info_fwd, ancestral_seq, 104, "T")
        self.assertEqual(result["mutation_type"], "NS")
        self.assertEqual(result["codon_before"], "CGT")
        self.assertEqual(result["codon_after"], "CTT")
        self.assertEqual(result["aa_before"], "R")
        self.assertEqual(result["aa_after"], "L")

        # Reverse strand, non-synonymous (Leu -> Leu)
        # This tests that the derived allele is correctly complemented.
        ancestral_seq_rev = "GTACTAACA"
        result = analyze_mutation_effect(
            self.gene_info_rev, ancestral_seq_rev, 203, "G"
        )
        self.assertEqual(result["mutation_type"], "S")
        self.assertEqual(result["codon_before"], "CTA")
        self.assertEqual(result["codon_after"], "CTC")  # G complements to C
        self.assertEqual(result["aa_before"], "L")
        self.assertEqual(result["aa_after"], "L")


class TestGlobalDndsCalculation(unittest.TestCase):
    """
    Unit tests for the calculate_global_dnds_for_sites function.
    """

    def setUp(self):
        """Set up mock data for testing the global dN/dS calculation."""

        # Mock prodigal records and ancestral sequences
        self.ancestral_seqs = {"gene1": "ATGCTTACCCGG"}  # Codons: ATG, CTT, ACC, CGG
        self.prodigal_records = {
            "gene1": {
                "record": SeqRecord(Seq(self.ancestral_seqs["gene1"]), id="gene1"),
                "start": 101,
                "end": 112,
                "strand": 1,
            }
        }

        # Mock substitution results DataFrame
        # Sequence: ATGCTTACCCGG (12 bases, positions 101-112 in 1-based coordinates)
        # Codons:   ATG(101-103) CTT(104-106) ACC(107-109) CGG(110-112)
        # Using get_codon_from_site formula: pos_in_gene = position - (start - 1) = position - 100
        # So position 101 → index 1, position 102 → index 2, etc.
        # To target index 0: need position 100
        # To target index 1: use position 101
        # To target index 3: use position 103
        # To target index 6: use position 106
        substitutions_data = {
            "mag_id": ["mag1", "mag1", "mag1", "mag1"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1"],
            "contig": ["contig1", "contig1", "contig1", "contig1"],
            "position": [101, 103, 104, 107],  # These will map to indices 1, 3, 4, 7
            "mutation_type": ["NS", "NS", "S", "S"],
            "strand": [1, 1, 1, 1],
            "major_allele_ancestral": ["T", "C", "T", "C"],
            "major_allele_derived": ["A", "A", "A", "G"],
            "codon_before": ["ATG", "CTT", "CTT", "ACC"],
            "codon_after": ["AAG", "ATT", "CAT", "AGC"],  # Mutated codons
            "aa_before": ["M", "L", "L", "T"],  # Original amino acids
            "aa_after": ["K", "I", "H", "S"],  # Mutated amino acids
        }
        self.substitutions_df = pd.DataFrame(substitutions_data)

    def test_global_dnds_calculation(self):
        """
        Tests the main dN/dS calculation logic, ensuring potential sites for a
        codon are counted only once, even with multiple mutations.
        """
        summary_df, codon_specific_df = calculate_global_dnds_for_sites(
            self.substitutions_df, self.ancestral_seqs, self.prodigal_records
        )

        # Expected values based on manual calculation:
        # Observed: 2 NS, 2 S
        # Potential sites are from codons with substitutions: ATG, CTT, ACC
        # ATG (Met): S=0, N=3
        # CTT (Leu): S=1, N=2
        # ACC (Thr): S=1, N=2
        # Total Potential S = 0 + 1 + 1 = 2.0
        # Total Potential N = 3 + 2 + 2 = 7.0
        # pN = 2 / 7.0 = 0.2857
        # pS = 2 / 2.0 = 1.0
        # dN/dS = pN / pS = 0.2857
        expected_values = {
            "Total Potential Non-Synonymous Sites (n)": "7.00",
            "Total Potential Synonymous Sites (s)": "2.00",
            "Observed Non-Synonymous Substitutions (NS)": 2,
            "Observed Synonymous Substitutions (S)": 2,
            "pN (NS / n)": "0.2857",
            "pS (S / s)": "1.0000",
            "Global dN/dS (pN/pS)": "0.2857",
        }

        # Convert the summary DataFrame into a dictionary for easy comparison
        result_dict = summary_df.set_index("Metric")["Value"].to_dict()

        self.assertDictEqual(result_dict, expected_values)

        # Test that codon_specific_df is also returned and has expected structure
        self.assertIsNotNone(codon_specific_df)
        self.assertIsInstance(codon_specific_df, pd.DataFrame)

        # Check that codon_specific_df has the expected columns
        expected_columns = [
            "mag_id",
            "contig",
            "position",
            "gene_id",
            "pos_in_gene",
            "codon_start_index",
            "position_in_codon",
            "strand",
            "major_allele_ancestral",
            "major_allele_derived",
            "codon_before",
            "codon_after",
            "aa_before",
            "aa_after",
            "mutation_type",
            "potential_S_sites_codon",
            "potential_N_sites_codon",
        ]
        for col in expected_columns:
            self.assertIn(col, codon_specific_df.columns, f"Missing column: {col}")

    def test_global_dnds_with_no_substitutions(self):
        """
        Tests that the function returns None when given an empty DataFrame.
        """
        empty_df = pd.DataFrame()
        summary_result, codon_specific_result = calculate_global_dnds_for_sites(
            empty_df, self.ancestral_seqs, self.prodigal_records
        )
        self.assertIsNone(summary_result)
        self.assertIsNone(codon_specific_result)


# This allows running the tests directly from the command line
if __name__ == "__main__":
    unittest.main()
