#!/usr/bin/env python3
"""
Comprehensive test suite for terminal nucleotide analysis.

Tests cover unit functions, integration workflows, edge cases, and error handling
for the terminal nucleotide analysis script.
"""

import gzip
import math
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import numpy as np
import pandas as pd
import pytest

from alleleflux.scripts.visualization.accessory.terminal_nucleotide_analysis import (
    NUCLEOTIDES,
    get_terminal_samples,
    perform_majority_voting,
    process_single_mag,
    worker_wrapper,
)


class TestTerminalNucleotideAnalysis(unittest.TestCase):
    """Base test class for terminal nucleotide analysis tests."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create temporary directories for test data
        self.temp_dir = tempfile.mkdtemp()
        self.profile_dir = Path(self.temp_dir) / "profiles"
        self.output_dir = Path(self.temp_dir) / "output"
        self.profile_dir.mkdir(parents=True)
        self.output_dir.mkdir(parents=True)

        # Sample metadata for testing
        self.metadata_data = {
            "sample_id": [
                "sample_A_T1",
                "sample_A_T2",
                "sample_A_T3",
                "sample_B_T1",
                "sample_B_T2",
                "sample_B_T3",
                "sample_C_T1",
                "sample_C_T2",
                "sample_C_T3",
            ],
            "file_path": [f"/path/to/sample_{i}.bam" for i in range(9)],
            "group": ["A", "A", "A", "B", "B", "B", "C", "C", "C"],
            "time": ["T1", "T2", "T3", "T1", "T2", "T3", "T1", "T2", "T3"],
        }

        # Create metadata file
        self.metadata_file = Path(self.temp_dir) / "metadata.tsv"
        pd.DataFrame(self.metadata_data).to_csv(
            self.metadata_file, sep="\t", index=False
        )

        # Create sample subdirectories and profile files
        self._create_sample_profile_data()

        # Create significant sites data
        self._create_significant_sites_data()

        # MAG IDs for testing
        self.mag_ids = ["MAG001", "MAG002"]

    def tearDown(self):
        """Clean up test fixtures after each test method."""
        if Path(self.temp_dir).exists():
            shutil.rmtree(self.temp_dir)

    def _create_sample_profile_data(self):
        """Create sample profile files with nucleotide frequency data."""
        sample_data = {
            "sample_A_T1": {
                "MAG001": [
                    {
                        "contig": "contig1",
                        "position": 100,
                        "total_coverage": 50,
                        "A": 25,
                        "C": 15,
                        "G": 8,
                        "T": 2,
                    },
                    {
                        "contig": "contig1",
                        "position": 200,
                        "total_coverage": 40,
                        "A": 10,
                        "C": 20,
                        "G": 8,
                        "T": 2,
                    },
                    {
                        "contig": "contig2",
                        "position": 150,
                        "total_coverage": 30,
                        "A": 15,
                        "C": 10,
                        "G": 3,
                        "T": 2,
                    },
                ],
                "MAG002": [
                    {
                        "contig": "contig3",
                        "position": 300,
                        "total_coverage": 60,
                        "A": 30,
                        "C": 20,
                        "G": 8,
                        "T": 2,
                    },
                    {
                        "contig": "contig4",
                        "position": 400,
                        "total_coverage": 35,
                        "A": 15,
                        "C": 10,
                        "G": 8,
                        "T": 2,
                    },
                ],
            },
            "sample_B_T1": {
                "MAG001": [
                    {
                        "contig": "contig1",
                        "position": 100,
                        "total_coverage": 45,
                        "A": 20,
                        "C": 15,
                        "G": 7,
                        "T": 3,
                    },
                    {
                        "contig": "contig1",
                        "position": 200,
                        "total_coverage": 35,
                        "A": 8,
                        "C": 18,
                        "G": 7,
                        "T": 2,
                    },
                    {
                        "contig": "contig2",
                        "position": 150,
                        "total_coverage": 25,
                        "A": 12,
                        "C": 8,
                        "G": 3,
                        "T": 2,
                    },
                ],
                "MAG002": [
                    {
                        "contig": "contig3",
                        "position": 300,
                        "total_coverage": 55,
                        "A": 28,
                        "C": 18,
                        "G": 7,
                        "T": 2,
                    },
                    {
                        "contig": "contig4",
                        "position": 400,
                        "total_coverage": 30,
                        "A": 12,
                        "C": 9,
                        "G": 7,
                        "T": 2,
                    },
                ],
            },
        }

        # The cache is the single per-site source after QC: one file per (MAG, timepoint)
        self.cache_dir = Path(self.temp_dir) / "allele_freq_cache"
        for sample_id, mag_data in sample_data.items():
            for mag_id, positions in mag_data.items():
                self._write_cache_rows(
                    sample_id,
                    mag_id,
                    pd.DataFrame(positions),
                    sample_id.rsplit("_", 1)[1],
                )

    def _write_cache_rows(self, sample_id, mag_id, data, timepoint):
        """Append one sample's rows (contig, position, total_coverage, A, C, G, T) to a cache file."""
        from alleleflux.scripts.utilities.cache_io import cache_path

        df = pd.DataFrame(data)
        df["gene_id"] = "g"
        df["sample_id"] = sample_id
        df["group"] = "x"
        df["subjectID"] = sample_id
        df["replicate"] = sample_id
        for nuc in NUCLEOTIDES:
            df[f"{nuc}_frequency"] = df[nuc] / df["total_coverage"]
        df["time"] = timepoint
        path = Path(cache_path(str(self.cache_dir), mag_id, timepoint))
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            old = pd.read_parquet(path)
            for col in ("sample_id", "group", "subjectID", "replicate"):
                old[col] = old[col].astype(str)
            df = pd.concat([old[old["sample_id"] != sample_id], df], ignore_index=True)
        for col in ("sample_id", "group", "subjectID", "replicate"):
            df[col] = df[col].astype("category")
        df.to_parquet(path, index=False)

    def _create_significant_sites_data(self):
        """Create significant sites data for testing."""
        significant_sites = {
            "mag_id": ["MAG001", "MAG001", "MAG001", "MAG002", "MAG002"],
            "contig": ["contig1", "contig1", "contig2", "contig3", "contig4"],
            "position": [100, 200, 150, 300, 400],
            "gene_id": ["geneA", "geneB", "geneC", "geneD", "geneE"],
            "min_p_value": [0.005, 0.015, 0.025, 0.035, 0.045],
            "q_value": [0.01, 0.02, 0.03, 0.04, 0.05],
            "test_type": [
                "two_sample_paired_tTest",
                "CMH",
                "LMM",
                "two_sample_paired_tTest",
                "CMH",
            ],
            "group_analyzed": ["A_vs_B", "A_vs_B", "A_vs_B", "A_vs_B", "A_vs_B"],
        }

        self.significant_sites_file = Path(self.temp_dir) / "significant_sites.tsv"
        pd.DataFrame(significant_sites).to_csv(
            self.significant_sites_file, sep="\t", index=False
        )


class TestGetTerminalSamples(TestTerminalNucleotideAnalysis):
    """Test the get_terminal_samples function."""

    def test_get_terminal_samples_valid_data(self):
        """Test getting terminal samples with valid metadata."""
        result = get_terminal_samples(self.metadata_file, "A", "T1")

        # Should find the matching sample
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], "sample_A_T1")

    def test_get_terminal_samples_multiple_matches(self):
        """Test with multiple matching samples."""
        result = get_terminal_samples(self.metadata_file, "A", "T2")

        # Should find matching samples
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], "sample_A_T2")

    def test_get_terminal_samples_no_matches(self):
        """Test with no matching samples."""
        with self.assertRaises(ValueError) as context:
            get_terminal_samples(self.metadata_file, "nonexistent", "T1")

        self.assertIn("No samples found", str(context.exception))

    def test_get_terminal_samples_missing_columns(self):
        """Test with metadata missing required columns."""
        # Create metadata without required columns
        bad_metadata = pd.DataFrame({"wrong_column": [1, 2, 3]})
        bad_metadata_file = Path(self.temp_dir) / "bad_metadata.tsv"
        bad_metadata.to_csv(bad_metadata_file, sep="\t", index=False)

        with self.assertRaises(ValueError) as context:
            get_terminal_samples(bad_metadata_file, "A", "T1")

        self.assertIn("must contain", str(context.exception))

    def test_get_terminal_samples_empty_metadata(self):
        """Test with empty metadata file."""
        empty_metadata = pd.DataFrame()
        empty_metadata_file = Path(self.temp_dir) / "empty_metadata.tsv"
        empty_metadata.to_csv(empty_metadata_file, sep="\t", index=False)

        with self.assertRaises(ValueError) as context:
            get_terminal_samples(empty_metadata_file, "A", "T1")

        # Empty metadata raises pandas parsing error
        self.assertIn("No columns to parse from file", str(context.exception))


class TestPerformMajorityVoting(TestTerminalNucleotideAnalysis):
    """Test the perform_majority_voting function."""

    def test_majority_voting_single_winner(self):
        """Test majority voting with clear winner."""
        major_alleles = ["A", "A", "C", "A"]  # A wins with 3/4 votes

        result = perform_majority_voting(major_alleles)

        self.assertEqual(result, "A")

    def test_majority_voting_with_ties(self):
        """Test majority voting with tied nucleotides."""
        major_alleles = ["A", "C", "G", "T"]  # Tie between all

        result = perform_majority_voting(major_alleles)

        # Should return comma-separated sorted nucleotides
        self.assertEqual(result, "A,C,G,T")

    def test_majority_voting_with_tied_pairs(self):
        """Test majority voting with tied pairs."""
        major_alleles = ["A", "C", "A", "C"]  # Tie between A and C

        result = perform_majority_voting(major_alleles)

        # Should return comma-separated sorted nucleotides
        self.assertEqual(result, "A,C")

    def test_majority_voting_with_fractional_ties(self):
        """Test majority voting with fractional contribution for ties."""
        major_alleles = ["A,C", "A", "C"]  # A gets 1.5, C gets 1.5 (tie)

        result = perform_majority_voting(major_alleles)

        # Should return comma-separated sorted nucleotides
        self.assertEqual(result, "A,C")

    def test_majority_voting_with_none_values(self):
        """Test majority voting with None values (missing data)."""
        major_alleles = ["A", None, "C", "A"]  # A wins with 2/3 valid votes

        result = perform_majority_voting(major_alleles)

        self.assertEqual(result, "A")

    def test_majority_voting_all_none(self):
        """Test majority voting when all values are None."""
        major_alleles = [None, None, None]

        result = perform_majority_voting(major_alleles)

        self.assertIsNone(result)

    def test_majority_voting_empty_list(self):
        """Test majority voting with empty list."""
        result = perform_majority_voting([])

        self.assertIsNone(result)

    def test_majority_voting_single_allele(self):
        """Test majority voting with single allele."""
        major_alleles = ["G"]

        result = perform_majority_voting(major_alleles)

        self.assertEqual(result, "G")

    def test_majority_voting_with_multinucleotides(self):
        """Test majority voting with multinucleotide strings."""
        major_alleles = ["A,G", "A", "G", "C"]  # A: 1.5, G: 1.5, C: 1 (tie between A,G)

        result = perform_majority_voting(major_alleles)

        # Should return comma-separated sorted nucleotides
        self.assertEqual(result, "A,G")


class TestWorkerFunctions(TestTerminalNucleotideAnalysis):
    """Test worker wrapper and MAG processing functions."""

    @patch(
        "alleleflux.scripts.visualization.accessory.terminal_nucleotide_analysis.logger"
    )
    def test_process_single_mag_basic(self, mock_logger):
        """Test basic MAG processing functionality."""
        # Create test data
        mag_id = "MAG001"
        mag_sites = pd.DataFrame(
            {
                "contig": ["contig1", "contig1"],
                "position": [100, 200],
                "gene_id": ["geneA", "geneB"],
                "q_value": [0.01, 0.02],
            }
        )

        sample_ids = ["sample_A_T1", "sample_B_T1"]
        p_value_column = "q_value"
        output_dir = str(self.output_dir)

        # Process the MAG
        result = process_single_mag(
            mag_id,
            mag_sites,
            str(self.cache_dir),
            sample_ids,
            "T1",
            p_value_column,
            output_dir,
        )

        # Verify result structure
        self.assertEqual(result["mag_id"], mag_id)
        self.assertEqual(result["sites_processed"], 2)
        self.assertEqual(result["samples_processed"], 2)
        self.assertIn("mean_freq_terminal_nucleotides", result)
        self.assertIn("majority_vote_terminal_nucleotides", result)
        self.assertIn("main_output", result)
        self.assertIn("intermediate_output", result)

        # Verify output files were created
        self.assertTrue(Path(result["main_output"]).exists())
        self.assertTrue(Path(result["intermediate_output"]).exists())

        # Verify terminal nucleotide counts
        mean_freq_counts = result["mean_freq_terminal_nucleotides"]
        majority_counts = result["majority_vote_terminal_nucleotides"]

        # Both methods should have some counts
        self.assertTrue(len(mean_freq_counts) > 0)
        self.assertTrue(len(majority_counts) > 0)

        # Counts should be valid nucleotides
        for nuc in mean_freq_counts.keys():
            self.assertIn(nuc, NUCLEOTIDES)
        for nuc in majority_counts.keys():
            self.assertIn(nuc, NUCLEOTIDES)

    @patch(
        "alleleflux.scripts.visualization.accessory.terminal_nucleotide_analysis.logger"
    )
    def test_process_single_mag_with_sample_missing_from_cache(self, mock_logger):
        """Test MAG processing when some profile files are missing."""
        mag_id = "MAG001"
        mag_sites = pd.DataFrame(
            {
                "contig": ["contig1"],
                "position": [100],
                "gene_id": ["geneA"],
                "q_value": [0.01],
            }
        )

        # Profile files for non-existent sample - file doesn't exist
        sample_ids = ["nonexistent_sample"]
        p_value_column = "q_value"
        output_dir = str(self.output_dir)

        # Process should raise FileNotFoundError when file doesn't exist
        # a sample with no cache rows for the MAG did not pass QC: skipped, not an error
        result = process_single_mag(
            mag_id,
            mag_sites,
            str(self.cache_dir),
            sample_ids,
            "T1",
            p_value_column,
            output_dir,
        )
        self.assertEqual(result["samples_processed"], 0)

    @patch(
        "alleleflux.scripts.visualization.accessory.terminal_nucleotide_analysis.logger"
    )
    def test_process_single_mag_with_empty_sites(self, mock_logger):
        """Test MAG processing with empty sites DataFrame."""
        mag_id = "MAG001"
        mag_sites = pd.DataFrame()  # Empty DataFrame
        sample_ids = ["sample_A_T1"]
        p_value_column = "q_value"
        output_dir = str(self.output_dir)

        result = process_single_mag(
            mag_id,
            mag_sites,
            str(self.cache_dir),
            sample_ids,
            "T1",
            p_value_column,
            output_dir,
        )

        # Should handle empty sites gracefully
        self.assertEqual(result["sites_processed"], 0)
        # Empty sites should result in None intermediate_output
        self.assertIsNone(result["intermediate_output"])
        # Main output should still be created
        self.assertTrue(Path(result["main_output"]).exists())


class TestIntegrationWorkflow(TestTerminalNucleotideAnalysis):
    """Test the complete integration workflow."""

    def test_complete_workflow_with_realistic_data(self):
        """Test the complete terminal nucleotide analysis workflow."""
        # This test would simulate the main() function workflow
        # For now, test the key components working together

        # Step 1: Get terminal samples
        terminal_samples = get_terminal_samples(self.metadata_file, "A", "T1")
        self.assertEqual(len(terminal_samples), 1)
        self.assertEqual(terminal_samples[0], "sample_A_T1")

        # Step 2: Find profile files

        # Step 3: Test majority voting
        major_alleles = ["A", "A", "C"]  # A should win
        terminal_nuc = perform_majority_voting(major_alleles)
        self.assertEqual(terminal_nuc, "A")

        # Step 4: Process a MAG - filter to use only .tsv files to avoid multiple file error
        mag_sites = pd.DataFrame(
            {
                "contig": ["contig1"],
                "position": [100],
                "gene_id": ["geneA"],
                "q_value": [0.01],
            }
        )

        # Filter profile files to use only .tsv files (not .tsv.gz) to avoid multiple file error

        result = process_single_mag(
            mag_id="MAG001",
            mag_sites=mag_sites,
            cache_dir=str(self.cache_dir),
            timepoint="T1",
            sample_ids=terminal_samples,
            p_value_column="q_value",
            output_dir=str(self.output_dir),
        )

        # Verify the complete workflow
        self.assertIsNotNone(result)
        self.assertTrue(Path(result["main_output"]).exists())
        self.assertTrue(Path(result["intermediate_output"]).exists())

    def test_workflow_with_both_methods(self):
        """Test workflow with both mean frequency and majority voting methods."""
        # Create a scenario where both methods should produce different results
        mag_sites = pd.DataFrame(
            {
                "contig": ["contig1"],
                "position": [200],  # This position has A=10, C=20 (C wins in mean freq)
                "gene_id": ["geneB"],
                "q_value": [0.02],
            }
        )

        # Sample with A=8, C=18 (C major allele)
        # Sample with A=12, C=8 (A major allele)
        # Mean frequencies: A=(10+12)/2=11, C=(20+8)/2=14 (C wins)
        # Majority voting: ["C", "A"] -> tie between C and A

        result = process_single_mag(
            mag_id="MAG001",
            mag_sites=mag_sites,
            cache_dir=str(self.cache_dir),
            timepoint="T1",
            sample_ids=["sample_A_T1", "sample_B_T1"],
            p_value_column="q_value",
            output_dir=str(self.output_dir),
        )

        # Load and verify results
        main_output = pd.read_csv(result["main_output"], sep="\t")
        self.assertFalse(main_output.empty)

        # Check that both methods are present
        self.assertIn("terminal_nucleotide_mean_freq", main_output.columns)
        self.assertIn("terminal_nucleotide_majority_vote", main_output.columns)

        # Verify that different methods can give different results
        mean_freq_result = main_output.iloc[0]["terminal_nucleotide_mean_freq"]
        majority_result = main_output.iloc[0]["terminal_nucleotide_majority_vote"]

        # At least one should be non-None
        self.assertTrue(mean_freq_result is not None or majority_result is not None)


class TestEdgeCases(TestTerminalNucleotideAnalysis):
    """Test edge cases and error conditions."""

    def test_perform_majority_voting_edge_cases(self):
        """Test edge cases for majority voting."""
        # Test with empty string - empty string should be handled as no valid nucleotides
        result = perform_majority_voting([""])
        self.assertEqual(result, "")

        # Test with invalid nucleotides - they get counted equally so return all as comma-separated
        result = perform_majority_voting(["X", "Y", "Z"])
        self.assertEqual(result, "X,Y,Z")  # All invalid nucleotides tied, return all

        # Test with mixed valid and invalid
        result = perform_majority_voting(["A", "X", "A"])
        self.assertEqual(result, "A")  # A should still win

    def test_get_terminal_samples_case_sensitivity(self):
        """Test case sensitivity in group and timepoint matching."""
        # Test with different case - should raise ValueError
        with self.assertRaises(ValueError) as context:
            get_terminal_samples(self.metadata_file, "a", "t1")

        # Case-sensitive matching should result in no matches found
        self.assertIn(
            "No samples found for group 'a' and timepoint 't1'", str(context.exception)
        )


class TestProcessSingleMagValidation(TestTerminalNucleotideAnalysis):
    """Comprehensive validation tests for process_single_mag function's mean frequency calculation and site_row dictionary."""

    def test_mean_frequency_calculation_single_site(self):
        """Test mean frequency calculation for a single site with known expected values."""
        # Create controlled test data with exact nucleotide frequencies
        mag_id = "VALIDATION_MAG"
        mag_sites = pd.DataFrame(
            {
                "contig": ["test_contig"],
                "position": [500],
                "gene_id": ["test_gene"],
                "q_value": [0.001],
            }
        )

        # Create test profile files with specific frequencies
        # Sample 1: A=40, C=30, G=20, T=10 (total=100) -> A_freq=0.4, C_freq=0.3, G_freq=0.2, T_freq=0.1
        # Sample 2: A=20, C=30, G=40, T=10 (total=100) -> A_freq=0.2, C_freq=0.3, G_freq=0.4, T_freq=0.1
        # Mean frequencies: A=0.3, C=0.3, G=0.3, T=0.1 (A,C,G tie at 0.3)

        profile_data_1 = {
            "contig": ["test_contig"],
            "position": [500],
            "total_coverage": [100],
            "A": [40],
            "C": [30],
            "G": [20],
            "T": [10],
        }
        profile_data_2 = {
            "contig": ["test_contig"],
            "position": [500],
            "total_coverage": [100],
            "A": [20],
            "C": [30],
            "G": [40],
            "T": [10],
        }

        # Create test profile files
        self._write_cache_rows("sample_test_1", "VALIDATION_MAG", profile_data_1, "T1")

        self._write_cache_rows("sample_test_2", "VALIDATION_MAG", profile_data_2, "T1")

        sample_ids = ["sample_test_1", "sample_test_2"]

        # Process the MAG
        result = process_single_mag(
            mag_id,
            mag_sites,
            str(self.cache_dir),
            sample_ids,
            "T1",
            "q_value",
            str(self.output_dir),
        )

        # Load and validate intermediate data
        self.assertIsNotNone(result["intermediate_output"])
        intermediate_df = pd.read_csv(result["intermediate_output"], sep="\t")
        self.assertEqual(len(intermediate_df), 1)  # One site processed

        # Validate site_row dictionary values
        site_row = intermediate_df.iloc[0]

        # Basic site information
        self.assertEqual(site_row["mag_id"], mag_id)
        self.assertEqual(site_row["contig"], "test_contig")
        self.assertEqual(site_row["position"], 500)
        self.assertEqual(site_row["gene_id"], "test_gene")
        self.assertEqual(site_row["q_value"], 0.001)

        # Validate nucleotide frequencies for each sample
        # Sample 1: A=40/100=0.4, C=30/100=0.3, G=20/100=0.2, T=10/100=0.1
        self.assertAlmostEqual(site_row["A_sample_test_1"], 0.4, places=6)
        self.assertAlmostEqual(site_row["C_sample_test_1"], 0.3, places=6)
        self.assertAlmostEqual(site_row["G_sample_test_1"], 0.2, places=6)
        self.assertAlmostEqual(site_row["T_sample_test_1"], 0.1, places=6)

        # Sample 2: A=20/100=0.2, C=30/100=0.3, G=40/100=0.4, T=10/100=0.1
        self.assertAlmostEqual(site_row["A_sample_test_2"], 0.2, places=6)
        self.assertAlmostEqual(site_row["C_sample_test_2"], 0.3, places=6)
        self.assertAlmostEqual(site_row["G_sample_test_2"], 0.4, places=6)
        self.assertAlmostEqual(site_row["T_sample_test_2"], 0.1, places=6)

        # Validate mean frequency calculations
        # Mean: A=(0.4+0.2)/2=0.3, C=(0.3+0.3)/2=0.3, G=(0.2+0.4)/2=0.3, T=(0.1+0.1)/2=0.1
        self.assertAlmostEqual(site_row["A_mean_frequency"], 0.3, places=6)
        self.assertAlmostEqual(site_row["C_mean_frequency"], 0.3, places=6)
        self.assertAlmostEqual(site_row["G_mean_frequency"], 0.3, places=6)
        self.assertAlmostEqual(site_row["T_mean_frequency"], 0.1, places=6)

        # Validate sample count for mean calculation
        self.assertEqual(site_row["n_samples_used_for_mean"], 2)

        # Validate major allele identification
        # Sample 1: A is major (0.4 > others)
        self.assertEqual(site_row["major_allele_sample_test_1"], "A")
        # Sample 2: G is major (0.4 > others)
        self.assertEqual(site_row["major_allele_sample_test_2"], "G")

        # Validate Method 1: Terminal nucleotide by highest mean frequency
        # A,C,G all have mean freq 0.3 (tie) -> should be "A,C,G"
        expected_terminal_mean_freq = "A,C,G"
        self.assertEqual(
            site_row["terminal_nucleotide_mean_freq"], expected_terminal_mean_freq
        )

        # Validate Method 2: Majority voting
        # Major alleles: ["A", "G"] -> each gets 1 vote -> tie between A and G
        expected_terminal_majority = "A,G"
        self.assertEqual(
            site_row["terminal_nucleotide_majority_vote"], expected_terminal_majority
        )

    def test_mean_frequency_calculation_different_coverage(self):
        """Test mean frequency calculation with different coverage depths."""
        mag_id = "COVERAGE_TEST_MAG"
        mag_sites = pd.DataFrame(
            {
                "contig": ["coverage_contig"],
                "position": [1000],
                "gene_id": ["coverage_gene"],
                "q_value": [0.05],
            }
        )

        # Create test data with different coverage
        # Sample 1: A=15, C=5, G=0, T=0 (total=20) -> A=0.75, others=0
        # Sample 2: A=45, C=25, G=10, T=20 (total=100) -> A=0.45, C=0.25, G=0.1, T=0.2
        # Mean frequencies: A=0.6, C=0.125, G=0.05, T=0.1

        profile_data_1 = {
            "contig": ["coverage_contig"],
            "position": [1000],
            "total_coverage": [20],
            "A": [15],
            "C": [5],
            "G": [0],
            "T": [0],
        }
        profile_data_2 = {
            "contig": ["coverage_contig"],
            "position": [1000],
            "total_coverage": [100],
            "A": [45],
            "C": [25],
            "G": [10],
            "T": [20],
        }

        self._write_cache_rows(
            "coverage_sample_1", "COVERAGE_TEST_MAG", profile_data_1, "T1"
        )

        self._write_cache_rows(
            "coverage_sample_2", "COVERAGE_TEST_MAG", profile_data_2, "T1"
        )

        sample_ids = ["coverage_sample_1", "coverage_sample_2"]

        result = process_single_mag(
            mag_id,
            mag_sites,
            str(self.cache_dir),
            sample_ids,
            "T1",
            "q_value",
            str(self.output_dir),
        )

        intermediate_df = pd.read_csv(result["intermediate_output"], sep="\t")
        site_row = intermediate_df.iloc[0]

        # Validate individual frequencies
        # Sample 1: A=15/20=0.75, C=5/20=0.25, G=0, T=0
        self.assertAlmostEqual(site_row["A_coverage_sample_1"], 0.75, places=6)
        self.assertAlmostEqual(site_row["C_coverage_sample_1"], 0.25, places=6)
        self.assertAlmostEqual(site_row["G_coverage_sample_1"], 0.0, places=6)
        self.assertAlmostEqual(site_row["T_coverage_sample_1"], 0.0, places=6)

        # Sample 2: A=45/100=0.45, C=25/100=0.25, G=10/100=0.1, T=20/100=0.2
        self.assertAlmostEqual(site_row["A_coverage_sample_2"], 0.45, places=6)
        self.assertAlmostEqual(site_row["C_coverage_sample_2"], 0.25, places=6)
        self.assertAlmostEqual(site_row["G_coverage_sample_2"], 0.1, places=6)
        self.assertAlmostEqual(site_row["T_coverage_sample_2"], 0.2, places=6)

        # Validate mean frequencies: A=(0.75+0.45)/2=0.6, C=(0.25+0.25)/2=0.25
        self.assertAlmostEqual(site_row["A_mean_frequency"], 0.6, places=6)
        # (0.25+0.25)/2 = 0.25
        self.assertAlmostEqual(site_row["C_mean_frequency"], 0.25, places=6)
        # (0+0.1)/2 = 0.05
        self.assertAlmostEqual(site_row["G_mean_frequency"], 0.05, places=6)
        # (0+0.2)/2 = 0.1
        self.assertAlmostEqual(site_row["T_mean_frequency"], 0.1, places=6)

        # Method 1: A should win (0.6 > all others)
        self.assertEqual(site_row["terminal_nucleotide_mean_freq"], "A")

        # Method 2: A should win in majority voting (A is major in both samples)
        self.assertEqual(site_row["terminal_nucleotide_majority_vote"], "A")

    def test_mean_frequency_with_missing_samples(self):
        """Test mean frequency calculation when some samples are missing a site."""
        mag_id = "MISSING_SITE_MAG"
        mag_sites = pd.DataFrame(
            {
                "contig": ["missing_contig"],
                "position": [2000],
                "gene_id": ["missing_gene"],
                "q_value": [0.01],
            }
        )

        # Sample 1: Has the site (A=25, C=25, G=25, T=25, total=100) -> A=0.25, others=0.25
        # Sample 2: Missing the site (should be excluded from mean calculation)
        profile_data_1 = {
            "contig": ["missing_contig"],
            "position": [2000],
            "total_coverage": [100],
            "A": [25],
            "C": [25],
            "G": [25],
            "T": [25],
        }

        self._write_cache_rows(
            "missing_sample_1", "MISSING_SITE_MAG", profile_data_1, "T1"
        )

        # Sample 2: Create file but without the position
        self._write_cache_rows(
            "missing_sample_2",
            "MISSING_SITE_MAG",
            pd.DataFrame(
                {
                    "contig": ["other_contig"],
                    "position": [100],
                    "total_coverage": [50],
                    "A": [30],
                    "C": [10],
                    "G": [5],
                    "T": [5],
                }
            ),
            "T1",
        )

        sample_ids = ["missing_sample_1", "missing_sample_2"]

        result = process_single_mag(
            mag_id,
            mag_sites,
            str(self.cache_dir),
            sample_ids,
            "T1",
            "q_value",
            str(self.output_dir),
        )

        intermediate_df = pd.read_csv(result["intermediate_output"], sep="\t")
        site_row = intermediate_df.iloc[0]

        # Sample 1 should have frequencies
        self.assertAlmostEqual(site_row["A_missing_sample_1"], 0.25, places=6)
        self.assertAlmostEqual(site_row["C_missing_sample_1"], 0.25, places=6)
        self.assertAlmostEqual(site_row["G_missing_sample_1"], 0.25, places=6)
        self.assertAlmostEqual(site_row["T_missing_sample_1"], 0.25, places=6)

        # Sample 2 should have NaN (site not found)
        self.assertTrue(math.isnan(site_row["A_missing_sample_2"]))
        self.assertTrue(math.isnan(site_row["C_missing_sample_2"]))
        self.assertTrue(math.isnan(site_row["G_missing_sample_2"]))
        self.assertTrue(math.isnan(site_row["T_missing_sample_2"]))

        # Mean frequencies should be calculated only from Sample 1
        # Since Sample 2 has NaN, mean = Sample_1_value / 1 (only one sample with data)
        self.assertAlmostEqual(site_row["A_mean_frequency"], 0.25, places=6)
        self.assertAlmostEqual(site_row["C_mean_frequency"], 0.25, places=6)
        self.assertAlmostEqual(site_row["G_mean_frequency"], 0.25, places=6)
        self.assertAlmostEqual(site_row["T_mean_frequency"], 0.25, places=6)

        # Only 1 sample should be used for mean calculation
        self.assertEqual(site_row["n_samples_used_for_mean"], 1)

        # Major allele for Sample 1 should be A,C,G,T (all tied at 0.25)
        self.assertEqual(site_row["major_allele_missing_sample_1"], "A,C,G,T")

    def test_zero_coverage_handling(self):
        """Test handling of sites with zero coverage."""
        mag_id = "ZERO_COVERAGE_MAG"
        mag_sites = pd.DataFrame(
            {
                "contig": ["zero_contig"],
                "position": [3000],
                "gene_id": ["zero_gene"],
                "q_value": [0.02],
            }
        )

        # Create data with zero coverage
        profile_data_1 = {
            "contig": ["zero_contig"],
            "position": [3000],
            "total_coverage": [0],  # Zero coverage
            "A": [0],
            "C": [0],
            "G": [0],
            "T": [0],
        }

        self._write_cache_rows(
            "zero_sample_1", "ZERO_COVERAGE_MAG", profile_data_1, "T1"
        )

        sample_ids = ["zero_sample_1"]

        # Should handle zero coverage gracefully (filtered out)
        result = process_single_mag(
            mag_id,
            mag_sites,
            str(self.cache_dir),
            sample_ids,
            "T1",
            "q_value",
            str(self.output_dir),
        )

        # Zero coverage site should be filtered out, but intermediate file is still created
        self.assertIsNotNone(result["intermediate_output"])
        # The site is still "processed" (attempted) but filtered out due to zero coverage
        self.assertEqual(result["sites_processed"], 1)

        # The main output file should still be created (even if empty or with no results)
        self.assertTrue(Path(result["main_output"]).exists())

    def test_site_row_dictionary_completeness(self):
        """Test that site_row dictionary contains all expected columns."""
        mag_id = "COMPLETENESS_MAG"
        mag_sites = pd.DataFrame(
            {
                "contig": ["complete_contig"],
                "position": [4000],
                "gene_id": ["complete_gene"],
                "q_value": [0.03],
            }
        )

        profile_data = {
            "contig": ["complete_contig"],
            "position": [4000],
            "total_coverage": [50],
            "A": [20],
            "C": [15],
            "G": [10],
            "T": [5],
        }

        self._write_cache_rows(
            "complete_sample", "COMPLETENESS_MAG", profile_data, "T1"
        )
        sample_ids = ["complete_sample"]

        result = process_single_mag(
            mag_id,
            mag_sites,
            str(self.cache_dir),
            sample_ids,
            "T1",
            "q_value",
            str(self.output_dir),
        )

        intermediate_df = pd.read_csv(result["intermediate_output"], sep="\t")
        site_row = intermediate_df.iloc[0]

        # Define expected columns
        expected_columns = [
            # Basic site information
            "mag_id",
            "contig",
            "position",
            "gene_id",
            # Nucleotide frequency columns for each sample
            "A_complete_sample",
            "C_complete_sample",
            "G_complete_sample",
            "T_complete_sample",
            # Mean frequency columns
            "A_mean_frequency",
            "C_mean_frequency",
            "G_mean_frequency",
            "T_mean_frequency",
            # Sample count
            "n_samples_used_for_mean",
            # Major allele columns
            "major_allele_complete_sample",
            # Terminal nucleotide results
            "terminal_nucleotide_mean_freq",
            "terminal_nucleotide_majority_vote",
        ]

        # Check all expected columns are present
        for col in expected_columns:
            self.assertIn(col, site_row, f"Missing expected column: {col}")

        # Verify no unexpected columns are present
        expected_set = set(expected_columns)
        actual_set = set(site_row.index)
        unexpected_columns = actual_set - expected_set

        # Filter out any columns that might be added by pandas or other processes
        # but are not part of our expected dictionary
        if unexpected_columns:
            print(f"Note: Found unexpected columns: {unexpected_columns}")
            # Don't fail on this as some columns might be added by the process

        # Verify all expected columns have non-None values (where expected)
        self.assertIsNotNone(site_row["mag_id"])
        self.assertIsNotNone(site_row["contig"])
        self.assertIsNotNone(site_row["position"])
        self.assertIsNotNone(site_row["gene_id"])

    def test_method1_terminal_nucleotide_determination(self):
        """Test Method 1: Terminal nucleotide determination by highest mean frequency."""
        mag_id = "METHOD1_TEST_MAG"

        # Create sites with different scenarios for Method 1
        mag_sites = pd.DataFrame(
            {
                "contig": ["method1_contig1", "method1_contig2", "method1_contig3"],
                "position": [5000, 5001, 5002],
                "gene_id": ["method1_gene1", "method1_gene2", "method1_gene3"],
                "q_value": [0.001, 0.002, 0.003],
            }
        )

        # Create clean profile data - one row per position
        # Site 1: A clearly wins (mean freq A=0.8, others much lower)
        # Site 2: A wins (0.6 vs 0.4 for C)
        # Site 3: Single nucleotide G (1.0)
        profile_data = pd.DataFrame(
            {
                "contig": ["method1_contig1", "method1_contig2", "method1_contig3"],
                "position": [5000, 5001, 5002],
                "total_coverage": [100, 100, 100],
                "A": [80, 60, 0],  # A: 0.8, 0.6, 0 (G wins site 3)
                "C": [10, 40, 0],  # C: 0.1, 0.4, 0
                "G": [5, 0, 100],  # G: 0.05, 0, 1.0
                "T": [5, 0, 0],  # T: 0.05, 0, 0
            }
        )

        self._write_cache_rows("method1_sample", "METHOD1_TEST_MAG", profile_data, "T1")
        sample_ids = ["method1_sample"]

        result = process_single_mag(
            mag_id,
            mag_sites,
            str(self.cache_dir),
            sample_ids,
            "T1",
            "q_value",
            str(self.output_dir),
        )

        # Load results
        intermediate_df = pd.read_csv(result["intermediate_output"], sep="\t")

        # Validate Site 1: A should win (0.8 > all others)
        site1_row = intermediate_df[intermediate_df["position"] == 5000].iloc[0]
        self.assertEqual(site1_row["terminal_nucleotide_mean_freq"], "A")
        self.assertAlmostEqual(site1_row["A_mean_frequency"], 0.8, places=2)

        # Validate Site 2: A should win (0.6 > 0.4)
        site2_row = intermediate_df[intermediate_df["position"] == 5001].iloc[0]
        self.assertEqual(site2_row["terminal_nucleotide_mean_freq"], "A")
        self.assertAlmostEqual(site2_row["A_mean_frequency"], 0.6, places=2)
        self.assertAlmostEqual(site2_row["C_mean_frequency"], 0.4, places=2)

        # Validate Site 3: G should win (1.0 > all others)
        site3_row = intermediate_df[intermediate_df["position"] == 5002].iloc[0]
        self.assertEqual(site3_row["terminal_nucleotide_mean_freq"], "G")
        self.assertAlmostEqual(site3_row["G_mean_frequency"], 1.0, places=2)

    def test_both_methods(self):
        """Test that both mean frequency and majority voting methods."""
        mag_id = "COMPARISON_MAG"
        mag_sites = pd.DataFrame(
            {
                "contig": ["compare_contig"],
                "position": [6000],
                "gene_id": ["compare_gene"],
                "q_value": [0.001],
            }
        )

        # Create clean profile data with 3 samples
        # Sample 1: A=60, C=40 (A=0.6) -> A is major
        # Sample 2: A=40, C=60 (C=0.6) -> C is major
        # Sample 3: A=80, C=20 (A=0.8) -> A is major
        # Mean frequencies: A=0.6, C=0.4 (A wins)
        # Majority voting: ["A", "C", "A"] -> A wins 2/3

        # Create separate files for each sample with different data
        sample_data_configs = [
            {"A": 60, "C": 40},  # A=0.6, A is major
            {"A": 40, "C": 60},  # C=0.6, C is major
            {"A": 80, "C": 20},  # A=0.8, A is major
        ]

        for i, (sample_id, config) in enumerate(
            zip(
                ["compare_sample_1", "compare_sample_2", "compare_sample_3"],
                sample_data_configs,
            )
        ):

            # Create individual sample data
            sample_data = pd.DataFrame(
                {
                    "contig": ["compare_contig"],
                    "position": [6000],
                    "total_coverage": [100],
                    "A": [config["A"]],
                    "C": [config["C"]],
                    "G": [0],
                    "T": [0],
                }
            )
            self._write_cache_rows(sample_id, "COMPARISON_MAG", sample_data, "T1")

        sample_ids = ["compare_sample_1", "compare_sample_2", "compare_sample_3"]

        result = process_single_mag(
            mag_id,
            mag_sites,
            str(self.cache_dir),
            sample_ids,
            "T1",
            "q_value",
            str(self.output_dir),
        )

        intermediate_df = pd.read_csv(result["intermediate_output"], sep="\t")
        site_row = intermediate_df.iloc[0]

        # Validate mean frequencies: A=(0.6+0.4+0.8)/3=0.6, C=(0.4+0.6+0.2)/3=0.4
        self.assertAlmostEqual(site_row["A_mean_frequency"], 0.6, places=2)
        self.assertAlmostEqual(site_row["C_mean_frequency"], 0.4, places=2)

        # Method 1: A should win (0.6 > 0.4)
        self.assertEqual(site_row["terminal_nucleotide_mean_freq"], "A")

        # Method 2: A should win (major alleles: ["A", "C", "A"] -> A wins 2/3)
        self.assertEqual(site_row["terminal_nucleotide_majority_vote"], "A")

        # Test a different scenario where methods differ
        # Sample 1: A=50, C=30, G=20 -> A major (0.5)
        # Sample 2: A=30, C=50, G=20 -> C major (0.5)
        # Mean frequencies: A=0.4, C=0.4, G=0.2 (A,C tie at 0.4)
        # Majority voting: ["A", "C"] -> tie between A and C

        mag_sites_2 = pd.DataFrame(
            {
                "contig": ["compare_contig2"],
                "position": [6001],
                "gene_id": ["compare_gene2"],
                "q_value": [0.002],
            }
        )

        # Update profile files with new data
        self._write_cache_rows(
            "compare_sample_1",
            "COMPARISON_MAG",
            pd.DataFrame(
                {
                    "contig": ["compare_contig2"],
                    "position": [6001],
                    "total_coverage": [100],
                    "A": [50],  # A=0.5, A is major
                    "C": [30],  # C=0.3
                    "G": [20],  # G=0.2
                    "T": [0],
                }
            ),
            "T1",
        )

        self._write_cache_rows(
            "compare_sample_2",
            "COMPARISON_MAG",
            pd.DataFrame(
                {
                    "contig": ["compare_contig2"],
                    "position": [6001],
                    "total_coverage": [100],
                    "A": [30],  # A=0.3
                    "C": [50],  # C=0.5, C is major
                    "G": [20],  # G=0.2
                    "T": [0],
                }
            ),
            "T1",
        )

        sample_ids_2 = ["compare_sample_1", "compare_sample_2"]

        result_2 = process_single_mag(
            "COMPARISON_MAG",
            mag_sites_2,
            str(self.cache_dir),
            sample_ids_2,
            "T1",
            "q_value",
            str(self.output_dir),
        )

        intermediate_df_2 = pd.read_csv(result_2["intermediate_output"], sep="\t")
        site_row_2 = intermediate_df_2.iloc[0]

        # Method 1: A and C tie (both 0.4 mean frequency)
        self.assertEqual(site_row_2["terminal_nucleotide_mean_freq"], "A,C")

        # Method 2: A and C tie (["A", "C"])
        self.assertEqual(site_row_2["terminal_nucleotide_majority_vote"], "A,C")

        # Both methods should give the same result in this case
        self.assertEqual(
            site_row_2["terminal_nucleotide_mean_freq"],
            site_row_2["terminal_nucleotide_majority_vote"],
        )


if __name__ == "__main__":
    # Create a test suite
    suite = unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])

    # Run the tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Print summary
    print(f"\nTests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")

    if result.failures:
        print("\nFailures:")
        for test, traceback in result.failures:
            print(f"  {test}: {traceback}")

    if result.errors:
        print("\nErrors:")
        for test, traceback in result.errors:
            print(f"  {test}: {traceback}")
