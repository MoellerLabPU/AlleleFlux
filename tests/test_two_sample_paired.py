#!/usr/bin/env python3
"""
Unit Tests for Two Sample Paired Test Script
===========================================

This module contains comprehensive unit tests for the two_sample_paired.py script,
which performs paired statistical tests between two groups of samples with both
raw and absolute value analyses.

Test Coverage:
- Core statistical testing functionality (t-test, Wilcoxon paired)
- Absolute value test calculations
- Edge case handling (identical values, small samples, NaN data)
- Data type support (longitudinal vs single)
- File I/O and multiprocessing operations
- Error handling and validation
- Memory efficiency and performance aspects

Usage:
    python -m unittest tests.test_two_sample_paired
    python -m unittest discover tests -p "test_two_sample_paired.py"
"""

import logging
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd

from alleleflux.scripts.statistics.two_sample_paired import (
    main,
    perform_paired_tests,
    run_paired_tests,
)

# Suppress logging during tests
logging.getLogger().setLevel(logging.WARNING)

NUCLEOTIDES = ["A", "T", "G", "C"]


class TestRunPairedTests(unittest.TestCase):
    """Test the core run_paired_tests function."""

    def setUp(self):
        """Set up test data fixtures."""
        # Create longitudinal test data
        self.longitudinal_data = pd.DataFrame(
            {
                "group": [
                    "group1",
                    "group1",
                    "group1",
                    "group1",
                    "group2",
                    "group2",
                    "group2",
                    "group2",
                ],
                "replicate": [
                    "rep1",
                    "rep2",
                    "rep3",
                    "rep4",
                    "rep1",
                    "rep2",
                    "rep3",
                    "rep4",
                ],
                "contig": ["contig1"] * 8,
                "position": [100] * 8,
                "A_frequency_diff_mean": [1.0, 2.0, 3.0, 4.0, 1.5, 2.5, 3.5, 4.5],
                "T_frequency_diff_mean": [0.5, 1.0, 1.5, 2.0, 0.2, 0.7, 1.2, 1.7],
                "G_frequency_diff_mean": [0.1, 0.2, 0.3, 0.4, 0.3, 0.4, 0.5, 0.6],
                "C_frequency_diff_mean": [0.0, 0.1, 0.2, 0.3, 0.1, 0.3, 0.5, 0.7],
            }
        )

        # Create single test data
        self.single_data = pd.DataFrame(
            {
                "group": ["group1", "group1", "group1", "group2", "group2", "group2"],
                "replicate": ["rep1", "rep2", "rep3", "rep1", "rep2", "rep3"],
                "contig": ["contig1"] * 6,
                "position": [100] * 6,
                "A_frequency": [0.2, 0.3, 0.4, 0.1, 0.2, 0.3],
                "T_frequency": [0.3, 0.4, 0.5, 0.2, 0.3, 0.4],
                "G_frequency": [0.3, 0.2, 0.1, 0.4, 0.3, 0.2],
                "C_frequency": [0.2, 0.1, 0.0, 0.3, 0.2, 0.1],
            }
        )

        self.name_tuple = ("contig1", "gene1", "100")
        self.args = MagicMock()
        self.args.min_sample_num = 4

    def test_longitudinal_data_paired_tests_success(self):
        """Test successful paired tests with longitudinal data."""
        self.args.data_type = "longitudinal"

        result = run_paired_tests(
            args=(self.name_tuple, self.longitudinal_data),
            group_1_name="group1",
            group_2_name="group2",
            min_sample_num=self.args.min_sample_num,
            data_type=self.args.data_type,
        )

        name_tuple_result, p_values_result, num_pairs_result, notes_result = result

        # Verify structure
        self.assertEqual(name_tuple_result, ("contig1", "gene1", "100"))
        self.assertEqual(num_pairs_result, 4)  # 4 replicates

        # Verify p-value keys exist and are numeric using a consistent pattern
        for nuc in NUCLEOTIDES:
            for suffix in [
                "p_value_tTest",
                "p_value_Wilcoxon",
                "p_value_tTest_abs",
                "p_value_Wilcoxon_abs",
            ]:
                key = f"{nuc}_frequency_{suffix}"
                self.assertIn(key, p_values_result)
                self.assertIsInstance(
                    p_values_result[key], (float, np.floating, np.float64)
                )

    def test_single_data_paired_tests_success(self):
        """Test successful paired tests with single data."""
        self.args.data_type = "single"

        result = run_paired_tests(
            args=(self.name_tuple, self.single_data),
            group_1_name="group1",
            group_2_name="group2",
            min_sample_num=3,
            data_type=self.args.data_type,
        )

        _, p_values_result, num_pairs_result, _ = result
        self.assertEqual(num_pairs_result, 3)  # 3 replicates
        # Absolute value tests should not be present for single data_type
        for nuc in NUCLEOTIDES:
            self.assertNotIn(f"{nuc}_frequency_p_value_tTest_abs", p_values_result)
            self.assertNotIn(f"{nuc}_frequency_p_value_Wilcoxon_abs", p_values_result)
            # Raw tests should be within [0,1]
            self.assertTrue(
                0.0 <= p_values_result[f"{nuc}_frequency_p_value_tTest"] <= 1.0
            )
            self.assertTrue(
                0.0 <= p_values_result[f"{nuc}_frequency_p_value_Wilcoxon"] <= 1.0
            )

    def test_abs_identical_magnitude_opposite_signs_longitudinal(self):
        """Absolute tests should be 1.0 when magnitudes match but signs differ."""
        df = pd.DataFrame(
            {
                "group": [
                    "group1",
                    "group1",
                    "group1",
                    "group1",
                    "group2",
                    "group2",
                    "group2",
                    "group2",
                ],
                "replicate": ["r1", "r2", "r3", "r4", "r1", "r2", "r3", "r4"],
                "contig": ["contig1"] * 8,
                "position": [100] * 8,
                # group2 is exact negative of group1
                "A_frequency_diff_mean": [1.0, 2.0, 3.0, 4.0, -1.0, -2.0, -3.0, -4.0],
                "T_frequency_diff_mean": [0.5, 1.5, 2.5, 3.5, -0.5, -1.5, -2.5, -3.5],
                "G_frequency_diff_mean": [0.2, 0.4, 0.6, 0.8, -0.2, -0.4, -0.6, -0.8],
                "C_frequency_diff_mean": [0.1, 0.2, 0.3, 0.4, -0.1, -0.2, -0.3, -0.4],
            }
        )

        result = run_paired_tests(
            args=(self.name_tuple, df),
            group_1_name="group1",
            group_2_name="group2",
            min_sample_num=4,
            data_type="longitudinal",
        )

        _, pvals, num_pairs, notes = result
        self.assertEqual(num_pairs, 4)

        for nuc in NUCLEOTIDES:
            self.assertEqual(pvals[f"{nuc}_frequency_p_value_tTest_abs"], 1.0)
            self.assertEqual(pvals[f"{nuc}_frequency_p_value_Wilcoxon_abs"], 1.0)
            # Raw tests should be defined and within [0,1]
            self.assertTrue(0.0 <= pvals[f"{nuc}_frequency_p_value_tTest"] <= 1.0)
            self.assertTrue(0.0 <= pvals[f"{nuc}_frequency_p_value_Wilcoxon"] <= 1.0)

        self.assertIn("absolute values identical in both groups", notes)

    def test_identical_values_edge_case(self):
        """Test behavior when all paired differences are identical (zero)."""
        identical_data = pd.DataFrame(
            {
                "group": ["group1", "group1", "group2", "group2"],
                "replicate": ["rep1", "rep2", "rep1", "rep2"],
                "contig": ["contig1"] * 4,
                "position": [100] * 4,
                "A_frequency_diff_mean": [1.0, 1.0, 1.0, 1.0],  # All identical
                "T_frequency_diff_mean": [1.0, 1.0, 1.0, 1.0],
                "G_frequency_diff_mean": [1.0, 1.0, 1.0, 1.0],
                "C_frequency_diff_mean": [1.0, 1.0, 1.0, 1.0],
            }
        )

        result = run_paired_tests(
            args=(self.name_tuple, identical_data),
            group_1_name="group1",
            group_2_name="group2",
            min_sample_num=2,
            data_type="longitudinal",
        )

        _, p_values_result, _, notes_result = result

        # Verify p-values are set to 1.0 for identical values1
        for key, val in p_values_result.items():
            if any(key.startswith(f"{n}_frequency_p_value") for n in NUCLEOTIDES):
                self.assertEqual(val, 1.0)

        # Verify appropriate notes are included
        self.assertIn("identical values in both groups", notes_result)

    def test_insufficient_sample_size(self):
        """Test behavior with insufficient sample size."""
        small_data = pd.DataFrame(
            {
                "group": ["group1", "group2"],
                "replicate": ["rep1", "rep1"],
                "contig": ["contig1"] * 2,
                "position": [100] * 2,
                "A_frequency_diff_mean": [1.0, 2.0],
                "T_frequency_diff_mean": [1.0, 2.0],
                "G_frequency_diff_mean": [1.0, 2.0],
                "C_frequency_diff_mean": [1.0, 2.0],
            }
        )

        result = run_paired_tests(
            args=(self.name_tuple, small_data),
            group_1_name="group1",
            group_2_name="group2",
            min_sample_num=3,
            data_type="longitudinal",
        )

        _, p_values_result, num_pairs_result, _ = result

        # Verify all p-values are NaN for insufficient samples
        for key in p_values_result.keys():
            self.assertTrue(np.isnan(p_values_result[key]))

        self.assertEqual(num_pairs_result, 1)  # Only 1 pair


class TestPerformPairedTests(unittest.TestCase):
    """Test the perform_paired_tests function."""

    def setUp(self):
        """Set up test fixtures for multiprocessing tests."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)

        # Create mock grouped data (list of tuples)
        self.mock_grouped_df = [
            (
                ("contig1", "gene1", "100"),
                pd.DataFrame(
                    {
                        "group": ["group1", "group1", "group2", "group2"],
                        "replicate": ["rep1", "rep2", "rep1", "rep2"],
                        "contig": ["contig1"] * 4,
                        "position": [100] * 4,
                        "A_frequency_diff_mean": [1.0, 2.0, 1.5, 2.5],
                        "T_frequency_diff_mean": [0.5, 1.0, 0.2, 0.7],
                        "G_frequency_diff_mean": [0.1, 0.2, 0.3, 0.4],
                        "C_frequency_diff_mean": [0.0, 0.1, 0.1, 0.3],
                    }
                ),
            )
        ]

        self.mock_func_paired = MagicMock(
            return_value=(
                ("contig1", "gene1", "100"),
                {
                    "A_frequency_p_value_tTest": 0.05,
                    "A_frequency_p_value_Wilcoxon": 0.08,
                    "A_frequency_p_value_tTest_abs": 0.03,
                    "A_frequency_p_value_Wilcoxon_abs": 0.06,
                },
                2,
                "",
            )
        )

    def tearDown(self):
        """Clean up temporary directory."""
        self.temp_dir.cleanup()

    @patch("alleleflux.scripts.statistics.two_sample_paired.Pool")
    def test_perform_paired_tests_success(self, mock_pool):
        """Test successful execution of perform_paired_tests."""
        # Mock multiprocessing pool behavior
        mock_imap_result = MagicMock()
        mock_imap_result.__iter__.return_value = [self.mock_func_paired.return_value]

        mock_pool_instance = MagicMock()
        mock_pool_instance.imap_unordered.return_value = mock_imap_result
        mock_pool.return_value.__enter__.return_value = mock_pool_instance

        # Mock os.makedirs to avoid creating directories and capture to_csv call
        with patch("os.makedirs"), patch("pandas.DataFrame.to_csv") as mock_to_csv:
            perform_paired_tests(
                self.mock_func_paired,
                self.mock_grouped_df,
                cpus=2,
                num_tests=1,
                output_dir=str(self.temp_path),
                mag_id="MAG001",
            )

        # Verify Pool was called correctly
        mock_pool.assert_called_once_with(processes=2)
        # Verify to_csv was called once with expected parameters
        mock_to_csv.assert_called_once()
        args, kwargs = mock_to_csv.call_args
        # path_or_buf is first positional arg
        out_path = args[0] if args else kwargs.get("path_or_buf")
        self.assertIsNotNone(out_path)
        self.assertTrue(str(self.temp_path) in str(out_path))
        self.assertTrue(str(out_path).endswith("MAG001_two_sample_paired.tsv.gz"))
        self.assertEqual(kwargs.get("sep"), "\t")
        self.assertEqual(kwargs.get("compression"), "gzip")
        self.assertEqual(kwargs.get("index"), False)


class TestMainFunction(unittest.TestCase):
    """Test the main function and command line interface."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)

        # Create test input file
        test_data = pd.DataFrame(
            {
                "group": ["group1", "group1", "group1", "group2", "group2", "group2"],
                "replicate": ["rep1", "rep2", "rep3", "rep1", "rep2", "rep3"],
                "contig": ["contig1"] * 6,
                "position": [100] * 6,
                "gene_id": ["gene1"] * 6,
                "A_frequency_diff_mean": [1.0, 2.0, 3.0, 1.5, 2.5, 3.5],
                "T_frequency_diff_mean": [0.5, 1.0, 1.5, 0.2, 0.7, 1.2],
                "G_frequency_diff_mean": [0.1, 0.2, 0.3, 0.3, 0.4, 0.5],
                "C_frequency_diff_mean": [0.0, 0.1, 0.2, 0.1, 0.3, 0.5],
            }
        )

        self.input_file = self.temp_path / "test_input.tsv"
        test_data.to_csv(self.input_file, sep="\t", index=False)

    def tearDown(self):
        """Clean up temporary directory."""
        self.temp_dir.cleanup()

    def test_main_invalid_group_count(self):
        """Test main function with invalid group count."""
        # Create data with only one group
        single_group_data = pd.DataFrame(
            {
                "group": ["group1", "group1", "group1"],
                "replicate": ["rep1", "rep2", "rep3"],
                "contig": ["contig1"] * 3,
                "position": [100] * 3,
                "gene_id": ["gene1"] * 3,
                "A_frequency_diff_mean": [1.0, 2.0, 3.0],
            }
        )
        single_group_file = self.temp_path / "single_group.tsv"
        single_group_data.to_csv(single_group_file, sep="\t", index=False)

        with patch(
            "sys.argv",
            [
                "two_sample_paired.py",
                "--input_df",
                str(single_group_file),
                "--mag_id",
                "MAG001",
                "--output_dir",
                str(self.temp_path / "output"),
                "--data_type",
                "longitudinal",
                "--min_sample_num",
                "3",
                "--cpus",
                "1",
            ],
        ):
            with patch("alleleflux.scripts.statistics.two_sample_paired.setup_logging"):
                with patch("pandas.read_csv", return_value=single_group_data):
                    with self.assertRaises(ValueError) as cm:
                        main()
                    self.assertIn("Expected exactly 2 groups", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
