#!/usr/bin/env python3
"""
Unit Tests for Two Sample Unpaired Test Script
=============================================

Covers:
- Core unpaired statistical testing (t-test, Mann-Whitney)
- Absolute-value test behavior including negative values
- Edge cases: identical values, insufficient sample size
- Multiprocessing + I/O pipeline in perform_unpaired_tests
- CLI validation errors
"""

import logging
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd

from alleleflux.scripts.statistics.two_sample_unpaired import (
    main,
    perform_unpaired_tests,
    run_unpaired_tests,
)

# Suppress logging during tests
logging.getLogger().setLevel(logging.WARNING)


NUCLEOTIDES = ["A", "T", "G", "C"]


class TestRunUnpairedTests(unittest.TestCase):
    """Test the core run_unpaired_tests function."""

    def setUp(self):
        self.name_tuple = ("contig1", "gene1", "100")

        # Longitudinal data: use *_diff_mean columns
        self.longitudinal_data = pd.DataFrame(
            {
                "group": ["group1"] * 5 + ["group2"] * 5,
                "contig": ["contig1"] * 10,
                "gene_id": ["gene1"] * 10,
                "position": [100] * 10,
                "A_frequency_diff_mean": [1, 2, 3, 4, 5, 1.5, 2.5, 3.5, 4.5, 5.5],
                "T_frequency_diff_mean": [0.5, 1, 1.5, 2, 2.5, 0.2, 0.7, 1.2, 1.7, 2.2],
                "G_frequency_diff_mean": [
                    0.1,
                    0.2,
                    0.3,
                    0.4,
                    0.5,
                    0.3,
                    0.4,
                    0.5,
                    0.6,
                    0.7,
                ],
                "C_frequency_diff_mean": [
                    0.0,
                    0.1,
                    0.2,
                    0.3,
                    0.4,
                    0.1,
                    0.3,
                    0.5,
                    0.7,
                    0.9,
                ],
            }
        )

        # Single data: use raw nucleotide columns
        self.single_data = pd.DataFrame(
            {
                "group": ["group1"] * 5 + ["group2"] * 5,
                "contig": ["contig1"] * 10,
                "gene_id": ["gene1"] * 10,
                "position": [100] * 10,
                "A_frequency": [0.2, 0.3, 0.4, 0.5, 0.6, 0.1, 0.2, 0.3, 0.4, 0.5],
                "T_frequency": [0.3, 0.4, 0.5, 0.6, 0.7, 0.2, 0.3, 0.4, 0.5, 0.6],
                "G_frequency": [0.3, 0.2, 0.1, 0.0, -0.1, 0.4, 0.3, 0.2, 0.1, 0.0],
                "C_frequency": [0.2, 0.1, 0.0, -0.1, -0.2, 0.3, 0.2, 0.1, 0.0, -0.1],
            }
        )

    def test_longitudinal_success(self):
        result = run_unpaired_tests(
            args=(self.name_tuple, self.longitudinal_data),
            group_1_name="group1",
            group_2_name="group2",
            min_sample_num=3,
            data_type="longitudinal",
        )

        name_tuple_result, pvals, n1, n2, notes = result
        self.assertEqual(name_tuple_result, ("contig1", "gene1", "100"))
        self.assertEqual(n1, 5)
        self.assertEqual(n2, 5)
        # Keys present and numeric
        for nuc in NUCLEOTIDES:
            for suffix in [
                "p_value_tTest",
                "p_value_MannWhitney",
                "p_value_tTest_abs",
                "p_value_MannWhitney_abs",
            ]:
                key = f"{nuc}_frequency_{suffix}"
                self.assertIn(key, pvals)
                self.assertIsInstance(pvals[key], (float, np.floating, np.float64))

    def test_single_success(self):
        result = run_unpaired_tests(
            args=(self.name_tuple, self.single_data),
            group_1_name="group1",
            group_2_name="group2",
            min_sample_num=3,
            data_type="single",
        )
        _, pvals, n1, n2, _ = result
        self.assertEqual(n1, 5)
        self.assertEqual(n2, 5)
        # Absolute value tests should not be present for single data_type
        for nuc in NUCLEOTIDES:
            self.assertNotIn(f"{nuc}_frequency_p_value_tTest_abs", pvals)
            self.assertNotIn(f"{nuc}_frequency_p_value_MannWhitney_abs", pvals)

    def test_identical_values_edge_case(self):
        # Both groups constant and identical -> tTest p set to 1; abs tTest also set to 1
        df = pd.DataFrame(
            {
                "group": ["group1"] * 4 + ["group2"] * 4,
                "contig": ["c1"] * 8,
                "gene_id": ["g1"] * 8,
                "position": [100] * 8,
                "A_frequency_diff_mean": [5, 5, 5, 5, 5, 5, 5, 5],
                "T_frequency_diff_mean": [5, 5, 5, 5, 5, 5, 5, 5],
                "G_frequency_diff_mean": [5, 5, 5, 5, 5, 5, 5, 5],
                "C_frequency_diff_mean": [5, 5, 5, 5, 5, 5, 5, 5],
            }
        )
        _, pvals, n1, n2, notes = run_unpaired_tests(
            args=(self.name_tuple, df),
            group_1_name="group1",
            group_2_name="group2",
            min_sample_num=4,
            data_type="longitudinal",
        )
        self.assertEqual(n1, 4)
        self.assertEqual(n2, 4)
        # Verify p-values are set to 1.0 for identical values
        for key, val in pvals.items():
            if any(key.startswith(f"{n}_frequency_p_value") for n in NUCLEOTIDES):
                self.assertEqual(val, 1.0)

        self.assertIn("identical values in both groups", notes)

    def test_abs_with_negative_values(self):
        # Same absolute distributions but different signs -> abs tests ~1, raw tests likely significant
        vals = list(range(1, 11))
        df = pd.DataFrame(
            {
                "group": ["group1"] * 10 + ["group2"] * 10,
                "contig": ["c1"] * 20,
                "gene_id": ["g1"] * 20,
                "position": [100] * 20,
                "A_frequency": vals + ([-v for v in vals]),
                "T_frequency": [v * 0.5 for v in vals] + ([-v * 0.5 for v in vals]),
                "G_frequency": [v * 0.2 for v in vals] + ([-v * 0.2 for v in vals]),
                "C_frequency": [v * 0.1 for v in vals] + ([-v * 0.1 for v in vals]),
            }
        )
        _, pvals, n1, n2, notes = run_unpaired_tests(
            args=(self.name_tuple, df),
            group_1_name="group1",
            group_2_name="group2",
            min_sample_num=5,
            data_type="single",
        )
        self.assertEqual(n1, 10)
        self.assertEqual(n2, 10)
        for nuc in NUCLEOTIDES:
            # Abs tests not present for single
            self.assertNotIn(f"{nuc}_frequency_p_value_tTest_abs", pvals)
            self.assertNotIn(f"{nuc}_frequency_p_value_MannWhitney_abs", pvals)
            # Raw tests should be within [0,1]
            self.assertTrue(0.0 <= pvals[f"{nuc}_frequency_p_value_tTest"] <= 1.0)
            self.assertTrue(0.0 <= pvals[f"{nuc}_frequency_p_value_MannWhitney"] <= 1.0)

    def test_insufficient_sample_size(self):
        # Fewer than min samples in one group -> all NaNs
        df = pd.DataFrame(
            {
                "group": ["group1", "group2"],
                "contig": ["c1", "c1"],
                "gene_id": ["g1", "g1"],
                "position": [100, 100],
                "A_frequency": [0.1, 0.2],
                "T_frequency": [0.1, 0.2],
                "G_frequency": [0.1, 0.2],
                "C_frequency": [0.1, 0.2],
            }
        )
        _, pvals, n1, n2, _ = run_unpaired_tests(
            args=(self.name_tuple, df),
            group_1_name="group1",
            group_2_name="group2",
            min_sample_num=3,
            data_type="single",
        )
        self.assertEqual(n1, 1)
        self.assertEqual(n2, 1)
        for v in pvals.values():
            self.assertTrue(np.isnan(v))


class TestPerformUnpairedTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)

        self.mock_grouped = [
            (
                ("contig1", "gene1", "100"),
                pd.DataFrame(
                    {
                        "group": ["group1", "group1", "group2", "group2"],
                        "contig": ["contig1"] * 4,
                        "gene_id": ["gene1"] * 4,
                        "position": [100] * 4,
                        "A_frequency_diff_mean": [1.0, 2.0, 1.5, 2.5],
                        "T_frequency_diff_mean": [0.5, 1.0, 0.2, 0.7],
                        "G_frequency_diff_mean": [0.1, 0.2, 0.3, 0.4],
                        "C_frequency_diff_mean": [0.0, 0.1, 0.1, 0.3],
                    }
                ),
            )
        ]

        # Mock result of run_unpaired_tests
        self.mock_func = MagicMock(
            return_value=(
                ("contig1", "gene1", "100"),
                {
                    "A_frequency_p_value_tTest": 0.05,
                    "A_frequency_p_value_MannWhitney": 0.08,
                    "A_frequency_p_value_tTest_abs": 0.97,
                    "A_frequency_p_value_MannWhitney_abs": 0.99,
                },
                2,
                2,
                "",
            )
        )

    def tearDown(self):
        self.temp_dir.cleanup()

    @patch("alleleflux.scripts.statistics.two_sample_unpaired.Pool")
    def test_perform_unpaired_tests_success(self, mock_pool):
        mock_imap = MagicMock()
        mock_imap.__iter__.return_value = [self.mock_func.return_value]
        mock_pool_instance = MagicMock()
        mock_pool_instance.imap_unordered.return_value = mock_imap
        mock_pool.return_value.__enter__.return_value = mock_pool_instance

        with patch("os.makedirs"), patch("pandas.DataFrame.to_csv") as mock_to_csv:
            perform_unpaired_tests(
                self.mock_func,
                self.mock_grouped,
                group_1="group1",
                group_2="group2",
                cpus=2,
                num_tests=1,
                output_dir=str(self.temp_path),
                mag_id="MAG001",
            )

        mock_pool.assert_called_once_with(processes=2)
        mock_to_csv.assert_called_once()
        args, kwargs = mock_to_csv.call_args
        out_path = args[0] if args else kwargs.get("path_or_buf")
        self.assertIsNotNone(out_path)
        self.assertTrue(str(self.temp_path) in str(out_path))
        self.assertTrue(str(out_path).endswith("MAG001_two_sample_unpaired.tsv.gz"))
        self.assertEqual(kwargs.get("sep"), "\t")
        self.assertEqual(kwargs.get("compression"), "gzip")
        self.assertEqual(kwargs.get("index"), False)


class TestMainFunction(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_main_invalid_group_count(self):
        single_group = pd.DataFrame(
            {
                "group": ["group1", "group1", "group1"],
                "contig": ["c1"] * 3,
                "gene_id": ["g1"] * 3,
                "position": [100] * 3,
                "A_frequency_diff_mean": [1.0, 2.0, 3.0],
            }
        )

        with patch(
            "sys.argv",
            [
                "two_sample_unpaired.py",
                "--input_df",
                str(self.temp_path / "input.tsv"),
                "--mag_id",
                "MAG001",
                "--output_dir",
                str(self.temp_path / "out"),
                "--data_type",
                "longitudinal",
                "--min_sample_num",
                "3",
                "--cpus",
                "1",
            ],
        ):
            with patch(
                "alleleflux.scripts.statistics.two_sample_unpaired.setup_logging"
            ), patch("pandas.read_csv", return_value=single_group):
                with self.assertRaises(ValueError) as cm:
                    main()
                self.assertIn("Expected exactly 2 groups", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
