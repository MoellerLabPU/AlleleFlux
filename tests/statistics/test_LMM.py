#!/usr/bin/env python3
"""
Unit tests for alleleflux.scripts.statistics.LMM

Covers:
- Data prep helpers: plan_csv_loading, convert_suffix_to_long, detect_frequency_columns
- Perfect separation detection
- LMM core fit wrapper behavior on LinAlgError/constant response (mocked)
- Per-position aggregation in run_lmm_for_position (mocked fit)
- Task preparation filtering logic
- CLI pipeline happy path and error handling (mocked Pool, I/O)
"""

import io
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd

from alleleflux.scripts.statistics import LMM as lmm


class TestPlanCsvLoading(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.tmp_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_column_mode_adds_time_and_base_freqs(self):
        # Create a minimal TSV to satisfy header read
        # Added subjectID column to ensure it is present
        df = pd.DataFrame(
            columns=["subjectID", "contig", "time", "A_frequency", "misc"]
        )
        f = self.tmp_path / "input.tsv"
        df.to_csv(f, sep="\t", index=False)

        base_dtype = {
            "subjectID": str,
            "contig": str,
            "gene_id": str,
            "position": int,
            "group": str,
            "replicate": str,
        }
        dtype = lmm.plan_csv_loading(base_dtype, "column", str(f))
        self.assertEqual(dtype["subjectID"], str)
        self.assertEqual(dtype["contig"], str)
        self.assertEqual(dtype["gene_id"], str)
        self.assertEqual(dtype["position"], int)
        # group/replicate remain as provided; time now set to category in column mode
        self.assertEqual(dtype["group"], str)
        self.assertEqual(dtype["replicate"], str)
        self.assertEqual(dtype["time"], "category")
        for col in ["A_frequency", "T_frequency", "G_frequency", "C_frequency"]:
            self.assertEqual(dtype[col], "float32")

    def test_suffix_mode_detects_only_suffixed_freqs(self):
        # Header contains suffixed and diff columns; only suffixed should be added
        cols = [
            "subjectID",
            "contig",
            "gene_id",
            "position",
            "A_frequency_t1",
            "G_frequency_t2",
            "T_frequency_diff",  # should be excluded by regex
            "misc",
        ]
        df = pd.DataFrame(columns=cols)
        f = self.tmp_path / "input2.tsv"
        df.to_csv(f, sep="\t", index=False)

        base_dtype = {
            "subjectID": str,
            "contig": str,
            "gene_id": str,
            "position": int,
            "group": str,
            "replicate": str,
        }
        dtype = lmm.plan_csv_loading(base_dtype, "suffix", str(f))
        self.assertEqual(dtype["subjectID"], str)
        self.assertEqual(dtype["contig"], str)
        self.assertEqual(dtype["gene_id"], str)
        self.assertEqual(dtype["position"], int)
        self.assertEqual(dtype["group"], str)
        self.assertEqual(dtype["replicate"], str)
        self.assertIn("A_frequency_t1", dtype)
        self.assertIn("G_frequency_t2", dtype)
        self.assertEqual(dtype["A_frequency_t1"], "float32")
        self.assertNotIn("T_frequency_diff", dtype)


class TestConvertSuffixToLong(unittest.TestCase):
    def test_convert_success(self):
        wide = pd.DataFrame(
            {
                "subjectID": ["S1", "S2"],
                "contig": ["c1", "c1"],
                "gene_id": ["g1", "g1"],
                "position": [1, 1],
                "replicate": ["r1", "r2"],
                # A/T have two timepoints; G/C omitted intentionally
                "A_frequency_t1": [0.1, 0.2],
                "A_frequency_t2": [0.3, 0.4],
                "T_frequency_t1": [0.5, 0.6],
                "T_frequency_t2": [0.7, 0.8],
            }
        )
        long = lmm.convert_suffix_to_long(wide)
        # 2 replicates x 2 timepoints = 4 rows
        self.assertEqual(len(long), 4)
        self.assertIn("time", long.columns)
        self.assertSetEqual(set(long["time"].unique()), {"t1", "t2"})
        # Ensure values preserved for A at one replicate
        r1 = long[long["replicate"] == "r1"].sort_values("time")
        self.assertListEqual(list(r1["A_frequency"]), [0.1, 0.3])

    def test_convert_suffix_to_long_insufficient_suffixes(self):
        wide = pd.DataFrame(
            {
                "subjectID": ["S1", "S2"],
                "contig": ["c1", "c1"],
                "gene_id": ["g1", "g1"],
                "position": [1, 1],
                "replicate": ["r1", "r2"],
                "A_frequency_t1": [0.1, 0.2],
                # Only one timepoint, insufficient
            }
        )
        with self.assertRaises(ValueError) as cm:
            lmm.convert_suffix_to_long(wide)
        self.assertIn("Less than two timepoint suffixes", str(cm.exception))


class TestHelperFunctions(unittest.TestCase):
    def test_detect_frequency_columns(self):
        self.assertEqual(
            lmm.detect_frequency_columns("longitudinal"),
            [
                "A_frequency_diff_mean",
                "T_frequency_diff_mean",
                "G_frequency_diff_mean",
                "C_frequency_diff_mean",
            ],
        )
        self.assertEqual(
            lmm.detect_frequency_columns("single"),
            ["A_frequency", "T_frequency", "G_frequency", "C_frequency"],
        )
        self.assertEqual(
            lmm.detect_frequency_columns("across_time"),
            ["A_frequency", "T_frequency", "G_frequency", "C_frequency"],
        )

    def test_has_perfect_separation_true(self):
        df = pd.DataFrame(
            {
                "subjectID": ["S1"] * 6,
                "group": ["g1"] * 3 + ["g2"] * 3,
                "y": [1, 1, 1, 2, 2, 2],
            }
        )
        self.assertTrue(lmm.has_perfect_separation(df, "y", "group"))

    def test_has_perfect_separation_false(self):
        df = pd.DataFrame(
            {
                "subjectID": ["S1"] * 4,
                "group": ["g1", "g1", "g2", "g2"],
                "y": [1.1, 1.2, 2.1, 2.2],  # variance > 0
            }
        )
        self.assertFalse(lmm.has_perfect_separation(df, "y", "group"))

    def test_detect_frequency_columns_raises(self):
        with self.assertRaises(ValueError) as cm:
            lmm.detect_frequency_columns("invalid")
        self.assertIn("Unknown data_type", str(cm.exception))


class TestFitWrapper(unittest.TestCase):
    @patch("alleleflux.scripts.statistics.LMM.smf.mixedlm")
    def test_linAlgError_constant_values_sets_p_to_one(self, mock_mixedlm):
        # Configure mixedlm().fit() to raise LinAlgError (e.g., singular matrix)
        instance = MagicMock()

        def raise_fit(*args, **kwargs):
            raise np.linalg.LinAlgError("Singular matrix")

        instance.fit.side_effect = raise_fit
        mock_mixedlm.return_value = instance

        df = pd.DataFrame(
            {
                "subjectID": ["S1", "S2", "S3", "S4"],
                "group": ["g1", "g1", "g2", "g2"],
                "replicate": ["r1", "r2", "r1", "r2"],
                "y": [5.0, 5.0, 5.0, 5.0],  # constant response
            }
        )
        res = lmm._fit_lmm_and_process_results(df, "y", "group")
        self.assertEqual(res["p_value"], 1.0)
        self.assertIn("p-value set to 1", res["notes"])  # explanatory note

    def test_fit_wrapper_predictor_already_categorical(self):
        # Ensure no error when predictor dtype is already category
        df = pd.DataFrame(
            {
                "subjectID": ["S1", "S2", "S3", "S4"],
                "group": pd.Series(["g1", "g1", "g2", "g2"], dtype="category"),
                "replicate": ["r1", "r2", "r1", "r2"],
                "y": [1.0, 1.1, 2.0, 2.1],
            }
        )
        with patch("alleleflux.scripts.statistics.LMM.smf.mixedlm") as mock_mixedlm:
            instance = MagicMock()
            # Minimal params structure to satisfy downstream logic
            instance.fit.return_value = MagicMock(
                params=pd.Series({"Intercept": 0.0, "group[T.g2]": 1.0}),
                pvalues=pd.Series({"Intercept": 0.05, "group[T.g2]": 0.01}),
                tvalues=pd.Series({"Intercept": 1.2, "group[T.g2]": 2.5}),
                converged=True,
            )
            mock_mixedlm.return_value = instance
            res = lmm._fit_lmm_and_process_results(df, "y", "group")
            self.assertAlmostEqual(res["p_value"], 0.01)
            self.assertTrue(res["converged"])  # Should reflect mocked converged flag


class TestEnsureCategorical(unittest.TestCase):
    def test_no_copy_when_all_categorical(self):
        df = pd.DataFrame(
            {
                "group": pd.Series(["g1", "g2"], dtype="category"),
                "replicate": pd.Series(["r1", "r2"], dtype="category"),
                "time": pd.Series(["t1", "t2"], dtype="category"),
            }
        )
        same_df = lmm._ensure_categorical(df, ["group", "replicate", "time"])
        self.assertIs(same_df, df)

    def test_copy_when_cast_needed(self):
        df = pd.DataFrame(
            {
                "group": ["g1", "g2"],
                "replicate": ["r1", "r2"],
                "time": ["t1", "t2"],
            }
        )
        new_df = lmm._ensure_categorical(df, ["group", "replicate", "time"])
        self.assertIsNot(new_df, df)
        self.assertTrue(
            all(
                new_df[c].dtype.name == "category"
                for c in ["group", "replicate", "time"]
            )
        )


class TestPrepareTasksForLMM(unittest.TestCase):
    def test_filters_on_groups_replicates_and_nans(self):
        # Two positions, only one valid
        df = pd.DataFrame(
            {
                "subjectID": ["S1"] * 10,
                "contig": ["c1"] * 8 + ["c1"] * 2,
                "gene_id": ["g1"] * 10,
                "position": [1] * 8 + [2] * 2,
                "group": ["g1", "g1", "g2", "g2"] * 2 + ["g1", "g1"],
                "replicate": ["r1", "r2", "r1", "r2"] * 2 + ["r1", "r1"],
                # Valid freq data for first position, NaNs for second
                "A_frequency": [0.1, 0.2, 0.3, 0.4] * 2 + [np.nan, np.nan],
                "T_frequency": [0.2, 0.3, 0.4, 0.5] * 2 + [np.nan, np.nan],
                "G_frequency": [0.3, 0.2, 0.1, 0.0] * 2 + [np.nan, np.nan],
                "C_frequency": [0.0, 0.1, 0.2, 0.3] * 2 + [np.nan, np.nan],
            }
        )

        args = type("Args", (), {"data_type": "single", "min_sample_num": 2})()
        tasks = lmm.prepare_tasks_for_lmm(df, args)
        # Only the first position passes: 1 task
        self.assertEqual(len(tasks), 1)
        contig, gene_id, pos, sub_df, data_type = tasks[0]
        self.assertEqual((contig, gene_id, pos), ("c1", "g1", 1))
        self.assertEqual(data_type, "single")
        # Sub_df has no NaNs in required columns and two groups
        self.assertFalse(
            sub_df[["A_frequency", "T_frequency", "G_frequency", "C_frequency"]]
            .isna()
            .any()
            .any()
        )
        self.assertEqual(sub_df["group"].nunique(), 2)

    def test_filters_insufficient_replicates_per_group(self):
        df = pd.DataFrame(
            {
                "subjectID": ["S1"] * 6,
                "contig": ["c1"] * 6,
                "gene_id": ["g1"] * 6,
                "position": [1] * 6,
                "group": ["g1", "g1", "g1", "g2", "g2", "g2"],
                "replicate": [
                    "r1",
                    "r2",
                    "r3",
                    "r1",
                    "r1",
                    "r1",
                ],  # g1 has 3, g2 has only r1 (1 replicate)
                "A_frequency": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
                "T_frequency": [0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
                "G_frequency": [0.3, 0.2, 0.1, 0.0, 0.1, 0.2],
                "C_frequency": [0.0, 0.1, 0.2, 0.3, 0.2, 0.1],
            }
        )

        args = type("Args", (), {"data_type": "single", "min_sample_num": 2})()
        tasks = lmm.prepare_tasks_for_lmm(df, args)
        self.assertEqual(len(tasks), 0)  # No tasks because g2 has only 1 replicate

    def test_filters_less_than_two_replicates_total(self):
        df = pd.DataFrame(
            {
                "subjectID": ["S1"] * 4,
                "contig": ["c1"] * 4,
                "gene_id": ["g1"] * 4,
                "position": [1] * 4,
                "group": ["g1", "g1", "g2", "g2"],
                "replicate": ["r1", "r1", "r1", "r1"],  # All same replicate
                "A_frequency": [0.1, 0.2, 0.3, 0.4],
                "T_frequency": [0.2, 0.3, 0.4, 0.5],
                "G_frequency": [0.3, 0.2, 0.1, 0.0],
                "C_frequency": [0.0, 0.1, 0.2, 0.3],
            }
        )

        args = type("Args", (), {"data_type": "single", "min_sample_num": 1})()
        tasks = lmm.prepare_tasks_for_lmm(df, args)
        self.assertEqual(
            len(tasks), 0
        )  # No tasks because only 1 unique replicate total


class TestMainCLI(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    @patch("alleleflux.scripts.statistics.LMM.setup_logging")
    def test_main_suffix_without_across_time_raises(self, mock_log):
        with patch(
            "sys.argv",
            [
                "LMM.py",
                "--input_df",
                "dummy.tsv",
                "--output_dir",
                str(self.temp_path / "out"),
                "--mag_id",
                "MAG001",
                "--data_type",
                "single",
                "--input_time_format",
                "suffix",
                "--cpus",
                "1",
            ],
        ):
            with self.assertRaises(SystemExit):  # parser.error raises SystemExit
                lmm.main()

    @patch("alleleflux.scripts.statistics.LMM.setup_logging")
    def test_main_across_time_without_group_raises(self, mock_log):
        with patch(
            "sys.argv",
            [
                "LMM.py",
                "--input_df",
                "dummy.tsv",
                "--output_dir",
                str(self.temp_path / "out"),
                "--mag_id",
                "MAG001",
                "--data_type",
                "across_time",
                "--cpus",
                "1",
            ],
        ):
            with self.assertRaises(SystemExit):  # parser.error raises SystemExit
                lmm.main()

    def test_main_group_not_found_raises(self):
        df = pd.DataFrame(
            {
                "subjectID": ["S1", "S2"],
                "group": ["G1", "G1"],
                "time": ["t1", "t2"],
                "A_frequency": [0.1, 0.2],
            }
        )

        with patch(
            "sys.argv",
            [
                "LMM.py",
                "--input_df",
                "dummy.tsv",
                "--output_dir",
                str(self.temp_path / "out"),
                "--mag_id",
                "MAG001",
                "--data_type",
                "across_time",
                "--group_to_analyze",
                "G2",  # not in df
                "--cpus",
                "1",
            ],
        ), patch("alleleflux.scripts.statistics.LMM.setup_logging"), patch(
            "pandas.read_csv", return_value=df
        ):
            with self.assertRaises(ValueError) as cm:
                lmm.main()
            self.assertIn("Group 'G2' not found", str(cm.exception))

    def test_main_time_column_not_found_raises(self):
        df = pd.DataFrame(
            {
                "subjectID": ["S1", "S2"],
                "group": ["G1", "G1"],
                "A_frequency": [0.1, 0.2],  # no 'time' column
            }
        )

        with patch(
            "sys.argv",
            [
                "LMM.py",
                "--input_df",
                "dummy.tsv",
                "--output_dir",
                str(self.temp_path / "out"),
                "--mag_id",
                "MAG001",
                "--data_type",
                "across_time",
                "--group_to_analyze",
                "G1",
                "--cpus",
                "1",
            ],
        ), patch("alleleflux.scripts.statistics.LMM.setup_logging"), patch(
            "pandas.read_csv", return_value=df
        ), patch(
            "alleleflux.scripts.statistics.LMM.load_and_filter_data", return_value=df
        ):
            with self.assertRaises(ValueError) as cm:
                lmm.main()
            self.assertIn("'time' column not found", str(cm.exception))

    def test_main_less_than_two_timepoints_raises(self):
        df = pd.DataFrame(
            {
                "subjectID": ["S1", "S2"],
                "group": ["G1", "G1"],
                "time": ["t1", "t1"],  # only one unique timepoint
                "A_frequency": [0.1, 0.2],
            }
        )

        with patch(
            "sys.argv",
            [
                "LMM.py",
                "--input_df",
                "dummy.tsv",
                "--output_dir",
                str(self.temp_path / "out"),
                "--mag_id",
                "MAG001",
                "--data_type",
                "across_time",
                "--group_to_analyze",
                "G1",
                "--cpus",
                "1",
            ],
        ), patch("alleleflux.scripts.statistics.LMM.setup_logging"), patch(
            "pandas.read_csv", return_value=df
        ), patch(
            "alleleflux.scripts.statistics.LMM.load_and_filter_data", return_value=df
        ):
            with self.assertRaises(ValueError) as cm:
                lmm.main()
            self.assertIn("Only two unique timepoints are required", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
