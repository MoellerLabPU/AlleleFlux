#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import numpy as np
import pandas as pd

from alleleflux.scripts.preprocessing.quality_control import (
    QC_SAMPLE_COLUMNS,
    init_worker,
    process_mag_files,
    qc_table_columns,
)
from alleleflux.scripts.utilities.qc_metrics import aggregate_mag_stats


class TestQualityControl(unittest.TestCase):

    def setUp(self):
        """Set up test fixtures before each test method."""
        # Sample metadata for testing
        self.metadata_dict = {
            "sample1": {
                "group": "treatment",
                "subjectID": "subject1",
                "replicate": "rep1",
                "time": "T1",
            },
            "sample2": {
                "group": "control",
                "subjectID": "subject1",
                "replicate": "rep1",
                "time": "T2",
            },
            "sample3": {
                "group": "treatment",
                "subjectID": "subject2",
                "replicate": "rep2",
                "time": "T1",
            },
            "sample4": {
                "group": "control",
                "subjectID": "subject2",
                "replicate": "rep2",
                "time": "T2",
            },
        }

        # MAG size dictionary
        self.mag_size_dict = {
            "MAG001": 1000000,  # 1M bp
            "MAG002": 2000000,  # 2M bp
        }

        # Contig length dictionary
        self.contig_length_dict = {
            "contig1": 500000,
            "contig2": 300000,
            "contig3": 200000,
        }

        # Contig to MAG mapping
        self.contig_to_mag = {
            "contig1": "MAG001",
            "contig2": "MAG001",
            "contig3": "MAG002",
        }

        # Pre-inverted MAG -> contigs index (mirrors what main() builds)
        from collections import defaultdict

        mag_to_contigs = defaultdict(set)
        for contig, mid in self.contig_to_mag.items():
            mag_to_contigs[mid].add(contig)

        # Initialize worker globals
        init_worker(
            self.metadata_dict,
            self.mag_size_dict,
            self.contig_length_dict,
            self.contig_to_mag,
            mag_to_contigs,
        )

    def create_temp_profile(self, data):
        """Create a temporary profile file with given data."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tsv")
        df = pd.DataFrame(data)
        df.to_csv(temp_file.name, sep="\t", index=False)
        temp_file.close()
        return temp_file.name


class TestProcessMagFiles(TestQualityControl):

    def test_normal_coverage_calculation(self):
        """Test normal coverage calculation with valid data."""
        # Create test profile data
        profile_data = {
            "contig": ["contig1", "contig1", "contig2", "contig2"],
            "position": [1, 2, 1, 2],
            "gene_id": ["gene1", "gene2", "gene3", "gene4"],
            "total_coverage": [10.0, 20.0, 5.0, 15.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            args = ("sample1", profile_file, "MAG001", 0.1, 1.0)
            result = process_mag_files(args)

            # Check basic structure
            self.assertEqual(result["sample_id"], "sample1")
            self.assertEqual(result["MAG_ID"], "MAG001")
            self.assertEqual(result["genome_size"], 1000000)

            # Check coverage calculations
            self.assertEqual(
                result["breadth"], 4 / 1000000
            )  # 4 positions with coverage >= 1
            self.assertEqual(
                result["average_coverage_including_zeros"], 50.0 / 1000000
            )  # sum of coverage / genome size: the QC threshold's number
            self.assertEqual(
                result["average_coverage"], 12.5
            )  # mean over the 4 covered positions: the site filter's reference
            self.assertEqual(
                result["median_coverage"], 12.5
            )  # median of [10, 20, 5, 15]

            # Check median including zeros
            expected_median_with_zeros = np.median(
                np.concatenate([[10.0, 20.0, 5.0, 15.0], np.zeros(1000000 - 4)])
            )
            self.assertEqual(
                result["median_coverage_including_zeros"], expected_median_with_zeros
            )

            # Verify that either threshold passed or a failure reason is provided
            self.assertTrue(
                result["breadth_threshold_passed"] or result["breadth_fail_reason"]
            )

        finally:
            os.unlink(profile_file)

    def test_coverage_threshold_uses_the_including_zeros_number(self):
        """One deep position in a 1 Mb genome: the covered-sites mean (100x)
        would sail past a 1x threshold, the genome-size mean (0.0001x) must
        not.  The threshold tests the latter; both are reported."""
        profile_data = {
            "contig": ["contig1"],
            "position": [1],
            "gene_id": ["gene1"],
            "total_coverage": [100.0],
        }
        profile_file = self.create_temp_profile(profile_data)
        try:
            result = process_mag_files(("sample1", profile_file, "MAG001", 0.0, 1.0))
            self.assertEqual(result["average_coverage"], 100.0)
            self.assertEqual(
                result["average_coverage_including_zeros"], 100.0 / 1000000
            )
            self.assertTrue(result["breadth_threshold_passed"])
            self.assertFalse(result["coverage_threshold_passed"])
            self.assertIn("including zeros", result["coverage_fail_reason"])
            # std over the whole genome is built from the including-zeros mean
            expected_std = np.sqrt(100.0**2 / 1000000 - (100.0 / 1000000) ** 2)
            self.assertAlmostEqual(result["coverage_std_including_zeros"], expected_std)
        finally:
            os.unlink(profile_file)

    def test_zero_coverage_warning(self):
        """Test that zero coverage values trigger a warning."""
        profile_data = {
            "contig": ["contig1", "contig1"],
            "position": [1, 2],
            "gene_id": ["gene1", "gene2"],
            "total_coverage": [10.0, 0.0],  # One zero coverage
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            with patch(
                "alleleflux.scripts.preprocessing.quality_control.logger"
            ) as mock_logger:
                args = ("sample1", profile_file, "MAG001", 0.1, 1.0)
                result = process_mag_files(args)

                # Check that warning was logged
                mock_logger.warning.assert_called()
                warning_call = mock_logger.warning.call_args[0][0]
                self.assertIn("positions with zero coverage", warning_call)

                # Check that median_coverage excludes zeros
                self.assertEqual(
                    result["median_coverage"], 10.0
                )  # Only the non-zero value

        finally:
            os.unlink(profile_file)

    def test_mag_filtering(self):
        """Test MAG filtering functionality."""
        profile_data = {
            "contig": ["contig1", "contig2", "contig3"],  # contig3 belongs to MAG002
            "position": [1, 1, 1],
            "gene_id": ["gene1", "gene2", "gene3"],
            "total_coverage": [10.0, 20.0, 30.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
                args = (
                    "sample1",
                    profile_file,
                    "MAG001",
                    0.1,
                    1.0,
                )  # Processing MAG001
                result = process_mag_files(args)

                # Should warn about filtering out contig3
                mock_logger.warning.assert_called()
                warning_call = mock_logger.warning.call_args[0][0]
                self.assertIn("1 rows not belonging to MAG MAG001", warning_call)

                # Should only include contig1 and contig2 in calculations
                self.assertEqual(result["median_coverage"], 15.0)  # median of [10, 20]

        finally:
            os.unlink(profile_file)

    def test_invalid_mag_size(self):
        """Test handling of invalid MAG sizes."""
        profile_data = {
            "contig": ["contig1"],
            "position": [1],
            "gene_id": ["gene1"],
            "total_coverage": [10.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            args = ("sample1", profile_file, "INVALID_MAG", 0.1, 1.0)
            result = process_mag_files(args)

            # Should fail with appropriate error message
            self.assertFalse(result["breadth_threshold_passed"])
            self.assertIn(
                "MAG size for INVALID_MAG not found", result["breadth_fail_reason"]
            )

        finally:
            os.unlink(profile_file)

    def test_length_weighted_coverage(self):
        """Test length-weighted coverage calculation."""
        profile_data = {
            "contig": ["contig1", "contig1", "contig2", "contig2"],
            "position": [1, 2, 1, 2],
            "gene_id": ["gene1", "gene2", "gene3", "gene4"],
            "total_coverage": [
                10.0,
                20.0,
                30.0,
                40.0,
            ],  # contig1: 30 total, contig2: 70 total
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            args = ("sample1", profile_file, "MAG001", 0.1, 1.0)
            result = process_mag_files(args)

            # contig1: mean_coverage = 30/500000, contig2: mean_coverage = 70/300000
            # length_weighted = (30 + 70) / mag_size = 100 / 1000000
            expected_lw_coverage = (30 + 70) / 1000000
            self.assertAlmostEqual(
                result["length_weighted_coverage"], expected_lw_coverage
            )

        finally:
            os.unlink(profile_file)


class TestAggregation(TestQualityControl):

    def test_aggregate_mag_stats_group(self):
        """Test MAG statistics aggregation by group."""
        df_results = pd.DataFrame(
            {
                "sample_id": [
                    "sample1",
                    "sample2",
                    "sample3",
                    "sample4",
                ],  # Add missing column
                "breadth_threshold_passed": [True, True, True, True],
                "group": ["treatment", "treatment", "control", "control"],
                "breadth": [0.8, 0.9, 0.7, 0.6],
                "average_coverage": [50.0, 60.0, 40.0, 30.0],
                "median_coverage": [45.0, 55.0, 35.0, 25.0],
                "median_coverage_including_zeros": [0.1, 0.2, 0.05, 0.03],
                "coverage_std": [15.0, 20.0, 12.0, 8.0],
                "coverage_std_including_zeros": [25.0, 30.0, 20.0, 15.0],
                "length_weighted_coverage": [48.0, 58.0, 38.0, 28.0],
                "subjects_per_group": [2, 2, 2, 2],
                "replicates_per_group": [2, 2, 2, 2],
                "paired_replicates_per_group": [1, 1, 1, 1],
            }
        )

        result = aggregate_mag_stats(df_results, "MAG001", ["group"], overall=False)

        # Should have 2 rows (treatment and control)
        self.assertEqual(len(result), 2)
        self.assertIn("MAG_ID", result.columns)
        self.assertTrue(all(result["MAG_ID"] == "MAG001"))

        # Check treatment group means
        treatment_row = result[result["group"] == "treatment"].iloc[0]
        self.assertAlmostEqual(treatment_row["breadth_mean"], 0.85)  # (0.8 + 0.9) / 2
        self.assertAlmostEqual(
            treatment_row["average_coverage_mean"], 55.0
        )  # (50 + 60) / 2

    def test_aggregate_includes_all_metrics(self):
        """Test that all metrics including breadth_genome are aggregated."""
        # Create sample data with breadth_genome field (from positions-based QC)
        df_results = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2"],
                "coverage_threshold_passed": [True, True],
                "group": ["treatment", "treatment"],
                "breadth": [0.8, 0.9],
                "breadth_genome": [
                    0.0001,
                    0.0002,
                ],  # Should be included in aggregation
                "average_coverage": [50.0, 60.0],
                "average_coverage_including_zeros": [0.05, 0.06],
                "median_coverage": [45.0, 55.0],
                "median_coverage_including_zeros": [0.1, 0.2],
                "coverage_std": [15.0, 20.0],
                "coverage_std_including_zeros": [25.0, 30.0],
                "subjects_per_group": [2, 2],
                "replicates_per_group": [2, 2],
                "paired_replicates_per_group": [1, 1],
            }
        )

        # Test group aggregation
        result = aggregate_mag_stats(df_results, "MAG001", ["group"], overall=False)

        # Both breadth metrics should be aggregated
        self.assertIn("breadth_mean", result.columns)
        self.assertIn("breadth_genome_mean", result.columns)

        # Verify breadth_genome aggregation
        self.assertAlmostEqual(result.iloc[0]["breadth_genome_mean"], 0.00015)

        # Both coverage numbers are aggregated (decision 2026-10-15)
        self.assertIn("average_coverage_including_zeros_mean", result.columns)
        self.assertAlmostEqual(
            result.iloc[0]["average_coverage_including_zeros_mean"], 0.055
        )

        # Test overall aggregation
        result_overall = aggregate_mag_stats(df_results, "MAG001", [], overall=True)

        # Both breadth metrics should be in overall results too
        self.assertIn("breadth_mean", result_overall.columns)
        self.assertIn("breadth_genome_mean", result_overall.columns)

    def test_aggregate_mag_stats_overall(self):
        """Test overall MAG statistics aggregation."""
        df_results = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2"],  # Add missing column
                "breadth_threshold_passed": [True, True],
                "breadth": [0.8, 0.6],
                "average_coverage": [50.0, 30.0],
                "median_coverage": [45.0, 25.0],
                "median_coverage_including_zeros": [0.1, 0.03],
                "coverage_std": [15.0, 8.0],
                "coverage_std_including_zeros": [25.0, 15.0],
                "length_weighted_coverage": [48.0, 28.0],
                "subjects_per_group": [2, 2],
                "replicates_per_group": [2, 2],
                "paired_replicates_per_group": [1, 1],
            }
        )

        result = aggregate_mag_stats(df_results, "MAG001", [], overall=True)

        # Should have 1 row with _mean suffix columns
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["MAG_ID"], "MAG001")
        self.assertEqual(result.iloc[0]["num_samples"], 2)
        self.assertAlmostEqual(result.iloc[0]["breadth_mean"], 0.7)  # (0.8 + 0.6) / 2
        self.assertAlmostEqual(
            result.iloc[0]["average_coverage_mean"], 40.0
        )  # (50 + 30) / 2

    def test_aggregate_empty_results(self):
        """Test aggregation with empty results."""
        df_results = pd.DataFrame()

        result = aggregate_mag_stats(df_results, "MAG001", ["group"], overall=False)

        # Should return None for empty input
        self.assertIsNone(result)


class TestEdgeCases(TestQualityControl):

    def test_empty_profile_after_filtering(self):
        """Test handling when profile is empty after MAG filtering."""
        profile_data = {
            "contig": ["contig3"],  # Only contig3, which belongs to MAG002
            "position": [1],
            "gene_id": ["gene1"],
            "total_coverage": [10.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            args = ("sample1", profile_file, "MAG001", 0.1, 1.0)  # Processing MAG001
            result = process_mag_files(args)

            # Should fail because no rows remain after filtering
            self.assertFalse(result["breadth_threshold_passed"])
            self.assertIn(
                "No rows remain after filtering", result["breadth_fail_reason"]
            )

        finally:
            os.unlink(profile_file)

    def test_all_zero_coverage(self):
        """Test profile with all zero coverage values."""
        profile_data = {
            "contig": ["contig1", "contig1"],
            "position": [1, 2],
            "gene_id": ["gene1", "gene2"],
            "total_coverage": [0.0, 0.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            args = ("sample1", profile_file, "MAG001", 0.1, 1.0)
            result = process_mag_files(args)

            # median_coverage should be NaN (no positive coverage)
            self.assertTrue(np.isnan(result["median_coverage"]))
            self.assertTrue(np.isnan(result["coverage_std"]))

            # But median_coverage_including_zeros should be 0.0
            self.assertEqual(result["median_coverage_including_zeros"], 0.0)

        finally:
            os.unlink(profile_file)

    def test_single_data_type_aggregation(self):
        """Test that aggregation works correctly when data_type='single'."""
        # Create sample data without 'time' column (single timepoint)
        df_results = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2", "sample3"],
                "breadth_threshold_passed": [True, True, True],
                "group": ["treatment", "treatment", "control"],
                "breadth": [0.8, 0.9, 0.7],
                "average_coverage": [50.0, 60.0, 40.0],
                "median_coverage": [45.0, 55.0, 35.0],
                "median_coverage_including_zeros": [0.1, 0.2, 0.05],
                "coverage_std": [15.0, 20.0, 12.0],
                "coverage_std_including_zeros": [25.0, 30.0, 20.0],
                "length_weighted_coverage": [48.0, 58.0, 38.0],
                # Note: NO 'time' column for single data type
            }
        )

        # Test group summary
        group_summary = aggregate_mag_stats(
            df_results, "MAG001", ["group"], overall=False
        )
        self.assertIsNotNone(group_summary)
        self.assertEqual(len(group_summary), 2)  # treatment and control
        # Verify _mean suffix is added to metric columns
        self.assertIn("breadth_mean", group_summary.columns)
        self.assertIn("average_coverage_mean", group_summary.columns)

        # Test overall summary (should work)
        overall_summary = aggregate_mag_stats(df_results, "MAG001", [], overall=True)
        self.assertIsNotNone(overall_summary)
        self.assertEqual(len(overall_summary), 1)

        # Test group-time summary (should return None since no 'time' column)
        # This simulates the check: if "time" in df_results.columns
        if "time" in df_results.columns:
            group_time_summary = aggregate_mag_stats(
                df_results, "MAG001", ["group", "time"], overall=False
            )
        else:
            group_time_summary = None

        self.assertIsNone(group_time_summary)  # Should be None for single data type


class TestPerSampleQcTable(unittest.TestCase):
    """The QC step writes per-sample facts only, one table per genome.

    Which subjects are readable at both timepoints and how many replicates a
    group has are questions about a COMPARISON; the eligibility step answers
    them from two of these tables.  Three contracts are pinned here, on the
    real command with tiny real inputs:

    * the table has exactly the per-sample columns (plus ``time`` for
      longitudinal data): no pairing or count column;
    * a genome whose metadata file is header-only (no profiled sample at this
      timepoint) gets a header-only QC table with the same header;
    * the summary files carry no count column, and list only genomes that
      have samples.

    Layout: genome MAG_A is one 100-base contig; sample s1 covers all 100
    positions at depth 10 (passes), sample s2 covers 5 positions (breadth
    0.05, fails the 0.1 breadth threshold).  Genome MAG_B has a header-only
    metadata file.
    """

    PAIRING_COLUMNS = {
        "two_timepoints_passed",
        "subjects_per_group",
        "replicates_per_group",
        "paired_replicates_per_group",
    }

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.meta = self.tmp / "inputMetadata_pre"
        self.out = self.tmp / "QC_pre"
        self.meta.mkdir()
        (self.tmp / "ref.fa").write_text(
            ">c1\n" + "A" * 100 + "\n>d1\n" + "C" * 100 + "\n"
        )
        (self.tmp / "mapping.tsv").write_text(
            "mag_id\tcontig_id\nMAG_A\tc1\nMAG_B\td1\n"
        )
        self.header = (
            "sample_id\tsubjectID\ttime\tgroup\treplicate\tfile_path\tbam_path\n"
        )
        rows = []
        for sample, subject, group, n_positions in (
            ("s1", "m1", "A", 100),
            ("s2", "m2", "B", 5),
        ):
            profile = self.tmp / f"{sample}_MAG_A_profiled.tsv"
            pd.DataFrame(
                {
                    "contig": ["c1"] * n_positions,
                    "position": range(n_positions),
                    "gene_id": ["g1"] * n_positions,
                    "total_coverage": [10.0] * n_positions,
                }
            ).to_csv(profile, sep="\t", index=False)
            rows.append(
                f"{sample}\t{subject}\tpre\t{group}\tcage_{subject}\t{profile}\t/bams/{sample}.bam\n"
            )
        (self.meta / "MAG_A_metadata.tsv").write_text(self.header + "".join(rows))
        # MAG_B: profiled somewhere in the run, but no sample at this timepoint.
        (self.meta / "MAG_B_metadata.tsv").write_text(self.header)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self, data_type="longitudinal"):
        return subprocess.run(
            [
                "alleleflux-qc",
                "--rootDir",
                str(self.meta),
                "--fasta",
                str(self.tmp / "ref.fa"),
                "--mag_mapping_file",
                str(self.tmp / "mapping.tsv"),
                "--output_dir",
                str(self.out),
                "--breadth_threshold",
                "0.1",
                "--coverage_threshold",
                "1",
                "--cpus",
                "2",
                "--data_type",
                data_type,
            ],
            capture_output=True,
            text=True,
        )

    def _table(self, name):
        return pd.read_csv(self.out / name, sep="\t")

    def test_qc_table_has_no_pairing_columns(self):
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        table = self._table("MAG_A_QC.tsv")
        # Exactly the per-sample columns, in the fixed order, then time.
        self.assertEqual(list(table.columns), list(QC_SAMPLE_COLUMNS) + ["time"])
        self.assertFalse(self.PAIRING_COLUMNS & set(table.columns))
        by_sample = table.set_index("sample_id")
        self.assertEqual(by_sample.loc["s1", "breadth"], 1.0)
        self.assertTrue(by_sample.loc["s1", "coverage_threshold_passed"])
        self.assertEqual(by_sample.loc["s2", "breadth"], 0.05)
        self.assertFalse(by_sample.loc["s2", "coverage_threshold_passed"])
        self.assertEqual(set(table["time"]), {"pre"})

    def test_header_only_metadata_gives_a_header_only_qc_table(self):
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        empty = self._table("MAG_B_QC.tsv")
        self.assertEqual(len(empty), 0)
        # Same header as a table with rows, so every reader finds its columns.
        self.assertEqual(list(empty.columns), list(self._table("MAG_A_QC.tsv").columns))

    def test_summaries_are_written_without_the_pairing_columns(self):
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        for name in (
            "ALL_MAGs_QC_group_summary.tsv",
            "ALL_MAGs_QC_group_time_summary.tsv",
            "ALL_MAGs_QC_overall_summary.tsv",
        ):
            summary = self._table(name)
            self.assertFalse(
                [c for c in summary.columns if "per_group" in c],
                f"{name}: {list(summary.columns)}",
            )
            # The genome without samples has nothing to summarise.
            self.assertEqual(set(summary["MAG_ID"]), {"MAG_A"})
            self.assertIn("breadth_mean", summary.columns)

    def test_no_sample_in_any_metadata_file_is_an_error(self):
        # Every file header-only: nothing to read, and a folder of header-only
        # QC tables would let the run carry on with no data.
        (self.meta / "MAG_A_metadata.tsv").write_text(self.header)
        done = self._run()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("no sample", done.stderr)

    def test_single_data_type_table_has_no_time_column(self):
        done = self._run(data_type="single")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual(
            list(self._table("MAG_A_QC.tsv").columns), list(QC_SAMPLE_COLUMNS)
        )
        self.assertEqual(
            list(self._table("MAG_B_QC.tsv").columns), list(QC_SAMPLE_COLUMNS)
        )
        self.assertFalse((self.out / "ALL_MAGs_QC_group_time_summary.tsv").exists())

    def test_column_list_matches_what_a_sample_worker_returns(self):
        # The constant is the header of an EMPTY table; the worker's result
        # keys are the header of a table with rows.  They must be the same
        # list, or the two kinds of table would disagree.
        self.assertEqual(
            qc_table_columns("longitudinal"), list(QC_SAMPLE_COLUMNS) + ["time"]
        )
        self.assertEqual(qc_table_columns("single"), list(QC_SAMPLE_COLUMNS))
        init_worker(
            {"s1": {"group": "A", "subjectID": "m1", "replicate": "c1", "time": "pre"}},
            {"MAG_A": 100},
            {"c1": 100},
            {"c1": "MAG_A"},
            {"MAG_A": {"c1"}},
        )
        profile = self.tmp / "s1_MAG_A_profiled.tsv"
        result = process_mag_files(("s1", str(profile), "MAG_A", 0.1, 1.0))
        self.assertEqual(list(result), qc_table_columns("longitudinal"))


if __name__ == "__main__":
    # Create a test suite
    suite = unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])

    # Run the tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
