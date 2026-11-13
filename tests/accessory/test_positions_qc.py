#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import numpy as np
import pandas as pd

from alleleflux.scripts.accessory.positions_qc import (
    apply_positions_filter,
    init_worker,
    load_positions_file,
    process_mag,
    process_positions_mag_sample,
)
from alleleflux.scripts.utilities.qc_metrics import aggregate_mag_stats


class TestPositionsQC(unittest.TestCase):
    """Base test class for positions QC tests."""

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
        }

        # MAG size dictionary
        self.mag_size_dict = {
            "MAG001": 1000000,  # 1M bp
            "MAG002": 2000000,  # 2M bp
        }

        # Contig to MAG mapping
        self.contig_to_mag = {
            "contig1": "MAG001",
            "contig2": "MAG001",
            "contig3": "MAG002",
        }

        # Positions filter map (MAG -> set of (contig, position) tuples)
        self.positions_filter = {
            "MAG001": {("contig1", 1), ("contig1", 2), ("contig2", 1)},
            "MAG002": {("contig3", 1), ("contig3", 2)},
        }

    def create_temp_profile(self, data):
        """Create a temporary profile file with given data."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tsv")
        df = pd.DataFrame(data)
        df.to_csv(temp_file.name, sep="\t", index=False)
        temp_file.close()
        return temp_file.name

    def create_temp_positions_file(self, data):
        """Create a temporary positions file with given data."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tsv")
        df = pd.DataFrame(data)
        df.to_csv(temp_file.name, sep="\t", index=False)
        temp_file.close()
        return temp_file.name


class TestLoadPositionsFile(TestPositionsQC):
    """Test the load_positions_file function."""

    def test_load_valid_positions_file(self):
        """Test loading a valid positions file."""
        positions_data = {
            "MAG": ["MAG001", "MAG001", "MAG001", "MAG002"],
            "contig": ["contig1", "contig1", "contig2", "contig3"],
            "position": [1, 2, 1, 5],
        }

        positions_file = self.create_temp_positions_file(positions_data)

        try:
            result = load_positions_file(Path(positions_file))

            # Check structure
            self.assertIsInstance(result, dict)
            self.assertIn("MAG001", result)
            self.assertIn("MAG002", result)

            # Check MAG001 positions (set of tuples)
            self.assertIsInstance(result["MAG001"], set)
            self.assertEqual(len(result["MAG001"]), 3)
            self.assertIn(("contig1", 1), result["MAG001"])
            self.assertIn(("contig1", 2), result["MAG001"])
            self.assertIn(("contig2", 1), result["MAG001"])

            # Check MAG002 positions
            self.assertEqual(len(result["MAG002"]), 1)
            self.assertIn(("contig3", 5), result["MAG002"])

        finally:
            os.unlink(positions_file)

    def test_load_positions_file_with_duplicates(self):
        """Test that duplicate positions raise an error."""
        positions_data = {
            "MAG": ["MAG001", "MAG001", "MAG001"],
            "contig": ["contig1", "contig1", "contig1"],
            "position": [1, 2, 1],  # Position 1 duplicated
        }

        positions_file = self.create_temp_positions_file(positions_data)

        try:
            with self.assertRaises(ValueError) as context:
                load_positions_file(Path(positions_file))

            self.assertIn("duplicate", str(context.exception).lower())
            self.assertIn("2", str(context.exception))  # 2 duplicate rows

        finally:
            os.unlink(positions_file)

    def test_load_positions_file_path_object(self):
        """Test that load_positions_file accepts Path objects."""
        positions_data = {
            "MAG": ["MAG001"],
            "contig": ["contig1"],
            "position": [1],
        }

        positions_file = self.create_temp_positions_file(positions_data)

        try:
            # Should accept Path object
            result = load_positions_file(Path(positions_file))
            self.assertIn("MAG001", result)

        finally:
            os.unlink(positions_file)

    def test_load_positions_file_string_path(self):
        """Test that load_positions_file accepts string paths."""
        positions_data = {
            "MAG": ["MAG001"],
            "contig": ["contig1"],
            "position": [1],
        }

        positions_file = self.create_temp_positions_file(positions_data)

        try:
            # Should accept string path
            result = load_positions_file(positions_file)
            self.assertIn("MAG001", result)

        finally:
            os.unlink(positions_file)

    def test_load_positions_empty_file(self):
        """Test loading an empty positions file."""
        positions_data = {"MAG": [], "contig": [], "position": []}

        positions_file = self.create_temp_positions_file(positions_data)

        try:
            result = load_positions_file(Path(positions_file))
            self.assertEqual(len(result), 0)

        finally:
            os.unlink(positions_file)


class TestApplyPositionsFilter(TestPositionsQC):
    """Test the apply_positions_filter function."""

    def test_filter_basic(self):
        """Test basic position filtering."""
        df = pd.DataFrame(
            {
                "contig": ["contig1", "contig1", "contig1", "contig2"],
                "position": [1, 2, 3, 1],
                "total_coverage": [10.0, 20.0, 30.0, 40.0],
            }
        )

        wanted = {("contig1", 1), ("contig1", 2), ("contig2", 1)}

        filtered_df, universe_n = apply_positions_filter(
            df, wanted, "MAG001", "sample1"
        )

        # Should only include positions 1, 2, and 4 (contig1:1, contig1:2, contig2:1)
        self.assertEqual(len(filtered_df), 3)
        self.assertEqual(universe_n, 3)
        self.assertTrue(all(filtered_df["position"].isin([1, 2])))

    def test_filter_no_matches(self):
        """Test filtering when no positions match."""
        df = pd.DataFrame(
            {
                "contig": ["contig1", "contig1"],
                "position": [10, 20],
                "total_coverage": [10.0, 20.0],
            }
        )

        wanted = {("contig1", 1), ("contig1", 2)}

        with patch("alleleflux.scripts.accessory.positions_qc.logger") as mock_logger:
            filtered_df, universe_n = apply_positions_filter(
                df, wanted, "MAG001", "sample1"
            )

            # Should return empty dataframe
            self.assertTrue(filtered_df.empty)
            self.assertEqual(universe_n, 2)

            # Should log info message
            mock_logger.info.assert_called()

    def test_filter_missing_columns(self):
        """Test filtering with missing required columns."""
        df = pd.DataFrame({"some_column": [1, 2]})

        wanted = {("contig1", 1)}

        with self.assertRaises(ValueError) as context:
            apply_positions_filter(df, wanted, "MAG001", "sample1")

        self.assertIn("contig", str(context.exception).lower())
        self.assertIn("position", str(context.exception).lower())

    def test_filter_type_conversion(self):
        """Test that filtering handles type conversions correctly."""
        df = pd.DataFrame(
            {
                "contig": [1, 2],  # Integer contigs
                "position": ["1", "2"],  # String positions
                "total_coverage": [10.0, 20.0],
            }
        )

        wanted = {("1", 1), ("2", 2)}  # String contig, int position

        filtered_df, universe_n = apply_positions_filter(
            df, wanted, "MAG001", "sample1"
        )

        # Should handle type conversions and match both
        self.assertEqual(len(filtered_df), 2)


class TestInitWorker(TestPositionsQC):
    """Test the init_worker function."""

    def test_init_worker_sets_globals(self):
        """Test that init_worker sets up global variables correctly."""
        init_worker(
            self.metadata_dict,
            self.mag_size_dict,
            self.contig_to_mag,
            self.positions_filter,
            "positions",
        )

        # Import the globals from the module
        from alleleflux.scripts.accessory import positions_qc

        self.assertEqual(positions_qc.metadata_dict, self.metadata_dict)
        self.assertEqual(positions_qc.mag_size_dict, self.mag_size_dict)
        self.assertEqual(positions_qc.contig_to_mag, self.contig_to_mag)
        self.assertEqual(positions_qc.positions_filter, self.positions_filter)
        self.assertEqual(positions_qc.positions_den_mode, "positions")

    def test_init_worker_genome_denominator(self):
        """Test init_worker with genome denominator mode."""
        init_worker(
            self.metadata_dict,
            self.mag_size_dict,
            self.contig_to_mag,
            self.positions_filter,
            "genome",
        )

        from alleleflux.scripts.accessory import positions_qc

        self.assertEqual(positions_qc.positions_den_mode, "genome")


class TestProcessPositionsMagSample(TestPositionsQC):
    """Test the process_positions_mag_sample function."""

    def test_positions_mode_basic(self):
        """Test basic processing with positions denominator."""
        # Create profile with 4 positions, filter to 3
        profile_data = {
            "contig": ["contig1", "contig1", "contig1", "contig2"],
            "position": [1, 2, 3, 1],
            "gene_id": ["gene1", "gene2", "gene3", "gene4"],
            "total_coverage": [10.0, 20.0, 30.0, 40.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            # Initialize worker
            init_worker(
                self.metadata_dict,
                self.mag_size_dict,
                self.contig_to_mag,
                self.positions_filter,
                "positions",
            )

            args = ("sample1", profile_file, "MAG001", 0.000001, 1.0)
            result = process_positions_mag_sample(args)

            # Check basic structure
            self.assertEqual(result["sample_id"], "sample1")
            self.assertEqual(result["MAG_ID"], "MAG001")
            self.assertEqual(result["positions_considered"], 3)

            # Check breadth (filtered): 3 positions with coverage / 3 total = 1.0
            self.assertAlmostEqual(result["breadth"], 1.0)

            # Check breadth_genome: 4 positions with coverage / 1M genome = 0.000004
            self.assertAlmostEqual(result["breadth_genome"], 4 / 1000000)

            # Check average_coverage (covered-sites mean of the filtered rows)
            # and the including-zeros mean over the 3 wanted positions: all
            # three are covered, so both are (10 + 20 + 40) / 3 = 23.33...
            self.assertAlmostEqual(result["average_coverage"], 70.0 / 3)
            self.assertAlmostEqual(result["average_coverage_including_zeros"], 70.0 / 3)

            # Check average_coverage_genome: 100 / 1M
            self.assertAlmostEqual(result["average_coverage_genome"], 100.0 / 1000000)

            # Should pass thresholds (genome-wide checked)
            self.assertTrue(result["breadth_threshold_passed"])
            # Coverage will fail because 100/1M = 0.0001 < 1.0 threshold
            self.assertFalse(result["coverage_threshold_passed"])

        finally:
            os.unlink(profile_file)

    def test_genome_mode_basic(self):
        """Test basic processing with genome denominator."""
        profile_data = {
            "contig": ["contig1", "contig1"],
            "position": [1, 2],
            "gene_id": ["gene1", "gene2"],
            "total_coverage": [10.0, 20.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            # Initialize worker with genome mode
            init_worker(
                self.metadata_dict,
                self.mag_size_dict,
                self.contig_to_mag,
                self.positions_filter,
                "genome",
            )

            args = ("sample1", profile_file, "MAG001", 0.000001, 1.0)
            result = process_positions_mag_sample(args)

            # Check positions_considered
            self.assertEqual(result["positions_considered"], 3)

            # Check breadth: 2 positions / 1M genome
            self.assertAlmostEqual(result["breadth"], 2.0 / 1000000)

            # Check breadth_genome (should still be calculated)
            self.assertAlmostEqual(result["breadth_genome"], 2.0 / 1000000)

            # Check average_coverage_including_zeros: 30 / 1M; the covered-sites
            # mean is 30 / 2 (decision 2026-10-15: two numbers, two jobs)
            self.assertAlmostEqual(
                result["average_coverage_including_zeros"], 30.0 / 1000000
            )
            self.assertAlmostEqual(result["average_coverage"], 15.0)

        finally:
            os.unlink(profile_file)

    def test_no_matching_positions(self):
        """Test processing when no positions match the filter."""
        profile_data = {
            "contig": ["contig1"],
            "position": [100],  # Not in filter
            "gene_id": ["gene1"],
            "total_coverage": [10.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            init_worker(
                self.metadata_dict,
                self.mag_size_dict,
                self.contig_to_mag,
                self.positions_filter,
                "positions",
            )

            args = ("sample1", profile_file, "MAG001", 0.1, 1.0)
            result = process_positions_mag_sample(args)

            # Should set positions_considered
            self.assertEqual(result["positions_considered"], 3)

            # Should not pass
            self.assertIn("No matching positions", result["breadth_fail_reason"])

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
            init_worker(
                self.metadata_dict,
                {"INVALID_MAG": 0},  # Invalid size
                self.contig_to_mag,
                self.positions_filter,
                "positions",
            )

            args = ("sample1", profile_file, "INVALID_MAG", 0.1, 1.0)
            result = process_positions_mag_sample(args)

            # Should fail with appropriate error
            self.assertIn("MAG size", result["breadth_fail_reason"])

        finally:
            os.unlink(profile_file)

    def test_no_positions_for_mag(self):
        """Test handling when no positions exist for MAG in filter."""
        profile_data = {
            "contig": ["contig1"],
            "position": [1],
            "gene_id": ["gene1"],
            "total_coverage": [10.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            init_worker(
                self.metadata_dict,
                self.mag_size_dict,
                self.contig_to_mag,
                {"MAG999": {("contig1", 1)}},  # Different MAG
                "positions",
            )

            args = ("sample1", profile_file, "MAG001", 0.1, 1.0)
            result = process_positions_mag_sample(args)

            # Should fail with appropriate error
            self.assertIn("No positions found in filter", result["breadth_fail_reason"])

        finally:
            os.unlink(profile_file)

    def test_zero_coverage_warning(self):
        """Test that zero coverage positions trigger warning."""
        profile_data = {
            "contig": ["contig1", "contig1", "contig1"],
            "position": [1, 2, 3],
            "gene_id": ["gene1", "gene2", "gene3"],
            "total_coverage": [10.0, 0.0, 20.0],  # One zero
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            init_worker(
                self.metadata_dict,
                self.mag_size_dict,
                self.contig_to_mag,
                self.positions_filter,
                "positions",
            )

            with patch(
                "alleleflux.scripts.accessory.positions_qc.logger"
            ) as mock_logger:
                args = ("sample1", profile_file, "MAG001", 0.000001, 1.0)
                result = process_positions_mag_sample(args)

                # Should log warning
                mock_logger.warning.assert_called()
                warning_call = str(mock_logger.warning.call_args)
                self.assertIn("zero coverage", warning_call.lower())

        finally:
            os.unlink(profile_file)

    def test_breadth_threshold_genome_check(self):
        """Test that breadth threshold uses genome-wide metric."""
        # 100 positions total in profile, only 10 in filter
        profile_data = {
            "contig": ["contig1"] * 100,
            "position": list(range(1, 101)),
            "gene_id": [f"gene{i}" for i in range(1, 101)],
            "total_coverage": [10.0] * 100,
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            # Filter to only first 10 positions
            positions_filter = {"MAG001": {("contig1", i) for i in range(1, 11)}}

            init_worker(
                self.metadata_dict,
                self.mag_size_dict,
                self.contig_to_mag,
                positions_filter,
                "positions",
            )

            # Set threshold that genome-wide breadth passes
            # breadth_genome = 100/1M = 0.0001
            args = ("sample1", profile_file, "MAG001", 0.00005, 1.0)
            result = process_positions_mag_sample(args)

            # Should pass based on genome-wide breadth
            self.assertTrue(result["breadth_threshold_passed"])
            self.assertAlmostEqual(result["breadth_genome"], 100 / 1000000)

            # Now test with higher threshold that fails
            args = ("sample1", profile_file, "MAG001", 0.0002, 1.0)
            result = process_positions_mag_sample(args)

            self.assertFalse(result["breadth_threshold_passed"])
            self.assertIn("Breadth (genome)", result["breadth_fail_reason"])

        finally:
            os.unlink(profile_file)

    def test_coverage_threshold_genome_check(self):
        """Test that coverage threshold uses genome-wide metric."""
        profile_data = {
            "contig": ["contig1", "contig1"],
            "position": [1, 2],
            "gene_id": ["gene1", "gene2"],
            "total_coverage": [100.0, 200.0],  # High coverage in filtered positions
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            init_worker(
                self.metadata_dict,
                self.mag_size_dict,
                self.contig_to_mag,
                self.positions_filter,
                "positions",
            )

            # average_coverage_genome = 300/1M = 0.0003
            # Set threshold higher than genome-wide but lower than filtered
            args = ("sample1", profile_file, "MAG001", 0.000001, 0.0005)
            result = process_positions_mag_sample(args)

            # Should fail coverage threshold (genome-wide is 0.0003 < 0.0005)
            self.assertFalse(result["coverage_threshold_passed"])
            # Check that the failure reason mentions coverage and threshold
            self.assertIn("coverage", result["coverage_fail_reason"].lower())
            self.assertIn("threshold", result["coverage_fail_reason"].lower())

        finally:
            os.unlink(profile_file)

    def test_metadata_with_time(self):
        """Test that time column is included when present in metadata."""
        profile_data = {
            "contig": ["contig1"],
            "position": [1],
            "gene_id": ["gene1"],
            "total_coverage": [10.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            init_worker(
                self.metadata_dict,
                self.mag_size_dict,
                self.contig_to_mag,
                self.positions_filter,
                "positions",
            )

            args = ("sample1", profile_file, "MAG001", 0.000001, 1.0)
            result = process_positions_mag_sample(args)

            # Should include time
            self.assertIn("time", result)
            self.assertEqual(result["time"], "T1")

        finally:
            os.unlink(profile_file)

    def test_metadata_without_time(self):
        """Test handling when time is not in metadata."""
        # Create metadata without time
        metadata_no_time = {
            "sample1": {
                "group": "treatment",
                "subjectID": "subject1",
                "replicate": "rep1",
            }
        }

        profile_data = {
            "contig": ["contig1"],
            "position": [1],
            "gene_id": ["gene1"],
            "total_coverage": [10.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            init_worker(
                metadata_no_time,
                self.mag_size_dict,
                self.contig_to_mag,
                self.positions_filter,
                "positions",
            )

            args = ("sample1", profile_file, "MAG001", 0.000001, 1.0)
            result = process_positions_mag_sample(args)

            # Should not include time
            self.assertNotIn("time", result)

        finally:
            os.unlink(profile_file)

    def test_median_and_std_calculations(self):
        """Test median and standard deviation calculations."""
        profile_data = {
            "contig": ["contig1", "contig1", "contig2"],
            "position": [1, 2, 1],
            "gene_id": ["gene1", "gene2", "gene3"],
            "total_coverage": [10.0, 20.0, 30.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            init_worker(
                self.metadata_dict,
                self.mag_size_dict,
                self.contig_to_mag,
                self.positions_filter,
                "positions",
            )

            args = ("sample1", profile_file, "MAG001", 0.000001, 1.0)
            result = process_positions_mag_sample(args)

            # All 3 positions are in the filter
            # Median of [10, 20, 30] = 20
            self.assertAlmostEqual(result["median_coverage"], 20.0)

            # Std of [10, 20, 30] with ddof=0
            expected_std = np.std([10.0, 20.0, 30.0], ddof=0)
            self.assertAlmostEqual(result["coverage_std"], expected_std)

            # Should have median_including_zeros
            self.assertIsNotNone(result["median_coverage_including_zeros"])

        finally:
            os.unlink(profile_file)


class TestAggregation(TestPositionsQC):
    """Test aggregation of positions QC results."""

    def test_aggregate_includes_breadth_genome(self):
        """Test that breadth_genome IS included in aggregation for positions QC."""
        df_results = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2"],
                "breadth_threshold_passed": [True, True],
                "group": ["treatment", "treatment"],
                "breadth": [0.8, 0.9],
                "breadth_genome": [0.0001, 0.0002],
                "average_coverage": [50.0, 60.0],
                "average_coverage_genome": [0.00005, 0.00006],
                "median_coverage": [45.0, 55.0],
                "median_coverage_including_zeros": [0.1, 0.2],
                "coverage_std": [15.0, 20.0],
                "coverage_std_including_zeros": [25.0, 30.0],
                "positions_considered": [10, 10],
            }
        )

        result = aggregate_mag_stats(df_results, "MAG001", ["group"], overall=False)

        # For positions QC, breadth_genome and average_coverage_genome ARE aggregated
        # (unlike quality_control.py where they are excluded)
        self.assertIn("breadth_genome_mean", result.columns)
        self.assertIn("average_coverage_genome_mean", result.columns)

        # Check values are correctly averaged
        self.assertAlmostEqual(result.iloc[0]["breadth_genome_mean"], 0.00015)
        self.assertAlmostEqual(result.iloc[0]["average_coverage_genome_mean"], 0.000055)

        # Regular breadth and average_coverage should also be aggregated
        self.assertIn("breadth_mean", result.columns)
        self.assertIn("average_coverage_mean", result.columns)

    def test_aggregate_group_summary(self):
        """Test group-level aggregation."""
        df_results = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2", "sample3"],
                "breadth_threshold_passed": [True, True, True],
                "group": ["treatment", "treatment", "control"],
                "breadth": [0.8, 0.9, 0.7],
                "breadth_genome": [0.0001, 0.0002, 0.00015],
                "average_coverage": [50.0, 60.0, 40.0],
                "average_coverage_genome": [0.00005, 0.00006, 0.00004],
                "median_coverage": [45.0, 55.0, 35.0],
                "median_coverage_including_zeros": [0.1, 0.2, 0.05],
                "coverage_std": [15.0, 20.0, 12.0],
                "coverage_std_including_zeros": [25.0, 30.0, 20.0],
                "positions_considered": [10, 10, 10],
            }
        )

        result = aggregate_mag_stats(df_results, "MAG001", ["group"], overall=False)

        # Should have 2 rows (treatment and control)
        self.assertEqual(len(result), 2)
        self.assertTrue(all(result["MAG_ID"] == "MAG001"))

        # Check treatment group means
        treatment_row = result[result["group"] == "treatment"].iloc[0]
        self.assertAlmostEqual(treatment_row["breadth_mean"], 0.85)
        self.assertAlmostEqual(treatment_row["average_coverage_mean"], 55.0)

    def test_aggregate_overall_summary(self):
        """Test overall aggregation."""
        df_results = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2"],
                "breadth_threshold_passed": [True, True],
                "breadth": [0.8, 0.6],
                "breadth_genome": [0.0001, 0.0002],
                "average_coverage": [50.0, 30.0],
                "average_coverage_genome": [0.00005, 0.00003],
                "median_coverage": [45.0, 25.0],
                "median_coverage_including_zeros": [0.1, 0.03],
                "coverage_std": [15.0, 8.0],
                "coverage_std_including_zeros": [25.0, 15.0],
                "positions_considered": [10, 10],
            }
        )

        result = aggregate_mag_stats(df_results, "MAG001", [], overall=True)

        # Should have 1 row
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["MAG_ID"], "MAG001")
        self.assertEqual(result.iloc[0]["num_samples"], 2)
        self.assertAlmostEqual(result.iloc[0]["breadth_mean"], 0.7)
        self.assertAlmostEqual(result.iloc[0]["average_coverage_mean"], 40.0)

    def test_aggregate_with_time_column(self):
        """Test aggregation when time column is present."""
        df_results = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2", "sample3", "sample4"],
                "breadth_threshold_passed": [True, True, True, True],
                "group": ["treatment", "treatment", "control", "control"],
                "time": ["T1", "T2", "T1", "T2"],
                "breadth": [0.8, 0.9, 0.7, 0.6],
                "breadth_genome": [0.0001, 0.0002, 0.00015, 0.00012],
                "average_coverage": [50.0, 60.0, 40.0, 30.0],
                "average_coverage_genome": [0.00005, 0.00006, 0.00004, 0.00003],
                "median_coverage": [45.0, 55.0, 35.0, 25.0],
                "median_coverage_including_zeros": [0.1, 0.2, 0.05, 0.03],
                "coverage_std": [15.0, 20.0, 12.0, 8.0],
                "coverage_std_including_zeros": [25.0, 30.0, 20.0, 15.0],
                "positions_considered": [10, 10, 10, 10],
            }
        )

        # Test group-time aggregation
        result = aggregate_mag_stats(
            df_results, "MAG001", ["group", "time"], overall=False
        )

        # Should have 4 rows (2 groups × 2 timepoints)
        self.assertEqual(len(result), 4)
        self.assertIn("group", result.columns)
        self.assertIn("time", result.columns)


class TestEdgeCases(TestPositionsQC):
    """Test edge cases and error handling."""

    def test_empty_results_aggregation(self):
        """Test aggregation with empty results."""
        df_results = pd.DataFrame()

        result = aggregate_mag_stats(df_results, "MAG001", ["group"], overall=False)

        # Should return None for empty input
        self.assertIsNone(result)

    def test_partial_coverage_in_positions(self):
        """Test with some positions having zero coverage."""
        profile_data = {
            "contig": ["contig1", "contig1", "contig2"],
            "position": [1, 2, 1],
            "gene_id": ["gene1", "gene2", "gene3"],
            "total_coverage": [10.0, 0.0, 20.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            init_worker(
                self.metadata_dict,
                self.mag_size_dict,
                self.contig_to_mag,
                self.positions_filter,
                "positions",
            )

            args = ("sample1", profile_file, "MAG001", 0.000001, 1.0)
            result = process_positions_mag_sample(args)

            # Breadth: 2 positions with coverage >= 1 / 3 total
            self.assertAlmostEqual(result["breadth"], 2.0 / 3)

            # Breadth_genome: 2 positions / 1M
            self.assertAlmostEqual(result["breadth_genome"], 2.0 / 1000000)

            # Median should exclude zero
            self.assertAlmostEqual(result["median_coverage"], 15.0)

        finally:
            os.unlink(profile_file)

    def test_all_positions_filtered_out(self):
        """Test when all positions are filtered by MAG assignment."""
        profile_data = {
            "contig": ["contig3"],  # Belongs to MAG002
            "position": [1],
            "gene_id": ["gene1"],
            "total_coverage": [10.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            init_worker(
                self.metadata_dict,
                self.mag_size_dict,
                self.contig_to_mag,
                self.positions_filter,
                "positions",
            )

            args = ("sample1", profile_file, "MAG001", 0.1, 1.0)

            with patch(
                "alleleflux.scripts.accessory.positions_qc.logger"
            ) as mock_logger:
                result = process_positions_mag_sample(args)

                # Should fail with appropriate error
                self.assertIn(
                    "No rows remain after filtering", result["breadth_fail_reason"]
                )

        finally:
            os.unlink(profile_file)


class TestGenomeWithoutSamples(unittest.TestCase):
    """A metadata folder may hold header-only files, and the tool copes.

    The metadata step writes one ``{mag}_metadata.tsv`` per genome of the run
    in every timepoint's folder: header-only when the genome has no sample at
    that timepoint.  Such a genome has nothing to report at these positions,
    so it is skipped with a warning and the other genomes are reported.  A
    folder in which NO genome has a sample is still an error: there is
    nothing to write, and silence would look like a finished run.

    Layout: genome MAG_A is one 100-base contig; sample s1 covers all 100
    positions at depth 10, sample s2 covers 5.  Genome MAG_B has a header-only
    metadata file.  The positions file names both genomes.
    """

    HEADER = "sample_id\tsubjectID\ttime\tgroup\treplicate\tfile_path\tbam_path\n"

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.meta = self.tmp / "inputMetadata_pre"
        self.out = self.tmp / "positions_qc"
        self.meta.mkdir()
        (self.tmp / "ref.fa").write_text(
            ">c1\n" + "A" * 100 + "\n>d1\n" + "C" * 100 + "\n"
        )
        (self.tmp / "mapping.tsv").write_text(
            "mag_id\tcontig_id\nMAG_A\tc1\nMAG_B\td1\n"
        )
        (self.tmp / "positions.tsv").write_text(
            "MAG\tcontig\tposition\nMAG_A\tc1\t1\nMAG_A\tc1\t2\nMAG_B\td1\t1\n"
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
        (self.meta / "MAG_A_metadata.tsv").write_text(self.HEADER + "".join(rows))
        # MAG_B: profiled somewhere in the run, but no sample at this timepoint.
        (self.meta / "MAG_B_metadata.tsv").write_text(self.HEADER)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self):
        """The real command on the fixture folder (``python -m``, so the test
        runs the code of this checkout whatever is on PATH)."""
        return subprocess.run(
            [
                sys.executable,
                "-m",
                "alleleflux.scripts.accessory.positions_qc",
                "--metadata_dir",
                str(self.meta),
                "--fasta",
                str(self.tmp / "ref.fa"),
                "--mag_mapping_file",
                str(self.tmp / "mapping.tsv"),
                "--positions_file",
                str(self.tmp / "positions.tsv"),
                "--output_dir",
                str(self.out),
                "--cpus",
                "2",
            ],
            capture_output=True,
            text=True,
        )

    def test_process_mag_has_nothing_to_report_for_a_genome_without_samples(self):
        self.out.mkdir()
        summaries = process_mag(
            (
                "MAG_B",
                self.meta / "MAG_B_metadata.tsv",
                {"MAG_A": 100, "MAG_B": 100},
                {"c1": "MAG_A", "d1": "MAG_B"},
                {"MAG_B": {("d1", 1)}},
                self.out,
                0.1,
                1.0,
                2,
                "positions",
                "longitudinal",
            )
        )
        # Same shape as a genome with samples (three summaries), all absent.
        self.assertEqual(summaries, (None, None, None))
        # No table for a genome that was not read.
        self.assertEqual(list(self.out.iterdir()), [])

    def test_genome_without_samples_is_skipped_and_the_others_are_reported(self):
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        # The genome with samples: one row per sample.
        table = pd.read_csv(self.out / "MAG_A_positions_QC.tsv", sep="\t")
        self.assertEqual(sorted(table["sample_id"]), ["s1", "s2"])
        # The genome without samples: no table, absent from the summary, and
        # named in the log so the reader knows why.
        self.assertFalse((self.out / "MAG_B_positions_QC.tsv").exists())
        overall = pd.read_csv(
            self.out / "ALL_MAGs_positions_QC_overall_summary.tsv", sep="\t"
        )
        self.assertEqual(list(overall["MAG_ID"]), ["MAG_A"])
        self.assertIn("MAG_B", done.stderr + done.stdout)
        self.assertIn("skipped", done.stderr + done.stdout)

    def test_folder_where_no_genome_has_a_sample_is_an_error(self):
        # Every metadata file header-only: nothing to report at all.
        (self.meta / "MAG_A_metadata.tsv").write_text(self.HEADER)
        done = self._run()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("no sample in any of the 2 metadata file(s)", done.stderr)
        # And no summary file that would look like a result.
        self.assertFalse(
            (self.out / "ALL_MAGs_positions_QC_overall_summary.tsv").exists()
        )

    def test_metadata_file_missing_a_column_is_still_an_error(self):
        # Guard: only "no rows" is tolerated.  A file that is not a metadata
        # table (here: no file_path column) must keep stopping the tool.
        (self.meta / "MAG_B_metadata.tsv").write_text(
            "sample_id\tsubjectID\ttime\tgroup\treplicate\n"
        )
        done = self._run()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("file_path", done.stderr)


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
