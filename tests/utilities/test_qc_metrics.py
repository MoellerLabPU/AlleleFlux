#!/usr/bin/env python3

import os
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd

from alleleflux.scripts.utilities.qc_metrics import (
    aggregate_mag_stats,
    calculate_breadth_metrics,
    calculate_coverage_metrics,
    calculate_length_weighted_coverage,
    calculate_median_and_std_metrics,
    check_breadth_threshold,
    check_coverage_threshold,
    load_and_validate_profile,
)


class TestQCMetrics(unittest.TestCase):
    """Base test class for QC metrics tests."""

    def setUp(self):
        """Set up test fixtures."""
        self.contig_to_mag = {
            "contig1": "MAG001",
            "contig2": "MAG001",
            "contig3": "MAG002",
        }

        self.contig_length_dict = {
            "contig1": 500000,
            "contig2": 300000,
            "contig3": 200000,
        }

    def create_temp_profile(self, data):
        """Create a temporary profile file with given data."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tsv")
        df = pd.DataFrame(data)
        df.to_csv(temp_file.name, sep="\t", index=False)
        temp_file.close()
        return temp_file.name


class TestLoadAndValidateProfile(TestQCMetrics):
    """Test load_and_validate_profile function."""

    def test_load_valid_profile(self):
        """Test loading a valid profile."""
        profile_data = {
            "contig": ["contig1", "contig1", "contig2"],
            "position": [1, 2, 1],
            "gene_id": ["gene1", "gene2", "gene3"],
            "total_coverage": [10.0, 20.0, 30.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            df = load_and_validate_profile(
                profile_file, "MAG001", "sample1", self.contig_to_mag
            )

            # Should only include contig1 and contig2 (MAG001)
            self.assertEqual(len(df), 3)
            self.assertTrue(all(df["contig"].isin(["contig1", "contig2"])))

        finally:
            os.unlink(profile_file)

    def test_filter_wrong_mag_contigs(self):
        """Test that contigs from other MAGs are filtered out."""
        profile_data = {
            "contig": ["contig1", "contig2", "contig3"],
            "position": [1, 1, 1],
            "gene_id": ["gene1", "gene2", "gene3"],
            "total_coverage": [10.0, 20.0, 30.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
                df = load_and_validate_profile(
                    profile_file, "MAG001", "sample1", self.contig_to_mag
                )

                # Should only include contig1 and contig2
                self.assertEqual(len(df), 2)
                self.assertNotIn("contig3", df["contig"].values)

                # Should log warning about filtering
                mock_logger.warning.assert_called()

        finally:
            os.unlink(profile_file)

    def test_empty_after_filtering(self):
        """Test error when no rows remain after filtering."""
        profile_data = {
            "contig": ["contig3"],  # Only MAG002 contigs
            "position": [1],
            "gene_id": ["gene1"],
            "total_coverage": [10.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            with self.assertRaises(ValueError) as context:
                load_and_validate_profile(
                    profile_file, "MAG001", "sample1", self.contig_to_mag
                )

            self.assertIn("No rows remain", str(context.exception))

        finally:
            os.unlink(profile_file)

    def test_gene_id_as_string(self):
        """Test that gene_id is preserved as string."""
        profile_data = {
            "contig": ["contig1"],
            "position": [1],
            "gene_id": ["001"],  # Could be parsed as int
            "total_coverage": [10.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            df = load_and_validate_profile(
                profile_file, "MAG001", "sample1", self.contig_to_mag
            )

            # gene_id should be string, not int
            self.assertEqual(df["gene_id"].dtype, object)
            self.assertEqual(df["gene_id"].iloc[0], "001")

        finally:
            os.unlink(profile_file)

    def test_profile_without_contig_column(self):
        """Test handling when contig column is missing."""
        profile_data = {
            "position": [1, 2],
            "gene_id": ["gene1", "gene2"],
            "total_coverage": [10.0, 20.0],
        }

        profile_file = self.create_temp_profile(profile_data)

        try:
            # Should return dataframe as-is without filtering
            df = load_and_validate_profile(
                profile_file, "MAG001", "sample1", self.contig_to_mag
            )

            self.assertEqual(len(df), 2)
            self.assertNotIn("contig", df.columns)

        finally:
            os.unlink(profile_file)


class TestCalculateBreadthMetrics(TestQCMetrics):
    """Test calculate_breadth_metrics function."""

    def test_breadth_no_filter(self):
        """Test breadth calculation without filter."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0, 0.0, 30.0]})

        result = calculate_breadth_metrics(
            df, denom=1000000, mag_size=1000000, using_filter=False
        )

        # 3 positions with coverage >= 1 out of 1M
        self.assertAlmostEqual(result["breadth"], 3 / 1000000)
        self.assertIsNone(result["breadth_genome"])

    def test_breadth_with_filter(self):
        """Test breadth calculation with filter."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0, 30.0]})

        result = calculate_breadth_metrics(
            df,
            denom=10,  # Filtered to 10 positions
            mag_size=1000000,
            using_filter=True,
            positions_with_coverage_genome=50,  # 50 positions genome-wide
        )

        # Filtered: 3 positions with coverage / 10 positions = 0.3
        self.assertAlmostEqual(result["breadth"], 3 / 10)
        # Genome-wide: 50 positions / 1M = 0.00005
        self.assertAlmostEqual(result["breadth_genome"], 50 / 1000000)

    def test_breadth_with_zeros(self):
        """Test breadth calculation with zero coverage positions."""
        df = pd.DataFrame({"total_coverage": [10.0, 0.0, 0.0, 30.0]})

        result = calculate_breadth_metrics(
            df, denom=100, mag_size=100, using_filter=False
        )

        # Only 2 positions with coverage >= 1
        self.assertAlmostEqual(result["breadth"], 2 / 100)

    def test_breadth_zero_denominator(self):
        """Test breadth with zero denominator returns NaN."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0]})

        result = calculate_breadth_metrics(
            df, denom=0, mag_size=1000000, using_filter=False
        )

        self.assertTrue(np.isnan(result["breadth"]))

    def test_breadth_missing_genome_param(self):
        """Test error when using_filter=True but genome param missing."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0]})

        with self.assertRaises(ValueError) as context:
            calculate_breadth_metrics(
                df,
                denom=10,
                mag_size=1000000,
                using_filter=True,
                positions_with_coverage_genome=None,  # Missing!
            )

        self.assertIn("positions_with_coverage_genome", str(context.exception))


class TestCalculateCoverageMetrics(TestQCMetrics):
    """Test calculate_coverage_metrics function."""

    def test_coverage_no_filter(self):
        """Test coverage calculation without filter."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0, 30.0]})

        result = calculate_coverage_metrics(
            df, denom=1000000, mag_size=1000000, using_filter=False
        )

        # Sum = 60, denom = 1M
        self.assertAlmostEqual(result["average_coverage"], 60.0 / 1000000)
        self.assertIsNone(result["average_coverage_genome"])

    def test_coverage_with_filter(self):
        """Test coverage calculation with filter."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0, 30.0]})

        result = calculate_coverage_metrics(
            df,
            denom=10,  # Filtered denominator
            mag_size=1000000,
            using_filter=True,
            total_coverage_sum_genome=1000.0,  # Genome-wide sum
        )

        # Filtered: 60 / 10 = 6.0
        self.assertAlmostEqual(result["average_coverage"], 6.0)
        # Genome-wide: 1000 / 1M = 0.001
        self.assertAlmostEqual(result["average_coverage_genome"], 0.001)

    def test_coverage_zero_denominator(self):
        """Test coverage with zero denominator returns NaN."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0]})

        result = calculate_coverage_metrics(
            df, denom=0, mag_size=1000000, using_filter=False
        )

        self.assertTrue(np.isnan(result["average_coverage"]))

    def test_coverage_missing_genome_param(self):
        """Test error when using_filter=True but genome param missing."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0]})

        with self.assertRaises(ValueError) as context:
            calculate_coverage_metrics(
                df,
                denom=10,
                mag_size=1000000,
                using_filter=True,
                total_coverage_sum_genome=None,  # Missing!
            )

        self.assertIn("total_coverage_sum_genome", str(context.exception))


class TestCalculateMedianAndStdMetrics(TestQCMetrics):
    """Test calculate_median_and_std_metrics function."""

    def test_median_and_std_basic(self):
        """Test basic median and std calculations."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0, 30.0]})

        result = calculate_median_and_std_metrics(df, denom=1000, avg_cov=0.06)

        # Median of [10, 20, 30] = 20
        self.assertAlmostEqual(result["median_coverage"], 20.0)

        # Std with ddof=0
        expected_std = np.std([10.0, 20.0, 30.0], ddof=0)
        self.assertAlmostEqual(result["coverage_std"], expected_std)

    def test_median_including_zeros(self):
        """Test median calculation including zeros for absent positions."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0, 30.0]})

        result = calculate_median_and_std_metrics(df, denom=10, avg_cov=6.0)

        # 3 observed + 7 absent (zeros) = [10, 20, 30, 0, 0, 0, 0, 0, 0, 0]
        # Median = 5.0
        expected_median = np.median([10.0, 20.0, 30.0] + [0.0] * 7)
        self.assertAlmostEqual(
            result["median_coverage_including_zeros"], expected_median
        )

    def test_std_including_zeros(self):
        """Test std calculation including zeros using variance formula."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0, 30.0]})
        avg_cov = 6.0  # (10+20+30) / 10 = 6.0
        denom = 10

        result = calculate_median_and_std_metrics(df, denom, avg_cov)

        # Var(X) = E[X²] - E[X]²
        # E[X²] = (100 + 400 + 900) / 10 = 140
        # E[X]² = 6.0² = 36
        # Var = 140 - 36 = 104
        # Std = sqrt(104) ≈ 10.198
        ex2 = (100 + 400 + 900) / 10
        var = ex2 - (avg_cov**2)
        expected_std = np.sqrt(var)
        self.assertAlmostEqual(result["coverage_std_including_zeros"], expected_std)

    def test_with_zero_coverage(self):
        """Test with some zero coverage values."""
        df = pd.DataFrame({"total_coverage": [10.0, 0.0, 20.0, 0.0, 30.0]})

        result = calculate_median_and_std_metrics(df, denom=100, avg_cov=0.6)

        # Median of non-zero: [10, 20, 30] = 20
        self.assertAlmostEqual(result["median_coverage"], 20.0)

        # Median including zeros: 5 observed + 95 absent
        all_values = [10.0, 0.0, 20.0, 0.0, 30.0] + [0.0] * 95
        expected_median = np.median(all_values)
        self.assertAlmostEqual(
            result["median_coverage_including_zeros"], expected_median
        )

    def test_empty_dataframe(self):
        """Test with empty dataframe."""
        df = pd.DataFrame({"total_coverage": []})

        result = calculate_median_and_std_metrics(df, denom=100, avg_cov=0.0)

        # All should be NaN except median_including_zeros (which is 0)
        self.assertTrue(np.isnan(result["median_coverage"]))
        self.assertTrue(np.isnan(result["coverage_std"]))
        self.assertAlmostEqual(result["median_coverage_including_zeros"], 0.0)

    def test_all_zeros(self):
        """Test with all zero coverage."""
        df = pd.DataFrame({"total_coverage": [0.0, 0.0, 0.0]})

        result = calculate_median_and_std_metrics(df, denom=100, avg_cov=0.0)

        # No positions with coverage > 0
        self.assertTrue(np.isnan(result["median_coverage"]))
        self.assertTrue(np.isnan(result["coverage_std"]))

        # Median including zeros should be 0.0
        self.assertAlmostEqual(result["median_coverage_including_zeros"], 0.0)

    def test_more_positions_than_denominator(self):
        """Test handling when profile has more positions than denominator."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0, 30.0, 40.0, 50.0]})

        with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
            result = calculate_median_and_std_metrics(df, denom=3, avg_cov=50.0)

            # Should log warning
            mock_logger.warning.assert_called()

            # Should fall back to median_coverage (no zeros added)
            self.assertAlmostEqual(result["median_coverage_including_zeros"], 30.0)


class TestCalculateLengthWeightedCoverage(TestQCMetrics):
    """Test calculate_length_weighted_coverage function."""

    def test_length_weighted_basic(self):
        """Test basic length-weighted coverage calculation."""
        df = pd.DataFrame(
            {
                "contig": ["contig1", "contig1", "contig2", "contig2"],
                "total_coverage": [100.0, 200.0, 300.0, 150.0],
            }
        )

        # contig1: sum=300, length=500000, mean=300/500000=0.0006
        # contig2: sum=450, length=300000, mean=450/300000=0.0015
        # weighted = (0.0006*500000 + 0.0015*300000) / 800000
        # weighted = (300 + 450) / 800000 = 750/800000 = 0.0009375
        mag_size = 800000

        result = calculate_length_weighted_coverage(
            df, self.contig_length_dict, "MAG001", "sample1", mag_size
        )

        expected = 750.0 / mag_size
        self.assertAlmostEqual(result, expected)

    def test_missing_contig_column(self):
        """Test handling when contig column is missing."""
        df = pd.DataFrame({"total_coverage": [10.0, 20.0]})

        with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
            result = calculate_length_weighted_coverage(
                df, self.contig_length_dict, "MAG001", "sample1", 1000000
            )

            self.assertIsNone(result)
            mock_logger.warning.assert_called()

    def test_invalid_mag_size(self):
        """Test handling with invalid MAG size."""
        df = pd.DataFrame({"contig": ["contig1"], "total_coverage": [10.0]})

        with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
            result = calculate_length_weighted_coverage(
                df, self.contig_length_dict, "MAG001", "sample1", 0
            )

            self.assertIsNone(result)
            mock_logger.warning.assert_called()

    def test_contig_without_length(self):
        """Test handling contigs without valid lengths."""
        df = pd.DataFrame(
            {
                "contig": ["contig1", "contig2", "contig_unknown"],
                "total_coverage": [100.0, 200.0, 300.0],
            }
        )

        with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
            result = calculate_length_weighted_coverage(
                df, self.contig_length_dict, "MAG001", "sample1", 800000
            )

            # Should drop contig_unknown and calculate with contig1 and contig2
            self.assertIsNotNone(result)
            mock_logger.warning.assert_called()

    def test_empty_after_filtering(self):
        """Test when no contigs remain after filtering invalid lengths."""
        df = pd.DataFrame({"contig": ["contig_unknown"], "total_coverage": [100.0]})

        result = calculate_length_weighted_coverage(
            df, self.contig_length_dict, "MAG001", "sample1", 1000000
        )

        self.assertIsNone(result)

    def test_empty_input(self):
        """Test with empty dataframe."""
        df = pd.DataFrame({"contig": [], "total_coverage": []})

        result = calculate_length_weighted_coverage(
            df, self.contig_length_dict, "MAG001", "sample1", 1000000
        )

        self.assertIsNone(result)


class TestCheckBreadthThreshold(TestQCMetrics):
    """Test check_breadth_threshold function."""

    def test_breadth_pass_no_filter(self):
        """Test breadth passes threshold without filter."""
        passed, reason = check_breadth_threshold(
            breadth=0.15,
            breadth_genome=None,
            threshold=0.1,
            using_filter=False,
            sample_id="sample1",
            mag_id="MAG001",
        )

        self.assertTrue(passed)
        self.assertEqual(reason, "")

    def test_breadth_fail_no_filter(self):
        """Test breadth fails threshold without filter."""
        with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
            passed, reason = check_breadth_threshold(
                breadth=0.05,
                breadth_genome=None,
                threshold=0.1,
                using_filter=False,
                sample_id="sample1",
                mag_id="MAG001",
            )

            self.assertFalse(passed)
            self.assertIn("Breadth", reason)
            self.assertIn("5.00%", reason)
            self.assertIn("10.00%", reason)
            mock_logger.info.assert_called()

    def test_breadth_pass_with_filter(self):
        """Test breadth passes when using filter (checks genome-wide)."""
        passed, reason = check_breadth_threshold(
            breadth=0.8,  # Filtered breadth (high)
            breadth_genome=0.15,  # Genome-wide breadth
            threshold=0.1,
            using_filter=True,
            sample_id="sample1",
            mag_id="MAG001",
        )

        # Should check breadth_genome (0.15 > 0.1), so passes
        self.assertTrue(passed)
        self.assertEqual(reason, "")

    def test_breadth_fail_with_filter(self):
        """Test breadth fails when using filter (checks genome-wide)."""
        with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
            passed, reason = check_breadth_threshold(
                breadth=0.9,  # Filtered breadth (high, but irrelevant)
                breadth_genome=0.05,  # Genome-wide breadth (low)
                threshold=0.1,
                using_filter=True,
                sample_id="sample1",
                mag_id="MAG001",
            )

            # Should check breadth_genome (0.05 < 0.1), so fails
            self.assertFalse(passed)
            self.assertIn("Breadth (genome)", reason)
            self.assertIn("5.00%", reason)
            mock_logger.info.assert_called()


class TestCheckCoverageThreshold(TestQCMetrics):
    """Test check_coverage_threshold function."""

    def test_coverage_pass_no_filter(self):
        """Test coverage passes threshold without filter."""
        passed, reason = check_coverage_threshold(
            avg_cov=15.0,
            avg_cov_genome=None,
            threshold=10.0,
            breadth_passed=True,
            using_filter=False,
            sample_id="sample1",
            mag_id="MAG001",
        )

        self.assertTrue(passed)
        self.assertEqual(reason, "")

    def test_coverage_fail_no_filter(self):
        """Test coverage fails threshold without filter."""
        with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
            passed, reason = check_coverage_threshold(
                avg_cov=5.0,
                avg_cov_genome=None,
                threshold=10.0,
                breadth_passed=True,
                using_filter=False,
                sample_id="sample1",
                mag_id="MAG001",
            )

            self.assertFalse(passed)
            self.assertIn("Average coverage", reason)
            self.assertIn("5.00x", reason)
            self.assertIn("10.00x", reason)
            mock_logger.info.assert_called()

    def test_coverage_cascading_failure(self):
        """Test coverage fails when breadth check failed (cascading)."""
        passed, reason = check_coverage_threshold(
            avg_cov=50.0,  # High coverage (irrelevant)
            avg_cov_genome=None,
            threshold=10.0,
            breadth_passed=False,  # Breadth failed
            using_filter=False,
            sample_id="sample1",
            mag_id="MAG001",
        )

        self.assertFalse(passed)
        self.assertEqual(reason, "Breadth check failed")

    def test_coverage_pass_with_filter(self):
        """Test coverage passes when using filter (checks genome-wide)."""
        passed, reason = check_coverage_threshold(
            avg_cov=50.0,  # Filtered coverage
            avg_cov_genome=15.0,  # Genome-wide coverage
            threshold=10.0,
            breadth_passed=True,
            using_filter=True,
            sample_id="sample1",
            mag_id="MAG001",
        )

        # Should check avg_cov_genome (15.0 > 10.0), so passes
        self.assertTrue(passed)
        self.assertEqual(reason, "")

    def test_coverage_fail_with_filter(self):
        """Test coverage fails when using filter (checks genome-wide)."""
        with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
            passed, reason = check_coverage_threshold(
                avg_cov=100.0,  # Filtered coverage (high, but irrelevant)
                avg_cov_genome=5.0,  # Genome-wide coverage (low)
                threshold=10.0,
                breadth_passed=True,
                using_filter=True,
                sample_id="sample1",
                mag_id="MAG001",
            )

            # Should check avg_cov_genome (5.0 < 10.0), so fails
            self.assertFalse(passed)
            self.assertIn("Average coverage (genome)", reason)
            self.assertIn("5.00x", reason)
            mock_logger.info.assert_called()


class TestAggregateMagStats(TestQCMetrics):
    """Test aggregate_mag_stats function."""

    def test_aggregate_by_group(self):
        """Test aggregation by group."""
        df = pd.DataFrame(
            {
                "sample_id": ["s1", "s2", "s3", "s4"],
                "coverage_threshold_passed": [True, True, True, True],
                "group": ["treatment", "treatment", "control", "control"],
                "breadth": [0.8, 0.9, 0.7, 0.6],
                "average_coverage": [50.0, 60.0, 40.0, 30.0],
                "median_coverage": [45.0, 55.0, 35.0, 25.0],
            }
        )

        result = aggregate_mag_stats(df, "MAG001", ["group"], overall=False)

        # Should have 2 rows
        self.assertEqual(len(result), 2)
        self.assertTrue(all(result["MAG_ID"] == "MAG001"))

        # Check treatment group
        treatment = result[result["group"] == "treatment"].iloc[0]
        self.assertAlmostEqual(treatment["breadth_mean"], 0.85)
        self.assertAlmostEqual(treatment["average_coverage_mean"], 55.0)
        self.assertEqual(treatment["num_samples"], 2)

    def test_aggregate_overall(self):
        """Test overall aggregation."""
        df = pd.DataFrame(
            {
                "sample_id": ["s1", "s2", "s3"],
                "coverage_threshold_passed": [True, True, True],
                "breadth": [0.8, 0.9, 0.7],
                "average_coverage": [50.0, 60.0, 40.0],
            }
        )

        result = aggregate_mag_stats(df, "MAG001", [], overall=True)

        # Should have 1 row
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["MAG_ID"], "MAG001")
        self.assertAlmostEqual(result.iloc[0]["breadth_mean"], 0.8)
        self.assertAlmostEqual(result.iloc[0]["average_coverage_mean"], 50.0)
        self.assertEqual(result.iloc[0]["num_samples"], 3)

    def test_aggregate_filters_passing_samples(self):
        """Test that only passing samples are aggregated."""
        df = pd.DataFrame(
            {
                "sample_id": ["s1", "s2", "s3"],
                "coverage_threshold_passed": [True, True, False],
                "breadth": [0.8, 0.9, 0.3],  # s3 should be excluded
                "average_coverage": [50.0, 60.0, 10.0],
            }
        )

        result = aggregate_mag_stats(df, "MAG001", [], overall=True)

        # Should only include s1 and s2
        self.assertEqual(result.iloc[0]["num_samples"], 2)
        self.assertAlmostEqual(result.iloc[0]["breadth_mean"], 0.85)

    def test_aggregate_fallback_all_samples(self):
        """Test fallback to all samples when none pass."""
        df = pd.DataFrame(
            {
                "sample_id": ["s1", "s2"],
                "coverage_threshold_passed": [False, False],
                "breadth": [0.3, 0.4],
                "average_coverage": [5.0, 6.0],
            }
        )

        with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
            result = aggregate_mag_stats(df, "MAG001", [], overall=True)

            # Should include all samples
            self.assertEqual(result.iloc[0]["num_samples"], 2)
            mock_logger.warning.assert_called()

    def test_aggregate_empty_dataframe(self):
        """Test aggregation with empty dataframe."""
        df = pd.DataFrame()

        with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
            result = aggregate_mag_stats(df, "MAG001", ["group"], overall=False)

            self.assertIsNone(result)
            mock_logger.warning.assert_called()

    def test_aggregate_no_threshold_column(self):
        """Test aggregation when coverage_threshold_passed column missing."""
        df = pd.DataFrame(
            {
                "sample_id": ["s1", "s2"],
                "breadth": [0.8, 0.9],
                "average_coverage": [50.0, 60.0],
            }
        )

        result = aggregate_mag_stats(df, "MAG001", [], overall=True)

        # Should use all samples
        self.assertEqual(result.iloc[0]["num_samples"], 2)

    def test_aggregate_with_breadth_genome(self):
        """Test that breadth_genome is included in aggregation."""
        df = pd.DataFrame(
            {
                "sample_id": ["s1", "s2"],
                "coverage_threshold_passed": [True, True],
                "breadth": [0.8, 0.9],
                "breadth_genome": [0.0001, 0.0002],
                "average_coverage": [50.0, 60.0],
                "average_coverage_genome": [0.00005, 0.00006],
            }
        )

        result = aggregate_mag_stats(df, "MAG001", [], overall=True)

        # Should include both regular and genome metrics
        self.assertIn("breadth_mean", result.columns)
        self.assertIn("breadth_genome_mean", result.columns)
        self.assertIn("average_coverage_mean", result.columns)
        self.assertIn("average_coverage_genome_mean", result.columns)

        # Check values
        self.assertAlmostEqual(result.iloc[0]["breadth_genome_mean"], 0.00015)
        self.assertAlmostEqual(result.iloc[0]["average_coverage_genome_mean"], 0.000055)

    def test_aggregate_group_and_time(self):
        """Test aggregation by multiple grouping columns."""
        df = pd.DataFrame(
            {
                "sample_id": ["s1", "s2", "s3", "s4"],
                "coverage_threshold_passed": [True, True, True, True],
                "group": ["treatment", "treatment", "control", "control"],
                "time": ["T1", "T2", "T1", "T2"],
                "breadth": [0.8, 0.9, 0.7, 0.6],
            }
        )

        result = aggregate_mag_stats(df, "MAG001", ["group", "time"], overall=False)

        # Should have 4 rows (2 groups × 2 timepoints)
        self.assertEqual(len(result), 4)
        self.assertIn("group", result.columns)
        self.assertIn("time", result.columns)

    def test_aggregate_empty_group_cols_not_overall(self):
        """Test error when group_cols is empty but overall=False."""
        df = pd.DataFrame(
            {
                "sample_id": ["s1", "s2"],
                "coverage_threshold_passed": [True, True],
                "breadth": [0.8, 0.9],
            }
        )

        with patch("alleleflux.scripts.utilities.qc_metrics.logger") as mock_logger:
            result = aggregate_mag_stats(df, "MAG001", [], overall=False)

            self.assertIsNone(result)
            mock_logger.error.assert_called()

    def test_aggregate_all_metric_types(self):
        """Test that all metric types are correctly aggregated."""
        df = pd.DataFrame(
            {
                "sample_id": ["s1", "s2"],
                "coverage_threshold_passed": [True, True],
                "breadth": [0.8, 0.9],
                "average_coverage": [50.0, 60.0],
                "median_coverage": [45.0, 55.0],
                "median_coverage_including_zeros": [0.1, 0.2],
                "coverage_std": [15.0, 20.0],
                "coverage_std_including_zeros": [25.0, 30.0],
                "length_weighted_coverage": [48.0, 58.0],
                "subjects_per_group": [5, 5],
                "replicates_per_group": [3, 3],
                "paired_replicates_per_group": [2, 2],
            }
        )

        result = aggregate_mag_stats(df, "MAG001", [], overall=True)

        # All metrics should be present with _mean suffix
        expected_cols = [
            "breadth_mean",
            "average_coverage_mean",
            "median_coverage_mean",
            "median_coverage_including_zeros_mean",
            "coverage_std_mean",
            "coverage_std_including_zeros_mean",
            "length_weighted_coverage_mean",
            "subjects_per_group_mean",
            "replicates_per_group_mean",
            "paired_replicates_per_group_mean",
        ]

        for col in expected_cols:
            self.assertIn(col, result.columns)


if __name__ == "__main__":
    # Create a test suite
    suite = unittest.TestLoader().loadTestsFromModule(__import__(__name__))

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
