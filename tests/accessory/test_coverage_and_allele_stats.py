import numpy as np
import pandas as pd
import pytest

from alleleflux.scripts.accessory.coverage_and_allele_stats import (
    ALLELES,
    _aggregate_with_freqs,
    _std_from_moments,
    compute_allele_frequencies,
    compute_grouped_statistics,
    compute_overall_statistics,
)


def make_df(rows):
    return pd.DataFrame(rows)


def test_std_from_moments_basic_and_edge():
    # More complex test data with 4 positions
    sum_ = pd.Series([10.0, 6.0, 25.0, 0.0])
    sumsq = pd.Series([110.0, 20.0, 175.0, 0.0])
    denom = pd.Series([2.0, 2.0, 5.0, 3.0])
    mean, std = _std_from_moments(sum_, sumsq, denom)

    # Expected values: mean = sum/denom, var = sumsq/denom - mean^2, std = sqrt(var)
    expected_mean = np.array([5.0, 3.0, 5.0, 0.0])
    expected_var = np.array([110.0 / 2 - 25.0, 20.0 / 2 - 9.0, 175.0 / 5 - 25.0, 0.0])
    expected_std = np.sqrt(np.maximum(expected_var, 0.0))

    # Test all values
    np.testing.assert_allclose(mean.values, expected_mean, rtol=1e-12)
    np.testing.assert_allclose(std.values, expected_std, rtol=1e-12)

    # Zero denominator returns NaN or inf
    mean2, std2 = _std_from_moments(
        pd.Series([0.0, 5.0]), pd.Series([0.0, 30.0]), pd.Series([0.0, 0.0])
    )
    assert np.isnan(mean2.iloc[0])  # 0/0 = NaN
    assert np.isinf(mean2.iloc[1])  # 5/0 = inf
    assert np.isnan(std2.iloc[0]) and np.isnan(std2.iloc[1])  # std of inf/nan is nan


def test_compute_allele_frequencies_nan_and_has_cov():
    df = make_df(
        [
            {
                "contig": "c1",
                "position": 1,
                "total_coverage": 0.0,
                "A": 0.0,
                "C": 0.0,
                "G": 0.0,
                "T": 0.0,
            },
            {
                "contig": "c1",
                "position": 2,
                "total_coverage": 10.0,
                "A": 7.0,
                "C": 1.0,
                "G": 2.0,
                "T": 0.0,
            },
            {
                "contig": "c2",
                "position": 1,
                "total_coverage": 20.0,
                "A": 5.0,
                "C": 5.0,
                "G": 5.0,
                "T": 5.0,
            },
            {
                "contig": "c2",
                "position": 2,
                "total_coverage": 8.0,
                "A": 8.0,
                "C": 0.0,
                "G": 0.0,
                "T": 0.0,
            },
        ]
    )
    df = compute_allele_frequencies(df)

    # has_cov flags and squared coverage
    expected_has_cov = [False, True, True, True]
    expected_cov_sq = [0.0, 100.0, 400.0, 64.0]
    np.testing.assert_array_equal(df["has_cov"].values, expected_has_cov)
    np.testing.assert_allclose(
        df["total_coverage_sq"].values, expected_cov_sq, rtol=1e-12
    )

    # zero-coverage row has NaN freqs for all alleles
    zero_cov_mask = df["total_coverage"] == 0
    for allele in ALLELES:
        assert np.isnan(df.loc[zero_cov_mask, f"{allele}_freq"]).all()
        assert np.isnan(df.loc[zero_cov_mask, f"{allele}_freq_sq"]).all()

    # Check all frequency values for covered rows
    expected_freqs = {
        1: {"A": 0.7, "C": 0.1, "G": 0.2, "T": 0.0},  # row 1: 7/10, 1/10, 2/10, 0/10
        2: {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},  # row 2: 5/20 each
        3: {"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0},  # row 3: 8/8, 0/8, 0/8, 0/8
    }

    for i, (_, row) in enumerate(df[df["has_cov"]].iterrows()):
        expected = expected_freqs[i + 1]  # Use 1-based indexing for expected_freqs keys
        actual_freqs = [row[f"{allele}_freq"] for allele in ALLELES]
        expected_freqs_list = [expected[allele] for allele in ALLELES]
        np.testing.assert_allclose(actual_freqs, expected_freqs_list, rtol=1e-12)

        # Check squared frequencies
        expected_freq_sq = [freq**2 for freq in expected_freqs_list]
        actual_freq_sq = [row[f"{allele}_freq_sq"] for allele in ALLELES]
        np.testing.assert_allclose(actual_freq_sq, expected_freq_sq, rtol=1e-12)


def test_aggregate_with_freqs_counts_and_moments():
    # Multiple positions across contigs with different coverage patterns
    df = make_df(
        [
            {
                "contig": "c1",
                "position": 1,
                "total_coverage": 10.0,
                "A": 6.0,
                "C": 2.0,
                "G": 2.0,
                "T": 0.0,
            },
            {
                "contig": "c1",
                "position": 1,
                "total_coverage": 5.0,
                "A": 1.0,
                "C": 4.0,
                "G": 0.0,
                "T": 0.0,
            },
            {
                "contig": "c1",
                "position": 2,
                "total_coverage": 8.0,
                "A": 8.0,
                "C": 0.0,
                "G": 0.0,
                "T": 0.0,
            },
            {
                "contig": "c2",
                "position": 1,
                "total_coverage": 12.0,
                "A": 3.0,
                "C": 3.0,
                "G": 3.0,
                "T": 3.0,
            },
        ]
    )
    df = compute_allele_frequencies(df)
    agg = _aggregate_with_freqs(
        df, ["contig", "position"], present_col_name="n_present"
    )

    # Expected aggregated values for each position
    expected_results = {
        ("c1", 1): {
            "sum_cov": 15.0,
            "sumsq_cov": 125.0,  # 10^2 + 5^2 = 100 + 25
            "n_present": 2,
            "A_freq_sum": 6.0 / 10.0 + 1.0 / 5.0,  # 0.6 + 0.2 = 0.8
            "C_freq_sum": 2.0 / 10.0 + 4.0 / 5.0,  # 0.2 + 0.8 = 1.0
        },
        ("c1", 2): {
            "sum_cov": 8.0,
            "sumsq_cov": 64.0,
            "n_present": 1,
            "A_freq_sum": 1.0,  # 8/8
            "T_freq_sum": 0.0,
        },
        ("c2", 1): {
            "sum_cov": 12.0,
            "sumsq_cov": 144.0,
            "n_present": 1,
            "A_freq_sum": 0.25,  # 3/12
            "C_freq_sum": 0.25,
            "G_freq_sum": 0.25,
            "T_freq_sum": 0.25,
        },
    }

    # Test all aggregated values
    for pos_key, expected in expected_results.items():
        for col, expected_val in expected.items():
            actual_val = agg.loc[pos_key, col]
            np.testing.assert_allclose(
                actual_val,
                expected_val,
                rtol=1e-12,
                err_msg=f"Mismatch at {pos_key}[{col}]: got {actual_val}, expected {expected_val}",
            )

    # Verify all frequency moment columns exist and have correct sums
    for allele in ALLELES:
        freq_sum_col = f"{allele}_freq_sum"
        freq_sq_col = f"{allele}_freq_sumsq"
        assert freq_sum_col in agg.columns
        assert freq_sq_col in agg.columns

        # Check that frequency squared sums are computed correctly
        for pos_key in expected_results.keys():
            expected_freq_sum = df[
                df.set_index(["contig", "position"]).index == pos_key
            ][f"{allele}_freq"].sum()
            expected_freq_sq_sum = df[
                df.set_index(["contig", "position"]).index == pos_key
            ][f"{allele}_freq_sq"].sum()

            np.testing.assert_allclose(
                agg.loc[pos_key, freq_sum_col], expected_freq_sum, rtol=1e-12
            )
            np.testing.assert_allclose(
                agg.loc[pos_key, freq_sq_col], expected_freq_sq_sum, rtol=1e-12
            )


def test_overall_stats_both_present_and_with_zeros():
    # Multiple positions with different coverage patterns
    # Simulating 3 total samples, but positions have different presence patterns
    df = make_df(
        [
            # Position 1: present in 1 sample (out of 3 total)
            {
                "contig": "c1",
                "position": 1,
                "total_coverage": 10.0,
                "A": 7.0,
                "C": 3.0,
                "G": 0.0,
                "T": 0.0,
            },
            # Position 2: present in 2 samples (out of 3 total)
            {
                "contig": "c1",
                "position": 2,
                "total_coverage": 12.0,
                "A": 6.0,
                "C": 6.0,
                "G": 0.0,
                "T": 0.0,
            },
            {
                "contig": "c1",
                "position": 2,
                "total_coverage": 9.0,
                "A": 3.0,
                "C": 3.0,
                "G": 3.0,
                "T": 0.0,
            },
            # Position 3: present in all 3 samples
            {
                "contig": "c2",
                "position": 1,
                "total_coverage": 8.0,
                "A": 8.0,
                "C": 0.0,
                "G": 0.0,
                "T": 0.0,
            },
            {
                "contig": "c2",
                "position": 1,
                "total_coverage": 6.0,
                "A": 0.0,
                "C": 6.0,
                "G": 0.0,
                "T": 0.0,
            },
            {
                "contig": "c2",
                "position": 1,
                "total_coverage": 4.0,
                "A": 0.0,
                "C": 0.0,
                "G": 4.0,
                "T": 0.0,
            },
        ]
    )
    df = compute_allele_frequencies(df)

    # n_samples=3 total - now returns both present-only and with-zeros columns
    overall = compute_overall_statistics(df, n_samples=3)

    # Verify both column sets exist
    assert "mean_coverage" in overall.columns
    assert "std_coverage" in overall.columns
    assert "mean_coverage_with_zeros" in overall.columns
    assert "std_coverage_with_zeros" in overall.columns

    # Expected values for each position
    expected_coverage_means = {
        0: {
            "with_zeros": 10.0 / 3,
            "present_only": 10.0 / 1,
        },  # c1,1: 10 total, 1 present
        1: {
            "with_zeros": 21.0 / 3,
            "present_only": 21.0 / 2,
        },  # c1,2: 21 total, 2 present
        2: {
            "with_zeros": 18.0 / 3,
            "present_only": 18.0 / 3,
        },  # c2,1: 18 total, 3 present
    }

    expected_n_present = [1, 2, 3]
    expected_n_samples = [3, 3, 3]

    # Expected allele frequencies (always present-only denominators)
    expected_freq_means = {
        0: {"A": 7.0 / 10.0, "C": 3.0 / 10.0, "G": 0.0, "T": 0.0},  # c1,1
        1: {
            "A": (6.0 / 12.0 + 3.0 / 9.0) / 2,
            "C": (6.0 / 12.0 + 3.0 / 9.0) / 2,
            "G": 3.0 / 9.0 / 2,
            "T": 0.0,
        },  # c1,2
        2: {
            "A": (8.0 / 8.0) / 3,
            "C": (6.0 / 6.0) / 3,
            "G": (4.0 / 4.0) / 3,
            "T": 0.0,
        },  # c2,1
    }

    # Test both coverage means (present-only and with-zeros)
    for i in range(3):
        m_with_zeros = overall.loc[i, "mean_coverage_with_zeros"]
        m_present = overall.loc[i, "mean_coverage"]

        np.testing.assert_allclose(
            m_with_zeros, expected_coverage_means[i]["with_zeros"], rtol=1e-12
        )
        np.testing.assert_allclose(
            m_present, expected_coverage_means[i]["present_only"], rtol=1e-12
        )

        # Test n_present and n_samples
        assert overall.loc[i, "n_present"] == expected_n_present[i]
        assert overall.loc[i, "n_samples"] == expected_n_samples[i]

    # Test that allele frequencies use present-only denominators
    for i in range(3):
        for allele in ALLELES:
            actual_freq = overall.loc[i, f"mean_freq_{allele}"]
            expected_freq = expected_freq_means[i][allele]
            np.testing.assert_allclose(actual_freq, expected_freq, rtol=1e-12)


def test_grouped_stats_both_present_and_with_zeros():
    # Simpler scenario: 4 samples across 2 groups, 1 position, varying presence patterns
    rows = [
        # Position 1: G1 has 1/2 present, G2 has 2/2 present
        {
            "contig": "c1",
            "position": 1,
            "total_coverage": 10.0,
            "A": 8.0,
            "C": 2.0,
            "G": 0.0,
            "T": 0.0,
            "sample_id": "s1",
            "group": "G1",
        },
        # s2 from G1 is absent (no row)
        {
            "contig": "c1",
            "position": 1,
            "total_coverage": 12.0,
            "A": 6.0,
            "C": 6.0,
            "G": 0.0,
            "T": 0.0,
            "sample_id": "s3",
            "group": "G2",
        },
        {
            "contig": "c1",
            "position": 1,
            "total_coverage": 8.0,
            "A": 4.0,
            "C": 2.0,
            "G": 2.0,
            "T": 0.0,
            "sample_id": "s4",
            "group": "G2",
        },
    ]
    df = compute_allele_frequencies(make_df(rows))

    # Metadata lists all samples and groups
    meta = pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "s3", "s4"],
            "file_path": ["/f1", "/f2", "/f3", "/f4"],
            "group": ["G1", "G1", "G2", "G2"],
        }
    )

    # overall first
    result_base = compute_overall_statistics(df, n_samples=4)

    # grouped - now returns both present-only and with-zeros columns
    grouped = compute_grouped_statistics(
        df.copy(),
        meta.copy(),
        result_base.copy(),
        mag_id="test_mag",
    )

    # Verify both column types exist for each group
    for group in ["G1", "G2"]:
        assert f"mean_coverage_group_{group}" in grouped.columns
        assert f"mean_coverage_with_zeros_group_{group}" in grouped.columns
        assert f"std_coverage_group_{group}" in grouped.columns
        assert f"std_coverage_with_zeros_group_{group}" in grouped.columns
        assert f"n_present_group_{group}" in grouped.columns
        assert f"n_samples_group_{group}" in grouped.columns

    # Expected n_present and n_samples values per group
    expected_n_present = {"G1": 1, "G2": 2}  # Position 1 only
    expected_n_samples = {"G1": 2, "G2": 2}  # Total samples per group

    # Expected coverage means (with_zeros uses group size 2, present_only uses actual present count)
    expected_coverage_means = {
        "G1": {"with_zeros": 10.0 / 2, "present_only": 10.0 / 1},
        "G2": {"with_zeros": 20.0 / 2, "present_only": 20.0 / 2},  # (12+8)/2 = 20/2
    }

    # Expected allele frequencies (always present-only denominators)
    expected_allele_freqs = {
        "G1": {"A": 8.0 / 10.0, "C": 2.0 / 10.0, "G": 0.0, "T": 0.0},
        "G2": {
            "A": (6.0 / 12.0 + 4.0 / 8.0) / 2,
            "C": (6.0 / 12.0 + 2.0 / 8.0) / 2,
            "G": 2.0 / 8.0 / 2,
            "T": 0.0,
        },
    }

    # Test all n_present and n_samples values (single position index 0)
    for group in ["G1", "G2"]:
        n_pres_col = f"n_present_group_{group}"
        n_samp_col = f"n_samples_group_{group}"

        expected_n_pres = expected_n_present[group]
        expected_n_samp = expected_n_samples[group]
        actual_n_pres = grouped.loc[0, n_pres_col]
        actual_n_samp = grouped.loc[0, n_samp_col]

        # Check values and int typing
        assert isinstance(
            actual_n_pres, (int, np.integer)
        ), f"{n_pres_col} should be int, got {type(actual_n_pres)}"
        assert isinstance(
            actual_n_samp, (int, np.integer)
        ), f"{n_samp_col} should be int, got {type(actual_n_samp)}"
        assert actual_n_pres == expected_n_pres
        assert actual_n_samp == expected_n_samp

    # Test both coverage means (present-only and with-zeros) for all groups
    for group in ["G1", "G2"]:
        mean_present_col = f"mean_coverage_group_{group}"
        mean_with_zeros_col = f"mean_coverage_with_zeros_group_{group}"

        expected_with_zeros = expected_coverage_means[group]["with_zeros"]
        expected_present = expected_coverage_means[group]["present_only"]
        actual_with_zeros = grouped.loc[0, mean_with_zeros_col]
        actual_present = grouped.loc[0, mean_present_col]

        np.testing.assert_allclose(actual_with_zeros, expected_with_zeros, rtol=1e-12)
        np.testing.assert_allclose(actual_present, expected_present, rtol=1e-12)

    # Test that allele frequencies use present-only denominators for all groups
    for group in ["G1", "G2"]:
        for allele in ALLELES:
            col = f"mean_freq_{allele}_group_{group}"
            assert col in grouped.columns

            expected_freq = expected_allele_freqs[group][allele]
            actual = grouped.loc[0, col]
            np.testing.assert_allclose(actual, expected_freq, rtol=1e-12)


def test_comprehensive_workflow_with_zeros_and_grouping():
    """
    Comprehensive test validating the complete workflow with:
    - Both present-only and include-zeros statistics
    - Group and time dimensions
    - Proper handling of absent samples (zero coverage)
    """
    # Toy metadata with 4 samples across 2 groups and 2 timepoints
    meta = pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "s3", "s4"],
            "file_path": ["f1", "f2", "f3", "f4"],
            "group": ["A", "A", "B", "B"],
            "time": ["1", "2", "1", "2"],
        }
    )

    # Toy combined dataframe (as if after load_and_combine_sample_data)
    # Position 100: s1(A/1)=10, s2(A/2)=0, s3(B/1)=20 present
    # Position 101: s1(A/1)=0, s2(A/2)=0, s4(B/2)=30 present
    df = pd.DataFrame(
        {
            "contig": ["c1", "c1", "c1", "c1", "c1", "c1"],
            "position": [100, 100, 100, 101, 101, 101],
            "total_coverage": [10, 0, 20, 0, 0, 30],
            "A": [6, 0, 12, 0, 0, 20],
            "C": [2, 0, 5, 0, 0, 5],
            "G": [2, 0, 3, 0, 0, 5],
            "T": [0, 0, 0, 0, 0, 0],
            "sample_id": ["s1", "s2", "s3", "s1", "s2", "s4"],
            "group": ["A", "A", "B", "A", "A", "B"],
            "time": ["1", "2", "1", "1", "2", "2"],
        }
    )

    # Compute allele frequencies
    df = compute_allele_frequencies(df)

    # Test 1: Overall statistics with both present-only and include-zeros
    agg = _aggregate_with_freqs(
        df, ["contig", "position"], present_col_name="n_present"
    )

    # Present-only statistics (denominator = n_present)
    mean_p, std_p = _std_from_moments(
        agg["sum_cov"], agg["sumsq_cov"], agg["n_present"].astype(float)
    )

    # Include-zeros statistics (denominator = total samples)
    mean_z, std_z = _std_from_moments(
        agg["sum_cov"], agg["sumsq_cov"], float(len(meta))
    )

    # Position 100: s1=10, s3=20 present (s2=0, s4 absent)
    assert np.isclose(
        mean_p.loc[("c1", 100)], (10 + 20) / 2
    ), "Present-only mean should be 15"
    assert np.isclose(
        mean_z.loc[("c1", 100)], (10 + 0 + 20 + 0) / 4
    ), "Include-zeros mean should be 7.5"
    assert agg.loc[("c1", 100), "n_present"] == 2, "Should have 2 samples with coverage"

    # Position 101: only s4=30 present (s1=0, s2=0, s3 absent)
    assert np.isclose(mean_p.loc[("c1", 101)], 30 / 1), "Present-only mean should be 30"
    assert np.isclose(
        mean_z.loc[("c1", 101)], (0 + 0 + 30) / 4
    ), "Include-zeros mean should be 7.5"
    assert agg.loc[("c1", 101), "n_present"] == 1, "Should have 1 sample with coverage"

    # Test 2: Grouped statistics with proper group size mapping
    group_sizes = meta.groupby(["group", "time"], observed=True).size()
    idx = ["contig", "position", "group", "time"]
    grouped = _aggregate_with_freqs(df, idx, present_col_name="n_present_group")

    # Verify group sizes are correctly mapped
    keys_df = grouped.index.to_frame(index=False)[["group", "time"]]
    mi = pd.MultiIndex.from_frame(keys_df)
    denom = pd.Series(group_sizes.reindex(mi).values, index=grouped.index, dtype=float)

    # Group A/time 1 has only s1 in metadata
    assert denom.loc[("c1", 100, "A", "1")] == 1, "Group A/time 1 should have 1 sample"
    # Group B/time 1 has only s3 in metadata
    assert denom.loc[("c1", 100, "B", "1")] == 1, "Group B/time 1 should have 1 sample"
    # Group B/time 2 has only s4 in metadata
    assert denom.loc[("c1", 101, "B", "2")] == 1, "Group B/time 2 should have 1 sample"

    # Test 3: Verify overall statistics via compute_overall_statistics
    result_overall = compute_overall_statistics(df, n_samples=len(meta))

    # Check that both column types exist
    assert "mean_coverage" in result_overall.columns
    assert "mean_coverage_with_zeros" in result_overall.columns
    assert "std_coverage" in result_overall.columns
    assert "std_coverage_with_zeros" in result_overall.columns

    # Verify values for position 100
    pos100 = result_overall[result_overall["position"] == 100].iloc[0]
    assert np.isclose(pos100["mean_coverage"], 15.0), "Present-only mean should be 15"
    assert np.isclose(
        pos100["mean_coverage_with_zeros"], 7.5
    ), "Include-zeros mean should be 7.5"
    assert pos100["n_present"] == 2
    assert pos100["n_samples"] == 4

    # Test 4: Verify grouped statistics via compute_grouped_statistics
    result_grouped = compute_grouped_statistics(
        df.copy(), meta.copy(), result_overall.copy(), mag_id="test_mag"
    )

    # Check that grouped columns exist for both modes
    assert "mean_coverage_group_A" in result_grouped.columns
    assert "mean_coverage_with_zeros_group_A" in result_grouped.columns
    assert "mean_coverage_group_A_time_1" in result_grouped.columns
    assert "mean_coverage_with_zeros_group_A_time_1" in result_grouped.columns

    # Verify n_present and n_samples columns exist
    assert "n_present_group_A" in result_grouped.columns
    assert "n_samples_group_A" in result_grouped.columns
    assert "n_present_group_A_time_1" in result_grouped.columns
    assert "n_samples_group_A_time_1" in result_grouped.columns
