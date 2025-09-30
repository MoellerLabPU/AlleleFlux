import numpy as np
import pandas as pd
import pytest

from alleleflux.scripts.accessory.mean_coverage_per_position import (
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


def test_overall_stats_freq_denominator_present_only():
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

    # n_samples=3 total
    overall_inc_zeros = compute_overall_statistics(
        df, n_samples=3, treat_absent_as_zero=True
    )
    overall_present_only = compute_overall_statistics(
        df, n_samples=3, treat_absent_as_zero=False
    )

    # Expected values for each position
    expected_coverage_means = {
        0: {
            "inc_zeros": 10.0 / 3,
            "present_only": 10.0 / 1,
        },  # c1,1: 10 total, 1 present
        1: {
            "inc_zeros": 21.0 / 3,
            "present_only": 21.0 / 2,
        },  # c1,2: 21 total, 2 present
        2: {
            "inc_zeros": 18.0 / 3,
            "present_only": 18.0 / 3,
        },  # c2,1: 18 total, 3 present
    }

    expected_n_present = [1, 2, 3]
    expected_n_samples = [3, 3, 3]

    # Expected allele frequencies (same for both modes - present-only denominators)
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

    # Test all coverage means
    for i in range(3):
        m_inc = overall_inc_zeros.loc[i, "mean_coverage"]
        m_pres = overall_present_only.loc[i, "mean_coverage"]

        np.testing.assert_allclose(
            m_inc, expected_coverage_means[i]["inc_zeros"], rtol=1e-12
        )
        np.testing.assert_allclose(
            m_pres, expected_coverage_means[i]["present_only"], rtol=1e-12
        )

        # Test n_present and n_samples
        assert overall_inc_zeros.loc[i, "n_present"] == expected_n_present[i]
        assert overall_present_only.loc[i, "n_present"] == expected_n_present[i]
        assert overall_inc_zeros.loc[i, "n_samples"] == expected_n_samples[i]
        assert overall_present_only.loc[i, "n_samples"] == expected_n_samples[i]

    # Test that allele frequencies are identical across modes (present-only denominators)
    for i in range(3):
        for allele in ALLELES:
            c1 = overall_inc_zeros.loc[i, f"mean_freq_{allele}"]
            c2 = overall_present_only.loc[i, f"mean_freq_{allele}"]
            expected_freq = expected_freq_means[i][allele]

            # Both should equal the expected frequency
            np.testing.assert_allclose(c1, expected_freq, rtol=1e-12)
            np.testing.assert_allclose(c2, expected_freq, rtol=1e-12)

            # And they should be identical to each other
            np.testing.assert_allclose(c1, c2, rtol=1e-12)


def test_grouped_stats_freq_denominator_present_only_and_counts():
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
    result_base = compute_overall_statistics(df, n_samples=4, treat_absent_as_zero=True)

    # grouped with include_zeros True/False
    grouped_inc = compute_grouped_statistics(
        df.copy(), meta.copy(), result_base.copy(), treat_absent_as_zero=True
    )
    grouped_pres = compute_grouped_statistics(
        df.copy(), meta.copy(), result_base.copy(), treat_absent_as_zero=False
    )

    # Expected n_present values per group
    expected_n_present = {"G1": 1, "G2": 2}  # Position 1 only

    # Expected coverage means (include_zeros uses group size 2, present_only uses actual present count)
    expected_coverage_means = {
        "G1": {"inc": 10.0 / 2, "pres": 10.0 / 1},
        "G2": {"inc": 20.0 / 2, "pres": 20.0 / 2},  # (12+8)/2 = 20/2
    }

    # Expected allele frequencies (always present-only denominators, same across modes)
    expected_allele_freqs = {
        "G1": {"A": 8.0 / 10.0, "C": 2.0 / 10.0, "G": 0.0, "T": 0.0},
        "G2": {
            "A": (6.0 / 12.0 + 4.0 / 8.0) / 2,
            "C": (6.0 / 12.0 + 2.0 / 8.0) / 2,
            "G": 2.0 / 8.0 / 2,
            "T": 0.0,
        },
    }

    # Test all n_present values and their typing (single position index 0)
    for group in ["G1", "G2"]:
        col = f"n_present_group_{group}"
        assert col in grouped_inc.columns and col in grouped_pres.columns

        expected_n_pres = expected_n_present[group]
        actual_inc = grouped_inc.loc[0, col]
        actual_pres = grouped_pres.loc[0, col]

        # Check values and int typing
        assert isinstance(
            actual_inc, (int, np.integer)
        ), f"{col} should be int, got {type(actual_inc)}"
        assert isinstance(
            actual_pres, (int, np.integer)
        ), f"{col} should be int, got {type(actual_pres)}"
        assert actual_inc == expected_n_pres
        assert actual_pres == expected_n_pres

    # Test coverage means for all groups (single position index 0)
    for group in ["G1", "G2"]:
        mean_col = f"mean_coverage_group_{group}"
        assert mean_col in grouped_inc.columns and mean_col in grouped_pres.columns

        expected_inc = expected_coverage_means[group]["inc"]
        expected_pres = expected_coverage_means[group]["pres"]
        actual_inc = grouped_inc.loc[0, mean_col]
        actual_pres = grouped_pres.loc[0, mean_col]

        np.testing.assert_allclose(actual_inc, expected_inc, rtol=1e-12)
        np.testing.assert_allclose(actual_pres, expected_pres, rtol=1e-12)

    # Test that allele frequencies are identical across modes for all groups (single position index 0)
    for group in ["G1", "G2"]:
        for allele in ALLELES:
            col = f"mean_freq_{allele}_group_{group}"
            assert col in grouped_inc.columns and col in grouped_pres.columns

            expected_freq = expected_allele_freqs[group][allele]
            actual_inc = grouped_inc.loc[0, col]
            actual_pres = grouped_pres.loc[0, col]

            # Both modes should give same result (present-only denominators)
            np.testing.assert_allclose(actual_inc, expected_freq, rtol=1e-12)
            np.testing.assert_allclose(actual_pres, expected_freq, rtol=1e-12)
            np.testing.assert_allclose(actual_inc, actual_pres, rtol=1e-12)
