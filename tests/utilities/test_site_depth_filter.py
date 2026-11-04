"""Tests for the MIDAS2 per-site depth filter (alleleflux.scripts.utilities.site_depth_filter).

The reference implementation inside this file (`_midas_reference`) is a literal
transcription of midas2 1.1.3 merge_snps.py lines 386-396 and run_snps.py lines
595-607 / 751, kept deliberately naive so the vectorised code is proven against
it by brute force on random tables.

House variations on MIDAS2 (spec revision 2026-10-15): the default floor is 0
(only the ratio ceiling removes positions), the mean coverage is passed IN (QC
computes it once; see Task 1b), and a mean of 0 makes the ratio test vacuous
rather than dropping everything.
"""

import argparse

import numpy as np
import pandas as pd
import pytest

from alleleflux.scripts.utilities.site_depth_filter import (
    MIDAS2_MIN_ALLELE_READS,
    MIDAS2_SITE_DEPTH,
    SiteDepthFilter,
    add_site_depth_filter_args,
    apply_site_depth_filter,
    mean_coverage,
    site_depth_filter_from_args,
    site_depth_mask,
)


def _midas_reference(rows, mean_cov, site_depth=5, site_ratio=3.0, min_allele_reads=2):
    """Naive MIDAS2 transcription: rows = list of (A, C, G, T)."""
    keep, low, high = [], 0, 0
    for r in rows:
        zeroed = [c if c > min_allele_reads else 0 for c in r]
        depth = sum(zeroed)
        if depth < site_depth:
            keep.append(False)
            low += 1
            continue
        if mean_cov > 0 and depth / mean_cov > site_ratio:
            keep.append(False)
            high += 1
            continue
        keep.append(True)
    return keep, low, high


def _frame(rows):
    return pd.DataFrame(rows, columns=["A", "C", "G", "T"])


MIDAS = SiteDepthFilter(
    site_depth=MIDAS2_SITE_DEPTH,
    min_allele_reads=MIDAS2_MIN_ALLELE_READS,
    keep_original_counts=False,
)  # MIDAS2's own rule, for the parity tests


class TestMeanCoverage:
    def test_ignores_zero_rows(self):
        # depths 10, 0, 20 -> covered sites are the 10 and the 20 -> 15, not 10
        assert (
            mean_coverage(_frame([(10, 0, 0, 0), (0, 0, 0, 0), (5, 5, 5, 5)])) == 15.0
        )

    def test_all_zero_profile_is_zero(self):
        assert mean_coverage(_frame([(0, 0, 0, 0), (0, 0, 0, 0)])) == 0.0

    def test_empty_frame_is_zero(self):
        assert mean_coverage(_frame([])) == 0.0


class TestDefaults:
    def test_default_floor_is_off(self):
        assert SiteDepthFilter() == SiteDepthFilter(True, 0, 3.0, 0, True)

    def test_default_keeps_thin_rows(self):
        # mean 10: a 1-read row and an all-2s row (depth' 0) both stay; MIDAS2
        # would drop both as "low".  Only the ceiling removes positions.
        df = _frame([(1, 0, 0, 0), (2, 2, 2, 2), (10, 0, 0, 0), (40, 0, 0, 0)])
        keep, s = site_depth_mask(df, SiteDepthFilter(), mean_cov=10.0)
        assert keep.tolist() == [True, True, True, False]
        assert s["n_removed_low"] == 0 and s["n_removed_high"] == 1 and s["n_kept"] == 3

    def test_zero_mean_removes_nothing_by_ratio(self):
        # no covered site means no pile to measure against: the ratio test is
        # vacuous, nothing is removed with the floor off
        keep, s = site_depth_mask(
            _frame([(0, 0, 0, 0), (3, 0, 0, 0)]), SiteDepthFilter(), mean_cov=0.0
        )
        assert keep.all()
        assert s == {
            "mean_coverage": 0.0,
            "mean_coverage_after_filter": 3.0,
            "n_sites_covered": 1,
            "n_removed_low": 0,
            "n_removed_high": 0,
            "n_kept": 2,
        }


class TestMidasParity:
    """With site_depth=5 the function is MIDAS2's merge_snps rule."""

    def test_counts_of_exactly_two_are_zeroed(self):
        # a base counts only when > 2; four 2s -> depth' 0 -> dropped as low,
        # even though the raw depth (8) clears the floor of 5
        df = _frame([(2, 2, 2, 2), (10, 0, 0, 0)])
        keep, s = site_depth_mask(df, MIDAS, mean_cov=mean_coverage(df))
        assert keep.tolist() == [False, True]
        assert s["n_removed_low"] == 1 and s["n_kept"] == 1

    def test_ratio_boundary_is_inclusive(self):
        # rows 3,3,3,3 and 18: mean (12 + 18) / 5 = 6; 18 / 6 = 3.0 exactly.
        # The 3-read rows fail the floor; the 18-read row sits AT 3.0 and MIDAS2
        # uses '>' not '>=', so it is kept.
        rows = [(3, 0, 0, 0)] * 4 + [(18, 0, 0, 0)]
        keep, s = site_depth_mask(_frame(rows), MIDAS, mean_cov=6.0)
        assert keep.tolist() == [False, False, False, False, True]
        assert s["n_removed_low"] == 4 and s["n_removed_high"] == 0

    def test_ratio_just_above_three_is_dropped(self):
        rows = [(3, 0, 0, 0)] * 4 + [(19, 0, 0, 0)]  # mean 6.2; 19 / 6.2 = 3.06
        keep, s = site_depth_mask(_frame(rows), MIDAS, mean_cov=6.2)
        assert keep.tolist()[-1] is False and s["n_removed_high"] == 1

    def test_low_wins_over_high(self):
        # a row failing both tests is counted once, as low (MIDAS2 continues at
        # the first failing test): mean 1 with min_allele_reads 0, floor 5,
        # row depth 4 -> 4 < 5 (low) and 4/1 > 3 (high) -> low
        keep, s = site_depth_mask(
            _frame([(4, 0, 0, 0)]),
            SiteDepthFilter(site_depth=5, min_allele_reads=0),
            mean_cov=1.0,
        )
        assert not keep.any() and (s["n_removed_low"], s["n_removed_high"]) == (1, 0)

    def test_all_zero_profile_drops_everything_under_a_floor(self):
        keep, s = site_depth_mask(
            _frame([(0, 0, 0, 0), (0, 0, 0, 0)]), MIDAS, mean_cov=0.0
        )
        assert not keep.any()
        assert s == {
            "mean_coverage": 0.0,
            "mean_coverage_after_filter": 0.0,
            "n_sites_covered": 0,
            "n_removed_low": 2,
            "n_removed_high": 0,
            "n_kept": 0,
        }


class TestMaskMechanics:
    def test_disabled_keeps_every_row_but_still_summarises(self):
        df = _frame([(2, 2, 2, 2), (100, 0, 0, 0), (10, 0, 0, 0)])
        keep, s = site_depth_mask(df, SiteDepthFilter(enabled=False), mean_cov=10.0)
        assert keep.all() and s["n_kept"] == 3 and s["n_sites_covered"] == 3
        assert s["mean_coverage"] == 10.0 and s["n_removed_high"] == 0

    def test_mask_is_aligned_to_a_non_default_index(self):
        df = _frame([(10, 0, 0, 0), (40, 0, 0, 0)])
        df.index = [7, 3]
        keep, _ = site_depth_mask(df, SiteDepthFilter(), mean_cov=10.0)
        assert df[keep].index.tolist() == [7]

    def test_mask_does_not_modify_the_input_frame(self):
        df = _frame([(2, 2, 2, 2), (10, 0, 0, 0)])
        before = df.copy()
        site_depth_mask(df, SiteDepthFilter(min_allele_reads=2), mean_cov=10.0)
        pd.testing.assert_frame_equal(df, before)


class TestMeanCoverageAfterFilter:
    """The summary also reports the covered-sites mean of the rows that continue."""

    def test_default_drops_the_pile_and_the_mean_falls(self):
        # rows 10 x 4, 40, 80 against QC mean 24: the 80 goes, the rest stay
        df = _frame([(10, 0, 0, 0)] * 4 + [(40, 0, 0, 0), (80, 0, 0, 0)])
        _, s = site_depth_mask(df, SiteDepthFilter(), mean_cov=24.0)
        assert s["mean_coverage"] == 24.0  # the reference, unchanged
        assert s["mean_coverage_after_filter"] == pytest.approx(80 / 5)  # (40 + 40) / 5

    def test_midas_zeroing_uses_the_zeroed_depths(self):
        # MIDAS2 settings: (18, 2, 0, 0) continues as depth 18, not 20
        df = _frame([(18, 2, 0, 0), (10, 0, 0, 0)])
        _, s = site_depth_mask(df, MIDAS, mean_cov=16.0)
        assert s["mean_coverage_after_filter"] == pytest.approx(14.0)

    def test_keep_original_counts_uses_the_original_depths(self):
        df = _frame([(18, 2, 0, 0), (10, 0, 0, 0)])
        flt = SiteDepthFilter(
            site_depth=5, min_allele_reads=2
        )  # keep_original_counts default
        _, s = site_depth_mask(df, flt, mean_cov=16.0)
        assert s["mean_coverage_after_filter"] == pytest.approx(15.0)

    def test_disabled_reports_the_unfiltered_mean(self):
        df = _frame([(10, 0, 0, 0), (0, 0, 0, 0), (30, 0, 0, 0)])
        _, s = site_depth_mask(df, SiteDepthFilter(enabled=False), mean_cov=99.0)
        assert s["mean_coverage_after_filter"] == 20.0


class TestApplySiteDepthFilter:
    """MIDAS2 semantics: the zeroed counts are what continues downstream."""

    def test_default_leaves_kept_counts_untouched(self):
        # min_allele_reads 0: nothing is zeroed, kept rows are byte-identical
        df = _frame([(18, 2, 0, 0), (60, 0, 0, 0)]).assign(total_coverage=[20, 60])
        out, s = apply_site_depth_filter(df, SiteDepthFilter(), mean_cov=16.0)
        assert out[["A", "C", "G", "T", "total_coverage"]].values.tolist() == [
            [18, 2, 0, 0, 20]
        ]
        assert s["n_removed_high"] == 1 and s["n_kept"] == 1

    def test_midas_zeroing_is_carried_downstream(self):
        # MIDAS2 merge_snps: a kept (18, 2, 0, 0) continues as (18, 0, 0, 0)
        # with depth 18; the 2-read minor allele disappears
        df = _frame([(18, 2, 0, 0), (2, 2, 2, 2)]).assign(total_coverage=[20, 8])
        out, s = apply_site_depth_filter(df, MIDAS, mean_cov=16.0)
        assert out[["A", "C", "G", "T", "total_coverage"]].values.tolist() == [
            [18, 0, 0, 0, 18]
        ]
        assert s["n_removed_low"] == 1

    def test_without_total_coverage_column(self):
        # MIDAS2 zeroing on a frame that has no total_coverage column (ANI loader)
        df = _frame([(18, 2, 0, 0)])
        out, _ = apply_site_depth_filter(
            df,
            SiteDepthFilter(min_allele_reads=2, keep_original_counts=False),
            mean_cov=16.0,
        )
        assert out.values.tolist() == [[18, 0, 0, 0]]
        assert "total_coverage" not in out.columns

    def test_other_columns_and_dtypes_survive(self):
        df = pd.DataFrame(
            {
                "contig": ["c1", "c1"],
                "position": [3, 9],
                "A": [18, 60],
                "C": [2, 0],
                "G": [0, 0],
                "T": [0, 0],
                "gene_id": ["g", "g"],
            }
        )
        df[["A", "C", "G", "T"]] = df[["A", "C", "G", "T"]].astype("int32")
        out, _ = apply_site_depth_filter(
            df, SiteDepthFilter(min_allele_reads=2), mean_cov=16.0
        )
        assert out["position"].tolist() == [3] and out["gene_id"].tolist() == ["g"]
        assert str(out["A"].dtype) == "int32"
        assert out.index.tolist() == [0]

    def test_keep_original_counts_default_skips_the_zeroing(self):
        # house default: zeroing decides keep/drop on a scratch copy, the kept
        # row continues with its original counts
        df = _frame([(18, 2, 0, 0), (2, 2, 2, 2)]).assign(total_coverage=[20, 8])
        flt = SiteDepthFilter(site_depth=5, min_allele_reads=2)
        assert flt.keep_original_counts is True
        out, s = apply_site_depth_filter(df, flt, mean_cov=16.0)
        assert out[["A", "C", "G", "T", "total_coverage"]].values.tolist() == [
            [18, 2, 0, 0, 20]
        ]
        assert s["n_removed_low"] == 1  # the (2,2,2,2) row is still dropped

    def test_disabled_returns_the_frame_as_is(self):
        df = _frame([(18, 2, 0, 0), (60, 0, 0, 0)])
        out, s = apply_site_depth_filter(
            df, SiteDepthFilter(enabled=False, min_allele_reads=2), mean_cov=16.0
        )
        pd.testing.assert_frame_equal(out, df)
        assert s["n_kept"] == 2

    @pytest.mark.parametrize("bad", [float("nan"), float("inf"), -1.0])
    def test_rejects_a_non_finite_or_negative_mean(self, bad):
        # the mean comes from a QC table on trust; NaN would silently keep
        # every row and a negative number is a corrupt table, so fail loud
        with pytest.raises(ValueError, match="mean coverage"):
            site_depth_mask(_frame([(10, 0, 0, 0)]), SiteDepthFilter(), mean_cov=bad)

    def test_empty_frame(self):
        keep, s = site_depth_mask(_frame([]), SiteDepthFilter(), mean_cov=0.0)
        assert len(keep) == 0 and s["n_kept"] == 0

    @pytest.mark.parametrize("seed", range(40))
    def test_brute_force_against_midas_transcription(self, seed):
        rng = np.random.default_rng(seed)
        n = int(rng.integers(0, 60))
        rows = [tuple(int(v) for v in rng.integers(0, 12, size=4)) for _ in range(n)]
        # sprinkle a few piles and a few all-zero rows
        rows += [
            (int(rng.integers(40, 200)), 0, 0, 0)
            for _ in range(int(rng.integers(0, 3)))
        ]
        rows += [(0, 0, 0, 0)] * int(rng.integers(0, 3))
        flt = SiteDepthFilter(
            site_depth=int(rng.integers(0, 8)),
            site_ratio=float(rng.choice([1.5, 2.0, 3.0, 5.0])),
            min_allele_reads=int(rng.integers(0, 4)),
        )
        df = _frame(rows)
        # the mean is QC's job; here it is the MIDAS2 mean of this very table,
        # or sometimes an arbitrary external value, since the function must
        # simply use what it is given
        mean_cov = (
            mean_coverage(df) if rng.random() < 0.7 else float(rng.integers(0, 30))
        )
        keep, s = site_depth_mask(df, flt, mean_cov=mean_cov)
        ref_keep, ref_low, ref_high = _midas_reference(
            rows, mean_cov, flt.site_depth, flt.site_ratio, flt.min_allele_reads
        )
        assert keep.tolist() == ref_keep
        assert s["mean_coverage"] == mean_cov
        assert s["n_sites_covered"] == sum(1 for r in rows if sum(r) > 0)
        assert (s["n_removed_low"], s["n_removed_high"], s["n_kept"]) == (
            ref_low,
            ref_high,
            sum(ref_keep),
        )


class TestCliHelpers:
    def test_defaults_are_the_house_defaults(self):
        p = argparse.ArgumentParser()
        add_site_depth_filter_args(p)
        flt = site_depth_filter_from_args(p.parse_args([]))
        assert flt == SiteDepthFilter(True, 0, 3.0, 0, True)

    def test_off_and_overrides(self):
        p = argparse.ArgumentParser()
        add_site_depth_filter_args(p)
        flt = site_depth_filter_from_args(
            p.parse_args(
                [
                    "--disable_site_depth_filter",
                    "--site_depth",
                    "7",
                    "--site_ratio",
                    "2.5",
                    "--min_allele_reads",
                    "0",
                ]
            )
        )
        assert flt == SiteDepthFilter(False, 7, 2.5, 0, True)

    def test_zero_low_counts_flag(self):
        p = argparse.ArgumentParser()
        add_site_depth_filter_args(p)
        flt = site_depth_filter_from_args(p.parse_args(["--zero_low_counts"]))
        assert flt.keep_original_counts is False
        assert (
            site_depth_filter_from_args(p.parse_args([])).keep_original_counts is True
        )
