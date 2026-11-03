"""Tests for the two Poisson rate tests (statistics/poisson_rate_tests.py).

Every number below was computed by hand and checked against scipy /
statsmodels on 2026-11-01 (plan Task 4).  Both tests are likelihood-ratio
tests of Poisson counts with an exposure offset, scaled by an overdispersion
factor ``max(1, Pearson X2 / df_resid)`` and read off chi-square(1).

Test A, gene vs genome within a group (``gene_vs_genome``): observed events
``o_c`` per cage against expected ``e_c`` = the cage's genome-wide rate times
the gene's readable bases.  One rate ratio ``rr = sum(o) / sum(e)``;
LRT = 2 * (sum(o) * ln rr - (rr - 1) * sum(e)).

Test B, gene across groups (``group_vs_group``): per-cage counts and
readable bases in two groups; rates r_a, r_b, pooled r_0;
LRT = 2 * (O_a ln(r_a / r_0) + O_b ln(r_b / r_0)) with 0 ln 0 = 0.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from alleleflux.scripts.statistics.poisson_rate_tests import (
    RateTest,
    gene_vs_genome,
    group_vs_group,
)


class TestGeneVsGenome:
    def test_hand_case_three_cages(self):
        # o = [3, 0, 1], e = [0.2, 0.14, 0.5]; sum o 4, sum e 0.84, rr 4.7619
        # LRT = 2*[4 ln 4.7619 - 3.7619*0.84] = 2*[6.2426 - 3.1600] = 6.1652
        # mu = rr*e = [0.952, 0.667, 2.381]
        # X2 = (3-.952)^2/.952 + .667 + (1-2.381)^2/2.381 = 4.403 + 0.667 + 0.801 = 5.870
        # phi = 5.870 / 2 = 2.935; p = chi2.sf(6.1652 / 2.935, 1) = 0.1472
        r = gene_vs_genome(np.array([3, 0, 1]), np.array([0.2, 0.14, 0.5]))
        assert isinstance(r, RateTest)
        assert r.rate_ratio == pytest.approx(4.7619, rel=1e-3)
        assert r.deviance == pytest.approx(6.1652, rel=1e-3)
        assert r.dispersion == pytest.approx(2.935, rel=1e-3)
        assert r.p_value == pytest.approx(0.1472, abs=1e-3)
        assert r.n == 3 and r.note == ""

    def test_underdispersed_case_keeps_dispersion_at_one(self):
        # o = [1,1,1,1], e = [0.1]*4: rr 10, LRT = 2*(4 ln 10 - 9*0.4) = 11.2207,
        # X2 = 0 -> phi floored at 1, p = chi2.sf(11.2207, 1) = 8.09e-4
        r = gene_vs_genome(np.ones(4), np.full(4, 0.1))
        assert r.rate_ratio == pytest.approx(10.0)
        assert r.deviance == pytest.approx(11.2207, rel=1e-4)
        assert r.dispersion == 1.0
        assert r.p_value == pytest.approx(8.09e-4, rel=1e-2)

    def test_observed_equal_to_expected_is_null(self):
        r = gene_vs_genome(np.array([1.0, 2.0]), np.array([1.0, 2.0]))
        assert r.rate_ratio == pytest.approx(1.0)
        assert r.deviance == pytest.approx(0.0, abs=1e-12)
        assert r.p_value == pytest.approx(1.0)

    def test_no_events_is_untestable_not_p_one(self):
        r = gene_vs_genome(np.zeros(4), np.full(4, 0.3))
        assert r.rate_ratio == 0.0
        assert math.isnan(r.p_value) and math.isnan(r.deviance)
        assert r.note == "no_events"

    def test_no_exposure_is_untestable(self):
        r = gene_vs_genome(np.array([1.0, 0.0]), np.zeros(2))
        assert math.isnan(r.p_value) and math.isnan(r.rate_ratio)
        assert r.note == "no_exposure"

    def test_cages_with_zero_expected_are_dropped_first(self):
        # the third cage searched nothing in this gene: it cannot contribute
        full = gene_vs_genome(np.array([3, 0]), np.array([0.2, 0.14]))
        padded = gene_vs_genome(np.array([3, 0, 5]), np.array([0.2, 0.14, 0.0]))
        assert padded.n == 2
        assert padded.deviance == pytest.approx(full.deviance)
        assert padded.p_value == pytest.approx(full.p_value)

    def test_one_cage_has_no_dispersion_estimate_and_uses_one(self):
        r = gene_vs_genome(np.array([3]), np.array([0.5]))
        assert r.n == 1 and r.dispersion == 1.0
        assert r.deviance == pytest.approx(2 * (3 * math.log(6) - 5 * 0.5), rel=1e-6)

    def test_mismatched_lengths_raise(self):
        with pytest.raises(ValueError, match="length"):
            gene_vs_genome(np.array([1, 2]), np.array([0.5]))


class TestGroupVsGroup:
    def test_all_events_in_one_group(self):
        # o_a = [1,1,1,1], L_a = [100]*4; o_b = [0]*4, L_b = [100]*4
        # r_a = 0.01, r_b = 0, r_0 = 0.005; LRT = 2 * 4 ln 2 = 5.5452
        # X2 = 0 (every cage sits on its group's rate) -> phi 1; p = 0.01853
        r = group_vs_group(
            np.ones(4), np.full(4, 100.0), np.zeros(4), np.full(4, 100.0)
        )
        assert r.rate_ratio == math.inf
        assert r.deviance == pytest.approx(5.5452, rel=1e-4)
        assert r.dispersion == 1.0
        assert r.p_value == pytest.approx(0.01853, rel=1e-3)
        assert r.n == 8 and r.note == ""

    def test_all_events_in_the_other_group_gives_ratio_zero_same_p(self):
        r = group_vs_group(
            np.zeros(4), np.full(4, 100.0), np.ones(4), np.full(4, 100.0)
        )
        assert r.rate_ratio == 0.0
        assert r.p_value == pytest.approx(0.01853, rel=1e-3)

    def test_regular_case_matches_statsmodels_glm(self):
        # o_a=[2,0,1,3], L_a=[120,80,100,150]; o_b=[0,1,0,1], L_b=[100,100,90,110]
        # r_a = 6/450 = 0.013333, r_b = 2/400 = 0.005, LRT = 1.64959,
        # X2 = 3.568 -> 3.568/6 < 1 so phi = 1, p = 0.1990
        oa, la = np.array([2, 0, 1, 3.0]), np.array([120, 80, 100, 150.0])
        ob, lb = np.array([0, 1, 0, 1.0]), np.array([100, 100, 90, 110.0])
        r = group_vs_group(oa, la, ob, lb)
        assert r.rate_ratio == pytest.approx(0.013333 / 0.005, rel=1e-4)
        assert r.deviance == pytest.approx(1.64959, rel=1e-4)
        assert r.dispersion == 1.0
        assert r.p_value == pytest.approx(0.1990, rel=1e-3)
        sm = pytest.importorskip("statsmodels.api")
        o = np.concatenate([oa, ob])
        x = sm.add_constant(np.concatenate([np.ones(4), np.zeros(4)]))
        fit = sm.GLM(
            o, x, family=sm.families.Poisson(), offset=np.log(np.concatenate([la, lb]))
        ).fit()
        assert r.deviance == pytest.approx(fit.null_deviance - fit.deviance, abs=1e-6)

    def test_overdispersed_case_inflates_the_scale(self):
        # one burst cage in group a: counts [10, 0, 0, 0] over equal exposure
        r = group_vs_group(
            np.array([10, 0, 0, 0.0]), np.full(4, 100.0), np.zeros(4), np.full(4, 100.0)
        )
        # r_a = 0.025, mu_a = 2.5 each: X2 = (7.5^2 + 3*2.5^2)/2.5 = (56.25+18.75)/2.5 = 30; group b adds 0
        # phi = 30 / (8 - 2) = 5
        assert r.dispersion == pytest.approx(5.0)
        plain = 2 * 10 * math.log(0.025 / 0.0125)  # 13.863
        assert r.deviance == pytest.approx(plain, rel=1e-4)
        from scipy.stats import chi2

        assert r.p_value == pytest.approx(chi2.sf(plain / 5.0, 1), rel=1e-6)

    def test_no_events_anywhere_is_untestable(self):
        r = group_vs_group(np.zeros(3), np.full(3, 50.0), np.zeros(3), np.full(3, 50.0))
        assert math.isnan(r.p_value) and math.isnan(r.rate_ratio)
        assert r.note == "no_events"

    def test_a_group_with_no_exposure_is_untestable(self):
        r = group_vs_group(
            np.array([1.0]), np.array([10.0]), np.array([0.0]), np.array([0.0])
        )
        assert math.isnan(r.p_value)
        assert r.note == "no_exposure"

    def test_mismatched_lengths_raise(self):
        with pytest.raises(ValueError, match="length"):
            group_vs_group(
                np.array([1.0, 2.0]),
                np.array([10.0]),
                np.array([0.0]),
                np.array([10.0]),
            )
