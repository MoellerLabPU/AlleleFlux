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
LRT = 2 * (O_a ln(r_a / r_0) + O_b ln(r_b / r_0)) with 0 ln 0 = 0.  Its
scaled statistic is referred to F(1, df_resid), not chi-square(1) (decision
2026-11-03 after the calibration simulation: chi-square raised false alarms
at 12 % with heterogeneous cages, F holds 5 %); Test A stays on chi-square,
where F halved the power for no calibration gain.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from alleleflux.scripts.statistics.poisson_rate_tests import (
    RateTest,
    gene_vs_genome,
    group_vs_group,
    result_columns,
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
        # the same surprise without the agreement discount: chi2.sf(6.1652, 1)
        assert r.p_value_raw == pytest.approx(0.0130, abs=1e-3)
        assert r.n == 3 and r.note == ""
        # Both curves for both statistics (collaborator request 2026-11-04):
        # Test A's house pick is chi-square, so p_value_chi2 IS p_value; the F
        # twin uses df_resid = cages - 1 = 2.
        from scipy.stats import f as f_dist

        assert r.df_resid == 2
        assert r.p_value_chi2 == r.p_value
        assert r.p_value_f == pytest.approx(f_dist.sf(6.1652 / 2.935, 1, 2), rel=1e-3)
        assert r.p_value_raw_chi2 == r.p_value_raw
        assert r.p_value_raw_f == pytest.approx(f_dist.sf(6.1652, 1, 2), rel=1e-3)
        assert r.p_value_f > r.p_value and r.p_value_raw_f > r.p_value_raw

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
        assert math.isnan(r.p_value_raw)
        assert r.note == "no_events"

    def test_no_exposure_is_untestable(self):
        # no cage read any of the gene, so no cage can have hit it either
        r = gene_vs_genome(np.zeros(2), np.zeros(2))
        assert math.isnan(r.p_value) and math.isnan(r.rate_ratio)
        assert r.note == "no_exposure"

    def test_cages_with_zero_expected_and_zero_observed_are_dropped_first(self):
        # the third cage searched nothing in this gene: it cannot contribute
        full = gene_vs_genome(np.array([3, 0]), np.array([0.2, 0.14]))
        padded = gene_vs_genome(np.array([3, 0, 0]), np.array([0.2, 0.14, 0.0]))
        assert padded.n == 2
        assert padded.deviance == pytest.approx(full.deviance)
        assert padded.p_value == pytest.approx(full.p_value)

    def test_an_event_where_nothing_was_expected_is_a_contract_violation(self):
        # expected is 0 only when the cage read none of the gene or had no
        # events anywhere; an observed event there means the tables disagree
        with pytest.raises(ValueError, match="expected 0"):
            gene_vs_genome(np.array([3, 0, 5]), np.array([0.2, 0.14, 0.0]))

    def test_one_cage_has_no_dispersion_estimate_and_uses_one(self):
        r = gene_vs_genome(np.array([3]), np.array([0.5]))
        assert r.n == 1 and r.dispersion == 1.0
        assert r.deviance == pytest.approx(2 * (3 * math.log(6) - 5 * 0.5), rel=1e-6)
        # no spare cage: the chi-square twins exist, the F twins cannot
        assert r.df_resid == 0
        assert not math.isnan(r.p_value_chi2) and not math.isnan(r.p_value_raw_chi2)
        assert math.isnan(r.p_value_f) and math.isnan(r.p_value_raw_f)

    def test_mismatched_lengths_raise(self):
        with pytest.raises(ValueError, match="length"):
            gene_vs_genome(np.array([1, 2]), np.array([0.5]))


class TestGroupVsGroup:
    def test_all_events_in_one_group(self):
        # o_a = [1,1,1,1], L_a = [100]*4; o_b = [0]*4, L_b = [100]*4
        # r_a = 0.01, r_b = 0, r_0 = 0.005; LRT = 2 * 4 ln 2 = 5.5452
        # X2 = 0 (every cage sits on its group's rate) -> phi 1
        # p = F.sf(5.5452, 1, 8 - 2) = 0.0567 (chi-square would say 0.0185)
        r = group_vs_group(
            np.ones(4), np.full(4, 100.0), np.zeros(4), np.full(4, 100.0)
        )
        assert r.rate_ratio == math.inf
        assert r.deviance == pytest.approx(5.5452, rel=1e-4)
        assert r.dispersion == 1.0
        assert r.p_value == pytest.approx(0.0567, rel=1e-2)
        assert r.p_value_raw == pytest.approx(0.0185, rel=1e-2)  # chi2, no discount
        assert r.n == 8 and r.note == ""
        assert r.df_resid == 6
        # Both curves: Test B's house pick is F, so p_value_f IS p_value; with
        # dispersion 1 the raw and discounted statistics coincide, so the raw F
        # twin equals the discounted p and the chi-square twin equals the raw p.
        assert r.p_value_f == r.p_value
        assert r.p_value_chi2 == pytest.approx(0.0185, rel=1e-2)
        assert r.p_value_raw_chi2 == r.p_value_raw
        assert r.p_value_raw_f == pytest.approx(0.0567, rel=1e-2)

    def test_eight_cages_per_group_one_event_each_is_significant(self):
        # LRT = 2 * 8 ln 2 = 11.09; F.sf(11.09, 1, 14) = 0.0050
        r = group_vs_group(
            np.ones(8), np.full(8, 100.0), np.zeros(8), np.full(8, 100.0)
        )
        assert r.p_value == pytest.approx(0.0050, rel=2e-2)
        assert r.df_resid == 14

    def test_all_events_in_the_other_group_gives_ratio_zero_same_p(self):
        r = group_vs_group(
            np.zeros(4), np.full(4, 100.0), np.ones(4), np.full(4, 100.0)
        )
        assert r.rate_ratio == 0.0
        assert r.p_value == pytest.approx(0.0567, rel=1e-2)

    def test_regular_case_matches_statsmodels_glm(self):
        # o_a=[2,0,1,3], L_a=[120,80,100,150]; o_b=[0,1,0,1], L_b=[100,100,90,110]
        # r_a = 6/450 = 0.013333, r_b = 2/400 = 0.005, LRT = 1.64959,
        # X2 = 3.568 -> 3.568/6 < 1 so phi = 1, p = F.sf(1.64959, 1, 6) = 0.2464
        oa, la = np.array([2, 0, 1, 3.0]), np.array([120, 80, 100, 150.0])
        ob, lb = np.array([0, 1, 0, 1.0]), np.array([100, 100, 90, 110.0])
        r = group_vs_group(oa, la, ob, lb)
        assert r.rate_ratio == pytest.approx(0.013333 / 0.005, rel=1e-4)
        assert r.deviance == pytest.approx(1.64959, rel=1e-4)
        assert r.dispersion == 1.0
        assert r.p_value == pytest.approx(0.2464, rel=1e-3)
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
        from scipy.stats import f as fdist

        assert r.p_value == pytest.approx(fdist.sf(plain / 5.0, 1, 6), rel=1e-6)
        from scipy.stats import chi2

        # raw: the burst looks overwhelming (0.0002) once cage agreement is ignored
        assert r.p_value_raw == pytest.approx(chi2.sf(plain, 1), rel=1e-6)
        assert r.p_value_raw < 1e-3 < r.p_value

    def test_one_cage_per_group_is_untestable_there_is_no_replicate(self):
        # 2 cages - 2 rates = 0 spare cages: the dispersion cannot be estimated
        # and a cage-as-replicate test has no replicate, so no p-value at all
        r = group_vs_group(
            np.array([3.0]), np.array([100.0]), np.array([0.0]), np.array([100.0])
        )
        assert r.df_resid == 0
        assert math.isnan(r.p_value) and math.isnan(r.dispersion)
        assert r.note == "no_residual_df"
        assert r.rate_ratio == math.inf  # the effect size is still reported
        # the raw p needs no dispersion estimate, so it IS reported: chi2 on 2*3 ln 2
        from scipy.stats import chi2

        assert r.p_value_raw == pytest.approx(chi2.sf(2 * 3 * math.log(2), 1), rel=1e-6)
        assert r.p_value_raw_chi2 == r.p_value_raw
        # no F curve without a spare cage, discounted or raw
        assert math.isnan(r.p_value_chi2) and math.isnan(r.p_value_f)
        assert math.isnan(r.p_value_raw_f)

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


# ---------------------------------------------------------------------------
# Task 5: the per-gene driver (floors + BH)
# ---------------------------------------------------------------------------

import pandas as pd
from statsmodels.stats.multitest import multipletests

from alleleflux.scripts.statistics.poisson_rate_tests import (
    test_every_gene as every_gene,  # aliased: pytest would collect the bare name
)

GROUPS = ("treatment", "control")
CAGES = [("treatment", c) for c in "ABCD"] + [("control", c) for c in "EFGH"]


def _gene_cage(rows, mag="M"):
    """rows: (gene_id, group, replicate, n_de_novo, covered_bp, gene_length) -> Task 3's gene x cage frame."""
    df = pd.DataFrame(
        rows,
        columns=[
            "gene_id",
            "group",
            "replicate",
            "n_de_novo",
            "covered_bp",
            "gene_length",
        ],
    )
    df["mag_id"] = mag
    df["entered"] = df["covered_bp"] > 0
    df["covered_fraction"] = df["covered_bp"] / df["gene_length"]
    return df


def _cage_totals(rows, mag="M"):
    """rows: (group, replicate, n_de_novo_total, covered_bp_total)."""
    df = pd.DataFrame(
        rows, columns=["group", "replicate", "n_de_novo_total", "covered_bp_total"]
    )
    df["mag_id"] = mag
    return df


def _plant_fixture(hot_events=(1, 1, 1, 1), gene_bp=100, genome_bp=1000):
    """Three genes x eight cages.  gX: one event per treatment cage (the
    regression plant's shape), none in control; gY: nothing anywhere;
    gZ: one event in one control cage.  Every gene is fully readable in every
    cage; each cage's genome-wide total is exactly its events in these genes."""
    gene_rows, totals = [], []
    for group, rep in CAGES:
        x = hot_events["ABCD".index(rep)] if group == "treatment" else 0
        z = 1 if (group, rep) == ("control", "E") else 0
        gene_rows += [
            ("gX", group, rep, x, gene_bp, gene_bp),
            ("gY", group, rep, 0, gene_bp, gene_bp),
            ("gZ", group, rep, z, gene_bp, gene_bp),
        ]
        totals.append((group, rep, x + z, genome_bp))
    return _gene_cage(gene_rows), _cage_totals(totals)


class TestTestEveryGene:
    def test_hot_gene_wins_both_tests_and_empty_gene_is_untested(self):
        gc, tot = _plant_fixture()
        out = every_gene(gc, tot, GROUPS, min_cages=4, fdr_by_mag=False)
        assert list(out.gene_id) == ["gX", "gY", "gZ"]
        x = out.set_index("gene_id").loc["gX"]
        # Test A, treatment: each treatment cage has 1 event over 1,000 bases,
        # so gX (100 bases) is owed 0.1 per cage: o = [1,1,1,1], e = [0.1]*4
        assert x["n_cages_treatment"] == 4
        assert x["n_de_novo_treatment"] == 4
        assert x["expected_treatment"] == pytest.approx(0.4)
        assert x["rate_ratio_vs_genome_treatment"] == pytest.approx(10.0)
        assert x["p_value_vs_genome_treatment"] == pytest.approx(
            gene_vs_genome(np.ones(4), np.full(4, 0.1)).p_value
        )
        assert x["note_treatment"] == ""
        # Test A, control: gX has no event in control -> untestable, not p = 1
        assert x["n_de_novo_control"] == 0
        assert math.isnan(x["p_value_vs_genome_control"])
        assert x["note_control"] == "no_events"
        # Test B: 4 cages x 1 event over 100 bases vs 4 cages x 0 over 100
        assert x["rate_treatment"] == pytest.approx(4 / 400)
        assert x["rate_control"] == 0.0
        assert x["rate_ratio_treatment_vs_control"] == math.inf
        assert x["p_value_between"] == pytest.approx(0.0567, rel=1e-2)  # F(1, 6)
        assert x["p_value_raw_between"] == pytest.approx(0.0185, rel=1e-2)  # chi2
        # two genes are tested between the groups (gX and gZ), so BH doubles
        # the smaller p: q = min(0.0185 * 2 / 1, gZ's larger raw p) = 0.037
        assert x["q_value_raw_between"] == pytest.approx(2 * x["p_value_raw_between"])
        assert x["p_value_raw_vs_genome_treatment"] == pytest.approx(
            gene_vs_genome(np.ones(4), np.full(4, 0.1)).p_value_raw
        )
        assert x["df_resid_between"] == 6
        assert x["note_between"] == ""
        # The statistics and both curves, per block (collaborator request).
        from scipy.stats import f as f_dist

        ra = gene_vs_genome(np.ones(4), np.full(4, 0.1))
        assert x["lrt_vs_genome_treatment"] == pytest.approx(ra.deviance)
        assert x["df_resid_vs_genome_treatment"] == 3
        assert x["p_value_chi2_vs_genome_treatment"] == x["p_value_vs_genome_treatment"]
        assert x["p_value_f_vs_genome_treatment"] == pytest.approx(
            f_dist.sf(ra.deviance / ra.dispersion, 1, 3)
        )
        assert (
            x["p_value_raw_chi2_vs_genome_treatment"]
            == x["p_value_raw_vs_genome_treatment"]
        )
        assert x["p_value_raw_f_vs_genome_treatment"] == pytest.approx(
            f_dist.sf(ra.deviance, 1, 3)
        )
        assert x["lrt_between"] == pytest.approx(5.5452, rel=1e-4)
        assert x["p_value_f_between"] == x["p_value_between"]
        assert x["p_value_chi2_between"] == pytest.approx(0.0185, rel=1e-2)
        assert x["p_value_raw_chi2_between"] == x["p_value_raw_between"]
        assert x["p_value_raw_f_between"] == pytest.approx(0.0567, rel=1e-2)
        # Every twin is BH-corrected as its own family (user 2026-11-04): the
        # headline twin's q equals the headline q, the other twin gets its own.
        # Two genes are tested between the groups, so the smaller p doubles.
        assert x["q_value_f_between"] == x["q_value_between"]
        assert x["q_value_chi2_between"] == pytest.approx(2 * x["p_value_chi2_between"])
        assert x["q_value_raw_chi2_between"] == x["q_value_raw_between"]
        assert x["q_value_raw_f_between"] == pytest.approx(
            2 * x["p_value_raw_f_between"]
        )
        assert x["q_value_chi2_vs_genome_treatment"] == x["q_value_vs_genome_treatment"]
        assert x["q_value_f_vs_genome_treatment"] == pytest.approx(
            x["p_value_f_vs_genome_treatment"]  # gX is the only treatment gene tested
        )
        # gY: no events anywhere -> every test untestable, every q NaN
        y = out.set_index("gene_id").loc["gY"]
        for col in (
            "p_value_vs_genome_treatment",
            "p_value_vs_genome_control",
            "p_value_between",
            "q_value_vs_genome_treatment",
            "q_value_vs_genome_control",
            "q_value_between",
            "p_value_raw_vs_genome_treatment",
            "q_value_raw_vs_genome_treatment",
            "p_value_raw_between",
            "q_value_raw_between",
            "lrt_vs_genome_treatment",
            "p_value_chi2_vs_genome_treatment",
            "p_value_f_vs_genome_treatment",
            "p_value_raw_chi2_vs_genome_treatment",
            "p_value_raw_f_vs_genome_treatment",
            "lrt_between",
            "p_value_chi2_between",
            "p_value_f_between",
            "p_value_raw_chi2_between",
            "p_value_raw_f_between",
            "q_value_chi2_between",
            "q_value_f_between",
            "q_value_raw_chi2_between",
            "q_value_raw_f_between",
            "q_value_f_vs_genome_treatment",
        ):
            assert math.isnan(y[col]), col
        assert (y["note_treatment"], y["note_control"], y["note_between"]) == (
            "no_events",
            "no_events",
            "no_events",
        )
        # gZ: one control event in cage E; control's genome total there is 1, so
        # gZ's expected in E is 0.1 and 0 elsewhere -> rr 10 on one informative cage
        z = out.set_index("gene_id").loc["gZ"]
        assert z["n_de_novo_control"] == 1 and z["note_control"] == ""
        assert z["rate_ratio_treatment_vs_control"] == 0.0
        assert z["min_covered_fraction"] == 1.0

    def test_bh_runs_per_family_over_the_tested_genes_only(self):
        gc, tot = _plant_fixture()
        out = every_gene(gc, tot, GROUPS, min_cages=4, fdr_by_mag=False)
        p = out["p_value_between"]
        tested = p.notna()
        expected_q = multipletests(p[tested], method="fdr_bh")[1]
        assert out.loc[tested, "q_value_between"].to_numpy() == pytest.approx(
            expected_q
        )
        assert out.loc[~tested, "q_value_between"].isna().all()
        # the within-group family is corrected on its own
        pa = out["p_value_vs_genome_treatment"]
        assert out.loc[
            pa.notna(), "q_value_vs_genome_treatment"
        ].to_numpy() == pytest.approx(multipletests(pa[pa.notna()], method="fdr_bh")[1])

    def test_cage_floor_blocks_a_group_with_too_few_entered_cages(self):
        gc, tot = _plant_fixture()
        # control cage H never read gX: three entered control cages for gX
        gc.loc[
            (gc.gene_id == "gX") & (gc.replicate == "H"), ["covered_bp", "entered"]
        ] = [0, False]
        out = every_gene(gc, tot, GROUPS, min_cages=4, fdr_by_mag=False)
        x = out.set_index("gene_id").loc["gX"]
        assert x["n_cages_control"] == 3
        assert x["note_control"] == "below_cage_floor"
        assert math.isnan(x["p_value_vs_genome_control"])
        # Test B needs both groups at the floor
        assert x["note_between"] == "below_cage_floor"
        assert math.isnan(x["p_value_between"]) and math.isnan(x["q_value_between"])
        # treatment still has its four cages
        assert x["note_treatment"] == ""
        # a floor of 3 lets it through again
        out3 = every_gene(gc, tot, GROUPS, min_cages=3, fdr_by_mag=False)
        assert out3.set_index("gene_id").loc["gX", "note_between"] == ""

    def test_fdr_by_mag_corrects_each_genome_separately(self):
        gc1, tot1 = _plant_fixture()
        gc2, tot2 = _plant_fixture(
            hot_events=(1, 1, 0, 0)
        )  # 2 vs 0: p 0.096, not M's 0.0185
        gc2["mag_id"] = tot2["mag_id"] = "N"
        gc2["gene_id"] = gc2["gene_id"] + "2"
        gc = pd.concat([gc1, gc2], ignore_index=True)
        tot = pd.concat([tot1, tot2], ignore_index=True)
        pooled = every_gene(gc, tot, GROUPS, min_cages=4, fdr_by_mag=False)
        by_mag = every_gene(gc, tot, GROUPS, min_cages=4, fdr_by_mag=True)
        p = pooled["p_value_between"]
        assert pooled.loc[p.notna(), "q_value_between"].to_numpy() == pytest.approx(
            multipletests(p[p.notna()], method="fdr_bh")[1]
        )
        for mag in ("M", "N"):
            sub = by_mag[by_mag.mag_id == mag]
            ps = sub["p_value_between"]
            assert sub.loc[ps.notna(), "q_value_between"].to_numpy() == pytest.approx(
                multipletests(ps[ps.notna()], method="fdr_bh")[1]
            )
        assert not np.allclose(
            pooled["q_value_between"].fillna(-1), by_mag["q_value_between"].fillna(-1)
        )

    def test_min_covered_fraction_is_over_entered_cages(self):
        gc, tot = _plant_fixture()
        gc.loc[
            (gc.gene_id == "gX") & (gc.replicate == "A"),
            ["covered_bp", "covered_fraction"],
        ] = [40, 0.4]
        gc.loc[
            (gc.gene_id == "gX") & (gc.replicate == "H"),
            ["covered_bp", "entered", "covered_fraction"],
        ] = [0, False, 0.0]
        out = every_gene(gc, tot, GROUPS, min_cages=3, fdr_by_mag=False)
        assert out.set_index("gene_id").loc[
            "gX", "min_covered_fraction"
        ] == pytest.approx(0.4)

    def test_unknown_group_name_raises(self):
        gc, tot = _plant_fixture()
        with pytest.raises(ValueError, match="fat"):
            every_gene(gc, tot, ("fat", "control"), min_cages=4, fdr_by_mag=False)

    def test_no_genes_gives_an_empty_frame_with_the_columns(self):
        gc, tot = _plant_fixture()
        out = every_gene(gc.iloc[0:0], tot, GROUPS, min_cages=4, fdr_by_mag=False)
        assert out.empty and "p_value_between" in out.columns


class TestResultColumns:
    def test_raw_columns_follow_the_discounted_ones_in_every_block(self):
        cols = result_columns(("a", "b"))
        i = cols.index("q_value_vs_genome_a")
        assert cols[i + 1 : i + 3] == [
            "p_value_raw_vs_genome_a",
            "q_value_raw_vs_genome_a",
        ]
        j = cols.index("q_value_between")
        assert cols[j + 1 : j + 3] == ["p_value_raw_between", "q_value_raw_between"]

    def test_statistic_and_both_curves_sit_after_the_q_values_in_every_block(self):
        """House picks (with q) first, then the diagnostic twins: the statistic
        on both curves, discounted and raw, then the note."""
        cols = result_columns(("a", "b"))
        k = cols.index("dispersion_a")
        assert cols[k + 1 : k + 3] == ["df_resid_vs_genome_a", "lrt_vs_genome_a"]
        i = cols.index("q_value_raw_vs_genome_a")
        assert cols[i + 1 : i + 10] == [
            "p_value_chi2_vs_genome_a",
            "q_value_chi2_vs_genome_a",
            "p_value_f_vs_genome_a",
            "q_value_f_vs_genome_a",
            "p_value_raw_chi2_vs_genome_a",
            "q_value_raw_chi2_vs_genome_a",
            "p_value_raw_f_vs_genome_a",
            "q_value_raw_f_vs_genome_a",
            "note_a",
        ]
        m = cols.index("df_resid_between")
        assert cols[m + 1] == "lrt_between"
        j = cols.index("q_value_raw_between")
        assert cols[j + 1 : j + 10] == [
            "p_value_chi2_between",
            "q_value_chi2_between",
            "p_value_f_between",
            "q_value_f_between",
            "p_value_raw_chi2_between",
            "q_value_raw_chi2_between",
            "p_value_raw_f_between",
            "q_value_raw_f_between",
            "note_between",
        ]
