"""Tests for the gene-level de novo accumulation counter (gene_de_novo.py).

Everything here is a pure function over small hand-built frames; the CLI
(Task 6) only streams real files through the same functions.  Fixture sizes
are chosen so every expected number can be checked by eye:

* ``cage_rises``: a mean-change table of three positions -> the letters whose
  cage-average share went UP, nothing else;
* ``paired_samples``: which samples belong to mice sampled at both timepoints
  (the mice the pipeline's mean-change table is built from, hence the ONLY
  mice whose reads count as "searched");
* ``CoverageMask``: per cage and timepoint, a bit per contig position set
  where some sample has at least ``min_cov`` reads; readable = set at both
  timepoints;
* ``covered_bp_by_gene`` / ``covered_bp_total``: readable positions inside a
  gene's span, and over the whole genome.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from alleleflux.scripts.analysis.gene_de_novo import (
    CoverageMask,
    cage_rises,
    covered_bp_by_gene,
    covered_bp_total,
    paired_samples,
)

BASES = ("A", "T", "G", "C")


def _mean_changes(rows):
    """rows: (contig, gene_id, position, replicate, group, dA, dT, dG, dC)."""
    cols = ["contig", "gene_id", "position", "replicate", "group"] + [
        f"{b}_frequency_diff_mean" for b in BASES
    ]
    df = pd.DataFrame(rows, columns=cols)
    df["subjectID_count"] = 1
    # the pipeline writes these two as pandas categoricals; the function must
    # hand back plain text whatever the category order was
    df["replicate"] = df["replicate"].astype("category")
    df["group"] = df["group"].astype("category")
    return df


class TestCageRises:
    def test_only_positive_mean_changes_become_rises(self):
        mc = _mean_changes(
            [
                ("c1", "g1", 10, "cageA", "fat", -0.3, 0.3, 0.0, 0.0),
                ("c1", "g1", 20, "cageA", "fat", 0.0, 0.0, 0.0, 0.0),
                ("c1", None, 30, "cageB", "ctl", 0.1, -0.1, 0.0, 0.0),
            ]
        )
        out = cage_rises(mc, ("fat", "ctl"))
        assert list(zip(out.position, out.replicate, out.allele)) == [
            (10, "cageA", "T"),
            (30, "cageB", "A"),
        ]
        assert out.loc[out.position == 30, "gene_id"].item() == ""
        assert out.loc[out.position == 10, "diff_mean"].item() == pytest.approx(0.3)
        assert out["replicate"].dtype == object and out["group"].dtype == object

    def test_two_letters_rising_at_one_position_give_two_rows(self):
        # a three-letter site: A fell, both G and T rose
        mc = _mean_changes([("c1", "g1", 10, "cageA", "fat", -0.5, 0.2, 0.3, 0.0)])
        out = cage_rises(mc, ("fat", "ctl"))
        assert sorted(out.allele) == ["G", "T"]

    def test_rows_of_other_groups_are_dropped(self):
        mc = _mean_changes(
            [
                ("c1", "g1", 10, "cageA", "fat", -0.3, 0.3, 0.0, 0.0),
                ("c1", "g1", 10, "cageZ", "other", -0.3, 0.3, 0.0, 0.0),
            ]
        )
        out = cage_rises(mc, ("fat", "ctl"))
        assert list(out.replicate) == ["cageA"]

    def test_group_names_compare_as_text(self):
        mc = _mean_changes([("c1", "g1", 10, "cageA", 40, -0.3, 0.3, 0.0, 0.0)])
        out = cage_rises(mc, ("40", "AL"))
        assert len(out) == 1 and out.group.item() == "40"

    def test_missing_diff_column_raises(self):
        mc = _mean_changes([("c1", "g1", 10, "cageA", "fat", -0.3, 0.3, 0.0, 0.0)])
        with pytest.raises(ValueError, match="G_frequency_diff_mean"):
            cage_rises(mc.drop(columns=["G_frequency_diff_mean"]), ("fat", "ctl"))


def _sheet(rows):
    """rows: (sample_id, subjectID, replicate, group, time)."""
    return pd.DataFrame(
        rows, columns=["sample_id", "subjectID", "replicate", "group", "time"]
    )


class TestPairedSamples:
    def test_only_mice_with_both_timepoints_are_kept(self):
        sheet = _sheet(
            [
                ("m1_pre", "m1", "cageA", "fat", "pre"),
                ("m1_post", "m1", "cageA", "fat", "post"),
                ("m2_pre", "m2", "cageA", "fat", "pre"),  # no post sample
                ("m3_post", "m3", "cageB", "fat", "post"),  # no pre sample
                ("m4_pre", "m4", "cageB", "ctl", "pre"),
                ("m4_mid", "m4", "cageB", "ctl", "mid"),  # other timepoint, ignored
                ("m4_post", "m4", "cageB", "ctl", "post"),
            ]
        )
        out = paired_samples(sheet, "pre", "post")
        assert sorted(out.sample_id) == ["m1_post", "m1_pre", "m4_post", "m4_pre"]
        assert list(out.columns) == [
            "sample_id",
            "subjectID",
            "replicate",
            "group",
            "time",
        ]

    def test_a_mouse_sampled_twice_at_one_timepoint_raises(self):
        sheet = _sheet(
            [
                ("m1_pre", "m1", "cageA", "fat", "pre"),
                ("m1_pre_b", "m1", "cageA", "fat", "pre"),
                ("m1_post", "m1", "cageA", "fat", "post"),
            ]
        )
        with pytest.raises(ValueError, match="m1"):
            paired_samples(sheet, "pre", "post")


def _cache(rows, group="fat"):
    """rows: (sample_id, replicate, timepoint, contig, position, total_coverage)."""
    df = pd.DataFrame(
        rows,
        columns=[
            "sample_id",
            "replicate",
            "timepoint",
            "contig",
            "position",
            "total_coverage",
        ],
    )
    df["group"] = group
    return df


class TestCoverageMask:
    def test_readable_needs_a_deep_sample_at_both_timepoints(self):
        mask = CoverageMask({"c1": 12}, min_cov=5)
        mask.add(
            _cache(
                [
                    ("m1_pre", "cageA", "pre", "c1", 10, 5),
                    ("m1_post", "cageA", "post", "c1", 10, 4),
                    (
                        "m2_post",
                        "cageA",
                        "post",
                        "c1",
                        10,
                        7,
                    ),  # another mouse rescues post
                    ("m1_pre", "cageA", "pre", "c1", 11, 9),  # no post row at all
                    (
                        "m3_pre",
                        "cageB",
                        "pre",
                        "c1",
                        10,
                        50,
                    ),  # cage B has no post sample
                ]
            )
        )
        readable = mask.readable("pre", "post")
        assert set(readable) == {("fat", "cageA"), ("fat", "cageB")}
        assert np.flatnonzero(readable[("fat", "cageA")]["c1"]).tolist() == [10]
        assert not readable[("fat", "cageB")]["c1"].any()

    def test_the_same_replicate_letter_in_two_groups_is_two_cages(self):
        # the bundled data names the cages of both groups A..D
        mask = CoverageMask({"c1": 4}, min_cov=5)
        mask.add(_cache([("f1_pre", "A", "pre", "c1", 1, 6)], group="fat"))
        mask.add(_cache([("f1_post", "A", "post", "c1", 1, 6)], group="fat"))
        mask.add(_cache([("c1_pre", "A", "pre", "c1", 2, 6)], group="ctl"))
        mask.add(_cache([("c1_post", "A", "post", "c1", 2, 6)], group="ctl"))
        readable = mask.readable("pre", "post")
        assert np.flatnonzero(readable[("fat", "A")]["c1"]).tolist() == [1]
        assert np.flatnonzero(readable[("ctl", "A")]["c1"]).tolist() == [2]

    def test_rows_arrive_in_batches_and_accumulate(self):
        mask = CoverageMask({"c1": 4}, min_cov=5)
        mask.add(_cache([("m1_pre", "cageA", "pre", "c1", 2, 6)]))
        mask.add(_cache([("m1_post", "cageA", "post", "c1", 2, 6)]))
        mask.add(_cache([("m1_post", "cageA", "post", "c1", 3, 6)]))
        readable = mask.readable("pre", "post")
        assert np.flatnonzero(readable[("fat", "cageA")]["c1"]).tolist() == [2]

    def test_position_beyond_the_contig_raises(self):
        mask = CoverageMask({"c1": 4}, min_cov=5)
        with pytest.raises(ValueError, match="c1"):
            mask.add(_cache([("m1_pre", "cageA", "pre", "c1", 4, 6)]))

    def test_unknown_contig_raises(self):
        mask = CoverageMask({"c1": 4}, min_cov=5)
        with pytest.raises(ValueError, match="c9"):
            mask.add(_cache([("m1_pre", "cageA", "pre", "c9", 0, 6)]))


def _genes(rows):
    """rows: (gene_id, contig, start, end) with 0-based inclusive coordinates."""
    return pd.DataFrame(rows, columns=["gene_id", "contig", "start", "end"])


class TestCoveredBp:
    def test_covered_bp_counts_readable_positions_inside_the_gene(self):
        readable = {("fat", "cageA"): {"c1": np.zeros(300, dtype=bool)}}
        readable[("fat", "cageA")]["c1"][[99, 100, 150, 201]] = True
        genes = _genes([("g1", "c1", 100, 200)])
        out = covered_bp_by_gene(readable, genes)
        assert out.to_dict("records") == [
            {"group": "fat", "replicate": "cageA", "gene_id": "g1", "covered_bp": 2}
        ]
        assert covered_bp_total(readable).to_dict("records") == [
            {"group": "fat", "replicate": "cageA", "covered_bp_total": 4}
        ]

    def test_every_gene_and_cage_gets_a_row_even_at_zero(self):
        readable = {
            ("fat", "A"): {
                "c1": np.ones(10, dtype=bool),
                "c2": np.zeros(10, dtype=bool),
            },
            ("ctl", "A"): {
                "c1": np.zeros(10, dtype=bool),
                "c2": np.zeros(10, dtype=bool),
            },
        }
        genes = _genes([("g1", "c1", 0, 9), ("g2", "c2", 0, 4)])
        out = covered_bp_by_gene(readable, genes)
        key = ["group", "replicate", "gene_id"]
        assert out.set_index(key)["covered_bp"].to_dict() == {
            ("fat", "A", "g1"): 10,
            ("fat", "A", "g2"): 0,
            ("ctl", "A", "g1"): 0,
            ("ctl", "A", "g2"): 0,
        }
        totals = covered_bp_total(readable).set_index(["group", "replicate"])
        assert totals["covered_bp_total"].to_dict() == {
            ("fat", "A"): 10,
            ("ctl", "A"): 0,
        }

    def test_gene_on_a_contig_the_mask_does_not_know_raises(self):
        readable = {("fat", "cageA"): {"c1": np.ones(10, dtype=bool)}}
        with pytest.raises(ValueError, match="c7"):
            covered_bp_by_gene(readable, _genes([("g1", "c7", 0, 4)]))


# ---------------------------------------------------------------------------
# Task 2: per-cage classification through the Phase 2 rules
# ---------------------------------------------------------------------------

from alleleflux.scripts.analysis.ani.null_model import build_error_model
from alleleflux.scripts.analysis.baseline_presence import Timepoints
from alleleflux.scripts.analysis.gene_de_novo import classify_cage_rises, is_event

MODEL = build_error_model(min_base_quality=30, fdr=1e-6)
TPS = Timepoints("pre", "post")
PERIOD = "pre_post"
# Three cages: fat/cageA holds m1 + m2, fat/cageB holds m3, ctl/cageC holds m9.
MICE = {
    "m1": ("cageA", "fat"),
    "m2": ("cageA", "fat"),
    "m3": ("cageB", "fat"),
    "m9": ("cageC", "ctl"),
}
DATES = {"pre": "2024-01-15", "post": "2024-04-15"}


def _sheet_rows():
    rows = []
    for mouse, (rep, group) in MICE.items():
        for tp in ("pre", "post"):
            rows.append((f"{mouse}_{tp}", mouse, rep, group, tp, DATES[tp]))
    return pd.DataFrame(
        rows,
        columns=[
            "sample_id",
            "subjectID",
            "replicate",
            "group",
            "time",
            "collection_date",
        ],
    )


def _rows(counts):
    """counts: {(mouse, tp, position): (A, C, G, T)} -> cache rows for contig c1."""
    recs = []
    for (mouse, tp, pos), (a, c, g, t) in counts.items():
        rep, group = MICE[mouse]
        recs.append(
            {
                "contig": "c1",
                "position": pos,
                "sample_id": f"{mouse}_{tp}",
                "subjectID": mouse,
                "replicate": rep,
                "group": group,
                "timepoint": tp,
                "A": a,
                "C": c,
                "G": g,
                "T": t,
                "total_coverage": a + c + g + t,
            }
        )
    return pd.DataFrame(recs)


def _background(positions, mice=MICE, tps=("pre", "post"), depth=40):
    """Every mouse at every timepoint: ``depth`` reads of A, nothing else."""
    return {(m, tp, p): (depth, 0, 0, 0) for m in mice for tp in tps for p in positions}


def _rise(rep, group, pos, allele="T"):
    return {
        "contig": "c1",
        "position": pos,
        "gene_id": "g1",
        "replicate": rep,
        "group": group,
        "allele": allele,
        "diff_mean": 0.5,
        "subjectID_count": 1,
    }


def _turnover(verdicts):
    """verdicts: {mouse: True/False/None} -> the per-mouse turnover table."""
    return pd.DataFrame(
        {
            "subjectID": list(verdicts),
            "transition": PERIOD,
            "sample_t2": [f"{m}_post" for m in verdicts],
            "dominant_strain_change": [
                pd.NA if v is None else v for v in verdicts.values()
            ],
        }
    )


def _classify(rises, counts, **over):
    kwargs = dict(
        mag_id="M",
        tps=TPS,
        scope="cohort",
        model=MODEL,
        min_freq=0.05,
        min_cov=5,
        window_bp=1000,
        max_neighbours=0,
        strain_gate="none",
        min_determined_carriers=1,
        turnover=None,
    )
    kwargs.update(over)
    return classify_cage_rises(
        pd.DataFrame(rises), _rows(counts), _sheet_rows(), **kwargs
    )


class TestClassifyCageRises:
    def test_standing_allele_is_not_an_event(self):
        counts = _background([2000])
        counts[("m9", "pre", 2000)] = (
            34,
            0,
            0,
            6,
        )  # another group's mouse had T before
        counts[("m1", "post", 2000)] = (20, 0, 0, 20)
        out = _classify([_rise("cageA", "fat", 2000)], counts)
        assert out["origin"].tolist() == ["standing_variation"]
        assert out["carriers"].item() == "m1"
        assert out["n_baseline_samples_present"].item() == 1
        assert not is_event(out["origin"], "exclude").any()

    def test_lone_de_novo_in_one_cage_is_an_event(self):
        counts = _background([3000])
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        out = _classify([_rise("cageA", "fat", 3000)], counts)
        assert out["origin"].tolist() == ["de_novo"]
        assert out["block_id"].item() == ""
        assert out["n_baseline_samples_covered"].item() == 4  # every pre row (cohort)
        assert is_event(out["origin"], "exclude").tolist() == [True]

    def test_a_whisper_at_baseline_is_de_novo_below_detection(self):
        counts = _background([3000])
        counts[("m2", "pre", 3000)] = (39, 0, 0, 1)  # 1 read of 40: under the bar
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        out = _classify([_rise("cageA", "fat", 3000)], counts)
        assert out["origin"].item() == "de_novo_below_detection"
        assert is_event(out["origin"], "exclude").item()

    def test_two_rises_within_the_window_in_one_cage_are_a_block(self):
        counts = _background([1172, 1472])
        counts[("m1", "post", 1172)] = (20, 0, 0, 20)
        counts[("m1", "post", 1472)] = (20, 0, 0, 20)
        out = _classify(
            [_rise("cageA", "fat", 1172), _rise("cageA", "fat", 1472)], counts
        )
        assert out["origin"].tolist() == ["block_event", "block_event"]
        assert out["block_id"].nunique() == 1 and out["block_id"].iloc[0] != ""
        assert out["n_de_novo_neighbours"].tolist() == [1, 1]
        assert is_event(out["origin"], "exclude").tolist() == [False, False]
        assert is_event(out["origin"], "once").tolist() == [True, True]

    def test_the_same_positions_in_different_cages_are_not_neighbours(self):
        counts = _background([1172, 1472])
        counts[("m1", "post", 1172)] = (20, 0, 0, 20)
        counts[("m3", "post", 1472)] = (20, 0, 0, 20)
        out = _classify(
            [_rise("cageA", "fat", 1172), _rise("cageB", "fat", 1472)], counts
        )
        assert out["origin"].tolist() == ["de_novo", "de_novo"]
        assert out["block_id"].tolist() == ["", ""]

    def test_rise_without_a_present_carrier_is_not_present_later(self):
        counts = _background([3000])
        counts[("m1", "post", 3000)] = (39, 0, 0, 1)  # one read: below the bar
        out = _classify([_rise("cageA", "fat", 3000)], counts)
        assert out["origin"].item() == "not_present_later"
        assert out["carriers"].item() == ""
        assert not is_event(out["origin"], "once").item()

    def test_baseline_nobody_covered_is_not_an_event(self):
        counts = _background([3000], tps=("post",))
        counts.update(
            {(m, "pre", 3000): (2, 0, 0, 0) for m in MICE}
        )  # all under min_cov
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        out = _classify([_rise("cageA", "fat", 3000)], counts)
        assert out["origin"].item() == "baseline_not_covered"
        assert out["n_baseline_samples_covered"].item() == 0

    def test_carrier_gate_one_clean_carrier_keeps_the_event(self):
        counts = _background([3000])
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        counts[("m2", "post", 3000)] = (20, 0, 0, 20)
        out = _classify(
            [_rise("cageA", "fat", 3000)],
            counts,
            strain_gate="carrier",
            turnover=_turnover({"m1": True, "m2": False}),
        )
        assert out["origin"].item() == "de_novo"
        assert out["carriers"].item() == "m1,m2"
        assert (
            out["n_carriers_determined"].item(),
            out["n_carriers_strain_changed"].item(),
        ) == (2, 1)

    def test_carrier_gate_every_carrier_changed_strain_is_new_strain(self):
        counts = _background([3000])
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        out = _classify(
            [_rise("cageA", "fat", 3000)],
            counts,
            strain_gate="carrier",
            turnover=_turnover({"m1": True, "m2": False}),
        )
        assert out["origin"].item() == "new_strain"
        assert not is_event(out["origin"], "once").item()

    def test_carrier_gate_without_a_verdict_is_undetermined_not_counted(self):
        counts = _background([3000])
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        for turnover in (None, _turnover({"m1": None}), _turnover({"m2": False})):
            out = _classify(
                [_rise("cageA", "fat", 3000)],
                counts,
                strain_gate="carrier",
                turnover=turnover,
            )
            assert out["origin"].item() == "strain_undetermined"
            assert not is_event(out["origin"], "once").item()

    def test_gate_none_leaves_a_de_novo_event_alone(self):
        counts = _background([3000])
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        out = _classify([_rise("cageA", "fat", 3000)], counts, strain_gate="none")
        assert out["origin"].item() == "de_novo"
        assert out["carrier_strain_status"].item() == ""

    @pytest.mark.parametrize("gate", ["group", "both"])
    def test_group_gates_are_refused(self, gate):
        counts = _background([3000])
        with pytest.raises(ValueError, match="strain_gate"):
            _classify([_rise("cageA", "fat", 3000)], counts, strain_gate=gate)

    def test_scope_decides_which_earlier_samples_count(self):
        # m3 (same group, other cage) had T before; m9 (other group) did not.
        counts = _background([3000])
        counts[("m3", "pre", 3000)] = (34, 0, 0, 6)
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        rise = [_rise("cageA", "fat", 3000)]
        assert (
            _classify(rise, counts, scope="cohort")["origin"].item()
            == "standing_variation"
        )
        assert (
            _classify(rise, counts, scope="comparison")["origin"].item()
            == "standing_variation"
        )
        assert (
            _classify(rise, counts, scope="own_group")["origin"].item()
            == "standing_variation"
        )
        assert _classify(rise, counts, scope="individual")["origin"].item() == "de_novo"
        # the carrier's own earlier sample had a whisper: individual says below detection
        counts[("m1", "pre", 3000)] = (39, 0, 0, 1)
        assert (
            _classify(rise, counts, scope="individual")["origin"].item()
            == "de_novo_below_detection"
        )

    def test_output_keeps_the_rise_columns_and_adds_the_verdict_columns(self):
        counts = _background([3000])
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        out = _classify([_rise("cageA", "fat", 3000)], counts)
        assert list(out.columns) == [
            "contig",
            "position",
            "gene_id",
            "replicate",
            "group",
            "allele",
            "diff_mean",
            "subjectID_count",
            "origin",
            "block_id",
            "n_de_novo_neighbours",
            "carriers",
            "n_carriers",
            "carrier_strain_status",
            "n_carriers_determined",
            "n_carriers_strain_changed",
            "n_baseline_samples_covered",
            "n_baseline_samples_present",
        ]

    def test_no_rises_gives_an_empty_frame_with_the_columns(self):
        out = _classify([], _background([3000]))
        assert out.empty and "origin" in out.columns


class TestIsEvent:
    def test_only_de_novo_flavours_count_and_blocks_follow_the_knob(self):
        origins = pd.Series(
            [
                "de_novo",
                "de_novo_below_detection",
                "block_event",
                "standing_variation",
                "new_strain",
                "strain_undetermined",
                "not_present_later",
                "baseline_not_covered",
            ]
        )
        assert is_event(origins, "exclude").tolist() == [
            True,
            True,
            False,
            False,
            False,
            False,
            False,
            False,
        ]
        assert is_event(origins, "once").tolist() == [
            True,
            True,
            True,
            False,
            False,
            False,
            False,
            False,
        ]

    def test_unknown_block_rule_raises(self):
        with pytest.raises(ValueError, match="block_rule"):
            is_event(pd.Series(["de_novo"]), "twice")


class TestCohortWithoutDates:
    """A sheet with no collection_date cannot order a third timepoint against
    the comparison, and the sheet's row order says nothing (user 2026-11-02).
    Under ``cohort`` the baseline then narrows to the comparison's EARLIER
    timepoint, exactly like ``comparison``; a 'mid' sample listed first in the
    sheet must not become baseline just because it is listed first."""

    def _undated_sheet_mid_first(self):
        sheet = _sheet_rows().drop(columns=["collection_date"])
        mid = pd.DataFrame(
            [("m7_mid", "m7", "cageZ", "other", "mid")],
            columns=["sample_id", "subjectID", "replicate", "group", "time"],
        )
        return pd.concat([mid, sheet], ignore_index=True)  # mid listed FIRST

    def test_undated_cohort_uses_only_the_earlier_timepoint(self):
        counts = _background([3000])
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        rows = _rows(counts)
        # a 'mid' sample of another group carries T; without dates nobody can
        # say whether mid is before post, so it must not be baseline
        mid_row = rows.iloc[[0]].copy()
        mid_row.loc[
            :, ["sample_id", "subjectID", "replicate", "group", "timepoint"]
        ] = ["m7_mid", "m7", "cageZ", "other", "mid"]
        mid_row.loc[:, ["A", "C", "G", "T", "total_coverage"]] = [30, 0, 0, 10, 40]
        rows = pd.concat([rows, mid_row], ignore_index=True)
        out = classify_cage_rises(
            pd.DataFrame([_rise("cageA", "fat", 3000)]),
            rows,
            self._undated_sheet_mid_first(),
            mag_id="M",
            tps=TPS,
            scope="cohort",
            model=MODEL,
            min_freq=0.05,
            min_cov=5,
            window_bp=1000,
            max_neighbours=0,
            strain_gate="none",
            min_determined_carriers=1,
            turnover=None,
        )
        assert out["origin"].item() == "de_novo"
        assert out["n_baseline_samples_covered"].item() == 4  # the four pre rows only
