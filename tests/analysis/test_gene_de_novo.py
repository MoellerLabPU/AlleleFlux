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
import pyarrow.parquet as pq
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


FIRST_SEEN = "first_seen_after_treatment_began"
MID_DATE = "2024-03-01"  # between pre (2024-01-15) and post (2024-04-15)


def _sheet_with_mid():
    """``_sheet_rows`` plus a sample of every mouse at a third timepoint
    ``mid``, dated between pre and post."""
    mid = pd.DataFrame(
        [
            (f"{m}_mid", m, rep, group, "mid", MID_DATE)
            for m, (rep, group) in MICE.items()
        ],
        columns=[
            "sample_id",
            "subjectID",
            "replicate",
            "group",
            "time",
            "collection_date",
        ],
    )
    return pd.concat([_sheet_rows(), mid], ignore_index=True)


def _classify_three_timepoints(rises, counts, pre_treatment):
    """``_classify`` on the three-timepoint sheet, cohort scope, blocks and gate off."""
    return classify_cage_rises(
        pd.DataFrame(rises),
        _rows(counts),
        _sheet_with_mid(),
        mag_id="M",
        tps=TPS,
        scope="cohort",
        model=MODEL,
        min_freq=0.05,
        min_cov=5,
        window_bp=1000,
        max_neighbours=10**6,
        strain_gate="none",
        min_determined_carriers=1,
        turnover=None,
        pre_treatment=pre_treatment,
    )


class TestFirstSeenAfterTreatment:
    """Pre-treatment timepoint: ``pre``.  ``mid`` (dated between pre and
    post) was sampled after treatment began.  T rises in cage A at 3000; m1
    carries it at post."""

    def _counts(self, mid_carrier=None, pre_carrier=None):
        counts = _background([3000], tps=("pre", "mid", "post"))
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        if mid_carrier:
            counts[(mid_carrier, "mid", 3000)] = (30, 0, 0, 10)
        if pre_carrier:
            counts[(pre_carrier, "pre", 3000)] = (30, 0, 0, 10)
        return counts

    def test_only_an_on_treatment_witness_carries_it(self):
        """m9 (another group) has T at mid and nobody has it at pre."""
        out = _classify_three_timepoints(
            [_rise("cageA", "fat", 3000)], self._counts(mid_carrier="m9"), ["pre"]
        )
        assert out["origin"].tolist() == [FIRST_SEEN]
        assert out["n_baseline_samples_present"].item() == 1
        assert out["n_baseline_samples_covered"].item() == 8  # 4 pre + 4 mid

    def test_a_pre_treatment_witness_carries_it_too(self):
        out = _classify_three_timepoints(
            [_rise("cageA", "fat", 3000)],
            self._counts(mid_carrier="m9", pre_carrier="m3"),
            ["pre"],
        )
        assert out["origin"].tolist() == ["standing_variation"]

    def test_no_witness_carries_it_is_still_de_novo(self):
        out = _classify_three_timepoints(
            [_rise("cageA", "fat", 3000)], self._counts(), ["pre"]
        )
        assert out["origin"].tolist() == ["de_novo"]

    def test_without_the_list_it_is_standing_variation(self):
        out = _classify_three_timepoints(
            [_rise("cageA", "fat", 3000)], self._counts(mid_carrier="m9"), None
        )
        assert out["origin"].tolist() == ["standing_variation"]

    def test_both_timepoints_pre_treatment_is_standing_variation(self):
        out = _classify_three_timepoints(
            [_rise("cageA", "fat", 3000)],
            self._counts(mid_carrier="m9"),
            ["pre", "mid"],
        )
        assert out["origin"].tolist() == ["standing_variation"]

    @pytest.mark.parametrize("block_rule", ["exclude", "once"])
    def test_it_is_never_an_event(self, block_rule):
        """Design D12: counting it would count one mutation again at every
        later timepoint and cage it spreads to."""
        assert not is_event(pd.Series([FIRST_SEEN]), block_rule).any()

    def test_it_never_joins_a_block_or_passes_the_strain_gate(self):
        """A de novo rise 100 bp away in the same cage: with the cluster rule
        ON (no neighbour tolerated) the de novo one stays a lone de novo and
        the first-seen one keeps its label, also under the carrier gate with
        a carrier that changed strain."""
        counts = _background([3000, 3100], tps=("pre", "mid", "post"))
        counts[("m1", "post", 3000)] = (20, 0, 0, 20)
        counts[("m9", "mid", 3000)] = (30, 0, 0, 10)  # witness for 3000 only
        counts[("m1", "post", 3100)] = (20, 0, 0, 20)
        out = classify_cage_rises(
            pd.DataFrame([_rise("cageA", "fat", 3000), _rise("cageA", "fat", 3100)]),
            _rows(counts),
            _sheet_with_mid(),
            mag_id="M",
            tps=TPS,
            scope="cohort",
            model=MODEL,
            min_freq=0.05,
            min_cov=5,
            window_bp=1000,
            max_neighbours=0,
            strain_gate="carrier",
            min_determined_carriers=1,
            turnover=_turnover({"m1": True}),
            pre_treatment=["pre"],
        )
        assert out["origin"].tolist() == [FIRST_SEEN, "new_strain"]
        assert out["block_id"].tolist() == ["", ""]


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


# ---------------------------------------------------------------------------
# Task 3: gene tables
# ---------------------------------------------------------------------------

from alleleflux.scripts.analysis.gene_de_novo import gene_tables, load_genes

# A 9-base gene: ATG GGA TGG = M G W; repeated to build the three test genes.
CODING = "ATGGGATGG"


class TestLoadGenes:
    def test_prodigal_headers_become_zero_based_rows_with_sequences(self, tmp_path):
        fasta = tmp_path / "genes.fna"
        fasta.write_text(
            ">c1_1 # 101 # 109 # 1 # ID=1_1;partial=00\n"
            + CODING
            + "\n>c1_2 # 201 # 209 # -1 # ID=1_2;partial=00\n"
            + CODING
            + "\n"
        )
        genes = load_genes(str(fasta))
        assert genes[["gene_id", "contig", "start", "end", "strand", "length"]].to_dict(
            "records"
        ) == [
            {
                "gene_id": "c1_1",
                "contig": "c1",
                "start": 100,
                "end": 108,
                "strand": 1,
                "length": 9,
            },
            {
                "gene_id": "c1_2",
                "contig": "c1",
                "start": 200,
                "end": 208,
                "strand": -1,
                "length": 9,
            },
        ]
        assert genes["sequence"].tolist() == [CODING, CODING]

    def test_a_sequence_whose_length_disagrees_with_the_header_raises(self, tmp_path):
        fasta = tmp_path / "genes.fna"
        fasta.write_text(">c1_1 # 101 # 112 # 1 # ID=1_1\n" + CODING + "\n")
        with pytest.raises(ValueError, match="c1_1"):
            load_genes(str(fasta))


def _genes_frame():
    # three genes on a 1,200 bp contig: g1 100..399 (+), g2 400..699 (-), g3 700..999 (+)
    seq = (CODING * 34)[:300]  # 300 bases = 100 codons; codon k starts at 3k
    return pd.DataFrame(
        [
            ("g1", "c1", 100, 399, 1, seq),
            ("g2", "c1", 400, 699, -1, seq),
            ("g3", "c1", 700, 999, 1, seq),
        ],
        columns=["gene_id", "contig", "start", "end", "strand", "sequence"],
    ).assign(length=300)


def _event(group, rep, pos, gene, alt, origin="de_novo"):
    return {
        "contig": "c1",
        "position": pos,
        "gene_id": gene,
        "replicate": rep,
        "group": group,
        "allele": alt,
        "origin": origin,
    }


def _readable(spec):
    """spec: {(group, rep): slice-or-None} -> readable dict over the 1,200 bp contig."""
    out = {}
    for cage, readable in spec.items():
        arr = np.zeros(1200, dtype=bool)
        if readable is not None:
            arr[readable] = True
        out[cage] = {"c1": arr}
    return out


CAGES = pd.DataFrame(
    [("fat", "A"), ("fat", "B"), ("ctl", "C")], columns=["group", "replicate"]
)


class TestGeneTables:
    def _run(self, events, readable, **over):
        kwargs = dict(mag_id="M", block_rule="exclude", min_gene_covered_fraction=0.0)
        kwargs.update(over)
        return gene_tables(
            pd.DataFrame(events), readable, _genes_frame(), CAGES, **kwargs
        )

    def test_counts_rates_and_totals_pin_the_arithmetic(self):
        # two events in g1 (150, 160), one between genes (1100), one in g2 (450)
        events = [
            _event("fat", "A", 150, "g1", "G"),
            _event("fat", "A", 160, "g1", "C"),
            _event("fat", "A", 1100, "", "T"),
            _event("fat", "B", 450, "g2", "G"),
        ]
        readable = _readable(
            {
                ("fat", "A"): slice(0, 1200),
                ("fat", "B"): slice(0, 600),
                ("ctl", "C"): slice(0, 1200),
            }
        )
        gene_cage, totals = self._run(events, readable)
        gc = gene_cage.set_index(["group", "replicate", "gene_id"])
        row = gc.loc[("fat", "A", "g1")]
        assert row.n_de_novo == 2
        assert (row.covered_bp, row.gene_length, row.covered_fraction) == (
            300,
            300,
            1.0,
        )
        assert row.rate == pytest.approx(2 / 300) and row.entered
        row = gc.loc[("fat", "B", "g2")]
        assert (row.n_de_novo, row.covered_bp) == (1, 200)  # 400..599 readable
        assert row.covered_fraction == pytest.approx(200 / 300)
        row = gc.loc[("fat", "B", "g3")]
        assert row.covered_bp == 0 and np.isnan(row.rate) and not row.entered
        assert len(gene_cage) == 9  # 3 genes x 3 cages, zeros kept
        assert gene_cage["n_de_novo"].sum() == 3  # the intergenic one is in no gene
        tot = totals.set_index(["group", "replicate"])
        assert tot.loc[("fat", "A")].to_dict() == {
            "mag_id": "M",
            "n_de_novo_total": 3,
            "covered_bp_total": 1200,
        }
        assert tot.loc[("fat", "B")].to_dict() == {
            "mag_id": "M",
            "n_de_novo_total": 1,
            "covered_bp_total": 600,
        }
        assert tot.loc[("ctl", "C")].to_dict() == {
            "mag_id": "M",
            "n_de_novo_total": 0,
            "covered_bp_total": 1200,
        }

    def test_block_rule_decides_whether_blocks_count(self):
        events = [
            _event("fat", "A", 150, "g1", "G", origin="block_event"),
            _event("fat", "A", 160, "g1", "C", origin="block_event"),
            _event("fat", "A", 180, "g1", "G", origin="standing_variation"),
        ]
        readable = _readable(
            {("fat", "A"): slice(0, 1200), ("fat", "B"): None, ("ctl", "C"): None}
        )
        gc, tot = self._run(events, readable)
        row = gc.set_index(["group", "replicate", "gene_id"]).loc[("fat", "A", "g1")]
        assert (row.n_de_novo, row.n_block_excluded) == (0, 2)
        assert (
            tot.set_index(["group", "replicate"]).loc[("fat", "A"), "n_de_novo_total"]
            == 0
        )
        gc, tot = self._run(events, readable, block_rule="once")
        row = gc.set_index(["group", "replicate", "gene_id"]).loc[("fat", "A", "g1")]
        assert (row.n_de_novo, row.n_block_excluded) == (2, 0)
        assert (
            tot.set_index(["group", "replicate"]).loc[("fat", "A"), "n_de_novo_total"]
            == 2
        )

    def test_covered_fraction_knob_drops_a_sliver(self):
        events = [_event("fat", "B", 450, "g2", "G")]
        readable = _readable(
            {("fat", "A"): None, ("fat", "B"): slice(0, 600), ("ctl", "C"): None}
        )
        gc, _ = self._run(events, readable, min_gene_covered_fraction=0.8)
        row = gc.set_index(["group", "replicate", "gene_id"]).loc[("fat", "B", "g2")]
        # 200 of 300 bases searched is under the 0.8 bar: the cage does not
        # enter for this gene, and its event is not counted
        assert not row.entered and row.covered_bp == 0 and row.n_de_novo == 0
        assert row.covered_fraction == pytest.approx(200 / 300)  # still reported
        assert np.isnan(row.rate)

    def test_cage_with_nothing_readable_appears_with_zeros(self):
        readable = _readable(
            {("fat", "A"): slice(0, 1200)}
        )  # B and C not in the mask at all
        gc, tot = self._run([], readable)
        assert len(gc) == 9
        sub = gc[(gc.group == "ctl") & (gc.replicate == "C")]
        assert (
            (sub.covered_bp == 0).all()
            and (~sub.entered).all()
            and sub.rate.isna().all()
        )
        assert tot.set_index(["group", "replicate"]).loc[("ctl", "C")].to_dict() == {
            "mag_id": "M",
            "n_de_novo_total": 0,
            "covered_bp_total": 0,
        }

    def test_an_event_in_an_unknown_gene_raises(self):
        with pytest.raises(ValueError, match="g9"):
            self._run(
                [_event("fat", "A", 150, "g9", "G")],
                _readable({("fat", "A"): slice(0, 1200)}),
            )


# ---------------------------------------------------------------------------
# Task 6: the classifier in bulk, and the CLI end to end
# ---------------------------------------------------------------------------

import os
import subprocess
import sys

from alleleflux.scripts.analysis.baseline_presence import (
    assess_rows,
    baseline_cutoff,
    baseline_rows,
)
from alleleflux.scripts.analysis.gene_de_novo import (
    CAGE_TOTAL_COLUMNS,
    EVENT_COLUMNS,
    GENE_CAGE_COLUMNS,
)
from alleleflux.scripts.statistics.poisson_rate_tests import result_columns


def _reference_origin(rises, rows, sheet, scope, pre_treatment=None):
    """Stages 1-3 the slow way: Phase 2's functions called once per candidate.

    With ``pre_treatment`` (timepoints sampled before treatment began) an
    allele some baseline row carries is standing only when one of those rows
    is at a pre-treatment timepoint; otherwise it is first seen after
    treatment began.

    The oracle for the bulk implementation; written out in full here so a
    change to ``baseline_rows`` / ``baseline_cutoff`` / ``assess_rows`` that
    the bulk path does not mirror shows up as a red test, not a silent drift.
    """
    rows = rows.copy()
    rows["time"] = rows["timepoint"]
    dates = sheet.set_index("sample_id")["collection_date"]
    rows["collection_date"] = dates.reindex(rows["sample_id"]).to_numpy()
    out = []
    for r in rises.itertuples(index=False):
        site = rows[(rows["contig"] == r.contig) & (rows["position"] == r.position)]
        site = pd.concat([site, assess_rows(site, r.allele, MODEL, 0.05, 5)], axis=1)
        later = site[
            (site["time"] == "post")
            & (site["group"] == r.group)
            & (site["replicate"] == r.replicate)
        ]
        carriers = set(later.loc[later["allele_present"], "subjectID"])
        if not carriers:
            out.append("not_present_later")
            continue
        pool = site[site["group"] == r.group] if scope == "own_group" else site
        t0 = baseline_rows(pool, scope, TPS, carriers, baseline_cutoff(later))
        t0 = t0[t0["allele_status"] != "not_covered"]
        if t0.empty:
            out.append("baseline_not_covered")
        elif t0["allele_present"].any():
            carrying = t0[t0["allele_present"]]
            if pre_treatment is None or carrying["time"].isin(pre_treatment).any():
                out.append("standing_variation")
            else:
                out.append("first_seen_after_treatment_began")
        elif (t0["allele_status"] == "below_detection").any():
            out.append("de_novo_below_detection")
        else:
            out.append("de_novo")
    return out


class TestClassifyInBulk:
    """The bulk classifier gives the per-candidate answer at every site."""

    @pytest.mark.parametrize(
        "scope", ["cohort", "comparison", "own_group", "individual"]
    )
    def test_bulk_verdicts_match_the_per_candidate_oracle(self, scope):
        rng = np.random.default_rng(7)
        positions = list(range(100, 160))
        counts = {}
        for m in MICE:
            for tp in ("pre", "post"):
                for p in positions:
                    # mostly A with a sprinkle of every other letter, so every
                    # verdict word (present / below / absent / not_covered)
                    # occurs somewhere in the grid
                    depth = int(rng.choice([2, 8, 20, 40, 60]))
                    other = rng.integers(0, 4, size=3)
                    a = max(depth - int(other.sum()), 0)
                    counts[(m, tp, p)] = (
                        a,
                        int(other[0]),
                        int(other[1]),
                        int(other[2]),
                    )
        rises = [
            _rise(rep, group, p, allele)
            for p in positions
            for (rep, group) in {("cageA", "fat"), ("cageB", "fat"), ("cageC", "ctl")}
            for allele in ("C", "G", "T")
        ]
        # blocks off (no cap on neighbours), gate off: stages 1-3 alone
        out = _classify(
            rises, counts, scope=scope, max_neighbours=10**6, strain_gate="none"
        )
        expected = _reference_origin(
            pd.DataFrame(rises), _rows(counts), _sheet_rows(), scope
        )
        assert out["origin"].tolist() == expected
        # the grid exercises more than one verdict word (not a trivially uniform case)
        assert len(set(out["origin"])) >= 3

    @pytest.mark.parametrize("pre_treatment", [["pre"], ["mid"], ["pre", "mid"]])
    def test_three_labels_match_the_oracle_on_a_three_timepoint_grid(
        self, pre_treatment
    ):
        """Random reads for every mouse at pre, mid and post (mid dated
        between the two), cohort scope.  The bulk classifier and the
        per-candidate oracle must give the same label at every candidate,
        whichever timepoints are called pre-treatment."""
        rng = np.random.default_rng(11)
        positions = list(range(100, 160))
        counts = {}
        for m in MICE:
            for tp in ("pre", "mid", "post"):
                for p in positions:
                    depth = int(rng.choice([2, 8, 20, 40, 60]))
                    # mostly clean piles, so that "nobody carried it" and
                    # "one class carried it" both occur often
                    other = rng.integers(0, 4, size=3) * (rng.random(3) < 0.2)
                    a = max(depth - int(other.sum()), 0)
                    counts[(m, tp, p)] = (
                        a,
                        int(other[0]),
                        int(other[1]),
                        int(other[2]),
                    )
        rises = [
            _rise(rep, group, p, allele)
            for p in positions
            for (rep, group) in {("cageA", "fat"), ("cageB", "fat"), ("cageC", "ctl")}
            for allele in ("C", "G", "T")
        ]
        out = _classify_three_timepoints(rises, counts, pre_treatment)
        expected = _reference_origin(
            pd.DataFrame(rises),
            _rows(counts),
            _sheet_with_mid(),
            "cohort",
            pre_treatment,
        )
        assert out["origin"].tolist() == expected
        labels = set(expected)
        assert "standing_variation" in labels
        if pre_treatment != ["pre", "mid"]:
            assert "first_seen_after_treatment_began" in labels
        else:  # every witness is pre-treatment: the third label cannot occur
            assert "first_seen_after_treatment_began" not in labels


# --- the CLI on a tiny run directory ---------------------------------------
GENOME = 1000  # one contig of this many bases
GENE_START1, GENE_END1 = 101, 400  # Prodigal's 1-based inclusive span
CLI_MICE = {  # mouse -> (cage, group); one mouse per cage, two cages per group
    "t1": ("A", "trt"),
    "t2": ("B", "trt"),
    "c1": ("C", "ctl"),
    "c2": ("D", "ctl"),
}
CLI_DATES = {"pre": "2024-01-15", "post": "2024-04-15"}


def _write_run(tmp, plant_control_standing=True):
    """A run directory with one MAG, one gene, and one planted event per treatment cage.

    Every sample reads 40 A at every base (readable everywhere); at 150 mouse
    t1's post sample and at 250 mouse t2's post sample carry 20 T (the
    plants); at 600 (between genes) cage C's post sample carries G which c1's
    pre sample already had (standing variation, not an event).
    """
    tmp = str(tmp)
    rng = np.random.default_rng(3)
    ref = "".join(rng.choice(list("ACGT"), size=GENOME))
    os.makedirs(os.path.join(tmp, "ref"))
    with open(os.path.join(tmp, "ref", "mags.fasta"), "w") as fh:
        fh.write(">M1_c1\n" + ref + "\n")
    with open(os.path.join(tmp, "ref", "genes.fna"), "w") as fh:
        fh.write(
            f">M1_c1_1 # {GENE_START1} # {GENE_END1} # 1 # ID=1_1;partial=00\n"
            + ref[GENE_START1 - 1 : GENE_END1]
            + "\n"
        )
    pd.DataFrame({"mag_id": ["M1"], "contig_id": ["M1_c1"]}).to_csv(
        os.path.join(tmp, "ref", "mag_mapping.tsv"), sep="\t", index=False
    )
    meta = pd.DataFrame(
        [
            (f"{m}_{tp}", m, group, cage, tp, CLI_DATES[tp])
            for m, (cage, group) in CLI_MICE.items()
            for tp in ("pre", "post")
        ],
        columns=[
            "sample_id",
            "subjectID",
            "group",
            "replicate",
            "time",
            "collection_date",
        ],
    )
    meta.to_csv(os.path.join(tmp, "metadata.tsv"), sep="\t", index=False)
    # the cache: one parquet per timepoint, every sample at every base
    plants = {("t1", "post", 150): "T", ("t2", "post", 250): "T"}
    if plant_control_standing:
        plants[("c1", "post", 600)] = "G"
        plants[("c1", "pre", 600)] = "G"
    for tp in ("pre", "post"):
        recs = []
        for m, (cage, group) in CLI_MICE.items():
            for p in range(GENOME):
                a, c, g, t = 40, 0, 0, 0
                letter = plants.get((m, tp, p))
                if letter == "T":
                    a, t = 20, 20
                elif letter == "G":
                    a, g = 20, 20
                recs.append(
                    dict(
                        contig="M1_c1",
                        position=p,
                        ref_base=ref[p],
                        gene_id=(
                            "M1_c1_1" if GENE_START1 - 1 <= p <= GENE_END1 - 1 else ""
                        ),
                        total_coverage=40,
                        A=a,
                        C=c,
                        G=g,
                        T=t,
                        sample_id=f"{m}_{tp}",
                        group=group,
                        subjectID=m,
                        replicate=cage,
                        A_frequency=a / 40,
                        T_frequency=t / 40,
                        G_frequency=g / 40,
                        C_frequency=c / 40,
                        time=tp,
                    )
                )
        df = pd.DataFrame(recs)
        for col in ("sample_id", "group", "subjectID", "replicate", "time"):
            df[col] = df[col].astype("category")
        d = os.path.join(tmp, "allele_freq_cache", tp)
        os.makedirs(d)
        df.to_parquet(os.path.join(d, f"M1_{tp}_allele_frequency.parquet"), index=False)
    # the mean-change table: the letters that rose per cage
    rows = [
        ("M1_c1", "M1_c1_1", 150, "A", "trt", -0.5, 0.5, 0.0, 0.0, 1),
        ("M1_c1", "M1_c1_1", 250, "B", "trt", -0.5, 0.5, 0.0, 0.0, 1),
        ("M1_c1", "", 600, "C", "ctl", -0.5, 0.0, 0.5, 0.0, 1),
        ("M1_c1", "M1_c1_1", 150, "C", "ctl", 0.0, 0.0, 0.0, 0.0, 1),  # flat: no rise
    ]
    mc = pd.DataFrame(
        rows,
        columns=["contig", "gene_id", "position", "replicate", "group"]
        + [f"{b}_frequency_diff_mean" for b in ("A", "T", "G", "C")]
        + ["subjectID_count"],
    )
    aa = os.path.join(tmp, "allele_analysis")
    os.makedirs(aa)
    mc.to_parquet(
        os.path.join(aa, "M1_allele_frequency_changes_mean.parquet"), index=False
    )
    # a turnover table: t1 changed strain, t2 kept it, controls kept it
    os.makedirs(os.path.join(tmp, "strain_turnover"))
    pd.DataFrame(
        {
            "MAG_ID": "M1",
            "subjectID": list(CLI_MICE),
            "transition": "pre_post",
            "sample_t2": [f"{m}_post" for m in CLI_MICE],
            "dominant_strain_change": [True, False, False, False],
        }
    ).to_csv(
        os.path.join(tmp, "strain_turnover", "M1_strain_turnover.tsv"),
        sep="\t",
        index=False,
    )
    return tmp


def _add_witness_sample(run, with_cache=True, empty_cache=False):
    """Add one sample OUTSIDE the comparison to a ``_write_run`` directory.

    Mouse w1 of a third group ``oth`` at a third timepoint ``mid``, collected
    before every other sample.  With ``with_cache`` its cache file is written:
    40 A at every base and 20 T at 150, the letter that later rises in cage A,
    so under the cohort scope it testifies that T was already there (or, with
    ``empty_cache``, the zero-row file the workflow writes for a genome that
    has no QC-passing sample at that timepoint).
    """
    from alleleflux.scripts.analysis.allele_frequency._allele_freq_common import (
        empty_cache_table,
    )

    sheet = os.path.join(run, "metadata.tsv")
    meta = pd.read_csv(sheet, sep="\t", dtype=str)
    witness = {
        "sample_id": "w1_mid",
        "subjectID": "w1",
        "group": "oth",
        "replicate": "W",
        "time": "mid",
        "collection_date": "2023-06-01",
    }
    pd.concat([meta, pd.DataFrame([witness])], ignore_index=True).to_csv(
        sheet, sep="\t", index=False
    )
    if not with_cache:
        return
    folder = os.path.join(run, "allele_freq_cache", "mid")
    os.makedirs(folder)
    path = os.path.join(folder, "M1_mid_allele_frequency.parquet")
    if empty_cache:
        pq.write_table(empty_cache_table("longitudinal"), path)
        return
    template = pd.read_parquet(
        os.path.join(run, "allele_freq_cache", "pre", "M1_pre_allele_frequency.parquet")
    )
    rows = template[template["sample_id"].astype(str) == "t1_pre"].copy()
    for col, value in (
        ("sample_id", "w1_mid"),
        ("group", "oth"),
        ("subjectID", "w1"),
        ("replicate", "W"),
        ("time", "mid"),
    ):
        rows[col] = pd.Categorical([value] * len(rows))
    at_150 = rows["position"] == 150
    rows.loc[at_150, ["A", "T"]] = [20, 20]
    rows.loc[at_150, ["A_frequency", "T_frequency"]] = [0.5, 0.5]
    rows.to_parquet(path, index=False)


def _origin_at(out, position):
    """The verdict of the candidate at ``position`` in the per-MAG events table."""
    events = pd.read_csv(os.path.join(out, "M1_gene_de_novo_events.tsv"), sep="\t")
    return events.loc[events["position"] == position, "origin"].item()


def _tables(out):
    """Every table the command wrote, as text keyed by file name."""
    return {
        name: pd.read_csv(os.path.join(out, name), sep="\t", dtype=str).to_csv(
            index=False
        )
        for name in sorted(os.listdir(out))
        if name.endswith(".tsv")
    }


def _run_cli(run, out, *extra):
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "alleleflux.scripts.analysis.gene_de_novo",
            "--comparison",
            "pre_post-trt_ctl",
            "--cache_dir",
            os.path.join(run, "allele_freq_cache"),
            "--allele_analysis_dir",
            os.path.join(run, "allele_analysis"),
            "--metadata",
            os.path.join(run, "metadata.tsv"),
            "--prodigal_fasta",
            os.path.join(run, "ref", "genes.fna"),
            "--fasta",
            os.path.join(run, "ref", "mags.fasta"),
            "--mag_mapping",
            os.path.join(run, "ref", "mag_mapping.tsv"),
            "--output_dir",
            out,
            "--min_cages",
            "2",
            "--cpus",
            "1",
            *extra,
        ],
        capture_output=True,
        text=True,
    )


class TestGeneDeNovoCLI:
    def test_planted_events_reach_the_tests_table(self, tmp_path):
        run = _write_run(tmp_path / "run")
        out = str(tmp_path / "out")
        done = _run_cli(run, out, "--strain_gate", "none")
        assert done.returncode == 0, done.stderr
        tests = pd.read_csv(
            os.path.join(out, "gene_de_novo_tests_pre_post-trt_ctl.tsv"), sep="\t"
        )
        assert list(tests.columns) == result_columns(("trt", "ctl"))
        assert len(tests) == 1
        row = tests.iloc[0]
        assert (row.mag_id, row.gene_id, row.gene_length) == ("M1", "M1_c1_1", 300)
        # two treatment cages, one event each in 300 readable bases; the whole
        # genome read: 1,000 bases per cage, one event per treatment cage
        assert (
            row.n_cages_trt,
            row.n_de_novo_trt,
            row.n_cages_ctl,
            row.n_de_novo_ctl,
        ) == (
            2,
            2,
            2,
            0,
        )
        assert row.expected_trt == pytest.approx(2 * 300 / 1000)
        # the control cages caught nothing anywhere in the genome, so the
        # gene's fair share there is 0: no expectation to test against
        assert row.expected_ctl == 0 and row.note_ctl == "no_exposure"
        assert np.isfinite(row.p_value_between) and 0 < row.p_value_between < 1
        assert row.df_resid_between == 2 and np.isinf(row.rate_ratio_trt_vs_ctl)
        # the three per-MAG files, with the pinned columns
        events = pd.read_csv(os.path.join(out, "M1_gene_de_novo_events.tsv"), sep="\t")
        assert list(events.columns) == EVENT_COLUMNS
        assert sorted(events["origin"]) == ["de_novo", "de_novo", "standing_variation"]
        genes = pd.read_csv(os.path.join(out, "M1_gene_de_novo_genes.tsv"), sep="\t")
        assert list(genes.columns) == GENE_CAGE_COLUMNS and len(genes) == 4
        assert (genes["covered_bp"] == 300).all()
        cages = pd.read_csv(os.path.join(out, "M1_gene_de_novo_cages.tsv"), sep="\t")
        assert list(cages.columns) == CAGE_TOTAL_COLUMNS
        assert cages.set_index("replicate")["n_de_novo_total"].to_dict() == {
            "A": 1,
            "B": 1,
            "C": 0,
            "D": 0,
        }
        assert (cages["covered_bp_total"] == GENOME).all()

    # ---- Phase 3b: witnesses outside the comparison, and the two lists that
    # say which of them the workflow built caches for.
    def test_cohort_witness_at_a_third_timepoint_makes_the_allele_standing(
        self, tmp_path
    ):
        """The reference case for the tests below: w1_mid (another group,
        another timepoint, collected first) already carries the T that rises
        in cage A, so that rise is standing variation, not an event."""
        run = _write_run(tmp_path / "run")
        _add_witness_sample(run)
        out = str(tmp_path / "out")
        done = _run_cli(run, out, "--strain_gate", "none")
        assert done.returncode == 0, done.stderr
        assert _origin_at(out, 150) == "standing_variation"
        assert _origin_at(out, 250) == "de_novo"  # cage B's T at 250: no witness

    # ---- three labels: --pre_treatment_timepoints
    def test_the_label_is_reported_and_no_count_or_test_moves(self, tmp_path):
        """w1_mid carries the T that rises in cage A at 150.  With 'pre' as
        the pre-treatment timepoint that witness was sampled after treatment
        began, so the rise reads first_seen_after_treatment_began in the
        events table.  It was not an event before (standing) and is not one
        now: the gene, cage and test tables are identical text."""
        run = _write_run(tmp_path / "run")
        _add_witness_sample(run)
        plain_out, split_out = str(tmp_path / "plain"), str(tmp_path / "split")
        assert _run_cli(run, plain_out, "--strain_gate", "none").returncode == 0
        done = _run_cli(
            run, split_out, "--strain_gate", "none", "--pre_treatment_timepoints", "pre"
        )
        assert done.returncode == 0, done.stderr
        assert _origin_at(plain_out, 150) == "standing_variation"
        assert _origin_at(split_out, 150) == "first_seen_after_treatment_began"
        assert _origin_at(split_out, 250) == "de_novo"
        assert _origin_at(split_out, 600) == "standing_variation"  # c1's own pre sample
        plain, split = _tables(plain_out), _tables(split_out)
        assert sorted(plain) == sorted(split)
        events_file = "M1_gene_de_novo_events.tsv"
        for name in plain:
            if name != events_file:
                assert plain[name] == split[name], name
        # the events table differs in that one cell only
        a = pd.read_csv(os.path.join(plain_out, events_file), sep="\t", dtype=str)
        b = pd.read_csv(os.path.join(split_out, events_file), sep="\t", dtype=str)
        assert list(a.columns) == list(b.columns) == EVENT_COLUMNS
        differs = (a.fillna("") != b.fillna("")).any(axis=1)
        assert a.loc[differs, "position"].tolist() == ["150"]
        pd.testing.assert_frame_equal(
            a.drop(columns="origin"), b.drop(columns="origin")
        )

    def test_a_witness_at_a_pre_treatment_timepoint_keeps_it_standing(self, tmp_path):
        run = _write_run(tmp_path / "run")
        _add_witness_sample(run)
        plain_out, split_out = str(tmp_path / "plain"), str(tmp_path / "split")
        assert _run_cli(run, plain_out, "--strain_gate", "none").returncode == 0
        done = _run_cli(
            run,
            split_out,
            "--strain_gate",
            "none",
            "--pre_treatment_timepoints",
            "pre",
            "mid",
        )
        assert done.returncode == 0, done.stderr
        assert _tables(split_out) == _tables(plain_out)

    def test_pre_treatment_flag_under_another_scope_is_refused(self, tmp_path):
        run = _write_run(tmp_path / "run")
        done = _run_cli(
            run,
            str(tmp_path / "out"),
            "--strain_gate",
            "none",
            "--baseline_scope",
            "comparison",
            "--pre_treatment_timepoints",
            "pre",
        )
        assert done.returncode != 0
        assert "pre_treatment_timepoints" in done.stderr

    def test_cohort_missing_cache_is_an_error_naming_the_file(self, tmp_path):
        run = _write_run(tmp_path / "run")
        _add_witness_sample(run, with_cache=False)
        done = _run_cli(run, str(tmp_path / "out"), "--strain_gate", "none")
        assert done.returncode != 0
        assert "FileNotFoundError" in done.stderr
        assert "M1_mid_allele_frequency.parquet" in done.stderr

    def test_cohort_zero_row_cache_at_a_third_timepoint_changes_nothing(self, tmp_path):
        """The genome has no QC-passing sample at ``mid``: the workflow wrote a
        zero-row cache there.  Every table equals the run on the sheet without
        that sample."""
        plain_run = _write_run(tmp_path / "plain")
        plain_out = str(tmp_path / "plain_out")
        assert _run_cli(plain_run, plain_out, "--strain_gate", "none").returncode == 0
        run = _write_run(tmp_path / "run")
        _add_witness_sample(run, empty_cache=True)
        out = str(tmp_path / "out")
        done = _run_cli(run, out, "--strain_gate", "none")
        assert done.returncode == 0, done.stderr
        assert _tables(out) == _tables(plain_out)

    def test_witness_lists_covering_the_sheet_change_nothing(self, tmp_path):
        run = _write_run(tmp_path / "run")
        _add_witness_sample(run)
        plain_out, listed_out = str(tmp_path / "plain"), str(tmp_path / "listed")
        assert _run_cli(run, plain_out, "--strain_gate", "none").returncode == 0
        done = _run_cli(
            run,
            listed_out,
            "--strain_gate",
            "none",
            "--witness_timepoints",
            "pre",
            "post",
            "mid",
            "--witness_groups",
            "trt",
            "ctl",
            "oth",
        )
        assert done.returncode == 0, done.stderr
        assert _tables(listed_out) == _tables(plain_out)
        assert "not read as baseline witnesses" not in done.stderr

    def test_a_sample_outside_the_witness_timepoints_does_not_testify(self, tmp_path):
        """``mid`` is left out of the witness timepoints: w1_mid is not read (its
        cache is not even there), so cage A's T at 150 is an event again and
        the run says one sample was left out."""
        run = _write_run(tmp_path / "run")
        _add_witness_sample(run, with_cache=False)
        out = str(tmp_path / "out")
        done = _run_cli(
            run, out, "--strain_gate", "none", "--witness_timepoints", "pre", "post"
        )
        assert done.returncode == 0, done.stderr
        assert _origin_at(out, 150) == "de_novo"
        assert "1 sample(s) of the sheet not read as baseline witnesses" in done.stderr

    def test_a_sample_outside_the_witness_groups_does_not_testify(self, tmp_path):
        run = _write_run(tmp_path / "run")
        _add_witness_sample(run)
        out = str(tmp_path / "out")
        done = _run_cli(
            run,
            out,
            "--strain_gate",
            "none",
            "--witness_timepoints",
            "pre",
            "post",
            "mid",
            "--witness_groups",
            "trt",
            "ctl",
        )
        assert done.returncode == 0, done.stderr
        assert _origin_at(out, 150) == "de_novo"
        assert "1 sample(s) of the sheet not read as baseline witnesses" in done.stderr

    def test_the_comparisons_own_samples_are_never_filtered_out(self, tmp_path):
        """Lists naming only the witness's timepoint and group: the comparison's
        own samples are read all the same, so the tables equal the plain run."""
        run = _write_run(tmp_path / "run")
        _add_witness_sample(run)
        plain_out, listed_out = str(tmp_path / "plain"), str(tmp_path / "listed")
        assert _run_cli(run, plain_out, "--strain_gate", "none").returncode == 0
        done = _run_cli(
            run,
            listed_out,
            "--strain_gate",
            "none",
            "--witness_timepoints",
            "mid",
            "--witness_groups",
            "oth",
        )
        assert done.returncode == 0, done.stderr
        assert _tables(listed_out) == _tables(plain_out)

    def test_carrier_gate_reads_the_turnover_table(self, tmp_path):
        run = _write_run(tmp_path / "run")
        out = str(tmp_path / "out")
        done = _run_cli(
            run,
            out,
            "--strain_gate",
            "carrier",
            "--turnover_dir",
            os.path.join(run, "strain_turnover"),
            "--min_determined_carriers",
            "1",
        )
        assert done.returncode == 0, done.stderr
        events = pd.read_csv(os.path.join(out, "M1_gene_de_novo_events.tsv"), sep="\t")
        by_cage = events.set_index("replicate")["origin"].to_dict()
        # t1 (cage A) changed strain -> its event is a new strain, not a mutation
        assert (by_cage["A"], by_cage["B"]) == ("new_strain", "de_novo")
        tests = pd.read_csv(
            os.path.join(out, "gene_de_novo_tests_pre_post-trt_ctl.tsv"), sep="\t"
        )
        assert tests.iloc[0].n_de_novo_trt == 1

    def test_group_gate_is_refused_at_the_command_line(self, tmp_path):
        run = _write_run(tmp_path / "run")
        done = _run_cli(run, str(tmp_path / "out"), "--strain_gate", "group")
        assert done.returncode != 0
        assert "invalid choice" in done.stderr

    def test_carrier_gate_needs_the_turnover_dir(self, tmp_path):
        run = _write_run(tmp_path / "run")
        done = _run_cli(run, str(tmp_path / "out"), "--strain_gate", "carrier")
        assert done.returncode != 0
        assert "--turnover_dir" in done.stderr
