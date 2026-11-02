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
