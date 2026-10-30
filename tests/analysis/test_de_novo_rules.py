"""Pure de novo rules shared by baseline presence and rising origin (de_novo_rules.py).

No I/O here: every function takes and returns DataFrames or numbers, so these
tests and both CLIs exercise exactly the same code.
"""

import numpy as np
import pandas as pd
import pytest

from alleleflux.scripts.analysis.de_novo_rules import (
    BLOCK_EVENT,
    assign_blocks,
    carrier_strain_verdicts,
)


def _alleles(rows, allele="A"):
    """rows: (mag, contig, position, rose_in, origin) -> the columns assign_blocks reads."""
    return pd.DataFrame(
        rows, columns=["mag_id", "contig", "position", "rose_in", "origin"]
    ).assign(allele=allele)


class TestAssignBlocks:
    def test_lone_allele_keeps_its_origin(self):
        out = assign_blocks(_alleles([("M", "c1", 100, "fat", "de_novo")]), 1000, 0)
        assert out.loc[0, "origin"] == "de_novo"
        assert out.loc[0, "block_id"] == ""
        assert out.loc[0, "n_de_novo_neighbours"] == 0

    def test_two_within_window_become_one_block(self):
        out = assign_blocks(
            _alleles(
                [
                    ("M", "c1", 100, "fat", "de_novo"),
                    ("M", "c1", 900, "fat", "de_novo_below_detection"),
                ]
            ),
            1000,
            0,
        )
        assert list(out["origin"]) == [BLOCK_EVENT, BLOCK_EVENT]
        assert out["block_id"].nunique() == 1
        assert out.loc[0, "block_id"].startswith("M:c1:fat:")
        assert list(out["n_de_novo_neighbours"]) == [1, 1]

    def test_single_linkage_chains(self):
        # 0 and 1600 are outside the window, but 800 links them: one block
        out = assign_blocks(
            _alleles([("M", "c1", p, "fat", "de_novo") for p in (0, 800, 1600)]),
            1000,
            0,
        )
        assert out["block_id"].nunique() == 1
        assert list(out["n_de_novo_neighbours"]) == [1, 2, 1]

    def test_window_is_inclusive(self):
        out = assign_blocks(
            _alleles(
                [("M", "c1", 0, "fat", "de_novo"), ("M", "c1", 1000, "fat", "de_novo")]
            ),
            1000,
            0,
        )
        assert out["block_id"].nunique() == 1
        out = assign_blocks(
            _alleles(
                [("M", "c1", 0, "fat", "de_novo"), ("M", "c1", 1001, "fat", "de_novo")]
            ),
            1000,
            0,
        )
        assert list(out["origin"]) == ["de_novo", "de_novo"]

    def test_other_group_contig_or_mag_does_not_link(self):
        out = assign_blocks(
            _alleles(
                [
                    ("M", "c1", 100, "fat", "de_novo"),
                    ("M", "c1", 200, "control", "de_novo"),
                    ("M", "c2", 150, "fat", "de_novo"),
                    ("N", "c1", 150, "fat", "de_novo"),
                ]
            ),
            1000,
            0,
        )
        assert list(out["origin"]) == ["de_novo"] * 4

    def test_standing_and_uncalled_alleles_are_not_neighbours(self):
        out = assign_blocks(
            _alleles(
                [
                    ("M", "c1", 100, "fat", "de_novo"),
                    ("M", "c1", 200, "fat", "standing_variation"),
                    ("M", "c1", 300, "fat", "not_present_later"),
                ]
            ),
            1000,
            0,
        )
        assert list(out["origin"]) == [
            "de_novo",
            "standing_variation",
            "not_present_later",
        ]
        assert list(out["n_de_novo_neighbours"]) == [0, 0, 0]

    def test_max_neighbours_tolerates_small_clusters(self):
        rows = [("M", "c1", p, "fat", "de_novo") for p in (100, 200, 300)]
        assert list(assign_blocks(_alleles(rows), 1000, 2)["origin"]) == ["de_novo"] * 3
        assert (
            list(assign_blocks(_alleles(rows), 1000, 1)["origin"]) == [BLOCK_EVENT] * 3
        )

    def test_same_allele_across_comparisons_counts_once(self):
        # the same (mag, contig, position, allele) rising in two comparisons is ONE
        # allele: it is never its own neighbour, and both rows get the same labels
        rows = _alleles(
            [("M", "c1", 100, "fat", "de_novo"), ("M", "c1", 100, "fat", "de_novo")]
        )
        out = assign_blocks(rows, 1000, 0)
        assert list(out["origin"]) == ["de_novo", "de_novo"]
        rows = _alleles(
            [
                ("M", "c1", 100, "fat", "de_novo"),
                ("M", "c1", 100, "fat", "de_novo"),
                ("M", "c1", 600, "fat", "de_novo"),
            ]
        )
        out = assign_blocks(rows, 1000, 0)
        assert list(out["origin"]) == [BLOCK_EVENT] * 3
        assert list(out["n_de_novo_neighbours"]) == [1, 1, 1]
        assert out["block_id"].nunique() == 1

    def test_two_alleles_at_one_position_are_two_nodes(self):
        # A and T at the same position: distinct alleles, 0 bp apart -> neighbours
        rows = pd.concat(
            [
                _alleles([("M", "c1", 100, "fat", "de_novo")], allele="A"),
                _alleles([("M", "c1", 100, "fat", "de_novo")], allele="T"),
            ],
            ignore_index=True,
        )
        out = assign_blocks(rows, 1000, 0)
        assert list(out["origin"]) == [BLOCK_EVENT] * 2

    def test_blocks_are_numbered_by_position_per_group(self):
        rows = [("M", "c1", p, "fat", "de_novo") for p in (5000, 5300, 100, 900)]
        out = assign_blocks(_alleles(rows), 1000, 0)
        assert list(out["block_id"]) == [
            "M:c1:fat:2",
            "M:c1:fat:2",
            "M:c1:fat:1",
            "M:c1:fat:1",
        ]

    def test_input_is_not_modified_and_order_is_kept(self):
        rows = _alleles(
            [("M", "c1", 900, "fat", "de_novo"), ("M", "c1", 100, "fat", "de_novo")]
        )
        before = rows.copy()
        out = assign_blocks(rows, 1000, 0)
        pd.testing.assert_frame_equal(rows, before)
        assert list(out["position"]) == [900, 100]

    def test_empty_table(self):
        out = assign_blocks(_alleles([]), 1000, 0)
        assert (
            list(out.columns)[-2:] == ["block_id", "n_de_novo_neighbours"] and out.empty
        )


def _turnover(rows):
    """rows: (subject, transition, sample_t2, dominant_strain_change) -> the slice of a
    {mag}_strain_turnover.tsv the carrier gate reads.  The verdict is True / False /
    None (None = undetermined, written blank on disk, read back as NaN)."""
    return pd.DataFrame(
        [
            {
                "subjectID": s,
                "transition": t,
                "sample_t2": s2,
                "dominant_strain_change": v,
            }
            for s, t, s2, v in rows
        ],
        columns=["subjectID", "transition", "sample_t2", "dominant_strain_change"],
    )


class TestCarrierStrainVerdicts:
    def test_majority_changed_is_new_strain(self):
        tbl = _turnover(
            [
                ("m1", "pre_end", "m1_end", True),
                ("m2", "pre_end", "m2_end", True),
                ("m3", "pre_end", "m3_end", False),
            ]
        )
        assert carrier_strain_verdicts(["m1", "m2", "m3"], tbl, "pre_end", 3) == (
            "new_strain_carriers",
            3,
            2,
        )

    def test_half_is_not_more_than_half(self):
        tbl = _turnover(
            [
                ("m1", "pre_end", "m1_end", True),
                ("m2", "pre_end", "m2_end", False),
                ("m3", "pre_end", "m3_end", False),
                ("m4", "pre_end", "m4_end", True),
            ]
        )
        assert carrier_strain_verdicts(["m1", "m2", "m3", "m4"], tbl, "pre_end", 3) == (
            "same_strain_carriers",
            4,
            2,
        )

    def test_too_few_determined(self):
        # m2 undetermined (blank verdict), m3 has no row at all: 1 of 3 determined
        tbl = _turnover(
            [("m1", "pre_end", "m1_end", True), ("m2", "pre_end", "m2_end", None)]
        )
        assert carrier_strain_verdicts(["m1", "m2", "m3"], tbl, "pre_end", 3) == (
            "strain_undetermined",
            1,
            1,
        )

    def test_floor_of_one(self):
        tbl = _turnover([("m1", "pre_end", "m1_end", True)])
        assert carrier_strain_verdicts(["m1"], tbl, "pre_end", 1) == (
            "new_strain_carriers",
            1,
            1,
        )

    def test_text_verdicts_as_read_from_disk(self):
        # a TSV read with dtype=str gives "True" / "False" / NaN
        tbl = _turnover(
            [
                ("m1", "pre_end", "m1_end", "True"),
                ("m2", "pre_end", "m2_end", "False"),
                ("m3", "pre_end", "m3_end", np.nan),
            ]
        )
        assert carrier_strain_verdicts(["m1", "m2", "m3"], tbl, "pre_end", 2) == (
            "same_strain_carriers",
            2,
            1,
        )

    def test_mouse_with_two_rows_for_one_transition_raises(self):
        """A mouse sampled twice at the later timepoint would have two verdicts
        for one transition.  That never happens in the designs this is built for
        (user decision 2026-10-26), so it is a data error: stop and name the mouse
        rather than silently pick one of its rows."""
        tbl = _turnover(
            [
                ("m1", "pre_end", "m1_end_b", False),
                ("m1", "pre_end", "m1_end_a", True),
                ("m2", "pre_end", "m2_end", True),
            ]
        )
        with pytest.raises(ValueError, match="m1"):
            carrier_strain_verdicts(["m1", "m2"], tbl, "pre_end", 1)
        # a duplicate in ANOTHER transition does not matter for this one
        tbl = _turnover(
            [
                ("m1", "mid_end", "m1_end_b", False),
                ("m1", "mid_end", "m1_end_a", True),
                ("m1", "pre_end", "m1_end", True),
            ]
        )
        assert carrier_strain_verdicts(["m1"], tbl, "pre_end", 1) == (
            "new_strain_carriers",
            1,
            1,
        )

    def test_other_transitions_or_mice_are_ignored(self):
        tbl = _turnover(
            [("m1", "mid_end", "m1_end", True), ("m9", "pre_end", "m9_end", True)]
        )
        assert carrier_strain_verdicts(["m1"], tbl, "pre_end", 1) == (
            "strain_undetermined",
            0,
            0,
        )

    def test_duplicate_carrier_name_raises(self):
        # a carrier list is a set by construction; a repeated name means a
        # hand-edited or concatenated rising table: stop, do not count twice
        tbl = _turnover([("m1", "pre_end", "m1_end", True)])
        with pytest.raises(ValueError, match="m1"):
            carrier_strain_verdicts(["m1", "m1"], tbl, "pre_end", 1)

    def test_no_carriers_or_empty_table(self):
        assert carrier_strain_verdicts([], _turnover([]), "pre_end", 3) == (
            "strain_undetermined",
            0,
            0,
        )
        assert carrier_strain_verdicts(["m1"], _turnover([]), "pre_end", 1) == (
            "strain_undetermined",
            0,
            0,
        )

    def test_unknown_verdict_word_raises(self):
        tbl = _turnover([("m1", "pre_end", "m1_end", "maybe")])
        with pytest.raises(ValueError):
            carrier_strain_verdicts(["m1"], tbl, "pre_end", 1)
