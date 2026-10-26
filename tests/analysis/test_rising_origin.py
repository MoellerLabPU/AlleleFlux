"""Tests for alleleflux-rising-origin (alleleflux/scripts/analysis/rising_origin.py).

The fixture is small enough to check by hand.  One comparison
(``5mo_16mo-40_AL``: group 40 against AL), within-group test:

    rising alleles                              replicate verdict (floor 8 | floor 4)
    MAG_A G  rose in 40  standing_variation     MAG_A 40  not_replaced   | not_replaced
    MAG_A T  rose in 40  de_novo_candidate      MAG_B 40  replaced       | not_replaced
    MAG_B C  rose in 40  de novo, below det.    MAG_C 40  too_few_voters | too_few_voters
    MAG_B A  rose in 40  not present at 16mo    MAG_A AL  replaced       | replaced
    MAG_C G  rose in 40  5mo not covered        MAG_D     (no row: never classified)
    MAG_A C  rose in AL  de_novo_candidate
    MAG_D T  rose in 40  de_novo_candidate

Expected under --strain_gate group (the across-cage vote):
    floor 8: MAG_B C (de novo, MAG_B replaced in 40)  -> new_strain
             MAG_A C in AL (MAG_A replaced in AL)      -> new_strain
             MAG_B A stays not_present_later (only de novo can be relabelled)
             MAG_D T -> strain_undetermined (no row in the table: unverified)
    floor 4: MAG_B is kept in 40, so MAG_B C stays de_novo_below_detection

Carriers (``16mo_mice_allele_present``) and per-mouse turnover verdicts
(``{mag}_strain_turnover.tsv``, dominant_strain_change for 5mo_16mo):
    MAG_A T  carriers m1,m2,m3   m1 True, m2 True, m3 False  -> 3 determined, 2 changed
    MAG_A C  carriers m7,m8,m9   m7 True, m8 False, m9 False -> 3 determined, 1 changed
    MAG_B C  carriers m1,m2,m3   m1 blank, m2 True, m3 True  -> 2 determined (< floor 3)
    MAG_D T  carriers m1,m2      no turnover table for MAG_D
Expected under --strain_gate carrier (floor 3):
    MAG_A T -> new_strain; MAG_A C stays de_novo; MAG_B C and MAG_D T -> strain_undetermined
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

from alleleflux.scripts.analysis.rising_origin import (
    CALLED,
    ORIGINS,
    apply_strain_gate,
    attach_carrier_gate,
    attach_strain,
    combine_arms,
    load_classification,
    load_rising,
    normalise_origin,
    parse_rising_name,
    summarise,
)

COMP = "5mo_16mo-40_AL"
WITHIN = f"{COMP}_single_sample_tTest_baseline_presence_rising.tsv"
BETWEEN = f"{COMP}_two_sample_unpaired_tTest_baseline_presence_rising.tsv"


def _rising_rows(within=True):
    """(mag, allele, rose_in, origin_any_mouse, carriers) -> rising-table rows."""
    rows = [
        ("MAG_A", "G", "40", "standing_variation", "m1,m2"),
        ("MAG_A", "T", "40", "de_novo_candidate", "m1,m2,m3"),
        ("MAG_B", "C", "40", "de_novo_candidate_below_detection_at_5mo", "m1,m2,m3"),
        ("MAG_B", "A", "40", "allele_not_present_at_16mo", ""),
        ("MAG_C", "G", "40", "5mo_not_covered", "m1"),
        ("MAG_A", "C", "AL", "de_novo_candidate", "m7,m8,m9"),
        ("MAG_D", "T", "40", "de_novo_candidate", "m1,m2"),
    ]
    return pd.DataFrame(
        [
            {
                "mag_id": mag,
                "contig": f"{mag}_c1",
                "position": i,
                "gene_id": f"{mag}_g1",
                "group_analyzed": rose if within else "",
                "allele": allele,
                "rose_in": rose,
                "frequency_change": 0.5,
                "n_replicates": 3,
                "n_mice": 4,
                "origin_any_mouse": origin,
                "16mo_mice_allele_present": carriers,
            }
            for i, (mag, allele, rose, origin, carriers) in enumerate(rows)
        ]
    )


def _write_turnover(directory):
    """Per-MAG strain-turnover tables (the columns the carrier gate reads) for
    MAG_A, MAG_B and MAG_C; none for MAG_D.  Verdicts for transition 5mo_16mo."""
    verdicts = {
        "MAG_A": [
            ("m1", True),
            ("m2", True),
            ("m3", False),
            ("m7", True),
            ("m8", False),
            ("m9", False),
        ],
        "MAG_B": [("m1", None), ("m2", True), ("m3", True)],
        "MAG_C": [("m1", True)],
    }
    os.makedirs(directory, exist_ok=True)
    for mag, rows in verdicts.items():
        pd.DataFrame(
            [
                {
                    "MAG_ID": mag,
                    "subjectID": m,
                    "transition": "5mo_16mo",
                    "sample_t1": f"{m}_5mo",
                    "sample_t2": f"{m}_16mo",
                    "conANI": 0.998 if v else 0.9999,
                    "dominant_strain_change": v,
                }
                for m, v in rows
            ]
        ).to_csv(
            os.path.join(directory, f"{mag}_strain_turnover.tsv"), sep="\t", index=False
        )


def _classification(floor):
    """A replacement_classification table for one floor (the columns read here)."""
    b40 = "replaced" if floor == 8 else "not_replaced"
    verdicts = [
        ("MAG_A", "40", "not_replaced"),
        ("MAG_B", "40", b40),
        ("MAG_C", "40", "too_few_voters"),
        ("MAG_A", "AL", "replaced"),
    ]
    rows = [
        {
            "MAG_ID": mag,
            "group": grp,
            "transition": "5mo_16mo",
            "metric": "dominant_strain_change",
            "min_voters": floor,
            "strain_status_replicates": status,
        }
        for mag, grp, status in verdicts
    ]
    # another metric's row must be ignored
    rows.append(
        {
            **rows[0],
            "metric": "strain_replacement",
            "strain_status_replicates": "replaced",
        }
    )
    return pd.DataFrame(rows)


def _tables():
    return [
        load_classification(
            _classification(f), "dominant_strain_change", "strain_status_replicates"
        )
        for f in (8, 4)
    ]


class TestParseRisingName(unittest.TestCase):
    def test_within_and_between(self):
        self.assertEqual(parse_rising_name(WITHIN), (COMP, "single_sample_tTest"))
        self.assertEqual(
            parse_rising_name("/x/" + BETWEEN), (COMP, "two_sample_unpaired_tTest")
        )

    def test_not_a_rising_file_raises(self):
        with self.assertRaises(ValueError):
            parse_rising_name(
                f"{COMP}_single_sample_tTest_baseline_presence_summary.tsv"
            )


class TestNormaliseOrigin(unittest.TestCase):
    def test_every_label(self):
        got = [
            normalise_origin(o, "5mo", "16mo")
            for o in (
                "standing_variation",
                "de_novo_candidate",
                "de_novo_candidate_below_detection_at_5mo",
                "allele_not_present_at_16mo",
                "5mo_not_covered",
            )
        ]
        self.assertEqual(
            got,
            [
                "standing_variation",
                "de_novo",
                "de_novo_below_detection",
                "not_present_later",
                "baseline_not_covered",
            ],
        )

    def test_unknown_label_raises(self):
        with self.assertRaises(ValueError):
            normalise_origin("de_novo_candidate_below_detection_at_9mo", "5mo", "16mo")


class TestLoadRising(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_reads_comparison_test_and_ages_from_the_file_name(self):
        _rising_rows().to_csv(os.path.join(self.tmp, WITHIN), sep="\t", index=False)
        got = load_rising(self.tmp)
        self.assertEqual(len(got), 7)
        row = got.iloc[0]
        self.assertEqual(
            (row.comparison, row.test_type, row.period, row.age),
            (COMP, "single_sample_tTest", "5mo_16mo", "16mo"),
        )
        self.assertEqual(
            got.origin.tolist()[:3],
            ["standing_variation", "de_novo", "de_novo_below_detection"],
        )
        self.assertEqual(got.rose_in.tolist()[5], "AL")

    def test_carriers_column_is_read_by_the_later_label(self):
        _rising_rows().to_csv(os.path.join(self.tmp, WITHIN), sep="\t", index=False)
        got = load_rising(self.tmp)
        row = got[(got.mag_id == "MAG_A") & (got.allele == "T")].iloc[0]
        self.assertEqual(row["carriers"], "m1,m2,m3")
        self.assertEqual(
            got.loc[got.mag_id == "MAG_B", "carriers"].tolist(), ["m1,m2,m3", ""]
        )

    def test_table_without_carriers_column_raises(self):
        # the carrier column is required like every other: a table without it
        # was written by an older baseline presence and must be regenerated
        _rising_rows().drop(columns="16mo_mice_allele_present").to_csv(
            os.path.join(self.tmp, WITHIN), sep="\t", index=False
        )
        with self.assertRaisesRegex(ValueError, "16mo_mice_allele_present"):
            load_rising(self.tmp)

    def test_numeric_group_names_stay_text(self):
        rows = _rising_rows().assign(rose_in="40", group_analyzed="40")
        rows.to_csv(os.path.join(self.tmp, WITHIN), sep="\t", index=False)
        got = load_rising(self.tmp)
        self.assertTrue(all(isinstance(g, str) for g in got.rose_in))

    def test_header_only_file_adds_no_rows(self):
        _rising_rows().iloc[0:0].to_csv(
            os.path.join(self.tmp, WITHIN), sep="\t", index=False
        )
        self.assertTrue(load_rising(self.tmp).empty)

    def test_no_rising_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            load_rising(self.tmp)

    def test_missing_column_raises(self):
        _rising_rows().drop(columns="rose_in").to_csv(
            os.path.join(self.tmp, WITHIN), sep="\t", index=False
        )
        with self.assertRaisesRegex(ValueError, "rose_in"):
            load_rising(self.tmp)


class TestLoadClassification(unittest.TestCase):
    def test_floor_and_metric(self):
        floor, v = load_classification(
            _classification(8), "dominant_strain_change", "strain_status_replicates"
        )
        self.assertEqual(floor, 8)
        self.assertEqual(len(v), 4)  # the strain_replacement row is not read
        self.assertEqual(list(v.columns), ["mag_id", "group", "period", "status"])

    def test_two_floors_in_one_table_raises(self):
        t = pd.concat([_classification(8), _classification(4)], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "one voter floor"):
            load_classification(t, "dominant_strain_change", "strain_status_replicates")

    def test_unknown_status_raises(self):
        t = _classification(8)
        t.loc[0, "strain_status_replicates"] = "maybe"
        with self.assertRaisesRegex(ValueError, "maybe"):
            load_classification(t, "dominant_strain_change", "strain_status_replicates")

    def test_missing_metric_raises(self):
        with self.assertRaisesRegex(ValueError, "popANI"):
            load_classification(
                _classification(8), "popANI", "strain_status_replicates"
            )


class TestCombineArms(unittest.TestCase):
    def test_rules(self):
        self.assertEqual(combine_arms("not_replaced", "replaced"), "replaced")
        self.assertEqual(combine_arms("not_replaced", "not_replaced"), "not_replaced")
        self.assertEqual(combine_arms("tie", "not_replaced"), "tie")
        self.assertEqual(combine_arms("not_replaced", "no_voters"), "unverified")
        self.assertEqual(combine_arms("not_in_strain_table", "replaced"), "replaced")

    def test_unknown_status_raises(self):
        with self.assertRaises(ValueError):
            combine_arms("maybe", "replaced")


class TestAttachStrain(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        _rising_rows().to_csv(os.path.join(self.tmp, WITHIN), sep="\t", index=False)
        self.rising = load_rising(self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _row(self, got, floor, mag, allele, rose):
        r = got[
            (got.floor == floor)
            & (got.mag_id == mag)
            & (got.allele == allele)
            & (got.rose_in == rose)
        ]
        self.assertEqual(len(r), 1)
        return r.iloc[0]

    def test_one_copy_per_floor(self):
        got = attach_strain(self.rising, _tables())
        self.assertEqual(sorted(got.floor.unique()), [4, 8])
        self.assertEqual(len(got), 2 * 7)

    def test_attach_strain_reports_and_never_relabels(self):
        # the group verdict is written beside every row; origins are untouched
        # (apply_strain_gate is the only place an origin changes)
        got = attach_strain(self.rising, _tables())
        b = self._row(got, 8, "MAG_B", "C", "40")
        self.assertEqual(
            (b.strain_class, b.origin), ("replaced", "de_novo_below_detection")
        )
        al = self._row(got, 8, "MAG_A", "C", "AL")  # the AL arm has its own verdict
        self.assertEqual((al.strain_class, al.origin), ("replaced", "de_novo"))
        a = self._row(got, 8, "MAG_A", "T", "40")  # same letter, kept MAG
        self.assertEqual((a.strain_class, a.origin), ("not_replaced", "de_novo"))
        self.assertEqual(
            self._row(got, 8, "MAG_B", "A", "40").origin, "not_present_later"
        )

    def test_floors_differ(self):
        got = attach_strain(self.rising, _tables())
        b4 = self._row(got, 4, "MAG_B", "C", "40")
        self.assertEqual(b4.strain_class, "not_replaced")
        self.assertEqual(self._row(got, 8, "MAG_B", "C", "40").strain_class, "replaced")
        self.assertEqual(
            self._row(got, 8, "MAG_C", "G", "40").strain_class, "unverified"
        )

    def test_never_classified_mag_is_unverified_not_guessed(self):
        got = attach_strain(self.rising, _tables())
        d = self._row(got, 8, "MAG_D", "T", "40")
        self.assertEqual(
            (d.strain_status, d.strain_class, d.origin),
            ("not_in_strain_table", "unverified", "de_novo"),
        )

    def test_raw_verdict_is_kept(self):
        got = attach_strain(self.rising, _tables())
        b = self._row(got, 8, "MAG_B", "C", "40")
        self.assertEqual(b.origin_any_mouse, "de_novo_candidate_below_detection_at_5mo")
        self.assertEqual(b.strain_status, "replaced")

    def test_between_group_rows_take_the_dirtier_arm(self):
        # MAG_A: 40 not_replaced, AL replaced -> replaced, whichever group rose
        rising = self.rising.assign(
            group_analyzed="", test_type="two_sample_unpaired_tTest"
        )
        got = attach_strain(rising, _tables())
        a40 = self._row(got, 8, "MAG_A", "T", "40")
        self.assertEqual((a40.strain_class, a40.origin), ("replaced", "de_novo"))
        self.assertEqual(a40.strain_status, "40:not_replaced|AL:replaced")
        # MAG_B at floor 4: 40 not_replaced, AL never classified -> unverified
        b4 = self._row(got, 4, "MAG_B", "C", "40")
        self.assertEqual(b4.strain_class, "unverified")


class TestAttachCarrierGate(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        _rising_rows().to_csv(os.path.join(self.tmp, WITHIN), sep="\t", index=False)
        self.turnover = os.path.join(self.tmp, "strain_turnover")
        _write_turnover(self.turnover)
        self.rising = load_rising(self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _row(self, got, mag, allele):
        return got[(got.mag_id == mag) & (got.allele == allele)].iloc[0]

    def test_four_columns_from_the_carriers_own_verdicts(self):
        got = attach_carrier_gate(self.rising, self.turnover, 3)
        t = self._row(got, "MAG_A", "T")
        self.assertEqual(
            (
                t.carrier_strain_status,
                t.n_carriers,
                t.n_carriers_determined,
                t.n_carriers_strain_changed,
            ),
            ("new_strain_carriers", 3, 3, 2),
        )
        c = self._row(got, "MAG_A", "C")
        self.assertEqual(
            (c.carrier_strain_status, c.n_carriers_strain_changed),
            ("same_strain_carriers", 1),
        )
        b = self._row(
            got, "MAG_B", "C"
        )  # m1 undetermined: 2 of 3 determined, under the floor
        self.assertEqual(
            (b.carrier_strain_status, b.n_carriers_determined),
            ("strain_undetermined", 2),
        )
        a = self._row(got, "MAG_B", "A")  # no carriers at all
        self.assertEqual(
            (a.carrier_strain_status, a.n_carriers), ("strain_undetermined", 0)
        )
        self.assertEqual(
            list(got["origin"]), list(self.rising["origin"])
        )  # no relabel here

    def test_mag_without_turnover_table_is_undetermined_with_a_warning(self):
        with self.assertLogs(
            "alleleflux.scripts.analysis.rising_origin", level="WARNING"
        ) as logs:
            got = attach_carrier_gate(self.rising, self.turnover, 3)
        d = self._row(got, "MAG_D", "T")
        self.assertEqual(
            (d.carrier_strain_status, d.n_carriers, d.n_carriers_determined),
            ("strain_undetermined", 2, 0),
        )
        self.assertTrue(any("MAG_D" in m for m in logs.output))

    def test_floor_of_one(self):
        got = attach_carrier_gate(self.rising, self.turnover, 1)
        self.assertEqual(
            self._row(got, "MAG_B", "C").carrier_strain_status, "new_strain_carriers"
        )


class TestApplyStrainGate(unittest.TestCase):
    """The ONE place origins change.  A small frame: every (origin x group class
    x carrier status) combination that matters."""

    def _frame(self):
        rows = [
            # origin,                  strain_class,   carrier_strain_status
            (
                "de_novo",
                "replaced",
                "same_strain_carriers",
            ),  # 0 group says swap, carriers say kept
            ("de_novo", "not_replaced", "new_strain_carriers"),  # 1 the reverse
            ("de_novo", "not_replaced", "same_strain_carriers"),  # 2 clean both ways
            ("de_novo", "unverified", "strain_undetermined"),  # 3 nobody can say
            ("de_novo", "tie", "same_strain_carriers"),  # 4 tie in the group vote
            ("de_novo_below_detection", "replaced", "strain_undetermined"),  # 5
            (
                "de_novo",
                "not_replaced",
                "strain_undetermined",
            ),  # 6 group clean, carriers unknown
            (
                "standing_variation",
                "replaced",
                "new_strain_carriers",
            ),  # 7 never relabelled
            ("block_event", "replaced", "new_strain_carriers"),  # 8 never relabelled
            (
                "not_present_later",
                "unverified",
                "strain_undetermined",
            ),  # 9 never relabelled
        ]
        return pd.DataFrame(
            rows, columns=["origin", "strain_class", "carrier_strain_status"]
        )

    def test_group(self):
        got = apply_strain_gate(self._frame(), "group")
        self.assertEqual(
            list(got["origin"]),
            [
                "new_strain",
                "de_novo",
                "de_novo",
                "strain_undetermined",
                "strain_undetermined",
                "new_strain",
                "de_novo",
                "standing_variation",
                "block_event",
                "not_present_later",
            ],
        )

    def test_carrier(self):
        got = apply_strain_gate(self._frame(), "carrier")
        self.assertEqual(
            list(got["origin"]),
            [
                "de_novo",
                "new_strain",
                "de_novo",
                "strain_undetermined",
                "de_novo",
                "strain_undetermined",
                "strain_undetermined",
                "standing_variation",
                "block_event",
                "not_present_later",
            ],
        )

    def test_both_needs_both_gates_to_agree(self):
        # user ruling 2026-10-26: under `both` an allele is only EXPLAINED when
        # both gates could speak; one undetermined verdict makes the row
        # undetermined even if the other gate says replaced (row 5)
        got = apply_strain_gate(self._frame(), "both")
        self.assertEqual(
            list(got["origin"]),
            [
                "new_strain",
                "new_strain",
                "de_novo",
                "strain_undetermined",
                "strain_undetermined",
                "strain_undetermined",
                "strain_undetermined",
                "standing_variation",
                "block_event",
                "not_present_later",
            ],
        )

    def test_none_changes_nothing(self):
        f = self._frame()
        self.assertEqual(
            list(apply_strain_gate(f, "none")["origin"]), list(f["origin"])
        )

    def test_missing_columns_raise(self):
        f = self._frame()
        with self.assertRaisesRegex(ValueError, "strain_class"):
            apply_strain_gate(f.drop(columns="strain_class"), "group")
        with self.assertRaisesRegex(ValueError, "carrier_strain_status"):
            apply_strain_gate(f.drop(columns="carrier_strain_status"), "carrier")
        with self.assertRaisesRegex(ValueError, "gate"):
            apply_strain_gate(f, "maybe")

    def test_input_untouched(self):
        f = self._frame()
        before = f.copy()
        apply_strain_gate(f, "both")
        pd.testing.assert_frame_equal(f, before)


class TestSummarise(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        rows = _rising_rows()
        rows = rows[rows.mag_id != "MAG_D"]  # the 6-allele hand count below
        rows.to_csv(os.path.join(self.tmp, WITHIN), sep="\t", index=False)
        self.rising = load_rising(self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _cell(self, summ, floor, rose, cls, origin):
        floor_match = summ.floor.isna() if floor is None else summ.floor == floor
        r = summ[
            floor_match
            & (summ.rose_in == rose)
            & (summ.strain_class == cls)
            & (summ.origin == origin)
        ]
        self.assertEqual(len(r), 1, (floor, rose, cls, origin))
        return r.iloc[0]

    def test_without_strain_every_origin_but_the_gate_outcomes_has_a_row(self):
        # no strain input at all: new_strain / strain_undetermined cannot occur
        summ = summarise(self.rising)
        diet = summ[summ.rose_in == "40"]
        self.assertEqual(
            sorted(diet.origin),
            sorted(
                o for o in ORIGINS if o not in ("new_strain", "strain_undetermined")
            ),
        )
        self.assertTrue((diet.strain_class == "all").all())
        self.assertTrue(diet.floor.isna().all())

    def test_percentages_of_all_and_of_called(self):
        summ = summarise(self.rising)
        # group 40: 5 rising alleles; called = standing + de novo + below det. = 3
        st = self._cell(summ, None, "40", "all", "standing_variation")
        self.assertEqual((int(st.n_alleles), int(st.n_mags)), (1, 1))
        self.assertAlmostEqual(st.pct_of_all, 100 * 1 / 5)
        self.assertAlmostEqual(st.pct_of_called, 100 * 1 / 3)
        nc = self._cell(summ, None, "40", "all", "baseline_not_covered")
        self.assertAlmostEqual(nc.pct_of_all, 100 * 1 / 5)
        self.assertTrue(np.isnan(nc.pct_of_called))

    def test_zero_counts_are_rows_not_gaps(self):
        summ = summarise(self.rising)
        al = self._cell(summ, None, "AL", "all", "standing_variation")
        self.assertEqual((int(al.n_alleles), al.pct_of_all), (0, 0.0))

    def test_with_strain_per_floor_and_class_plus_all(self):
        summ = summarise(
            apply_strain_gate(attach_strain(self.rising, _tables()), "group")
        )
        ns = self._cell(summ, 8, "40", "all", "new_strain")
        self.assertEqual(int(ns.n_alleles), 1)
        self.assertAlmostEqual(ns.pct_of_called, 100 * 1 / 3)
        rep = self._cell(summ, 8, "40", "replaced", "new_strain")
        self.assertAlmostEqual(
            rep.pct_of_all, 50.0
        )  # MAG_B: new_strain + not_present_later
        self.assertAlmostEqual(rep.pct_of_called, 100.0)
        self.assertEqual(
            int(self._cell(summ, 4, "40", "all", "new_strain").n_alleles), 0
        )

    def test_percentages_add_up(self):
        summ = summarise(
            apply_strain_gate(attach_strain(self.rising, _tables()), "group")
        )
        for _, g in summ.groupby(["floor", "comparison", "rose_in", "strain_class"]):
            if g.n_alleles.sum():
                self.assertAlmostEqual(g.pct_of_all.sum(), 100.0)
            if g[g.origin.isin(CALLED)].n_alleles.sum():
                self.assertAlmostEqual(g.pct_of_called.sum(), 100.0)

    def test_comparisons_are_never_pooled(self):
        other = self.rising.assign(
            comparison="5mo_22mo-40_AL", period="5mo_22mo", age="22mo"
        )
        summ = summarise(pd.concat([self.rising, other], ignore_index=True))
        st = summ[(summ.rose_in == "40") & (summ.origin == "standing_variation")]
        self.assertEqual(sorted(st.comparison), ["5mo_16mo-40_AL", "5mo_22mo-40_AL"])
        self.assertTrue((st.n_alleles == 1).all())


class TestRisingOriginCLI(unittest.TestCase):
    """The command end to end: rising tables + two floors' classification tables."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.rising_dir = os.path.join(self.tmp, "baseline_presence")
        os.makedirs(self.rising_dir)
        _rising_rows().to_csv(
            os.path.join(self.rising_dir, WITHIN), sep="\t", index=False
        )
        self.tables = []
        for f in (8, 4):
            path = os.path.join(self.tmp, f"classification_{f}.tsv")
            _classification(f).to_csv(path, sep="\t", index=False)
            self.tables.append(path)
        self.turnover = os.path.join(self.tmp, "strain_turnover")
        _write_turnover(self.turnover)
        self.out = os.path.join(self.tmp, "out")

    def _alleles(self):
        return pd.read_csv(
            os.path.join(self.out, "rising_origin_alleles.tsv"), sep="\t", dtype=str
        )

    def _origin(self, alleles, mag, allele, floor=None):
        r = alleles[(alleles.mag_id == mag) & (alleles.allele == allele)]
        if floor is not None:
            r = r[r.floor == str(floor)]
        self.assertEqual(len(r), 1)
        return r.iloc[0]

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self, *extra):
        return subprocess.run(
            [
                sys.executable,
                "-m",
                "alleleflux.scripts.analysis.rising_origin",
                "--rising_dir",
                self.rising_dir,
                "--output_dir",
                self.out,
                *extra,
            ],
            capture_output=True,
            text=True,
        )

    def test_group_gate_with_classification_tables(self):
        done = self._run("--classification", *self.tables, "--strain_gate", "group")
        self.assertEqual(done.returncode, 0, done.stderr)
        alleles = self._alleles()
        summ = pd.read_csv(
            os.path.join(self.out, "rising_origin_summary.tsv"), sep="\t"
        )
        self.assertEqual(len(alleles), 2 * 7)
        self.assertEqual(sorted(summ.floor.unique()), [4, 8])
        self.assertEqual(
            int((alleles.origin == "new_strain").sum()), 3
        )  # 2 at 8, 1 at 4
        # the ruling: an unclassified genome's de novo allele is undetermined, not de novo
        self.assertEqual(
            int((alleles.origin == "strain_undetermined").sum()), 2
        )  # MAG_D, both floors
        self.assertEqual(
            self._origin(alleles, "MAG_D", "T", 8).origin, "strain_undetermined"
        )
        self.assertNotIn(
            "carrier_strain_status", alleles.columns
        )  # no turnover dir given

    def test_carrier_gate(self):
        done = self._run("--turnover_dir", self.turnover, "--strain_gate", "carrier")
        self.assertEqual(done.returncode, 0, done.stderr)
        alleles = self._alleles()
        self.assertEqual(len(alleles), 7)  # no floors without classification tables
        t = self._origin(alleles, "MAG_A", "T")
        self.assertEqual(
            (t.origin, t.carrier_strain_status, t.n_carriers_strain_changed),
            ("new_strain", "new_strain_carriers", "2"),
        )
        self.assertEqual(
            self._origin(alleles, "MAG_A", "C").origin, "de_novo"
        )  # 1 of 3 changed
        self.assertEqual(
            self._origin(alleles, "MAG_B", "C").origin, "strain_undetermined"
        )
        self.assertEqual(
            self._origin(alleles, "MAG_D", "T").origin, "strain_undetermined"
        )
        self.assertEqual(
            self._origin(alleles, "MAG_A", "G").origin, "standing_variation"
        )
        self.assertIn("MAG_D", done.stderr)  # the missing-table warning
        summ = pd.read_csv(
            os.path.join(self.out, "rising_origin_summary.tsv"), sep="\t", dtype=str
        ).set_index(["rose_in", "origin"])
        self.assertEqual(summ.loc[("40", "new_strain"), "n_alleles"], "1")
        self.assertEqual(summ.loc[("40", "strain_undetermined"), "n_alleles"], "2")

    def test_carrier_gate_is_the_default(self):
        done = self._run("--turnover_dir", self.turnover)
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual(
            self._origin(self._alleles(), "MAG_A", "T").origin, "new_strain"
        )

    def test_both_takes_the_stricter(self):
        done = self._run(
            "--classification",
            *self.tables,
            "--turnover_dir",
            self.turnover,
            "--strain_gate",
            "both",
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        alleles = self._alleles()
        # floor 4: MAG_B kept in 40 (group says keep) but carriers undetermined -> undetermined
        self.assertEqual(
            self._origin(alleles, "MAG_B", "C", 4).origin, "strain_undetermined"
        )
        # floor 8: group says replaced, carriers undetermined -> undetermined (one
        # gate could not speak, so the allele is not explained, only removed)
        self.assertEqual(
            self._origin(alleles, "MAG_B", "C", 8).origin, "strain_undetermined"
        )
        # MAG_A T: group keep, carriers new strain -> new_strain at both floors
        self.assertEqual(self._origin(alleles, "MAG_A", "T", 8).origin, "new_strain")
        # both verdict columns travel with every row
        self.assertEqual(
            self._origin(alleles, "MAG_A", "T", 8).strain_status, "not_replaced"
        )

    def test_gate_needs_its_inputs(self):
        self.assertNotEqual(
            self._run("--strain_gate", "carrier").returncode, 0
        )  # no turnover dir
        self.assertNotEqual(
            self._run("--strain_gate", "group").returncode, 0
        )  # no classification
        self.assertNotEqual(
            self._run(
                "--turnover_dir", self.turnover, "--strain_gate", "both"
            ).returncode,
            0,
        )
        self.assertEqual(self._run("--strain_gate", "none").returncode, 0)

    def test_cluster_rule_makes_block_events(self):
        """Two de novo alleles of MAG_A in group 40, 500 bp apart (the fixture's T at
        position 1 plus an added A at 501): one block, both block_event, counted in
        the summary as block_event and no longer as de novo.  With one tolerated
        neighbour they stay de novo."""
        extra = _rising_rows().iloc[[1]].assign(allele="A", position=501)
        pd.concat([_rising_rows(), extra]).to_csv(
            os.path.join(self.rising_dir, WITHIN), sep="\t", index=False
        )
        done = self._run("--strain_gate", "none")
        self.assertEqual(done.returncode, 0, done.stderr)
        alleles = pd.read_csv(
            os.path.join(self.out, "rising_origin_alleles.tsv"), sep="\t", dtype=str
        )
        pair = alleles[(alleles.mag_id == "MAG_A") & (alleles.rose_in == "40")]
        pair = pair[pair.allele.isin(["T", "A"])].set_index("allele")
        self.assertEqual(set(pair["origin"]), {"block_event"})
        self.assertEqual(pair.loc["T", "block_id"], pair.loc["A", "block_id"])
        self.assertTrue(pair.loc["T", "block_id"].startswith("MAG_A:MAG_A_c1:40:"))
        self.assertEqual(list(pair["n_de_novo_neighbours"]), ["1", "1"])
        self.assertEqual(
            alleles.loc[alleles.mag_id == "MAG_D", "origin"].iloc[0], "de_novo"
        )
        summ = pd.read_csv(
            os.path.join(self.out, "rising_origin_summary.tsv"), sep="\t", dtype=str
        ).set_index(["rose_in", "origin"])
        self.assertEqual(summ.loc[("40", "block_event"), "n_alleles"], "2")
        self.assertEqual(summ.loc[("40", "de_novo"), "n_alleles"], "1")  # MAG_D only
        self.assertTrue(pd.isna(summ.loc[("40", "block_event"), "pct_of_called"]))
        self.assertIn("1 blocks", done.stderr)
        done = self._run("--strain_gate", "none", "--max_de_novo_neighbours", "1")
        self.assertEqual(done.returncode, 0, done.stderr)
        alleles = pd.read_csv(
            os.path.join(self.out, "rising_origin_alleles.tsv"), sep="\t", dtype=str
        )
        self.assertNotIn("block_event", set(alleles["origin"]))

    def test_without_classification(self):
        done = self._run("--strain_gate", "none")
        self.assertEqual(done.returncode, 0, done.stderr)
        summ = pd.read_csv(
            os.path.join(self.out, "rising_origin_summary.tsv"), sep="\t"
        )
        self.assertTrue(summ.floor.isna().all())
        self.assertNotIn("new_strain", set(summ.origin))


if __name__ == "__main__":
    unittest.main()
