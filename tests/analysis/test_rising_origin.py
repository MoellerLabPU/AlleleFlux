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

Expected with the strain tables:
    floor 8: MAG_B C (de novo, MAG_B replaced in 40)  -> new_strain
             MAG_A C in AL (MAG_A replaced in AL)      -> new_strain
             MAG_B A stays not_present_later (only de novo can become new_strain)
             MAG_D T stays de_novo, class unverified (no row in the table)
    floor 4: MAG_B is kept in 40, so MAG_B C stays de_novo_below_detection
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
    """(mag, allele, rose_in, origin_any_mouse) -> rising-table rows."""
    rows = [
        ("MAG_A", "G", "40", "standing_variation"),
        ("MAG_A", "T", "40", "de_novo_candidate"),
        ("MAG_B", "C", "40", "de_novo_candidate_below_detection_at_5mo"),
        ("MAG_B", "A", "40", "allele_not_present_at_16mo"),
        ("MAG_C", "G", "40", "5mo_not_covered"),
        ("MAG_A", "C", "AL", "de_novo_candidate"),
        ("MAG_D", "T", "40", "de_novo_candidate"),
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
            }
            for i, (mag, allele, rose, origin) in enumerate(rows)
        ]
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

    def test_de_novo_in_a_replaced_mag_becomes_new_strain(self):
        got = attach_strain(self.rising, _tables())
        b = self._row(got, 8, "MAG_B", "C", "40")
        self.assertEqual((b.strain_class, b.origin), ("replaced", "new_strain"))
        al = self._row(got, 8, "MAG_A", "C", "AL")  # the AL arm has its own verdict
        self.assertEqual((al.strain_class, al.origin), ("replaced", "new_strain"))
        a = self._row(got, 8, "MAG_A", "T", "40")  # same letter, kept MAG
        self.assertEqual((a.strain_class, a.origin), ("not_replaced", "de_novo"))

    def test_only_de_novo_is_relabelled(self):
        got = attach_strain(self.rising, _tables())
        self.assertEqual(
            self._row(got, 8, "MAG_B", "A", "40").origin, "not_present_later"
        )

    def test_floors_differ(self):
        got = attach_strain(self.rising, _tables())
        b4 = self._row(got, 4, "MAG_B", "C", "40")
        self.assertEqual(
            (b4.strain_class, b4.origin), ("not_replaced", "de_novo_below_detection")
        )
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

    def test_raw_verdict_is_kept_next_to_the_relabel(self):
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
        self.assertEqual((a40.strain_class, a40.origin), ("replaced", "new_strain"))
        self.assertEqual(a40.strain_status, "40:not_replaced|AL:replaced")
        # MAG_B at floor 4: 40 not_replaced, AL never classified -> unverified
        b4 = self._row(got, 4, "MAG_B", "C", "40")
        self.assertEqual(b4.strain_class, "unverified")


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

    def test_without_strain_every_origin_but_new_strain_has_a_row(self):
        summ = summarise(self.rising)
        diet = summ[summ.rose_in == "40"]
        self.assertEqual(
            sorted(diet.origin), sorted(o for o in ORIGINS if o != "new_strain")
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
        summ = summarise(attach_strain(self.rising, _tables()))
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
        summ = summarise(attach_strain(self.rising, _tables()))
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
        self.out = os.path.join(self.tmp, "out")

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

    def test_with_classification_tables(self):
        done = self._run("--classification", *self.tables)
        self.assertEqual(done.returncode, 0, done.stderr)
        alleles = pd.read_csv(
            os.path.join(self.out, "rising_origin_alleles.tsv"), sep="\t"
        )
        summ = pd.read_csv(
            os.path.join(self.out, "rising_origin_summary.tsv"), sep="\t"
        )
        self.assertEqual(len(alleles), 2 * 7)
        self.assertEqual(sorted(summ.floor.unique()), [4, 8])
        self.assertEqual(
            int((alleles.origin == "new_strain").sum()), 3
        )  # 2 at 8, 1 at 4

    def test_without_classification(self):
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        summ = pd.read_csv(
            os.path.join(self.out, "rising_origin_summary.tsv"), sep="\t"
        )
        self.assertTrue(summ.floor.isna().all())
        self.assertNotIn("new_strain", set(summ.origin))


if __name__ == "__main__":
    unittest.main()
