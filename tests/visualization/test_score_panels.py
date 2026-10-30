"""Tests for score_panels: a test family -> its three cells -> one MAG x cell grid.

Synthetic score tables shaped like the pipeline's combined MAG tables
(``MAG_ID``, taxonomy, ``score_<X> (%)``, ``total_sites_per_group_<X>``) and a
synthetic cell summary (``significant_sites_mag_cell_stats_long.tsv``).
Example groups are ``treatment`` / ``control``.
"""

import os
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd

from alleleflux.scripts.visualization.score_panels import (
    load_family,
    q_stars,
    resolve_family,
    taxonomy_order,
    total_column,
    tree_order,
)

TAX = ["domain", "phylum", "class", "order", "family", "genus", "species"]
CMP = "pre_end-treatment_control"


def _write_scores(directory, name, rows, score_col):
    """rows: (mag, phylum, score, total_sites); taxonomy below phylum is 'x'."""
    recs = []
    for mag, phylum, score, total in rows:
        rec = {"MAG_ID": mag, "domain": "Bacteria", "phylum": phylum}
        rec.update({c: "x" for c in TAX[2:]})
        rec[total_column(score_col)] = total
        rec[score_col] = score
        recs.append(rec)
    pd.DataFrame(recs).to_csv(os.path.join(directory, name), sep="\t", index=False)


def _cell_stats(rows):
    """rows: (test_type, group_analyzed or None, mag, min_q)."""
    return pd.DataFrame(
        [
            {
                "comparison": CMP,
                "test_type": t,
                "group_analyzed": g,
                "mag_id": m,
                "min_q_value": q,
            }
            for t, g, m, q in rows
        ]
    )


class TestResolveFamily(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        _write_scores(
            self.dir,
            f"scores_single_sample-{CMP}_treatment-MAGs.tsv",
            [("M1", "P1", 10, 100)],
            "score_tTest_treatment (%)",
        )
        _write_scores(
            self.dir,
            f"scores_single_sample-{CMP}_control-MAGs.tsv",
            [("M1", "P1", 5, 100)],
            "score_tTest_control (%)",
        )
        _write_scores(
            self.dir,
            f"scores_two_sample_paired-{CMP}-MAGs.tsv",
            [("M1", "P1", 20, 100)],
            "score_tTest (%)",
        )

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_ttest_resolves_three_cells(self):
        fam = resolve_family(self.dir, CMP, ["treatment", "control"], "ttest")
        self.assertEqual(
            [c.label for c in fam.cells],
            ["treatment parallelism", "control parallelism", "divergence"],
        )
        self.assertEqual(fam.cells[0].score_column, "score_tTest_treatment (%)")
        self.assertEqual(
            (fam.cells[0].test_type, fam.cells[0].group_analyzed),
            ("single_sample_tTest", "treatment"),
        )
        self.assertEqual(
            (fam.cells[2].test_type, fam.cells[2].group_analyzed),
            ("two_sample_paired_tTest", ""),
        )
        self.assertTrue(all(c.score_file for c in fam.cells))

    def test_group_label_with_underscore_resolves_files(self):
        cmp_ = "pre_end-high_fat_control"
        _write_scores(
            self.dir,
            f"scores_single_sample-{cmp_}_high_fat-MAGs.tsv",
            [("M1", "P1", 1, 9)],
            "score_tTest_high_fat (%)",
        )
        fam = resolve_family(self.dir, cmp_, ["high_fat", "control"], "ttest")
        self.assertTrue(fam.cells[0].score_file.endswith("_high_fat-MAGs.tsv"))
        self.assertEqual(fam.cells[0].score_column, "score_tTest_high_fat (%)")
        self.assertIsNone(
            fam.cells[1].score_file
        )  # no control file for this comparison

    def test_missing_within_file_gives_a_cell_without_file(self):
        os.remove(
            os.path.join(self.dir, f"scores_single_sample-{CMP}_treatment-MAGs.tsv")
        )
        fam = resolve_family(self.dir, CMP, ["treatment", "control"], "ttest")
        self.assertIsNone(fam.cells[0].score_file)
        self.assertIsNotNone(fam.cells[1].score_file)

    def test_unpaired_family_uses_the_unpaired_table(self):
        _write_scores(
            self.dir,
            f"scores_two_sample_unpaired-{CMP}-MAGs.tsv",
            [("M1", "P1", 7, 100)],
            "score_MannWhitney (%)",
        )
        fam = resolve_family(self.dir, CMP, ["treatment", "control"], "mannwhitney")
        div = fam.cells[2]
        self.assertTrue(
            div.score_file.endswith(f"scores_two_sample_unpaired-{CMP}-MAGs.tsv")
        )
        self.assertEqual(div.score_column, "score_MannWhitney (%)")
        self.assertEqual(div.test_type, "two_sample_unpaired_MannWhitney")
        # the within-group cells are the same single-sample Wilcoxon cells
        self.assertEqual(fam.cells[0].test_type, "single_sample_Wilcoxon")
        # no unpaired t-test column in that table -> ttest_unpaired still resolves
        # (file exists); loading it is what would fail, loudly
        self.assertIsNotNone(
            resolve_family(self.dir, CMP, ["treatment", "control"], "ttest_unpaired")
        )

    def test_no_file_at_all_is_none(self):
        self.assertIsNone(
            resolve_family(self.dir, CMP, ["treatment", "control"], "lmm")
        )

    def test_single_data_type_has_divergence_only(self):
        fam = resolve_family(
            self.dir, CMP, ["treatment", "control"], "ttest", data_type="single"
        )
        self.assertEqual([c.kind for c in fam.cells], ["divergence"])

    def test_cmh_one_divergence_cell_per_focus_file(self):
        _write_scores(
            self.dir,
            f"scores_cmh-{CMP}-MAGs-end.tsv",
            [("M1", "P1", 1, 9)],
            "score_CMH (%)",
        )
        _write_scores(
            self.dir,
            f"scores_cmh-{CMP}-MAGs-pre.tsv",
            [("M1", "P1", 1, 9)],
            "score_CMH (%)",
        )
        fam = resolve_family(self.dir, CMP, ["treatment", "control"], "cmh")
        div = [c for c in fam.cells if c.kind == "divergence"]
        self.assertEqual(
            [c.label for c in div], ["divergence (end)", "divergence (pre)"]
        )
        self.assertEqual(
            [(c.test_type, c.group_analyzed) for c in div],
            [("CMH", "end"), ("CMH", "pre")],
        )
        self.assertEqual(fam.cells[0].test_type, "CMH_across_time")


class TestLoadFamily(unittest.TestCase):
    """M1 in all three files; M2 only in the treatment file; M3 in the paired file with
    zero sites (the CMH-style 'insufficient data' row)."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        _write_scores(
            self.dir,
            f"scores_single_sample-{CMP}_treatment-MAGs.tsv",
            [("M1", "P1", 10, 100), ("M2", "P2", 30, 50)],
            "score_tTest_treatment (%)",
        )
        _write_scores(
            self.dir,
            f"scores_single_sample-{CMP}_control-MAGs.tsv",
            [("M1", "P1", 5, 100)],
            "score_tTest_control (%)",
        )
        _write_scores(
            self.dir,
            f"scores_two_sample_paired-{CMP}-MAGs.tsv",
            [("M1", "P1", 20, 100), ("M3", "P1", 0, 0)],
            "score_tTest (%)",
        )
        self.fam = resolve_family(self.dir, CMP, ["treatment", "control"], "ttest")
        self.cs = _cell_stats(
            [
                ("single_sample_tTest", "treatment", "M1", 0.004),
                ("single_sample_tTest", "treatment", "M2", 0.5),
                ("single_sample_tTest", "control", "M1", 0.9),
                ("two_sample_paired_tTest", None, "M1", 0.2),
            ]
        )
        self.data = load_family(self.fam, self.cs, CMP)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_union_keeps_a_mag_tested_in_one_cell(self):
        self.assertEqual(sorted(self.data.scores.index), ["M1", "M2"])
        self.assertEqual(self.data.scores.loc["M2", "treatment parallelism"], 30)
        self.assertTrue(np.isnan(self.data.scores.loc["M2", "control parallelism"]))
        self.assertTrue(np.isnan(self.data.min_q.loc["M2", "divergence"]))

    def test_zero_total_sites_is_untested(self):
        self.assertNotIn("M3", self.data.scores.index)

    def test_min_q_comes_from_the_matching_cell(self):
        self.assertAlmostEqual(
            self.data.min_q.loc["M1", "treatment parallelism"], 0.004
        )
        self.assertAlmostEqual(self.data.min_q.loc["M1", "control parallelism"], 0.9)
        self.assertAlmostEqual(self.data.min_q.loc["M1", "divergence"], 0.2)

    def test_score_without_cell_stats_row_raises(self):
        # a tested MAG always has a summary row (both come from the same p-value
        # file); a missing one is a broken run, not a blank star
        with self.assertRaisesRegex(ValueError, "M2"):
            load_family(self.fam, self.cs[self.cs.mag_id != "M2"], CMP)

    def test_cell_stats_row_without_score_raises(self):
        extra = _cell_stats([("single_sample_tTest", "control", "M2", 0.1)])
        with self.assertRaisesRegex(ValueError, "M2"):
            load_family(self.fam, pd.concat([self.cs, extra]), CMP)

    def test_duplicate_cell_stats_row_raises(self):
        dup = _cell_stats([("single_sample_tTest", "treatment", "M1", 0.5)])
        with self.assertRaisesRegex(ValueError, "M1"):
            load_family(self.fam, pd.concat([self.cs, dup]), CMP)

    def test_missing_total_sites_column_raises(self):
        path = os.path.join(self.dir, f"scores_single_sample-{CMP}_control-MAGs.tsv")
        df = pd.read_csv(path, sep="\t")
        df.drop(columns=[total_column("score_tTest_control (%)")]).to_csv(
            path, sep="\t", index=False
        )
        with self.assertRaisesRegex(ValueError, "total_sites_per_group_tTest_control"):
            load_family(self.fam, self.cs, CMP)

    def test_cell_without_file_is_all_nan_column(self):
        os.remove(
            os.path.join(self.dir, f"scores_single_sample-{CMP}_control-MAGs.tsv")
        )
        fam = resolve_family(self.dir, CMP, ["treatment", "control"], "ttest")
        # no control score table means no control p-value files, so no summary rows
        cs = self.cs[self.cs.group_analyzed != "control"]
        data = load_family(fam, cs, CMP)
        self.assertTrue(data.scores["control parallelism"].isna().all())
        self.assertEqual(list(data.scores.columns), [c.label for c in fam.cells])

    def test_taxonomy_is_taken_from_any_file(self):
        self.assertEqual(self.data.taxonomy.loc["M2", "phylum"], "P2")
        self.assertEqual(list(self.data.taxonomy.columns), TAX)


class TestOrdering(unittest.TestCase):
    def test_taxonomy_order_sorts_phylum_then_family_then_id(self):
        tax = pd.DataFrame(
            {
                "phylum": ["B", "A", "A", "A"],
                "class": "c",
                "order": "o",
                "family": ["f", "f2", "f1", "f1"],
                "genus": "g",
                "species": "s",
            },
            index=["M4", "M3", "M2", "M1"],
        )
        self.assertEqual(taxonomy_order(tax), ["M1", "M2", "M3", "M4"])

    def test_taxonomy_order_puts_missing_taxonomy_last(self):
        tax = pd.DataFrame(
            {
                "phylum": [np.nan, "A"],
                "class": "c",
                "order": "o",
                "family": "f",
                "genus": "g",
                "species": "s",
            },
            index=["M9", "M1"],
        )
        self.assertEqual(taxonomy_order(tax), ["M1", "M9"])

    def test_tree_order_follows_tips_and_lists_absent(self):
        d = tempfile.mkdtemp()
        path = os.path.join(d, "t.nwk")
        open(path, "w").write("((M2:1,M1:1):1,M3:2);\n")
        in_tree, absent = tree_order(["M1", "M2", "M3", "M4"], path)
        shutil.rmtree(d)
        self.assertEqual(set(in_tree), {"M1", "M2", "M3"})
        self.assertEqual(in_tree.index("M3") in (0, 2), True)  # M3 is the outer tip
        self.assertEqual(absent, ["M4"])


class TestQStars(unittest.TestCase):
    def test_boundaries(self):
        self.assertEqual(q_stars(0.0009, 0.05), "***")
        self.assertEqual(q_stars(0.001, 0.05), "**")
        self.assertEqual(q_stars(0.0099, 0.05), "**")
        self.assertEqual(q_stars(0.01, 0.05), "*")
        self.assertEqual(q_stars(0.049, 0.05), "*")
        self.assertEqual(q_stars(0.05, 0.05), "")
        self.assertEqual(q_stars(float("nan"), 0.05), "")


class TestTotalColumn(unittest.TestCase):
    def test_mapping(self):
        self.assertEqual(
            total_column("score_tTest_control (%)"),
            "total_sites_per_group_tTest_control",
        )
        self.assertEqual(total_column("score_CMH (%)"), "total_sites_per_group_CMH")


if __name__ == "__main__":
    unittest.main()
