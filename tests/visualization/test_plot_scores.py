"""Tests for plot_scores: the drawing functions and the alleleflux-plot-scores command.

A synthetic run directory with four MAGs: M1, M2, M3 (phylum P1) tested in every
ttest and lmm cell, M4 (phylum P2) tested only in the treatment parallelism cell.
"""

import os
import shutil
import subprocess
import tempfile
import unittest

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

from alleleflux.scripts.visualization.plot_scores import (
    _long_scores,
    draw_score_boxplot,
    draw_score_boxplot_by_phylum,
    draw_score_heatmap,
)
from alleleflux.scripts.visualization.score_panels import (
    load_family,
    resolve_family,
    total_column,
)

CMP = "pre_end-treatment_control"
TAX = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def _scores(directory, name, rows, score_col):
    recs = []
    for mag, phylum, score, total in rows:
        rec = {"MAG_ID": mag, "domain": "Bacteria", "phylum": phylum}
        rec.update({c: f"{c}_of_{mag}" for c in TAX[2:]})
        rec[total_column(score_col)] = total
        rec[score_col] = score
        recs.append(rec)
    pd.DataFrame(recs).to_csv(os.path.join(directory, name), sep="\t", index=False)


def make_run(run):
    """scores/processed/combined/MAG + the cell summary, ttest and lmm families."""
    sd = os.path.join(run, "scores", "processed", "combined", "MAG")
    os.makedirs(sd)
    three = [("M1", "P1", 10, 100), ("M2", "P1", 50, 100), ("M3", "P1", 90, 100)]
    four = three + [("M4", "P2", 30, 40)]
    _scores(
        sd,
        f"scores_single_sample-{CMP}_treatment-MAGs.tsv",
        four,
        "score_tTest_treatment (%)",
    )
    _scores(
        sd,
        f"scores_single_sample-{CMP}_control-MAGs.tsv",
        three,
        "score_tTest_control (%)",
    )
    _scores(sd, f"scores_two_sample_paired-{CMP}-MAGs.tsv", three, "score_tTest (%)")
    _scores(
        sd, f"scores_lmm_across_time-{CMP}_treatment-MAGs.tsv", four, "score_LMM (%)"
    )
    _scores(
        sd, f"scores_lmm_across_time-{CMP}_control-MAGs.tsv", three, "score_LMM (%)"
    )
    _scores(sd, f"scores_lmm-{CMP}-MAGs.tsv", three, "score_LMM (%)")
    rows = []
    for tt, g, mags in (
        ("single_sample_tTest", "treatment", ["M1", "M2", "M3", "M4"]),
        ("single_sample_tTest", "control", ["M1", "M2", "M3"]),
        ("two_sample_paired_tTest", None, ["M1", "M2", "M3"]),
        ("LMM_across_time", "treatment", ["M1", "M2", "M3", "M4"]),
        ("LMM_across_time", "control", ["M1", "M2", "M3"]),
        ("LMM", None, ["M1", "M2", "M3"]),
    ):
        for i, m in enumerate(mags):
            rows.append(
                {
                    "comparison": CMP,
                    "test_type": tt,
                    "group_analyzed": g,
                    "mag_id": m,
                    "n_sig_q": i,
                    "min_q_value": [0.0004, 0.004, 0.04, 0.4][i],
                }
            )
    cs = os.path.join(run, "p_value_summary", "significant_sites_summary")
    os.makedirs(cs)
    pd.DataFrame(rows).to_csv(
        os.path.join(cs, "significant_sites_mag_cell_stats_long.tsv"),
        sep="\t",
        index=False,
    )
    return sd, os.path.join(cs, "significant_sites_mag_cell_stats_long.tsv")


class TestDrawing(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.run_dir = os.path.join(self.tmp, "run")
        sd, cs = make_run(self.run_dir)
        fam = resolve_family(sd, CMP, ["treatment", "control"], "ttest")
        self.data = load_family(fam, pd.read_csv(cs, sep="\t"), CMP)
        self.out = os.path.join(self.tmp, "out")
        os.makedirs(self.out)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_heatmap_writes_a_file(self):
        path = os.path.join(self.out, "h.png")
        draw_score_heatmap(self.data, ["M1", "M2", "M3", "M4"], [], "t", path, 0.05)
        self.assertGreater(os.path.getsize(path), 1000)

    def test_heatmap_legend_lists_the_phyla_and_labels_follow_order(self):
        path = os.path.join(self.out, "h.png")
        fig = draw_score_heatmap(self.data, ["M3", "M1", "M2"], ["M4"], "t", path, 0.05)
        legend = fig.axes[1].get_legend()
        self.assertEqual([t.get_text() for t in legend.get_texts()], ["P1", "P2"])
        labels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
        self.assertEqual([l.split()[0] for l in labels], ["M3", "M1", "M2", "M4"])

    def test_heatmap_with_tree_order_and_absent_rows(self):
        path = os.path.join(self.out, "h.png")
        draw_score_heatmap(self.data, ["M3", "M1", "M2"], ["M4"], "t", path, 0.05)
        self.assertGreater(os.path.getsize(path), 1000)

    def test_long_frame_marks_pairs_and_keeps_one_jitter_per_mag(self):
        long = _long_scores(self.data, np.random.default_rng(0))
        # M1, M2, M3 have both parallelism scores; M4 only treatment -> no pair
        paired = long[long["paired"]]
        self.assertEqual(sorted(paired["mag"].unique()), ["M1", "M2", "M3"])
        self.assertEqual(len(paired), 6)  # 3 MAGs x 2 parallelism rows
        self.assertEqual(int((long["mag"] == "M4").sum()), 1)
        # one jitter per MAG: x - cell_idx is identical across a MAG's rows
        off = (long["x"] - long["cell_idx"]).round(9).groupby(long["mag"]).nunique()
        self.assertTrue((off == 1).all())
        self.assertEqual(int((long["cell"] == "treatment parallelism").sum()), 4)
        self.assertEqual(int((long["cell"] == "divergence").sum()), 3)

    def test_boxplot_writes_a_file(self):
        path = os.path.join(self.out, "b.png")
        draw_score_boxplot(self.data, "t", path)
        self.assertGreater(os.path.getsize(path), 1000)

    def test_by_phylum_pools_small_phyla(self):
        path = os.path.join(self.out, "p.png")
        p = draw_score_boxplot_by_phylum(self.data, "t", path, 3)
        panels = list(p.data["panel"].cat.categories)
        self.assertEqual(panels, ["P1 (n=3)", "other (1 phyla, n=1)"])
        self.assertGreater(os.path.getsize(path), 1000)


class TestPlotScoresCLI(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.run_dir = os.path.join(self.tmp, "run")
        self.sd, _ = make_run(self.run_dir)
        self.out = os.path.join(self.tmp, "figs")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self, *extra):
        return subprocess.run(
            [
                "alleleflux-plot-scores",
                "--run_dir",
                self.run_dir,
                "--comparison",
                CMP,
                "--groups",
                "treatment",
                "control",
                "--output_dir",
                self.out,
                *extra,
            ],
            capture_output=True,
            text=True,
        )

    def test_writes_figures_and_sentinel_for_present_families_only(self):
        done = self._run("--families", "ttest", "lmm", "cmh")
        self.assertEqual(done.returncode, 0, done.stderr)
        files = sorted(os.listdir(self.out))
        self.assertEqual(
            files,
            [
                "score_boxplot_lmm.png",
                "score_boxplot_lmm_by_phylum.png",
                "score_boxplot_ttest.png",
                "score_boxplot_ttest_by_phylum.png",
                "score_figures.done",
                "score_heatmap_lmm.png",
                "score_heatmap_ttest.png",
            ],
        )
        listed = open(os.path.join(self.out, "score_figures.done")).read().split()
        self.assertEqual(
            sorted(os.path.basename(f) for f in listed),
            [f for f in files if f != "score_figures.done"],
        )
        self.assertIn("cmh", done.stderr)  # skipped, named

    def _cell_stats_path(self):
        return os.path.join(
            self.run_dir,
            "p_value_summary",
            "significant_sites_summary",
            "significant_sites_mag_cell_stats_long.tsv",
        )

    def test_missing_within_file_is_a_grey_column_not_a_skip(self):
        os.remove(os.path.join(self.sd, f"scores_single_sample-{CMP}_control-MAGs.tsv"))
        # no control score table means no control summary rows either
        cs = pd.read_csv(self._cell_stats_path(), sep="\t")
        cs[cs.group_analyzed != "control"].to_csv(
            self._cell_stats_path(), sep="\t", index=False
        )
        done = self._run("--families", "ttest")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertTrue(
            os.path.exists(os.path.join(self.out, "score_heatmap_ttest.png"))
        )

    def test_single_data_type_draws_divergence_only(self):
        for f in os.listdir(self.sd):
            if "single_sample" in f or "across_time" in f:
                os.remove(os.path.join(self.sd, f))
        done = self._run("--families", "ttest", "--data_type", "single")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertTrue(
            os.path.exists(os.path.join(self.out, "score_heatmap_ttest.png"))
        )
        self.assertTrue(
            os.path.exists(os.path.join(self.out, "score_boxplot_ttest.png"))
        )

    def test_unpaired_family_draws_when_its_table_exists(self):
        _scores(
            self.sd,
            f"scores_two_sample_unpaired-{CMP}-MAGs.tsv",
            [("M1", "P1", 10, 100), ("M2", "P1", 50, 100), ("M3", "P1", 90, 100)],
            "score_tTest (%)",
        )
        cs = pd.read_csv(self._cell_stats_path(), sep="\t")
        extra = cs[cs.test_type == "two_sample_paired_tTest"].assign(
            test_type="two_sample_unpaired_tTest"
        )
        pd.concat([cs, extra]).to_csv(self._cell_stats_path(), sep="\t", index=False)
        done = self._run("--families", "ttest_unpaired")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertTrue(
            os.path.exists(os.path.join(self.out, "score_heatmap_ttest_unpaired.png"))
        )

    def test_group_labels_with_underscore(self):
        cmp_ = "pre_end-high_fat_control"
        _scores(
            self.sd,
            f"scores_single_sample-{cmp_}_high_fat-MAGs.tsv",
            [("M1", "P1", 1, 9)],
            "score_tTest_high_fat (%)",
        )
        _scores(
            self.sd,
            f"scores_two_sample_paired-{cmp_}-MAGs.tsv",
            [("M1", "P1", 1, 9)],
            "score_tTest (%)",
        )
        cs = pd.read_csv(self._cell_stats_path(), sep="\t")
        extra = pd.DataFrame(
            [
                {
                    "comparison": cmp_,
                    "test_type": "single_sample_tTest",
                    "group_analyzed": "high_fat",
                    "mag_id": "M1",
                    "n_sig_q": 1,
                    "min_q_value": 0.01,
                },
                {
                    "comparison": cmp_,
                    "test_type": "two_sample_paired_tTest",
                    "group_analyzed": None,
                    "mag_id": "M1",
                    "n_sig_q": 1,
                    "min_q_value": 0.01,
                },
            ]
        )
        pd.concat([cs, extra]).to_csv(self._cell_stats_path(), sep="\t", index=False)
        done = subprocess.run(
            [
                "alleleflux-plot-scores",
                "--run_dir",
                self.run_dir,
                "--comparison",
                cmp_,
                "--groups",
                "high_fat",
                "control",
                "--output_dir",
                self.out,
                "--families",
                "ttest",
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertTrue(
            os.path.exists(os.path.join(self.out, "score_heatmap_ttest.png"))
        )

    def test_pdf_format_still_writes_the_png_the_report_needs(self):
        done = self._run("--families", "ttest", "--format", "pdf")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertTrue(
            os.path.exists(os.path.join(self.out, "score_heatmap_ttest.pdf"))
        )
        self.assertTrue(
            os.path.exists(os.path.join(self.out, "score_heatmap_ttest.png"))
        )


if __name__ == "__main__":
    unittest.main()
