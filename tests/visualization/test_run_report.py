"""Tests for run_report: the funnel, the per-cell significance table, the HTML.

A synthetic run directory (the data-type directory) with two MAGs, four samples
(two mice, treatment / control, pre / end), the bookkeeping tables the pipeline
writes, a cell summary and one figure.  ``config_runtime_*.yml`` sits in the
run root, one level up, as the workflow writes it.
"""

import base64
import os
import shutil
import subprocess
import tempfile
import unittest

import pandas as pd
import yaml

from alleleflux.scripts.visualization.run_report import (
    cell_significance,
    collect_figures,
    eligibility_counts,
    mags_with_profiles,
    newest_runtime_config,
    preprocessing_counts,
    render_report,
    sample_counts,
    site_filter_section,
)

CMP = "pre_end-treatment_control"
TP = "pre_end"
PNG_1PX = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=="
)


def _tsv(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def make_run(root, reuse=None, data_type="longitudinal"):
    """``root/longitudinal`` (or ``single``) plus the sheet and runtime config in ``root``.
    QC and inputMetadata go under ``reuse`` when given (a profile-reusing run)."""
    run = os.path.join(root, data_type)
    side = reuse or run
    samples = [
        ("s1", "m1", "treatment", "pre"),
        ("s2", "m1", "treatment", "end"),
        ("s3", "m2", "control", "pre"),
        ("s4", "m2", "control", "end"),
    ]
    sheet = os.path.join(root, "metadata.tsv")
    rows = [
        {"sample_id": s, "subjectID": m, "group": g, "time": t, "bam_path": "x"}
        for s, m, g, t in samples
    ]
    if data_type == "single":
        for r in rows:
            del r["time"]
    _tsv(sheet, rows)
    _tsv(
        os.path.join(
            side, "inputMetadata", f"inputMetadata_{TP}", "summary_mag_index.tsv"
        ),
        [{"mag_id": mag, "sample_id": s} for mag in ("M1", "M2") for s, *_ in samples],
    )
    # QC: M1 passes in 3 of 4 samples, M2 in 1 -> 2 MAGs with a passing sample, 4 passing pairs
    for mag, passes in (("M1", 3), ("M2", 1)):
        _tsv(
            os.path.join(side, "QC", f"QC_{TP}", f"{mag}_QC.tsv"),
            [
                {
                    "sample_id": s,
                    "MAG_ID": mag,
                    "breadth_threshold_passed": i < passes,
                    "coverage_threshold_passed": True,
                }
                for i, (s, *_) in enumerate(samples)
            ],
        )
    _tsv(
        os.path.join(run, f"eligibility_table_{CMP}.tsv"),
        [
            {
                "MAG_ID": "M1",
                "unpaired_test_eligible": True,
                "paired_test_eligible": True,
                "single_sample_eligible_treatment": True,
                "single_sample_eligible_control": False,
            },
            {
                "MAG_ID": "M2",
                "unpaired_test_eligible": False,
                "paired_test_eligible": False,
                "single_sample_eligible_treatment": True,
                "single_sample_eligible_control": False,
            },
        ],
    )
    _tsv(
        os.path.join(
            run,
            "preprocessing_eligibility",
            f"preprocessing_eligibility_between_groups_{CMP}.tsv",
        ),
        [
            {
                "MAG_ID": "M1",
                "positions_eligible_unpaired": 172,  # a count, not a flag
                "two_sample_unpaired_eligible": True,
                "two_sample_paired_eligible": True,
            }
        ],
    )
    if data_type == "longitudinal":
        _tsv(
            os.path.join(
                run,
                "preprocessing_eligibility",
                f"preprocessing_eligibility_within_groups_{CMP}.tsv",
            ),
            [
                {
                    "MAG_ID": "M1",
                    "single_sample_eligible_treatment": True,
                    "single_sample_eligible_control": False,
                },
                {
                    "MAG_ID": "M2",
                    "single_sample_eligible_treatment": True,
                    "single_sample_eligible_control": False,
                },
            ],
        )
    cells = [
        ("single_sample_tTest", "treatment", m, q, p)
        for m, q, p in (("M1", 0.001, 30.0), ("M2", 0.2, 0.0), ("M3", 0.049, 10.0))
    ]
    cells.append(("two_sample_paired_tTest", None, "M1", 0.5, 0.0))
    _tsv(
        os.path.join(
            run,
            "p_value_summary",
            "significant_sites_summary",
            "significant_sites_mag_cell_stats_long.tsv",
        ),
        [
            {
                "comparison": CMP,
                "test_type": t,
                "group_analyzed": g,
                "mag_id": m,
                "n_sig_q": int(q < 0.05),
                "min_q_value": q,
                "pct_sig_q": p,
            }
            for t, g, m, q, p in cells
        ],
    )
    _tsv(
        os.path.join(run, "site_depth_filter", "site_depth_filter_summary.tsv"),
        [
            {
                "sample_id": s,
                "mag_id": "M1",
                "timepoint": TP,
                "mean_coverage": 2.0,
                "mean_coverage_after_filter": 1.9,
                "n_sites_covered": cov,
                "n_removed_low": low,
                "n_removed_high": high,
                "n_kept": cov - low - high,
            }
            for s, cov, low, high in [("s1", 100, 5, 0), ("s2", 200, 0, 15)]
        ],
    )
    _tsv(
        os.path.join(run, "site_depth_filter", "site_depth_filter_summary_by_mag.tsv"),
        [
            {
                "mag_id": m,
                "timepoint": TP,
                "n_samples": 2,
                "sites_removed_high_total": 15,
                "sites_covered_total": 300,
                "fraction_removed_high": f,
            }
            for m, f in [("M1", 0.05), ("M2", 0.2)]
        ],
    )
    fig_dir = os.path.join(run, "figures", CMP)
    os.makedirs(fig_dir)
    open(os.path.join(fig_dir, "score_heatmap_ttest.png"), "wb").write(PNG_1PX)
    cfg = {
        "run_name": "demo_run",
        "input": {"metadata_path": sheet, "reuse_from": reuse or ""},
        "analysis": {"data_type": data_type},
        "permutation": {"enabled": False},
        "quality_control": {"site_depth_filter": {"enabled": True}},
    }
    yaml.safe_dump(
        cfg, open(os.path.join(root, "config_runtime_2026-01-01_00-00-00.yml"), "w")
    )
    return run, side, sheet


class TestTables(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.run_dir, self.side, self.sheet = make_run(self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_cell_significance_count_and_percent(self):
        cs = pd.read_csv(
            os.path.join(
                self.run_dir,
                "p_value_summary",
                "significant_sites_summary",
                "significant_sites_mag_cell_stats_long.tsv",
            ),
            sep="\t",
        )
        t = cell_significance(cs, CMP, 0.05).set_index(["test_type", "group_analyzed"])
        row = t.loc[("single_sample_tTest", "treatment")]
        self.assertEqual((int(row.n_tested), int(row.n_significant)), (3, 2))
        self.assertEqual(row.pct_significant, 66.7)
        self.assertNotIn("median_pct_sig_q", t.columns)
        self.assertEqual(int(t.loc[("two_sample_paired_tTest", "")].n_significant), 0)

    def test_sample_counts(self):
        t = sample_counts(self.sheet)
        self.assertEqual(len(t), 4)
        self.assertEqual(
            int(t.set_index(["group", "time"]).loc[("treatment", "pre"), "n_samples"]),
            1,
        )

    def test_sample_counts_single_sheet_without_time(self):
        run, _, sheet = make_run(
            os.path.join(self.tmp, "single_root"), data_type="single"
        )
        t = sample_counts(sheet)
        self.assertEqual(set(t["time"]), {""})
        self.assertEqual(int(t.set_index("group").loc["control", "n_samples"]), 2)

    def test_profiles_eligibility_preprocessing(self):
        self.assertEqual(
            mags_with_profiles(os.path.join(self.side, "inputMetadata"), TP), 2
        )
        self.assertEqual(
            eligibility_counts(
                os.path.join(self.run_dir, f"eligibility_table_{CMP}.tsv")
            ),
            {
                "unpaired_test_eligible": 1,
                "paired_test_eligible": 1,
                "single_sample_eligible_treatment": 2,
                "single_sample_eligible_control": 0,
            },
        )
        self.assertEqual(
            preprocessing_counts(
                os.path.join(self.run_dir, "preprocessing_eligibility"), CMP
            ),
            {
                "two_sample_unpaired_eligible": 1,
                "two_sample_paired_eligible": 1,
                "single_sample_eligible_treatment": 2,
                "single_sample_eligible_control": 0,
            },
        )

    def test_site_filter_section_none_when_absent(self):
        self.assertIsNone(site_filter_section(os.path.join(self.run_dir, "nope")))
        section = site_filter_section(os.path.join(self.run_dir, "site_depth_filter"))
        row = section["summary"].set_index("timepoint").loc[TP]
        self.assertEqual(row["sample_mag_pairs_checked"], 2)
        self.assertEqual(row["mags"], 1)
        self.assertEqual(row["sites_covered"], 300)
        self.assertEqual(row["sites_removed"], 20)
        self.assertAlmostEqual(row["pct_removed"], 100 * 20 / 300)
        self.assertAlmostEqual(row["median_row_pct"], 100 * (0.05 + 0.075) / 2)
        self.assertAlmostEqual(row["max_row_pct"], 7.5)
        self.assertTrue(section["histogram_png"].startswith(b"\x89PNG"))

    def test_collect_figures(self):
        figs = collect_figures(os.path.join(self.run_dir, "figures"))
        self.assertEqual(
            [(f["comparison"], f["family"], f["kind"]) for f in figs],
            [(CMP, "ttest", "heatmap")],
        )
        self.assertTrue(figs[0]["png_b64"])

    def test_newest_runtime_config(self):
        yaml.safe_dump(
            {"run_name": "newer"},
            open(os.path.join(self.tmp, "config_runtime_2026-02-01_00-00-00.yml"), "w"),
        )
        self.assertEqual(newest_runtime_config(self.tmp)["run_name"], "newer")
        self.assertIsNone(newest_runtime_config(self.run_dir))


class TestRender(unittest.TestCase):
    def _ctx(self, **over):
        ctx = {
            "run_name": "demo_run",
            "run_dir": "/r",
            "date": "2026-01-01",
            "version": "1.0",
            "commit": "",
            "config_text": "run_name: demo_run",
            "permuted": False,
            "data_type": "longitudinal",
            "comparisons": [
                {
                    "name": CMP,
                    "funnel": [("samples in the metadata sheet (treatment, pre)", 1)],
                    "significance": pd.DataFrame(
                        [
                            {
                                "test_type": "single_sample_tTest",
                                "group_analyzed": "treatment",
                                "n_tested": 3,
                                "n_significant": 2,
                                "pct_significant": 66.7,
                            }
                        ]
                    ),
                    "preprocessing": None,
                }
            ],
            "site_filter": None,
            "figures": [
                {
                    "comparison": CMP,
                    "family": "ttest",
                    "kind": "heatmap",
                    "caption": "colour = score",
                    "png_b64": base64.b64encode(PNG_1PX).decode(),
                }
            ],
            "q_threshold": 0.05,
        }
        ctx.update(over)
        return ctx

    def test_sections_and_embedded_image(self):
        html = render_report(self._ctx())
        for h in ("Run", "Funnel", "Significant sites", "Site-depth filter", "Figures"):
            self.assertIn(h, html)
        self.assertIn("data:image/png;base64,", html)
        self.assertIn("<td>66.7", html)
        self.assertNotIn("<script", html)

    def test_missing_sections_are_announced_not_fatal(self):
        html = render_report(self._ctx())
        self.assertIn("not available in this run", html)

    def test_integers_get_thousands_separators(self):
        ctx = self._ctx()
        ctx["site_filter"] = {
            "summary": pd.DataFrame(
                {
                    "timepoint": ["end"],
                    "sites_covered": [5279499787],
                    "pct_removed": [1.04],
                }
            ),
            "histogram_b64": "",
        }
        ctx["comparisons"][0]["funnel"] = [
            ("MAGs with profiles (pre_end)", 6675),
            ("missing", None),
        ]
        html = render_report(ctx)
        self.assertIn("5,279,499,787", html)
        self.assertIn("<td>1.0</td>", html)
        self.assertIn("6,675", html)
        self.assertIn("not available in this run", html)

    def test_permuted_run_banner(self):
        self.assertIn("permuted (null) run", render_report(self._ctx(permuted=True)))
        self.assertNotIn("permuted (null) run", render_report(self._ctx()))


class TestRunReportCLI(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self, run, *extra):
        out = os.path.join(self.tmp, "report.html")
        done = subprocess.run(
            [
                "alleleflux-run-report",
                "--run_dir",
                run,
                "--figures_dir",
                os.path.join(run, "figures"),
                "--output",
                out,
                "--comparisons",
                CMP,
                *extra,
            ],
            capture_output=True,
            text=True,
        )
        return done, out

    def test_funnel_reads_explicit_metadata_dir_and_names_the_tests(self):
        reuse = os.path.join(self.tmp, "reuse")
        run, side, sheet = make_run(os.path.join(self.tmp, "root"), reuse=reuse)
        done, out = self._run(
            run,
            "--metadata_dir",
            os.path.join(reuse, "inputMetadata"),
            "--metadata_path",
            sheet,
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        html = open(out).read()
        self.assertIn("<th>MAGs with profiles (pre_end)</th><td>2</td>", html)
        self.assertIn(
            "<th>MAGs with enough samples for the paired between-group tests and CMH (divergence)</th><td>1</td>",
            html,
        )
        self.assertIn(
            "<th>MAGs with enough sites after preprocessing for the within-group tests in treatment (parallelism)</th><td>2</td>",
            html,
        )
        self.assertNotIn("QC-passing", html)
        self.assertIn("<td>66.7", html)

    def test_defaults_from_runtime_config(self):
        reuse = os.path.join(self.tmp, "reuse")
        run, _, _ = make_run(os.path.join(self.tmp, "root"), reuse=reuse)
        done, out = self._run(
            run
        )  # no dirs given: the runtime config's reuse_from decides
        self.assertEqual(done.returncode, 0, done.stderr)
        html = open(out).read()
        # the reuse dir's inputMetadata was found through the runtime config
        self.assertIn("<th>MAGs with profiles (pre_end)</th><td>2</td>", html)
        self.assertIn("demo_run", html)

    def test_report_renders_for_a_single_run(self):
        run, _, _ = make_run(os.path.join(self.tmp, "root"), data_type="single")
        done, out = self._run(run, "--data_type", "single")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("Funnel", open(out).read())


if __name__ == "__main__":
    unittest.main()
