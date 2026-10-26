"""Tests for ``alleleflux-baseline-presence``: was the significant allele already
there at the baseline timepoint?

Unit tests cover each pure step; the end-to-end test builds a miniature run
directory (summary, per-MAG test file, profiles, metadata, reference) and runs
the real console script.
"""

import gzip
import os
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
import pandas as pd

from alleleflux.scripts.analysis.ani.null_model import build_error_model
from alleleflux.scripts.analysis.baseline_presence import (
    BASELINE_SCOPE_DEFAULT,
    EVIDENCE_ABSENT,
    EVIDENCE_BELOW_DETECTION,
    EVIDENCE_NOT_COVERED,
    EVIDENCE_PRESENT,
    ORIGIN_STANDING,
    Timepoints,
    label_origin_in_own_mouse,
    allele_p_column,
    assess_allele,
    find_summary_file,
    output_label,
    candidate_alleles,
    load_significant_sites,
    restrict_to_own_group,
    frequency_change_by_group,
    summarise_rising,
    rising_columns,
    RISING_COLUMNS,
    check_baseline_scope,
    SUMMARY_FAMILIES,
    WITHIN_GROUP_FAMILIES,
    parse_comparison,
    summarise_sites,
    summary_columns,
)

MODEL = build_error_model(min_base_quality=30, fdr=1e-6)  # bar = 3 reads for 5..99x


class TestParseComparison(unittest.TestCase):
    def test_splits_period_and_groups(self):
        self.assertEqual(
            parse_comparison("pre_end-fat_control"), ("pre", "end", "fat", "control")
        )
        self.assertEqual(
            parse_comparison("5mo_22mo-40_AL"), ("5mo", "22mo", "40", "AL")
        )

    def test_rejects_malformed_label(self):
        for bad in ("pre_end", "pre-end-fat_control", "pre_mid_end-fat_control"):
            with self.assertRaises(ValueError):
                parse_comparison(bad)


class TestFindSummaryFile(unittest.TestCase):
    def test_family_prefix_does_not_swallow_across_time(self):
        tmp = tempfile.mkdtemp()
        for name in (
            "p_value_summary_lmm_pre_end.tsv",
            "p_value_summary_lmm_across_time_pre_end.tsv",
            "p_value_summary_two_sample_paired_5mo_22mo-40_AL.tsv",
        ):
            open(os.path.join(tmp, name), "w").write("x\n")
        self.assertTrue(
            find_summary_file(tmp, "lmm", "pre_end").endswith(
                "p_value_summary_lmm_pre_end.tsv"
            )
        )
        self.assertTrue(
            find_summary_file(tmp, "lmm_across_time", "pre_end").endswith(
                "across_time_pre_end.tsv"
            )
        )
        self.assertTrue(
            find_summary_file(tmp, "two_sample_paired", "5mo_22mo").endswith(
                "-40_AL.tsv"
            )
        )  # DRiDO naming
        with self.assertRaises(FileNotFoundError):
            find_summary_file(tmp, "two_sample_unpaired", "pre_end")
        shutil.rmtree(tmp)


class TestOutputLabel(unittest.TestCase):
    def test_family_is_always_in_the_name_and_never_repeated(self):
        self.assertEqual(
            output_label("two_sample_paired", "two_sample_paired_tTest"),
            "two_sample_paired_tTest",
        )
        self.assertEqual(output_label("lmm", "LMM_abs"), "lmm_LMM_abs")
        self.assertEqual(
            output_label("lmm_across_time", "LMM"), "lmm_across_time_LMM"
        )  # distinct from lmm's


class TestLoadSignificantSites(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.path = os.path.join(
            self.tmp, "p_value_summary_two_sample_paired_pre_end.tsv"
        )
        pd.DataFrame(
            {
                "period": ["pre_end"] * 4,
                "mag_id": ["MAG_A"] * 4,
                "contig": ["c1"] * 4,
                "position": [10, 10, 20, 30],
                "gene_id": [
                    "g1 ",
                    "g1 ",
                    "",
                    "g2 ",
                ],  # real files carry a trailing space
                "test_type": [
                    "two_sample_paired_tTest",
                    "two_sample_paired_Wilcoxon",
                    "two_sample_paired_tTest",
                    "two_sample_paired_tTest",
                ],
                "min_p_value": [0.001, 0.002, 0.03, 0.2],
                "source_file": ["MAG_A_two_sample_paired.tsv.gz"] * 4,
                "q_value": [0.01, 0.02, 0.06, 0.5],
            }
        ).to_csv(self.path, sep="\t", index=False)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_filters_to_the_test_type_and_the_threshold_column(self):
        by_q = load_significant_sites(
            self.path, "two_sample_paired_tTest", "q_value", 0.05
        )
        self.assertEqual(
            by_q.position.tolist(), [10]
        )  # Wilcoxon row and q=0.06 row dropped
        self.assertEqual(by_q.gene_id.tolist(), ["g1"])  # stripped
        by_p = load_significant_sites(
            self.path, "two_sample_paired_tTest", "min_p_value", 0.05
        )
        self.assertEqual(by_p.position.tolist(), [10, 20])  # p=0.03 passes on p
        self.assertEqual(
            by_p.group_analyzed.tolist(), ["", ""]
        )  # two-sample: no group column

    def test_any_test_type_in_the_file_is_selectable(self):
        wilcoxon = load_significant_sites(
            self.path, "two_sample_paired_Wilcoxon", "q_value", 0.05
        )
        self.assertEqual(wilcoxon.position.tolist(), [10])

    def test_wrong_test_type_finds_nothing_and_says_so(self):
        with self.assertRaisesRegex(ValueError, "single_sample_tTest"):
            load_significant_sites(self.path, "single_sample_tTest", "q_value", 0.05)

    def test_test_type_is_required(self):
        with self.assertRaisesRegex(
            ValueError, "Wilcoxon"
        ):  # None -> error listing what the file holds
            load_significant_sites(self.path, None, "q_value", 0.05)


class TestCandidateAlleles(unittest.TestCase):
    """The candidate rule: every base with a real test (p < 1) at a two-letter
    site; at a site with three or four tested bases, only those significant on
    their own (p <= the comparison's BH line for that group).  No tie tolerance."""

    # staticmethod: a plain function stored on a class would be bound as a method
    TWO = staticmethod(allele_p_column("two_sample_paired", "two_sample_paired_tTest"))

    def test_two_letter_site_keeps_both_bases_even_when_float32_noise_breaks_the_tie(
        self,
    ):
        # position 10: exact twins; position 20: twins that differ in the 9th digit
        # (the cache stores float32, so the seesaw partners' p never tie bit-for-bit);
        # position 30: a single tested base.
        sites = pd.DataFrame(
            {
                "contig": ["c1"] * 3,
                "position": [10, 20, 30],
                "gene_id": ["g1", "", ""],
                "min_p_value": [0.001, 0.001455849740, 0.03],
                "group_analyzed": [""] * 3,
            }
        )
        source = pd.DataFrame(
            {
                "contig": ["c1"] * 3,
                "position": [10, 20, 30],
                "gene_id": ["g1 ", np.nan, np.nan],  # real-file quirks
                "A_frequency_p_value_tTest": [0.001, 1.0, 0.03],
                "C_frequency_p_value_tTest": [1.0, 0.001455849740, 1.0],
                "G_frequency_p_value_tTest": [0.001, 1.0, 1.0],
                "T_frequency_p_value_tTest": [1.0, 0.001455849790, 1.0],
            }
        )
        got = candidate_alleles(sites, source, self.TWO, {"": 0.03})
        self.assertEqual(
            got[got.position == 10].allele.tolist(), ["A", "G"]
        )  # complementary pair
        self.assertEqual(
            got[got.position == 20].allele.tolist(), ["C", "T"]
        )  # the old rtol=1e-9 rule kept only C
        self.assertEqual(got[got.position == 30].allele.tolist(), ["A"])
        self.assertEqual(got.n_candidate_alleles.tolist(), [2, 2, 2, 2, 1])

    def test_multi_letter_site_keeps_only_bases_significant_on_their_own(self):
        # Three tested bases: A 1e-4 (the site's min), C 0.002, G 0.5; T untested.
        sites = pd.DataFrame(
            {
                "contig": ["c1"],
                "position": [10],
                "gene_id": ["g1"],
                "min_p_value": [1e-4],
                "group_analyzed": [""],
            }
        )
        source = pd.DataFrame(
            {
                "contig": ["c1"],
                "position": [10],
                "gene_id": ["g1"],
                "A_frequency_p_value_tTest": [1e-4],
                "C_frequency_p_value_tTest": [0.002],
                "G_frequency_p_value_tTest": [0.5],
                "T_frequency_p_value_tTest": [1.0],
            }
        )
        wide = candidate_alleles(
            sites, source, self.TWO, {"": 0.0025}
        )  # line above C's p: A and C
        self.assertEqual(wide.allele.tolist(), ["A", "C"])
        self.assertEqual(wide.n_candidate_alleles.tolist(), [2, 2])
        narrow = candidate_alleles(
            sites, source, self.TWO, {"": 0.0015}
        )  # line between them: A only
        self.assertEqual(narrow.allele.tolist(), ["A"])
        self.assertEqual(narrow.n_candidate_alleles.tolist(), [1])

    def test_single_sample_reads_the_group_suffixed_column_and_that_groups_line(self):
        sites = pd.DataFrame(
            {
                "contig": ["c1"],
                "position": [10],
                "gene_id": ["g1"],
                "min_p_value": [0.004],
                "group_analyzed": ["fat"],
            }
        )
        source = pd.DataFrame(
            {
                "contig": ["c1"],
                "position": [10],
                "gene_id": ["g1"],
                "A_frequency_p_value_tTest_fat": [0.9],
                "C_frequency_p_value_tTest_fat": [0.004],
                "G_frequency_p_value_tTest_fat": [0.9],
                "T_frequency_p_value_tTest_fat": [1.0],
            }
        )
        col = allele_p_column("single_sample", "single_sample_tTest")
        self.assertEqual(
            candidate_alleles(sites, source, col, {"fat": 0.01}).allele.tolist(), ["C"]
        )  # 3 tested: A, G above the line
        with self.assertRaisesRegex(
            ValueError, "fat"
        ):  # a group with no line is a broken call
            candidate_alleles(sites, source, col, {"control": 0.01})

    def test_significant_site_with_no_tested_base_raises(self):
        sites = pd.DataFrame(
            {
                "contig": ["c1"],
                "position": [10],
                "gene_id": ["g1"],
                "min_p_value": [0.001],
                "group_analyzed": [""],
            }
        )
        source = pd.DataFrame(
            {
                "contig": ["c1"],
                "position": [10],
                "gene_id": ["g1"],
                **{f"{b}_frequency_p_value_tTest": [1.0] for b in "ACGT"},
            }
        )
        with self.assertRaisesRegex(ValueError, "no tested base"):
            candidate_alleles(sites, source, self.TWO, {"": 0.05})

    def test_site_missing_from_source_raises_and_names_it(self):
        sites = pd.DataFrame(
            {
                "contig": ["c1", "c1"],
                "position": [10, 99],
                "gene_id": ["g1", "g1"],
                "min_p_value": [1.0, 0.001],
                "group_analyzed": ["", ""],
            }
        )
        source = pd.DataFrame(
            {
                "contig": ["c1"],
                "position": [10],
                "gene_id": ["g1"],
                **{f"{b}_frequency_p_value_tTest": [1.0] for b in "ACGT"},
            }
        )
        with self.assertRaisesRegex(ValueError, "1 of 2 sites.*99"):
            candidate_alleles(sites, source, self.TWO, {"": 0.05})

    def test_gene_id_disagreement_is_a_mismatch(self):
        sites = pd.DataFrame(
            {
                "contig": ["c1"],
                "position": [10],
                "gene_id": ["g1"],
                "min_p_value": [1.0],
                "group_analyzed": [""],
            }
        )
        source = pd.DataFrame(
            {
                "contig": ["c1"],
                "position": [10],
                "gene_id": ["g2"],
                **{f"{b}_frequency_p_value_tTest": [1.0] for b in "ACGT"},
            }
        )
        with self.assertRaises(ValueError):
            candidate_alleles(sites, source, self.TWO, {"": 0.05})

    def test_column_registry_covers_every_family_that_reports_per_base_p(self):
        col = allele_p_column(
            "two_sample_unpaired", "two_sample_unpaired_MannWhitney_abs"
        )
        self.assertEqual(col("G", ""), "G_frequency_p_value_MannWhitney_abs")
        col = allele_p_column("single_sample", "single_sample_Wilcoxon")
        self.assertEqual(col("A", "fat"), "A_frequency_p_value_Wilcoxon_fat")
        col = allele_p_column("lmm", "LMM_abs")
        self.assertEqual(col("T", ""), "T_p_value_LMM_abs")
        col = allele_p_column("lmm_across_time", "LMM")
        self.assertEqual(
            col("C", "control"), "C_p_value_LMM"
        )  # one FILE per group, plain columns


class TestAssessAllele(unittest.TestCase):
    def test_four_evidence_tiers(self):
        counts = np.zeros((5, 4), dtype=np.uint16)
        counts[0] = (25, 0, 5, 0)  # G 5 of 30: present
        counts[1] = (29, 0, 1, 0)  # G 1 of 30: below the bar of 3
        counts[2] = (30, 0, 0, 0)  # G 0 of 30: absent
        counts[3] = (2, 0, 0, 0)  # 2 reads total: not covered at min_cov 5
        counts[4] = (96, 0, 4, 0)  # G 4 of 100: clears bar 4 but under 5 %
        got = assess_allele(
            counts,
            np.array([0, 1, 2, 3, 4]),
            base_idx=2,
            model=MODEL,
            min_freq=0.05,
            min_cov=5,
        )
        self.assertEqual(
            got["allele_status"].tolist(),
            [
                EVIDENCE_PRESENT,
                EVIDENCE_BELOW_DETECTION,
                EVIDENCE_ABSENT,
                EVIDENCE_NOT_COVERED,
                EVIDENCE_BELOW_DETECTION,
            ],
        )
        self.assertEqual(got["allele_reads"].tolist(), [5, 1, 0, 0, 4])
        self.assertEqual(got["total_reads"].tolist(), [30, 30, 30, 2, 100])
        self.assertEqual(
            got["detection_threshold_reads"].tolist(), [3, 3, 3, 2, 4]
        )  # bar is 2 reads at 2x, 4 at 100x
        self.assertEqual(
            got["allele_present"].tolist(), [True, False, False, False, False]
        )
        self.assertTrue(
            np.isnan(got["allele_frequency"][3])
        )  # no frequency without coverage

    def test_min_cov_toggle_off(self):
        counts = np.array([[1, 0, 1, 0]], dtype=np.uint16)  # 2 reads, G 1 of 2
        got = assess_allele(counts, np.array([0]), 2, MODEL, 0.05, min_cov=1)
        self.assertEqual(
            got["allele_status"].tolist(), [EVIDENCE_BELOW_DETECTION]
        )  # covered now, but 1 < bar


class TestLabelOriginInOwnMouse(unittest.TestCase):
    """Later-timepoint rows get a verdict from BOTH their own status and their own
    mouse's earlier-timepoint status; earlier rows stay blank.  One site-allele,
    eight mice covering every case.  Labels are spelled with the comparison's real
    timepoint names (pre/end here), never a fixed t0/t1 vocabulary."""

    CASES = {  # mouse: (pre status or None for no pre sample, end status)
        "m1": ("present", "present"),  # standing_variation
        "m2": (
            "below_detection",
            "present",
        ),  # de_novo_candidate_below_detection_at_pre
        "m3": ("absent", "present"),  # de_novo_candidate
        "m4": ("not_covered", "present"),  # pre_not_covered
        "m5": (None, "present"),  # no_pre_sample
        "m6": (
            "absent",
            "below_detection",
        ),  # allele_below_detection_at_end (pre irrelevant)
        "m7": ("present", "not_covered"),  # end_not_covered
        "m8": ("present", "absent"),  # allele_absent_at_end (pre irrelevant)
    }
    TPS = Timepoints("pre", "end")

    def _long(self, earlier="pre", later="end"):
        site = {
            "mag_id": "MAG_A",
            "contig": "c1",
            "position": 10,
            "gene_id": "g1",
            "allele": "G",
            "test_type": "t",
            "group_analyzed": "",
        }
        rows = []
        for mouse, (t0, t1) in self.CASES.items():
            if t0 is not None:
                rows.append(
                    {
                        **site,
                        "sample_id": f"{mouse}_t0",
                        "subjectID": mouse,
                        "time": earlier,
                        "allele_status": t0,
                    }
                )
            rows.append(
                {
                    **site,
                    "sample_id": f"{mouse}_t1",
                    "subjectID": mouse,
                    "time": later,
                    "allele_status": t1,
                }
            )
        return pd.DataFrame(rows)

    def test_every_label(self):
        got = (
            label_origin_in_own_mouse(self._long(), self.TPS)
            .set_index("sample_id")
            .origin_in_own_mouse
        )
        self.assertEqual(got["m1_t1"], "standing_variation")
        self.assertEqual(got["m2_t1"], "de_novo_candidate_below_detection_at_pre")
        self.assertEqual(got["m3_t1"], "de_novo_candidate")
        self.assertEqual(got["m4_t1"], "pre_not_covered")
        self.assertEqual(got["m5_t1"], "no_pre_sample")
        self.assertEqual(got["m6_t1"], "allele_below_detection_at_end")
        self.assertEqual(got["m8_t1"], "allele_absent_at_end")
        self.assertEqual(got["m7_t1"], "end_not_covered")
        self.assertTrue(got[[s for s in got.index if s.endswith("_t0")]].isna().all())

    def test_labels_follow_the_comparison_not_a_fixed_vocabulary(self):
        # DRiDO-style names: the same table, told it is a 5mo -> 22mo comparison
        got = (
            label_origin_in_own_mouse(
                self._long("5mo", "22mo"), Timepoints("5mo", "22mo")
            )
            .set_index("sample_id")
            .origin_in_own_mouse
        )
        self.assertEqual(got["m2_t1"], "de_novo_candidate_below_detection_at_5mo")
        self.assertEqual(got["m7_t1"], "22mo_not_covered")


class TestSummariseSites(unittest.TestCase):
    """One site, allele G.  pre: m1 present, m2 absent, m3 below, m4 not covered.
    end: m1..m4 all present.  m1+m2 share replicate r1 (fat), m3+m4 replicate r2 (control).
    Reads (allele/total): pre m1 5/30, m2 0/30, m3 1/30, m4 0/2; end all 24/30."""

    TPS = Timepoints("pre", "end")

    def _long(self):
        rows = []
        pre = {
            "m1": EVIDENCE_PRESENT,
            "m2": EVIDENCE_ABSENT,
            "m3": EVIDENCE_BELOW_DETECTION,
            "m4": EVIDENCE_NOT_COVERED,
        }
        meta = {
            "m1": ("r1", "fat"),
            "m2": ("r1", "fat"),
            "m3": ("r2", "control"),
            "m4": ("r2", "control"),
        }
        reads = {"m1": (5, 30), "m2": (0, 30), "m3": (1, 30), "m4": (0, 2)}
        for mouse, ev in pre.items():
            rep, grp = meta[mouse]
            for time, evidence, freq, (ar, tr) in (
                (
                    "pre",
                    ev,
                    {"allele_present": 0.2, "absent": 0.0}.get(ev, 0.02),
                    reads[mouse],
                ),
                ("end", EVIDENCE_PRESENT, 0.8, (24, 30)),
            ):
                rows.append(
                    {
                        "mag_id": "MAG_A",
                        "contig": "c1",
                        "position": 10,
                        "gene_id": "g1",
                        "allele": "G",
                        "test_type": "two_sample_paired_tTest",
                        "group_analyzed": "",
                        "sample_id": f"{mouse}_{time}",
                        "subjectID": mouse,
                        "replicate": rep,
                        "group": grp,
                        "time": time,
                        "allele_status": evidence,
                        "allele_reads": ar,
                        "total_reads": tr,
                        "allele_present": evidence == EVIDENCE_PRESENT,
                        "allele_frequency": (
                            np.nan if evidence == EVIDENCE_NOT_COVERED else freq
                        ),
                    }
                )
        return pd.DataFrame(rows)

    def test_site_verdict_requires_the_allele_to_be_seen_at_the_later_timepoint(self):
        long = label_origin_in_own_mouse(self._long(), self.TPS)
        long.loc[long.time == "end", "allele_status"] = (
            "absent"  # nobody shows it at end
        )
        long.loc[long.time == "end", "allele_present"] = False
        got = summarise_sites(long, self.TPS).iloc[0]
        self.assertEqual(got.origin_any_mouse, "allele_not_present_at_end")
        self.assertEqual(
            int(got.n_pre_samples_allele_present), 1
        )  # the pre fact itself is unchanged

    def test_lean_summary_row(self):
        got = summarise_sites(
            label_origin_in_own_mouse(self._long(), self.TPS), self.TPS
        ).iloc[0]
        self.assertEqual(list(got.index), summary_columns(self.TPS))
        self.assertEqual(
            list(got.index),
            [
                "mag_id",
                "contig",
                "position",
                "gene_id",
                "group_analyzed",
                "allele",
                "n_candidate_alleles",
                "q_value",
                "origin_any_mouse",
                "first_carrier_date",
                "first_carrier_cohort",
                "first_detection_date",
                "first_detection_cohort",
                "n_pre_samples_allele_present",
                "n_pre_samples_covered",
                "n_replicates_with_allele_at_pre",
                "pre_mice_allele_present",
                "end_mice_allele_present",
                "n_mice_standing_variation",
                "n_mice_de_novo_candidate",
                "n_mice_de_novo_candidate_below_detection_at_pre",
                "total_reads_pre",
                "allele_reads_pre",
                "allele_frequency_pre",
                "total_reads_end",
                "allele_reads_end",
                "allele_frequency_end",
            ],
        )
        self.assertEqual(got.origin_any_mouse, ORIGIN_STANDING)  # m1 had it at pre
        # pre: m1 present, m2 absent, m3 below, m4 not covered -> 1 present of 3 covered, 1 replicate (r1)
        self.assertEqual(
            (
                int(got.n_pre_samples_allele_present),
                int(got.n_pre_samples_covered),
                int(got.n_replicates_with_allele_at_pre),
                got.pre_mice_allele_present,
            ),
            (1, 3, 1, "m1"),
        )
        # per-mouse: m1 standing, m2 de novo, m3 de novo (below det.); m4's pre not covered is not a summary column
        self.assertEqual(
            (
                int(got.n_mice_standing_variation),
                int(got.n_mice_de_novo_candidate),
                int(got.n_mice_de_novo_candidate_below_detection_at_pre),
            ),
            (1, 1, 1),
        )

    def test_read_columns_are_unfiltered_sums_over_every_sample(self):
        # The requested bound -- "1000 reads at pre, allele never seen -> frequency < 1/1000" --
        # only holds if EVERY read counts, so the not-covered m4 (2 reads) is in the total.
        got = summarise_sites(
            label_origin_in_own_mouse(self._long(), self.TPS), self.TPS
        ).iloc[0]
        self.assertEqual(
            (int(got.total_reads_pre), int(got.allele_reads_pre)),
            (30 + 30 + 30 + 2, 5 + 0 + 1 + 0),
        )
        self.assertAlmostEqual(got.allele_frequency_pre, 6 / 92)
        self.assertEqual(
            (int(got.total_reads_end), int(got.allele_reads_end)), (120, 96)
        )
        self.assertAlmostEqual(got.allele_frequency_end, 0.8)

    def test_read_frequency_is_nan_when_no_reads_at_all(self):
        long = self._long()
        long.loc[long.time == "pre", ["allele_reads", "total_reads"]] = 0
        got = summarise_sites(label_origin_in_own_mouse(long, self.TPS), self.TPS).iloc[
            0
        ]
        self.assertEqual(int(got.total_reads_pre), 0)
        self.assertTrue(np.isnan(got.allele_frequency_pre))


class TestFrequencyChangeByGroup(unittest.TestCase):
    """The change the within-group test is built on: each mouse later minus
    earlier, mean per replicate, mean over replicates."""

    TPS = Timepoints("pre", "end")

    @staticmethod
    def _long(rows):
        """rows: (mouse, replicate, group, time, allele, allele_reads, total_reads)."""
        return pd.DataFrame(
            [
                {
                    "mag_id": "MAG_A",
                    "contig": "c1",
                    "position": 1,
                    "gene_id": "g1",
                    "test_type": "two_sample_paired_tTest",
                    "group_analyzed": "",
                    "allele": allele,
                    "subjectID": mouse,
                    "replicate": rep,
                    "group": grp,
                    "time": time,
                    "allele_reads": ar,
                    "total_reads": tr,
                }
                for mouse, rep, grp, time, allele, ar, tr in rows
            ]
        )

    def test_seesaw_pair_gets_opposite_changes(self):
        # G 5/30 -> 24/30 and 0/30 -> 21/30 in r1; A is the other side of the seesaw
        long = self._long(
            [
                ("m1", "r1", "fat", "pre", "G", 5, 30),
                ("m1", "r1", "fat", "end", "G", 24, 30),
                ("m2", "r1", "fat", "pre", "G", 0, 30),
                ("m2", "r1", "fat", "end", "G", 21, 30),
                ("m1", "r1", "fat", "pre", "A", 25, 30),
                ("m1", "r1", "fat", "end", "A", 6, 30),
                ("m2", "r1", "fat", "pre", "A", 30, 30),
                ("m2", "r1", "fat", "end", "A", 9, 30),
            ]
        )
        got = frequency_change_by_group(long, self.TPS).set_index("allele")
        # brute force: m1 24/30-5/30 = 0.6333, m2 21/30-0 = 0.7, r1 mean = 0.6667
        self.assertAlmostEqual(
            got.loc["G", "frequency_change"], (19 / 30 + 21 / 30) / 2
        )
        self.assertAlmostEqual(
            got.loc["A", "frequency_change"], -(19 / 30 + 21 / 30) / 2
        )
        self.assertEqual(int(got.loc["G", "n_replicates"]), 1)
        self.assertEqual(int(got.loc["G", "n_mice"]), 2)  # m1, m2 in one cage
        self.assertEqual(got.loc["G", "group"], "fat")

    def test_replicates_weigh_equally(self):
        # r1: two mice (+0.6, +0.2) -> +0.4 ; r2: one mouse (+0.1) -> +0.1
        # group = (0.4 + 0.1) / 2 = 0.25, NOT the mouse mean (0.6+0.2+0.1)/3 = 0.3
        long = self._long(
            [
                ("m1", "r1", "fat", "pre", "G", 0, 10),
                ("m1", "r1", "fat", "end", "G", 6, 10),
                ("m2", "r1", "fat", "pre", "G", 0, 10),
                ("m2", "r1", "fat", "end", "G", 2, 10),
                ("m3", "r2", "fat", "pre", "G", 0, 10),
                ("m3", "r2", "fat", "end", "G", 1, 10),
            ]
        )
        got = frequency_change_by_group(long, self.TPS).iloc[0]
        self.assertAlmostEqual(got.frequency_change, 0.25)
        self.assertEqual(int(got.n_replicates), 2)
        self.assertEqual(int(got.n_mice), 3)  # 3 mice in 2 cages

    def test_unpaired_mouse_and_zero_read_sample_are_ignored(self):
        long = self._long(
            [
                ("m1", "r1", "fat", "pre", "G", 0, 10),
                ("m1", "r1", "fat", "end", "G", 5, 10),  # +0.5, the only pair
                ("m2", "r1", "fat", "end", "G", 10, 10),  # no pre sample: ignored
                ("m3", "r1", "fat", "pre", "G", 0, 0),  # 0 reads at pre: no frequency
                ("m3", "r1", "fat", "end", "G", 10, 10),
            ]
        )
        got = frequency_change_by_group(long, self.TPS)
        self.assertEqual(len(got), 1)
        self.assertAlmostEqual(got.iloc[0].frequency_change, 0.5)
        self.assertEqual(int(got.iloc[0].n_mice), 1)  # only m1 is paired

    def test_site_with_no_paired_mouse_has_no_row(self):
        long = self._long(
            [
                ("m1", "r1", "fat", "pre", "G", 0, 10),
                ("m2", "r1", "fat", "end", "G", 5, 10),
            ]
        )
        self.assertTrue(frequency_change_by_group(long, self.TPS).empty)

    def test_groups_are_separate(self):
        # fat rises, control falls: one row each, never averaged together
        long = self._long(
            [
                ("m1", "r1", "fat", "pre", "G", 0, 10),
                ("m1", "r1", "fat", "end", "G", 5, 10),
                ("m3", "r2", "control", "pre", "G", 5, 10),
                ("m3", "r2", "control", "end", "G", 0, 10),
            ]
        )
        got = frequency_change_by_group(long, self.TPS).set_index("group")
        self.assertAlmostEqual(got.loc["fat", "frequency_change"], 0.5)
        self.assertAlmostEqual(got.loc["control", "frequency_change"], -0.5)

    def test_numeric_group_names_stay_text(self):
        long = self._long(
            [
                ("m1", "r1", "40", "pre", "G", 0, 10),
                ("m1", "r1", "40", "end", "G", 5, 10),
            ]
        )
        got = frequency_change_by_group(long, self.TPS)
        self.assertEqual(got.iloc[0].group, "40")
        self.assertIsInstance(got.iloc[0].group, str)

    def test_two_samples_of_one_mouse_at_one_timepoint_raises(self):
        long = self._long(
            [
                ("m1", "r1", "fat", "pre", "G", 0, 10),
                ("m1", "r1", "fat", "pre", "G", 1, 10),
                ("m1", "r1", "fat", "end", "G", 5, 10),
            ]
        )
        with self.assertRaises(pd.errors.MergeError):
            frequency_change_by_group(long, self.TPS)

    def test_blank_replicate_raises_instead_of_being_lumped(self):
        # a blank cage label must not silently become one shared "blank" cage
        long = self._long(
            [
                ("m1", "r1", "fat", "pre", "G", 0, 10),
                ("m1", "r1", "fat", "end", "G", 5, 10),
            ]
        )
        long["replicate"] = np.nan
        with self.assertRaisesRegex(ValueError, "replicate"):
            frequency_change_by_group(long, self.TPS)


class TestSummariseRising(unittest.TestCase):
    """Site c1:1, G and A tie (seesaw).  fat m1, m2 (r1); control m3, m4 (r2).
    G reads pre/end: m1 5/30 -> 24/30, m2 0/30 -> 21/30, m3 1/30 -> 2/30,
    m4 0/2 -> 0/30.  A is the rest of each pile.
    Changes: G fat +0.667, G control +0.0167 (m3 +0.033, m4 0); A the negatives."""

    TPS = Timepoints("pre", "end")
    # mouse -> (replicate, group, G at pre, G at end, total at pre, total at end)
    MICE = {
        "m1": ("r1", "fat", 5, 24, 30, 30),
        "m2": ("r1", "fat", 0, 21, 30, 30),
        "m3": ("r2", "control", 1, 2, 30, 30),
        "m4": ("r2", "control", 0, 0, 2, 30),
    }
    # (allele, mouse, time) -> status, by the presence rule (min_cov 5, bar 3, 5%)
    STATUS = {
        ("G", "m1", "pre"): "present",
        ("G", "m1", "end"): "present",
        ("G", "m2", "pre"): "absent",
        ("G", "m2", "end"): "present",
        ("G", "m3", "pre"): "below_detection",
        ("G", "m3", "end"): "below_detection",
        ("G", "m4", "pre"): "not_covered",
        ("G", "m4", "end"): "absent",
        ("A", "m1", "pre"): "present",
        ("A", "m1", "end"): "present",
        ("A", "m2", "pre"): "present",
        ("A", "m2", "end"): "present",
        ("A", "m3", "pre"): "present",
        ("A", "m3", "end"): "present",
        ("A", "m4", "pre"): "not_covered",
        ("A", "m4", "end"): "present",
    }

    def _long(self, group_analyzed=""):
        rows = []
        for mouse, (rep, grp, g_pre, g_end, n_pre, n_end) in self.MICE.items():
            for time, g, n in (("pre", g_pre, n_pre), ("end", g_end, n_end)):
                for allele, reads in (("G", g), ("A", n - g)):
                    status = self.STATUS[(allele, mouse, time)]
                    rows.append(
                        {
                            "mag_id": "MAG_A",
                            "contig": "c1",
                            "position": 1,
                            "gene_id": "g1",
                            "test_type": "two_sample_paired_tTest",
                            "group_analyzed": group_analyzed,
                            "allele": allele,
                            "q_value": 0.01,
                            "n_candidate_alleles": 2,
                            "sample_id": f"{mouse}_{time}",
                            "subjectID": mouse,
                            "replicate": rep,
                            "group": grp,
                            "time": time,
                            "allele_reads": reads,
                            "total_reads": n,
                            "allele_status": status,
                            "allele_present": status == "present",
                        }
                    )
        return label_origin_in_own_mouse(pd.DataFrame(rows), self.TPS)

    def _rising(self, long):
        change = frequency_change_by_group(long, self.TPS)
        return summarise_rising(long, change, self.TPS)

    def test_fallen_seesaw_partner_is_left_out(self):
        got = self._rising(self._long())
        self.assertEqual(set(got.allele), {"G"})  # A fell in both groups

    def test_rose_in_both_groups_gives_one_row_per_group_with_its_own_end_check(
        self,
    ):
        got = self._rising(self._long()).set_index("rose_in")
        self.assertEqual(sorted(got.index), ["control", "fat"])
        fat, control = got.loc["fat"], got.loc["control"]
        # fat: m1 had G at pre (baseline = every pre sample) -> standing
        self.assertEqual(fat.origin_any_mouse, "standing_variation")
        self.assertAlmostEqual(fat.frequency_change, (19 / 30 + 21 / 30) / 2)
        self.assertEqual((int(fat.n_replicates), int(fat.n_mice)), (1, 2))
        self.assertEqual(
            (int(fat.n_mice_standing_variation), int(fat.n_mice_de_novo_candidate)),
            (1, 1),
        )
        # baseline reads: every pre sample; later reads: fat's own end samples only
        self.assertEqual((int(fat.total_reads_pre), int(fat.allele_reads_pre)), (92, 6))
        self.assertEqual(
            (int(fat.total_reads_end), int(fat.allele_reads_end)), (60, 45)
        )
        # control: G rose a little, but no CONTROL end mouse has it present (m3 2/30
        # is under the bar, m4 0/30) -> nothing to explain, judged on control alone
        self.assertEqual(control.origin_any_mouse, "allele_not_present_at_end")
        self.assertEqual(
            (int(control.total_reads_end), int(control.allele_reads_end)), (60, 2)
        )

    def test_full_summary_is_unchanged_by_the_refactor(self):
        # the pooled summary still reports both alleles, G standing via fat
        summ = summarise_sites(self._long(), self.TPS).set_index("allele")
        self.assertEqual(sorted(summ.index), ["A", "G"])
        self.assertEqual(summ.loc["G", "origin_any_mouse"], "standing_variation")

    def test_within_group_site_considers_only_its_own_group(self):
        # a site significant within fat: control's small rise is not a candidate
        got = self._rising(self._long(group_analyzed="fat"))
        self.assertEqual(got.rose_in.tolist(), ["fat"])

    def test_flat_or_unpaired_is_not_rising(self):
        long = self._long()
        flat = long.copy()
        # make G identical at pre and end in every mouse: change 0 -> no G row
        for mouse, (_, _, g_pre, _, n_pre, _) in self.MICE.items():
            end_g = (
                (flat.subjectID == mouse) & (flat.time == "end") & (flat.allele == "G")
            )
            flat.loc[end_g, ["allele_reads", "total_reads"]] = (g_pre, n_pre)
        self.assertNotIn("G", set(self._rising(flat).allele))
        # drop every pre sample: no mouse is paired -> no row at all
        self.assertTrue(self._rising(long[long.time == "end"]).empty)

    def test_unknown_change_is_counted_in_the_log(self):
        # a within-fat site whose fat mice are unpaired (no pre sample) while a
        # control mouse is paired: fat's change is unknown, so no rising row, and
        # the log must count that (site, allele, group) candidate
        long = self._long(group_analyzed="fat")
        long = long[~((long.group == "fat") & (long.time == "pre"))]
        with self.assertLogs(
            "alleleflux.scripts.analysis.baseline_presence", level="INFO"
        ) as logs:
            got = self._rising(long)
        self.assertTrue(got.empty)
        self.assertTrue(
            any(
                "2 (site, allele, group) candidates had no mouse" in m
                for m in logs.output
            ),
            logs.output,
        )  # G-fat and A-fat

    def test_columns_and_empty_input(self):
        got = self._rising(self._long().iloc[0:0])
        self.assertTrue(got.empty)
        self.assertEqual(list(got.columns), rising_columns(self.TPS))
        cols = rising_columns(self.TPS)
        start = cols.index("q_value") + 1
        self.assertEqual(cols[start : start + len(RISING_COLUMNS)], RISING_COLUMNS)


# ---------------------------------------------------------------------------
# End-to-end
# ---------------------------------------------------------------------------
def _write_profile(directory, sample, mag, rows):
    """rows: (position, A, C, G, T) on contig c1, gene g1.

    Writes the sample's rows into the run's allele-frequency CACHE (the single
    per-site source after QC) at ``{run}/allele_freq_cache/{tp}/{mag}_{tp}_allele_frequency.parquet``;
    ``directory`` is the run's ``profiles`` folder as before, the timepoint is the
    sample id's suffix, and calling again for the same sample replaces its rows.
    """
    from alleleflux.scripts.utilities.cache_io import cache_path

    run = os.path.dirname(directory.rstrip("/"))
    mouse, tp = sample.rsplit("_", 1)
    path = cache_path(os.path.join(run, "allele_freq_cache"), mag, tp)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    recs = []
    for p, a, c, g, tt in rows:
        tot = a + c + g + tt
        recs.append(
            dict(
                contig="c1",
                position=p,
                gene_id="g1",
                total_coverage=tot,
                A=a,
                C=c,
                G=g,
                T=tt,
                sample_id=sample,
                group="x",
                subjectID=mouse,
                replicate=mouse,
                A_frequency=a / tot,
                T_frequency=tt / tot,
                G_frequency=g / tot,
                C_frequency=c / tot,
                time=tp,
            )
        )
    new = pd.DataFrame(recs)
    if os.path.exists(path):
        old = pd.read_parquet(path)
        for col in ("sample_id", "group", "subjectID", "replicate"):
            old[col] = old[col].astype(str)
        new = pd.concat([old[old["sample_id"] != sample], new], ignore_index=True)
    for col in ("sample_id", "group", "subjectID", "replicate"):
        new[col] = new[col].astype("category")
    new.to_parquet(path, index=False)


def _drop_sample_from_cache(directory, sample, mag):
    """Remove one sample's rows from the cache (it failed QC for this MAG)."""
    run = os.path.dirname(directory.rstrip("/"))
    tp = sample.rsplit("_", 1)[1]
    path = os.path.join(
        run, "allele_freq_cache", tp, f"{mag}_{tp}_allele_frequency.parquet"
    )
    df = pd.read_parquet(path)
    df["sample_id"] = df["sample_id"].astype(str)
    df = df[df["sample_id"] != sample]
    for col in ("sample_id", "group", "subjectID", "replicate"):
        df[col] = df[col].astype("category")
    df.to_parquet(path, index=False)


class TestBaselinePresenceCLI(unittest.TestCase):
    """8 samples: mice m1,m2 (fat, replicate r1) and m3,m4 (control, r2), each at pre and end.
    One significant site, c1:10, where G rose at end in the fat mice.  t0 G reads:
    m1 5/30 present, m2 0/30 absent, m3 1/30 below detection, m4 2 reads total (not covered).
    """

    MAG = "MAG_A"

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        run = self.run = os.path.join(self.tmp, "run")
        prof = os.path.join(run, "profiles")
        base = [(p, 30, 0, 0, 0) for p in range(6)]

        def with_site(a, c, g, t):
            return (
                [(p, 30, 0, 0, 0) if p != 3 else (3, a, c, g, t) for p in range(6)]
                if False
                else [(p, 30, 0, 0, 0) for p in range(6) if p != 1] + [(1, a, c, g, t)]
            )

        pre = {
            "m1": (25, 0, 5, 0),
            "m2": (30, 0, 0, 0),
            "m3": (29, 0, 1, 0),
            "m4": (2, 0, 0, 0),
        }
        end = {
            "m1": (6, 0, 24, 0),
            "m2": (9, 0, 21, 0),
            "m3": (28, 0, 2, 0),
            "m4": (30, 0, 0, 0),
        }
        for mouse in pre:
            _write_profile(prof, f"{mouse}_pre", self.MAG, with_site(*pre[mouse]))
            _write_profile(prof, f"{mouse}_end", self.MAG, with_site(*end[mouse]))
        meta = pd.DataFrame(
            [
                {
                    "sample_id": f"{m}_{t}",
                    "subjectID": m,
                    "group": g,
                    "time": t,
                    "replicate": r,
                    "bam_path": "x",
                }
                for m, g, r in (
                    ("m1", "fat", "r1"),
                    ("m2", "fat", "r1"),
                    ("m3", "control", "r2"),
                    ("m4", "control", "r2"),
                )
                for t in ("pre", "end")
            ]
        )
        self.metadata = os.path.join(self.tmp, "metadata.tsv")
        meta.to_csv(self.metadata, sep="\t", index=False)
        # p_value_summary + the per-MAG test file the summary points at (site c1:1, G and A tie)
        sd = os.path.join(run, "p_value_summary", "pre_end-fat_control")
        os.makedirs(sd)
        pd.DataFrame(
            {
                "period": ["pre_end"],
                "mag_id": [self.MAG],
                "contig": ["c1"],
                "position": [1],
                "gene_id": ["g1 "],
                "test_type": ["two_sample_paired_tTest"],
                "min_p_value": [0.001],
                "source_file": [f"{self.MAG}_two_sample_paired.tsv.gz"],
                "q_value": [0.01],
            }
        ).to_csv(
            os.path.join(sd, "p_value_summary_two_sample_paired_pre_end.tsv"),
            sep="\t",
            index=False,
        )
        td = os.path.join(
            run, "significance_tests", "two_sample_paired_pre_end-fat_control"
        )
        os.makedirs(td)
        pd.DataFrame(
            {
                "contig": ["c1"],
                "gene_id": ["g1"],
                "position": [1],
                "A_frequency_p_value_tTest": [0.001],
                "C_frequency_p_value_tTest": [1.0],
                "G_frequency_p_value_tTest": [0.001],
                "T_frequency_p_value_tTest": [1.0],
                "num_pairs": [4],
            }
        ).to_csv(
            os.path.join(td, f"{self.MAG}_two_sample_paired.tsv.gz"),
            sep="\t",
            index=False,
            compression="gzip",
        )
        self.fasta = os.path.join(self.tmp, "ref.fa")
        open(self.fasta, "w").write(">c1\nAAAAAA\n")
        open(self.fasta + ".fai", "w").write("c1\t6\t4\t6\t7\n")
        self.mapping = os.path.join(self.tmp, "map.tsv")
        pd.DataFrame({"mag_id": [self.MAG], "contig_id": ["c1"]}).to_csv(
            self.mapping, sep="\t", index=False
        )
        self.out = os.path.join(self.tmp, "out")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self, *extra):
        return subprocess.run(
            [
                "alleleflux-baseline-presence",
                "--run_dir",
                self.run,
                "--comparison",
                "pre_end-fat_control",
                "--summary",
                "two_sample_paired",
                "--test_type",
                "two_sample_paired_tTest",
                "--cache_dir",
                os.path.join(self.run, "allele_freq_cache"),
                "--metadata",
                self.metadata,
                "--output_dir",
                self.out,
                "--cpus",
                "2",
                *extra,
            ],
            capture_output=True,
            text=True,
        )

    def test_site_filter_and_profile_flags_are_gone(self):
        """The cache is already site-filtered and QC-gated: the profile, FASTA, mapping and
        site-filter flags no longer exist on this command."""
        for flag in (
            ["--site_depth", "5"],
            ["--disable_site_depth_filter"],
            ["--profiles_dir", self.tmp],
        ):
            done = self._run(*flag)
            self.assertEqual(done.returncode, 2, flag)

    def test_long_table_one_row_per_site_allele_sample(self):
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        long = pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_two_sample_paired_tTest_baseline_presence.tsv.gz",
            ),
            sep="\t",
        )
        self.assertEqual(len(long), 2 * 8)  # 2 candidate alleles x 8 samples
        g = long[(long.allele == "G")].set_index("sample_id")
        self.assertEqual(g.loc["m1_pre"].allele_status, "present")
        self.assertEqual(int(g.loc["m1_pre"].allele_reads), 5)
        self.assertEqual(g.loc["m2_pre"].allele_status, "absent")
        self.assertEqual(g.loc["m3_pre"].allele_status, "below_detection")
        self.assertEqual(int(g.loc["m3_pre"].detection_threshold_reads), 3)
        self.assertEqual(g.loc["m4_pre"].allele_status, "not_covered")
        self.assertTrue(pd.isna(g.loc["m4_pre"].allele_frequency))
        row = g.loc["m1_pre"]
        self.assertEqual(
            (
                row.mag_id,
                row.contig,
                int(row.position),
                row.gene_id,
                row.subjectID,
                row.replicate,
                row.group,
                row.time,
                int(row.n_candidate_alleles),
                int(row.min_cov),
            ),
            (self.MAG, "c1", 1, "g1", "m1", "r1", "fat", "pre", 2, 5),
        )
        self.assertNotIn(
            "timepoint_role", long.columns
        )  # the real label is the only timepoint column
        self.assertEqual(sorted(long.allele.unique()), ["A", "G"])
        self.assertEqual(g.loc["m1_end"].origin_in_own_mouse, "standing_variation")
        self.assertEqual(g.loc["m2_end"].origin_in_own_mouse, "de_novo_candidate")
        self.assertEqual(
            g.loc["m3_end"].origin_in_own_mouse, "allele_below_detection_at_end"
        )  # G 2/30 at end
        self.assertEqual(
            g.loc["m4_end"].origin_in_own_mouse, "allele_absent_at_end"
        )  # G 0/30 at end
        self.assertTrue(pd.isna(g.loc["m1_pre"].origin_in_own_mouse))
        self.assertTrue(long.strain_background.isna().all())  # no --turnover_dir given

    def test_summary_counts(self):
        self._run()
        summ = pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_two_sample_paired_tTest_baseline_presence_summary.tsv",
            ),
            sep="\t",
        )
        g = summ[summ.allele == "G"].iloc[0]
        self.assertEqual(
            (
                int(g.n_pre_samples_allele_present),
                int(g.n_pre_samples_covered),
                int(g.n_replicates_with_allele_at_pre),
            ),
            (1, 3, 1),
        )
        self.assertEqual(g.pre_mice_allele_present, "m1")
        self.assertEqual(g.origin_any_mouse, "standing_variation")
        self.assertEqual(
            (int(g.n_mice_standing_variation), int(g.n_mice_de_novo_candidate)), (1, 1)
        )  # m1; m2
        # reads pooled over all four mice, thin m4 included: pre 5+0+1+0 of 30+30+30+2; end 24+21+2+0 of 120
        self.assertEqual((int(g.total_reads_pre), int(g.allele_reads_pre)), (92, 6))
        self.assertEqual((int(g.total_reads_end), int(g.allele_reads_end)), (120, 47))
        self.assertAlmostEqual(g.allele_frequency_end, 47 / 120)
        self.assertEqual(
            len(summ.columns), 27
        )  # 22 + first_carrier_date, first_detection_date, first_detection_cohort

    def test_strain_background_joins_when_given(self):
        tdir = os.path.join(self.tmp, "turnover")
        os.makedirs(tdir)
        pd.DataFrame(
            {
                "MAG_ID": [self.MAG] * 2,
                "subjectID": ["m1", "m2"],
                "transition": ["pre_end"] * 2,
                "background": ["stable", "strain_replacement"],
            }
        ).to_csv(
            os.path.join(tdir, f"{self.MAG}_strain_turnover.tsv"), sep="\t", index=False
        )
        done = self._run("--turnover_dir", tdir)
        self.assertEqual(done.returncode, 0, done.stderr)
        long = pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_two_sample_paired_tTest_baseline_presence.tsv.gz",
            ),
            sep="\t",
        )
        bg = long.drop_duplicates("subjectID").set_index("subjectID").strain_background
        self.assertEqual((bg["m1"], bg["m2"]), ("stable", "strain_replacement"))
        m2_end = long[(long.sample_id == "m2_end") & (long.allele == "G")].iloc[0]
        self.assertEqual(
            m2_end.origin_in_own_mouse, "de_novo_candidate"
        )  # label does NOT depend on the strain column
        self.assertTrue(pd.isna(bg["m3"]))  # mouse not in the turnover table

    def test_sample_absent_from_the_cache_gets_no_rows(self):
        """A sample with no rows in the MAG's cache did not pass QC for that MAG (spec 2.1):
        it is not assessed and appears nowhere in the long table."""
        _drop_sample_from_cache(os.path.join(self.run, "profiles"), "m4_end", self.MAG)
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        long = pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_two_sample_paired_tTest_baseline_presence.tsv.gz",
            ),
            sep="\t",
        )
        self.assertNotIn("m4_end", set(long["sample_id"]))
        self.assertEqual(len(long), 2 * 7)

    def test_no_significant_sites_gives_header_only_outputs(self):
        done = self._run("--threshold", "0.0001")  # q=0.01 site no longer passes
        self.assertEqual(done.returncode, 0, done.stderr)
        long = pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_two_sample_paired_tTest_baseline_presence.tsv.gz",
            ),
            sep="\t",
        )
        self.assertEqual(len(long), 0)
        self.assertIn("allele_status", long.columns)
        rise = pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_two_sample_paired_tTest_baseline_presence_rising.tsv",
            ),
            sep="\t",
        )
        self.assertEqual(len(rise), 0)
        self.assertIn("rose_in", rise.columns)  # header-only, never a missing file

    def test_missing_summary_for_the_chosen_family_fails_loud(self):
        done = self._run(
            "--summary",
            "two_sample_unpaired",
            "--test_type",
            "two_sample_unpaired_tTest",
        )
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("two_sample_unpaired", done.stderr)

    def test_own_group_with_a_between_group_test_fails_before_writing(self):
        # the same run succeeds with the default scope; own_group must stop it
        # outright instead of quietly judging the site on both groups
        done = self._run("--baseline_scope", "own_group")
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("own_group", done.stderr)
        self.assertIn("two_sample_paired", done.stderr)
        self.assertEqual(os.listdir(self.out) if os.path.isdir(self.out) else [], [])

    def test_rising_table_keeps_the_risen_allele_per_group(self):
        # G rose in fat (+0.667) and a little in control (+0.017); A fell in both
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        rise = pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_two_sample_paired_tTest_baseline_presence_rising.tsv",
            ),
            sep="\t",
        ).set_index("rose_in")
        self.assertEqual(set(rise.allele), {"G"})  # A fell: left out
        self.assertEqual(sorted(rise.index), ["control", "fat"])  # G rose in both
        fat, control = rise.loc["fat"], rise.loc["control"]
        self.assertEqual(fat.origin_any_mouse, "standing_variation")
        self.assertAlmostEqual(fat.frequency_change, (19 / 30 + 21 / 30) / 2)
        self.assertEqual((int(fat.n_replicates), int(fat.n_mice)), (1, 2))
        # judged on control's end mice only: m3 2/30 under the bar, m4 0/30
        self.assertEqual(control.origin_any_mouse, "allele_not_present_at_end")
        # the carriers (later-timepoint mice with the allele present), for the
        # carrier-level strain gate downstream: fat m1 and m2, control nobody
        self.assertEqual(fat.end_mice_allele_present, "m1,m2")
        self.assertTrue(pd.isna(control.end_mice_allele_present))
        self.assertIn("rose in at least one group", done.stderr)  # counts logged

    # ---- Task 5 (spec 2.2, 2.5): optional metadata columns ride into the long table.
    def _long(self):
        return pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_two_sample_paired_tTest_baseline_presence.tsv.gz",
            ),
            sep="\t",
            dtype=str,
        )

    def test_optional_metadata_columns_ride_along(self):
        """collection_date (ISO YYYY-MM-DD) and cohort copied from the metadata onto
        every long row of that sample, verbatim."""
        meta = pd.read_csv(self.metadata, sep="\t", dtype=str)
        # 8 samples -> dates 2016-01-01 .. 2016-01-08 in metadata order (m1_pre first)
        meta["collection_date"] = [f"2016-01-0{i + 1}" for i in range(len(meta))]
        meta["cohort"] = "batch1"
        meta.to_csv(self.metadata, sep="\t", index=False)
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        long = self._long()
        self.assertEqual(
            long.loc[long.sample_id == "m1_pre", "collection_date"].iloc[0],
            "2016-01-01",
        )
        self.assertEqual(
            long.loc[long.sample_id == "m4_end", "collection_date"].iloc[0],
            "2016-01-08",
        )
        self.assertEqual(set(long["cohort"]), {"batch1"})

    def test_missing_optional_columns_are_blank(self):
        """Metadata without the two columns: they exist in the long table, all blank."""
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        long = self._long()
        self.assertIn("collection_date", long.columns)
        self.assertIn("cohort", long.columns)
        self.assertTrue(long["collection_date"].isna().all())
        self.assertTrue(long["cohort"].isna().all())
        # placed right after `time`, so the sample-context block stays together
        cols = list(long.columns)
        self.assertEqual(
            cols[cols.index("time") : cols.index("time") + 3],
            ["time", "collection_date", "cohort"],
        )

    def test_non_iso_date_is_refused_before_any_output(self):
        """A date in any format but YYYY-MM-DD (here DD/MM/YYYY) is a malformed
        metadata file: the job stops and names the format, no guessing, no output."""
        meta = pd.read_csv(self.metadata, sep="\t", dtype=str)
        meta["collection_date"] = "10/01/2017"
        meta.loc[meta.sample_id == "m1_pre", "collection_date"] = (
            "2017-01-10"  # one good
        )
        meta.to_csv(self.metadata, sep="\t", index=False)
        done = self._run()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("YYYY-MM-DD", done.stderr)
        self.assertIn("10/01/2017", done.stderr)
        self.assertEqual(os.listdir(self.out) if os.path.isdir(self.out) else [], [])

    def test_blank_dates_are_allowed(self):
        """An empty collection_date cell is "unknown", not malformed: the job runs and
        the row's date is blank while the dated rows keep theirs."""
        meta = pd.read_csv(self.metadata, sep="\t", dtype=str)
        meta["collection_date"] = ""
        meta.loc[meta.sample_id == "m1_pre", "collection_date"] = "2017-01-10"
        meta.to_csv(self.metadata, sep="\t", index=False)
        # explicit scope: under the default (cohort) a partly dated sheet is refused
        done = self._run("--baseline_scope", "comparison")
        self.assertEqual(done.returncode, 0, done.stderr)
        long = self._long()
        self.assertEqual(
            long.loc[long.sample_id == "m1_pre", "collection_date"].iloc[0],
            "2017-01-10",
        )
        self.assertTrue(
            long.loc[long.sample_id == "m2_pre", "collection_date"].isna().all()
        )

    # ---- Task 6 (spec 2.2): baseline scopes individual and cohort.
    def _summary(self, **read_kw):
        return pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_two_sample_paired_tTest_baseline_presence_summary.tsv",
            ),
            sep="\t",
            **read_kw,
        )

    def test_individual_scope_uses_carriers_own_baseline_only(self):
        """Carriers at end are m1 and m2 (m3's 2/30 is under the bar).  Under
        `individual` only THEIR pre samples are the baseline, so m3's pre with its
        one stray G read is not consulted.  With m1's pre G reads removed, both
        carriers' pre samples are covered and absent: de_novo_candidate, judged on
        2 covered baseline samples (the default scope would say below_detection
        from m3's pre, on 3)."""
        prof = os.path.join(self.run, "profiles")
        _write_profile(prof, "m1_pre", self.MAG, [(p, 30, 0, 0, 0) for p in range(6)])
        done = self._run("--baseline_scope", "individual")
        self.assertEqual(done.returncode, 0, done.stderr)
        g = self._summary().set_index("allele").loc["G"]
        self.assertEqual(g["origin_any_mouse"], "de_novo_candidate")
        self.assertEqual(int(g["n_pre_samples_covered"]), 2)
        self.assertEqual(
            (
                g["pre_mice_allele_present"]
                if isinstance(g["pre_mice_allele_present"], str)
                else ""
            ),
            "",
        )
        # and the default scope on the same data still consults m3's pre
        done = self._run("--baseline_scope", "comparison")
        self.assertEqual(done.returncode, 0, done.stderr)
        g = self._summary().set_index("allele").loc["G"]
        self.assertEqual(
            g["origin_any_mouse"], "de_novo_candidate_below_detection_at_pre"
        )
        self.assertEqual(int(g["n_pre_samples_covered"]), 3)

    def _add_cohort_samples(self, dates=None):
        """Two mice of a third group 'other' at a timepoint 'mid' (replicate r3),
        G present in both (20/30).  With ``dates`` (sample -> 'YYYY-MM-DD' or ''),
        every sample is dated from it and given a cohort label; with None the
        sheet gets NO collection_date / cohort column at all."""
        prof = os.path.join(self.run, "profiles")
        meta = pd.read_csv(self.metadata, sep="\t", dtype=str)
        extra = []
        for mouse in ("m5", "m6"):
            s = f"{mouse}_mid"
            _write_profile(
                prof,
                s,
                self.MAG,
                [(p, 30, 0, 0, 0) for p in range(6) if p != 1] + [(1, 10, 0, 20, 0)],
            )
            extra.append(
                {
                    "sample_id": s,
                    "subjectID": mouse,
                    "group": "other",
                    "time": "mid",
                    "replicate": "r3",
                    "bam_path": "x",
                }
            )
        meta = pd.concat([meta, pd.DataFrame(extra)], ignore_index=True)
        if dates is not None:
            meta["collection_date"] = [dates[s] for s in meta["sample_id"]]
            meta["cohort"] = [
                "wave1" if s.endswith("_mid") else "wave2" for s in meta["sample_id"]
            ]
        meta.to_csv(self.metadata, sep="\t", index=False)

    def test_cohort_scope_finds_the_allele_in_another_group_before_the_first_carrier(
        self,
    ):
        """The fat carriers' end samples are dated 2017-06-10.  A third-group sample
        dated 2016-09-01 already carries G: standing under `cohort` (the allele was
        in the building), de novo under `comparison`.  The 2018 sample is after the
        first carrier and is not baseline."""
        dates = {f"{m}_pre": "2017-01-10" for m in ("m1", "m2", "m3", "m4")}
        dates.update({f"{m}_end": "2017-06-10" for m in ("m1", "m2", "m3", "m4")})
        dates.update({"m5_mid": "2016-09-01", "m6_mid": "2018-01-01"})
        self._add_cohort_samples(dates)
        prof = os.path.join(self.run, "profiles")
        _write_profile(prof, "m1_pre", self.MAG, [(p, 30, 0, 0, 0) for p in range(6)])
        done = self._run("--baseline_scope", "cohort")
        self.assertEqual(done.returncode, 0, done.stderr)
        g = (
            self._summary(dtype={"first_detection_cohort": str})
            .set_index("allele")
            .loc["G"]
        )
        self.assertEqual(g["origin_any_mouse"], "standing_variation")
        self.assertEqual(g["first_detection_date"], "2016-09-01")
        self.assertEqual(g["first_detection_cohort"], "wave1")
        self.assertEqual(g["first_carrier_date"], "2017-06-10")
        # covered baseline: m1, m2, m3 pre (m4 pre has 2 reads: not covered) + m5_mid
        self.assertEqual(int(g["n_pre_samples_covered"]), 4)
        self.assertEqual(g["pre_mice_allele_present"], "m5")
        # the comparison-level numbers stay the comparison's: end reads from fat+control only
        self.assertEqual(int(g["total_reads_end"]), 120)
        rising = pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_two_sample_paired_tTest_baseline_presence_rising.tsv",
            ),
            sep="\t",
        ).set_index("rose_in")
        self.assertEqual(sorted(rising.index), ["control", "fat"])  # no row for 'other'
        self.assertEqual(rising.loc["fat", "origin_any_mouse"], "standing_variation")
        long = self._long()
        self.assertEqual(
            set(long.loc[long.group == "other", "is_comparison_sample"]), {"False"}
        )
        self.assertEqual(
            set(long.loc[long.group != "other", "is_comparison_sample"]), {"True"}
        )
        # same data, default scope: the third group is never loaded
        done = self._run("--baseline_scope", "comparison")
        self.assertEqual(done.returncode, 0, done.stderr)
        g = self._summary().set_index("allele").loc["G"]
        self.assertTrue(g["origin_any_mouse"].startswith("de_novo_candidate"))
        self.assertNotIn("other", set(self._long()["group"]))

    def test_cohort_is_the_default_scope(self):
        """No --baseline_scope on the command line means cohort (user decision
        2026-10-22): on the dated sheet the default run finds G in the third-group
        sample from 2016 and calls it standing, exactly as --baseline_scope cohort."""
        dates = {f"{m}_pre": "2017-01-10" for m in ("m1", "m2", "m3", "m4")}
        dates.update({f"{m}_end": "2017-06-10" for m in ("m1", "m2", "m3", "m4")})
        dates.update({"m5_mid": "2016-09-01", "m6_mid": "2018-01-01"})
        self._add_cohort_samples(dates)
        prof = os.path.join(self.run, "profiles")
        _write_profile(prof, "m1_pre", self.MAG, [(p, 30, 0, 0, 0) for p in range(6)])
        done = self._run()  # no scope flag
        self.assertEqual(done.returncode, 0, done.stderr)
        g = self._summary().set_index("allele").loc["G"]
        self.assertEqual(g["origin_any_mouse"], "standing_variation")
        self.assertEqual(g["pre_mice_allele_present"], "m5")
        self.assertEqual(BASELINE_SCOPE_DEFAULT, "cohort")

    def test_no_carrier_row_is_cut_at_the_groups_sampling_day_and_carrier_wave_is_reported(
        self,
    ):
        """Control rose (1 -> 2 reads) but no control mouse clears the bar: no carrier.
        Its baseline is cut at the day control was sampled at the later timepoint
        (2017-06-10), the same line as the fat row, so both rows describe the same
        4 covered baseline samples with m5 present; its first_carrier_* stay blank.
        The fat row reports the carriers' wave (wave2) beside the first detection's
        (wave1).  Dates are never replaced by label order on a dated sheet."""
        dates = {f"{m}_pre": "2017-01-10" for m in ("m1", "m2", "m3", "m4")}
        dates.update({f"{m}_end": "2017-06-10" for m in ("m1", "m2", "m3", "m4")})
        dates.update({"m5_mid": "2016-09-01", "m6_mid": "2018-01-01"})
        self._add_cohort_samples(dates)
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        rising = pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_two_sample_paired_tTest_baseline_presence_rising.tsv",
            ),
            sep="\t",
            dtype=str,
        ).set_index("rose_in")
        fat, control = rising.loc["fat"], rising.loc["control"]
        self.assertEqual(int(fat["n_pre_samples_covered"]), 4)
        self.assertEqual(int(control["n_pre_samples_covered"]), 4)
        self.assertEqual(
            control["pre_mice_allele_present"], "m1,m5"
        )  # m1 pre has 5 G reads here
        self.assertEqual(control["origin_any_mouse"], "allele_not_present_at_end")
        self.assertTrue(pd.isna(control["first_carrier_date"]))
        self.assertTrue(pd.isna(control["first_carrier_cohort"]))
        self.assertEqual(control["first_detection_date"], "2016-09-01")
        self.assertEqual(fat["first_carrier_date"], "2017-06-10")
        self.assertEqual(fat["first_carrier_cohort"], "wave2")
        self.assertEqual(fat["first_detection_cohort"], "wave1")

    def test_cohort_partly_dated_sheet_is_refused(self):
        """A collection_date column with blanks under `cohort` is a data error, not
        a design: ordering the baseline by label while some samples ARE dated
        would be silent damage.  Stop before any output, name every undated sample."""
        dates = {
            f"{m}_{t}": "" for m in ("m1", "m2", "m3", "m4") for t in ("pre", "end")
        }
        dates.update({"m5_mid": "2016-09-01", "m6_mid": "2018-01-01"})
        self._add_cohort_samples(dates)
        done = self._run("--baseline_scope", "cohort")
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("collection_date", done.stderr)
        for s in ("m1_pre", "m4_end"):
            self.assertIn(s, done.stderr)
        self.assertEqual(os.listdir(self.out) if os.path.isdir(self.out) else [], [])
        # the same sheet is fine under every other scope: dates are only carried along
        done = self._run("--baseline_scope", "comparison")
        self.assertEqual(done.returncode, 0, done.stderr)

    def test_cohort_without_a_date_column_falls_back_to_label_order_with_warning(self):
        """No collection_date column anywhere: a study that never recorded dates.
        'before' means an earlier timepoint LABEL in the sheet's order (pre, end,
        mid): 'mid' sorts after 'end', neither extra sample is baseline, and the
        run says so once."""
        self._add_cohort_samples(dates=None)
        done = self._run("--baseline_scope", "cohort")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("timepoint order", done.stderr)
        g = self._summary().set_index("allele").loc["G"]
        self.assertEqual(int(g["n_pre_samples_covered"]), 3)  # m1, m2, m3 pre only
        self.assertEqual(g["origin_any_mouse"], "standing_variation")  # m1 pre has G
        self.assertTrue(pd.isna(g["first_carrier_date"]))

    def test_cohort_skips_timepoints_without_a_cache_for_the_mag(self):
        """A third timepoint in the sheet whose cache has no file for this genome
        (the genome was not QC-eligible there): skipped with a log line, no crash."""
        meta = pd.read_csv(self.metadata, sep="\t", dtype=str)
        meta = pd.concat(
            [
                meta,
                pd.DataFrame(
                    [
                        {
                            "sample_id": "m7_late",
                            "subjectID": "m7",
                            "group": "other",
                            "time": "late",
                            "replicate": "r4",
                            "bam_path": "x",
                        }
                    ]
                ),
            ]
        )
        meta.to_csv(self.metadata, sep="\t", index=False)
        done = self._run("--baseline_scope", "cohort")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("late", done.stderr)
        self.assertNotIn("m7_late", set(self._long()["sample_id"]))


class TestCheckBaselineScope(unittest.TestCase):
    """own_group needs a site that belongs to ONE group: only the within-group
    families (single_sample, lmm_across_time) have that.  The between-group
    families compare the two groups, so they are refused, not silently widened."""

    def test_comparison_scope_accepts_every_family(self):
        for family in SUMMARY_FAMILIES:
            check_baseline_scope(family, "comparison")  # no raise

    def test_own_group_accepts_within_group_families(self):
        for family in ("single_sample", "lmm_across_time"):
            check_baseline_scope(family, "own_group")  # no raise

    def test_own_group_refuses_between_group_families(self):
        # two-sample tests and the two-way (group x time) LMM compare the groups
        for family in ("two_sample_paired", "two_sample_unpaired", "lmm"):
            with self.assertRaisesRegex(ValueError, "own_group") as caught:
                check_baseline_scope(family, "own_group")
            self.assertIn(family, str(caught.exception))  # names the culprit

    def test_every_family_is_classified(self):
        # a family added to SUMMARY_FAMILIES later must be placed on one side
        self.assertTrue(set(WITHIN_GROUP_FAMILIES) <= set(SUMMARY_FAMILIES))

    def test_cmh_is_not_offered_at_all(self):
        # CMH gives one p per site, never per base, so it cannot name the allele:
        # it is refused for every scope by --summary's choices, before any check
        self.assertNotIn("cmh", SUMMARY_FAMILIES)


class TestRestrictToOwnGroup(unittest.TestCase):
    """--baseline_scope own_group: a within-group site keeps only its own group's
    samples, at both timepoints."""

    def test_keeps_own_group_at_both_timepoints(self):
        rows = [
            {"group_analyzed": "fat", "group": group, "time": time}
            for group in ("fat", "control")
            for time in ("pre", "end")
        ]
        kept = restrict_to_own_group(pd.DataFrame(rows))
        self.assertEqual(sorted(kept.group.unique()), ["fat"])  # control rows dropped
        self.assertEqual(sorted(kept.time.unique()), ["end", "pre"])  # both timepoints

    def test_between_group_row_raises(self):
        # a blank group_analyzed (a between-group site) must never get here:
        # check_baseline_scope refuses those families first
        df = pd.DataFrame({"group_analyzed": ["fat", ""], "group": ["fat", "fat"]})
        with self.assertRaisesRegex(ValueError, "between-group"):
            restrict_to_own_group(df)

    def test_numeric_group_names_compare_as_text(self):
        # group names like 40 are text throughout ("40" == "40"), never 40 vs "40"
        df = pd.DataFrame({"group_analyzed": ["40", "40"], "group": ["40", "AL"]})
        self.assertEqual(restrict_to_own_group(df).group.tolist(), ["40"])


class TestBaselineScopeWithinGroupCLI(unittest.TestCase):
    """A WITHIN-group site (single_sample, significant in fat) where only the
    CONTROL mice carry G at baseline.  Fat m1, m2: G 0/30 at pre, 24/30 and
    21/30 at end.  Control m3, m4: G 5/30 at pre (present), 0/30 at end.

    comparison (default): control's m3/m4 count too -> standing.  own_group:
    the baseline is fat's own mice -> G never seen at pre -> de_novo_candidate."""

    MAG = "MAG_A"

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        run = self.run = os.path.join(self.tmp, "run")
        prof = os.path.join(run, "profiles")

        def with_site(a, c, g, t):
            return [(p, 30, 0, 0, 0) for p in range(6) if p != 1] + [(1, a, c, g, t)]

        pre = {
            "m1": (30, 0, 0, 0),
            "m2": (30, 0, 0, 0),
            "m3": (25, 0, 5, 0),
            "m4": (25, 0, 5, 0),
        }
        end = {
            "m1": (6, 0, 24, 0),
            "m2": (9, 0, 21, 0),
            "m3": (30, 0, 0, 0),
            "m4": (30, 0, 0, 0),
        }
        for mouse in pre:
            _write_profile(prof, f"{mouse}_pre", self.MAG, with_site(*pre[mouse]))
            _write_profile(prof, f"{mouse}_end", self.MAG, with_site(*end[mouse]))
        meta = pd.DataFrame(
            [
                {
                    "sample_id": f"{m}_{t}",
                    "subjectID": m,
                    "group": g,
                    "time": t,
                    "replicate": r,
                    "bam_path": "x",
                }
                for m, g, r in (
                    ("m1", "fat", "r1"),
                    ("m2", "fat", "r2"),
                    ("m3", "control", "r3"),
                    ("m4", "control", "r4"),
                )
                for t in ("pre", "end")
            ]
        )
        self.metadata = os.path.join(self.tmp, "metadata.tsv")
        meta.to_csv(self.metadata, sep="\t", index=False)
        sd = os.path.join(run, "p_value_summary", "pre_end-fat_control")
        os.makedirs(sd)
        pd.DataFrame(
            {
                "period": ["pre_end"],
                "mag_id": [self.MAG],
                "contig": ["c1"],
                "position": [1],
                "gene_id": ["g1"],
                "test_type": ["single_sample_tTest"],
                "group_analyzed": ["fat"],
                "min_p_value": [0.001],
                "source_file": [f"{self.MAG}_single_sample_fat.tsv.gz"],
                "q_value": [0.01],
            }
        ).to_csv(
            os.path.join(sd, "p_value_summary_single_sample_pre_end.tsv"),
            sep="\t",
            index=False,
        )
        td = os.path.join(
            run, "significance_tests", "single_sample_pre_end-fat_control"
        )
        os.makedirs(td)
        # A and G are the seesaw pair (A falls as G rises): one test seen twice,
        # p-values equal up to float32 noise.  Both are candidates; the rising
        # table, not the p-values, decides that G is the one that rose.
        pd.DataFrame(
            {
                "contig": ["c1"],
                "gene_id": ["g1"],
                "position": [1],
                "A_frequency_p_value_tTest_fat": [0.001000001],
                "C_frequency_p_value_tTest_fat": [1.0],
                "G_frequency_p_value_tTest_fat": [0.001],
                "T_frequency_p_value_tTest_fat": [1.0],
            }
        ).to_csv(
            os.path.join(td, f"{self.MAG}_single_sample_fat.tsv.gz"),
            sep="\t",
            index=False,
            compression="gzip",
        )
        self.fasta = os.path.join(self.tmp, "ref.fa")
        open(self.fasta, "w").write(">c1\nAAAAAA\n")
        open(self.fasta + ".fai", "w").write("c1\t6\t4\t6\t7\n")
        self.mapping = os.path.join(self.tmp, "map.tsv")
        pd.DataFrame({"mag_id": [self.MAG], "contig_id": ["c1"]}).to_csv(
            self.mapping, sep="\t", index=False
        )
        self.out = os.path.join(self.tmp, "out")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self, *extra):
        done = subprocess.run(
            [
                "alleleflux-baseline-presence",
                "--run_dir",
                self.run,
                "--comparison",
                "pre_end-fat_control",
                "--summary",
                "single_sample",
                "--test_type",
                "single_sample_tTest",
                "--cache_dir",
                os.path.join(self.run, "allele_freq_cache"),
                "--metadata",
                self.metadata,
                "--output_dir",
                self.out,
                "--cpus",
                "2",
                *extra,
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        stem = os.path.join(
            self.out, "pre_end-fat_control_single_sample_tTest_baseline_presence"
        )
        summary = pd.read_csv(stem + "_summary.tsv", sep="\t")
        self.assertEqual(sorted(summary.allele), ["A", "G"])  # both seesaw partners
        return (
            pd.read_csv(stem + ".tsv.gz", sep="\t"),
            summary[summary.allele == "G"].iloc[0],  # the riser's row
        )

    def test_own_group_scope_uses_only_the_sites_group(self):
        long, s = self._run("--baseline_scope", "own_group")
        self.assertEqual(
            sorted(long.subjectID.unique()), ["m1", "m2"]
        )  # control mice gone, both timepoints
        self.assertEqual(len(long), 8)  # 2 alleles x 2 fat mice x 2 timepoints
        self.assertEqual(s.allele, "G")
        self.assertEqual(s.origin_any_mouse, "de_novo_candidate")
        self.assertEqual(
            (int(s.n_pre_samples_allele_present), int(s.n_pre_samples_covered)), (0, 2)
        )
        self.assertEqual(
            (int(s.total_reads_pre), int(s.allele_reads_pre)), (60, 0)
        )  # fat mice only
        self.assertEqual(int(s.n_mice_de_novo_candidate), 2)

    def test_default_is_comparison_both_groups(self):
        long, s = self._run()
        self.assertEqual(len(long), 16)  # 2 alleles x all four mice x 2 timepoints
        self.assertEqual(
            s.origin_any_mouse, "standing_variation"
        )  # m3, m4 carry G at pre
        self.assertEqual(
            (int(s.n_pre_samples_allele_present), int(s.n_pre_samples_covered)), (2, 4)
        )
        self.assertEqual((int(s.total_reads_pre), int(s.allele_reads_pre)), (120, 10))

    def test_unknown_scope_is_rejected(self):
        done = subprocess.run(
            ["alleleflux-baseline-presence", "--baseline_scope", "everyone", "--help"],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(done.returncode, 0)

    def _rising(self):
        return pd.read_csv(
            os.path.join(
                self.out,
                "pre_end-fat_control_single_sample_tTest_baseline_presence_rising.tsv",
            ),
            sep="\t",
        )

    def test_rising_table_under_each_scope(self):
        # comparison: control is in scope, so its change is reported alongside
        self._run()
        rise = self._rising()
        self.assertEqual(len(rise), 1)  # the site's own group only
        r = rise.iloc[0]
        self.assertEqual((r.allele, r.rose_in), ("G", "fat"))
        self.assertAlmostEqual(r.frequency_change, (0.8 + 0.7) / 2)  # r1, r2
        self.assertEqual((int(r.n_replicates), int(r.n_mice)), (2, 2))
        self.assertEqual(r.origin_any_mouse, "standing_variation")  # control had G
        # own_group: fat alone; baseline = fat's mice -> de novo
        shutil.rmtree(self.out)
        self._run("--baseline_scope", "own_group")
        r = self._rising().iloc[0]
        self.assertEqual(r.origin_any_mouse, "de_novo_candidate")

    def test_rising_table_has_no_other_group_columns(self):
        # dropped: under own_group the other group's mice are gone before counting,
        # so the columns were always blank there, and nothing downstream reads them
        self._run()
        cols = set(self._rising().columns)
        self.assertNotIn("other_group", cols)
        self.assertNotIn("other_group_frequency_change", cols)


if __name__ == "__main__":
    unittest.main()
