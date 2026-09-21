"""Tests for rolling per-mouse strain-background calls up to MAGs.

Every MAG-level row carries two blocks: counts over MICE, and counts over
REPLICATES (a replicate "changed" if ANY of its mice changed).  With no
replicate column in the metadata replicate == mouse and the blocks agree.
"""
import unittest

import pandas as pd

from alleleflux.scripts.analysis.ani.replacement_classification import METRICS, classify_mags


def _call_row(mag, mouse, group, transition, replacement, dominant, replicate=None):
    """One turnover-table row (Task 6 schema).  ``None`` for a flag means undetermined (pd.NA)."""
    return {
        "MAG_ID": mag, "subjectID": mouse, "group": group,
        "replicate": replicate or mouse, "transition": transition,
        "strain_replacement": pd.NA if replacement is None else replacement,
        "dominant_strain_change": pd.NA if dominant is None else dominant,
    }


def _frame(rows):
    df = pd.DataFrame(rows)
    for col in METRICS:
        df[col] = df[col].astype("boolean")
    return df


def _classify(df, **kwargs):
    """``classify_mags`` with the voter floor pinned to 1.

    The fixtures below have 2-4 mice per key, far under the production default
    of 8, so with the default every verdict would be blank.  These classes test
    the VOTING ARITHMETIC (majority / any / all / none); the floor itself has
    its own class, ``TestVoterFloorAndStatus``.
    """
    return classify_mags(df, **{"min_voters": 1, **kwargs})


class TestMouseBlock(unittest.TestCase):
    ROWS = [
        # MAG_A, group 40, 3 mice: replacement T/T/F -> majority; dominant F/F/F -> none.
        _call_row("MAG_A", "m1", "40", "5mo_22mo", True, False),
        _call_row("MAG_A", "m2", "40", "5mo_22mo", True, False),
        _call_row("MAG_A", "m3", "40", "5mo_22mo", False, False),
        # MAG_B, group 40, 2 called mice (T/F) + 1 undetermined -> denominator 2, no majority.
        _call_row("MAG_B", "m4", "40", "5mo_22mo", True, True),
        _call_row("MAG_B", "m5", "40", "5mo_22mo", False, False),
        _call_row("MAG_B", "m6", "40", "5mo_22mo", None, None),
        # MAG_B, group 20, everything undetermined -> no verdict of any kind.
        _call_row("MAG_B", "m7", "20", "5mo_22mo", None, None),
        # MAG_C, group 40, 2 mice both changed -> all_changed.
        _call_row("MAG_C", "m8", "40", "5mo_22mo", True, True),
        _call_row("MAG_C", "m9", "40", "5mo_22mo", True, False),
    ]

    def _mags(self, metric):
        both = _classify(_frame(self.ROWS))
        return both[both.metric == metric].set_index(["MAG_ID", "group"])

    def test_majority_uses_only_called_mice(self):
        got = self._mags("strain_replacement")
        a = got.loc[("MAG_A", "40")]
        self.assertEqual((a.n_mice_with_call, a.n_mice_undetermined, a.n_mice_changed), (3, 0, 2))
        self.assertTrue(a.majority_mice_changed and a.any_mouse_changed)
        self.assertFalse(a.all_mice_changed or a.no_mouse_changed)
        b = got.loc[("MAG_B", "40")]
        self.assertEqual((b.n_mice_with_call, b.n_mice_undetermined, b.n_mice_changed), (2, 1, 1))
        self.assertFalse(b.majority_mice_changed)      # 1 of 2 is NOT a majority
        self.assertTrue(b.any_mouse_changed)

    def test_all_changed(self):
        got = self._mags("strain_replacement")
        self.assertTrue(got.loc[("MAG_C", "40")].all_mice_changed)
        self.assertFalse(got.loc[("MAG_A", "40")].all_mice_changed)   # 2 of 3
        # dominant: MAG_C is T/F -> not all
        self.assertFalse(self._mags("dominant_strain_change").loc[("MAG_C", "40")].all_mice_changed)

    def test_no_change_and_all_changed_need_at_least_one_called_mouse(self):
        got = self._mags("dominant_strain_change")
        self.assertTrue(got.loc[("MAG_A", "40")].no_mouse_changed)
        b20 = got.loc[("MAG_B", "20")]
        self.assertEqual(b20.n_mice_with_call, 0)
        for col in ("majority_mice_changed", "any_mouse_changed", "all_mice_changed", "no_mouse_changed"):
            # zero evidence is not a verdict: the cell is BLANK, never False.
            # False would read as "checked, and it did not change".
            self.assertTrue(pd.isna(b20[col]), col)
        self.assertEqual(b20.strain_status, "no_voters")

    def test_metric_column_names_the_evidence(self):
        got = _classify(_frame(self.ROWS))
        self.assertEqual(got.metric.value_counts().to_dict(),
                         {"strain_replacement": 4, "dominant_strain_change": 4})
        self.assertEqual(list(got.columns[:4]), ["MAG_ID", "group", "transition", "metric"])

    def test_groups_are_never_pooled(self):
        got = self._mags("strain_replacement").reset_index()
        self.assertEqual(sorted(zip(got.MAG_ID, got.group)),
                         [("MAG_A", "40"), ("MAG_B", "20"), ("MAG_B", "40"), ("MAG_C", "40")])



class TestReplicateBlock(unittest.TestCase):
    ROWS = [
        # cage c1: m1 changed, m2 not  -> the cage counts as changed (ANY mouse).
        _call_row("MAG_A", "m1", "40", "5mo_22mo", True, False, replicate="c1"),
        _call_row("MAG_A", "m2", "40", "5mo_22mo", False, False, replicate="c1"),
        # cage c2: m3 alone, not changed.
        _call_row("MAG_A", "m3", "40", "5mo_22mo", False, False, replicate="c2"),
        # cage c3: only an undetermined mouse -> the cage has no call.
        _call_row("MAG_A", "m4", "40", "5mo_22mo", None, None, replicate="c3"),
    ]

    def test_replicate_counts_sit_beside_mouse_counts(self):
        both = _classify(_frame(self.ROWS))
        row = both[both.metric == "strain_replacement"].iloc[0]
        # mice: 3 called (m1,m2,m3), 1 changed -> 1 of 3 is no majority
        self.assertEqual((row.n_mice_with_call, row.n_mice_changed), (3, 1))
        self.assertFalse(row.majority_mice_changed)
        # replicates: c1 and c2 called, c3 not; c1 changed -> 1 of 2, still no majority
        self.assertEqual((row.n_replicates_with_call, row.n_replicates_changed), (2, 1))
        self.assertFalse(row.majority_replicates_changed)
        self.assertFalse(row.all_replicates_changed)

    def test_replicate_majority_can_differ_from_mouse_majority(self):
        # c1 = {m1 changed, m2 not, m3 not}; c2 = {m4 changed}.  Mice: 2 of 4 (no
        # majority).  Cages: both changed -> majority AND all.
        rows = [
            _call_row("MAG_A", "m1", "40", "5mo_22mo", True, False, replicate="c1"),
            _call_row("MAG_A", "m2", "40", "5mo_22mo", False, False, replicate="c1"),
            _call_row("MAG_A", "m3", "40", "5mo_22mo", False, False, replicate="c1"),
            _call_row("MAG_A", "m4", "40", "5mo_22mo", True, False, replicate="c2"),
        ]
        both = _classify(_frame(rows))
        row = both[both.metric == "strain_replacement"].iloc[0]
        self.assertFalse(row.majority_mice_changed)
        self.assertTrue(row.majority_replicates_changed)
        self.assertTrue(row.all_replicates_changed)

    def test_replicate_equal_to_mouse_makes_the_blocks_agree(self):
        # DRiDO: no replicate column -> replicate == subjectID.
        rows = [_call_row("MAG_A", m, "40", "5mo_22mo", flag, False)
                for m, flag in (("m1", True), ("m2", False), ("m3", None))]
        both = _classify(_frame(rows))
        row = both[both.metric == "strain_replacement"].iloc[0]
        self.assertEqual(row.n_replicates_with_call, row.n_mice_with_call)
        self.assertEqual(row.n_replicates_changed, row.n_mice_changed)
        self.assertEqual(row.majority_replicates_changed, row.majority_mice_changed)


def _key_rows(mag, n_changed, n_same, n_undetermined=0, group="40", transition="5mo_22mo"):
    """Mice for ONE (MAG, group, transition) key: so many changed, same, undetermined.

    Only ``strain_replacement`` varies; ``dominant_strain_change`` is False for
    every called mouse, so the two metrics give different statuses on purpose.
    """
    flags = [True] * n_changed + [False] * n_same + [None] * n_undetermined
    return [_call_row(mag, f"{mag}_m{i}", group, transition, flag, None if flag is None else False)
            for i, flag in enumerate(flags)]


class TestVoterFloorAndStatus(unittest.TestCase):
    """The voter floor, the blank verdicts under it, and the ``strain_status`` word.

    Five keys, one per situation (voters = called mice)::

        MAG_R  8 voters, 5 changed             -> replaced
        MAG_N  9 voters, 4 changed             -> not_replaced
        MAG_T  8 voters, 4 changed (a tie)     -> not_replaced by default
        MAG_F  5 voters, 3 changed             -> too_few_voters (a majority of 5, but under the floor)
        MAG_Z  0 voters, 2 undetermined        -> no_voters
    """

    ROWS = (_key_rows("MAG_R", 5, 3) + _key_rows("MAG_N", 4, 5) + _key_rows("MAG_T", 4, 4)
            + _key_rows("MAG_F", 3, 2, n_undetermined=6) + _key_rows("MAG_Z", 0, 0, n_undetermined=2))
    VERDICTS = ("majority_mice_changed", "any_mouse_changed", "all_mice_changed", "no_mouse_changed")

    def _mags(self, metric="strain_replacement", **kwargs):
        both = classify_mags(_frame(self.ROWS), **kwargs)
        return both[both.metric == metric].set_index("MAG_ID")

    def test_default_floor_is_eight_voters(self):
        from alleleflux.scripts.analysis.ani.replacement_classification import DEFAULT_MIN_VOTERS
        self.assertEqual(DEFAULT_MIN_VOTERS, 8)
        got = self._mags()                                   # no min_voters passed
        self.assertEqual(got.strain_status.to_dict(), {
            "MAG_R": "replaced", "MAG_N": "not_replaced", "MAG_T": "not_replaced",
            "MAG_F": "too_few_voters", "MAG_Z": "no_voters"})

    def test_below_the_floor_every_verdict_is_blank_but_counts_stay(self):
        f = self._mags().loc["MAG_F"]
        # the facts are still reported ...
        self.assertEqual((f.n_mice_with_call, f.n_mice_undetermined, f.n_mice_changed), (5, 6, 3))
        # ... but 3 of 5 is NOT allowed to read as "majority changed": blank, not True
        for col in self.VERDICTS:
            self.assertTrue(pd.isna(f[col]), col)

    def test_at_or_above_the_floor_verdicts_are_real_booleans(self):
        got = self._mags()
        self.assertTrue(got.loc["MAG_R"].majority_mice_changed)      # 5 of 8
        self.assertFalse(got.loc["MAG_N"].majority_mice_changed)     # 4 of 9
        for mag in ("MAG_R", "MAG_N", "MAG_T"):
            for col in self.VERDICTS:
                self.assertFalse(pd.isna(got.loc[mag][col]), (mag, col))

    def test_status_is_per_metric(self):
        # dominant_strain_change is False for every called mouse -> nobody replaced
        got = self._mags("dominant_strain_change")
        self.assertEqual(got.loc["MAG_R"].strain_status, "not_replaced")
        self.assertEqual(got.loc["MAG_F"].strain_status, "too_few_voters")

    def test_lowering_the_floor_lets_the_five_voter_key_speak(self):
        got = self._mags(min_voters=5)
        self.assertEqual(got.loc["MAG_F"].strain_status, "replaced")  # 3 of 5
        self.assertTrue(got.loc["MAG_F"].majority_mice_changed)
        self.assertEqual(got.loc["MAG_Z"].strain_status, "no_voters")  # zero voters never pass any floor

    def test_tie_rule_toggle(self):
        self.assertEqual(self._mags().loc["MAG_T"].strain_status, "not_replaced")            # default
        self.assertEqual(self._mags(tie="not_replaced").loc["MAG_T"].strain_status, "not_replaced")
        self.assertEqual(self._mags(tie="replaced").loc["MAG_T"].strain_status, "replaced")
        self.assertEqual(self._mags(tie="unresolved").loc["MAG_T"].strain_status, "tie")
        for tie in ("not_replaced", "replaced", "unresolved"):
            got = self._mags(tie=tie)
            # the majority column states a FACT (4 is not more than 8/2) and ignores the tie rule
            self.assertFalse(got.loc["MAG_T"].majority_mice_changed, tie)
            # ... and the tie rule touches nothing but exact ties
            self.assertEqual(got.loc["MAG_R"].strain_status, "replaced", tie)
            self.assertEqual(got.loc["MAG_N"].strain_status, "not_replaced", tie)

    def test_vote_rule_toggle(self):
        # Same five keys, three ways of turning votes into "replaced".
        #   MAG_R 5 of 8 changed, MAG_N 4 of 9, MAG_T 4 of 8 -- all pass the floor.
        majority = self._mags(vote_rule="majority").strain_status.to_dict()
        self.assertEqual((majority["MAG_R"], majority["MAG_N"], majority["MAG_T"]),
                         ("replaced", "not_replaced", "not_replaced"))
        any_rule = self._mags(vote_rule="any").strain_status.to_dict()       # one changed mouse is enough
        self.assertEqual((any_rule["MAG_R"], any_rule["MAG_N"], any_rule["MAG_T"]),
                         ("replaced", "replaced", "replaced"))
        all_rule = self._mags(vote_rule="all").strain_status.to_dict()       # every voter must have changed
        self.assertEqual((all_rule["MAG_R"], all_rule["MAG_N"], all_rule["MAG_T"]),
                         ("not_replaced", "not_replaced", "not_replaced"))
        # the floor holds under every rule: 3 changed of 5 voters is still no verdict
        for rule in ("majority", "any", "all"):
            got = self._mags(vote_rule=rule)
            self.assertEqual(got.loc["MAG_F"].strain_status, "too_few_voters", rule)
            self.assertEqual(got.loc["MAG_Z"].strain_status, "no_voters", rule)
            self.assertEqual(set(got.vote_rule), {rule})                     # stamped on every row

    def test_all_rule_flags_a_unanimous_key(self):
        rows = _key_rows("MAG_U", 8, 0) + _key_rows("MAG_V", 7, 1)
        got = classify_mags(_frame(rows), vote_rule="all")
        got = got[got.metric == "strain_replacement"].set_index("MAG_ID")
        self.assertEqual(got.strain_status.to_dict(), {"MAG_U": "replaced", "MAG_V": "not_replaced"})

    def test_status_never_disagrees_with_the_matching_yes_no_column(self):
        # The status word and the yes/no columns answer from the SAME raw
        # answers.  For every key above the floor (exact ties aside, which the
        # tie rule owns), "replaced" must mean exactly "the chosen rule's column
        # is True".  If the two were ever computed separately and one definition
        # drifted, this is the test that notices.
        column_of = {"majority": "majority_mice_changed", "any": "any_mouse_changed",
                     "all": "all_mice_changed"}
        for rule, column in column_of.items():
            got = self._mags(vote_rule=rule)
            judged = got[got.n_mice_with_call >= 8]
            if rule == "majority":
                judged = judged[judged.n_mice_changed * 2 != judged.n_mice_with_call]
            self.assertGreater(len(judged), 0, rule)
            for mag, row in judged.iterrows():
                self.assertEqual(row.strain_status == "replaced", bool(row[column]), (rule, mag))

    def test_tie_rule_only_matters_under_majority(self):
        # 4 of 8 is not a tie for "any" (4 > 0 -> replaced) nor for "all" (4 != 8 -> not replaced)
        self.assertEqual(self._mags(vote_rule="any", tie="unresolved").loc["MAG_T"].strain_status, "replaced")
        self.assertEqual(self._mags(vote_rule="all", tie="replaced").loc["MAG_T"].strain_status, "not_replaced")

    def test_bad_settings_raise(self):
        with self.assertRaises(ValueError):
            classify_mags(_frame(self.ROWS), vote_rule="plurality")
        with self.assertRaises(ValueError):
            classify_mags(_frame(self.ROWS), tie="coin_flip")
        with self.assertRaises(ValueError):
            classify_mags(_frame(self.ROWS), min_voters=0)

    def test_every_row_says_which_rule_made_it(self):
        got = classify_mags(_frame(self.ROWS), min_voters=5, tie="replaced")
        self.assertEqual(set(got.min_voters), {5})
        self.assertEqual(set(got.tie_rule), {"replaced"})

    def test_replicate_status_matches_mouse_status_when_replicate_is_the_mouse(self):
        got = self._mags()
        self.assertEqual(got.strain_status_replicates.to_dict(), got.strain_status.to_dict())
        self.assertTrue(pd.isna(got.loc["MAG_F"].majority_replicates_changed))


if __name__ == "__main__":
    unittest.main()


# ---------------------------------------------------------------------------
# End-to-end: the ``alleleflux-replacement-classification`` command
# ---------------------------------------------------------------------------
import os
import shutil
import subprocess
import tempfile

from alleleflux.scripts.analysis.ani.strain_turnover import TURNOVER_COLUMNS


def _turnover_row(mag, mouse, group, transition, background, replacement, dominant):
    """One row in the exact {mag}_strain_turnover.tsv schema (Task 6 output)."""
    row = {c: "" for c in TURNOVER_COLUMNS}
    row.update({
        "MAG_ID": mag, "subjectID": mouse, "replicate": mouse, "group": group,
        "transition": transition, "sample_t1": f"{mouse}_a", "sample_t2": f"{mouse}_b",
        "compared_bases_count": 1000, "percent_genome_compared": 0.5,
        "conANI": 0.9999, "popANI": 0.99999, "frequency_shift": 0.0,
        "strain_replacement": replacement, "dominant_strain_change": dominant,
        "background": background, "min_compared": 0.1, "pop_threshold": 0.99999,
        "con_threshold": 0.999, "min_cov": 5,
    })
    return row


class TestReplacementClassificationCLI(unittest.TestCase):
    """Two MAGs' turnover files plus one header-only file in a directory.

    MAG_A, group 40, 5mo_22mo: 3 mice, replacement T/T/F -> majority; dominant F/F/F.
    MAG_B, group 40, 5mo_22mo: 1 mouse undetermined -> no verdict of any kind.
    MAG_C: header only (a MAG with no matching transitions) -> contributes nothing.
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.turnover = os.path.join(self.tmp, "strain_turnover")
        os.makedirs(self.turnover)
        rows_a = [
            _turnover_row("MAG_A", "m1", "40", "5mo_22mo", "strain_replacement", True, False),
            _turnover_row("MAG_A", "m2", "40", "5mo_22mo", "strain_replacement", True, False),
            _turnover_row("MAG_A", "m3", "40", "5mo_22mo", "stable", False, False),
        ]
        rows_b = [_turnover_row("MAG_B", "m4", "40", "5mo_22mo", "undetermined", "", "")]
        for mag, rows in (("MAG_A", rows_a), ("MAG_B", rows_b), ("MAG_C", [])):
            pd.DataFrame(rows, columns=TURNOVER_COLUMNS).to_csv(
                os.path.join(self.turnover, f"{mag}_strain_turnover.tsv"), sep="\t", index=False)
        self.out = os.path.join(self.tmp, "replacement_classification.tsv")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self, turnover_dir=None, extra=()):
        return subprocess.run([
            "alleleflux-replacement-classification",
            "--turnover_dir", turnover_dir or self.turnover, "--output_path", self.out,
            *extra,
        ], capture_output=True, text=True)

    def test_default_floor_blanks_the_three_mouse_mag_on_disk(self):
        # No flags: the production default of 8 voters.  MAG_A has 3, so it gets
        # no verdict -- and the cells must be EMPTY in the file, not "False".
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        raw = pd.read_csv(self.out, sep="\t", dtype=str, keep_default_na=False)
        a = raw[(raw.MAG_ID == "MAG_A") & (raw.metric == "strain_replacement")].iloc[0]
        self.assertEqual(a.strain_status, "too_few_voters")
        self.assertEqual((a.n_mice_with_call, a.n_mice_changed), ("3", "2"))     # facts survive
        for col in ("majority_mice_changed", "any_mouse_changed", "all_mice_changed", "no_mouse_changed"):
            self.assertEqual(a[col], "", col)                                    # blank, never "False"
        b = raw[(raw.MAG_ID == "MAG_B") & (raw.metric == "strain_replacement")].iloc[0]
        self.assertEqual(b.strain_status, "no_voters")
        self.assertEqual((a.min_voters, a.tie_rule), ("8", "not_replaced"))

    def test_tie_flag_reaches_the_table(self):
        done = self._run(extra=("--min_voters", "1", "--tie", "unresolved"))
        self.assertEqual(done.returncode, 0, done.stderr)
        got = pd.read_csv(self.out, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(set(got.tie_rule), {"unresolved"})
        self.assertEqual(set(got.min_voters), {"1"})

    def test_classifies_every_mag_for_both_metrics(self):
        done = self._run(extra=("--min_voters", "1"))
        self.assertEqual(done.returncode, 0, done.stderr)
        got = pd.read_csv(self.out, sep="\t", dtype={"group": str})
        # 2 MAGs with rows x 2 metrics; the header-only MAG_C contributes nothing
        self.assertEqual(sorted(got.MAG_ID.unique()), ["MAG_A", "MAG_B"])
        self.assertEqual(got.metric.value_counts().to_dict(), {"strain_replacement": 2, "dominant_strain_change": 2})
        a = got[(got.MAG_ID == "MAG_A") & (got.metric == "strain_replacement")].iloc[0]
        self.assertEqual((int(a.n_mice_with_call), int(a.n_mice_changed)), (3, 2))
        self.assertTrue(bool(a.majority_mice_changed))
        b = got[(got.MAG_ID == "MAG_B") & (got.metric == "strain_replacement")].iloc[0]
        self.assertEqual((int(b.n_mice_with_call), int(b.n_mice_undetermined)), (0, 1))
        # zero voters: blank on disk -> NaN on read, never False
        self.assertTrue(pd.isna(b.majority_mice_changed) and pd.isna(b.no_mouse_changed))
        self.assertEqual((a.strain_status, b.strain_status), ("replaced", "no_voters"))
        self.assertEqual(list(got.columns[:4]), ["MAG_ID", "group", "transition", "metric"])

    def test_directory_without_turnover_files_fails_loud(self):
        empty = os.path.join(self.tmp, "nothing"); os.makedirs(empty)
        done = self._run(empty)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("strain_turnover", done.stderr)
