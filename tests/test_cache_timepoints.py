"""Unit tests for the parse-time "which timepoints get QC and caches" helpers in ``common.smk``.

Under a dated ``cohort`` baseline scope the de novo analyses judge "absent at
baseline" against every sample of the metadata sheet dated before the first
carrier, whatever its timepoint.  A timepoint can be in the sheet and in no
``timepoints_combinations`` entry (for example animals first sampled at an age
nobody compares); its samples are witnesses all the same, so the workflow must
build QC and allele-frequency caches there too.  The same holds for a group
that is in the sheet and in no ``groups_combinations`` entry.  Five helpers
decide this while the Snakefile is being read:

* ``sheet_fully_dated``        -- does every sample carry a collection date?
* ``cohort_witnesses_active``  -- do the de novo analyses read witnesses
  outside the comparison at all?  Only for longitudinal data,
  ``baseline_scope: cohort``, a fully dated sheet AND at least one de novo
  analysis switched on (user rule: no extra work otherwise).
* ``cohort_reads_whole_sheet`` -- are the witnesses the whole sheet
  (``cohort_witnesses: sheet``, the default) or only the samples at a compared
  timepoint in a compared group (``comparisons``)?
* ``extra_cache_timepoints``   -- the sheet timepoints outside every
  comparison, ordered by their first collection date.
* ``extra_witness_groups``     -- the sheet groups outside every comparison.

Sliced out of the ``.smk`` source and exec'd alone, as in
``test_timepoint_order.py``; pandas is handed in because the helpers use it.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

COMMON_SMK = (
    Path(__file__).resolve().parents[1]
    / "alleleflux"
    / "smk_workflow"
    / "alleleflux_pipeline"
    / "shared"
    / "common.smk"
)
FUNCS = (
    "sheet_fully_dated",
    "extra_cache_timepoints",
    "extra_witness_groups",
    "cohort_witnesses_active",
    "cohort_reads_whole_sheet",
)


def _slice_function(lines, name):
    """The source of one top-level function of the .smk file (def line to next top-level line)."""
    start = next(i for i, ln in enumerate(lines) if ln.startswith(f"def {name}("))
    end = start + 1
    while end < len(lines):
        ln = lines[end]
        if ln.strip() and not ln[0].isspace():
            break
        end += 1
    return "\n".join(lines[start:end])


def _load():
    """Exec the helpers alone and hand them back by name."""
    lines = COMMON_SMK.read_text().splitlines()
    src = "\n\n".join(_slice_function(lines, f) for f in FUNCS)
    namespace = {"pd": pd}
    exec(compile(src, str(COMMON_SMK), "exec"), namespace)
    return {f: namespace[f] for f in FUNCS}


def _sheet(times, dates=None, groups=None):
    """A metadata sheet as the workflow reads it (every cell text): one sample per row.

    ``dates`` None -> no ``collection_date`` column at all; a None entry -> a
    blank cell (pandas reads a blank TSV cell as NaN under ``dtype=str``).
    ``groups`` None -> no ``group`` column; a None entry -> a blank cell.
    """
    data = {
        "sample_id": [f"s{i}" for i in range(len(times))],
        "subjectID": [f"m{i}" for i in range(len(times))],
        "time": [str(t) for t in times],
    }
    if groups is not None:
        data["group"] = groups
    if dates is not None:
        data["collection_date"] = dates
    return pd.DataFrame(data, dtype=object)


class TestExtraCacheTimepoints:
    def test_extra_timepoints_are_the_sheet_ages_outside_every_comparison_in_date_order(
        self,
    ):
        # x is collected first of all, y last; a and b are compared, so only x
        # and y are extra, and x (earliest first date) comes before y.
        sheet = _sheet(
            times=["b", "a", "y", "x"],
            dates=["2020-02-01", "2020-01-01", "2021-01-01", "2019-06-01"],
        )
        assert _load()["extra_cache_timepoints"](sheet, ["a", "b"]) == ["x", "y"]

    def test_order_uses_the_first_collection_date_of_each_timepoint(self):
        # y has one late sample and one very early one: its FIRST date (2018)
        # precedes x's (2019), so y leads although its other sample is the latest.
        sheet = _sheet(
            times=["a", "x", "y", "y"],
            dates=["2020-01-01", "2019-06-01", "2022-01-01", "2018-03-01"],
        )
        assert _load()["extra_cache_timepoints"](sheet, ["a"]) == ["y", "x"]

    def test_timepoints_first_collected_on_the_same_day_are_ordered_by_label(self):
        # Same first date -> the label breaks the tie, so the list (and with it
        # the wildcard constraint and the DAG) never depends on row order.
        sheet = _sheet(
            times=["a", "z", "k"], dates=["2020-01-01", "2019-06-01", "2019-06-01"]
        )
        assert _load()["extra_cache_timepoints"](sheet, ["a"]) == ["k", "z"]

    def test_extra_timepoints_empty_when_every_sheet_age_is_in_a_comparison(self):
        sheet = _sheet(times=["a", "b", "a"], dates=["2020-01-01"] * 3)
        assert _load()["extra_cache_timepoints"](sheet, ["a", "b"]) == []

    def test_label_with_underscore_is_refused(self):
        # "x_1" could never be told apart from the comparison label "x" + "1"
        # (labels are joined with "_"), so parse time must stop, naming it.
        sheet = _sheet(times=["a", "x_1"], dates=["2020-01-01", "2019-01-01"])
        with pytest.raises(ValueError, match="x_1"):
            _load()["extra_cache_timepoints"](sheet, ["a"])

    @pytest.mark.parametrize("label", ["x 1", "x/1", "x|1", "é"])
    def test_label_that_is_not_a_plain_wildcard_token_is_refused(self, label):
        sheet = _sheet(times=["a", label], dates=["2020-01-01", "2019-01-01"])
        with pytest.raises(ValueError, match="timepoint label"):
            _load()["extra_cache_timepoints"](sheet, ["a"])

    def test_numeric_labels_compare_as_strings(self):
        # The sheet says 5 and 8 (text); the config may hand the comparison
        # timepoint over as the number 5.  "5" must still count as compared.
        sheet = _sheet(times=["5", "8"], dates=["2020-01-01", "2019-01-01"])
        fn = _load()["extra_cache_timepoints"]
        assert fn(sheet, [5]) == ["8"]
        assert fn(sheet, ["5"]) == ["8"]

    def test_blank_timepoint_cell_is_refused_naming_the_sample(self):
        # A blank ``time`` cell would otherwise turn into the label "nan" and
        # ask for QC and caches at a timepoint that does not exist.
        sheet = _sheet(times=["a", "x"], dates=["2020-01-01", "2019-01-01"])
        sheet.loc[1, "time"] = None
        with pytest.raises(ValueError, match="s1"):
            _load()["extra_cache_timepoints"](sheet, ["a"])

    def test_unreadable_date_at_an_extra_timepoint_is_refused_naming_it(self):
        # The order check only reads the comparison timepoints' dates; a bad
        # date at an extra timepoint must not slip through as a guessed order.
        sheet = _sheet(times=["a", "x"], dates=["2020-01-01", "01/06/2019"])
        with pytest.raises(ValueError, match="01/06/2019"):
            _load()["extra_cache_timepoints"](sheet, ["a"])


class TestCohortWitnessesActive:
    def test_flag_true_only_for_dated_longitudinal_cohort_with_a_de_novo_consumer(
        self,
    ):
        fn = _load()["cohort_witnesses_active"]
        # (data_type, baseline_scope, use_baseline_presence, use_gene_de_novo, sheet_dated)
        assert fn("longitudinal", "cohort", True, False, True)
        assert fn("longitudinal", "cohort", False, True, True)
        assert fn("longitudinal", "cohort", True, True, True)
        # undated sheet: the consumers narrow the scope themselves, nothing extra to build
        assert not fn("longitudinal", "cohort", True, True, False)
        # another scope never reads samples outside the comparison
        assert not fn("longitudinal", "comparison", True, True, True)
        # single-timepoint data has no de novo analysis
        assert not fn("single", "cohort", True, True, True)
        # user rule: no de novo analysis on -> no extra work
        assert not fn("longitudinal", "cohort", False, False, True)

    def test_flag_is_a_plain_bool(self):
        fn = _load()["cohort_witnesses_active"]
        assert fn("longitudinal", "cohort", True, False, True) is True
        assert fn("single", "cohort", True, False, True) is False


class TestCohortReadsWholeSheet:
    def test_whole_sheet_needs_the_setting_to_say_sheet(self):
        fn = _load()["cohort_reads_whole_sheet"]
        # (witnesses_active, cohort_witnesses)
        assert fn(True, "sheet") is True
        assert fn(True, "comparisons") is False

    def test_the_setting_is_ignored_when_no_witness_is_read(self):
        # Another scope, an undated sheet or no de novo analysis: the config
        # template always carries the line, so "sheet" there must change nothing.
        fn = _load()["cohort_reads_whole_sheet"]
        assert fn(False, "sheet") is False
        assert fn(False, "comparisons") is False

    def test_comparisons_setting_adds_no_timepoint_and_no_group(self):
        # What the module block does with the two answers: the extra lists are
        # asked for only when the whole sheet is read.
        fns = _load()
        sheet = _sheet(
            times=["a", "b", "x"],
            dates=["2020-01-01", "2020-02-01", "2019-06-01"],
            groups=["g1", "g2", "g3"],
        )
        for setting, timepoints, groups in (
            ("sheet", ["x"], ["g3"]),
            ("comparisons", [], []),
        ):
            whole = fns["cohort_reads_whole_sheet"](True, setting)
            assert (
                fns["extra_cache_timepoints"](sheet, ["a", "b"]) if whole else []
            ) == timepoints
            assert (
                fns["extra_witness_groups"](sheet, ["g1", "g2"]) if whole else []
            ) == groups


class TestExtraWitnessGroups:
    def test_extra_witness_groups_are_the_sheet_groups_in_no_comparison(self):
        sheet = _sheet(times=["a"] * 5, groups=["g1", "zeta", "g2", "alpha", "zeta"])
        assert _load()["extra_witness_groups"](sheet, ["g1", "g2"]) == [
            "alpha",
            "zeta",
        ]

    def test_empty_when_every_sheet_group_is_compared(self):
        sheet = _sheet(times=["a"] * 3, groups=["g1", "g2", "g1"])
        assert _load()["extra_witness_groups"](sheet, ["g1", "g2"]) == []

    def test_numeric_group_labels_are_compared_as_text(self):
        # The sheet is read as text ("20"); a config may hold the number 20.
        sheet = _sheet(times=["a"] * 3, groups=["20", "40", "60"])
        assert _load()["extra_witness_groups"](sheet, [20, "40"]) == ["60"]

    def test_surrounding_spaces_do_not_make_a_new_group(self):
        sheet = _sheet(times=["a"] * 2, groups=["g1 ", " g2"])
        assert _load()["extra_witness_groups"](sheet, ["g1", "g2"]) == []

    @pytest.mark.parametrize("blank", [None, "", "  "])
    def test_blank_group_cell_is_refused_naming_the_samples(self, blank):
        sheet = _sheet(times=["a"] * 3, groups=["g1", blank, "g2"])
        with pytest.raises(ValueError, match="blank group.*s1"):
            _load()["extra_witness_groups"](sheet, ["g1", "g2"])

    @pytest.mark.parametrize("label", ["old mice", "a;b", "a/b", "$x", "é"])
    def test_group_label_that_is_not_one_plain_word_is_refused(self, label):
        # An extra group is handed to the metadata command as one word of
        # --groups; a space would split it and a shell character would be read
        # by the shell.
        sheet = _sheet(times=["a"] * 2, groups=["g1", label])
        with pytest.raises(ValueError, match="group label"):
            _load()["extra_witness_groups"](sheet, ["g1"])

    def test_underscore_dot_and_dash_are_fine_in_an_extra_group(self):
        sheet = _sheet(times=["a"] * 2, groups=["g1", "old_mice-2.0"])
        assert _load()["extra_witness_groups"](sheet, ["g1"]) == ["old_mice-2.0"]


class TestSheetFullyDated:
    def test_true_when_every_row_has_a_date(self):
        sheet = _sheet(times=["a", "b"], dates=["2020-01-01", "2020-02-01"])
        assert _load()["sheet_fully_dated"](sheet) is True

    def test_false_without_the_column(self):
        assert _load()["sheet_fully_dated"](_sheet(times=["a", "b"])) is False

    @pytest.mark.parametrize("blank", [None, "", "   "])
    def test_false_on_a_blank_cell(self, blank):
        sheet = _sheet(times=["a", "b"], dates=["2020-01-01", blank])
        assert _load()["sheet_fully_dated"](sheet) is False
