"""Which samples of the sheet are read as baseline witnesses under the cohort scope.

``alleleflux-baseline-presence`` and ``alleleflux-gene-de-novo`` both read the
metadata sheet and, under a dated ``cohort`` scope, load samples outside
their comparison as potential baseline witnesses.  The workflow tells them
which timepoints and groups it built caches for (``--witness_timepoints``,
``--witness_groups``); ``cohort_witness_rows`` is the one function that turns
those two lists into "this row is read / is not read", so the two commands
cannot disagree.

Small sheet used below: the comparison is groups A and B at timepoints t0 and
t1; the sheet also holds a group C and a timepoint tx.
"""

import pandas as pd
import pytest

from alleleflux.scripts.utilities.sample_metadata import (
    checked_pre_treatment_timepoints,
    cohort_witness_rows,
)


def _sheet():
    rows = [
        ("a_t0", "A", "t0"),
        ("a_t1", "A", "t1"),
        ("b_t0", "B", "t0"),
        ("a_tx", "A", "tx"),  # compared group, timepoint in no comparison
        ("c_t0", "C", "t0"),  # compared timepoint, group in no comparison
        ("c_tx", "C", "tx"),  # neither
    ]
    meta = pd.DataFrame(rows, columns=["sample_id", "group", "time"])
    in_comparison = meta["group"].isin(["A", "B"]) & meta["time"].isin(["t0", "t1"])
    return meta, in_comparison


def _read(meta, mask):
    return sorted(meta.loc[mask, "sample_id"])


def test_no_list_means_every_sample_of_the_sheet():
    meta, in_comparison = _sheet()
    mask = cohort_witness_rows(meta, in_comparison)
    assert mask.all() and mask.dtype == bool


def test_a_sample_must_be_in_both_lists_to_testify():
    meta, in_comparison = _sheet()
    mask = cohort_witness_rows(
        meta, in_comparison, witness_timepoints=["t0", "t1"], witness_groups=["A", "B"]
    )
    assert _read(meta, mask) == ["a_t0", "a_t1", "b_t0"]
    # every timepoint, compared groups only -> a_tx comes in, the C samples do not
    mask = cohort_witness_rows(
        meta,
        in_comparison,
        witness_timepoints=["t0", "t1", "tx"],
        witness_groups=["A", "B"],
    )
    assert _read(meta, mask) == ["a_t0", "a_t1", "a_tx", "b_t0"]
    # every group, compared timepoints only -> c_t0 comes in, the tx samples do not
    mask = cohort_witness_rows(
        meta,
        in_comparison,
        witness_timepoints=["t0", "t1"],
        witness_groups=["A", "B", "C"],
    )
    assert _read(meta, mask) == ["a_t0", "a_t1", "b_t0", "c_t0"]


def test_one_list_alone_filters_on_that_column_only():
    meta, in_comparison = _sheet()
    mask = cohort_witness_rows(meta, in_comparison, witness_timepoints=["t0", "t1"])
    assert _read(meta, mask) == ["a_t0", "a_t1", "b_t0", "c_t0"]
    mask = cohort_witness_rows(meta, in_comparison, witness_groups=["C"])
    assert _read(meta, mask) == ["a_t0", "a_t1", "b_t0", "c_t0", "c_tx"]


def test_the_comparisons_own_samples_are_always_read():
    meta, in_comparison = _sheet()
    mask = cohort_witness_rows(
        meta, in_comparison, witness_timepoints=["tx"], witness_groups=["C"]
    )
    assert _read(meta, mask) == ["a_t0", "a_t1", "b_t0", "c_tx"]


def test_labels_are_compared_as_text():
    # The sheet is read as text; a list may arrive with numbers in it.
    meta = pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "s3"],
            "group": ["20", "40", "60"],
            "time": ["5"] * 3,
        }
    )
    none_in_comparison = pd.Series(False, index=meta.index)
    mask = cohort_witness_rows(
        meta, none_in_comparison, witness_timepoints=[5], witness_groups=[20, "40"]
    )
    assert _read(meta, mask) == ["s1", "s2"]


def test_the_mask_follows_the_sheets_own_index():
    # A sheet that was filtered before (non-contiguous index) keeps its labels.
    meta, in_comparison = _sheet()
    meta, in_comparison = meta.iloc[[5, 0, 3]], in_comparison.iloc[[5, 0, 3]]
    mask = cohort_witness_rows(
        meta, in_comparison, witness_timepoints=["t0"], witness_groups=["A"]
    )
    assert list(mask.index) == [5, 0, 3]
    assert _read(meta, mask) == ["a_t0"]


# ---------------------------------------------------------------------------
# --pre_treatment_timepoints: which timepoints were sampled before treatment
# began.  The two commands share one check, so they refuse the same mistakes.
# Read timepoints below: t0 and t1 (the comparison) and tx (a witness-only one).
# ---------------------------------------------------------------------------
READ = ["t0", "t1", "tx"]


def test_no_flag_means_no_split():
    """Flag left out -> None, under any scope and on any sheet: nothing to check."""
    assert checked_pre_treatment_timepoints(None, "cohort", True, READ) is None
    assert checked_pre_treatment_timepoints(None, "comparison", False, READ) is None


def test_labels_come_back_as_text_once_each_in_the_order_given():
    got = checked_pre_treatment_timepoints(["t0", "tx", "t0"], "cohort", True, READ)
    assert got == ["t0", "tx"]
    # a label typed as a number matches the text cell of the sheet
    assert checked_pre_treatment_timepoints([5], "cohort", True, ["5", "10"]) == ["5"]


@pytest.mark.parametrize("scope", ["individual", "own_group", "comparison"])
def test_the_flag_under_another_scope_is_refused(scope):
    """Only the cohort scope reads witnesses sampled after treatment began."""
    with pytest.raises(ValueError, match=f"pre_treatment_timepoints.*cohort.*{scope}"):
        checked_pre_treatment_timepoints(["t0"], scope, True, READ)


def test_the_flag_on_an_undated_sheet_is_refused():
    """Without dates the cohort baseline is the earlier timepoint only."""
    with pytest.raises(ValueError, match="collection_date"):
        checked_pre_treatment_timepoints(["t0"], "cohort", False, READ)


def test_a_label_no_read_sample_has_is_refused_naming_it():
    with pytest.raises(ValueError, match=r"\['t9'\].*t0.*t1.*tx"):
        checked_pre_treatment_timepoints(["t0", "t9"], "cohort", True, READ)
