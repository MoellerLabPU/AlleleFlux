"""Unit tests for the parse-time timepoint-order helpers in ``common.smk``.

Every longitudinal step of the pipeline reads the two timepoints of a
``timepoints_combinations`` entry as ``earlier`` and ``later`` (the config
spells them so since 2026-11-02).  ``check_timepoint_order`` compares that
order with the metadata sheet's ``collection_date`` column and refuses the
combination as soon as ONE mouse has its "later" sample dated before its
"earlier" one (user ruling 2026-11-02): a swapped pair or a wrong date stops
before any job instead of flipping a sign.

Sliced out of the ``.smk`` source and exec'd alone, as in
``test_data_type_checks.py``; pandas is handed in because the check uses it.
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
FUNCS = ("check_timepoint_order",)


def _slice_function(lines, name):
    start = next(i for i, ln in enumerate(lines) if ln.startswith(f"def {name}("))
    end = start + 1
    while end < len(lines):
        ln = lines[end]
        if ln.strip() and not ln[0].isspace():
            break
        end += 1
    return "\n".join(lines[start:end])


def _load():
    lines = COMMON_SMK.read_text().splitlines()
    src = "\n\n".join(_slice_function(lines, f) for f in FUNCS)
    namespace = {"pd": pd}
    exec(compile(src, str(COMMON_SMK), "exec"), namespace)
    return {f: namespace[f] for f in FUNCS}


def _sheet(rows):
    """rows: (subjectID, time, collection_date) -> the columns the check reads."""
    return pd.DataFrame(rows, columns=["subjectID", "time", "collection_date"])


class TestCheckTimepointOrder:
    def test_paired_mice_dated_forward_pass_and_are_counted(self):
        meta = _sheet(
            [
                ("m1", "pre", "2024-01-15"),
                ("m1", "post", "2024-04-15"),
                ("m2", "pre", "2024-01-15"),
                ("m2", "post", "2024-04-15"),
            ]
        )
        out = _load()["check_timepoint_order"](meta, "pre", "post")
        assert out == {"checked": True, "paired": 2}

    def test_swapped_order_is_refused_naming_both_timepoints(self):
        meta = _sheet(
            [
                ("m1", "pre", "2024-01-15"),
                ("m1", "post", "2024-04-15"),
            ]
        )
        with pytest.raises(ValueError, match="post.*pre|pre.*post"):
            _load()["check_timepoint_order"](meta, "post", "pre")

    def test_a_single_backward_mouse_stops_the_run_naming_it(self):
        meta = _sheet(
            [
                ("m1", "pre", "2024-01-15"),
                ("m1", "post", "2024-04-15"),
                ("m2", "pre", "2024-01-15"),
                ("m2", "post", "2024-04-15"),
                (
                    "m3",
                    "pre",
                    "2024-05-01",
                ),  # this mouse's "pre" is dated after its "post"
                ("m3", "post", "2024-04-15"),
            ]
        )
        with pytest.raises(ValueError, match="m3") as exc:
            _load()["check_timepoint_order"](meta, "pre", "post")
        assert "1 of 3" in str(exc.value)

    def test_no_mouse_sampled_at_both_timepoints_is_refused(self):
        # longitudinal AlleleFlux pairs each subject's two samples (the
        # allele-frequency step inner-joins on subjectID), so a sheet with no
        # subject at both timepoints has nothing to analyse
        meta = _sheet(
            [
                ("m1", "pre", "2024-01-15"),
                ("m2", "pre", "2024-02-15"),
                ("m3", "post", "2024-04-15"),
                ("m4", "post", "2024-05-15"),
            ]
        )
        with pytest.raises(ValueError, match="no subject") as exc:
            _load()["check_timepoint_order"](meta, "pre", "post")
        assert "pre" in str(exc.value) and "post" in str(exc.value)

    def test_sheet_without_dates_is_not_checked(self):
        meta = pd.DataFrame({"subjectID": ["m1", "m1"], "time": ["pre", "post"]})
        out = _load()["check_timepoint_order"](meta, "pre", "post")
        assert out == {"checked": False, "paired": 0}

    def test_a_blank_date_among_the_two_timepoints_skips_the_check(self):
        meta = _sheet(
            [
                ("m1", "pre", "2024-01-15"),
                ("m1", "post", ""),
                ("m9", "mid", "2024-03-01"),  # other timepoints are not consulted
            ]
        )
        out = _load()["check_timepoint_order"](meta, "pre", "post")
        assert out["checked"] is False

    def test_a_malformed_date_is_refused_by_value(self):
        meta = _sheet([("m1", "pre", "15/01/2024"), ("m1", "post", "2024-04-15")])
        with pytest.raises(ValueError, match="15/01/2024"):
            _load()["check_timepoint_order"](meta, "pre", "post")
