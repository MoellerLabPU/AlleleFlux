"""``transition_labels`` in ``common.smk``: the EARLIER:LATER list that the
pairwise-ANI rule (``pairs: transitions``) and the strain-turnover rule pass on
the command line, built from ``analysis.timepoints_combinations`` in the ONE
spelling the workflow accepts (``earlier:`` / ``later:``).

Found by the first cohort dry run after the spelling change (2026-11-04): the
pairwise-ANI rule still read ``tc["timepoint"][0]`` and failed at parse time
with ``pairs: transitions``; the bundled example uses ``within_subject`` and
never took that branch.  One helper, two callers, one test.

The helper is sliced out of the ``.smk`` source and exec'd alone, as in
``test_timepoint_order.py``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

COMMON_SMK = (
    Path(__file__).resolve().parents[1]
    / "alleleflux"
    / "smk_workflow"
    / "alleleflux_pipeline"
    / "shared"
    / "common.smk"
)


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
    namespace: dict = {}
    exec(
        compile(_slice_function(lines, "transition_labels"), str(COMMON_SMK), "exec"),
        namespace,
    )
    return namespace["transition_labels"]


class TestTransitionLabels:
    def test_longitudinal_entries_become_earlier_colon_later_in_order(self):
        fn = _load()
        combos = [
            {"earlier": "5mo", "later": "10mo"},
            {"earlier": "5mo", "later": "16mo"},
        ]
        assert fn(combos) == "5mo:10mo 5mo:16mo"

    def test_numeric_timepoint_names_are_written_as_given(self):
        fn = _load()
        assert fn([{"earlier": 0, "later": 7}]) == "0:7"

    def test_single_timepoint_entries_have_no_transition(self):
        fn = _load()
        assert fn([{"timepoint": "t1"}]) == ""

    def test_the_old_spelling_is_not_silently_accepted(self):
        # The parse-time shape check refuses it first; the helper must not
        # paper over a stray old-form entry that slipped past it.
        fn = _load()
        with pytest.raises(KeyError):
            fn([{"timepoint": ["pre", "post"], "focus": "post", "later": "post"}])
