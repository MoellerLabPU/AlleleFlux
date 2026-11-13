"""Unit tests for the rising-origin target generators in ``dynamic_targets.smk``.

``rising_origin`` is a run-once rollup over every comparison's rising table.  Its
summary must be requested only when at least one comparison actually produces a
rising table; otherwise the rule runs ``alleleflux-rising-origin`` on an empty
folder and the script raises ``FileNotFoundError``.  These tests pin that gate.

The ``.smk`` module body references Snakemake globals, so (as in
``test_eligibility_helpers.py``) the two generators are sliced out of the source
by indentation and exec'd with their free globals stubbed:

* ``generate_baseline_presence_targets(tp, gr)`` -- stubbed per test to say which
  comparisons have a rising table;
* ``timepoints_labels`` / ``groups_labels`` -- the comparison grid;
* ``OUTDIR`` / ``os``.
"""

from __future__ import annotations

import os
from pathlib import Path

DYNAMIC_TARGETS_SMK = (
    Path(__file__).resolve().parents[1]
    / "alleleflux"
    / "smk_workflow"
    / "alleleflux_pipeline"
    / "shared"
    / "dynamic_targets.smk"
)

FUNCS = ("generate_rising_table_targets", "generate_rising_origin_targets")


def _slice_function(lines, name):
    """Source of top-level ``def name`` up to the next column-0 statement."""
    start = next(i for i, ln in enumerate(lines) if ln.startswith(f"def {name}("))
    end = start + 1
    while end < len(lines):
        ln = lines[end]
        if ln.strip() and not ln[0].isspace():
            break
        end += 1
    return "\n".join(lines[start:end])


def _load_generators(outdir, baseline_targets, timepoints, groups):
    """Exec the two generators with stubbed globals; return them as a dict."""
    lines = DYNAMIC_TARGETS_SMK.read_text().splitlines()
    src = "\n\n".join(_slice_function(lines, f) for f in FUNCS)
    namespace = {
        "os": os,
        "OUTDIR": str(outdir),
        "timepoints_labels": timepoints,
        "groups_labels": groups,
        "generate_baseline_presence_targets": baseline_targets,
    }
    exec(compile(src, str(DYNAMIC_TARGETS_SMK), "exec"), namespace)
    return {f: namespace[f] for f in FUNCS}


def _stem(outdir, tp, gr):
    return os.path.join(str(outdir), "baseline_presence", f"{tp}-{gr}_baseline_presence")


def _summary_path(outdir):
    return os.path.join(
        str(outdir), "baseline_presence", "rising_origin", "rising_origin_summary.tsv"
    )


class TestRisingOriginTargets:
    TIMEPOINTS = ["t0_t1", "t0_t2"]
    GROUPS = ["a_b"]

    def test_no_rising_table_anywhere_means_no_summary_target(self, tmp_path):
        """Flag on but every comparison empty (no eligible MAG): nothing requested."""
        gens = _load_generators(
            tmp_path, lambda tp, gr: [], self.TIMEPOINTS, self.GROUPS
        )
        assert gens["generate_rising_table_targets"]() == []
        assert gens["generate_rising_origin_targets"]() == []

    def test_one_comparison_with_rising_table_requests_summary(self, tmp_path):
        """Only t0_t2 has a rising table: it is the sole input, summary requested."""

        def baseline(tp, gr):
            if tp != "t0_t2":
                return []
            stem = _stem(tmp_path, tp, gr)
            return [f"{stem}.tsv.gz", f"{stem}_summary.tsv", f"{stem}_rising.tsv"]

        gens = _load_generators(tmp_path, baseline, self.TIMEPOINTS, self.GROUPS)
        assert gens["generate_rising_table_targets"]() == [
            f"{_stem(tmp_path, 't0_t2', 'a_b')}_rising.tsv"
        ]
        assert gens["generate_rising_origin_targets"]() == [_summary_path(tmp_path)]

    def test_rising_tables_only_never_the_long_or_summary_files(self, tmp_path):
        """The rollup's inputs are the *_rising.tsv files alone."""

        def baseline(tp, gr):
            stem = _stem(tmp_path, tp, gr)
            return [f"{stem}.tsv.gz", f"{stem}_summary.tsv", f"{stem}_rising.tsv"]

        gens = _load_generators(tmp_path, baseline, self.TIMEPOINTS, self.GROUPS)
        rising = gens["generate_rising_table_targets"]()
        assert len(rising) == 2
        assert all(t.endswith("_rising.tsv") for t in rising)
