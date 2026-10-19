"""Unit tests for ``switch_off_label_blind_extras`` in ``common.smk``.

A permuted (null) run must not run dN/dS, baseline presence, pairwise ANI or
strain turnover: they read labels from the real metadata / QC, never from the
permuted sheet.  The helper and its constant are sliced out of the ``.smk``
source and exec'd with a stub ``config``, as in ``test_site_depth_filter_params.py``.
"""

from __future__ import annotations

from pathlib import Path

COMMON_SMK = (
    Path(__file__).resolve().parents[1]
    / "alleleflux"
    / "smk_workflow"
    / "alleleflux_pipeline"
    / "shared"
    / "common.smk"
)
FLAGS = ("use_dnds", "use_baseline_presence", "use_pairwise_ani", "use_strain_turnover")


def _slice(lines, start_marker):
    """Lines from ``start_marker`` up to the next top-level statement."""
    start = next(i for i, ln in enumerate(lines) if ln.startswith(start_marker))
    end = start + 1
    while end < len(lines):
        ln = lines[end]
        if ln.strip() and not ln[0].isspace() and not ln.startswith(")"):
            break
        end += 1
    return "\n".join(lines[start:end])


def _load():
    lines = COMMON_SMK.read_text().splitlines()
    src = "\n\n".join(
        [
            _slice(lines, "LABEL_BLIND_EXTRAS = ("),
            _slice(lines, "def switch_off_label_blind_extras("),
        ]
    )
    namespace: dict = {}
    exec(compile(src, str(COMMON_SMK), "exec"), namespace)
    return namespace["switch_off_label_blind_extras"]


class TestSwitchOffLabelBlindExtras:
    def test_permuted_run_forces_all_four_off_and_reports_what_was_on(self):
        fn = _load()
        # use_dnds absent (defaults ON), ANI explicitly on, the other two off
        config = {
            "analysis": {"use_pairwise_ani": True, "use_baseline_presence": False}
        }
        was_on = fn(config, is_permuted_run=True)
        assert was_on == ["use_dnds", "use_pairwise_ani"]
        assert all(config["analysis"][f] is False for f in FLAGS)

    def test_normal_run_is_left_untouched(self):
        fn = _load()
        config = {"analysis": {"use_pairwise_ani": True}}
        assert fn(config, is_permuted_run=False) == []
        assert config == {"analysis": {"use_pairwise_ani": True}}  # no keys added

    def test_permuted_run_with_everything_already_off_reports_nothing(self):
        fn = _load()
        config = {"analysis": {f: False for f in FLAGS}}
        assert fn(config, is_permuted_run=True) == []
        assert all(config["analysis"][f] is False for f in FLAGS)
