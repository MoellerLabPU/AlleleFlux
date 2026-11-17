"""Unit tests for the figure/report config and target helpers.

``figures_params`` and ``groups_for_label`` live in ``common.smk``;
``figure_sentinel``, ``run_report_path`` and ``generate_figure_targets`` in
``dynamic_targets.smk``.  As in ``test_site_depth_filter_params.py`` the
functions are sliced out of the ``.smk`` source by indentation and exec'd with
their free globals stubbed, so no Snakemake session is needed.

Example config: groups_combinations lists treatment ``high_fat`` vs
``control``, whose label is ``high_fat_control`` -- an underscore inside a
group name, which a naive split on "_" would get wrong.
"""

from __future__ import annotations

import os
from pathlib import Path

PIPE = (
    Path(__file__).resolve().parents[1]
    / "alleleflux"
    / "smk_workflow"
    / "alleleflux_pipeline"
)
COMMON_SMK = PIPE / "shared" / "common.smk"
TARGETS_SMK = PIPE / "shared" / "dynamic_targets.smk"


def _slice_function(lines, name):
    start = next(i for i, ln in enumerate(lines) if ln.startswith(f"def {name}("))
    end = start + 1
    while end < len(lines):
        ln = lines[end]
        if ln.strip() and not ln[0].isspace():
            break
        end += 1
    return "\n".join(lines[start:end])


def _load(path, funcs, namespace):
    lines = path.read_text().splitlines()
    src = "\n\n".join(_slice_function(lines, f) for f in funcs)
    exec(compile(src, str(path), "exec"), namespace)
    return {f: namespace[f] for f in funcs}


GROUPS = [{"treatment": "high_fat", "control": "control"}]


def _common(config):
    return _load(COMMON_SMK, ("figures_params", "groups_for_label"), {"config": config})


class TestFiguresParams:
    def test_defaults(self):
        fns = _common({"analysis": {"groups_combinations": GROUPS}})
        assert fns["figures_params"]() == {
            "enabled": True,
            "families": [],
            "tree_path": "",
            "format": "png",
            "min_mags_per_phylum": 3,
        }

    def test_overrides(self):
        cfg = {
            "analysis": {
                "groups_combinations": GROUPS,
                "figures": {"enabled": False, "families": ["ttest"], "format": "both"},
            }
        }
        p = _common(cfg)["figures_params"]()
        assert p["enabled"] is False
        assert p["families"] == ["ttest"]
        assert p["format"] == "both"
        assert p["min_mags_per_phylum"] == 3  # untouched default


class TestGroupsForLabel:
    def test_label_with_underscore_in_a_group_name(self):
        fns = _common({"analysis": {"groups_combinations": GROUPS}})
        assert fns["groups_for_label"]("high_fat_control") == ("high_fat", "control")

    def test_unknown_label_raises(self):
        fns = _common({"analysis": {"groups_combinations": GROUPS}})
        try:
            fns["groups_for_label"]("fat_control")
        except KeyError as e:
            assert "fat_control" in str(e)
        else:
            raise AssertionError("no KeyError")


def _targets(outdir, enabled, allele_only, scores_by_combo, timepoints, groups):
    ns = {
        "os": os,
        "OUTDIR": str(outdir),
        "FIGURES": {"enabled": enabled},
        "config": {"analysis": {"allele_analysis_only": allele_only}},
        "timepoints_labels": timepoints,
        "groups_labels": groups,
        "generate_mag_scores_targets": lambda tp, gr: scores_by_combo.get((tp, gr), []),
    }
    return _load(
        TARGETS_SMK,
        ("figure_sentinel", "run_report_path", "generate_figure_targets"),
        ns,
    )


class TestFigureTargets:
    def test_paths(self, tmp_path):
        fns = _targets(tmp_path, True, False, {}, [], [])
        assert fns["figure_sentinel"](
            str(tmp_path), "pre_end", "high_fat_control"
        ) == os.path.join(
            tmp_path, "figures", "pre_end-high_fat_control", "score_figures.done"
        )
        assert fns["run_report_path"](str(tmp_path)) == os.path.join(
            tmp_path, "figures", "run_report.html"
        )

    def test_one_sentinel_per_combo_with_scores_plus_the_report(self, tmp_path):
        scores = {("pre_end", "high_fat_control"): ["a.tsv"]}  # pre_post has no scores
        fns = _targets(
            tmp_path, True, False, scores, ["pre_end", "pre_post"], ["high_fat_control"]
        )
        assert fns["generate_figure_targets"]() == [
            fns["figure_sentinel"](str(tmp_path), "pre_end", "high_fat_control"),
            fns["run_report_path"](str(tmp_path)),
        ]

    def test_no_combo_with_scores_means_nothing(self, tmp_path):
        fns = _targets(tmp_path, True, False, {}, ["pre_end"], ["high_fat_control"])
        assert fns["generate_figure_targets"]() == []

    def test_disabled_or_allele_only_means_nothing(self, tmp_path):
        scores = {("pre_end", "high_fat_control"): ["a.tsv"]}
        assert (
            _targets(tmp_path, False, False, scores, ["pre_end"], ["high_fat_control"])[
                "generate_figure_targets"
            ]()
            == []
        )
        assert (
            _targets(tmp_path, True, True, scores, ["pre_end"], ["high_fat_control"])[
                "generate_figure_targets"
            ]()
            == []
        )
