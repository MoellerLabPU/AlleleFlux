"""Unit tests for the site-depth-filter summary targets in ``dynamic_targets.smk``.

The run-once rollup must be requested exactly when at least one cache sidecar
will exist, and its input must list every sidecar of every QC-eligible
(MAG, timepoint).  Same exec-by-indentation harness as the rising-origin tests:
the two generators are sliced out and run with their free globals stubbed.
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
FUNCS = (
    "generate_allele_freq_cache_sidecar_targets",
    "generate_site_depth_filter_summary_targets",
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


def _load(outdir, eligible, timepoints, groups, data_type="longitudinal"):
    lines = DYNAMIC_TARGETS_SMK.read_text().splitlines()
    src = "\n\n".join(_slice_function(lines, f) for f in FUNCS)

    def cache_path(mag_wildcard, timepoint_wildcard):
        return os.path.join(
            outdir,
            "allele_freq_cache",
            timepoint_wildcard,
            f"{mag_wildcard}_{timepoint_wildcard}_allele_frequency.parquet",
        )

    namespace = {
        "os": os,
        "OUTDIR": outdir,
        "DATA_TYPE": data_type,
        "timepoints_labels": timepoints,
        "groups_labels": groups,
        "_get_mags_by_eligibility": lambda tp, gr, eligibility_type: eligible.get(
            (tp, gr), []
        ),
        "get_allele_freq_cache_path": cache_path,
    }
    exec(compile(src, str(DYNAMIC_TARGETS_SMK), "exec"), namespace)
    return {f: namespace[f] for f in FUNCS}


class TestSidecarTargets:
    def test_union_over_comparisons_and_timepoints(self, tmp_path):
        out = str(tmp_path)
        fns = _load(
            out,
            {
                ("pre_end", "fat_control"): ["MAG_A", "MAG_B"],
                ("pre_post", "fat_control"): ["MAG_A"],
            },
            ["pre_end", "pre_post"],
            ["fat_control"],
        )
        got = fns["generate_allele_freq_cache_sidecar_targets"]()
        expect = sorted(
            os.path.join(
                out,
                "allele_freq_cache",
                tp,
                f"{mag}_{tp}_allele_frequency_site_depth_filter.tsv",
            )
            for mag, tp in [
                ("MAG_A", "pre"),
                ("MAG_A", "end"),
                ("MAG_A", "post"),
                ("MAG_B", "pre"),
                ("MAG_B", "end"),
            ]
        )
        assert got == expect  # de-duplicated: MAG_A at pre appears once

    def test_single_data_uses_the_label_as_the_timepoint(self, tmp_path):
        fns = _load(
            str(tmp_path),
            {("t1", "a_b"): ["MAG_A"]},
            ["t1"],
            ["a_b"],
            data_type="single",
        )
        got = fns["generate_allele_freq_cache_sidecar_targets"]()
        assert got == [
            os.path.join(
                str(tmp_path),
                "allele_freq_cache",
                "t1",
                "MAG_A_t1_allele_frequency_site_depth_filter.tsv",
            )
        ]


class TestSummaryTarget:
    def test_requested_when_a_sidecar_exists(self, tmp_path):
        fns = _load(
            str(tmp_path),
            {("pre_end", "fat_control"): ["MAG_A"]},
            ["pre_end"],
            ["fat_control"],
        )
        assert fns["generate_site_depth_filter_summary_targets"]() == [
            os.path.join(
                str(tmp_path), "site_depth_filter", "site_depth_filter_summary.tsv"
            )
        ]

    def test_not_requested_without_any_cache(self, tmp_path):
        fns = _load(str(tmp_path), {}, ["pre_end"], ["fat_control"])
        assert fns["generate_site_depth_filter_summary_targets"]() == []
