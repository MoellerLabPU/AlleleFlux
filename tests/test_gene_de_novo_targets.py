"""Unit tests for ``generate_gene_de_novo_targets`` in ``dynamic_targets.smk``.

The gene de novo rule runs once per comparison and reads every QC-eligible
MAG's mean-change table, so its one target is requested only when the
feature is on, the data is longitudinal and at least one MAG passed QC for
the comparison.  Sliced out of the ``.smk`` source and exec'd with its free
globals stubbed (``config``, ``DATA_TYPE``, ``_get_mags_by_eligibility``,
``get_gene_de_novo_tests_path``), as in ``test_rising_origin_targets.py``.
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
FUNC = "generate_gene_de_novo_targets"


def _slice_function(lines, name):
    start = next(i for i, ln in enumerate(lines) if ln.startswith(f"def {name}("))
    end = start + 1
    while end < len(lines):
        ln = lines[end]
        if ln.strip() and not ln[0].isspace():
            break
        end += 1
    return "\n".join(lines[start:end])


def _load(flag, data_type, eligible, allele_analysis_only=False):
    lines = DYNAMIC_TARGETS_SMK.read_text().splitlines()
    namespace = {
        "os": os,
        "config": {
            "analysis": {
                "use_gene_de_novo": flag,
                "allele_analysis_only": allele_analysis_only,
            }
        },
        "DATA_TYPE": data_type,
        "_get_mags_by_eligibility": lambda tp, gr, eligibility_type: list(eligible),
        "get_gene_de_novo_tests_path": lambda timepoints, groups: (
            f"OUT/gene_de_novo/{timepoints}-{groups}/gene_de_novo_tests_{timepoints}-{groups}.tsv"
        ),
    }
    exec(
        compile(_slice_function(lines, FUNC), str(DYNAMIC_TARGETS_SMK), "exec"),
        namespace,
    )
    return namespace[FUNC]


class TestGeneDeNovoTargets:
    def test_one_target_per_comparison_when_on_and_eligible(self):
        gen = _load(True, "longitudinal", ["MAG_1"])
        assert gen("pre_post", "a_b") == [
            "OUT/gene_de_novo/pre_post-a_b/gene_de_novo_tests_pre_post-a_b.tsv"
        ]

    def test_flag_off_means_nothing(self):
        assert _load(False, "longitudinal", ["MAG_1"])("pre_post", "a_b") == []

    def test_single_data_means_nothing(self):
        assert _load(True, "single", ["MAG_1"])("t1", "a_b") == []

    def test_allele_analysis_only_means_nothing(self):
        gen = _load(True, "longitudinal", ["MAG_1"], allele_analysis_only=True)
        assert gen("pre_post", "a_b") == []

    def test_no_eligible_mag_means_nothing(self):
        assert _load(True, "longitudinal", [])("pre_post", "a_b") == []
