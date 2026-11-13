"""Unit tests for the per-timepoint metadata / QC path helpers and their rules.

Metadata and QC run once per single timepoint.  Everything that needs a QC
table or folder asks one of the small helpers in ``common.smk`` for its path,
so the layout is spelled in one place:

    inputMetadata/inputMetadata_{timepoint}/.metadata_done
    QC/QC_{timepoint}/.qc_done
    QC/QC_{timepoint}/{mag}_QC.tsv

A comparison (``pre_post``) reads the folders of its two timepoints;
``comparison_timepoints`` is the helper that turns the label into them.

The helpers are sliced out of the ``.smk`` source and exec'd with their free
globals stubbed, as in ``test_cache_timepoints.py``.  The last tests read the
rule files as text: they pin that the old "one QC folder per comparison"
lookup is gone and that every flag a rule passes exists in the command it
calls (a renamed flag would otherwise surface only when a job runs).
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PIPELINE = REPO / "alleleflux" / "smk_workflow" / "alleleflux_pipeline"
COMMON_SMK = PIPELINE / "shared" / "common.smk"
RULES = PIPELINE / "rules"
SCRIPTS = REPO / "alleleflux" / "scripts"

HELPERS = (
    "metadata_sentinel",
    "qc_sentinel",
    "qc_dir_of",
    "comparison_timepoints",
    "get_timepoint_qc_file",
    "get_all_qc_files_for_mag",
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


def _load(data_type="longitudinal", unique_timepoints=("pre", "post")):
    """Exec the helpers with a stub run folder ``RUN`` and hand them back by name."""
    lines = COMMON_SMK.read_text().splitlines()
    src = "\n\n".join(_slice_function(lines, f) for f in HELPERS)
    namespace = {
        "os": os,
        # The folder QC and metadata live under (the run's own folder, or the
        # reused run's folder for a permuted run).
        "REUSE_DIR": "RUN",
        "SENTINEL_METADATA": ".metadata_done",
        "SENTINEL_QC": ".qc_done",
        "DATA_TYPE": data_type,
        "unique_timepoints": list(unique_timepoints),
    }
    exec(compile(src, str(COMMON_SMK), "exec"), namespace)
    return {f: namespace[f] for f in HELPERS}


class TestPerTimepointPaths:
    def test_metadata_and_qc_markers_are_per_single_timepoint(self):
        fn = _load()
        assert (
            fn["metadata_sentinel"]("t9")
            == "RUN/inputMetadata/inputMetadata_t9/.metadata_done"
        )
        assert fn["qc_sentinel"]("t9") == "RUN/QC/QC_t9/.qc_done"

    def test_qc_folder_and_table_of_a_timepoint(self):
        fn = _load()
        assert fn["qc_dir_of"]("t9") == "RUN/QC/QC_t9"
        assert fn["get_timepoint_qc_file"]("MAG_7", "t9") == "RUN/QC/QC_t9/MAG_7_QC.tsv"
        # The marker sits in the folder the helper names.
        assert os.path.dirname(fn["qc_sentinel"]("t9")) == fn["qc_dir_of"]("t9")

    def test_wildcard_patterns_pass_through(self):
        # Rules call the helpers with the wildcard itself to build a pattern.
        fn = _load()
        assert fn["qc_sentinel"]("{timepoint}") == "RUN/QC/QC_{timepoint}/.qc_done"
        assert (
            fn["get_timepoint_qc_file"]("{mag}", "{timepoint}")
            == "RUN/QC/QC_{timepoint}/{mag}_QC.tsv"
        )


class TestComparisonTimepoints:
    def test_longitudinal_label_is_its_two_timepoints_in_order(self):
        assert _load("longitudinal")["comparison_timepoints"]("pre_post") == [
            "pre",
            "post",
        ]

    def test_single_label_is_the_timepoint_itself(self):
        assert _load("single")["comparison_timepoints"]("t1") == ["t1"]

    @pytest.mark.parametrize("label", ["pre", "a_b_c"])
    def test_longitudinal_label_that_is_not_two_timepoints_is_refused(self, label):
        # Timepoint labels never hold "_" (parse time refuses it), so anything
        # but exactly two parts is a wiring mistake, not data.
        with pytest.raises(ValueError, match=label):
            _load("longitudinal")["comparison_timepoints"](label)


class TestQcTablesOfOneGenome:
    def test_one_table_per_comparison_timepoint_in_config_order(self):
        # Pairwise ANI unions the samples that passed QC at any comparison
        # timepoint: one table per timepoint, each sample in exactly one.
        fn = _load(unique_timepoints=("pre", "post", "end"))
        assert fn["get_all_qc_files_for_mag"]("MAG_7") == [
            "RUN/QC/QC_pre/MAG_7_QC.tsv",
            "RUN/QC/QC_post/MAG_7_QC.tsv",
            "RUN/QC/QC_end/MAG_7_QC.tsv",
        ]


class TestOldLayoutIsGone:
    def test_canonical_comparison_lookup_is_removed_everywhere(self):
        # QC used to live in one folder per COMPARISON, and a timepoint's cache
        # read the table of the first comparison naming it.  With one folder
        # per timepoint there is nothing to look up.
        sources = [COMMON_SMK, PIPELINE / "Snakefile", *sorted(RULES.glob("*.smk"))]
        sources += sorted((PIPELINE / "shared").glob("*.smk"))
        for path in sources:
            text = path.read_text()
            assert "timepoint_to_canonical_tp" not in text, path.name
            assert "get_canonical_qc_file" not in text, path.name

    def test_metadata_and_qc_rules_are_keyed_on_one_timepoint(self):
        metadata = (RULES / "metadata.smk").read_text()
        qc = (RULES / "quality_control.smk").read_text()
        assert 'metadata_sentinel("{timepoint}")' in metadata
        assert 'metadata_sentinel("{timepoint}")' in qc
        assert 'qc_sentinel("{timepoint}")' in qc
        # No comparison-label wildcard is left in either rule.
        assert "{timepoints}" not in metadata
        assert "{timepoints}" not in qc


def _shell_flags(rule_file, command):
    """The ``--flags`` a rule's shell block passes to ``command``."""
    text = (RULES / rule_file).read_text()
    block = text[text.index(command) :]
    block = block[: block.index('"""')]
    return set(re.findall(r"(--[a-zA-Z_]+)", block))


def _cli_flags(script):
    """The ``--flags`` a script's argument parser defines."""
    text = (SCRIPTS / script).read_text()
    return set(re.findall(r"""["'](--[a-zA-Z_]+)["']""", text))


class TestRuleFlagsExistInTheCommands:
    @pytest.mark.parametrize(
        "rule_file, command, script",
        [
            ("metadata.smk", "alleleflux-metadata", "preprocessing/mag_metadata.py"),
            (
                "quality_control.smk",
                "alleleflux-qc",
                "preprocessing/quality_control.py",
            ),
            (
                "eligibility.smk",
                "alleleflux-eligibility",
                "preprocessing/eligibility_table.py",
            ),
        ],
    )
    def test_every_flag_the_rule_passes_is_defined_by_the_command(
        self, rule_file, command, script
    ):
        passed = _shell_flags(rule_file, command)
        assert passed, f"no flag found in the shell block of {rule_file}"
        assert passed <= _cli_flags(script), passed - _cli_flags(script)

    def test_eligibility_rule_passes_the_per_timepoint_folders(self):
        flags = _shell_flags("eligibility.smk", "alleleflux-eligibility")
        assert "--qc_dirs" in flags
        assert "--qc_dir" not in flags
