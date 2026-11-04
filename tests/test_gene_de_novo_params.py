"""Unit tests for the gene de novo config helpers in ``common.smk``.

``analysis.gene_de_novo`` holds the three knobs that belong to the gene-level
analysis alone (``block_rule``, ``min_gene_covered_fraction``,
``min_determined_carriers``); the shared de novo rules (baseline scope,
cluster window, strain gate) come from ``analysis.baseline_presence`` through
``de_novo_rule_params``.  ``gene_de_novo_params`` reads them once with the
house defaults and refuses what the cage-level analysis cannot do (a group
strain vote); ``gene_de_novo_cli_args`` turns them into the rule's flag
string.  Sliced out of the ``.smk`` source and exec'd with a stub ``config``,
as in ``test_de_novo_rule_params.py``.
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
FUNCS = ("de_novo_rule_params", "gene_de_novo_params", "gene_de_novo_cli_args")


def _slice_function(lines, name):
    start = next(i for i, ln in enumerate(lines) if ln.startswith(f"def {name}("))
    end = start + 1
    while end < len(lines):
        ln = lines[end]
        if ln.strip() and not ln[0].isspace():
            break
        end += 1
    return "\n".join(lines[start:end])


def _load(config):
    lines = COMMON_SMK.read_text().splitlines()
    src = "\n\n".join(_slice_function(lines, f) for f in FUNCS)
    namespace = {"config": config}
    exec(compile(src, str(COMMON_SMK), "exec"), namespace)
    return {f: namespace[f] for f in FUNCS}


def _cfg(block=None, presence=None, turnover=True):
    analysis = {"use_strain_turnover": turnover, "use_gene_de_novo": True}
    if presence is not None:
        analysis["baseline_presence"] = presence
    if block is not None:
        analysis["gene_de_novo"] = block
    return {"analysis": analysis}


class TestGeneDeNovoParams:
    def test_defaults_without_a_block(self):
        fns = _load(_cfg())
        assert fns["gene_de_novo_params"]() == {
            "baseline_scope": "cohort",
            "cluster_window_bp": 1000,
            "max_de_novo_neighbours": 0,
            "strain_gate": "carrier",
            "min_determined_carriers": 1,  # the gene analysis's own floor, not 3
            "block_rule": "exclude",
            "min_gene_covered_fraction": 0.0,
        }

    def test_shared_rules_come_from_the_baseline_presence_block(self):
        fns = _load(
            _cfg(
                block={"block_rule": "once", "min_gene_covered_fraction": 0.5},
                presence={
                    "baseline_scope": "comparison",
                    "cluster_window_bp": 500,
                    "max_de_novo_neighbours": 2,
                    "strain_gate": "none",
                    "min_determined_carriers": 3,  # baseline presence's floor
                },
                turnover=False,
            )
        )
        p = fns["gene_de_novo_params"]()
        assert (
            p["baseline_scope"],
            p["cluster_window_bp"],
            p["max_de_novo_neighbours"],
        ) == (
            "comparison",
            500,
            2,
        )
        assert p["strain_gate"] == "none"
        # the gene block's own floor, defaulting to 1, NOT the presence block's 3
        assert p["min_determined_carriers"] == 1
        assert (p["block_rule"], p["min_gene_covered_fraction"]) == ("once", 0.5)

    @pytest.mark.parametrize("gate", ["group", "both"])
    def test_group_strain_votes_are_refused_naming_the_block(self, gate):
        fns = _load(_cfg(presence={"strain_gate": gate}))
        with pytest.raises(ValueError, match="gene_de_novo"):
            fns["gene_de_novo_params"]()

    def test_unknown_block_rule_is_refused(self):
        fns = _load(_cfg(block={"block_rule": "twice"}))
        with pytest.raises(ValueError, match="block_rule"):
            fns["gene_de_novo_params"]()

    @pytest.mark.parametrize("fraction", [-0.1, 1.5])
    def test_covered_fraction_outside_zero_one_is_refused(self, fraction):
        fns = _load(_cfg(block={"min_gene_covered_fraction": fraction}))
        with pytest.raises(ValueError, match="min_gene_covered_fraction"):
            fns["gene_de_novo_params"]()

    def test_carrier_floor_below_one_is_refused(self):
        fns = _load(_cfg(block={"min_determined_carriers": 0}))
        with pytest.raises(ValueError, match="min_determined_carriers"):
            fns["gene_de_novo_params"]()

    def test_carrier_gate_without_turnover_is_refused_by_the_shared_check(self):
        fns = _load(_cfg(turnover=False))
        with pytest.raises(ValueError, match="use_strain_turnover"):
            fns["gene_de_novo_params"]()


class TestGeneDeNovoCliArgs:
    def test_flag_string_carries_every_knob(self):
        fns = _load(_cfg(block={"block_rule": "once"}))
        assert fns["gene_de_novo_cli_args"]() == (
            "--baseline_scope cohort --cluster_window_bp 1000 "
            "--max_de_novo_neighbours 0 --strain_gate carrier "
            "--min_determined_carriers 1 --block_rule once "
            "--min_gene_covered_fraction 0.0"
        )
