"""Unit tests for the de novo rule config helpers in ``common.smk``.

``analysis.baseline_presence`` holds the knobs of the de novo analysis (the
baseline scope, the cluster rule, the strain gate); ``de_novo_rule_params``
reads them once with the house defaults and checks them against the rest of
the config, and ``de_novo_rule_cli_args`` turns the rising-origin ones into
the flag string both rules share.  Sliced out of the ``.smk`` source and
exec'd with a stub ``config``, as in ``test_site_depth_filter_params.py``.
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
FUNCS = ("de_novo_rule_params", "de_novo_rule_cli_args")


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


def _cfg(block=None, turnover=True):
    analysis = {"use_strain_turnover": turnover}
    if block is not None:
        analysis["baseline_presence"] = block
    return {"analysis": analysis}


class TestDeNovoRuleParams:
    def test_defaults_when_the_block_has_none_of_the_keys(self):
        fns = _load(_cfg({"summary": "two_sample_paired"}))
        assert fns["de_novo_rule_params"]() == {
            "baseline_scope": "cohort",
            "cluster_window_bp": 1000,
            "max_de_novo_neighbours": 0,
            "strain_gate": "carrier",
            "min_determined_carriers": 3,
        }
        assert fns["de_novo_rule_cli_args"]() == (
            "--cluster_window_bp 1000 --max_de_novo_neighbours 0 "
            "--strain_gate carrier --min_determined_carriers 3"
        )

    def test_no_block_at_all(self):
        fns = _load(_cfg(None))
        assert fns["de_novo_rule_params"]()["strain_gate"] == "carrier"

    def test_values_are_read_and_typed(self):
        fns = _load(
            _cfg(
                {
                    "baseline_scope": "individual",
                    "cluster_window_bp": "500",
                    "max_de_novo_neighbours": 2,
                    "strain_gate": "both",
                    "min_determined_carriers": "1",
                }
            )
        )
        p = fns["de_novo_rule_params"]()
        assert (
            p["cluster_window_bp"],
            p["max_de_novo_neighbours"],
            p["min_determined_carriers"],
        ) == (500, 2, 1)
        assert fns["de_novo_rule_cli_args"]() == (
            "--cluster_window_bp 500 --max_de_novo_neighbours 2 "
            "--strain_gate both --min_determined_carriers 1"
        )

    def test_carrier_or_both_need_strain_turnover(self):
        # the carrier gate reads the per-mouse turnover tables; without that
        # step the run must stop at parse time, not after the statistics
        for gate in ("carrier", "both"):
            with pytest.raises(ValueError, match="use_strain_turnover"):
                _load(_cfg({"strain_gate": gate}, turnover=False))[
                    "de_novo_rule_params"
                ]()
        _load(_cfg({"strain_gate": "none"}, turnover=False))["de_novo_rule_params"]()

    def test_group_needs_strain_turnover_too(self):
        # the group vote is the replacement-classification table, written by
        # the same feature switch
        with pytest.raises(ValueError, match="use_strain_turnover"):
            _load(_cfg({"strain_gate": "group"}, turnover=False))[
                "de_novo_rule_params"
            ]()

    def test_unknown_words_raise(self):
        with pytest.raises(ValueError, match="strain_gate"):
            _load(_cfg({"strain_gate": "maybe"}))["de_novo_rule_params"]()
        with pytest.raises(ValueError, match="baseline_scope"):
            _load(_cfg({"baseline_scope": "everyone"}))["de_novo_rule_params"]()

    def test_negative_numbers_raise(self):
        with pytest.raises(ValueError, match="cluster_window_bp"):
            _load(_cfg({"cluster_window_bp": -1}))["de_novo_rule_params"]()
        with pytest.raises(ValueError, match="min_determined_carriers"):
            _load(_cfg({"min_determined_carriers": 0}))["de_novo_rule_params"]()
