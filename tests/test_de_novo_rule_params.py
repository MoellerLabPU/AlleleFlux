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
            "cohort_witnesses": "sheet",
            "pre_treatment_timepoints": [],
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

    def test_cohort_witnesses_defaults_to_sheet_and_refuses_an_unknown_word(self):
        # Who may testify under the cohort scope: the whole sheet (default) or
        # only samples at a compared timepoint in a compared group.
        assert _load(_cfg({}))["de_novo_rule_params"]()["cohort_witnesses"] == "sheet"
        assert (
            _load(_cfg({"cohort_witnesses": "comparisons"}))["de_novo_rule_params"]()[
                "cohort_witnesses"
            ]
            == "comparisons"
        )
        with pytest.raises(ValueError, match="cohort_witnesses"):
            _load(_cfg({"cohort_witnesses": "everyone"}))["de_novo_rule_params"]()

    def test_cohort_witnesses_is_kept_but_not_refused_under_another_scope(self):
        # The config template always carries the line; under a scope that
        # reads no outside witness it is simply not used.
        p = _load(_cfg({"baseline_scope": "comparison", "cohort_witnesses": "sheet"}))[
            "de_novo_rule_params"
        ]()
        assert (p["baseline_scope"], p["cohort_witnesses"]) == ("comparison", "sheet")

    def test_cohort_witnesses_never_reaches_the_rising_origin_flags(self):
        # It decides which samples are read, upstream; rising_origin has no such flag.
        args = _load(_cfg({"cohort_witnesses": "comparisons"}))[
            "de_novo_rule_cli_args"
        ]()
        assert "cohort_witnesses" not in args and "comparisons" not in args

    def test_pre_treatment_timepoints_default_to_an_empty_list(self):
        # No list = no split: the two labels of before (de novo / standing).
        assert (
            _load(_cfg({}))["de_novo_rule_params"]()["pre_treatment_timepoints"] == []
        )
        assert (
            _load(_cfg(None))["de_novo_rule_params"]()["pre_treatment_timepoints"] == []
        )
        # An empty key in the YAML (``pre_treatment_timepoints:``) reads as None.
        assert (
            _load(_cfg({"pre_treatment_timepoints": None}))["de_novo_rule_params"]()[
                "pre_treatment_timepoints"
            ]
            == []
        )

    def test_pre_treatment_labels_are_kept_as_text(self):
        # The sheet is read as text, so a YAML number 5 and the text "5" are
        # the same timepoint; spaces around a label and a repeated label go.
        p = _load(_cfg({"pre_treatment_timepoints": [5, " t0 ", "5", "t0"]}))[
            "de_novo_rule_params"
        ]()
        assert p["pre_treatment_timepoints"] == ["5", "t0"]

    def test_pre_treatment_timepoints_must_be_a_list(self):
        # ``pre_treatment_timepoints: t0`` is one string; reading it letter by
        # letter would make "t" and "0" two timepoints.
        for not_a_list in ("t0", 5, {"t0": True}):
            with pytest.raises(ValueError, match="pre_treatment_timepoints.*list"):
                _load(_cfg({"pre_treatment_timepoints": not_a_list}))[
                    "de_novo_rule_params"
                ]()

    def test_a_blank_pre_treatment_label_is_refused(self):
        for blank in ("", "  ", None):
            with pytest.raises(ValueError, match="pre_treatment_timepoints.*blank"):
                _load(_cfg({"pre_treatment_timepoints": ["t0", blank]}))[
                    "de_novo_rule_params"
                ]()

    @pytest.mark.parametrize("scope", ["individual", "own_group", "comparison"])
    def test_pre_treatment_timepoints_under_another_scope_stop_the_run(self, scope):
        # Only the cohort scope reads witnesses outside the comparison's
        # earlier timepoint, so only there can a witness be "on treatment".
        # A filled list under any other scope is a mistake, not a leftover
        # (user ruling 2026-11-16): stop, naming both settings.
        with pytest.raises(
            ValueError, match=f"pre_treatment_timepoints.*baseline_scope.*{scope}"
        ):
            _load(_cfg({"baseline_scope": scope, "pre_treatment_timepoints": ["t0"]}))[
                "de_novo_rule_params"
            ]()
        # The template carries the line empty: that is accepted everywhere.
        p = _load(_cfg({"baseline_scope": scope, "pre_treatment_timepoints": []}))[
            "de_novo_rule_params"
        ]()
        assert p["pre_treatment_timepoints"] == []

    def test_pre_treatment_timepoints_never_reach_the_rising_origin_flags(self):
        # They are a flag of the two commands that read samples;
        # alleleflux-rising-origin only counts the labels it is handed.
        args = _load(_cfg({"pre_treatment_timepoints": ["t0"]}))[
            "de_novo_rule_cli_args"
        ]()
        assert "pre_treatment" not in args and "t0" not in args

    def test_negative_numbers_raise(self):
        with pytest.raises(ValueError, match="cluster_window_bp"):
            _load(_cfg({"cluster_window_bp": -1}))["de_novo_rule_params"]()
        with pytest.raises(ValueError, match="min_determined_carriers"):
            _load(_cfg({"min_determined_carriers": 0}))["de_novo_rule_params"]()
