"""Unit tests for the parse-time data-type guard in ``common.smk``.

``check_data_type_features`` runs while the Snakefile is still being read and
refuses every setting that cannot mean anything for the configured data type:
longitudinal-only features switched on for single-timepoint data, timepoint
combinations of the wrong length, an unknown ``data_type``, and strain
turnover without the ANI table it reads.  ``pairwise_ani_settings`` resolves
the two pairwise-ANI mode knobs with data-type-aware defaults (single data
has one sample per mouse, so same-mouse modes are refused there).

Both functions are sliced out of the ``.smk`` source and exec'd on their own,
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
FUNCS = ("check_data_type_features", "pairwise_ani_settings")


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
    src = "\n\n".join(_slice_function(lines, f) for f in FUNCS)
    namespace = {}
    exec(compile(src, str(COMMON_SMK), "exec"), namespace)
    return {f: namespace[f] for f in FUNCS}


def _cfg(data_type, combos=None, **flags):
    """A minimal config: one well-formed timepoint combination unless given."""
    if combos is None:
        combos = (
            [{"timepoint": ["t1"], "focus": "t1"}]
            if data_type == "single"
            else [{"timepoint": ["pre", "post"], "focus": "post"}]
        )
    analysis = {"data_type": data_type, "timepoints_combinations": combos}
    analysis.update(flags)
    return {"analysis": analysis}


class TestCheckDataTypeFeatures:
    def test_clean_single_config_passes_and_returns_the_data_type(self):
        fns = _load()
        assert fns["check_data_type_features"](_cfg("single")) == "single"

    def test_clean_longitudinal_config_passes(self):
        fns = _load()
        cfg = _cfg("longitudinal", use_dnds=True, use_baseline_presence=True)
        assert fns["check_data_type_features"](cfg) == "longitudinal"

    def test_unknown_data_type_is_refused_by_name(self):
        fns = _load()
        with pytest.raises(ValueError, match="cross_sectional"):
            fns["check_data_type_features"](_cfg("cross_sectional"))

    def test_single_with_dnds_on_is_refused(self):
        fns = _load()
        with pytest.raises(ValueError, match="use_dnds"):
            fns["check_data_type_features"](_cfg("single", use_dnds=True))

    def test_single_with_absent_longitudinal_flags_passes(self):
        # Absent means "default", and the default for single data is off; only
        # an explicit True is a user mistake worth stopping the run for.
        fns = _load()
        assert fns["check_data_type_features"](_cfg("single")) == "single"

    def test_single_with_flags_explicitly_false_passes(self):
        fns = _load()
        cfg = _cfg("single", use_dnds=False, use_baseline_presence=False)
        assert fns["check_data_type_features"](cfg) == "single"

    def test_single_names_every_offending_flag_at_once(self):
        fns = _load()
        cfg = _cfg(
            "single",
            use_regional_contrast=True,
            use_baseline_presence=True,
            use_strain_turnover=True,
            use_pairwise_ani=True,
        )
        with pytest.raises(ValueError) as exc:
            fns["check_data_type_features"](cfg)
        msg = str(exc.value)
        assert "use_regional_contrast" in msg
        assert "use_baseline_presence" in msg
        assert "use_strain_turnover" in msg
        assert "use_dnds" not in msg

    def test_single_with_a_two_timepoint_combination_is_refused(self):
        fns = _load()
        cfg = _cfg("single", combos=[{"timepoint": ["pre", "post"], "focus": "post"}])
        with pytest.raises(ValueError, match="exactly 1"):
            fns["check_data_type_features"](cfg)

    def test_longitudinal_with_a_one_timepoint_combination_is_refused(self):
        fns = _load()
        cfg = _cfg(
            "longitudinal",
            combos=[
                {"timepoint": ["pre", "post"], "focus": "post"},
                {"timepoint": ["mid"], "focus": "mid"},
            ],
        )
        with pytest.raises(ValueError, match="exactly 2") as exc:
            fns["check_data_type_features"](cfg)
        assert "mid" in str(exc.value)

    def test_turnover_without_pairwise_ani_is_refused(self):
        fns = _load()
        cfg = _cfg("longitudinal", use_strain_turnover=True)
        with pytest.raises(ValueError, match="use_pairwise_ani"):
            fns["check_data_type_features"](cfg)

    def test_turnover_with_pairwise_ani_passes(self):
        fns = _load()
        cfg = _cfg("longitudinal", use_strain_turnover=True, use_pairwise_ani=True)
        assert fns["check_data_type_features"](cfg) == "longitudinal"


class TestPairwiseAniSettings:
    def test_single_defaults_to_all_pairs_and_no_snp_locations(self):
        fns = _load()
        assert fns["pairwise_ani_settings"](_cfg("single"), "single") == {
            "pairs": "all",
            "store_snp_locations": "none",
        }

    def test_longitudinal_defaults_to_same_mouse_pairs(self):
        fns = _load()
        out = fns["pairwise_ani_settings"](_cfg("longitudinal"), "longitudinal")
        assert out == {
            "pairs": "within_subject",
            "store_snp_locations": "within_subject",
        }

    def test_explicit_values_are_kept_on_longitudinal(self):
        fns = _load()
        cfg = _cfg(
            "longitudinal",
            pairwise_ani={"pairs": "transitions", "store_snp_locations": "none"},
        )
        out = fns["pairwise_ani_settings"](cfg, "longitudinal")
        assert out == {"pairs": "transitions", "store_snp_locations": "none"}

    @pytest.mark.parametrize("mode", ["within_subject", "transitions"])
    def test_single_refuses_same_mouse_pair_modes(self, mode):
        fns = _load()
        cfg = _cfg("single", pairwise_ani={"pairs": mode})
        with pytest.raises(ValueError, match=mode):
            fns["pairwise_ani_settings"](cfg, "single")

    def test_single_refuses_same_mouse_snp_locations(self):
        fns = _load()
        cfg = _cfg("single", pairwise_ani={"store_snp_locations": "within_subject"})
        with pytest.raises(ValueError, match="store_snp_locations"):
            fns["pairwise_ani_settings"](cfg, "single")

    def test_single_explicit_all_and_all_passes(self):
        fns = _load()
        cfg = _cfg(
            "single", pairwise_ani={"pairs": "all", "store_snp_locations": "all"}
        )
        out = fns["pairwise_ani_settings"](cfg, "single")
        assert out == {"pairs": "all", "store_snp_locations": "all"}

    def test_other_knobs_in_the_block_are_ignored(self):
        fns = _load()
        cfg = _cfg("single", pairwise_ani={"min_cov": 5, "min_freq": 0.05})
        out = fns["pairwise_ani_settings"](cfg, "single")
        assert out == {"pairs": "all", "store_snp_locations": "none"}
