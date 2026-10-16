"""Unit tests for the site-depth-filter config helpers in ``common.smk``.

The three rules that read profiles (allele_freq_cache, pairwise_ani,
baseline_presence) must pass the SAME filter flags, built from one config
block, so the cache, the ANI tables and the baseline-presence tables can never
be filtered differently.  The two helpers are sliced out of the ``.smk`` source
by indentation and exec'd with a stub ``config``, as in
``test_rising_origin_targets.py``.
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
FUNCS = ("site_depth_filter_params", "site_depth_filter_cli_args")


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


class TestSiteDepthFilterParams:
    def test_params_default_to_the_house_values(self):
        fns = _load({"quality_control": {}})
        assert fns["site_depth_filter_params"]() == {
            "enabled": True,
            "site_depth": 0,
            "site_ratio": 3.0,
            "min_allele_reads": 0,
            "keep_original_counts": True,
        }
        assert (
            fns["site_depth_filter_cli_args"]()
            == "--site_depth 0 --site_ratio 3.0 --min_allele_reads 0"
        )

    def test_missing_block_entirely(self):
        fns = _load({})
        assert fns["site_depth_filter_params"]()["enabled"] is True

    def test_disabled_and_overridden_to_midas(self):
        fns = _load(
            {
                "quality_control": {
                    "site_depth_filter": {
                        "enabled": False,
                        "site_depth": 5,
                        "site_ratio": 2,
                        "min_allele_reads": 2,
                        "keep_original_counts": False,
                    }
                }
            }
        )
        assert (
            fns["site_depth_filter_cli_args"]()
            == "--disable_site_depth_filter --site_depth 5 --site_ratio 2.0 "
            "--min_allele_reads 2 --zero_low_counts"
        )
