"""Which cache files a run produces, and which of them each consumer waits for.

``_cache_universe`` in ``dynamic_targets.smk`` walks every (comparison, group
pair, eligible genome, timepoint) a cache file exists for.  Two kinds of
reader take slices of it:

* the statistics-side readers (pairwise ANI, dN/dS) want a genome only at the
  timepoints of the comparisons it is tested in;
* the de novo readers (baseline presence, gene de novo) under a dated cohort
  scope also read **witnesses**: samples of the genome at every witness
  timepoint (``CACHE_TIMEPOINTS``), whether or not the genome is tested there.

So the walk takes one argument, ``witness_timepoints``.  With it on, and the
cohort condition met (``COHORT_WITNESSES_ACTIVE``), every eligible genome is
listed at every witness timepoint; otherwise at its comparison's own.

Small example used below: comparisons ``a_b`` and ``a_c``; M1 is eligible in
both, M2 in ``a_b`` only; the sheet also holds timepoint ``x`` (in no
comparison).  M2 has samples at ``c`` and ``x`` that can testify for an allele
rising in ``a_b``, so under the cohort condition M2 needs caches at c and x.

Sliced out of the ``.smk`` source and exec'd with stubs, as in
``test_site_depth_filter_summary_targets.py``.
"""

from __future__ import annotations

import os
from pathlib import Path

PIPELINE = (
    Path(__file__).resolve().parents[1]
    / "alleleflux"
    / "smk_workflow"
    / "alleleflux_pipeline"
)
DYNAMIC_TARGETS_SMK = PIPELINE / "shared" / "dynamic_targets.smk"
RULES = PIPELINE / "rules"
FUNCS = (
    "_cache_universe",
    "generate_allele_freq_cache_sidecar_targets",
    "allele_freq_cache_files_for_mag",
    "allele_freq_cache_files_for_comparison",
)
ELIGIBLE = {("a_b", "g1_g2"): ["M1", "M2"], ("a_c", "g1_g2"): ["M1"]}


def _slice_function(lines, name):
    start = next(i for i, ln in enumerate(lines) if ln.startswith(f"def {name}("))
    end = start + 1
    while end < len(lines):
        ln = lines[end]
        if ln.strip() and not ln[0].isspace():
            break
        end += 1
    return "\n".join(lines[start:end])


def _cache(mag, timepoint):
    return os.path.join(
        "RUN",
        "allele_freq_cache",
        timepoint,
        f"{mag}_{timepoint}_allele_frequency.parquet",
    )


def _load(active, cache_tps, eligible=ELIGIBLE, data_type="longitudinal"):
    """Exec the helpers with a stub run: ``active`` is COHORT_WITNESSES_ACTIVE,
    ``cache_tps`` is CACHE_TIMEPOINTS, ``eligible`` maps (comparison, groups)
    to its eligible genomes."""
    lines = DYNAMIC_TARGETS_SMK.read_text().splitlines()
    src = "\n\n".join(_slice_function(lines, f) for f in FUNCS)
    labels = list(dict.fromkeys(tp for tp, _ in eligible))
    namespace = {
        "os": os,
        "DATA_TYPE": data_type,
        "timepoints_labels": labels,
        "groups_labels": list(dict.fromkeys(gr for _, gr in eligible)),
        "CACHE_TIMEPOINTS": list(cache_tps),
        "COHORT_WITNESSES_ACTIVE": active,
        "_get_mags_by_eligibility": lambda tp, gr, eligibility_type: eligible.get(
            (tp, gr), []
        ),
        "get_allele_freq_cache_path": lambda mag_wildcard, timepoint_wildcard: _cache(
            mag_wildcard, timepoint_wildcard
        ),
    }
    exec(compile(src, str(DYNAMIC_TARGETS_SMK), "exec"), namespace)
    return {f: namespace[f] for f in FUNCS}


class TestCacheUniverse:
    def test_cohort_universe_covers_every_cache_timepoint_for_every_eligible_mag(
        self,
    ):
        fns = _load(active=True, cache_tps=["a", "b", "c", "x"])
        rows = set(fns["_cache_universe"](True))
        assert rows == {
            ("a_b", "g1_g2", mag, tp)
            for mag in ("M1", "M2")
            for tp in ("a", "b", "c", "x")
        } | {("a_c", "g1_g2", "M1", tp) for tp in ("a", "b", "c", "x")}

    def test_default_universe_is_the_comparison_timepoints_only(self):
        # Condition not met (another scope, undated sheet, no de novo
        # analysis): asking for witness timepoints changes nothing.
        fns = _load(active=False, cache_tps=["a", "b", "c"])
        expect = {
            ("a_b", "g1_g2", "M1", "a"),
            ("a_b", "g1_g2", "M1", "b"),
            ("a_b", "g1_g2", "M2", "a"),
            ("a_b", "g1_g2", "M2", "b"),
            ("a_c", "g1_g2", "M1", "a"),
            ("a_c", "g1_g2", "M1", "c"),
        }
        assert set(fns["_cache_universe"](True)) == expect
        assert set(fns["_cache_universe"](False)) == expect

    def test_without_witness_timepoints_the_cohort_condition_changes_nothing(self):
        fns = _load(active=True, cache_tps=["a", "b", "c", "x"])
        assert ("a_b", "g1_g2", "M2", "c") not in set(fns["_cache_universe"](False))
        assert ("a_b", "g1_g2", "M2", "x") not in set(fns["_cache_universe"](False))

    def test_comparisons_setting_still_caches_a_tested_genome_at_every_compared_timepoint(
        self,
    ):
        # cohort_witnesses: comparisons -> no extra timepoint (no x), but M2
        # (tested in a_b only) still gets a cache at c: its samples there are
        # at a compared timepoint, so they may testify.
        fns = _load(active=True, cache_tps=["a", "b", "c"])
        rows = set(fns["_cache_universe"](True))
        assert ("a_b", "g1_g2", "M2", "c") in rows
        assert not any(tp == "x" for *_, tp in rows)


class TestSlices:
    def test_every_cache_file_of_the_run_has_its_sidecar_target(self):
        fns = _load(active=True, cache_tps=["a", "b", "c", "x"])
        got = fns["generate_allele_freq_cache_sidecar_targets"]()
        assert got == sorted(
            _cache(mag, tp)[: -len(".parquet")] + "_site_depth_filter.tsv"
            for mag in ("M1", "M2")
            for tp in ("a", "b", "c", "x")
        )

    def test_de_novo_readers_wait_for_the_witness_caches_of_their_comparison(self):
        fns = _load(active=True, cache_tps=["a", "b", "c", "x"])
        got = fns["allele_freq_cache_files_for_comparison"](
            "a_c", "g1_g2", witness_timepoints=True
        )
        assert got == sorted(_cache("M1", tp) for tp in ("a", "b", "c", "x"))

    def test_dnds_waits_for_its_own_two_timepoints_only(self):
        fns = _load(active=True, cache_tps=["a", "b", "c", "x"])
        got = fns["allele_freq_cache_files_for_comparison"](
            "a_c", "g1_g2", witness_timepoints=False
        )
        assert got == sorted(_cache("M1", tp) for tp in ("a", "c"))

    def test_ani_keeps_to_the_timepoints_a_genome_is_tested_at(self):
        # Pairwise ANI must not change with the de novo settings: M2 is tested
        # in a_b only, so ANI reads its a and b samples, whether or not a
        # cache of M2 at c or x exists for the de novo analysis.
        for active in (True, False):
            fns = _load(active=active, cache_tps=["a", "b", "c", "x"])
            assert fns["allele_freq_cache_files_for_mag"]("M2") == sorted(
                [_cache("M2", "a"), _cache("M2", "b")]
            )
            assert fns["allele_freq_cache_files_for_mag"]("M1") == sorted(
                _cache("M1", tp) for tp in ("a", "b", "c")
            )
            assert fns["allele_freq_cache_files_for_mag"]("M9") == []


class TestRulesAskForTheRightSlice:
    """The rule files say which slice each reader takes (text checks: a rule
    cannot be exec'd alone)."""

    def test_de_novo_rules_ask_for_witness_timepoints(self):
        for name in ("baseline_presence.smk", "gene_de_novo.smk"):
            text = (RULES / name).read_text()
            assert "witness_timepoints=True" in text, name

    def test_dnds_rule_does_not(self):
        text = (RULES / "dnds_analysis.smk").read_text()
        assert "witness_timepoints=False" in text

    def test_ani_rule_reads_only_the_tested_timepoints_and_the_config_groups(self):
        text = (RULES / "pairwise_ani.smk").read_text()
        # The QC tables follow the cache files the job waits for, so the
        # sample list can never name a timepoint whose cache is not an input.
        assert "qc_tables_for_cache_files(wildcards.mag, input.caches)" in text
        assert "--groups {params.groups}" in text

    def test_metadata_rule_passes_the_witness_groups(self):
        text = (RULES / "metadata.smk").read_text()
        assert '" ".join(WITNESS_GROUPS)' in text
        assert "_all_group_values" not in text

    def test_de_novo_rules_hand_both_lists_to_their_command(self):
        for name in ("baseline_presence.smk", "gene_de_novo.smk"):
            text = (RULES / name).read_text()
            assert "--witness_timepoints {params.witness_timepoints}" in text, name
            assert "--witness_groups {params.witness_groups}" in text, name
