"""Shared comparison scope + float tolerance for branch-vs-main regression.

Both regression tools diff a candidate branch against a ``main``-derived
reference, and both must apply the **same** contract or they disagree:

* ``test_pipeline_e2e.py`` (via ``conftest.SCENARIOS``) — snapshot-based.
* ``_diff_outputs.py`` (via ``diff_branches.sh``) — live cross-branch.

This module is that single source of truth.  It lives here — not in
``compare.py`` (which stays a generic, AlleleFlux-agnostic oracle) and not in
``conftest.py`` (which imports pytest, and ``_diff_outputs`` must not) — so both
importers stay clean.

The goal: verify the branch produces the same *science* as main while ignoring
intentional structural drift.  Rationale per prefix below.
"""

from __future__ import annotations

# --- What to compare ------------------------------------------------------
# Structurally-stable statistical outputs — identical filenames/format on both
# branches, and the integration point where any upstream numerical change
# surfaces:
#   * significance_tests/, p_value_summary/  — raw + summarised stats.
#   * scores/intermediate/                   — per-MAG scores; identical 1:1.
#   * scores/processed/combined/             — taxonomic aggregation; main
#     computes all 7 levels, the branch only the ``taxa_score_levels`` subset,
#     so the extra golden levels are allowed-missing (shared levels strict).
#   * preprocessing_eligibility/, eligibility_table — QC/filter decisions; 1:1.
STAT_COMPARE_ONLY = [
    "significance_tests",
    "p_value_summary",
    "scores/intermediate",
    "scores/processed/combined",
    "preprocessing_eligibility",
    "eligibility_table",
]

# Golden (main) files the branch intentionally does not emit — not MISSING:
#   * scores/processed/combined/<level>/ for taxa levels omitted by config.
# Excluded entirely (not even in COMPARE_ONLY): profiles/ (mapq_scores column
# dropped), allele_analysis/ (tsv.gz -> parquet, files dropped), QC/ +
# inputMetadata/ (renamed group-independent), outlier_genes/
# (use_outlier_detection=False), scores/processed/gene_scores_*
# (use_gene_scores=False).
STAT_ALLOW_GOLDEN_ONLY = ["scores/processed/combined"]

# --- How loosely to compare floats ----------------------------------------
# The golden is float64 (main), but the branch stores allele frequencies as
# float32 (~1e-7 precision — a deliberate memory optimization in
# _allele_freq_common.py), so the strict comparator default (rtol=1e-9) is
# unachievable by construction.  These bounds sit ~4 orders above the observed
# benign drift (max ~1.3e-4 on a q_value) yet still catch any scientifically
# meaningful change (a p-value moving 0.01, a frequency moving 0.001).  Tighten
# back toward the defaults if the golden and the branch are ever brought to the
# same float precision.
STAT_RTOL = 1e-5
STAT_ATOL = 1e-3

# --- Per-scenario carve-outs ----------------------------------------------
# Substrings matched anywhere in a golden file's relative path, applied AFTER
# STAT_COMPARE_ONLY.  Used to drop non-comparable subtrees that fall *inside* a
# broad include prefix (e.g. the preprocessed intermediates live under
# ``significance_tests/preprocessed_*/``).  Keyed by scenario name.
#
#   single:
#     * two_sample_paired — DETERMINISTIC but float32-sensitive (verified: two
#       same-branch runs are bit-identical; the paired-test code is identical
#       across branches — empty git diff — so only the input differs).  The
#       single paired test runs Wilcoxon on RAW frequencies (longitudinal uses
#       diff_means).  The branch stores frequencies as float32, so its values
#       differ from the float64 golden by ~1e-7, which flips a near-tie RANK in
#       the discrete small-n Wilcoxon → the p-value jumps a discrete step
#       (e.g. 0.375↔0.625) that the float tolerance can't absorb, cascading into
#       q-values.  Scientifically null (a non-significant p stays non-significant).
#       longitudinal's paired test uses diff_means (continuous, no rank flips),
#       so it is unaffected and stays in scope.
#     * preprocessed — the *_allele_frequency_preprocessed intermediate inputs
#       carry metadata columns (mapq_scores, coverage, breadth, …) the branch's
#       refactor prunes, so they column-set-mismatch.  Intermediate, not output.
#   longitudinal: none — its paired test (pre-vs-post) is deterministic and its
#     preprocessed files match, so both stay in scope for full coverage.
SCENARIO_EXCLUDE_SUBSTRINGS = {
    "longitudinal": [],
    "single": ["two_sample_paired", "preprocessed"],
}


def exclude_substrings_for(scenario: str) -> list[str]:
    """Return the substring carve-outs for ``scenario`` (empty if it has none)."""
    return SCENARIO_EXCLUDE_SUBSTRINGS.get(scenario, [])
