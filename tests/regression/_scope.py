"""The regression comparison contract: what is compared, what is skipped, what must be absent.

Both regression entry points apply this one contract so they can never
disagree:

* ``test_pipeline_e2e.py`` (via ``conftest.SCENARIOS``) -- a scenario's run
  against the committed golden snapshot;
* ``diff_runs.py`` -- any two output trees against each other.

This module imports nothing beyond ``__future__`` so every importer stays
clean (``conftest`` imports pytest; ``compare`` / ``diff_runs`` must not).
``compare.py`` imports the column rule from here rather than hard-coding
AlleleFlux column names.

The golden is captured from the branch under test, so there is no
cross-branch drift to tolerate: every tracked file is compared at the
comparator's strict default tolerance, minus the short exclusion list below.
"""

from __future__ import annotations

# Substrings identifying NON-DETERMINISTIC diagnostic columns whose *values*
# vary run-to-run and so must not gate equivalence.  The statsmodels LMM
# optimizer emits a run-dependent number of convergence warnings — and can flip
# the per-fit ``converged`` boolean — on near-degenerate positions (e.g. p≈1.0
# perfect ties), while the scientific outputs (p-value, coefficient, t-value) on
# those very rows are identical.  These columns are still required to be PRESENT
# (the comparator's schema gate catches a dropped/renamed column); only their
# *values* are skipped.  ``"_warnings"`` also matches ``"_n_warnings"`` as a
# substring, so one token covers both the warning text and the count.
IGNORE_VALUE_SUBSTRINGS = (
    "_warnings",  # {NUC}_warnings (text) and {NUC}_n_warnings (count)
    "_converged_LMM",  # per-fit convergence flag
)


# --- FULL contract: the whole output tree from THIS branch ------------------
# The golden is captured from the branch under test (``--update-golden``), so
# there is no cross-branch drift to carve out: every tracked file is compared
# at the comparator's strict default tolerance.  Only files that are inputs in
# disguise, large by construction, or pinned through what is derived from
# them are left out:
#   * ``_allele_frequency.parquet`` — the allele-frequency cache files
#     (``{mag}_{timepoint}_allele_frequency.parquet``).  Their content is
#     pinned by the site-filter sidecars beside them (compared) and by every
#     table read from them; the bytes themselves are not worth the golden's
#     size.  The pattern is the cache's file name, not ``.parquet``: the
#     allele-analysis step also writes parquet (``*_changes_mean.parquet``,
#     ``*_changes_no_zero-diff.parquet``) and those are derived outputs that
#     stay in.
#   * ``profiles/`` — the profiled TSVs are copies of the input profiles
#     (profiles mode) or the profiler's output (BAM mode); inputs either way.
# Logs, ``config_runtime_*`` and Snakemake sentinels are already untracked by
# the snapshot copier (``test_pipeline_e2e._copy_golden_payload``).
#   * ``inputMetadata/`` — the per-MAG sample sheets carry each profile's
#     ABSOLUTE path, so they differ between machines (CI vs a laptop) while
#     saying nothing the eligibility table does not.
#   * ``figures/`` — PNG images and an HTML report: not comparable byte-wise
#     and not what the regression pins.
FULL_EXCLUDE_SUBSTRINGS = (
    "_allele_frequency.parquet",
    "profiles/",
    "inputMetadata/",
    "figures/",
)

# Outputs a scenario must NOT produce, relative to the pipeline's data dir.
# The single-timepoint run has no earlier timepoint, so the four
# longitudinal-only trees must be absent (the workflow refuses the flags at
# parse time; this pins that nothing downstream writes them anyway).
EXPECTED_ABSENT = {
    "longitudinal": [],
    "single": [
        "dnds_analysis",
        "baseline_presence",
        "strain_turnover",
        "regional_contrast",
    ],
}


# --- The permutation-equivalence test's own scope ------------------------------
# ``test_permutation_equivalence.py`` compares a null run that REUSES the
# normal run's cache (and relabels in memory) against a null run REBUILT from
# the permuted sheet.  Two things differ between those by design, so the
# comparison is scoped to the statistics the relabelling feeds:
#   * the reuse path keeps the normal run's QC tables, whose per-group
#     summaries carry the unpermuted labels, while the rebuild recomputes them
#     from the permuted sheet;
#   * the rebuild is an ordinary run of a permuted sheet, not a flagged null
#     run, so it still writes the label-blind extras (dN/dS, ANI, turnover)
#     that a flagged null run switches off.
# The tolerance is the comparator's strict default: the test keeps LMM off
# (like the generated configs), and LMM was the one step whose numbers moved
# between runs (the statsmodels optimizer lands ~1e-7 apart; 115 LMM rows
# "mismatched" at rtol 1e-9 with identical printed values when it was on).
PERMUTATION_COMPARE_PREFIXES = (
    "significance_tests",
    "p_value_summary",
    "scores/intermediate",
    "scores/processed/combined",
    "preprocessing_eligibility",
    "eligibility_table",
)
