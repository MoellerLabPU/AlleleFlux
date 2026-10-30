"""
Unit tests for tests/regression/compare.py.

These keep the comparator itself trustworthy.  Without them, the regression
test would be a black box — if compare.py had a bug it would silently mask
or invent failures.

Each test follows the arrange / act / assert shape: build a golden and an
actual file on a temporary path, then either call ``compare_tables`` expecting
silence (a match) or wrap it in ``pytest.raises`` expecting a specific
complaint (a mismatch).  The ``match=`` regex pins down *which* failure fired,
so a test cannot pass by raising the wrong error for the wrong reason.

The most important pair is ``test_sign_flip_on_diff_column_is_tolerated`` and
``test_magnitude_change_on_diff_column_still_fails``: together they fence in
the one clever, slightly-dangerous feature in compare.py (magnitude-only
comparison of derivative columns) from both sides — it must forgive a sign
flip but never a magnitude change.
"""

from __future__ import annotations

import pandas as pd
import pytest

from tests.regression.compare import compare_tables, compare_trees


def _write_tsv(path, df):
    """Write ``df`` to ``path`` as a TSV, gzip-compressing iff the name ends in .gz.

    A tiny helper so each test reads as arrange/act/assert without repeating
    the ``to_csv`` boilerplate.  Compression is inferred from the suffix so the
    same helper serves ``.tsv`` and ``.tsv.gz`` fixtures.

    Parameters
    ----------
    path : Path
        Destination file.  A ``.gz`` suffix triggers gzip compression.
    df : pandas.DataFrame
        Data to write (no index column is written).
    """
    df.to_csv(
        path,
        sep="\t",
        index=False,
        compression="gzip" if str(path).endswith(".gz") else None,
    )


def test_identical_tables_match(tmp_path):
    """Baseline sanity: two byte-identical files must compare as a match.

    If this fails, the comparator is broken at the most basic level and no
    other test result can be trusted.
    """
    # Arrange: write the same frame to both actual and golden.
    df = pd.DataFrame(
        {"contig": ["c1", "c2"], "position": [10, 20], "A_change": [0.1, -0.2]}
    )
    a, g = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a, df)
    _write_tsv(g, df)
    # Act + Assert: a clean match returns None and must not raise.
    compare_tables(a, g)  # must not raise


def test_row_order_does_not_matter(tmp_path):
    """Row ordering is cosmetic: a reversed-row file must still match.

    Locks in the ``_sort_key_columns`` + stable-sort behaviour — parallel
    pipeline stages emit rows in arbitrary order, so the comparator canonicalises
    order before comparing.
    """
    # Arrange: actual is golden with its rows reversed.
    g = pd.DataFrame(
        {"contig": ["c1", "c2"], "position": [10, 20], "A_change": [0.1, -0.2]}
    )
    a = g.iloc[::-1].reset_index(drop=True)
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    # Act + Assert: sorting by (contig, position) realigns the rows -> match.
    compare_tables(a_path, g_path)


def test_column_order_does_not_matter(tmp_path):
    """Column ordering is cosmetic: a reordered-column file must still match.

    Locks in the set-level column check plus the reorder-to-golden step in
    ``compare_tables`` — column order varies run-to-run and is not semantic.
    """
    # Arrange: actual has the same columns as golden, shuffled.
    g = pd.DataFrame({"contig": ["c1"], "position": [10], "A_change": [0.1]})
    a = g[["A_change", "position", "contig"]]
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    # Act + Assert: column SETS match, so order is reconciled -> match.
    compare_tables(a_path, g_path)


def test_sign_flip_on_diff_column_fails(tmp_path):
    """A sign flip on a change / diff column is a real difference and must raise.

    History: the comparator once folded the sign away on ``_change`` /
    ``_diff`` / ``_diff_mean`` columns because two runs had been seen to
    "disagree in sign".  That was a misdiagnosis: the per-subject tables carry
    categorical keys whose category order differs between runs, and sorting by
    category code paired the control row with the treatment row of the same
    position -- which looks exactly like a flipped diff.  With keys sorted by
    value (``_sort_by_key_values``) two runs match sign and all, so the
    tolerance is gone and the sign is compared like any other digit.
    """
    g = pd.DataFrame(
        {
            "contig": ["c1", "c1"],
            "position": [10, 20],
            "A_frequency_diff_mean": [0.5, -0.3],
            "C_frequency_diff_mean": [-0.5, 0.3],
        }
    )
    a = g.copy()
    a["A_frequency_diff_mean"] *= -1
    a["C_frequency_diff_mean"] *= -1
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    with pytest.raises(AssertionError, match="mismatched"):
        compare_tables(a_path, g_path)


def test_magnitude_change_on_diff_column_still_fails(tmp_path):
    """A magnitude change on a change column raises (0.5 -> 0.6 is a regression)."""
    # Arrange: same column name pattern, but the magnitude actually changes.
    g = pd.DataFrame({"contig": ["c1"], "position": [10], "A_change": [0.5]})
    a = g.copy()
    a["A_change"] = [0.6]  # different magnitude — should fail
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    # Act + Assert: magnitudes differ beyond tolerance -> AssertionError.
    with pytest.raises(AssertionError, match="mismatched"):
        compare_tables(a_path, g_path)


def test_non_diff_column_sign_change_fails(tmp_path):
    """A sign flip on a plain frequency column raises (0.5 -> -0.5)."""
    # Arrange: negate a column that is NOT in the sign-insensitive set.
    g = pd.DataFrame({"contig": ["c1"], "position": [10], "A_frequency": [0.5]})
    a = g.copy()
    a["A_frequency"] = [-0.5]
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    # Act + Assert: signed comparison on a normal column -> AssertionError.
    with pytest.raises(AssertionError):
        compare_tables(a_path, g_path)


def test_nondeterministic_diag_column_value_is_ignored(tmp_path):
    """LMM convergence-diagnostic columns are compared on presence, not value.

    ``*_warnings`` / ``*_n_warnings`` / ``*_converged_LMM`` are non-reproducible
    run-to-run — the statsmodels optimizer's warning count and converged flag
    wobble on near-degenerate fits while the p-value on the same row is
    identical.  A pure diagnostic difference must therefore be forgiven.  This
    is the "forgive the wobble" half; the guardrails are the next two tests.
    """
    # Arrange: identical p-values, but the diagnostic columns differ.
    g = pd.DataFrame(
        {
            "contig": ["c1", "c1"],
            "position": [10, 20],
            "A_p_value_LMM": [0.374867, 1.0],
            "A_n_warnings": [4, 12],
            "A_converged_LMM": [True, False],
        }
    )
    a = g.copy()
    a["A_n_warnings"] = [8, 12]  # warning count wobbles on row 0
    a["A_converged_LMM"] = [False, False]  # converged flag flips on row 0
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    # Act + Assert: only diagnostic columns differ -> match.
    compare_tables(a_path, g_path)  # must not raise


def test_real_value_diff_beside_ignored_diag_still_fails(tmp_path):
    """Skipping diagnostic columns must not blind the comparator to the science.

    The guardrail: even when a non-deterministic ``A_n_warnings`` also differs,
    a genuine ``A_p_value_LMM`` change must still raise.  The value-skip is
    scoped to the diagnostic columns, not the row or the file.
    """
    # Arrange: a real p-value change AND a diagnostic wobble in the same row.
    g = pd.DataFrame(
        {
            "contig": ["c1"],
            "position": [10],
            "A_p_value_LMM": [0.374867],
            "A_n_warnings": [4],
        }
    )
    a = g.copy()
    a["A_p_value_LMM"] = [0.99]  # real regression
    a["A_n_warnings"] = [9]  # plus diagnostic noise
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    # Act + Assert: the p-value change is still caught despite the diag noise.
    with pytest.raises(AssertionError, match="mismatched"):
        compare_tables(a_path, g_path)


def test_dropping_ignored_diag_column_still_fails_schema(tmp_path):
    """Value-skip is not presence-skip: a dropped diagnostic column still fails.

    The diagnostic columns are exempt from VALUE comparison but NOT from the
    schema gate — removing ``A_converged_LMM`` entirely is a real schema
    regression and must be reported as a column-set mismatch.
    """
    # Arrange: golden has the diagnostic column; actual drops it.
    g = pd.DataFrame(
        {
            "contig": ["c1"],
            "position": [10],
            "A_p_value_LMM": [0.5],
            "A_converged_LMM": [True],
        }
    )
    a = g.drop(columns=["A_converged_LMM"])
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    # Act + Assert: schema gate fires before any value comparison.
    with pytest.raises(AssertionError, match="column-set mismatch"):
        compare_tables(a_path, g_path)


def test_missing_column_fails(tmp_path):
    """A differing column SET is a schema regression and must raise.

    Locks in the set-level column gate: a column present in golden but absent
    in actual (or vice versa) must be reported as a ``column-set mismatch``.
    """
    # Arrange: golden has ``extra_col``; actual does not.
    g = pd.DataFrame({"contig": ["c1"], "position": [10], "extra_col": [1.0]})
    a = pd.DataFrame({"contig": ["c1"], "position": [10]})
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    # Act + Assert: the column-set check fires before any value comparison.
    with pytest.raises(AssertionError, match="column-set mismatch"):
        compare_tables(a_path, g_path)


def test_row_count_mismatch_fails(tmp_path):
    """A differing row count must raise before any per-cell comparison.

    Locks in the shape gate: comparing files of different lengths is
    meaningless, so the comparator bails early with a ``row count`` message.
    """
    # Arrange: golden has two rows, actual has one.
    g = pd.DataFrame({"contig": ["c1", "c2"], "position": [10, 20]})
    a = pd.DataFrame({"contig": ["c1"], "position": [10]})
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    # Act + Assert: row-count gate fires -> AssertionError.
    with pytest.raises(AssertionError, match="row count"):
        compare_tables(a_path, g_path)


def test_tied_rows_distinguished_by_extra_keys_match(tmp_path):
    """Rows that tie on (mag_id, contig, position, gene_id) but differ only in
    test_type / group_analyzed must still canonicalise to the same order.

    Regression for the p_value_summary inflation: with only
    [mag_id, contig, position, gene_id] as sort keys, the four
    (test_type x group_analyzed) rows per position all tie, a stable sort
    cannot reorder them, and a shuffled actual mis-aligns row-for-row against
    golden -> thousands of *phantom* mismatches.  Widening the sort keys to
    include test_type/group_analyzed/source_file makes the row identity total,
    so a reordered-but-identical file matches.
    """
    # Arrange: four rows that share every current sort key, distinguished only
    # by test_type + group_analyzed (mirrors p_value_summary_single_sample).
    g = pd.DataFrame(
        {
            "mag_id": ["MAG_001"] * 4,
            "contig": ["c1"] * 4,
            "position": [14, 14, 14, 14],
            "gene_id": ["g1"] * 4,
            "test_type": ["tTest", "Wilcoxon", "tTest", "Wilcoxon"],
            "group_analyzed": ["treatment", "treatment", "control", "control"],
            "source_file": ["t.tsv.gz", "t.tsv.gz", "c.tsv.gz", "c.tsv.gz"],
            "min_p_value": [0.5, 1.0, 1.0, 1.0],
            "q_value": [0.72, 1.0, 1.0, 1.0],
        }
    )
    # actual = same rows in a different (parallel-emit) order.
    a = g.iloc[[2, 0, 3, 1]].reset_index(drop=True)
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    # Act + Assert: total identity key realigns the rows -> match.
    compare_tables(a_path, g_path)


def test_real_value_diff_among_tied_rows_still_fails(tmp_path):
    """Widening the sort keys must sharpen the signal, not swallow it.

    Same tied-row structure as the match test, but one row carries a genuine
    q_value change.  After canonical ordering the rows align by identity and
    the single real difference must still raise -- proving the fix does not
    over-tolerate.
    """
    # Arrange: identical to the match case, then change one real value.
    g = pd.DataFrame(
        {
            "mag_id": ["MAG_001"] * 4,
            "contig": ["c1"] * 4,
            "position": [14, 14, 14, 14],
            "gene_id": ["g1"] * 4,
            "test_type": ["tTest", "Wilcoxon", "tTest", "Wilcoxon"],
            "group_analyzed": ["treatment", "treatment", "control", "control"],
            "source_file": ["t.tsv.gz", "t.tsv.gz", "c.tsv.gz", "c.tsv.gz"],
            "min_p_value": [0.5, 1.0, 1.0, 1.0],
            "q_value": [0.72, 1.0, 1.0, 1.0],
        }
    )
    a = g.copy()
    a.loc[0, "q_value"] = 0.99  # genuine change on the (tTest, treatment) row
    a = a.iloc[[2, 0, 3, 1]].reset_index(drop=True)  # shuffle to prove alignment
    a_path, g_path = tmp_path / "a.tsv", tmp_path / "g.tsv"
    _write_tsv(a_path, a)
    _write_tsv(g_path, g)
    # Act + Assert: the one real diff survives canonical alignment -> raise.
    with pytest.raises(AssertionError, match="mismatched"):
        compare_tables(a_path, g_path)


def test_include_prefixes_scopes_the_sweep(tmp_path):
    """compare_trees(include_prefixes=...) diffs only the allowlisted subtrees.

    A golden file outside the allowlist (here ``profiles/``) is neither
    compared nor flagged MISSING, so intended structural divergence in excluded
    directories cannot fail the run.
    """
    # Arrange: a shared in-scope file, plus a golden-only OUT-of-scope file.
    golden, actual = tmp_path / "golden", tmp_path / "actual"
    (golden / "significance_tests").mkdir(parents=True)
    (actual / "significance_tests").mkdir(parents=True)
    df = pd.DataFrame({"contig": ["c1"], "position": [10], "p": [0.5]})
    _write_tsv(golden / "significance_tests" / "s.tsv", df)
    _write_tsv(actual / "significance_tests" / "s.tsv", df)
    (golden / "profiles").mkdir()
    _write_tsv(golden / "profiles" / "p.tsv", df)  # out of scope -> ignored
    # Act: restrict the sweep to significance_tests.
    fails = compare_trees(actual, golden, include_prefixes=["significance_tests"])
    # Assert: the out-of-scope golden-only file is not reported.
    assert fails == []


def test_allow_golden_only_tolerates_intended_subset(tmp_path):
    """allow_golden_only lets a golden file the branch never emits pass.

    Models the config-driven taxa-level subset: main computes 7 taxonomic
    levels, the branch computes 3, so golden carries levels the branch does not
    produce.  Inside an allow_golden_only prefix those absences are expected
    (not MISSING) -- while a *shared* file is still compared strictly, and
    without the tolerance the same absence IS flagged.
    """
    # Arrange: one shared level (genus) and one golden-only level (class).
    golden, actual = tmp_path / "golden", tmp_path / "actual"
    base = "scores/processed/combined"
    (golden / base / "genus").mkdir(parents=True)
    (actual / base / "genus").mkdir(parents=True)
    (golden / base / "class").mkdir(parents=True)
    df = pd.DataFrame({"contig": ["c1"], "position": [10], "p": [0.5]})
    _write_tsv(golden / base / "genus" / "x.tsv", df)
    _write_tsv(actual / base / "genus" / "x.tsv", df)
    _write_tsv(golden / base / "class" / "y.tsv", df)  # no actual twin
    # Act + Assert: tolerated when the prefix is allow_golden_only.
    fails = compare_trees(
        actual, golden, include_prefixes=[base], allow_golden_only=[base]
    )
    assert fails == []
    # And without the tolerance, the absent level IS reported MISSING.
    fails2 = compare_trees(actual, golden, include_prefixes=[base])
    assert any("MISSING" in f and "class" in f for f in fails2)


def test_compare_trees_threads_tolerance(tmp_path):
    """compare_trees forwards rtol/atol down to the table comparator.

    The golden is captured from float64 ``main`` while the branch under test
    stores frequencies as float32 (~1e-7 precision), so the strict default
    tolerance flags benign rounding.  A float drift below a loosened ``atol``
    must pass when the tolerance is threaded through, and still fail at the
    strict default -- proving the knob actually reaches ``np.isclose`` and is
    not silently dropped at a wrapper.
    """
    # Arrange: a ~4e-4 drift on a signed (not sign-folded) frequency column.
    golden, actual = tmp_path / "golden", tmp_path / "actual"
    golden.mkdir()
    actual.mkdir()
    g = pd.DataFrame(
        {"contig": ["c1", "c2"], "position": [10, 20], "A_frequency": [0.5, 0.25]}
    )
    a = g.copy()
    a["A_frequency"] = [0.5004, 0.2503]
    _write_tsv(golden / "f.tsv", g)
    _write_tsv(actual / "f.tsv", a)
    # Act + Assert: strict default flags it; loosened tolerance clears it.
    assert compare_trees(actual, golden)  # non-empty -> drift caught at 1e-9
    assert compare_trees(actual, golden, rtol=1e-5, atol=1e-3) == []


def test_exclude_substrings_skips_matching_files(tmp_path):
    """exclude_substrings drops in-scope files whose path contains a substring.

    Models the single scenario: ``two_sample_paired`` is non-deterministic on
    single-timepoint data and the ``preprocessed`` intermediates are
    column-pruned, so both are excluded by substring even though they fall
    under an included prefix (``significance_tests``).  Substring (not prefix)
    matching is required because ``two_sample_paired`` also appears mid-path in
    the ``p_value_summary`` filenames.  A kept in-scope file still compares.
    """
    # Arrange: a kept file (single_sample) and an excluded one (two_sample_paired
    # with no actual twin, which would otherwise be flagged MISSING).
    golden, actual = tmp_path / "golden", tmp_path / "actual"
    df = pd.DataFrame({"contig": ["c1"], "position": [10], "p": [0.5]})
    (golden / "significance_tests" / "single_sample_x").mkdir(parents=True)
    (actual / "significance_tests" / "single_sample_x").mkdir(parents=True)
    _write_tsv(golden / "significance_tests" / "single_sample_x" / "a.tsv", df)
    _write_tsv(actual / "significance_tests" / "single_sample_x" / "a.tsv", df)
    (golden / "significance_tests" / "two_sample_paired_x").mkdir(parents=True)
    _write_tsv(golden / "significance_tests" / "two_sample_paired_x" / "b.tsv", df)
    # Act + Assert: excluded substring -> the paired file is skipped, not MISSING.
    fails = compare_trees(
        actual,
        golden,
        include_prefixes=["significance_tests"],
        exclude_substrings=["two_sample_paired"],
    )
    assert fails == []
    # Without the exclude, the same paired file IS reported MISSING.
    fails2 = compare_trees(actual, golden, include_prefixes=["significance_tests"])
    assert any("MISSING" in f and "two_sample_paired" in f for f in fails2)


def test_compare_trees_reports_missing_file(tmp_path):
    """A golden file with no actual twin is reported as MISSING, not skipped.

    Locks in the "golden defines the contract" rule in ``compare_trees``: every
    golden file must have a counterpart under the actual tree, and the absence
    is collected as a failure string rather than silently ignored.
    """
    # Arrange: a golden file that has no counterpart in the actual tree.
    golden = tmp_path / "golden"
    actual = tmp_path / "actual"
    golden.mkdir()
    actual.mkdir()
    df = pd.DataFrame({"contig": ["c1"], "position": [10]})
    _write_tsv(golden / "a.tsv", df)
    # Don't write to actual — file is missing
    # Act: sweep the tree (compare_trees collects failures, does not raise).
    fails = compare_trees(actual, golden)
    # Assert: exactly one failure, flagged as MISSING.
    assert len(fails) == 1
    assert "MISSING" in fails[0]


def test_present_forbidden_paths_are_reported(tmp_path):
    """``forbidden_present`` names every expected-absent path that exists."""
    from tests.regression.compare import forbidden_present

    (tmp_path / "dnds_analysis").mkdir()
    (tmp_path / "dnds_analysis" / "x.tsv").write_text("a\n1\n")
    (tmp_path / "lonely.tsv").write_text("a\n1\n")
    out = forbidden_present(
        tmp_path, ["dnds_analysis", "baseline_presence", "lonely.tsv"]
    )
    assert out == [
        "FORBIDDEN (exists): dnds_analysis",
        "FORBIDDEN (exists): lonely.tsv",
    ]


def test_absent_forbidden_paths_pass(tmp_path):
    from tests.regression.compare import forbidden_present

    assert forbidden_present(tmp_path, ["dnds_analysis"]) == []


def test_full_contract_excludes_the_cache_but_not_the_allele_analysis_parquet():
    """Only the allele-frequency cache files are left out of the full-tree
    contract; the allele-analysis step's parquet tables are derived outputs
    and must stay pinned."""
    from tests.regression._scope import FULL_EXCLUDE_SUBSTRINGS

    def excluded(rel):
        return any(s in rel for s in FULL_EXCLUDE_SUBSTRINGS)

    assert excluded("allele_freq_cache/MAG_001_post_allele_frequency.parquet")
    assert not excluded(
        "allele_freq_cache/MAG_001_post_allele_frequency_site_depth_filter.tsv"
    )
    assert not excluded(
        "allele_analysis/allele_analysis_pre_post-treatment_control/"
        "MAG_001_allele_frequency_changes_mean.parquet"
    )
    assert not excluded(
        "allele_analysis/allele_analysis_pre_post-treatment_control/"
        "MAG_001_allele_frequency_changes_no_zero-diff.parquet"
    )


def test_categorical_keys_with_different_category_order_still_align(tmp_path):
    """Parquet keeps pandas categoricals, and a pipeline run records the
    categories in first-seen order, so two runs of the same table can carry
    the same rows under different category orders.  Sorting by category CODE
    would then pair row i of one file with a different subject in the other;
    the comparator must sort by the key's value instead."""
    from tests.regression.compare import compare_tables

    rows = pd.DataFrame(
        {
            "contig": ["c1"] * 4,
            "position": [0, 0, 1, 1],
            "subjectID": ["subj3", "subj8", "subj3", "subj8"],
            "A_frequency": [0.9, 0.1, 0.8, 0.2],
        }
    )
    golden = rows.copy()
    golden["subjectID"] = pd.Categorical(
        golden["subjectID"], categories=["subj3", "subj8"]
    )
    actual = rows.iloc[[1, 0, 3, 2]].copy()  # same rows, other order ...
    actual["subjectID"] = pd.Categorical(
        actual["subjectID"], categories=["subj8", "subj3"]
    )  # ... other category order
    g, a = tmp_path / "g.parquet", tmp_path / "a.parquet"
    golden.to_parquet(g, index=False)
    actual.to_parquet(a, index=False)
    compare_tables(a, g)  # must not raise
