"""Unit tests for ``alleleflux-eligibility`` reading one QC folder per timepoint.

QC writes one table per genome per timepoint, holding per-sample facts only
(was this sample readable for this genome?).  The eligibility step reads the
tables of the comparison's timepoints and works out, for the two tested
groups:

* which subjects are readable at BOTH timepoints (longitudinal) -- only their
  samples count;
* how many replicates (cages, blocks) each tested group has among them;
* how many replicates the two tested groups SHARE (the paired count).

The shared count is taken over the two tested groups only.  With three groups
that all reuse the same replicate IDs, counting over every group of the run
used to call no replicate paired; that case and its neighbours are pinned
below, and a seeded brute-force check compares the script with a plain set
definition on a few hundred random designs.

Test data are written as real ``{mag}_QC.tsv`` files in temporary folders, so
the text-reading rules (numeric-looking group names, header-only tables) are
exercised exactly as in a run.
"""

from __future__ import annotations

import random
import subprocess
import unittest

import numpy as np
import pandas as pd
import pytest

from alleleflux.scripts.preprocessing.eligibility_table import (
    add_subject_count_per_group,
    analyze_qc_dirs,
    check_timepoints,
    count_paired_replicates,
    load_genome_qc,
    pair_and_count,
)

# Columns of the QC tables written here.  A real table carries a dozen more
# per-sample statistics; ``file_path`` and ``breadth`` stand in for them to
# show that columns the eligibility step does not need are ignored.
QC_COLUMNS = [
    "sample_id",
    "MAG_ID",
    "file_path",
    "group",
    "subjectID",
    "replicate",
    "breadth",
    "coverage_threshold_passed",
]


def _subject(subject, group, replicate, pre=True, post=True):
    """The samples of one subject as ``(subject, group, replicate, time, readable)``.

    ``pre`` / ``post``: True = a sample that passed QC, False = a sample that
    failed QC, None = no sample at that timepoint.

    Example: ``_subject("m1", "A", "c1", post=False)`` -> a readable ``pre``
    sample and a ``post`` sample that failed QC.
    """
    return [
        (subject, group, replicate, time, state)
        for time, state in (("pre", pre), ("post", post))
        if state is not None
    ]


def _write_qc(folder, mag, samples, with_time=True):
    """Write ``{mag}_QC.tsv`` into ``folder``; no sample -> a header-only table."""
    folder.mkdir(parents=True, exist_ok=True)
    columns = QC_COLUMNS + (["time"] if with_time else [])
    rows = []
    for subject, group, replicate, time, readable in samples:
        sample_id = f"{subject}_{time}" if with_time else subject
        row = {
            "sample_id": sample_id,
            "MAG_ID": mag,
            "file_path": f"/profiles/{sample_id}.tsv.gz",
            "group": group,
            "subjectID": subject,
            "replicate": replicate,
            "breadth": 0.9 if readable else 0.01,
            "coverage_threshold_passed": readable,
        }
        if with_time:
            row["time"] = time
        rows.append(row)
    pd.DataFrame(rows, columns=columns).to_csv(
        folder / f"{mag}_QC.tsv", sep="\t", index=False
    )


def _folders(tmp_path, samples_by_mag):
    """One QC folder per timepoint (``pre``, ``post``) from ``{mag: [samples]}``.

    Every genome gets a table in BOTH folders, header-only where it has no
    sample, as the QC step guarantees.  Returns the two folder paths as text,
    earlier timepoint first.
    """
    out = []
    for time in ("pre", "post"):
        folder = tmp_path / f"QC_{time}"
        for mag, samples in samples_by_mag.items():
            _write_qc(folder, mag, [s for s in samples if s[3] == time])
        out.append(str(folder))
    return out


def _row(table, mag="MAG_A"):
    """The eligibility row of one genome as a plain dict."""
    rows = table[table["MAG_ID"] == mag]
    assert len(rows) == 1, f"expected one row for {mag}, got {len(rows)}"
    return rows.iloc[0].to_dict()


# ---------------------------------------------------------------------------
# Counting
# ---------------------------------------------------------------------------


def test_unpaired_design_counts_cages_of_subjects_readable_at_both_timepoints(tmp_path):
    samples = (
        # Group A: cage c1 holds two subjects, cage c2 one.  All readable at
        # both timepoints -> 2 cages.
        _subject("a1", "A", "c1")
        + _subject("a2", "A", "c1")
        + _subject("a3", "A", "c2")
        # Group B: b1 and b2 are readable at both timepoints (cages c3, c4).
        + _subject("b1", "B", "c3")
        + _subject("b2", "B", "c4")
        # b3 fails QC at post and b4 has no post sample: neither counts, so
        # cage c5 does not count.
        + _subject("b3", "B", "c5", post=False)
        + _subject("b4", "B", "c5", post=None)
        # Cage c6 is readable at pre through b5 and at post through b6, but no
        # single subject is readable at both -> the cage does not count.
        + _subject("b5", "B", "c6", post=None)
        + _subject("b6", "B", "c6", pre=None)
    )
    table = analyze_qc_dirs(
        _folders(tmp_path, {"MAG_A": samples}), 2, "longitudinal", ["A", "B"]
    )
    row = _row(table)
    assert row["replicates_per_group_A"] == 2
    assert row["replicates_per_group_B"] == 2
    # No replicate ID is shared between the groups -> nothing is paired.
    assert row["paired_replicates_per_group_A"] == 0
    assert row["paired_replicates_per_group_B"] == 0
    assert row["unpaired_test_eligible"]
    assert not row["paired_test_eligible"]
    assert row["single_sample_eligible_A"] and row["single_sample_eligible_B"]


def test_two_group_paired_design_counts_the_shared_replicates(tmp_path):
    # Four blocks r1..r4, each with one subject in A and one in B.
    samples = []
    for i in range(1, 5):
        samples += _subject(f"a{i}", "A", f"r{i}") + _subject(f"b{i}", "B", f"r{i}")
    # A fifth block exists in A only: it counts for A, but it is not paired.
    samples += _subject("a5", "A", "r5")
    table = analyze_qc_dirs(
        _folders(tmp_path, {"MAG_A": samples}), 4, "longitudinal", ["A", "B"]
    )
    row = _row(table)
    assert row["replicates_per_group_A"] == 5
    assert row["replicates_per_group_B"] == 4
    assert row["paired_replicates_per_group_A"] == 4
    assert row["paired_replicates_per_group_B"] == 4
    assert row["paired_test_eligible"]


def test_three_groups_sharing_replicates_are_paired_for_the_tested_pair(tmp_path):
    # A, B and C all use blocks r1..r4.  Counting "replicates seen in exactly
    # two groups" over all three groups finds none (each is seen in three);
    # over the tested pair A-B every block is shared.
    samples = []
    for i in range(1, 5):
        for group in ("A", "B", "C"):
            samples += _subject(f"{group.lower()}{i}", group, f"r{i}")
    table = analyze_qc_dirs(
        _folders(tmp_path, {"MAG_A": samples}), 4, "longitudinal", ["A", "B"]
    )
    row = _row(table)
    assert row["paired_replicates_per_group_A"] == 4
    assert row["paired_replicates_per_group_B"] == 4
    assert row["paired_test_eligible"]
    # The third group leaves no trace in the table.
    assert not [c for c in table.columns if c.endswith("_C")]


def test_pair_without_shared_replicates_is_zero_not_an_error_when_another_pair_shares(
    tmp_path,
):
    # A and B share blocks r1..r4; C has blocks of its own.  For the pair A-C
    # nothing is shared: the paired count is 0 for both, and the unpaired test
    # is still eligible.  (Counting over all three groups marked A's blocks as
    # paired through B, left C without any, and stopped the job.)
    samples = []
    for i in range(1, 5):
        samples += _subject(f"a{i}", "A", f"r{i}") + _subject(f"b{i}", "B", f"r{i}")
        samples += _subject(f"c{i}", "C", f"x{i}")
    table = analyze_qc_dirs(
        _folders(tmp_path, {"MAG_A": samples}), 4, "longitudinal", ["A", "C"]
    )
    row = _row(table)
    assert row["replicates_per_group_A"] == 4
    assert row["replicates_per_group_C"] == 4
    assert row["paired_replicates_per_group_A"] == 0
    assert row["paired_replicates_per_group_C"] == 0
    assert row["unpaired_test_eligible"]
    assert not row["paired_test_eligible"]


def test_pair_and_count_marks_subjects_over_all_groups_then_keeps_the_tested_pair():
    # The in-memory step alone.  s3 is in the third group: it is marked like
    # everyone else, then dropped because its group is not tested.
    frame = pd.DataFrame(
        {
            "sample_id": [
                "s1_pre",
                "s1_post",
                "s2_pre",
                "s2_post",
                "s3_pre",
                "s3_post",
            ],
            "MAG_ID": ["MAG_A"] * 6,
            "group": ["A", "A", "B", "B", "C", "C"],
            "subjectID": ["s1", "s1", "s2", "s2", "s3", "s3"],
            "replicate": ["r1", "r1", "r1", "r1", "r1", "r1"],
            "coverage_threshold_passed": [True, True, True, False, True, True],
            "time": ["pre", "post", "pre", "post", "pre", "post"],
        }
    )
    out = pair_and_count(frame, ["A", "B"], "longitudinal")
    assert set(out["group"]) == {"A", "B"}
    by_sample = out.set_index("sample_id")
    # s1 is readable at both timepoints; s2 failed at post, so neither of its
    # samples counts.
    assert by_sample.loc["s1_pre", "two_timepoints_passed"]
    assert not by_sample.loc["s2_pre", "two_timepoints_passed"]
    assert by_sample.loc["s1_pre", "replicates_per_group"] == 1
    assert np.isnan(by_sample.loc["s2_pre", "replicates_per_group"])
    # B has no counting subject, so block r1 is not shared.
    assert by_sample.loc["s1_pre", "paired_replicates_per_group"] == 0


def test_longitudinal_rows_without_a_time_column_are_refused():
    frame = pd.DataFrame(
        {
            "sample_id": ["s1"],
            "MAG_ID": ["MAG_A"],
            "group": ["A"],
            "subjectID": ["s1"],
            "replicate": ["r1"],
            "coverage_threshold_passed": [True],
        }
    )
    with pytest.raises(ValueError, match="time"):
        pair_and_count(frame, ["A", "B"], "longitudinal")


# ---------------------------------------------------------------------------
# Which genomes reach the table
# ---------------------------------------------------------------------------


def _four_cages(group, prefix):
    """Four subjects of one group in four cages, readable at both timepoints."""
    samples = []
    for i in range(1, 5):
        samples += _subject(f"{prefix}{i}", group, f"{prefix}c{i}")
    return samples


def test_genome_with_rows_at_one_timepoint_only_is_left_out(tmp_path):
    full = _four_cages("A", "a") + _four_cages("B", "b")
    # MAG_B has samples at pre only: its post table is header-only, so no
    # subject is readable at both timepoints.
    pre_only = [s for s in full if s[3] == "pre"]
    table = analyze_qc_dirs(
        _folders(tmp_path, {"MAG_A": full, "MAG_B": pre_only}),
        4,
        "longitudinal",
        ["A", "B"],
    )
    assert list(table["MAG_ID"]) == ["MAG_A"]


def test_no_genome_with_both_tested_groups_is_an_error(tmp_path):
    # Only group A has samples -> nothing to compare, and an empty eligibility
    # table would let the run carry on with no genome.
    with pytest.raises(ValueError, match="No MAGs passed"):
        analyze_qc_dirs(
            _folders(tmp_path, {"MAG_A": _four_cages("A", "a")}),
            4,
            "longitudinal",
            ["A", "B"],
        )


def test_rows_are_sorted_by_genome_and_the_tested_pair_keeps_its_order(tmp_path):
    samples = _four_cages("A", "a") + _four_cages("B", "b")
    folders = _folders(tmp_path, {"MAG_C": samples, "MAG_A": samples, "MAG_B": samples})
    table = analyze_qc_dirs(folders, 4, "longitudinal", ["B", "A"])
    # Genomes in name order whatever the folder listing order.
    assert list(table["MAG_ID"]) == ["MAG_A", "MAG_B", "MAG_C"]
    # The single-sample columns follow the pair as given (B then A), not the
    # order the samples happen to sit in the QC tables.
    single = [c for c in table.columns if c.startswith("single_sample_eligible_")]
    assert single == ["single_sample_eligible_B", "single_sample_eligible_A"]


# ---------------------------------------------------------------------------
# Reading the tables
# ---------------------------------------------------------------------------


def test_header_only_table_contributes_nothing(tmp_path):
    pre_only = [s for s in _four_cages("A", "a") if s[3] == "pre"]
    folders = _folders(tmp_path, {"MAG_A": pre_only, "MAG_B": []})
    rows = load_genome_qc(folders, "MAG_A")
    assert len(rows) == 4
    assert set(rows["time"]) == {"pre"}
    # The verdict is a real True/False column although one table had no row.
    assert rows["coverage_threshold_passed"].dtype == bool
    # A genome whose every table is header-only: no row, same columns.
    empty = load_genome_qc(folders, "MAG_B")
    assert len(empty) == 0
    assert list(empty.columns) == list(rows.columns)


def test_table_missing_from_one_folder_is_an_error(tmp_path):
    # Every genome has a table at every timepoint (header-only when it has no
    # sample there), so a missing table means a failed step, not "no sample".
    folders = _folders(tmp_path, {"MAG_A": _four_cages("A", "a")})
    (tmp_path / "QC_post" / "MAG_A_QC.tsv").unlink()
    with pytest.raises(FileNotFoundError, match="MAG_A_QC.tsv"):
        load_genome_qc(folders, "MAG_A")


def test_same_folder_given_twice_is_refused(tmp_path):
    # The same samples read twice would look like subjects with two timepoints.
    folders = _folders(tmp_path, {"MAG_A": _four_cages("A", "a")})
    with pytest.raises(ValueError, match="more than once"):
        load_genome_qc([folders[0], folders[0]], "MAG_A")


def test_verdict_that_is_not_true_or_false_is_refused(tmp_path):
    folders = _folders(tmp_path, {"MAG_A": _four_cages("A", "a")})
    path = tmp_path / "QC_pre" / "MAG_A_QC.tsv"
    path.write_text(path.read_text().replace("True", "yes", 1))
    with pytest.raises(ValueError, match="coverage_threshold_passed"):
        load_genome_qc(folders, "MAG_A")


def test_numeric_group_names_are_matched_as_text(tmp_path):
    # Groups named 20 and 40, cages named 1..8: read with pandas' defaults
    # they become numbers and "--groups 20 40" (text) matches no row.
    samples = []
    for i in range(1, 5):
        samples += _subject(f"m{i}", "20", str(i))
        samples += _subject(f"m{i + 4}", "40", str(i + 4))
    table = analyze_qc_dirs(
        _folders(tmp_path, {"MAG_A": samples}), 4, "longitudinal", ["20", "40"]
    )
    row = _row(table)
    assert row["replicates_per_group_20"] == 4
    assert row["replicates_per_group_40"] == 4
    assert row["unpaired_test_eligible"]


def test_single_data_type_uses_one_folder_and_the_verdict_as_is(tmp_path):
    # One timepoint, no time column: a sample counts when it passed QC.
    folder = tmp_path / "QC_t1"
    samples = [
        ("a1", "A", "c1", None, True),
        ("a2", "A", "c2", None, True),
        ("a3", "A", "c3", None, False),  # failed QC -> cage c3 does not count
        ("b1", "B", "c4", None, True),
        ("b2", "B", "c5", None, True),
    ]
    _write_qc(folder, "MAG_A", samples, with_time=False)
    table = analyze_qc_dirs([str(folder)], 2, "single", ["A", "B"])
    row = _row(table)
    assert row["replicates_per_group_A"] == 2
    assert row["replicates_per_group_B"] == 2
    assert row["unpaired_test_eligible"]
    # The single-sample test is a longitudinal test: blank for single data.
    assert np.isnan(row["single_sample_eligible_A"])
    assert np.isnan(row["single_sample_eligible_B"])


def test_wrong_number_of_folders_is_refused(tmp_path):
    folders = _folders(tmp_path, {"MAG_A": _four_cages("A", "a")})
    # Longitudinal: one folder per timepoint, so exactly two.
    with pytest.raises(ValueError, match="exactly 2"):
        analyze_qc_dirs(folders[:1], 4, "longitudinal", ["A", "B"])
    with pytest.raises(ValueError, match="exactly 2"):
        analyze_qc_dirs(folders + folders[:1], 4, "longitudinal", ["A", "B"])
    # Single: exactly one.
    with pytest.raises(ValueError, match="exactly 1"):
        analyze_qc_dirs(folders, 4, "single", ["A", "B"])


def test_permuted_metadata_relabels_groups_before_counting(tmp_path):
    # Real labels: A has cages c1..c3, B has cage c4 only -> 3 and 1.
    samples = (
        _subject("m1", "A", "c1")
        + _subject("m2", "A", "c2")
        + _subject("m3", "A", "c3")
        + _subject("m4", "B", "c4")
    )
    folders = _folders(tmp_path, {"MAG_A": samples})
    real = analyze_qc_dirs(folders, 2, "longitudinal", ["A", "B"])
    assert _row(real)["replicates_per_group_A"] == 3
    assert not _row(real)["unpaired_test_eligible"]
    # The permuted sheet moves subject m3 (both of its samples) to group B:
    # the counts must be taken on the permuted labels -> 2 and 2.
    sheet = pd.DataFrame(
        {
            "sample_id": [f"{s}_{t}" for s, _, _, t, _ in samples],
            "group": ["B" if s == "m3" else g for s, g, _, _, _ in samples],
        }
    )
    permuted = tmp_path / "permuted_metadata.tsv"
    sheet.to_csv(permuted, sep="\t", index=False)
    null = analyze_qc_dirs(
        folders, 2, "longitudinal", ["A", "B"], permuted_metadata=str(permuted)
    )
    row = _row(null)
    assert row["replicates_per_group_A"] == 2
    assert row["replicates_per_group_B"] == 2
    assert row["unpaired_test_eligible"]


# ---------------------------------------------------------------------------
# Brute force: the script against a plain definition on random designs
# ---------------------------------------------------------------------------


def _random_design(rng, n_groups):
    """Samples of one genome for a random design of ``n_groups`` groups.

    In 70 % of the designs the groups reuse the same replicate IDs (a paired
    design); subjects miss a sample at a timepoint 15 % of the time and a
    sample fails QC 25 % of the time.
    """
    shared = rng.random() < 0.7
    samples = []
    for g in range(n_groups):
        group = f"g{g}"
        for r in range(rng.randint(1, 6)):
            replicate = f"r{r}" if shared else f"{group}r{r}"
            for m in range(rng.randint(1, 3)):
                subject = f"{group}_{replicate}_m{m}"
                for time in ("pre", "post"):
                    if rng.random() < 0.85:
                        samples.append(
                            (subject, group, replicate, time, rng.random() < 0.75)
                        )
    return samples


def _plain_definition(samples, pair, floor):
    """The expected eligibility row, from sets; None when the genome is left out.

    A subject counts when it has a readable sample at both timepoints.  A
    group's replicates are those of its counting subjects.  Paired replicates
    are the ones both tested groups have.
    """
    times_readable = {}
    for subject, _group, _replicate, time, readable in samples:
        if readable:
            times_readable.setdefault(subject, set()).add(time)
    counting = {s for s, times in times_readable.items() if times == {"pre", "post"}}
    cages = {g: set() for g in pair}
    for subject, group, replicate, _time, _readable in samples:
        if subject in counting and group in cages:
            cages[group].add(replicate)
    if any(not c for c in cages.values()):
        return None
    g1, g2 = pair
    both = cages[g1] & cages[g2]
    expected = {
        "unpaired_test_eligible": len(cages[g1]) >= floor and len(cages[g2]) >= floor,
        "paired_test_eligible": len(both) >= floor,
    }
    for g in pair:
        expected[f"replicates_per_group_{g}"] = len(cages[g])
        expected[f"paired_replicates_per_group_{g}"] = len(both)
        expected[f"single_sample_eligible_{g}"] = len(cages[g]) >= floor
    return expected


def test_matches_the_plain_definition_on_random_designs(tmp_path):
    rng = random.Random(20261113)
    pair, floor = ["g0", "g1"], 3
    designs = {
        f"MAG_{i:03d}": _random_design(rng, rng.randint(2, 5)) for i in range(300)
    }
    expected = {
        mag: _plain_definition(samples, pair, floor) for mag, samples in designs.items()
    }
    kept = sorted(mag for mag, row in expected.items() if row is not None)
    # The seed must exercise every branch, or the comparison proves little.
    assert 0 < len(kept) < len(designs)
    assert any(expected[m]["paired_test_eligible"] for m in kept)
    assert any(expected[m]["paired_replicates_per_group_g0"] == 0 for m in kept)

    table = analyze_qc_dirs(_folders(tmp_path, designs), floor, "longitudinal", pair)

    assert list(table["MAG_ID"]) == kept
    for mag in kept:
        row = _row(table, mag)
        for column, value in expected[mag].items():
            assert row[column] == value, f"{mag} {column}: {row[column]} != {value}"


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def test_cli_reads_two_folders_and_the_old_single_folder_flag_is_gone(tmp_path):
    folders = _folders(
        tmp_path, {"MAG_A": _four_cages("A", "a") + _four_cages("B", "b")}
    )
    out = tmp_path / "eligibility.tsv"
    base = ["alleleflux-eligibility", "--min_sample_num", "4", "--groups", "A", "B"]
    done = subprocess.run(
        base + ["--qc_dirs", *folders, "--output_file", str(out)],
        capture_output=True,
        text=True,
    )
    assert done.returncode == 0, done.stderr
    table = pd.read_csv(out, sep="\t")
    assert list(table["MAG_ID"]) == ["MAG_A"]
    assert table.loc[0, "replicates_per_group_A"] == 4
    assert bool(table.loc[0, "unpaired_test_eligible"])
    # --qc_dir named ONE folder holding both timepoints; that layout is gone.
    old = subprocess.run(
        base + ["--qc_dir", folders[0], "--output_file", str(tmp_path / "old.tsv")],
        capture_output=True,
        text=True,
    )
    assert old.returncode != 0
    assert not (tmp_path / "old.tsv").exists()


# ---------------------------------------------------------------------------
# The three counting functions (moved here from the QC tests with the functions)
# ---------------------------------------------------------------------------


class TestCheckTimepoints(unittest.TestCase):

    def test_single_data_type(self):
        """Test single data type skips timepoint checks."""
        df = pd.DataFrame(
            {
                "subjectID": ["subject1", "subject2"],
                "coverage_threshold_passed": [True, False],
            }
        )

        result = check_timepoints(df, "single")

        # Should match coverage_threshold_passed
        self.assertTrue(result.loc[0, "two_timepoints_passed"])
        self.assertFalse(result.loc[1, "two_timepoints_passed"])

    def test_longitudinal_valid_timepoints(self):
        """Test longitudinal data with valid timepoints."""
        df = pd.DataFrame(
            {
                "subjectID": ["subject1", "subject1", "subject2", "subject2"],
                "time": ["T1", "T2", "T1", "T2"],
                "coverage_threshold_passed": [True, True, True, True],
            }
        )

        result = check_timepoints(df, "longitudinal")

        # All should pass since both subjects have 2 timepoints
        self.assertTrue(all(result["two_timepoints_passed"]))

    def test_longitudinal_missing_timepoint(self):
        """Test longitudinal data with missing timepoints."""
        df = pd.DataFrame(
            {
                "subjectID": [
                    "subject1",
                    "subject1",
                    "subject2",
                ],  # subject2 missing T2
                "time": ["T1", "T2", "T1"],
                "coverage_threshold_passed": [True, True, True],
            }
        )

        result = check_timepoints(df, "longitudinal")

        # Only subject1 should pass (has both timepoints)
        subject1_mask = result["subjectID"] == "subject1"
        subject2_mask = result["subjectID"] == "subject2"

        self.assertTrue(all(result.loc[subject1_mask, "two_timepoints_passed"]))
        self.assertFalse(any(result.loc[subject2_mask, "two_timepoints_passed"]))

    def test_too_many_timepoints(self):
        """Test error when more than 2 timepoints found."""
        df = pd.DataFrame(
            {
                "subjectID": ["subject1", "subject1", "subject1"],
                "time": ["T1", "T2", "T3"],  # 3 timepoints
                "coverage_threshold_passed": [True, True, True],
            }
        )

        with self.assertRaises(ValueError) as context:
            check_timepoints(df, "longitudinal")

        self.assertIn("More than 2 unique timepoints", str(context.exception))


class TestSubjectCounting(unittest.TestCase):

    def test_add_subject_count_per_group(self):
        """Test subject and replicate counting per group."""
        df = pd.DataFrame(
            {
                "two_timepoints_passed": [True, True, True, False],
                "group": ["treatment", "treatment", "control", "control"],
                "subjectID": ["subject1", "subject2", "subject1", "subject3"],
                "replicate": ["rep1", "rep2", "rep1", "rep3"],
            }
        )

        result = add_subject_count_per_group(df)

        # Treatment group: 2 subjects, 2 replicates
        treatment_mask = (result["group"] == "treatment") & result[
            "two_timepoints_passed"
        ]
        self.assertTrue(all(result.loc[treatment_mask, "subjects_per_group"] == 2))
        self.assertTrue(all(result.loc[treatment_mask, "replicates_per_group"] == 2))

        # Control group: 1 subject, 1 replicate (subject3 failed timepoints)
        control_mask = (result["group"] == "control") & result["two_timepoints_passed"]
        self.assertTrue(all(result.loc[control_mask, "subjects_per_group"] == 1))
        self.assertTrue(all(result.loc[control_mask, "replicates_per_group"] == 1))

        # Failed samples should have NaN
        failed_mask = ~result["two_timepoints_passed"]
        self.assertTrue(all(result.loc[failed_mask, "subjects_per_group"].isna()))

    def test_count_paired_replicates(self):
        """Test paired replicate counting."""
        df = pd.DataFrame(
            {
                "two_timepoints_passed": [True, True, True, True, True],
                "group": ["treatment", "control", "treatment", "control", "treatment"],
                "replicate": [
                    "rep1",
                    "rep1",
                    "rep2",
                    "rep2",
                    "rep3",
                ],  # rep1 and rep2 are paired
            }
        )

        result = count_paired_replicates(df)

        # rep1 and rep2 appear in both groups, so 2 paired replicates per group
        paired_mask = result["replicate"].isin(["rep1", "rep2"])
        self.assertTrue(
            all(result.loc[paired_mask, "paired_replicates_per_group"] == 2)
        )

        # rep3 only in treatment, should be NaN
        unpaired_mask = result["replicate"] == "rep3"
        self.assertTrue(
            all(result.loc[unpaired_mask, "paired_replicates_per_group"].isna())
        )


class TestCountingEdgeCases(unittest.TestCase):

    def test_no_valid_subjects(self):
        """Test subject counting with no valid subjects."""
        df = pd.DataFrame(
            {
                "two_timepoints_passed": [False, False],
                "group": ["treatment", "control"],
                "subjectID": ["subject1", "subject2"],
                "replicate": ["rep1", "rep2"],
            }
        )

        result = add_subject_count_per_group(df)

        # All should be NaN since no subjects passed
        self.assertTrue(all(result["subjects_per_group"].isna()))
        self.assertTrue(all(result["replicates_per_group"].isna()))

    def test_single_data_type_counts_the_samples_that_passed_qc(self):
        """Single data: the three steps in a row, with no time column at all."""
        frame = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2", "sample3"],
                "subjectID": ["subject1", "subject2", "subject3"],
                "group": ["treatment", "control", "control"],
                "replicate": ["rep1", "rep2", "rep3"],
                "coverage_threshold_passed": [True, True, False],
            }
        )

        marked = check_timepoints(frame, "single")
        # A sample counts exactly when it passed QC.
        self.assertEqual(list(marked["two_timepoints_passed"]), [True, True, False])

        counted = count_paired_replicates(add_subject_count_per_group(marked))
        by_sample = counted.set_index("sample_id")
        # One counting replicate in each group; sample3 failed QC, so its
        # replicate does not count and its own cells stay blank.
        self.assertEqual(by_sample.loc["sample1", "replicates_per_group"], 1)
        self.assertEqual(by_sample.loc["sample2", "replicates_per_group"], 1)
        self.assertTrue(np.isnan(by_sample.loc["sample3", "replicates_per_group"]))
        # rep1 and rep2 each sit in one group only: nothing is paired.
        self.assertEqual(by_sample.loc["sample1", "paired_replicates_per_group"], 0)
        self.assertEqual(by_sample.loc["sample2", "paired_replicates_per_group"], 0)
