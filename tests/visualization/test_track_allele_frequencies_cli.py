"""End-to-end test of the ``alleleflux-track-alleles`` console script.

The unit tests in ``test_track_allele_frequencies.py`` cover the worker
functions one at a time.  This file runs the whole CLI as a subprocess on a
synthetic allele-frequency cache and checks BOTH tables it writes
(``<mag>_frequency_table.tsv`` wide, ``<mag>_frequency_table.long.tsv`` long)
value by value.

Fixture geometry (one MAG, two timepoints, two mice, three anchor sites):

    metadata                cache rows at the three anchor sites
    --------                -----------------------------------------------
    s1  m1  pre   day 0     c1:100 A=8 C=2 (0.8) | c1:200 G=6 T=4 (0.6) | c2:50 T=7 A=3 (0.7)
    s2  m2  pre   day 0     c1:100 A=5 C=5 (0.5) | c1:200 G=3     (depth 3 < 5 -> NaN) | c2:50 T=10 (1.0)
    s3  m1  post  day 10    c1:100 A=2 C=8 (0.2) | c1:200 G=10    (1.0) | c2:50 ABSENT -> NaN
    s4  m2  post  day 10    c1:100 A=1 C=9 (0.1) | c1:200 T=10    (0.0) | c2:50 T=4 A=6 (0.4)
    s9  m3  post  day 10    NOT in the cache at all -> dropped with a warning, no column

Anchor file: the three sites above (anchor letters A, G, T) plus one extra row
(c2:75) whose anchor is NaN, which the loader must drop before any work.

So the wide table has 3 rows (sites) x 4 sample columns, and the long table has
3 x 4 = 12 rows with ``group``/``time``/``subjectID``/``day`` merged in from
the metadata.  Every NaN in the output has a distinct, deliberate cause:
depth floor (s2 @ c1:200), missing row (s3 @ c2:50), and the whole-sample case
(s9) which must not surface as a column at all.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tests.utilities.test_cache_io import write_cache

MAG = "MAG_A"
ANCHOR_COL = "terminal_nucleotide_mean_freq"
MIN_COV = 5  # --min-cov-per-site: s2's c1:200 row has depth 3, so it falls below

# (contig, position, anchor letter, q_value) -- one row per anchor site in the
# anchor file.  The q_value column is carried into both output tables, so we
# use it to prove the "p-value column passthrough" path too.
ANCHOR_SITES = [
    ("c1", 100, "A", 0.001),
    ("c1", 200, "G", 0.010),
    ("c2", 50, "T", 0.020),
]

# Metadata roster: sample_id, group, time, subjectID, day.
# s9 is the decoy: listed here, absent from the cache.
METADATA = [
    ("s1", "gA", "pre", "m1", 0),
    ("s2", "gB", "pre", "m2", 0),
    ("s3", "gA", "post", "m1", 10),
    ("s4", "gB", "post", "m2", 10),
    ("s9", "gB", "post", "m3", 10),
]

# Cache rows per timepoint in write_cache's tuple layout:
# (contig, position, sample_id, group, subjectID, A, C, G, T)
CACHE_ROWS = {
    "pre": [
        ("c1", 100, "s1", "gA", "m1", 8, 2, 0, 0),  # anchor A -> 8/10 = 0.8
        ("c1", 200, "s1", "gA", "m1", 0, 0, 6, 4),  # anchor G -> 6/10 = 0.6
        ("c2", 50, "s1", "gA", "m1", 3, 0, 0, 7),  # anchor T -> 7/10 = 0.7
        ("c1", 100, "s2", "gB", "m2", 5, 5, 0, 0),  # anchor A -> 5/10 = 0.5
        ("c1", 200, "s2", "gB", "m2", 0, 0, 3, 0),  # depth 3 < MIN_COV -> NaN
        ("c2", 50, "s2", "gB", "m2", 0, 0, 0, 10),  # anchor T -> 10/10 = 1.0
    ],
    "post": [
        ("c1", 100, "s3", "gA", "m1", 2, 8, 0, 0),  # anchor A -> 2/10 = 0.2
        ("c1", 200, "s3", "gA", "m1", 0, 0, 10, 0),  # anchor G -> 10/10 = 1.0
        # s3 has NO row at c2:50 -> NaN there (sample itself is in the cache)
        ("c1", 100, "s4", "gB", "m2", 1, 9, 0, 0),  # anchor A -> 1/10 = 0.1
        ("c1", 200, "s4", "gB", "m2", 0, 0, 0, 10),  # anchor G -> 0/10 = 0.0
        ("c2", 50, "s4", "gB", "m2", 6, 0, 0, 4),  # anchor T -> 4/10 = 0.4
        # s9: no rows at all -> dropped with a warning, never a column
    ],
}

# Ground truth, keyed by (contig, position) -> {sample_id: frequency}.
EXPECTED = {
    ("c1", 100): {"s1": 0.8, "s2": 0.5, "s3": 0.2, "s4": 0.1},
    ("c1", 200): {"s1": 0.6, "s2": np.nan, "s3": 1.0, "s4": 0.0},
    ("c2", 50): {"s1": 0.7, "s2": 1.0, "s3": np.nan, "s4": 0.4},
}
PRESENT_SAMPLES = ["s1", "s2", "s3", "s4"]
ABSENT_SAMPLE = "s9"


def _cli() -> list[str]:
    """The console script from the SAME env as pytest, else ``python -m`` fallback.

    Resolving next to ``sys.executable`` guards against a stale copy earlier on
    PATH (another interpreter's install) being picked up by accident.
    """
    exe = Path(sys.executable).parent / "alleleflux-track-alleles"
    if exe.exists():
        return [str(exe)]
    return [
        sys.executable,
        "-m",
        "alleleflux.scripts.visualization.accessory.track_allele_frequencies",
    ]


def _build_fixture(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    """Write cache, metadata and anchor file under ``tmp_path``; return their paths."""
    cache_dir = tmp_path / "allele_freq_cache"
    for timepoint, rows in CACHE_ROWS.items():
        # Layout the pipeline uses: {cache_dir}/{tp}/{mag}_{tp}_allele_frequency.parquet
        write_cache(str(cache_dir), MAG, timepoint, rows)

    metadata = tmp_path / "metadata.tsv"
    pd.DataFrame(
        METADATA, columns=["sample_id", "group", "time", "subjectID", "day"]
    ).to_csv(metadata, sep="\t", index=False)

    # Anchor file shaped like terminal_nucleotide_analysis output: mag_id,
    # contig, position, gene_id, the anchor column, q_value.  The last row has
    # a NaN anchor and must be dropped by load_and_prepare_anchor_sites.
    anchor_rows = [(MAG, c, p, f"gene_{c}", base, q) for c, p, base, q in ANCHOR_SITES]
    anchor_rows.append((MAG, "c2", 75, "gene_c2", np.nan, 0.030))
    anchor = tmp_path / "anchors.tsv"
    pd.DataFrame(
        anchor_rows,
        columns=["mag_id", "contig", "position", "gene_id", ANCHOR_COL, "q_value"],
    ).to_csv(anchor, sep="\t", index=False)

    out_dir = tmp_path / "track_out"
    return cache_dir, metadata, anchor, out_dir


@pytest.fixture(scope="module")
def tracker_run(tmp_path_factory):
    """Run the CLI once for the module; return (CompletedProcess, out_dir).

    ``--cpus 2`` with two timepoints means the per-timepoint multiprocessing
    pool really forks two workers (n_workers = min(cpus, n_timepoints)).
    """
    tmp_path = tmp_path_factory.mktemp("track_cli")
    cache_dir, metadata, anchor, out_dir = _build_fixture(tmp_path)
    cmd = _cli() + [
        "--mag-id", MAG,
        "--anchor-file", str(anchor),
        "--cache-dir", str(cache_dir),
        "--metadata", str(metadata),
        "--output-dir", str(out_dir),
        "--anchor-column", ANCHOR_COL,
        "--cpus", "2",
        "--min-cov-per-site", str(MIN_COV),
        "--log-level", "INFO",
    ]  # fmt: skip
    done = subprocess.run(cmd, capture_output=True, text=True)
    return done, out_dir


def _wide(out_dir: Path) -> pd.DataFrame:
    return pd.read_csv(out_dir / f"{MAG}_frequency_table.tsv", sep="\t")


def _long(out_dir: Path) -> pd.DataFrame:
    return pd.read_csv(out_dir / f"{MAG}_frequency_table.long.tsv", sep="\t")


class TestTrackAllelesCli:
    def test_exit_code_and_dropped_sample_warning(self, tracker_run):
        done, out_dir = tracker_run
        # stderr in the message so a failure shows the traceback, not just "1 != 0"
        assert done.returncode == 0, f"stderr:\n{done.stderr}"
        # The decoy sample is reported, by name, as having no cache rows.
        assert f"Sample {ABSENT_SAMPLE} has no rows in the cache" in done.stderr
        # The NaN-anchor row is counted and dropped before the workers run.
        assert f"Found 1 sites with NaN values in '{ANCHOR_COL}'" in done.stderr
        assert (out_dir / f"{MAG}_frequency_table.tsv").exists()
        assert (out_dir / f"{MAG}_frequency_table.long.tsv").exists()

    # ---- wide table ---------------------------------------------------------

    def test_wide_table_shape_and_columns(self, tracker_run):
        _, out_dir = tracker_run
        wide = _wide(out_dir)
        # 3 anchor sites survive (the NaN-anchor c2:75 is gone).
        assert len(wide) == 3
        assert set(zip(wide["contig"], wide["position"])) == set(EXPECTED)
        # Exact column set: 6 id columns, q_value (passthrough), the anchor
        # column itself (an "other" column the merge keeps), 4 sample columns,
        # and the n_samples_non_nan QC column.  Sample-column ORDER is not
        # asserted: it follows the order the timepoint workers return, which
        # imap_unordered does not fix.
        expected_cols = {
            "mag_id",
            "contig",
            "position",
            "gene_id",
            "anchor_allele",
            "anchor_criteria",
            "q_value",
            ANCHOR_COL,
            *PRESENT_SAMPLES,
            "n_samples_non_nan",
        }
        assert set(wide.columns) == expected_cols
        assert len(wide.columns) == 13
        # The absent sample is not a column -- it was dropped, not NaN-filled.
        assert ABSENT_SAMPLE not in wide.columns
        # id columns come first, in the documented order
        assert list(wide.columns[:7]) == [
            "mag_id",
            "contig",
            "position",
            "gene_id",
            "anchor_allele",
            "anchor_criteria",
            "q_value",
        ]
        assert list(wide.columns[-1:]) == ["n_samples_non_nan"]

    def test_wide_table_values(self, tracker_run):
        _, out_dir = tracker_run
        wide = _wide(out_dir).set_index(["contig", "position"])
        for contig, pos, base, q in ANCHOR_SITES:
            row = wide.loc[(contig, pos)]
            assert row["mag_id"] == MAG
            assert row["anchor_allele"] == base
            assert row["anchor_criteria"] == ANCHOR_COL
            assert row["gene_id"] == f"gene_{contig}"
            assert row["q_value"] == pytest.approx(q)
            for sample, want in EXPECTED[(contig, pos)].items():
                got = row[sample]
                if np.isnan(want):
                    assert np.isnan(
                        got
                    ), f"{contig}:{pos} {sample} should be NaN, got {got}"
                else:
                    assert got == pytest.approx(want), f"{contig}:{pos} {sample}"
        # QC column counts the non-NaN sample cells per site:
        #   c1:100 all four; c1:200 three (s2 below depth); c2:50 three (s3 no row)
        assert wide.loc[("c1", 100), "n_samples_non_nan"] == 4
        assert wide.loc[("c1", 200), "n_samples_non_nan"] == 3
        assert wide.loc[("c2", 50), "n_samples_non_nan"] == 3

    # ---- long table ---------------------------------------------------------

    def test_long_table_shape_and_columns(self, tracker_run):
        _, out_dir = tracker_run
        long = _long(out_dir)
        # 3 sites x 4 present samples; the absent sample contributes no rows.
        assert len(long) == 12
        assert set(long["sample_id"]) == set(PRESENT_SAMPLES)
        assert ABSENT_SAMPLE not in set(long["sample_id"])
        # Exact column list: id columns (q_value kept, the anchor column itself
        # is NOT carried into the long table), the melted pair, then the
        # metadata merge (group/time/subjectID always; day because the
        # metadata has it; no replicate because the metadata lacks it).
        assert list(long.columns) == [
            "mag_id",
            "contig",
            "position",
            "gene_id",
            "anchor_allele",
            "anchor_criteria",
            "q_value",
            "sample_id",
            "frequency",
            "group",
            "time",
            "subjectID",
            "day",
        ]
        assert ANCHOR_COL not in long.columns
        assert "replicate" not in long.columns

    def test_long_table_values_and_metadata_merge(self, tracker_run):
        _, out_dir = tracker_run
        long = _long(out_dir).set_index(["contig", "position", "sample_id"])
        meta = {s: (g, t, subj, d) for s, g, t, subj, d in METADATA}
        for (contig, pos), per_sample in EXPECTED.items():
            for sample, want in per_sample.items():
                row = long.loc[(contig, pos, sample)]
                got = row["frequency"]
                if np.isnan(want):
                    assert np.isnan(got), f"{contig}:{pos} {sample} should be NaN"
                else:
                    assert got == pytest.approx(want), f"{contig}:{pos} {sample}"
                # Metadata columns carried through the merge, per sample.
                group, time, subject, day = meta[sample]
                assert row["group"] == group
                assert row["time"] == time
                assert row["subjectID"] == subject
                assert row["day"] == day
        # Every (site, sample) pair appears exactly once.
        assert not long.index.duplicated().any()

    def test_wide_and_long_agree(self, tracker_run):
        """The long table is a pure melt of the wide one: same numbers, NaN for NaN."""
        _, out_dir = tracker_run
        wide = _wide(out_dir).set_index(["contig", "position"])
        long = _long(out_dir).set_index(["contig", "position", "sample_id"])
        for contig, pos in EXPECTED:
            for sample in PRESENT_SAMPLES:
                w = wide.loc[(contig, pos), sample]
                l = long.loc[(contig, pos, sample), "frequency"]
                assert (np.isnan(w) and np.isnan(l)) or w == pytest.approx(l)
