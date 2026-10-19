"""Tests for the ANI inputs: the QC sample gate, contig lengths, dense count arrays.

Fixtures mirror the production schema exactly: caches are written with the
columns the cache writer emits (one row group per sample), and QC fixtures use
the real column names from ``quality_control.py`` -- including a numeric group
label, to guard the "numeric labels silently become int64" trap.
"""

import os
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd

from alleleflux.scripts.analysis.ani.ani_io import (
    contig_lengths_for_mag,
    scatter_rows,
    dense_counts_from_cache,
    qc_passing_samples,
)


def _write_qc(path, mag_id, rows):
    """rows: (sample_id, group, subjectID, time, replicate, passed).

    Column names and their order match quality_control.py's real output (a subset).
    """
    pd.DataFrame(
        {
            "sample_id": [r[0] for r in rows],
            "MAG_ID": [mag_id] * len(rows),
            "file_path": ["/data/profiles/x.tsv.gz"] * len(rows),
            "group": [r[1] for r in rows],
            "subjectID": [r[2] for r in rows],
            "replicate": [r[4] for r in rows],
            "time": [r[3] for r in rows],
            # both averages, as quality_control.py writes them: the covered-sites
            # mean (the filter's reference) and the marker column the roster requires
            "average_coverage": [20.0 for _ in rows],
            "average_coverage_including_zeros": [20.0 for _ in rows],
            "coverage_threshold_passed": [r[5] for r in rows],
        }
    ).to_csv(path, sep="\t", index=False)


class TestQcPassingSamples(unittest.TestCase):
    MAG = "MAG_1"

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_keeps_only_passing_samples(self):
        qc = os.path.join(self.tmp, f"{self.MAG}_QC.tsv")
        _write_qc(
            qc,
            self.MAG,
            [
                ("S1", "treatment", "m1", "pre", 1, True),
                ("S2", "control", "m2", "pre", 1, False),
                ("S3", "treatment", "m3", "end", 1, True),
            ],
        )
        got = qc_passing_samples([qc], self.MAG)
        self.assertEqual(sorted(got["sample_id"]), ["S1", "S3"])

    def test_unions_across_files_without_duplicating(self):
        """A sample judged usable at two timepoint combinations appears once."""
        paths = []
        for name in ("qc_a.tsv", "qc_b.tsv"):
            path = os.path.join(self.tmp, name)
            _write_qc(path, self.MAG, [("S1", "treatment", "m1", "pre", 1, True)])
            paths.append(path)
        got = qc_passing_samples(paths, self.MAG)
        self.assertEqual(len(got), 1)

    def test_numeric_labels_stay_strings(self):
        """Some cohorts label groups and subjects with bare numbers; if pandas sniffs
        them as int64 every downstream string comparison silently matches nothing."""
        qc = os.path.join(self.tmp, "numeric_QC.tsv")
        _write_qc(
            qc,
            self.MAG,
            [
                ("S1", 10, 101, 4, 1, True),
                ("S2", 30, 102, 4, 1, True),
            ],
        )
        got = qc_passing_samples([qc], self.MAG)
        self.assertEqual(sorted(got["group"]), ["10", "30"])
        self.assertEqual(sorted(got["subjectID"]), ["101", "102"])
        self.assertTrue(all(isinstance(v, str) for v in got["time"]))

    def test_zero_passing_samples_is_an_empty_frame_not_an_error(self):
        """A low-abundance MAG where nothing passes QC is a legitimate outcome --
        the CLI writes header-only outputs, the pipeline stays green."""
        qc = os.path.join(self.tmp, "none_QC.tsv")
        _write_qc(qc, self.MAG, [("S1", "treatment", "m1", "pre", 1, False)])
        got = qc_passing_samples([qc], self.MAG)
        self.assertTrue(got.empty)
        self.assertIn("sample_id", got.columns)

    def test_raises_when_verdict_column_missing(self):
        qc = os.path.join(self.tmp, "bad_QC.tsv")
        pd.DataFrame({"sample_id": ["S1"]}).to_csv(qc, sep="\t", index=False)
        with self.assertRaises(ValueError):
            qc_passing_samples([qc], self.MAG)


def _dense_into(table, contig_lengths):
    """Allocate every contig at zero and scatter ``table`` in, as dense_counts_from_cache does."""
    dense = {c: np.zeros((n, 4), dtype=np.uint16) for c, n in contig_lengths.items()}
    scatter_rows(dense, table, contig_lengths)
    return dense


class TestScatterRows(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _load(self, rows):
        """rows: (contig, position, ref, A, C, G, T, n) -> the frame the cache gives the loader."""
        return pd.DataFrame(
            [
                {
                    "contig": c,
                    "position": pos,
                    "A": a,
                    "C": cc,
                    "G": g,
                    "T": t,
                    "coverage": a + cc + g + t,
                }
                for c, pos, _, a, cc, g, t, _n in rows
            ]
        )

    def test_places_rows_and_leaves_gaps_at_zero(self):
        """An absent position means zero coverage, not missing data."""
        table = self._load(
            [("ctg1", 0, "A", 4, 0, 0, 0, 0), ("ctg1", 3, "C", 0, 6, 0, 0, 0)]
        )
        dense = _dense_into(table, {"ctg1": 5})
        self.assertEqual(dense["ctg1"].shape, (5, 4))
        self.assertEqual(int(dense["ctg1"][0, 0]), 4)  # A at position 0
        self.assertEqual(int(dense["ctg1"][3, 1]), 6)  # C at position 3
        self.assertEqual(int(dense["ctg1"][1].sum()), 0)  # uncovered position
        self.assertEqual(dense["ctg1"].dtype, np.uint16)

    def test_every_mag_contig_gets_an_array_even_if_uncovered(self):
        """All samples must share one position index per contig, so a contig this
        sample never touched still gets a zero array."""
        table = self._load([("ctg1", 0, "A", 4, 0, 0, 0, 0)])
        dense = _dense_into(table, {"ctg1": 2, "ctg2": 3})
        self.assertEqual(dense["ctg2"].shape, (3, 4))
        self.assertEqual(int(dense["ctg2"].sum()), 0)

    def test_rejects_positions_beyond_contig_length(self):
        table = self._load([("ctg1", 99, "A", 4, 0, 0, 0, 0)])
        with self.assertRaises(ValueError):
            _dense_into(table, {"ctg1": 5})

    def test_rejects_contigs_the_mag_does_not_own(self):
        """A profile row on a contig outside the MAG's mapping means the mapping and
        the profiles are out of sync -- fail loud, don't silently drop data."""
        table = self._load(
            [("ctg1", 0, "A", 4, 0, 0, 0, 0), ("ghost", 0, "A", 4, 0, 0, 0, 0)]
        )
        with self.assertRaises(ValueError):
            _dense_into(table, {"ctg1": 5})


class TestContigLengthsForMag(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.mapping = os.path.join(self.tmp, "map.tsv")
        pd.DataFrame(
            {"mag_id": ["MAG_1", "MAG_1", "MAG_2"], "contig_id": ["c1", "c2", "other"]}
        ).to_csv(self.mapping, sep="\t", index=False)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _fasta(self, name, records, with_fai=False):
        path = os.path.join(self.tmp, name)
        with open(path, "w") as handle:
            for contig, seq in records:
                handle.write(f">{contig}\n{seq}\n")
        if with_fai:
            with open(path + ".fai", "w") as handle:
                offset = 0
                for contig, seq in records:
                    offset += len(contig) + 2  # ">" + name + newline
                    handle.write(
                        f"{contig}\t{len(seq)}\t{offset}\t{len(seq)}\t{len(seq) + 1}\n"
                    )
                    offset += len(seq) + 1
        return path

    def test_reads_lengths_from_the_fai_index(self):
        """The .fai beside the FASTA is a 1.5 MB text file versus a 0.5 GB parse."""
        fasta = self._fasta(
            "ref.fa", [("c1", "ACGT"), ("c2", "ACGTAC"), ("other", "AC")], with_fai=True
        )
        self.assertEqual(
            contig_lengths_for_mag(fasta, self.mapping, "MAG_1"), {"c1": 4, "c2": 6}
        )

    def test_falls_back_to_parsing_the_fasta(self):
        fasta = self._fasta("noindex.fa", [("c1", "ACGT"), ("c2", "ACGTAC")])
        self.assertEqual(
            contig_lengths_for_mag(fasta, self.mapping, "MAG_1"), {"c1": 4, "c2": 6}
        )

    def test_raises_when_mag_absent_from_mapping(self):
        fasta = self._fasta("r2.fa", [("c1", "ACGT")])
        with self.assertRaises(ValueError):
            contig_lengths_for_mag(fasta, self.mapping, "MAG_MISSING")

    def test_raises_when_mapping_contig_missing_from_reference(self):
        """Mapping names c2 but the reference lacks it: the two inputs disagree."""
        fasta = self._fasta("r3.fa", [("c1", "ACGT")])
        with self.assertRaises(ValueError):
            contig_lengths_for_mag(fasta, self.mapping, "MAG_1")


if __name__ == "__main__":
    unittest.main()


class TestDenseCountsFromCache(unittest.TestCase):
    """Pairwise ANI's per-sample dense arrays come from the cache, one timepoint file at a time."""

    def setUp(self):
        from tests.utilities.test_cache_io import write_cache_in_batches

        self.tmp = tempfile.mkdtemp()
        self.cache = os.path.join(self.tmp, "allele_freq_cache")
        # two samples at 'pre' in two row groups (as the writer produces), a 6 bp contig
        write_cache_in_batches(
            self.cache,
            "MAG_1",
            "pre",
            [
                [("c1", p, "S1", "treatment", "m1", 20, 0, 0, 0) for p in range(6)],
                [("c1", p, "S2", "control", "m2", 0, 0, 20, 0) for p in range(5)]
                + [("c1", 5, "S2", "control", "m2", 2, 0, 0, 0)],
            ],
        )

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_each_sample_gets_dense_arrays_and_stats(self):
        got = dense_counts_from_cache(
            self.cache, "MAG_1", "pre", ["S1", "S2"], {"c1": 6}, min_cov=5
        )
        self.assertEqual(sorted(got), ["S1", "S2"])
        dense, stats = got["S2"]
        self.assertEqual(dense["c1"][0].tolist(), [0, 0, 20, 0])
        self.assertEqual(dense["c1"][5].tolist(), [2, 0, 0, 0])
        self.assertEqual(stats, {"n_positions": 6, "breadth_minCov": 5 / 6})

    def test_only_requested_samples_are_loaded(self):
        got = dense_counts_from_cache(
            self.cache, "MAG_1", "pre", ["S1"], {"c1": 6}, min_cov=5
        )
        self.assertEqual(list(got), ["S1"])

    def test_requested_sample_absent_from_the_cache_raises(self):
        with self.assertRaises(ValueError) as ctx:
            dense_counts_from_cache(
                self.cache, "MAG_1", "pre", ["S1", "S9"], {"c1": 6}, min_cov=5
            )
        self.assertIn("S9", str(ctx.exception))

    def test_contig_outside_the_mapping_raises(self):
        with self.assertRaises(ValueError):
            dense_counts_from_cache(
                self.cache, "MAG_1", "pre", ["S1"], {"zzz": 6}, min_cov=5
            )
