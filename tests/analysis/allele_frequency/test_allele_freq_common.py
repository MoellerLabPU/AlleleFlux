"""Tests for the cache worker's profile read (process_mag_files).

The profile's ``total_coverage`` is A+C+G+T by construction (no N column since
2026-10-15); the cache carries it through unchanged and the four frequencies
sum to 1 at every position.
"""

import gzip
import os
import tempfile
import unittest

import numpy as np

from alleleflux.scripts.analysis.allele_frequency._allele_freq_common import (
    init_worker,
    process_mag_files,
)

HEADER = "contig\tposition\tref_base\ttotal_coverage\tA\tC\tG\tT\tgene_id\n"
# (contig, position, ref, A, C, G, T)
ROWS = [
    ("c1", 0, "A", 8, 0, 2, 0),
    ("c1", 1, "C", 0, 4, 0, 0),
    ("c2", 0, "T", 0, 0, 0, 60),
]


class TestProcessMagFilesCoverage(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        init_worker(
            {
                "s1": {
                    "group": "g",
                    "subjectID": "m1",
                    "replicate": "r1",
                    "average_coverage": 30.0,
                }
            },
            "single",
        )
        self.path = os.path.join(self.tmp, "s1_MAG1_profiled.tsv.gz")
        with gzip.open(self.path, "wt") as fh:
            fh.write(HEADER)
            for contig, pos, ref, a, c, g, t in ROWS:
                fh.write(
                    "\t".join(
                        map(str, [contig, pos, ref, a + c + g + t, a, c, g, t, "g1"])
                    )
                    + "\n"
                )

    def test_total_coverage_and_frequencies(self):
        df = process_mag_files(("s1", self.path, "MAG1"))
        df = df.sort_values(["contig", "position"]).reset_index(drop=True)
        self.assertEqual(df["total_coverage"].tolist(), [10, 4, 60])
        self.assertNotIn("N", df.columns)
        # ref_base is kept for the cache (dN/dS reads it), as a category so
        # Parquet dictionary-encodes it: one letter per row costs ~1 byte.
        self.assertEqual(df["ref_base"].astype(str).tolist(), [r[2] for r in ROWS])
        self.assertEqual(df["ref_base"].dtype.name, "category")
        freq = df[["A_frequency", "C_frequency", "G_frequency", "T_frequency"]].sum(
            axis=1
        )
        np.testing.assert_allclose(freq.to_numpy(), 1.0, rtol=1e-6)


if __name__ == "__main__":
    unittest.main()
