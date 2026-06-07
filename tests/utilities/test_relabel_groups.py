#!/usr/bin/env python3
"""Unit tests for ``relabel_groups_from_metadata`` (Phase 1 of the null-run feature).

This helper is the in-memory seam that lets a permuted run reuse the real run's
cache/QC: it re-derives permuted ``group``/``subjectID`` labels by joining the
loaded frame to a permuted metadata sheet on ``sample_id``.  These tests pin the
contract every cache consumer relies on — correct join, idempotency, identity
no-op, categorical preservation, graceful handling of a missing key, and the
numeric-group → str coercion guard.
"""

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from alleleflux.scripts.utilities.utilities import relabel_groups_from_metadata


def _cache_like():
    """Frame shaped like a loaded cache slice: per-sample rows + a value column."""
    return pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2", "s2", "s3", "s3"],
            "group": ["A", "A", "A", "A", "B", "B"],
            "subjectID": ["s1", "s1", "s2", "s2", "s3", "s3"],
            "A_frequency": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        }
    )


def _permuted_meta():
    """Permuted sheet: s2 moves A->B, s3 moves B->A (per-sample labels)."""
    return pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "s3"],
            "group": ["A", "B", "A"],
            "subjectID": ["s1", "s2", "s3"],
        }
    )


class TestRelabelGroups(unittest.TestCase):
    def test_relabels_group_by_sample_id(self):
        out = relabel_groups_from_metadata(_cache_like(), _permuted_meta())
        # s1->A (unchanged), s2->B, s3->A
        self.assertEqual(list(out["group"]), ["A", "A", "B", "B", "A", "A"])

    def test_frequency_column_untouched(self):
        df = _cache_like()
        before = df["A_frequency"].copy()
        out = relabel_groups_from_metadata(df, _permuted_meta())
        pd.testing.assert_series_equal(before, out["A_frequency"])

    def test_subjectid_never_touched(self):
        # Permutation moves ONLY the group label; subjectID must stay tied to its
        # original rows even if the permuted sheet carries different subjectIDs
        # (rewriting it would risk re-coercing the dtype of the merge key).
        meta = _permuted_meta()
        meta["subjectID"] = ["s1", "X2", "X3"]  # deliberately different
        out = relabel_groups_from_metadata(_cache_like(), meta)
        self.assertEqual(list(out["subjectID"]), ["s1", "s1", "s2", "s2", "s3", "s3"])
        # group still relabeled as expected
        self.assertEqual(list(out["group"]), ["A", "A", "B", "B", "A", "A"])

    def test_idempotent(self):
        once = relabel_groups_from_metadata(_cache_like(), _permuted_meta())
        twice = relabel_groups_from_metadata(once.copy(), _permuted_meta())
        pd.testing.assert_frame_equal(once, twice)

    def test_identity_metadata_is_noop(self):
        df = _cache_like()
        identity = df[["sample_id", "group", "subjectID"]].drop_duplicates()
        out = relabel_groups_from_metadata(df.copy(), identity)
        self.assertEqual(list(out["group"]), ["A", "A", "A", "A", "B", "B"])

    def test_preserves_categorical_dtype(self):
        df = _cache_like()
        df["group"] = df["group"].astype("category")
        out = relabel_groups_from_metadata(df, _permuted_meta())
        self.assertIsInstance(out["group"].dtype, pd.CategoricalDtype)
        self.assertEqual(list(out["group"]), ["A", "A", "B", "B", "A", "A"])

    def test_missing_key_in_df_returns_unchanged(self):
        # An already-aggregated frame that dropped sample_id is returned as-is.
        df = _cache_like().drop(columns=["sample_id"])
        out = relabel_groups_from_metadata(df.copy(), _permuted_meta())
        pd.testing.assert_frame_equal(out, df)

    def test_sample_absent_from_metadata_raises(self):
        meta = _permuted_meta()
        meta = meta[meta["sample_id"] != "s3"]  # drop s3 from the map
        with self.assertRaises(ValueError) as ctx:
            relabel_groups_from_metadata(_cache_like(), meta)
        self.assertIn("s3", str(ctx.exception))

    def test_numeric_groups_coerced_to_str(self):
        # DRIDO-style numeric group names must relabel as strings, not ints.
        df = pd.DataFrame(
            {"sample_id": ["s1", "s2"], "group": ["20", "40"], "v": [1.0, 2.0]}
        )
        meta = pd.DataFrame({"sample_id": ["s1", "s2"], "group": [40, 20]})  # int
        out = relabel_groups_from_metadata(df, meta)
        self.assertEqual(list(out["group"]), ["40", "20"])
        self.assertTrue(all(isinstance(v, str) for v in out["group"]))

    def test_missing_key_in_metadata_raises(self):
        bad_meta = _permuted_meta().rename(columns={"sample_id": "id"})
        with self.assertRaises(ValueError):
            relabel_groups_from_metadata(_cache_like(), bad_meta)

    def test_accepts_tsv_path(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "perm.tsv"
            _permuted_meta().to_csv(p, sep="\t", index=False)
            out = relabel_groups_from_metadata(_cache_like(), str(p))
        self.assertEqual(list(out["group"]), ["A", "A", "B", "B", "A", "A"])


if __name__ == "__main__":
    unittest.main()
