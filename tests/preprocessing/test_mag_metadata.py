#!/usr/bin/env python3
"""Tests for ``alleleflux-metadata`` (mag_metadata.py): profiles + run metadata
-> one ``{mag}_metadata.tsv`` per MAG.

The script is the FIRST consumer of the user's metadata file, so this is where
the metadata contract is pinned:

* every identifier column is text exactly as written (``007`` stays ``007``,
  a replicate column with one blank cell does not turn ``20`` into ``20.0``);
* required columns, duplicate sample ids and profiles missing from the
  metadata stop the job;
* the longitudinal filter keeps only the requested timepoints and groups;
* the per-MAG files carry a fixed column order and a stable sort.
"""

import os
import shutil
import subprocess
import tempfile
import unittest

import pandas as pd

from alleleflux.scripts.preprocessing.mag_metadata import (
    find_mag_files,
    merge_metadata,
    write_output_files,
)


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #
def _profiles(root, samples, mags=("MAG_A",)):
    """``{root}/{sample}/{sample}_{mag}_profiled.tsv.gz`` for every sample x MAG (empty files)."""
    for s in samples:
        os.makedirs(os.path.join(root, s), exist_ok=True)
        for m in mags:
            open(os.path.join(root, s, f"{s}_{m}_profiled.tsv.gz"), "w").close()


def _sheet(rows, path):
    """Write ``rows`` (dicts) as the run metadata TSV and return the path."""
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def _row(sample, subject, group, time, replicate=None, **extra):
    row = {
        "sample_id": sample,
        "bam_path": f"/bams/{sample}.bam",
        "subjectID": subject,
        "group": group,
        "time": time,
    }
    if replicate is not None:
        row["replicate"] = replicate
    row.update(extra)
    return row


class TestMergeMetadataTextColumns(unittest.TestCase):
    """Identifier columns are text exactly as the metadata spells them."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.root = os.path.join(self.tmp, "profiles")
        self.meta = os.path.join(self.tmp, "metadata.tsv")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _merged(self, rows, timepoints=("pre", "end")):
        samples = [r["sample_id"] for r in rows]
        groups = sorted(
            {str(r["group"]).strip() for r in rows}
        )  # as the rule passes them
        _profiles(self.root, samples)
        _sheet(rows, self.meta)
        merged = merge_metadata(
            find_mag_files(self.root),
            self.meta,
            list(timepoints),
            groups,
            "longitudinal",
        )
        return {s["sample_id"]: s for s in merged["MAG_A"]}

    def test_numeric_group_labels_stay_text(self):
        """Groups 20 and 40 (ages) sniff as int64; they must come out as '20' / '40'."""
        got = self._merged([_row("s1", "m1", 20, "pre"), _row("s2", "m1", 40, "end")])
        self.assertEqual(got["s1"]["group"], "20")
        self.assertEqual(got["s2"]["group"], "40")

    def test_leading_zeros_survive(self):
        """subjectID 007 and sample_id 0123 are text, not the integers 7 and 123."""
        got = self._merged(
            [_row("0123", "007", "20", "pre"), _row("0124", "007", "20", "end")]
        )
        self.assertEqual(set(got), {"0123", "0124"})
        self.assertEqual(got["0123"]["subjectID"], "007")

    def test_blank_identifier_cell_is_refused(self):
        """One blank replicate cell used to make pandas read the column as float,
        so '12' became '12.0' in every other row and matched no cage label
        downstream.  A blank identifier is a malformed file: stop, name the
        column and the sample, never guess."""
        with self.assertRaisesRegex(ValueError, r"replicate.*s2"):
            self._merged(
                [
                    _row("s1", "m1", "20", "pre", replicate=12),
                    _row("s2", "m2", "20", "pre", replicate=None),
                    _row("s3", "m1", "20", "end", replicate=12),
                ]
            )

    def test_na_words_count_as_blank(self):
        """'NA', 'null', 'None' in an identifier column are missing values, not
        labels (pandas convention): refused like a blank, with the sample named."""
        with self.assertRaisesRegex(ValueError, r"group.*s2"):
            self._merged([_row("s1", "m1", "20", "pre"), _row("s2", "m1", "NA", "end")])

    def test_integer_replicates_stay_integers_as_text(self):
        got = self._merged(
            [
                _row("s1", "m1", "20", "pre", replicate=12),
                _row("s2", "m1", "20", "end", replicate=12),
            ]
        )
        self.assertEqual(got["s1"]["replicate"], "12")

    def test_replicate_defaults_to_subject(self):
        got = self._merged(
            [_row("s1", "m1", "20", "pre"), _row("s2", "m1", "20", "end")]
        )
        self.assertEqual(got["s1"]["replicate"], "m1")

    def test_whitespace_is_stripped(self):
        got = self._merged(
            [_row("s1", " m1 ", " 20", "pre "), _row("s2", "m1", "20", "end")]
        )
        self.assertEqual(
            (got["s1"]["subjectID"], got["s1"]["group"], got["s1"]["time"]),
            ("m1", "20", "pre"),
        )


class TestMergeMetadataValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.root = os.path.join(self.tmp, "profiles")
        self.meta = os.path.join(self.tmp, "metadata.tsv")
        _profiles(self.root, ["s1", "s2"])

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _merge(self, rows, **kw):
        _sheet(rows, self.meta)
        args = dict(timepoints=["pre", "end"], groups=["A"], data_type="longitudinal")
        args.update(kw)
        return merge_metadata(find_mag_files(self.root), self.meta, **args)

    def test_missing_required_column_raises(self):
        rows = [_row("s1", "m1", "A", "pre"), _row("s2", "m1", "A", "end")]
        for r in rows:
            del r["subjectID"]
        with self.assertRaisesRegex(ValueError, "subjectID"):
            self._merge(rows)

    def test_duplicate_sample_id_raises(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            self._merge([_row("s1", "m1", "A", "pre"), _row("s1", "m1", "A", "end")])

    def test_profile_without_metadata_row_raises(self):
        # s2 has a profile on disk but no metadata row
        with self.assertRaisesRegex(ValueError, "not found"):
            self._merge([_row("s1", "m1", "A", "pre")])

    def test_unknown_group_or_timepoint_raises(self):
        rows = [_row("s1", "m1", "A", "pre"), _row("s2", "m1", "A", "end")]
        with self.assertRaisesRegex(ValueError, "groups"):
            self._merge(rows, groups=["B"])
        with self.assertRaisesRegex(ValueError, "timepoints"):
            self._merge(rows, timepoints=["pre", "mid"])

    def test_longitudinal_filter_keeps_requested_timepoints_and_groups_only(self):
        _profiles(self.root, ["s3", "s4"])
        rows = [
            _row("s1", "m1", "A", "pre"),
            _row("s2", "m1", "A", "end"),
            _row("s3", "m2", "B", "pre"),  # other group
            _row("s4", "m1", "A", "mid"),  # other timepoint
        ]
        merged = self._merge(rows)
        self.assertEqual(sorted(s["sample_id"] for s in merged["MAG_A"]), ["s1", "s2"])

    def test_single_with_two_timepoints_and_no_choice_raises(self):
        rows = [_row("s1", "m1", "A", "pre"), _row("s2", "m2", "A", "end")]
        with self.assertRaisesRegex(ValueError, "multiple unique time"):
            self._merge(rows, timepoints=None, data_type="single")


class TestWriteOutputFiles(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_per_mag_file_columns_and_sort(self):
        mag_dict = {
            "MAG_A": [
                dict(
                    sample_id="s2",
                    subjectID="m1",
                    group="A",
                    replicate="r1",
                    time="end",
                    file_path="/p/s2",
                    bam_path="/b/s2.bam",
                    extra="dropped",
                ),
                dict(
                    sample_id="s1",
                    subjectID="m1",
                    group="A",
                    replicate="r1",
                    time="pre",
                    file_path="/p/s1",
                    bam_path="/b/s1.bam",
                    extra="dropped",
                ),
            ]
        }
        write_output_files(mag_dict, self.tmp, "longitudinal")
        out = pd.read_csv(
            os.path.join(self.tmp, "MAG_A_metadata.tsv"), sep="\t", dtype=str
        )
        self.assertEqual(
            list(out.columns),
            [
                "sample_id",
                "subjectID",
                "time",
                "group",
                "replicate",
                "file_path",
                "bam_path",
            ],
        )
        self.assertEqual(list(out["sample_id"]), ["s1", "s2"])

    def test_single_has_no_time_column(self):
        mag_dict = {
            "MAG_A": [
                dict(
                    sample_id="s1",
                    subjectID="m1",
                    group="A",
                    replicate="r1",
                    file_path="/p/s1",
                    bam_path="/b/s1.bam",
                )
            ]
        }
        write_output_files(mag_dict, self.tmp, "single")
        out = pd.read_csv(
            os.path.join(self.tmp, "MAG_A_metadata.tsv"), sep="\t", dtype=str
        )
        self.assertNotIn("time", out.columns)


class TestMagMetadataCLI(unittest.TestCase):
    """End to end through the entry point, as the Snakemake rule calls it."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.root = os.path.join(self.tmp, "profiles")
        self.out = os.path.join(self.tmp, "inputMetadata")
        self.meta = os.path.join(self.tmp, "metadata.tsv")
        _profiles(self.root, ["s1", "s2", "s3", "s4"], mags=("MAG_A", "MAG_B"))
        _sheet(
            [
                _row("s1", "007", 20, "pre", replicate=1),
                _row("s2", "007", 20, "end", replicate=1),
                _row("s3", "008", 40, "pre", replicate=2),
                _row("s4", "008", 40, "end", replicate=2),
            ],
            self.meta,
        )

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_writes_one_file_per_mag_and_an_index_with_text_labels(self):
        done = subprocess.run(
            [
                "alleleflux-metadata",
                "--rootDir",
                self.root,
                "--metadata",
                self.meta,
                "--outDir",
                self.out,
                "--groups",
                "20",
                "40",
                "--data_type",
                "longitudinal",
                "--timepoints",
                "pre",
                "end",
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual(
            sorted(os.listdir(self.out)),
            ["MAG_A_metadata.tsv", "MAG_B_metadata.tsv", "summary_mag_index.tsv"],
        )
        a = pd.read_csv(
            os.path.join(self.out, "MAG_A_metadata.tsv"), sep="\t", dtype=str
        )
        self.assertEqual(list(a["sample_id"]), ["s1", "s2", "s3", "s4"])
        self.assertEqual(list(a["subjectID"]), ["007", "007", "008", "008"])
        self.assertEqual(list(a["group"]), ["20", "20", "40", "40"])
        self.assertEqual(list(a["replicate"]), ["1", "1", "2", "2"])  # text, never 1.0
        idx = pd.read_csv(
            os.path.join(self.out, "summary_mag_index.tsv"), sep="\t", dtype=str
        )
        self.assertEqual(set(idx["mag_id"]), {"MAG_A", "MAG_B"})
        self.assertEqual(len(idx), 8)

    def _run(self):
        return subprocess.run(
            [
                "alleleflux-metadata",
                "--rootDir",
                self.root,
                "--metadata",
                self.meta,
                "--outDir",
                self.out,
                "--groups",
                "20",
                "40",
                "--data_type",
                "longitudinal",
                "--timepoints",
                "pre",
                "end",
            ],
            capture_output=True,
            text=True,
        )

    def test_bad_collection_date_stops_the_run_at_the_metadata_step(self):
        """The metadata step is the first consumer of the sheet and runs in seconds:
        a date that is not YYYY-MM-DD must stop the run HERE, naming the sample,
        not hours later in baseline presence.  Nothing is written."""
        meta = pd.read_csv(self.meta, sep="\t", dtype=str)
        meta["collection_date"] = "2017-01-10"
        meta.loc[meta.sample_id == "s3", "collection_date"] = "10/01/2017"
        meta.to_csv(self.meta, sep="\t", index=False)
        done = self._run()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("YYYY-MM-DD", done.stderr)
        self.assertIn("s3", done.stderr)
        self.assertFalse(os.path.isdir(self.out))

    def test_good_and_blank_collection_dates_pass_and_stay_out_of_per_mag_files(self):
        """ISO dates and blanks are fine; the column is for baseline presence (which
        reads the run sheet itself) and is NOT copied into the per-MAG files."""
        meta = pd.read_csv(self.meta, sep="\t", dtype=str)
        meta["collection_date"] = ["2017-01-10", "2017-06-10", "", "2017-06-10"]
        meta.to_csv(self.meta, sep="\t", index=False)
        done = self._run()
        self.assertEqual(done.returncode, 0, done.stderr)
        a = pd.read_csv(
            os.path.join(self.out, "MAG_A_metadata.tsv"), sep="\t", dtype=str
        )
        self.assertNotIn("collection_date", a.columns)
        self.assertEqual(list(a["sample_id"]), ["s1", "s2", "s3", "s4"])


if __name__ == "__main__":
    unittest.main()
