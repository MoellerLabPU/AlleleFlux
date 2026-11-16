"""``alleleflux-cache-allele-freq --allow_empty``: a genome with no usable sample still gets a cache.

Under the cohort baseline scope the workflow asks for a cache of every tested
genome at every witness timepoint.  A genome can have no QC-passing sample at
one of them (its QC table there is header-only, or every sample failed).  The
workflow's list of cache files is fixed before any QC table is read, so that
job must still write its two outputs: a cache with the normal columns and zero
rows, and a header-only sidecar.  The readers then see "no sample here".

Without the flag nothing changes: no QC-passing sample is an error.  With the
flag a REAL failure is still an error: samples that passed QC but whose
profile cannot be read never produce an empty cache.

Each test drives ``main()`` on small files written in ``tmp_path``.
"""

import gzip

import pandas as pd
import pyarrow.parquet as pq
import pytest

from alleleflux.scripts.analysis.allele_frequency import _allele_freq_common as common
from alleleflux.scripts.analysis.allele_frequency import allele_freq_cache as cli

PROFILE_HEADER = "contig\tposition\tref_base\ttotal_coverage\tA\tC\tG\tT\tgene_id\n"
# The QC columns the cache writer reads (a real QC table has more).
QC_COLUMNS = [
    "sample_id",
    "MAG_ID",
    "file_path",
    "group",
    "subjectID",
    "replicate",
    "time",
    "breadth",
    "genome_size",
    "average_coverage",
    "average_coverage_including_zeros",
    "coverage_threshold_passed",
]
SIDECAR_COLUMNS = [
    "sample_id",
    "mag_id",
    "timepoint",
    "mean_coverage",
    "mean_coverage_after_filter",
    "n_sites_covered",
    "n_removed_low",
    "n_removed_high",
    "n_kept",
]


def _write_profile(path):
    """Four positions at 10 reads each, all base A."""
    with gzip.open(path, "wt") as fh:
        fh.write(PROFILE_HEADER)
        for pos in range(1, 5):
            fh.write(f"c1\t{pos}\tA\t10\t10\t0\t0\t0\tgene_1\n")
    return str(path)


def _write_qc(path, rows, with_time=True):
    """rows: (sample_id, profile_path, passed).  ``with_time`` False = a single-data table."""
    columns = [c for c in QC_COLUMNS if with_time or c != "time"]
    records = []
    for sample_id, profile, passed in rows:
        record = {
            "sample_id": sample_id,
            "MAG_ID": "MAG",
            "file_path": profile,
            "group": "a",
            "subjectID": f"m_{sample_id}",
            "replicate": "r1",
            "time": "t0",
            "breadth": 0.5,
            "genome_size": 7,
            "average_coverage": 10.0,
            "average_coverage_including_zeros": 5.0,
            "coverage_threshold_passed": passed,
        }
        records.append({c: record[c] for c in columns})
    pd.DataFrame(records, columns=columns).to_csv(path, sep="\t", index=False)
    return str(path)


def _run(monkeypatch, qc, out, *extra, data_type="longitudinal"):
    """Call the command's ``main()`` as the workflow would."""
    argv = [
        "alleleflux-cache-allele-freq",
        "--magID",
        "MAG",
        "--qc_files",
        qc,
        "--data_type",
        data_type,
        "--output_path",
        str(out),
        "--cpus",
        "1",
        *extra,
    ]
    if data_type == "longitudinal":
        argv += ["--timepoint", "t0"]
    monkeypatch.setattr("sys.argv", argv)
    cli.main()


def _sidecar_of(out):
    return out.with_name(out.name.replace(".parquet", "_site_depth_filter.tsv"))


class TestAllowEmpty:
    def test_allow_empty_writes_a_zero_row_cache_with_the_longitudinal_schema(
        self, tmp_path, monkeypatch
    ):
        """Every sample failed QC -> zero rows, and exactly the columns and
        types of a cache that has rows (compared against one written by the
        same command from a passing sample)."""
        profile = _write_profile(tmp_path / "s1_MAG_profiled.tsv.gz")
        failed = _write_qc(tmp_path / "failed_QC.tsv", [("s1", profile, False)])
        passed = _write_qc(tmp_path / "passed_QC.tsv", [("s1", profile, True)])
        empty_out = tmp_path / "empty" / "MAG_t0_allele_frequency.parquet"
        real_out = tmp_path / "real" / "MAG_t0_allele_frequency.parquet"

        _run(monkeypatch, failed, empty_out, "--allow_empty")
        _run(monkeypatch, passed, real_out)

        empty_schema = pq.read_schema(empty_out)
        assert pq.read_metadata(empty_out).num_rows == 0
        assert pq.read_metadata(real_out).num_rows == 4
        assert empty_schema.equals(pq.read_schema(real_out), check_metadata=False)
        assert empty_schema.equals(
            common.empty_cache_table("longitudinal").schema, check_metadata=False
        )
        assert empty_schema.names == common.CACHE_COLUMNS_LONGITUDINAL
        # No half-written file is left behind.
        assert not empty_out.with_suffix(".parquet.tmp").exists()

    def test_allow_empty_accepts_a_header_only_qc_table(self, tmp_path, monkeypatch):
        """The genome has no sample at all at this timepoint: QC wrote a header-only table."""
        qc = _write_qc(tmp_path / "MAG_QC.tsv", [])
        out = tmp_path / "MAG_t0_allele_frequency.parquet"
        _run(monkeypatch, qc, out, "--allow_empty")
        assert pq.read_metadata(out).num_rows == 0
        assert pq.read_schema(out).names == common.CACHE_COLUMNS_LONGITUDINAL

    def test_allow_empty_writes_a_header_only_sidecar(self, tmp_path, monkeypatch):
        qc = _write_qc(tmp_path / "MAG_QC.tsv", [])
        out = tmp_path / "MAG_t0_allele_frequency.parquet"
        _run(monkeypatch, qc, out, "--allow_empty")
        side = pd.read_csv(_sidecar_of(out), sep="\t")
        assert side.columns.tolist() == SIDECAR_COLUMNS
        assert len(side) == 0

    def test_allow_empty_on_single_data_has_no_time_column(self, tmp_path, monkeypatch):
        qc = _write_qc(tmp_path / "MAG_QC.tsv", [], with_time=False)
        out = tmp_path / "MAG_t0_allele_frequency.parquet"
        _run(monkeypatch, qc, out, "--allow_empty", data_type="single")
        assert pq.read_metadata(out).num_rows == 0
        assert pq.read_schema(out).names == common.CACHE_COLUMNS_BASE

    def test_a_cache_with_samples_is_the_same_with_and_without_the_flag(
        self, tmp_path, monkeypatch
    ):
        profile = _write_profile(tmp_path / "s1_MAG_profiled.tsv.gz")
        qc = _write_qc(tmp_path / "MAG_QC.tsv", [("s1", profile, True)])
        plain = tmp_path / "plain" / "MAG_t0_allele_frequency.parquet"
        flagged = tmp_path / "flagged" / "MAG_t0_allele_frequency.parquet"
        _run(monkeypatch, qc, plain)
        _run(monkeypatch, qc, flagged, "--allow_empty")
        pd.testing.assert_frame_equal(
            pd.read_parquet(plain), pd.read_parquet(flagged), check_exact=True
        )
        pd.testing.assert_frame_equal(
            pd.read_csv(_sidecar_of(plain), sep="\t"),
            pd.read_csv(_sidecar_of(flagged), sep="\t"),
        )


class TestStillAnError:
    @pytest.mark.parametrize("rows", [[], [("s1", "unused", False)]])
    def test_without_allow_empty_no_passing_sample_is_an_error(
        self, tmp_path, monkeypatch, rows
    ):
        """Unchanged behaviour: the command stops and writes neither output."""
        qc = _write_qc(tmp_path / "MAG_QC.tsv", rows)
        out = tmp_path / "MAG_t0_allele_frequency.parquet"
        with pytest.raises(ValueError, match="MAG.*t0"):
            _run(monkeypatch, qc, out)
        assert not out.exists() and not _sidecar_of(out).exists()

    def test_allow_empty_does_not_hide_an_unreadable_profile(
        self, tmp_path, monkeypatch
    ):
        """A sample passed QC but its profile is not on disk: a real failure,
        never an empty cache."""
        qc = _write_qc(
            tmp_path / "MAG_QC.tsv",
            [("s1", str(tmp_path / "not_there_profiled.tsv.gz"), True)],
        )
        out = tmp_path / "MAG_t0_allele_frequency.parquet"
        with pytest.raises(FileNotFoundError):
            _run(monkeypatch, qc, out, "--allow_empty")
        assert not out.exists()

    def test_allow_empty_does_not_hide_a_qc_table_without_the_verdict_column(
        self, tmp_path, monkeypatch
    ):
        qc = tmp_path / "MAG_QC.tsv"
        qc.write_text("sample_id\tfile_path\tgroup\ttime\n")
        out = tmp_path / "MAG_t0_allele_frequency.parquet"
        with pytest.raises(ValueError, match="coverage_threshold_passed"):
            _run(monkeypatch, str(qc), out, "--allow_empty")
        assert not out.exists()

    def test_allow_empty_does_not_hide_a_longitudinal_table_without_time(
        self, tmp_path, monkeypatch
    ):
        qc = _write_qc(tmp_path / "MAG_QC.tsv", [], with_time=False)
        out = tmp_path / "MAG_t0_allele_frequency.parquet"
        with pytest.raises(ValueError, match="time"):
            _run(monkeypatch, qc, out, "--allow_empty")
        assert not out.exists()


def test_empty_cache_table_refuses_an_unknown_data_type():
    with pytest.raises(ValueError, match="data_type"):
        common.empty_cache_table("paired")
