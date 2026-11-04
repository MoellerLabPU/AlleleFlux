"""The cache worker applies the per-site depth filter before frequencies are computed.

The mean each sample is judged against is QC's ``average_coverage`` (covered-sites
mean), carried into the worker through the metadata dict; the QC table must carry
the ``average_coverage_including_zeros`` column as proof it was written by a version
whose ``average_coverage`` means that.
"""

import gzip

import pandas as pd
import pytest

from alleleflux.scripts.analysis.allele_frequency import _allele_freq_common as common
from alleleflux.scripts.utilities.site_depth_filter import SiteDepthFilter

PROFILE_HEADER = "contig\tposition\tref_base\ttotal_coverage\tA\tC\tG\tT\tgene_id\n"


def _write_profile(path, rows):
    """rows: (contig, position, A, C, G, T); total_coverage = A+C+G+T."""
    with gzip.open(path, "wt") as fh:
        fh.write(PROFILE_HEADER)
        for contig, pos, a, c, g, t in rows:
            fh.write(
                f"{contig}\t{pos}\tA\t{a + c + g + t}\t{a}\t{c}\t{g}\t{t}\tgene_1\n"
            )


@pytest.fixture
def profile(tmp_path):
    # four rows at 10x, a 40x position and an 80x position; QC's covered-sites
    # mean for this sample is (40 + 40 + 80) / 6 = 26.67, so 80 / 26.67 = 3.0
    # exactly is KEPT ('>' strict) -- use mean 24.0 from the QC row instead:
    # 80 / 24 = 3.33 dropped, 40 / 24 = 1.67 kept.
    rows = [("c1", i, 10, 0, 0, 0) for i in range(1, 5)]
    rows += [("c1", 5, 40, 0, 0, 0), ("c1", 6, 80, 0, 0, 0)]
    p = tmp_path / "s1_MAG_profiled.tsv.gz"
    _write_profile(p, rows)
    return str(p)


META = {
    "s1": {
        "group": "a",
        "subjectID": "m1",
        "replicate": "r1",
        "time": "t0",
        "average_coverage": 24.0,
    }
}


class TestProcessMagFilesSiteFilter:
    def test_filter_drops_the_pile_before_frequencies(self, profile):
        common.init_worker(META, "longitudinal", SiteDepthFilter())
        df = common.process_mag_files(("s1", profile, "MAG"))
        assert sorted(df["position"].tolist()) == [1, 2, 3, 4, 5]
        s = df.attrs["site_depth_filter"]
        assert s["sample_id"] == "s1" and s["mag_id"] == "MAG"
        assert s["mean_coverage"] == 24.0
        assert (s["n_removed_low"], s["n_removed_high"], s["n_kept"]) == (0, 1, 5)
        assert df.loc[df["position"] == 5, "A_frequency"].item() == pytest.approx(1.0)

    def test_disabled_filter_keeps_every_row(self, profile):
        common.init_worker(META, "longitudinal", SiteDepthFilter(enabled=False))
        df = common.process_mag_files(("s1", profile, "MAG"))
        assert len(df) == 6
        assert df.attrs["site_depth_filter"]["n_kept"] == 6

    def test_two_argument_init_worker_means_default_filter(self, profile):
        # existing callers pass two arguments; the house default applies
        common.init_worker(META, "longitudinal")
        df = common.process_mag_files(("s1", profile, "MAG"))
        assert len(df) == 5


class TestQcFormatMarker:
    def _qc(self, with_marker):
        cols = {
            "sample_id": ["s1"],
            "file_path": ["x"],
            "group": ["a"],
            "subjectID": ["m1"],
            "replicate": ["r1"],
            "breadth": [0.5],
            "genome_size": [7],
            "average_coverage": [24.0],
            "coverage_threshold_passed": [True],
        }
        if with_marker:
            cols["average_coverage_including_zeros"] = [0.1]
        return pd.DataFrame(cols)

    def test_old_qc_table_is_refused(self):
        with pytest.raises(ValueError, match="average_coverage_including_zeros"):
            common.build_metadata_from_qc(self._qc(with_marker=False))

    def test_new_qc_table_carries_the_mean(self):
        metadata, _ = common.build_metadata_from_qc(self._qc(with_marker=True))
        assert metadata["s1"]["average_coverage"] == 24.0


class TestCacheCliSidecar:
    def test_sidecar_written_next_to_cache(self, tmp_path, profile, monkeypatch):
        """Drive main() with a one-sample QC file; expect the cache and its sidecar."""
        from alleleflux.scripts.analysis.allele_frequency import (
            allele_freq_cache as cli,
        )

        qc = tmp_path / "MAG_QC.tsv"
        qc.write_text(
            "sample_id\tMAG_ID\tfile_path\tgroup\tsubjectID\treplicate\ttime\tbreadth\t"
            "genome_size\taverage_coverage\taverage_coverage_including_zeros\tcoverage_threshold_passed\n"
            f"s1\tMAG\t{profile}\ta\tm1\tr1\tt0\t0.5\t7\t24.0\t0.1\tTrue\n"
        )
        out = tmp_path / "MAG_t0_allele_frequency.parquet"
        monkeypatch.setattr(
            "sys.argv",
            [
                "alleleflux-cache-allele-freq",
                "--magID",
                "MAG",
                "--qc_files",
                str(qc),
                "--timepoint",
                "t0",
                "--data_type",
                "longitudinal",
                "--output_path",
                str(out),
                "--cpus",
                "1",
            ],
        )
        cli.main()
        side = pd.read_csv(
            tmp_path / "MAG_t0_allele_frequency_site_depth_filter.tsv", sep="\t"
        )
        assert side.columns.tolist() == [
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
        assert side.iloc[0]["n_kept"] == 5 and side.iloc[0]["timepoint"] == "t0"
        # The reference the sample was judged against is the QC table's
        # average_coverage (24.0 above), passed through untouched: QC row ->
        # build_metadata_from_qc -> init_worker -> apply_site_depth_filter.
        assert side.iloc[0]["sample_id"] == "s1"
        assert side.iloc[0]["mean_coverage"] == 24.0
        assert side.iloc[0]["mean_coverage_after_filter"] == pytest.approx(80 / 5)
        assert len(pd.read_parquet(out)) == 5
