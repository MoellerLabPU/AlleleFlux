"""The allele-frequency cache as the single per-site source after QC."""

import os

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from alleleflux.scripts.utilities.cache_io import (
    cache_path,
    cache_samples,
    load_cache_sites,
)

COLS = [
    "contig",
    "position",
    "gene_id",
    "ref_base",
    "total_coverage",
    "A",
    "C",
    "G",
    "T",
    "sample_id",
    "group",
    "subjectID",
    "replicate",
    "A_frequency",
    "T_frequency",
    "G_frequency",
    "C_frequency",
    "time",
]


def _frame(rows, timepoint):
    """rows: (contig, position, sample, group, subject, A, C, G, T); frequencies derived,
    ref_base fixed at "A" (the cache stores the reference base, dictionary-encoded)."""
    recs = []
    for contig, pos, sample, group, subject, a, c, g, t in rows:
        tot = a + c + g + t
        recs.append(
            dict(
                contig=contig,
                position=pos,
                ref_base="A",
                gene_id="g1",
                total_coverage=tot,
                A=a,
                C=c,
                G=g,
                T=t,
                sample_id=sample,
                group=group,
                subjectID=subject,
                replicate=subject,
                A_frequency=a / tot,
                T_frequency=t / tot,
                G_frequency=g / tot,
                C_frequency=c / tot,
                time=timepoint,
            )
        )
    df = pd.DataFrame(recs, columns=COLS)
    for col in ("sample_id", "group", "subjectID", "replicate", "ref_base"):
        df[col] = df[col].astype("category")  # as the cache writer stores them
    return df


def write_cache(cache_dir, mag, timepoint, rows):
    """One cache file written like the pipeline does (one row group per call)."""
    path = cache_path(cache_dir, mag, timepoint)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    _frame(rows, timepoint).to_parquet(path, index=False)
    return path


def write_cache_in_batches(cache_dir, mag, timepoint, batches):
    """Several row groups with their own dictionaries, as the streaming writer produces."""
    path = cache_path(cache_dir, mag, timepoint)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    writer = None
    for rows in batches:
        table = pa.Table.from_pandas(_frame(rows, timepoint), preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(path, table.schema, use_dictionary=True)
        writer.write_table(table)
    writer.close()
    return path


@pytest.fixture
def cache(tmp_path):
    d = str(tmp_path / "allele_freq_cache")
    write_cache(
        d,
        "MAG_A",
        "pre",
        [
            ("c1", 1, "s1_pre", "20", "m1", 25, 0, 5, 0),
            ("c1", 2, "s1_pre", "20", "m1", 30, 0, 0, 0),
            ("c1", 1, "s2_pre", "20", "m2", 30, 0, 0, 0),
        ],
    )
    write_cache(d, "MAG_A", "end", [("c1", 1, "s1_end", "20", "m1", 6, 0, 24, 0)])
    return d


def test_cache_path_layout():
    assert (
        cache_path("/x", "MAG_A", "pre") == "/x/pre/MAG_A_pre_allele_frequency.parquet"
    )


def test_sites_and_timepoints_are_filtered(cache):
    sites = pd.DataFrame({"contig": ["c1"], "position": [1]})
    got = load_cache_sites(cache, "MAG_A", ["pre", "end"], sites=sites)
    assert sorted(zip(got["sample_id"], got["timepoint"])) == [
        ("s1_end", "end"),
        ("s1_pre", "pre"),
        ("s2_pre", "pre"),
    ]
    assert set(got["position"]) == {1}
    assert (
        got["sample_id"].dtype == object and got["group"].tolist()[0] == "20"
    )  # labels as str


def test_sample_filter_and_column_projection(cache):
    got = load_cache_sites(
        cache, "MAG_A", ["pre"], sample_ids=["s2_pre"], columns=["A", "total_coverage"]
    )
    assert got["sample_id"].tolist() == ["s2_pre"]
    assert set(got.columns) == {
        "contig",
        "position",
        "sample_id",
        "timepoint",
        "A",
        "total_coverage",
    }


def test_no_sites_means_every_position(cache):
    got = load_cache_sites(cache, "MAG_A", ["pre"])
    assert len(got) == 3


def test_missing_cache_file_raises(cache):
    with pytest.raises(FileNotFoundError, match="MAG_B_pre_allele_frequency.parquet"):
        load_cache_sites(cache, "MAG_B", ["pre"])


def test_nothing_matches_gives_an_empty_frame_with_keys(cache):
    got = load_cache_sites(
        cache, "MAG_A", ["pre"], sites=pd.DataFrame({"contig": ["c9"], "position": [1]})
    )
    assert len(got) == 0 and {"contig", "position", "sample_id", "timepoint"} <= set(
        got.columns
    )


def test_cache_samples_is_the_qc_passing_roster(cache):
    assert cache_samples(cache, "MAG_A", "pre") == ["s1_pre", "s2_pre"]


def test_batches_are_filtered_independently(tmp_path):
    """Row groups carry their own dictionaries; both samples must come back as plain strings."""
    d = str(tmp_path / "allele_freq_cache")
    write_cache_in_batches(
        d,
        "MAG_A",
        "pre",
        [
            [("c1", 1, "s1_pre", "20", "m1", 10, 0, 0, 0)],
            [("c1", 1, "s2_pre", "40", "m2", 0, 0, 10, 0)],
        ],
    )
    got = load_cache_sites(
        d, "MAG_A", ["pre"], sites=pd.DataFrame({"contig": ["c1"], "position": [1]})
    )
    assert sorted(got["sample_id"]) == ["s1_pre", "s2_pre"] and sorted(
        got["group"]
    ) == ["20", "40"]
    assert cache_samples(d, "MAG_A", "pre") == ["s1_pre", "s2_pre"]
