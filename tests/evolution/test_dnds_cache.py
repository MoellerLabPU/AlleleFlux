"""dN/dS reads its two samples from the allele-frequency cache, not from raw profiles."""

import pandas as pd
import pytest

from alleleflux.scripts.evolution.dnds_from_timepoints import load_sample_counts
from tests.utilities.test_cache_io import write_cache


@pytest.fixture
def cache(tmp_path):
    d = str(tmp_path / "allele_freq_cache")
    write_cache(
        d,
        "MAG_A",
        "pre",
        [
            ("c1", 0, "s1_pre", "g", "m1", 10, 0, 0, 0),
            ("c1", 1, "s1_pre", "g", "m1", 0, 8, 0, 0),
            ("c1", 0, "s2_pre", "g", "m2", 5, 0, 0, 0),
        ],
    )
    return d


def test_one_sample_one_timepoint_with_the_reconstruction_columns(cache):
    df = load_sample_counts(cache, "MAG_A", "pre", "s1_pre")
    # the reference base rides along as plain str: dN/dS's sanity check reads it
    assert df["ref_base"].tolist() == ["A", "A"] and df["ref_base"].dtype == object
    assert list(df.columns) == [
        "contig",
        "position",
        "ref_base",
        "gene_id",
        "total_coverage",
        "A",
        "C",
        "G",
        "T",
    ]
    assert df["position"].tolist() == [0, 1] and df["C"].tolist() == [0, 8]


def test_sample_absent_from_the_cache_is_none(cache):
    assert load_sample_counts(cache, "MAG_A", "pre", "s9_pre") is None
