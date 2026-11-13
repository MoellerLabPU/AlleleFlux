"""track_allele_frequencies reads per-sample anchor-site counts from the cache."""

import numpy as np
import pandas as pd

from alleleflux.scripts.visualization.accessory.track_allele_frequencies import (
    load_timepoint_rows,
    process_single_sample,
)
from tests.utilities.test_cache_io import write_cache

ANCHORS = {("c1", 0): ["A", "C"], ("c1", 5): ["A"], ("c1", 9): ["G"]}
SITES = pd.DataFrame([(c, p) for (c, p) in ANCHORS], columns=["contig", "position"])


def test_anchor_frequencies_from_the_cache(tmp_path):
    d = str(tmp_path / "allele_freq_cache")
    write_cache(
        d,
        "MAGX",
        "t1",
        [
            ("c1", 0, "s1", "g", "m1", 8, 2, 0, 0),
            ("c1", 5, "s1", "g", "m1", 1, 0, 0, 0),
        ],
    )  # below min_cov 5
    # one scan of the t1 file serves every sample at t1
    tp, rows, roster = load_timepoint_rows(("t1", ["s1"], SITES, "MAGX", d))
    assert tp == "t1" and roster == {"s1"}
    sample_id, freqs = process_single_sample(
        "s1", rows[rows["sample_id"] == "s1"], "s1" in roster, ANCHORS, 5
    )
    assert sample_id == "s1"
    assert freqs[("c1", 0, "A")] == 0.8 and freqs[("c1", 0, "C")] == 0.2
    assert np.isnan(freqs[("c1", 5, "A")])  # covered below min_cov -> NaN, as before
    assert np.isnan(freqs[("c1", 9, "G")])  # not in the cache -> NaN


def test_sample_absent_from_the_cache_gives_none(tmp_path):
    d = str(tmp_path / "allele_freq_cache")
    write_cache(d, "MAGX", "t1", [("c1", 0, "s1", "g", "m1", 8, 2, 0, 0)])
    _, rows, roster = load_timepoint_rows(("t1", ["s9"], SITES, "MAGX", d))
    assert rows.empty and "s9" not in roster
    assert process_single_sample("s9", rows, False, ANCHORS, 5) == ("s9", None)
