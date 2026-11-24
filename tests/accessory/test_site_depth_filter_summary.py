"""alleleflux-site-filter-summary: every cache sidecar -> one run table + a per-MAG rollup."""

import warnings

import pandas as pd
import pytest

from alleleflux.scripts.accessory.site_depth_filter_summary import main

COLS = [
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


def _sidecar(path, rows):
    pd.DataFrame(rows, columns=COLS).to_csv(path, sep="\t", index=False)


def test_run_table_and_rollup(tmp_path):
    cache = tmp_path / "allele_freq_cache"
    (cache / "t0").mkdir(parents=True)
    (cache / "t1").mkdir()
    _sidecar(
        cache / "t0" / "MAG_A_t0_allele_frequency_site_depth_filter.tsv",
        [
            ("s1", "MAG_A", "t0", 20.0, 18.0, 1000, 0, 10, 990),
            ("s2", "MAG_A", "t0", 30.0, 26.0, 2000, 0, 30, 1970),
        ],
    )
    _sidecar(
        cache / "t1" / "MAG_A_t1_allele_frequency_site_depth_filter.tsv",
        [
            ("s1", "MAG_A", "t1", 40.0, 40.0, 1000, 0, 0, 1000),
        ],
    )
    out = tmp_path / "site_depth_filter_summary.tsv"
    assert main(["--cache_dir", str(cache), "--output", str(out)]) == 0
    run = pd.read_csv(out, sep="\t")
    assert len(run) == 3 and run.columns.tolist() == COLS
    roll = pd.read_csv(tmp_path / "site_depth_filter_summary_by_mag.tsv", sep="\t")
    assert roll.columns.tolist() == [
        "mag_id",
        "timepoint",
        "n_samples",
        "mean_of_mean_coverage",
        "mean_of_mean_coverage_after_filter",
        "sites_removed_high_total",
        "sites_covered_total",
        "fraction_removed_high",
    ]
    t0 = roll.set_index("timepoint").loc["t0"]
    assert t0["n_samples"] == 2 and t0["mean_of_mean_coverage"] == 25.0
    assert t0["mean_of_mean_coverage_after_filter"] == 22.0
    assert t0["sites_removed_high_total"] == 40 and t0["sites_covered_total"] == 3000
    assert t0["fraction_removed_high"] == pytest.approx(40 / 3000)


def test_no_sidecars_fails_loud(tmp_path):
    (tmp_path / "allele_freq_cache" / "t0").mkdir(parents=True)
    with pytest.raises(FileNotFoundError, match="site_depth_filter"):
        main(
            [
                "--cache_dir",
                str(tmp_path / "allele_freq_cache"),
                "--output",
                str(tmp_path / "o.tsv"),
            ]
        )


def test_header_only_sidecar_adds_no_row_and_keeps_numbers_numeric(tmp_path):
    """A zero-row cache (a genome with no QC-passing sample at a timepoint)
    comes with a header-only sidecar.  It adds no row to the run table, and
    the count columns of the other sidecars stay numbers (an empty table read
    next to them must not turn them into text)."""
    cache = tmp_path / "allele_freq_cache"
    (cache / "t0").mkdir(parents=True)
    (cache / "t1").mkdir()
    _sidecar(
        cache / "t0" / "MAG_A_t0_allele_frequency_site_depth_filter.tsv",
        [("s1", "MAG_A", "t0", 20.0, 18.0, 1000, 0, 10, 990)],
    )
    _sidecar(cache / "t1" / "MAG_A_t1_allele_frequency_site_depth_filter.tsv", [])
    out = tmp_path / "site_depth_filter_summary.tsv"
    # pandas warns that stacking an empty table will change the result's types
    # in a later version; the command must not depend on that behaviour.
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        assert main(["--cache_dir", str(cache), "--output", str(out)]) == 0
    run = pd.read_csv(out, sep="\t")
    assert len(run) == 1 and run.columns.tolist() == COLS
    assert run["n_removed_high"].tolist() == [10]
    assert pd.api.types.is_integer_dtype(run["n_sites_covered"])
    roll = pd.read_csv(tmp_path / "site_depth_filter_summary_by_mag.tsv", sep="\t")
    assert len(roll) == 1 and roll.iloc[0]["sites_covered_total"] == 1000


def test_only_header_only_sidecars_give_header_only_tables(tmp_path):
    """Every cache of the folder is empty: both tables are written with their
    columns and no row (nothing was filtered because nothing was read)."""
    cache = tmp_path / "allele_freq_cache"
    (cache / "t0").mkdir(parents=True)
    _sidecar(cache / "t0" / "MAG_A_t0_allele_frequency_site_depth_filter.tsv", [])
    out = tmp_path / "site_depth_filter_summary.tsv"
    assert main(["--cache_dir", str(cache), "--output", str(out)]) == 0
    run = pd.read_csv(out, sep="\t")
    assert len(run) == 0 and run.columns.tolist() == COLS
    roll = pd.read_csv(tmp_path / "site_depth_filter_summary_by_mag.tsv", sep="\t")
    assert len(roll) == 0 and roll.columns.tolist() == [
        "mag_id",
        "timepoint",
        "n_samples",
        "mean_of_mean_coverage",
        "mean_of_mean_coverage_after_filter",
        "sites_removed_high_total",
        "sites_covered_total",
        "fraction_removed_high",
    ]
