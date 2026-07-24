#!/usr/bin/env python3
"""Regression tests for regional_contrast.aggregate_region_scores.

Guards against the categorical-groupby explosion: when a grouping column is a
pandas Categorical (as ``group`` is when loaded from parquet), a groupby that
does not pass ``observed=True`` reindexes the result to the full Cartesian
product of every grouping key's levels.  With the collinear region-metadata
columns (region_id/region_start/region_end/region_length) this exploded to a
23.7 PiB allocation in production.
"""

import pandas as pd
import pytest

from alleleflux.scripts.analysis.regional_contrast import (
    FISHER_PVAL_CONTROL_COL,
    FISHER_PVAL_TREATMENT_COL,
    _fisher_merge_keys,
    aggregate_region_scores,
    fisher_combine_empirical_pvalues,
)
from alleleflux.scripts.analysis.regional_contrast import (
    test_region_contrasts as run_region_contrasts,
)


def _make_inputs():
    """Two regions on one contig where each group is observed in only one region.

    Observed (host, group, region) combinations are exactly two.  Under
    ``observed=False`` a Categorical ``group`` grouper forces pandas to emit
    phantom rows for the unobserved combinations.
    """
    df = pd.DataFrame(
        {
            "replicate": [1, 1, 1, 1],
            # Categorical, mirroring how parquet loads the group column.
            "group": pd.Categorical(["A", "A", "B", "B"], categories=["A", "B"]),
            "contig": ["c1", "c1", "c1", "c1"],
            "position": [1, 2, 3, 4],
            "site_score": [0.1, 0.2, 0.3, 0.4],
        }
    )
    region_mapping = pd.DataFrame(
        {
            "contig": ["c1", "c1", "c1", "c1"],
            "position": [1, 2, 3, 4],
            "region_id": ["r1", "r1", "r2", "r2"],
            "region_type": ["gene", "gene", "gene", "gene"],
            "region_start": [1, 1, 3, 3],
            "region_end": [2, 2, 4, 4],
            "region_length": [2, 2, 2, 2],
        }
    )
    return df, region_mapping


def test_categorical_group_does_not_explode():
    df, region_mapping = _make_inputs()

    result = aggregate_region_scores(
        df,
        region_mapping,
        host_col="replicate",
        group_col="group",
        contig_col="contig",
        position_col="position",
        score_col="site_score",
        agg_method="median",
        min_sites=0,
        min_fraction=0.0,
    )

    # Only the two observed (host, group, region) combinations should appear.
    assert len(result) == 2
    # No phantom groups: every emitted region must have a real score.
    assert not result["region_score"].isna().any()
    assert set(zip(result["group"], result["region_id"])) == {("A", "r1"), ("B", "r2")}


def _make_paired_df():
    """Two regions, three replicates each, with treatment/control percentiles."""
    return pd.DataFrame(
        {
            "region_id": ["r1", "r1", "r1", "r2", "r2", "r2"],
            "region_type": ["gene"] * 6,
            "contig": ["c1"] * 6,
            "region_start": [1, 1, 1, 100, 100, 100],
            "region_end": [50, 50, 50, 150, 150, 150],
            "replicate": [1, 2, 3, 1, 2, 3],
            "contrast": [0.2, 0.3, 0.25, -0.1, -0.2, -0.15],
            "percentile_treatment": [80.0, 90.0, 85.0, 40.0, 30.0, 35.0],
            "percentile_control": [20.0, 10.0, 15.0, 60.0, 70.0, 65.0],
        }
    )


def test_summary_fisher_merge_joins_on_region_keys_only():
    """summary_df and fisher_df must merge on region identity, not replicate counts.

    fisher_df carries per-group ``n_replicates_treatment`` / ``n_replicates_control``
    columns that summary_df lacks.  Joining on every non-p-value fisher column (the
    old behavior) referenced a column absent from summary_df → KeyError in production.
    """
    paired = _make_paired_df()
    summary = run_region_contrasts(
        paired, host_col="replicate", contig_col="contig", min_replicates=2
    )
    fisher = fisher_combine_empirical_pvalues(
        paired, contig_col="contig", min_replicates=2
    )
    assert not fisher.empty
    assert "n_replicates_treatment" in fisher.columns
    assert "n_replicates_treatment" not in summary.columns

    # Reproduce the original (buggy) key selection: it raises the production KeyError.
    naive_keys = [
        c
        for c in fisher.columns
        if c not in (FISHER_PVAL_TREATMENT_COL, FISHER_PVAL_CONTROL_COL)
    ]
    with pytest.raises(KeyError):
        summary.merge(fisher, on=naive_keys, how="outer")

    # The fixed key selection joins on region-identity columns only and succeeds.
    keys = _fisher_merge_keys(fisher, summary)
    assert set(keys) == {
        "region_id",
        "region_type",
        "contig",
        "region_start",
        "region_end",
    }
    merged = summary.merge(fisher, on=keys, how="outer")
    assert len(merged) == 2  # one row per region, no duplication
    assert FISHER_PVAL_TREATMENT_COL in merged.columns
    assert "n_replicates" in merged.columns  # summary payload preserved
    assert "n_replicates_treatment" in merged.columns  # fisher payload carried along
