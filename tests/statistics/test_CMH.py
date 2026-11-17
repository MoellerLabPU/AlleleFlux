"""CMH pairs: a replicate counts only when BOTH groups have a row for it at the site.

The CMH test compares two groups replicate by replicate (the replicate is the
stratum).  ``prepare_cmh_input`` builds one site's table of pairs and returns
``None`` when fewer than ``min_sample_num`` replicates are present in both
groups.

The loader hands ``group`` and ``replicate`` over either as plain text or as
pandas categories, depending on whether the two timepoint files happen to
hold the same label sets.  A category remembers every label it has ever seen,
and a pandas groupby on categories fills in a zero row for every label
combination that has no data.  Left alone, that turned a replicate with NO
row in one group into a "pair" (reads on one side, zeros on the other), so
``min_sample_num`` was met by replicates that carry no comparison, and a
category label with no row at all reached R as an empty stratum, which
``mantelhaen.test`` refuses.

Contract pinned here: the table of pairs is the same whichever way the labels
are stored, and it is the one plain-text labels give.

Small example used throughout: one site, replicates 1-4; replicate 4 has
reads in ``fat`` and no row in ``control``.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
import pytest

from alleleflux.scripts.statistics.CMH import (
    NUCLEOTIDES,
    prepare_cmh_input,
    run_cmh_test_in_r,
)

COLUMNS = ["contig", "position", "group", "replicate", *NUCLEOTIDES]


def _site(rows):
    """One site's rows as CMH receives them.

    ``rows``: (group, replicate, A reads, G reads); the other two bases are 0.
    """
    records = []
    for group, replicate, a, g in rows:
        counts = dict.fromkeys(NUCLEOTIDES, 0)
        counts["A"], counts["G"] = a, g
        records.append(
            {"contig": "c1", "position": 10, "group": group, "replicate": replicate}
            | counts
        )
    return pd.DataFrame(records, columns=COLUMNS)


def _as_categories(df, groups=None, replicates=None):
    """The same rows with ``group`` / ``replicate`` as categories.

    ``groups`` / ``replicates`` name the full label sets the category knows
    (more than the rows use, as after a filter); default: the labels present.
    """
    out = df.copy()
    out["group"] = pd.Categorical(out["group"], categories=groups)
    out["replicate"] = pd.Categorical(out["replicate"], categories=replicates)
    return out


def _plain(table):
    """A pairs table with the replicate label as text, for comparing two results."""
    if table is None:
        return None
    return table.astype({"replicate": str}).reset_index(drop=True)


# Replicate 4: reads in fat, no row in control.
ONE_SIDED = [
    ("fat", "1", 20, 10),
    ("fat", "2", 18, 12),
    ("fat", "3", 25, 5),
    ("fat", "4", 15, 15),
    ("control", "1", 30, 0),
    ("control", "2", 28, 2),
    ("control", "3", 29, 1),
]


class TestPairsAreReplicatesInBothGroups:
    @pytest.mark.parametrize("as_categories", [False, True])
    def test_replicate_missing_from_one_group_is_not_a_pair(self, as_categories):
        site = _site(ONE_SIDED)
        if as_categories:
            site = _as_categories(site)
        pairs = prepare_cmh_input(site, "fat", "control", min_sample_num=2)
        assert sorted(pairs["replicate"].astype(str)) == ["1", "2", "3"]
        assert len(pairs) == 3  # what the CMH table reports as num_pairs

    @pytest.mark.parametrize("as_categories", [False, True])
    def test_the_floor_is_met_by_real_pairs_only(self, as_categories):
        """Four replicates, floor 4, one of them present in one group only:
        three real pairs, so the site is not tested."""
        site = _site(ONE_SIDED)
        if as_categories:
            site = _as_categories(site)
        assert prepare_cmh_input(site, "fat", "control", min_sample_num=4) is None

    @pytest.mark.parametrize("as_categories", [False, True])
    def test_complete_site_keeps_every_replicate(self, as_categories):
        site = _site(ONE_SIDED + [("control", "4", 27, 3)])
        if as_categories:
            site = _as_categories(site)
        pairs = prepare_cmh_input(site, "fat", "control", min_sample_num=4)
        assert sorted(pairs["replicate"].astype(str)) == ["1", "2", "3", "4"]
        row = pairs[pairs["replicate"].astype(str) == "4"].iloc[0]
        assert (row["fat_A"], row["fat_G"], row["control_A"], row["control_G"]) == (
            15,
            15,
            27,
            3,
        )

    def test_labels_the_rows_do_not_use_change_nothing(self):
        """After the group filter a category still lists the labels that were
        filtered out (another group, its replicates).  They add no pair."""
        site = _site(ONE_SIDED)
        wide = _as_categories(
            site,
            groups=["control", "fat", "other"],
            replicates=["1", "2", "3", "4", "5", "9"],
        )
        got = prepare_cmh_input(wide, "fat", "control", min_sample_num=2)
        want = prepare_cmh_input(site, "fat", "control", min_sample_num=2)
        pd.testing.assert_frame_equal(_plain(got), _plain(want))

    def test_across_time_shape_group_as_text_replicate_as_category(self):
        """The across-time test relabels ``group`` with the timepoint (text)
        and leaves ``replicate`` as it was loaded (a category)."""
        site = _site(
            [
                ("pre", "A", 30, 0),
                ("pre", "B", 28, 2),
                ("pre", "C", 29, 1),
                ("end", "A", 20, 10),
                ("end", "B", 18, 12),
            ]
        )
        mixed = site.copy()
        mixed["replicate"] = pd.Categorical(
            mixed["replicate"], categories=["A", "B", "C", "D"]
        )
        got = prepare_cmh_input(mixed, "pre", "end", min_sample_num=2)
        assert sorted(got["replicate"].astype(str)) == ["A", "B"]
        pd.testing.assert_frame_equal(
            _plain(got), _plain(prepare_cmh_input(site, "pre", "end", 2))
        )


class TestSameTableHoweverTheLabelsAreStored:
    def test_every_pattern_of_missing_rows_gives_the_text_result(self):
        """Brute force: 3 replicates x 2 groups = 6 possible rows; every one
        of the 64 subsets of rows, at floors 1, 2 and 3.  The category input
        (with labels the rows do not use) must give exactly the text result,
        ``None`` included."""
        cells = [(g, r) for g in ("fat", "control") for r in ("1", "2", "3")]
        reads = {
            ("fat", "1"): (20, 10),
            ("fat", "2"): (18, 12),
            ("fat", "3"): (25, 5),
            ("control", "1"): (30, 0),
            ("control", "2"): (28, 2),
            ("control", "3"): (29, 1),
        }
        checked = 0
        for present in itertools.product([False, True], repeat=len(cells)):
            rows = [(g, r, *reads[(g, r)]) for (g, r), on in zip(cells, present) if on]
            if not rows:
                continue
            site = _site(rows)
            wide = _as_categories(
                site,
                groups=["control", "fat", "other"],
                replicates=["1", "2", "3", "7"],
            )
            for floor in (1, 2, 3):
                want = prepare_cmh_input(site, "fat", "control", floor)
                got = prepare_cmh_input(wide, "fat", "control", floor)
                assert (got is None) == (want is None), (rows, floor)
                if want is not None:
                    pd.testing.assert_frame_equal(_plain(got), _plain(want))
                checked += 1
        assert checked == 63 * 3

    def test_number_of_pairs_never_exceeds_replicates_in_both_groups(self):
        """Random sites with holes: the pairs returned are always a subset of
        the replicates that have a row in both groups."""
        rng = np.random.default_rng(7)
        replicates = [str(i) for i in range(1, 9)]
        for _ in range(200):
            rows = [
                (g, r, int(rng.integers(0, 40)), int(rng.integers(0, 40)))
                for g in ("fat", "control")
                for r in replicates
                if rng.random() > 0.25
            ]
            site = _as_categories(
                _site(rows), groups=["control", "fat"], replicates=replicates
            )
            in_both = {r for g, r, *_ in rows if g == "fat"} & {
                r for g, r, *_ in rows if g == "control"
            }
            pairs = prepare_cmh_input(site, "fat", "control", min_sample_num=1)
            if pairs is not None:
                assert set(pairs["replicate"].astype(str)) <= in_both


class TestWhatReachesR:
    def test_category_site_with_a_one_sided_replicate_gets_the_text_p_value(self):
        """End to end for one site: the table built from category labels goes
        through R's ``mantelhaen.test`` and returns the p-value of the three
        real pairs.  A leftover category label would reach R as an empty
        stratum, which it refuses ("sample size in each stratum must be > 1")."""
        site = _site(ONE_SIDED)
        text_p = run_cmh_test_in_r(prepare_cmh_input(site, "fat", "control", 2))
        wide = _as_categories(
            site,
            groups=["control", "fat", "other"],
            replicates=["1", "2", "3", "4", "9"],
        )
        category_p = run_cmh_test_in_r(prepare_cmh_input(wide, "fat", "control", 2))
        assert category_p == pytest.approx(text_p, rel=1e-12)
        assert 0 < text_p < 1
