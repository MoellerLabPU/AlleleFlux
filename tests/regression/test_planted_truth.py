"""Planted-truth assertions: does the pipeline find what the generator planted?

The golden test (``test_pipeline_e2e.py``) proves the outputs did not CHANGE;
this file proves they are RIGHT, for the facts the synthetic-data generator
deliberately planted and recorded in ``planted_truth.json`` beside each
dataset (see ``docs/source/examples/generate_synthetic_data.py``,
``plan_the_plants``).  Every assertion reads an id from that card and looks
it up in the run's outputs, so a renamed file or column fails loudly here
while a passing run means the science is intact.

The pipeline run is shared with the golden test through the ``pipeline_run``
fixture in ``conftest.py`` (one run per scenario per session), so this file
adds seconds, not minutes.

What is asserted, per scenario:

* every planted significant site is significant (q < 0.05) in the family the
  baseline-presence step reads (paired t-test on longitudinal, unpaired on
  single), and its rising allele really rose in the treatment group;
* longitudinal only: the standing sites are ``standing_variation`` with the
  planted control mouse named at baseline; the lone de novo site is
  ``de_novo_candidate`` / ``de_novo`` with every carrier absent at its own
  baseline; the three block sites are one ``block_event``; the swapped mouse
  is the only strain change (turnover table, pairwise ANI, classification);
* the depth spikes are removed for exactly the spiked samples
  (per-sample sidecar + summary);
* longitudinal only: the hot gene (one de novo event per strain-keeping
  treatment mouse, none in control) is the only gene of its genome with
  events, its between-group test is significant and its events are
  ``de_novo`` in their cages, the block's positions are ``block_event`` in
  every treatment cage; the strain-swapped cage contributes no de novo
  event at all (its burst of new alleles is folded into blocks);
* single only: pairwise ANI compared EVERY pair, and only the pairs touching
  the divergent sample fall below the same-strain line.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

CON_THRESHOLD = 0.999  # strain_turnover.con_threshold in the generated configs
Q_THRESHOLD = 0.05  # statistics.p_value_threshold / baseline_presence.threshold


# --- helpers ------------------------------------------------------------------


def _card(run) -> dict:
    return json.loads((Path(run.example_dir) / "planted_truth.json").read_text())


def _comparison(card) -> str:
    """The one {timepoints}-{groups} label the generated configs define."""
    tps = "_".join(card["timepoints"])
    return f"{tps}-{card['groups']['treatment']}_{card['groups']['control']}"


def _site_key(site) -> tuple:
    return (site["mag"], site["contig"], site["position"])


def _subject_of(sample_id: str) -> str:
    """``treatment_subj5_post`` -> ``subj5``; ``treatment_sample5`` -> ``subj5``."""
    stem = sample_id.split("_")[1]
    return stem.replace("sample", "subj")


def _read(path: Path, **kwargs) -> pd.DataFrame:
    assert path.exists(), f"expected output missing: {path}"
    return pd.read_csv(path, sep="\t", **kwargs)


# --- significance ----------------------------------------------------------------


@pytest.mark.regression
def test_every_planted_site_is_significant_and_rose_in_treatment(pipeline_run):
    card = _card(pipeline_run)
    comp = _comparison(card)
    family = (
        "two_sample_paired"
        if card["data_type"] == "longitudinal"
        else "two_sample_unpaired"
    )
    summary = _read(
        pipeline_run.pipeline_data_dir
        / "p_value_summary"
        / comp
        / f"p_value_summary_{family}_{comp}.tsv"
    )
    sig = summary[summary["q_value"] < Q_THRESHOLD]
    found = {(r.mag_id, r.contig, int(r.position)) for r in sig.itertuples()}
    missed = [s for s in card["significant_sites"] if _site_key(s) not in found]
    assert not missed, f"planted sites not significant in {family}: {missed}"

    # Direction: the rising allele's frequency in the treatment carriers is
    # above its frequency everywhere else.  Read from the allele-analysis
    # tables, one row per (sample, position).
    aa_dir = (
        pipeline_run.pipeline_data_dir / "allele_analysis" / f"allele_analysis_{comp}"
    )
    for site in card["significant_sites"]:
        if card["data_type"] == "longitudinal":
            tbl = pd.read_parquet(
                aa_dir / f"{site['mag']}_allele_frequency_changes_no_zero-diff.parquet"
            )
            rows = tbl[
                (tbl["contig"] == site["contig"])
                & (tbl["position"] == site["position"])
            ]
            later = card["timepoints"][-1]
            col = f"{site['rising_allele']}_frequency_{later}"
            treat = rows[rows["group"].astype(str) == card["groups"]["treatment"]][col]
            ctrl = rows[rows["group"].astype(str) == card["groups"]["control"]][col]
        else:
            tbl = pd.read_parquet(
                aa_dir / f"{site['mag']}_allele_frequency_single.parquet"
            )
            rows = tbl[
                (tbl["contig"] == site["contig"])
                & (tbl["position"] == site["position"])
            ]
            col = f"{site['rising_allele']}_frequency"
            treat = rows[rows["group"].astype(str) == card["groups"]["treatment"]][col]
            ctrl = rows[rows["group"].astype(str) == card["groups"]["control"]][col]
        assert len(treat) and len(ctrl), f"no rows for {_site_key(site)}"
        assert treat.min() > ctrl.max(), (
            f"{_site_key(site)}: rising allele {site['rising_allele']} not higher in every "
            f"treatment sample ({treat.tolist()}) than in any control sample ({ctrl.tolist()})"
        )


# --- baseline presence / rising origin (longitudinal) --------------------------------


def _baseline_summary(run, card) -> pd.DataFrame:
    comp = _comparison(card)
    return _read(
        run.pipeline_data_dir
        / "baseline_presence"
        / f"{comp}_two_sample_paired_tTest_baseline_presence_summary.tsv"
    )


def _rising_rows(summary: pd.DataFrame, card) -> dict:
    """The summary row of each planted site for its RISING allele, indexed by site key."""
    out = {}
    for site in card["significant_sites"]:
        m = summary[
            (summary["mag_id"] == site["mag"])
            & (summary["contig"] == site["contig"])
            & (summary["position"] == site["position"])
            & (summary["allele"] == site["rising_allele"])
        ]
        assert (
            len(m) == 1
        ), f"{_site_key(site)}: expected one summary row for allele {site['rising_allele']}, got {len(m)}"
        out[_site_key(site)] = m.iloc[0]
    return out


@pytest.mark.regression
def test_standing_sites_are_standing_and_name_the_planted_mouse(pipeline_run):
    card = _card(pipeline_run)
    if card["data_type"] != "longitudinal":
        pytest.skip("baseline presence is longitudinal-only")
    rows = _rising_rows(_baseline_summary(pipeline_run, card), card)
    mouse = card["standing_mouse"]["subjectID"]
    for site in card["significant_sites"]:
        if site["kind"] != "standing":
            continue
        row = rows[_site_key(site)]
        assert row["origin_any_mouse"] == "standing_variation", (
            _site_key(site),
            row["origin_any_mouse"],
        )
        assert row["pre_mice_allele_present"] == mouse, (
            _site_key(site),
            row["pre_mice_allele_present"],
        )


@pytest.mark.regression
def test_lone_de_novo_site_is_de_novo_with_every_carrier_absent_before(pipeline_run):
    card = _card(pipeline_run)
    if card["data_type"] != "longitudinal":
        pytest.skip("baseline presence is longitudinal-only")
    lone = next(s for s in card["significant_sites"] if s["label"] == "lone")
    row = _rising_rows(_baseline_summary(pipeline_run, card), card)[_site_key(lone)]
    carriers = sorted(_subject_of(s) for s in card["carrier_samples"])
    assert row["origin_any_mouse"] == "de_novo_candidate"
    assert pd.isna(row["pre_mice_allele_present"])
    assert sorted(row["post_mice_allele_present"].split(",")) == carriers
    assert int(row["n_carriers_own_baseline_absent"]) == len(carriers)

    alleles = _read(
        pipeline_run.pipeline_data_dir
        / "baseline_presence"
        / "rising_origin"
        / "rising_origin_alleles.tsv"
    )
    a = alleles[
        (alleles["contig"] == lone["contig"])
        & (alleles["position"] == lone["position"])
    ]
    assert len(a) == 1 and a.iloc[0]["allele"] == lone["rising_allele"]
    assert a.iloc[0]["origin"] == "de_novo", a.iloc[0]["origin"]
    assert int(a.iloc[0]["n_de_novo_neighbours"]) == 0
    assert sorted(a.iloc[0]["carriers"].split(",")) == carriers


@pytest.mark.regression
def test_block_sites_are_one_block_event(pipeline_run):
    card = _card(pipeline_run)
    if card["data_type"] != "longitudinal":
        pytest.skip("rising origin is longitudinal-only")
    block = [s for s in card["significant_sites"] if s["label"] == "block"]
    assert len(block) == 3
    alleles = _read(
        pipeline_run.pipeline_data_dir
        / "baseline_presence"
        / "rising_origin"
        / "rising_origin_alleles.tsv"
    )
    rows = pd.concat(
        alleles[
            (alleles["contig"] == s["contig"])
            & (alleles["position"] == s["position"])
            & (alleles["allele"] == s["rising_allele"])
        ]
        for s in block
    )
    assert len(rows) == 3, rows
    assert set(rows["origin"]) == {"block_event"}, rows["origin"].tolist()
    assert rows["block_id"].nunique() == 1, rows["block_id"].tolist()
    # Three sites 300 bp apart inside a 1000 bp window: each sees the other two.
    assert rows["n_de_novo_neighbours"].tolist() == [2, 2, 2]

    summary = _read(
        pipeline_run.pipeline_data_dir
        / "baseline_presence"
        / "rising_origin"
        / "rising_origin_summary.tsv"
    )
    treat = summary[
        (summary["rose_in"] == card["groups"]["treatment"])
        & (summary["strain_class"] == "all")
    ]
    by_origin = treat.set_index("origin")["n_alleles"]
    assert int(by_origin["block_event"]) == 3
    assert int(by_origin["de_novo"]) == 1  # the lone one, and only it


# --- strain swap ------------------------------------------------------------------


@pytest.mark.regression
def test_swapped_mouse_is_the_only_strain_change(pipeline_run):
    card = _card(pipeline_run)
    swap = card["strain_swap"]
    ani_dir = pipeline_run.pipeline_data_dir / "pairwise_ani"
    mags = sorted(
        {s["mag"] for s in card["significant_sites"]} | set(swap["positions"])
    )
    for mag in mags:
        ani = _read(ani_dir / f"{mag}_pairwise_ani.tsv")
        touches = (ani["sample1"] == swap["sample_id"]) | (
            ani["sample2"] == swap["sample_id"]
        )
        assert (ani.loc[touches, "conANI"] < CON_THRESHOLD).all(), ani.loc[
            touches, ["sample1", "sample2", "conANI"]
        ]
        assert (ani.loc[~touches, "conANI"] >= CON_THRESHOLD).all(), ani.loc[
            ~touches & (ani["conANI"] < CON_THRESHOLD), ["sample1", "sample2", "conANI"]
        ]
        if card["data_type"] == "single":
            n = len({*ani["sample1"], *ani["sample2"]})
            assert len(ani) == n * (n - 1) // 2, "pairs: all must compare every pair"
            continue
        # Longitudinal: within-subject pairs only, and the turnover verdicts.
        assert len(ani) == len({*ani["sample1"], *ani["sample2"]}) // 2
        turnover = _read(
            pipeline_run.pipeline_data_dir
            / "strain_turnover"
            / f"{mag}_strain_turnover.tsv"
        )
        changed = turnover[turnover["dominant_strain_change"].astype(str) == "True"][
            "subjectID"
        ].tolist()
        assert changed == [swap["subjectID"]], changed

    if card["data_type"] == "longitudinal":
        cls = _read(
            pipeline_run.pipeline_data_dir
            / "strain_turnover"
            / "replacement_classification.tsv"
        )
        con = cls[cls["metric"] == "dominant_strain_change"]
        # One swapped mouse out of four voters in its group: a minority, so the
        # group is not_replaced; the other group has no change at all.
        assert set(con["strain_status"]) == {"not_replaced"}, con[
            ["MAG_ID", "group", "strain_status"]
        ]
        swapped_group = con[con["group"] == card["groups"]["treatment"]]
        assert (swapped_group["n_mice_changed"] == 1).all()
        assert (
            con[con["group"] == card["groups"]["control"]]["n_mice_changed"] == 0
        ).all()


# --- gene-level de novo accumulation (longitudinal) ---------------------------------


def _gene_de_novo_dir(run, card) -> Path:
    return run.pipeline_data_dir / "gene_de_novo" / _comparison(card)


def _replicate_of(run) -> dict:
    """sample_id -> replicate letter, from the dataset's sample sheet (the card
    names the carrier SAMPLES; the gene tables key cages by replicate)."""
    sheet = pd.read_csv(
        Path(run.example_dir) / "metadata" / "sample_metadata.tsv",
        sep="\t",
        dtype=str,
    )
    return dict(zip(sheet["sample_id"], sheet["replicate"]))


@pytest.mark.regression
def test_hot_gene_accumulates_in_treatment_and_no_other_gene_does(pipeline_run):
    card = _card(pipeline_run)
    if card["data_type"] != "longitudinal":
        pytest.skip("gene-level de novo accumulation is longitudinal-only")
    hot = card["hot_gene"]
    treat, ctrl = card["groups"]["treatment"], card["groups"]["control"]
    d = _gene_de_novo_dir(pipeline_run, card)

    tests = _read(d / f"gene_de_novo_tests_{_comparison(card)}.tsv")
    genome = tests[tests["mag_id"] == hot["mag"]]
    row = genome[genome["gene_id"] == hot["gene_id"]]
    assert len(row) == 1, row
    row = row.iloc[0]
    n_planted = len(hot["events"])
    assert int(row[f"n_de_novo_{treat}"]) == n_planted
    assert int(row[f"n_de_novo_{ctrl}"]) == 0
    # Control cages collected no event anywhere in the genome, so the gene's
    # expected count there is 0 for every cage: nothing to test against.
    assert row[f"note_{ctrl}"] == "no_exposure", row[f"note_{ctrl}"]
    # Test A, treatment: one event per strain-keeping cage where the gene's
    # share of the genome predicts ~0.13 each (chi-square curve).
    assert row[f"p_value_vs_genome_{treat}"] < 0.02, row[f"p_value_vs_genome_{treat}"]
    # Test B: n vs 0 events over matched readable bases.  Raw: chi-square of
    # 2 n ln 2 (3 events: 0.041).  Discounted: the same score on F(1, cages-2),
    # which charges for estimating the dispersion from six spare cages (0.087).
    assert np.isinf(row[f"rate_ratio_{treat}_vs_{ctrl}"])
    assert row["p_value_raw_between"] < 0.05, row["p_value_raw_between"]
    assert row["p_value_between"] < 0.1, row["p_value_between"]
    assert row["note_between"] == "" or pd.isna(row["note_between"])
    # No other gene of that genome collected a single event in either group:
    # the block plant is excluded by rule and the noise never clears the bar.
    others = genome[genome["gene_id"] != hot["gene_id"]]
    assert len(others) > 0
    assert (others[f"n_de_novo_{treat}"] == 0).all(), others
    assert (others[f"n_de_novo_{ctrl}"] == 0).all(), others

    # The events file: each planted position is a clean de novo event in its
    # own cage and nowhere else.
    # (low_memory=False: the strain-status column is blank on most rows and
    # pandas would otherwise warn about mixed types while chunking.)
    events = _read(d / f"{hot['mag']}_gene_de_novo_events.tsv", low_memory=False)
    for e in hot["events"]:
        m = events[
            (events["contig"] == hot["contig"])
            & (events["position"] == e["position"])
            & (events["allele"] == e["rising_allele"])
        ]
        assert len(m) == 1, m
        assert m.iloc[0]["replicate"] == e["replicate"]
        assert m.iloc[0]["group"] == treat
        assert m.iloc[0]["origin"] == "de_novo", m.iloc[0]["origin"]
        assert m.iloc[0]["gene_id"] == hot["gene_id"]
    # The block's three positions are block_event in EVERY treatment cage.
    replicate_of = _replicate_of(pipeline_run)
    treatment_cages = sorted(replicate_of[s] for s in card["carrier_samples"])
    block = [s for s in card["significant_sites"] if s["label"] == "block"]
    block_events = _read(
        d / f"{block[0]['mag']}_gene_de_novo_events.tsv", low_memory=False
    )
    for s in block:
        m = block_events[
            (block_events["contig"] == s["contig"])
            & (block_events["position"] == s["position"])
            & (block_events["allele"] == s["rising_allele"])
            & (block_events["group"] == treat)
        ]
        assert sorted(m["replicate"]) == treatment_cages, m
        assert set(m["origin"]) == {"block_event"}, m["origin"].tolist()


@pytest.mark.regression
def test_strain_swapped_cage_contributes_no_de_novo_event(pipeline_run):
    """A mouse whose strain was replaced shows new alleles at 2 % of every
    contig's positions at once.  Per cage, the cluster rule folds them into
    one block per contig, so that cage's genome-wide de novo total is 0 on
    every genome and the gene tests never see the burst."""
    card = _card(pipeline_run)
    if card["data_type"] != "longitudinal":
        pytest.skip("gene-level de novo accumulation is longitudinal-only")
    swap = card["strain_swap"]
    treat = card["groups"]["treatment"]
    rep = _replicate_of(pipeline_run)[swap["sample_id"]]
    d = _gene_de_novo_dir(pipeline_run, card)
    for mag in sorted(swap["positions"]):
        cages = _read(d / f"{mag}_gene_de_novo_cages.tsv")
        mine = cages[(cages["group"] == treat) & (cages["replicate"] == rep)]
        assert len(mine) == 1, mine
        assert int(mine.iloc[0]["n_de_novo_total"]) == 0, mine
        events = _read(d / f"{mag}_gene_de_novo_events.tsv", low_memory=False)
        ev = events[(events["group"] == treat) & (events["replicate"] == rep)]
        assert not ev["origin"].isin(["de_novo", "de_novo_below_detection"]).any(), ev[
            ev["origin"].isin(["de_novo", "de_novo_below_detection"])
        ]
        # Folded, not missed: the burst is there as block events.
        assert (ev["origin"] == "block_event").sum() >= 10, ev["origin"].value_counts()


# --- depth spikes -----------------------------------------------------------------


@pytest.mark.regression
def test_depth_spikes_are_removed_for_exactly_the_spiked_samples(pipeline_run):
    card = _card(pipeline_run)
    spikes = card["depth_spikes"]
    expected = {(s["sample_id"], s["mag"]): s["end"] - s["start"] for s in spikes}
    summary = _read(
        pipeline_run.pipeline_data_dir
        / "site_depth_filter"
        / "site_depth_filter_summary.tsv"
    )
    for row in summary.itertuples():
        want = expected.get((row.sample_id, row.mag_id), 0)
        assert int(row.n_removed_high) == want, (
            row.sample_id,
            row.mag_id,
            row.n_removed_high,
            want,
        )
        assert int(row.n_removed_low) == 0
    # Every spike appears in the summary (the dict above only checks rows that exist).
    seen = {(r.sample_id, r.mag_id) for r in summary.itertuples()}
    assert set(expected) <= seen, set(expected) - seen
    # The per-timepoint sidecar beside the cache tells the same story.
    for s in spikes:
        tp = (
            next(t for t in card["timepoints"] if s["sample_id"].endswith(f"_{t}"))
            if card["data_type"] == "longitudinal"
            else card["timepoints"][0]
        )
        side = _read(
            pipeline_run.pipeline_data_dir
            / "allele_freq_cache"
            / tp
            / f"{s['mag']}_{tp}_allele_frequency_site_depth_filter.tsv"
        )
        mine = side[side["sample_id"] == s["sample_id"]]
        assert (
            len(mine) == 1
            and int(mine.iloc[0]["n_removed_high"]) == s["end"] - s["start"]
        ), mine
