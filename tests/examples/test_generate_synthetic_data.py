"""Tests for the planted truths of ``docs/source/examples/generate_synthetic_data.py``.

The generator is the regression test's only source of ground truth: it decides
every significant site, standing allele, de novo allele, block, strain swap and
depth spike up front (``plan_the_plants``), writes them into the profiles
(``generate_profile_data``) and records them in ``planted_truth.json``
(``generate_dataset``).  These tests pin the contract between the three.
"""

from __future__ import annotations

import importlib.util
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

GENERATOR = (
    Path(__file__).resolve().parents[2]
    / "docs"
    / "source"
    / "examples"
    / "generate_synthetic_data.py"
)


def _load_generator():
    spec = importlib.util.spec_from_file_location("generate_synthetic_data", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


gen = _load_generator()

BASE_COLUMNS = ("A", "C", "G", "T")


def _rows_by_position(rows):
    """Profile row strings of one (sample, MAG) -> {(contig, pos): counts dict}."""
    out = {}
    for row in rows:
        contig, pos, ref, cov, a, c, g, t, gene = row.split("\t")
        out[(contig, int(pos))] = {
            "ref": ref,
            "cov": int(cov),
            "A": int(a),
            "C": int(c),
            "G": int(g),
            "T": int(t),
        }
    return out


def _dominant(counts):
    return max(BASE_COLUMNS, key=lambda b: counts[b])


@pytest.fixture
def longitudinal():
    """Two small MAGs, 8 mice (16 samples), the plan and every sample's profiles."""
    random.seed(7)
    np.random.seed(7)
    mags = [
        gen.create_mag_structure(
            f"MAG_{i:03d}", num_contigs=2, contig_length_range=(1200, 1500)
        )
        for i in (1, 2)
    ]
    samples = gen.generate_sample_structure(16, data_type="longitudinal")
    truth = gen.plan_the_plants(mags, samples, "longitudinal")
    profiles = {
        s["sample_id"]: gen.generate_profile_data(s, mags, truth=truth) for s in samples
    }
    return mags, samples, truth, profiles


class TestPlan:
    def test_longitudinal_plan_lists_standing_lone_and_block_sites(self, longitudinal):
        _, _, truth, _ = longitudinal
        kinds = {}
        for site in truth["significant_sites"]:
            kinds.setdefault((site["kind"], site["label"]), []).append(site)
        assert len(kinds[("standing", "standing")]) == truth["n_standing"]
        assert len(kinds[("de_novo", "lone")]) == 1
        block = sorted(kinds[("de_novo", "block")], key=lambda s: s["position"])
        assert len(block) == 3
        assert len({(s["mag"], s["contig"]) for s in block}) == 1
        spacing = truth["block_spacing_bp"]
        assert [s["position"] - block[0]["position"] for s in block] == [
            0,
            spacing,
            2 * spacing,
        ]
        # The lone de novo allele must never be clustered with the block: it
        # lives on a different contig, so the window cannot reach it.
        lone = kinds[("de_novo", "lone")][0]
        assert (lone["mag"], lone["contig"]) != (block[0]["mag"], block[0]["contig"])

    def test_every_site_names_the_rising_allele_as_the_next_base_in_the_cycle(
        self, longitudinal
    ):
        _, _, truth, _ = longitudinal
        for site in truth["significant_sites"]:
            ref_idx = gen.BASES.index(site["ref_base"])
            assert site["rising_allele"] == gen.BASES[(ref_idx + 1) % 4]

    def test_single_plan_has_shift_sites_and_a_divergent_treatment_sample_only(self):
        random.seed(3)
        np.random.seed(3)
        mags = [gen.create_mag_structure("MAG_001", contig_length_range=(1200, 1500))]
        samples = gen.generate_sample_structure(8, data_type="single")
        truth = gen.plan_the_plants(mags, samples, "single")
        assert {s["kind"] for s in truth["significant_sites"]} == {"shift"}
        assert "standing_mouse" not in truth
        assert truth["strain_swap"]["sample_id"].startswith("treatment_")
        assert len(truth["depth_spikes"]) == 2

    def test_hot_gene_has_one_event_per_strain_keeping_treatment_mouse(
        self, longitudinal
    ):
        """The hot gene is the longest gene on a contig that holds no other de
        novo plant; each treatment mouse that keeps its strain gets ONE event
        inside it.  The strain-swapped mouse gets none: its cage shows a burst
        of new alleles on every contig, which the per-cage cluster rule folds
        into blocks, and a single event planted there would be folded too."""
        mags, samples, truth, _ = longitudinal
        hot = truth["hot_gene"]
        contig = next(
            c
            for m in mags
            if m["mag_id"] == hot["mag"]
            for c in m["contigs"]
            if c["contig_id"] == hot["contig"]
        )
        gene = next(g for g in contig["genes"] if g["gene_id"] == hot["gene_id"])
        assert gene["end"] - gene["start"] == max(
            g["end"] - g["start"] for g in contig["genes"]
        )
        de_novo_contigs = {
            (s["mag"], s["contig"])
            for s in truth["significant_sites"]
            if s["kind"] == "de_novo"
        }
        assert (hot["mag"], hot["contig"]) not in de_novo_contigs
        events = hot["events"]
        positions = [e["position"] for e in events]
        assert all(gene["start"] <= p < gene["end"] for p in positions)
        assert len(set(positions)) == len(positions)
        swapped = truth["strain_swap"]["subjectID"]
        later = truth["timepoints"][-1]
        keepers = sorted(
            s["subjectID"]
            for s in samples
            if s["group"] == truth["groups"]["treatment"]
            and s["time"] == later
            and s["subjectID"] != swapped
        )
        assert sorted(e["subjectID"] for e in events) == keepers
        assert len({e["replicate"] for e in events}) == len(events)
        for e in events:
            assert e["sample_id"] in truth["carrier_samples"]
            assert e["ref_base"] == contig["sequence"][e["position"]]
            ref_idx = gen.BASES.index(e["ref_base"])
            assert e["rising_allele"] == gen.BASES[(ref_idx + 1) % 4]

    def test_plants_never_overlap(self, longitudinal):
        _, _, truth, _ = longitudinal
        sites = {
            (s["mag"], s["contig"], s["position"]) for s in truth["significant_sites"]
        }
        assert len(sites) == len(truth["significant_sites"])
        hot = truth["hot_gene"]
        hot_keys = {(hot["mag"], hot["contig"], e["position"]) for e in hot["events"]}
        assert not sites & hot_keys
        sites |= hot_keys
        for spike in truth["depth_spikes"]:
            for pos in range(spike["start"], spike["end"]):
                assert (spike["mag"], spike["contig"], pos) not in sites
        swap = truth["strain_swap"]
        for mag, contigs in swap["positions"].items():
            for contig, positions in contigs.items():
                assert not sites & {(mag, contig, p) for p in positions}
        assert swap["n_positions"] == sum(
            len(p) for contigs in swap["positions"].values() for p in contigs.values()
        )


class TestProfiles:
    def test_de_novo_site_has_no_rising_reads_before_and_about_half_after(
        self, longitudinal
    ):
        _, samples, truth, profiles = longitudinal
        lone = next(s for s in truth["significant_sites"] if s["label"] == "lone")
        key = (lone["contig"], lone["position"])
        for sample in samples:
            counts = _rows_by_position(profiles[sample["sample_id"]][lone["mag"]])[key]
            rising = counts[lone["rising_allele"]]
            if sample["sample_id"] in truth["carrier_samples"]:
                assert (
                    abs(rising / counts["cov"] - truth["shift_rising_fraction"]) < 0.03
                )
                # The consensus must stay the reference so the shift does not
                # count as a strain change in the ANI step.
                assert _dominant(counts) == lone["ref_base"]
            else:
                assert rising == 0

    def test_hot_gene_allele_has_zero_reads_outside_its_one_sample(self, longitudinal):
        _, samples, truth, profiles = longitudinal
        hot = truth["hot_gene"]
        for event in hot["events"]:
            key = (hot["contig"], event["position"])
            for sample in samples:
                counts = _rows_by_position(profiles[sample["sample_id"]][hot["mag"]])[
                    key
                ]
                rising = counts[event["rising_allele"]]
                if sample["sample_id"] == event["sample_id"]:
                    assert (
                        abs(rising / counts["cov"] - truth["shift_rising_fraction"])
                        < 0.03
                    )
                    assert _dominant(counts) == event["ref_base"]
                else:
                    assert rising == 0

    def test_background_noise_is_at_most_one_read_per_minor_base(self):
        """Sequencer noise is one stray read, not three: three reads of one
        letter at depth 20-40 clear the baseline-presence bar and read as a
        real minor allele, so every low-depth position became a fake de novo
        candidate in the gene-level analysis."""
        random.seed(1)
        seen = set()
        for _ in range(500):
            counts = gen._background_counts(
                "A", random.randint(20, 80), ["T", "G", "C"]
            )
            for b in ("T", "G", "C"):
                assert counts[b] <= 1, counts
                seen.add(counts[b])
        assert seen == {0, 1}  # the noise is still there, just one read deep

    def test_standing_allele_persists_in_one_control_mouse_at_both_timepoints(
        self, longitudinal
    ):
        """The standing allele is a stable minority in ONE control mouse, present
        at pre AND post at the same share, so that mouse's pre-to-post change is
        ~0 and the control group keeps zero variance for the paired tests."""
        _, samples, truth, profiles = longitudinal
        mouse = truth["standing_mouse"]
        assert mouse["subjectID"].startswith("subj")
        assert set(mouse["sample_ids"]) == {
            s["sample_id"] for s in samples if s["subjectID"] == mouse["subjectID"]
        }
        assert all(sid.startswith("control_") for sid in mouse["sample_ids"])
        for site in truth["significant_sites"]:
            if site["kind"] != "standing":
                continue
            key = (site["contig"], site["position"])
            for sample in samples:
                if sample["sample_id"] in truth["carrier_samples"]:
                    continue
                counts = _rows_by_position(profiles[sample["sample_id"]][site["mag"]])[
                    key
                ]
                rising = counts[site["rising_allele"]]
                if sample["sample_id"] in mouse["sample_ids"]:
                    assert abs(rising / counts["cov"] - mouse["fraction"]) < 0.03
                    assert rising >= 5  # unmistakably present, never a 1-read fluke
                else:
                    assert rising == 0

    def test_swapped_sample_changes_consensus_exactly_at_the_swap_positions(
        self, longitudinal
    ):
        _, samples, truth, profiles = longitudinal
        swap = truth["strain_swap"]
        pre_id = swap["sample_id"].replace(f"_{swap['timepoint']}", "_pre")
        for mag, contigs in swap["positions"].items():
            post = _rows_by_position(profiles[swap["sample_id"]][mag])
            pre = _rows_by_position(profiles[pre_id][mag])
            changed = {
                (contig, pos)
                for (contig, pos), counts in post.items()
                if _dominant(counts) != _dominant(pre[(contig, pos)])
            }
            planned = {(contig, pos) for contig, ps in contigs.items() for pos in ps}
            assert changed == planned

    def test_spike_rows_carry_tenfold_coverage_only_inside_the_stretch(
        self, longitudinal
    ):
        _, _, truth, profiles = longitudinal
        for spike in truth["depth_spikes"]:
            rows = _rows_by_position(profiles[spike["sample_id"]][spike["mag"]])
            inside = range(spike["start"], spike["end"])
            for (contig, pos), counts in rows.items():
                if contig == spike["contig"] and pos in inside:
                    assert counts["cov"] == spike["depth"]
                else:
                    assert counts["cov"] <= gen.DEFAULT_COVERAGE_RANGE[1]


class TestDataset:
    def test_card_ids_resolve_to_profile_rows_and_metadata_has_dates(self, tmp_path):
        truth = gen.generate_dataset(
            tmp_path / "ds",
            data_type="longitudinal",
            num_samples=16,
            contig_length_range=(1200, 1500),
            seed=11,
        )
        card = json.loads((tmp_path / "ds" / "planted_truth.json").read_text())
        assert card == truth
        carrier = card["carrier_samples"][0]
        for site in card["significant_sites"]:
            path = (
                tmp_path
                / "ds"
                / "profiles"
                / carrier
                / f"{carrier}_{site['mag']}_profiled.tsv.gz"
            )
            import gzip

            with gzip.open(path, "rt") as fh:
                next(fh)
                rows = _rows_by_position([ln.rstrip("\n") for ln in fh])
            counts = rows[(site["contig"], site["position"])]
            assert (
                abs(
                    counts[site["rising_allele"]] / counts["cov"]
                    - card["shift_rising_fraction"]
                )
                < 0.03
            )
        meta = (
            (tmp_path / "ds" / "metadata" / "sample_metadata.tsv")
            .read_text()
            .splitlines()
        )
        header = meta[0].split("\t")
        assert "collection_date" in header
        col = header.index("collection_date")
        dates = {ln.split("\t")[col] for ln in meta[1:]}
        assert dates == set(card["collection_dates"].values())
        assert all(len(d) == 10 and d[4] == "-" and d[7] == "-" for d in dates)

    def test_same_seed_twice_is_byte_identical(self, tmp_path):
        kwargs = dict(
            data_type="single", num_samples=8, contig_length_range=(1000, 1200), seed=5
        )
        gen.generate_dataset(tmp_path / "a", **kwargs)
        gen.generate_dataset(tmp_path / "b", **kwargs)
        files_a = sorted(
            p.relative_to(tmp_path / "a")
            for p in (tmp_path / "a").rglob("*")
            if p.is_file()
        )
        files_b = sorted(
            p.relative_to(tmp_path / "b")
            for p in (tmp_path / "b").rglob("*")
            if p.is_file()
        )
        assert files_a == files_b
        for rel in files_a:
            if rel.name == "README.md":  # carries the absolute output path
                continue
            a = (tmp_path / "a" / rel).read_bytes()
            b = (tmp_path / "b" / rel).read_bytes()
            if rel.suffix == ".gz":
                import gzip

                a, b = gzip.decompress(a), gzip.decompress(b)
            assert a == b, rel


def _key_paths(d, prefix=()):
    """Every nested key path of a dict as a tuple, e.g. ('analysis', 'use_lmm')."""
    out = set()
    for k, v in d.items():
        path = prefix + (k,)
        out.add(path)
        if isinstance(v, dict):
            out |= _key_paths(v, path)
    return out


TEMPLATE = (
    Path(__file__).resolve().parents[2]
    / "alleleflux"
    / "smk_workflow"
    / "config.template.yml"
)


class TestWriteConfig:
    """The emitted configs switch on every feature the data type supports."""

    def _emit(self, tmp_path, data_type, mode="profiles"):
        import yaml

        timepoints = ["t1"] if data_type == "single" else ["pre", "post"]
        path = gen.write_config(
            tmp_path, data_type, ["control", "treatment"], timepoints, mode=mode
        )
        return yaml.safe_load(path.read_text())

    @pytest.mark.parametrize("data_type", ["longitudinal", "single"])
    @pytest.mark.parametrize("mode", ["profiles", "bams"])
    def test_every_emitted_key_exists_in_the_template(self, tmp_path, data_type, mode):
        import yaml

        template = yaml.safe_load(TEMPLATE.read_text())
        cfg = self._emit(tmp_path, data_type, mode)
        # run_name / log_level / input.profiles_path are documented top-level
        # knobs the template does not spell out; everything else must match.
        known_extra = {("run_name",), ("log_level",), ("input", "profiles_path")}
        unknown = _key_paths(cfg) - _key_paths(template) - known_extra
        assert not unknown, sorted(unknown)

    def test_longitudinal_switches_on_the_whole_longitudinal_surface(self, tmp_path):
        a = self._emit(tmp_path, "longitudinal")["analysis"]
        assert a["use_lmm"] is False  # ~58 CPU-min on this design; unit-tested instead
        assert a["use_cmh"] is True
        assert a["figures"]["enabled"] is True
        assert a["use_dnds"] is True
        assert a["use_gene_scores"] is True
        assert a["use_outlier_detection"] is True
        assert a["use_pairwise_ani"] is True
        assert a["pairwise_ani"]["pairs"] == "within_subject"
        assert a["use_strain_turnover"] is True
        # 4 mice per group: the real-cohort floor of 8 voters would leave every
        # verdict blank, so the test data votes with 2.
        assert a["strain_turnover"]["min_voters"] == 2
        assert a["use_baseline_presence"] is True
        bp = a["baseline_presence"]
        assert bp["baseline_scope"] == "cohort"
        assert bp["strain_gate"] == "carrier"
        assert bp["min_determined_carriers"] == 1
        assert a["use_regional_contrast"] is False
        assert a["use_gene_de_novo"] is True
        assert a["gene_de_novo"] == {
            "block_rule": "exclude",
            "min_gene_covered_fraction": 0.0,
            "min_determined_carriers": 1,
        }

    def test_single_keeps_the_longitudinal_only_features_off(self, tmp_path):
        a = self._emit(tmp_path, "single")["analysis"]
        assert a["use_lmm"] is False
        assert a["use_cmh"] is True
        assert a["figures"]["enabled"] is False
        assert a["use_gene_scores"] is True
        assert a["use_outlier_detection"] is True
        assert a["use_pairwise_ani"] is True
        assert a["pairwise_ani"]["pairs"] == "all"
        assert a["pairwise_ani"]["store_snp_locations"] == "none"
        for flag in (
            "use_dnds",
            "use_regional_contrast",
            "use_baseline_presence",
            "use_strain_turnover",
            "use_gene_de_novo",
        ):
            assert a[flag] is False, flag

    @pytest.mark.parametrize("data_type", ["longitudinal", "single"])
    def test_site_depth_filter_is_on(self, tmp_path, data_type):
        qc = self._emit(tmp_path, data_type)["quality_control"]
        assert qc["site_depth_filter"]["enabled"] is True
        assert qc["site_depth_filter"]["site_ratio"] == 3.0


# ---------------------------------------------------------------------------
# The allele "first seen after treatment began" and the witness timepoint
# ---------------------------------------------------------------------------
def _witness_profiles(longitudinal):
    """Profiles of the witness-timepoint samples of the ``longitudinal`` fixture."""
    mags, _, truth, _ = longitudinal
    return {
        w["sample_id"]: gen.generate_profile_data(w, mags, truth=truth)
        for w in truth["witness_timepoint"]["samples"]
    }


class TestFirstSeenAfterTreatment:
    """One listed site is carried, before the comparison's later timepoint, only
    by a sample taken AFTER treatment began (the witness timepoint).  The
    pipeline's three-label split must call it first_seen_after_treatment_began;
    here the generator's side of that fact is pinned."""

    def test_plan_turns_one_drawn_site_into_the_first_seen_site(self, longitudinal):
        _, _, truth, _ = longitudinal
        fact = truth["first_seen_after_treatment"]
        listed = [
            s
            for s in truth["significant_sites"]
            if s["kind"] == "first_seen_after_treatment"
        ]
        assert len(listed) == 1
        assert (fact["mag"], fact["contig"], fact["position"]) == (
            listed[0]["mag"],
            listed[0]["contig"],
            listed[0]["position"],
        )
        assert fact["rising_allele"] == listed[0]["rising_allele"]
        # It is one of the N_STANDING drawn sites, so the card's standing
        # count is one lower and the number of listed sites does not change.
        assert truth["n_standing"] == gen.N_STANDING - 1
        assert len(truth["significant_sites"]) == gen.N_STANDING + 4

    def test_witness_timepoint_is_outside_the_comparison_and_dated_between(
        self, longitudinal
    ):
        _, samples, truth, _ = longitudinal
        witness = truth["witness_timepoint"]
        assert witness["label"] not in truth["timepoints"]
        earlier, later = truth["timepoints"][0], truth["timepoints"][-1]
        dates = truth["collection_dates"]
        assert dates[earlier] < witness["collection_date"] < dates[later]
        assert dates[witness["label"]] == witness["collection_date"]
        assert witness["pre_treatment_timepoints"] == [earlier]
        # Two animals that are in no comparison: new subjects, a third group
        # and the witness timepoint.  The group must be outside the compared
        # pair: a null run permutes labels within the matched pairs of the two
        # compared groups, and an unpaired extra sample there stops it.
        assert len(witness["samples"]) == 2
        assert witness["group"] == gen.WITNESS_GROUP
        assert witness["group"] not in truth["groups"].values()
        compared_subjects = {s["subjectID"] for s in samples}
        compared_replicates = {s["replicate"] for s in samples}
        for w in witness["samples"]:
            assert w["group"] == witness["group"]
            assert w["time"] == witness["label"]
            assert w["subjectID"] not in compared_subjects
            assert w["replicate"] not in compared_replicates
            assert w["sample_id"] == (
                f"{w['group']}_{w['subjectID']}_{witness['label']}"
            )
        ids = [w["sample_id"] for w in witness["samples"]]
        assert truth["first_seen_after_treatment"]["carried_by"] == ids[0]

    def test_only_the_carriers_and_one_witness_carry_the_allele(self, longitudinal):
        _, samples, truth, profiles = longitudinal
        fact = truth["first_seen_after_treatment"]
        key = (fact["contig"], fact["position"])
        every = {**profiles, **_witness_profiles(longitudinal)}
        carrying = set(truth["carrier_samples"]) | {fact["carried_by"]}
        assert len(every) == len(samples) + 2
        for sample_id, by_mag in every.items():
            counts = _rows_by_position(by_mag[fact["mag"]])[key]
            rising = counts[fact["rising_allele"]]
            if sample_id in carrying:
                assert abs(rising / counts["cov"] - fact["fraction"]) < 0.03
                assert rising >= 5
            else:
                # every pre sample, the standing mouse and the second witness
                assert rising == 0, sample_id

    def test_witnesses_carry_no_other_planted_allele(self, longitudinal):
        """The de novo alleles must stay de novo and the standing ones must stay
        the standing mouse's alone: a witness with reads of them would change
        both labels."""
        _, _, truth, _ = longitudinal
        fact = truth["first_seen_after_treatment"]
        for sample_id, by_mag in _witness_profiles(longitudinal).items():
            for site in truth["significant_sites"]:
                if site["kind"] == "first_seen_after_treatment":
                    continue
                counts = _rows_by_position(by_mag[site["mag"]])[
                    (site["contig"], site["position"])
                ]
                assert counts[site["rising_allele"]] == 0, (sample_id, site)
            hot = truth["hot_gene"]
            for event in hot["events"]:
                counts = _rows_by_position(by_mag[hot["mag"]])[
                    (hot["contig"], event["position"])
                ]
                assert counts[event["rising_allele"]] == 0
        assert fact["carried_by"] in _witness_profiles(longitudinal)

    def test_single_design_has_no_witness_timepoint(self):
        random.seed(3)
        mags = [gen.create_mag_structure("MAG_001", num_contigs=2)]
        samples = gen.generate_sample_structure(
            8, timepoints=["t1"], data_type="single"
        )
        truth = gen.plan_the_plants(mags, samples, "single")
        assert "witness_timepoint" not in truth
        assert "first_seen_after_treatment" not in truth

    def test_a_timepoint_named_like_the_witness_label_is_refused(self):
        random.seed(3)
        mags = [
            gen.create_mag_structure(
                f"MAG_{i:03d}", num_contigs=2, contig_length_range=(1200, 1500)
            )
            for i in (1, 2)
        ]
        samples = gen.generate_sample_structure(
            16, timepoints=["pre", gen.WITNESS_TIMEPOINT], data_type="longitudinal"
        )
        with pytest.raises(ValueError, match=gen.WITNESS_TIMEPOINT):
            gen.plan_the_plants(mags, samples, "longitudinal")

    def test_a_group_named_like_the_witness_group_is_refused(self):
        random.seed(3)
        mags = [
            gen.create_mag_structure(
                f"MAG_{i:03d}", num_contigs=2, contig_length_range=(1200, 1500)
            )
            for i in (1, 2)
        ]
        samples = gen.generate_sample_structure(
            16, groups=["control", gen.WITNESS_GROUP], data_type="longitudinal"
        )
        with pytest.raises(ValueError, match=gen.WITNESS_GROUP):
            gen.plan_the_plants(mags, samples, "longitudinal")

    def test_null_run_permutation_leaves_the_witnesses_alone(self, tmp_path):
        """The tool that builds a null run swaps the two compared groups'
        labels within matched pairs; it must accept the bundled sheet and
        return the witness rows unchanged."""
        from alleleflux.scripts.preprocessing.permute_metadata import permute_metadata

        truth = gen.generate_dataset(
            tmp_path / "ds",
            data_type="longitudinal",
            num_samples=16,
            contig_length_range=(1200, 1500),
            seed=11,
        )
        sheet = pd.read_csv(
            tmp_path / "ds" / "metadata" / "sample_metadata.tsv", sep="\t", dtype=str
        )
        permuted = permute_metadata(
            sheet,
            unit="subjectID",
            block="replicate",
            seed=5,
            swaps=1,
            groups=[truth["groups"]["treatment"], truth["groups"]["control"]],
        )
        ids = [w["sample_id"] for w in truth["witness_timepoint"]["samples"]]
        is_witness = sheet["sample_id"].isin(ids)
        assert permuted[is_witness].equals(sheet[is_witness])
        assert (permuted["group"] != sheet["group"]).sum() == 4  # one pair, 2 x 2 rows

    def test_dataset_writes_the_witnesses_in_profile_mode_only(self, tmp_path):
        import yaml

        truth = gen.generate_dataset(
            tmp_path / "ds",
            data_type="longitudinal",
            num_samples=16,
            contig_length_range=(1200, 1500),
            seed=11,
        )
        witness = truth["witness_timepoint"]
        ids = [w["sample_id"] for w in witness["samples"]]
        meta = pd.read_csv(
            tmp_path / "ds" / "metadata" / "sample_metadata.tsv", sep="\t", dtype=str
        )
        assert len(meta) == 18
        rows = meta[meta["sample_id"].isin(ids)]
        assert list(rows["time"]) == [witness["label"]] * 2
        assert list(rows["group"]) == [gen.WITNESS_GROUP] * 2
        assert set(rows["collection_date"]) == {witness["collection_date"]}
        for sample_id in ids:
            folder = tmp_path / "ds" / "profiles" / sample_id
            assert len(list(folder.glob("*_profiled.tsv.gz"))) == 2  # one per MAG
        # The BAMs are simulated from the unmodified reference and carry no
        # plant, so the BAM-mode sheet keeps to the compared samples.
        bams = pd.read_csv(
            tmp_path / "ds" / "metadata" / "sample_metadata_bams.tsv",
            sep="\t",
            dtype=str,
        )
        assert len(bams) == 16
        assert not bams["sample_id"].isin(ids).any()
        cfg = yaml.safe_load(
            (tmp_path / "ds" / "config_generated_longitudinal.yml").read_text()
        )
        block = cfg["analysis"]["baseline_presence"]
        assert block["pre_treatment_timepoints"] == ["pre"]
