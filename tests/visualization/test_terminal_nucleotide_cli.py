"""End-to-end run of ``alleleflux-terminal-nucleotide`` on a tiny synthetic cache.

Why this file exists
--------------------
``test_terminal_nucleotide_analysis.py`` drives ``process_single_mag`` in
process, so argparse, the significance / test-type filters, the metadata ->
sample roster step, the subject pairing for the frequency-change method and the
summary table were only ever exercised by the visualization pipeline.  These
tests run the REAL console script as a subprocess on a hand-built
allele-frequency cache + metadata + significant-sites TSV, then open every
table it wrote and check the planted winners of all three methods:

1. ``terminal_nucleotide_mean_freq``     - base with the highest MEAN frequency
   across the terminal samples that have the site.
2. ``terminal_nucleotide_majority_vote`` - each terminal sample votes for its
   own major allele (tied bases split the vote); ties in the tally come back
   as a sorted comma-joined string, e.g. ``"A,T"``.
3. ``terminal_nucleotide_freq_change``   - per subject (same mouse at both
   timepoints) terminal minus initial frequency, averaged over the subjects
   that have the site at BOTH timepoints; the base with the largest positive
   mean change wins.  Only computed with ``--initial_timepoint``.

Fixture geometry
----------------
One MAG (``MAG_T``), one contig (``c1``), two timepoints (``pre`` -> ``post``).
Group ``g1`` has two mice (``m1``, ``m2``); group ``g2`` has one mouse (``m3``)
that is pure ``A`` at every site and must never be read: the CLI is run with
``--group g1``.  Counts per (sample, site), ``A10 G10`` = 10 A reads + 10 G reads:

    site  pos | m1_pre    m2_pre  | m1_post   m2_post | mean_freq  majority  freq_change (n mice)
    ----  --- | --------  ------- | --------  ------- | ---------  --------  --------------------
    (a)    10 | A10 G10   A12 G8  | A8 G12    A8 G12  |    G          G         G   (2)
    (b)    20 | A20       A20     | A20       A8 T12  |    A         A,T        T   (2)
    (c)    30 | G20       G20     | A8 G12    A8 G12  |    G          G         A   (2)
    (d)    40 | C20       (none)  | C15 T5    C20     |    C          C         T   (1)
    decoy  50 | A20       A20     | A20       A20     |  q_value 0.5 > 0.05: filtered out

* (a) all three agree: G is 0.6 in both post samples, rose from 0.5 / 0.4.
  If the group filter leaked and m3_post (A 1.0) were pooled in, the mean A
  would become (0.4 + 0.4 + 1.0) / 3 = 0.6 > G 0.4 and flip the answer.
* (b) mean vs vote: m1_post is 100 % A, m2_post is 60 % T.  Mean: A 0.7,
  T 0.3 -> ``A``.  Votes: one for A, one for T -> tie ``"A,T"`` (with two
  samples the vote can only agree or tie; this is the code's tie output).
  Change: m1 0.0 for every base, m2 T +0.6 / A -0.6 -> mean T +0.3 -> ``T``.
* (c) change contradicts the terminal majority: G is major at post (0.6) but
  FELL from 1.0, while A rose 0.0 -> 0.4 in both mice -> change picks ``A``.
* (d) m2 has no ``pre`` row at pos 40, so the change method uses m1 only:
  C 1.0 -> 0.75 (-0.25), T 0.0 -> 0.25 (+0.25) -> ``T`` with
  ``n_subjects_used_for_freq_change`` = 1.  The terminal methods still see
  both post samples: mean C (0.75 + 1.0) / 2 = 0.875 -> ``C``, votes C, C.

Two decoy rows in the significant-sites TSV must not reach the output: the
``MAG_T`` pos 50 row above the q-value threshold, and a ``MAG_Z`` row whose
``test_type`` is ``CMH`` while the CLI asks for ``two_sample_paired_tTest``
(a tiny ``MAG_Z`` cache exists so a leak shows up as an unexpected output
directory rather than as a crash).
"""

import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from tests.utilities.test_cache_io import write_cache

# --- Fixture constants --------------------------------------------------------
MAG = "MAG_T"
DECOY_MAG = "MAG_Z"
CONTIG = "c1"
GROUP = "g1"
OTHER_GROUP = "g2"
PRE, POST = "pre", "post"
TEST_TYPE = "two_sample_paired_tTest"
Q_THRESHOLD = 0.05

# Mouse -> (group, sample id at pre, sample id at post).  Sample ids are
# arbitrary labels; they only have to match between metadata and cache.
MICE = {
    "m1": (GROUP, "m1_pre", "m1_post"),
    "m2": (GROUP, "m2_pre", "m2_post"),
    "m3": (OTHER_GROUP, "m3_pre", "m3_post"),
}
TERMINAL_SAMPLES = ["m1_post", "m2_post"]  # metadata order == column order

# Planted positions (see the module docstring table) and the decoy.
POS_AGREE, POS_MEAN_VS_VOTE, POS_CHANGE_VS_VOTE, POS_MISSING_PRE = 10, 20, 30, 40
POS_DECOY = 50
PLANTED = [POS_AGREE, POS_MEAN_VS_VOTE, POS_CHANGE_VS_VOTE, POS_MISSING_PRE]

# Counts as (A, C, G, T) per (sample, position).  A position absent from this
# table for a sample means "no row in the cache" (the site-depth filter dropped
# it, or the sample never covered it).
COUNTS = {
    # (a) agreement
    ("m1_pre", POS_AGREE): (10, 0, 10, 0),
    ("m2_pre", POS_AGREE): (12, 0, 8, 0),
    ("m1_post", POS_AGREE): (8, 0, 12, 0),
    ("m2_post", POS_AGREE): (8, 0, 12, 0),
    # (b) mean vs majority vote (and change picks a third answer)
    ("m1_pre", POS_MEAN_VS_VOTE): (20, 0, 0, 0),
    ("m2_pre", POS_MEAN_VS_VOTE): (20, 0, 0, 0),
    ("m1_post", POS_MEAN_VS_VOTE): (20, 0, 0, 0),
    ("m2_post", POS_MEAN_VS_VOTE): (8, 0, 0, 12),
    # (c) change contradicts the terminal majority
    ("m1_pre", POS_CHANGE_VS_VOTE): (0, 0, 20, 0),
    ("m2_pre", POS_CHANGE_VS_VOTE): (0, 0, 20, 0),
    ("m1_post", POS_CHANGE_VS_VOTE): (8, 0, 12, 0),
    ("m2_post", POS_CHANGE_VS_VOTE): (8, 0, 12, 0),
    # (d) m2 has no pre row here -> change uses m1 only
    ("m1_pre", POS_MISSING_PRE): (0, 20, 0, 0),
    ("m1_post", POS_MISSING_PRE): (0, 15, 0, 5),
    ("m2_post", POS_MISSING_PRE): (0, 20, 0, 0),
    # decoy site: covered everywhere, but its q_value is above the threshold
    ("m1_pre", POS_DECOY): (20, 0, 0, 0),
    ("m2_pre", POS_DECOY): (20, 0, 0, 0),
    ("m1_post", POS_DECOY): (20, 0, 0, 0),
    ("m2_post", POS_DECOY): (20, 0, 0, 0),
}
# The other-group mouse is pure A everywhere at both timepoints: pooling it in
# would flip site (a)'s mean-frequency winner from G to A.
for _pos in PLANTED + [POS_DECOY]:
    COUNTS[("m3_pre", _pos)] = (20, 0, 0, 0)
    COUNTS[("m3_post", _pos)] = (20, 0, 0, 0)

# Expected winners per method, keyed by position.
EXPECTED_MEAN = {
    POS_AGREE: "G",
    POS_MEAN_VS_VOTE: "A",
    POS_CHANGE_VS_VOTE: "G",
    POS_MISSING_PRE: "C",
}
EXPECTED_VOTE = {
    POS_AGREE: "G",
    POS_MEAN_VS_VOTE: "A,T",
    POS_CHANGE_VS_VOTE: "G",
    POS_MISSING_PRE: "C",
}
EXPECTED_CHANGE = {
    POS_AGREE: "G",
    POS_MEAN_VS_VOTE: "T",
    POS_CHANGE_VS_VOTE: "A",
    POS_MISSING_PRE: "T",
}
EXPECTED_N_SUBJECTS = {
    POS_AGREE: 2,
    POS_MEAN_VS_VOTE: 2,
    POS_CHANGE_VS_VOTE: 2,
    POS_MISSING_PRE: 1,
}


# --- Fixture builders ---------------------------------------------------------
def _cache_rows(timepoint):
    """``write_cache`` rows (contig, pos, sample, group, subject, A, C, G, T) for one timepoint.

    Example: ``_cache_rows("post")`` yields ``("c1", 10, "m1_post", "g1", "m1", 8, 0, 12, 0)``
    among others; m2 contributes no pos-40 row at ``pre``.
    """
    rows = []
    for mouse, (group, s_pre, s_post) in MICE.items():
        sample = s_pre if timepoint == PRE else s_post
        for (sid, pos), (a, c, g, t) in COUNTS.items():
            if sid == sample:
                rows.append((CONTIG, pos, sample, group, mouse, a, c, g, t))
    return rows


def _write_metadata(path):
    """Metadata TSV with the columns the CLI needs: sample_id, subjectID, group, time."""
    recs = []
    for mouse, (group, s_pre, s_post) in MICE.items():
        recs.append(dict(sample_id=s_pre, subjectID=mouse, group=group, time=PRE))
        recs.append(dict(sample_id=s_post, subjectID=mouse, group=group, time=POST))
    pd.DataFrame(recs).to_csv(path, sep="\t", index=False)


def _write_sig_sites(path):
    """Significant-sites TSV in the p_value_summary layout (period, mag_id, contig,
    position, gene_id, test_type, min_p_value, source_file, q_value) with the four
    planted sites plus the two decoys."""
    recs = [
        dict(
            period=f"{PRE}_{POST}",
            mag_id=MAG,
            contig=CONTIG,
            position=pos,
            gene_id=f"{MAG}_{CONTIG}_gene1",
            test_type=TEST_TYPE,
            min_p_value=0.001,
            source_file=f"{MAG}_two_sample_paired.tsv.gz",
            q_value=0.01,
        )
        for pos in PLANTED
    ]
    # Decoy 1: same MAG, significant p but q_value above the threshold.
    recs.append(
        dict(
            period=f"{PRE}_{POST}",
            mag_id=MAG,
            contig=CONTIG,
            position=POS_DECOY,
            gene_id=f"{MAG}_{CONTIG}_gene1",
            test_type=TEST_TYPE,
            min_p_value=0.2,
            source_file=f"{MAG}_two_sample_paired.tsv.gz",
            q_value=0.5,
        )
    )
    # Decoy 2: another MAG, significant, but a test_type the CLI is not asked for.
    recs.append(
        dict(
            period=f"{PRE}_{POST}",
            mag_id=DECOY_MAG,
            contig="c9",
            position=5,
            gene_id=f"{DECOY_MAG}_c9_gene1",
            test_type="CMH",
            min_p_value=0.001,
            source_file=f"{DECOY_MAG}_cmh.tsv.gz",
            q_value=0.001,
        )
    )
    pd.DataFrame(recs).to_csv(path, sep="\t", index=False)


def _cli():
    """Command prefix for the console script.

    Prefer the entry point next to ``sys.executable`` so a stale copy earlier on
    PATH (another interpreter) can never be picked up by accident.
    """
    exe = Path(sys.executable).parent / "alleleflux-terminal-nucleotide"
    if exe.exists():
        return [str(exe)]
    return [
        sys.executable,
        "-m",
        "alleleflux.scripts.visualization.accessory.terminal_nucleotide_analysis",
    ]


def _run_cli(tmp_path, with_initial):
    """Build the fixture under ``tmp_path``, run the CLI, return (CompletedProcess, outdir)."""
    cache_dir = tmp_path / "allele_freq_cache"
    for tp in (PRE, POST):
        write_cache(str(cache_dir), MAG, tp, _cache_rows(tp))
        # A tiny MAG_Z cache so a leaked decoy row would produce an output
        # directory (asserted absent) instead of a FileNotFoundError.
        write_cache(
            str(cache_dir),
            DECOY_MAG,
            tp,
            [("c9", 5, f"m1_{tp}", GROUP, "m1", 20, 0, 0, 0)],
        )
    metadata = tmp_path / "metadata.tsv"
    _write_metadata(metadata)
    sig = tmp_path / "p_value_summary.tsv"
    _write_sig_sites(sig)
    outdir = tmp_path / "terminal_out"
    cmd = _cli() + [
        "--significant_sites", str(sig),
        "--cache_dir", str(cache_dir),
        "--group", GROUP,
        "--timepoint", POST,
        "--metadata", str(metadata),
        "--output", str(outdir),
        "--p_value_column", "q_value",
        "--p_value_threshold", str(Q_THRESHOLD),
        "--test-type", TEST_TYPE,
        "--cpus", "1",
    ]  # fmt: skip
    if with_initial:
        cmd += ["--initial_timepoint", PRE]
    done = subprocess.run(cmd, capture_output=True, text=True)
    return done, outdir


def _terminal_table(outdir):
    """``<out>/MAG_T/MAG_T_terminal_nucleotides.tsv`` indexed by position."""
    path = outdir / MAG / f"{MAG}_terminal_nucleotides.tsv"
    assert path.exists(), f"missing {path}; outdir has {sorted(outdir.rglob('*'))}"
    return pd.read_csv(path, sep="\t").set_index("position")


def _frequencies_table(outdir):
    """``<out>/MAG_T/MAG_T_nucleotide_frequencies.tsv.gz`` indexed by position
    (opened through gzip explicitly so the compression is part of what is tested)."""
    path = outdir / MAG / f"{MAG}_nucleotide_frequencies.tsv.gz"
    assert path.exists(), f"missing {path}; outdir has {sorted(outdir.rglob('*'))}"
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, sep="\t").set_index("position")


def _summary(outdir):
    """``<out>/terminal_nucleotide_analysis_summary.tsv``, one row per MAG."""
    path = outdir / "terminal_nucleotide_analysis_summary.tsv"
    assert path.exists(), f"missing {path}; outdir has {sorted(outdir.rglob('*'))}"
    return pd.read_csv(path, sep="\t")


# --- Tests --------------------------------------------------------------------
class TestTerminalNucleotideCli:
    """The console script on the synthetic cache in the module docstring."""

    @pytest.fixture(scope="class")
    def with_initial(self, tmp_path_factory):
        """One run WITH ``--initial_timepoint`` (all three methods), shared by the class."""
        done, outdir = _run_cli(tmp_path_factory.mktemp("with_initial"), True)
        assert done.returncode == 0, f"stdout:\n{done.stdout}\nstderr:\n{done.stderr}"
        return outdir

    @pytest.fixture(scope="class")
    def without_initial(self, tmp_path_factory):
        """One run WITHOUT ``--initial_timepoint`` (mean + vote only), shared by the class."""
        done, outdir = _run_cli(tmp_path_factory.mktemp("without_initial"), False)
        assert done.returncode == 0, f"stdout:\n{done.stdout}\nstderr:\n{done.stderr}"
        return outdir

    def test_only_the_planted_sites_survive_the_filters(self, with_initial):
        """The q-value decoy (pos 50) and the MAG_Z CMH row are filtered before
        processing: exactly the four planted positions, one MAG directory,
        and ``source_file`` dropped from the per-MAG table."""
        table = _terminal_table(with_initial)
        assert sorted(table.index) == PLANTED
        assert POS_DECOY not in table.index
        assert "source_file" not in table.columns
        assert not (with_initial / DECOY_MAG).exists()
        # The sig-sites columns ride along (period, mag_id, gene_id, q_value ...).
        assert set(table["mag_id"]) == {MAG}
        assert (table["q_value"] <= Q_THRESHOLD).all()

    def test_three_methods_pick_the_planted_winners(self, with_initial):
        """Each planted site isolates one disagreement mode (docstring table).

        Worked example, site (c) pos 30: post counts A8 G12 in both mice ->
        G 0.6 wins the mean and both votes; but pre was G20 in both, so the
        per-mouse change is A +0.4 / G -0.4 -> the change method says ``A``.
        """
        table = _terminal_table(with_initial)
        for pos in PLANTED:
            row = table.loc[pos]
            assert row["terminal_nucleotide_mean_freq"] == EXPECTED_MEAN[pos], pos
            assert row["terminal_nucleotide_majority_vote"] == EXPECTED_VOTE[pos], pos
            assert row["terminal_nucleotide_freq_change"] == EXPECTED_CHANGE[pos], pos

    def test_per_sample_frequencies_and_mice_used(self, with_initial):
        """The frequencies table carries the raw per-sample fractions, the means,
        the per-sample major alleles, the mean changes and the n-mice counts.

        Worked example, site (d) pos 40: m1_post C15 T5 -> C 0.75 / T 0.25;
        m2_post C20 -> C 1.0; mean C 0.875.  m2 has no pre row, so the change
        is m1 alone: C 0.75 - 1.0 = -0.25, T 0.25 - 0.0 = +0.25, n = 1.
        """
        freqs = _frequencies_table(with_initial)
        assert sorted(freqs.index) == PLANTED

        # Group filter: only g1's post samples get columns; m3 never appears.
        for nuc in "ACGT":
            assert f"{nuc}_m3_post" not in freqs.columns
            assert (
                f"{nuc}_m1_post" in freqs.columns and f"{nuc}_m2_post" in freqs.columns
            )
        assert (freqs["n_samples_used_for_mean"] == 2).all()

        # (a) pos 10: A8 G12 in both post samples.
        a = freqs.loc[POS_AGREE]
        assert a["A_m1_post"] == pytest.approx(0.4) and a["G_m1_post"] == pytest.approx(
            0.6
        )
        assert a["A_m2_post"] == pytest.approx(0.4) and a["G_m2_post"] == pytest.approx(
            0.6
        )
        assert a["G_mean_frequency"] == pytest.approx(0.6)
        assert a["major_allele_m1_post"] == "G" and a["major_allele_m2_post"] == "G"
        # change: m1 G 0.6 - 0.5 = +0.1, m2 G 0.6 - 0.4 = +0.2 -> mean +0.15
        assert a["G_frequency_change"] == pytest.approx(0.15)
        assert a["A_frequency_change"] == pytest.approx(-0.15)

        # (b) pos 20: m1_post pure A, m2_post 60 % T.
        b = freqs.loc[POS_MEAN_VS_VOTE]
        assert b["A_m1_post"] == pytest.approx(1.0) and b["T_m1_post"] == pytest.approx(
            0.0
        )
        assert b["A_m2_post"] == pytest.approx(0.4) and b["T_m2_post"] == pytest.approx(
            0.6
        )
        assert b["A_mean_frequency"] == pytest.approx(0.7)
        assert b["T_mean_frequency"] == pytest.approx(0.3)
        assert b["major_allele_m1_post"] == "A" and b["major_allele_m2_post"] == "T"
        # change: m1 0 everywhere; m2 T +0.6, A -0.6 -> means +0.3 / -0.3
        assert b["T_frequency_change"] == pytest.approx(0.3)
        assert b["A_frequency_change"] == pytest.approx(-0.3)

        # (c) pos 30: G major at post, but A is what rose (0.0 -> 0.4 in both).
        c = freqs.loc[POS_CHANGE_VS_VOTE]
        assert c["G_mean_frequency"] == pytest.approx(0.6)
        assert c["A_frequency_change"] == pytest.approx(0.4)
        assert c["G_frequency_change"] == pytest.approx(-0.4)

        # (d) pos 40: m2 lacks a pre row -> one mouse used for the change.
        d = freqs.loc[POS_MISSING_PRE]
        assert d["C_m1_post"] == pytest.approx(0.75) and d[
            "T_m1_post"
        ] == pytest.approx(0.25)
        assert d["C_m2_post"] == pytest.approx(1.0)
        assert d["C_mean_frequency"] == pytest.approx(0.875)
        assert d["T_frequency_change"] == pytest.approx(0.25)
        assert d["C_frequency_change"] == pytest.approx(-0.25)

        # n mice used by the change method, per site.
        for pos in PLANTED:
            assert (
                freqs.loc[pos, "n_subjects_used_for_freq_change"]
                == EXPECTED_N_SUBJECTS[pos]
            ), pos

        # The winners are also written into the frequencies table.
        for pos in PLANTED:
            assert freqs.loc[pos, "terminal_nucleotide_mean_freq"] == EXPECTED_MEAN[pos]
            assert (
                freqs.loc[pos, "terminal_nucleotide_majority_vote"]
                == EXPECTED_VOTE[pos]
            )
            assert (
                freqs.loc[pos, "terminal_nucleotide_freq_change"]
                == EXPECTED_CHANGE[pos]
            )

    def test_summary_names_the_mag_with_its_site_and_winner_counts(self, with_initial):
        """One summary row for MAG_T: 4 sites, 2 terminal samples, 2 paired mice,
        and per-method tallies of the winners (a tie such as ``"A,T"`` counts
        once for each base)."""
        summary = _summary(with_initial)
        assert summary["mag_id"].tolist() == [MAG]
        row = summary.iloc[0]
        assert row["sites_processed"] == 4
        assert row["samples_processed"] == 2
        assert row["subjects_used"] == 2
        assert Path(row["main_output"]).exists()
        assert Path(row["intermediate_output"]).exists()
        # mean: G (a), A (b), G (c), C (d)
        assert row["mean_freq_terminal_G_count"] == 2
        assert row["mean_freq_terminal_A_count"] == 1
        assert row["mean_freq_terminal_C_count"] == 1
        # vote: G, "A,T", G, C -> the tie adds one to A and one to T
        assert row["majority_vote_terminal_G_count"] == 2
        assert row["majority_vote_terminal_A_count"] == 1
        assert row["majority_vote_terminal_T_count"] == 1
        assert row["majority_vote_terminal_C_count"] == 1
        # change: G, T, A, T
        assert row["freq_change_terminal_G_count"] == 1
        assert row["freq_change_terminal_T_count"] == 2
        assert row["freq_change_terminal_A_count"] == 1

    def test_without_initial_timepoint_only_two_methods_are_written(
        self, without_initial
    ):
        """No ``--initial_timepoint``: the mean and vote answers are unchanged and
        no frequency-change column exists in either table or in the summary."""
        table = _terminal_table(without_initial)
        assert sorted(table.index) == PLANTED
        assert "terminal_nucleotide_freq_change" not in table.columns
        for pos in PLANTED:
            assert (
                table.loc[pos, "terminal_nucleotide_mean_freq"] == EXPECTED_MEAN[pos]
            ), pos
            assert (
                table.loc[pos, "terminal_nucleotide_majority_vote"]
                == EXPECTED_VOTE[pos]
            ), pos

        freqs = _frequencies_table(without_initial)
        assert not [c for c in freqs.columns if c.endswith("_frequency_change")]
        assert "n_subjects_used_for_freq_change" not in freqs.columns
        assert "terminal_nucleotide_freq_change" not in freqs.columns
        # The per-sample fractions do not depend on the initial timepoint.
        assert freqs.loc[POS_MEAN_VS_VOTE, "T_m2_post"] == pytest.approx(0.6)
        assert freqs.loc[POS_MISSING_PRE, "C_mean_frequency"] == pytest.approx(0.875)

        summary = _summary(without_initial)
        assert summary["mag_id"].tolist() == [MAG]
        assert summary.iloc[0]["sites_processed"] == 4
        assert "subjects_used" not in summary.columns
        assert not [c for c in summary.columns if c.startswith("freq_change_terminal_")]
