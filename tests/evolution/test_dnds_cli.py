"""End-to-end run of ``alleleflux-dnds-from-timepoints`` on a tiny synthetic genome.

Why this file exists
--------------------
``test_dnds_cache.py`` only checks ``load_sample_counts``; the rest of the dN/dS
path (ancestral reconstruction -> Prodigal-vs-cache ``ref_base`` sanity check ->
NG86 codon events -> gene / MAG / global tables) was only ever driven by the full
pipeline regression test, which is how a ``KeyError: 'ref_base'`` in the cache
reader slipped past the unit suite.  These tests drive the console script as a
subprocess on a hand-built cache + Prodigal FASTA + significant-sites TSV, then
open the tables it wrote and check the planted answers.

Fixture geometry (one contig ``c1``, 27 bp, 0-based positions 0..26)
--------------------------------------------------------------------
::

    contig pos : 0  1  2  3  4  5  6  7  8  9 10 11 | 12 13 14 | 15 ... 26
    fwd base   : A  T  G  A  A  A  G  G  C  T  A  A |  C  C  C | T C A G A A C A G C A T
                 [---------- g1, + strand ---------]  intergenic  [------ g2, - strand ------]

* ``g1`` = Prodigal ``# 1 # 12 # 1``  (1-based inclusive), coding seq ``ATG AAA GGC TAA``
  (M K G *).  Gene index == contig position because the gene starts at 1-based 1.
* ``g2`` = Prodigal ``# 16 # 27 # -1``; coding seq ``ATG CTG TTC TGA`` (M L F *), which
  is the reverse complement of forward bases 15..26 (``TCAGAACAGCAT``).
  Gene index = (end - 1) - contig_pos = 26 - contig_pos, and every allele that
  enters the coding sequence is the COMPLEMENT of the forward-strand base.

Planted substitutions (ancestral sample fixed for the reference base, derived
sample carries the new base at 17 of 20 reads):

    pos  gene  strand  ref(fwd)  derived(fwd)  coding codon change    class
    ---  ----  ------  --------  ------------  --------------------   -----
      5   g1     +        A           G         AAA -> AAG  (K -> K)    S
      6   g1     +        G           A         GGC -> AGC  (G -> S)    NS
     18   g2     -        G           A         TTC -> TTT  (F -> F)    S
     23   g2     -        G           T         CTG -> ATG  (L -> M)    NS

The reverse-strand rows are the ones that matter most: at contig 18 the forward
reference is ``G``, Prodigal's coding base is ``C`` (gene index 8), and the
sanity check must complement Prodigal's ``C`` to ``G`` before comparing, or
every reverse-strand site would log a spurious mismatch.
"""

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from alleleflux.scripts.utilities.cache_io import cache_path

# --- The synthetic genome -----------------------------------------------------
MAG = "MAG_X"
CONTIG = "c1"
G1, G2 = f"{MAG}_{CONTIG}_g1", f"{MAG}_{CONTIG}_g2"
G1_CODING = "ATGAAAGGCTAA"  # + strand, contig 0..11
G2_CODING = "ATGCTGTTCTGA"  # - strand, contig 15..26 (coding = revcomp of forward)
INTERGENIC = "CCC"  # contig 12..14, no gene_id -> must be dropped, not crash
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _revcomp(seq):
    """Reverse complement, e.g. ``ATGCTGTTCTGA`` -> ``TCAGAACAGCAT``."""
    return "".join(COMPLEMENT[b] for b in reversed(seq))


# Forward-strand reference of the whole contig: g1 + intergenic + revcomp(g2).
CONTIG_FWD = G1_CODING + INTERGENIC + _revcomp(G2_CODING)
assert len(CONTIG_FWD) == 27

# contig position -> gene_id (None for the intergenic stretch)
GENE_AT = {p: G1 for p in range(0, 12)}
GENE_AT.update({p: None for p in range(12, 15)})
GENE_AT.update({p: G2 for p in range(15, 27)})

# The planted derived alleles, forward strand (see module docstring table).
PLANTED = {5: "G", 6: "A", 18: "A", 23: "T"}
EXPECTED_CLASS = {5: "S", 6: "NS", 18: "S", 23: "NS"}

ANC_SAMPLE, DER_SAMPLE = "s1_pre", "s1_post"
ANC_TP, DER_TP = "pre", "post"
TEST_TYPE = "two_sample_paired_tTest"
DEPTH = 20
MINOR = 3  # reads of the old base left in the derived sample (17 vs 3, no tie)


# --- Fixture builders ---------------------------------------------------------
def _cache_frame(sample_id, timepoint, derived, ref_base_override=None):
    """One sample's rows in the allele-frequency cache schema.

    ``derived=False`` -> every position is fixed for the reference base.
    ``derived=True``  -> the four PLANTED positions carry the new base at
    ``DEPTH - MINOR`` reads, the reference base at ``MINOR``.
    ``ref_base_override`` = {pos: base} lets a test write a WRONG reference base
    into the cache to provoke the sanity-check warning.
    """
    recs = []
    for pos, ref in enumerate(CONTIG_FWD):
        counts = {b: 0 for b in "ACGT"}
        if derived and pos in PLANTED:
            counts[PLANTED[pos]] = DEPTH - MINOR
            counts[ref] = MINOR
        else:
            counts[ref] = DEPTH
        stored_ref = (ref_base_override or {}).get(pos, ref)
        recs.append(
            dict(
                contig=CONTIG,
                position=pos,
                ref_base=stored_ref,
                gene_id=GENE_AT[pos],
                total_coverage=DEPTH,
                A=counts["A"],
                C=counts["C"],
                G=counts["G"],
                T=counts["T"],
                sample_id=sample_id,
                group="treatment",
                subjectID="m1",
                replicate="m1",
                A_frequency=counts["A"] / DEPTH,
                T_frequency=counts["T"] / DEPTH,
                G_frequency=counts["G"] / DEPTH,
                C_frequency=counts["C"] / DEPTH,
                time=timepoint,
            )
        )
    df = pd.DataFrame(recs)
    # The cache writer stores its label columns dictionary-encoded; ref_base too.
    for col in ("sample_id", "group", "subjectID", "replicate", "time", "ref_base"):
        df[col] = df[col].astype("category")
    return df


def _write_cache(cache_dir, ref_base_override=None):
    """``{cache_dir}/pre/MAG_X_pre_allele_frequency.parquet`` and the ``post`` twin."""
    for sample, tp, derived in (
        (ANC_SAMPLE, ANC_TP, False),
        (DER_SAMPLE, DER_TP, True),
    ):
        path = cache_path(str(cache_dir), MAG, tp)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        _cache_frame(sample, tp, derived, ref_base_override).to_parquet(
            path, index=False
        )


def _write_prodigal(path):
    """Prodigal-style ``>id # start # end # strand # ...`` headers, 1-based inclusive."""
    path.write_text(
        f">{G1} # 1 # 12 # 1 # ID=1_1;partial=00;start_type=ATG\n{G1_CODING}\n"
        f">{G2} # 16 # 27 # -1 # ID=1_2;partial=00;start_type=ATG\n{G2_CODING}\n"
    )


def _write_sig_sites(path):
    """A p_value_summary-shaped TSV: the four planted sites plus three decoys.

    Decoys that the CLI must drop: a non-significant row (q = 0.8), a row of
    another MAG, and a row of another test type.  One gene_id carries a trailing
    space, as real summaries do; the CLI strips it.
    """
    rows = [
        dict(mag_id=MAG, position=p, gene_id=GENE_AT[p], test_type=TEST_TYPE, q=0.01)
        for p in PLANTED
    ]
    rows[0]["gene_id"] = rows[0]["gene_id"] + " "
    rows += [
        dict(mag_id=MAG, position=0, gene_id=G1, test_type=TEST_TYPE, q=0.8),
        dict(mag_id="MAG_Y", position=5, gene_id=G1, test_type=TEST_TYPE, q=0.01),
        dict(
            mag_id=MAG,
            position=9,
            gene_id=G1,
            test_type="two_sample_paired_Wilcoxon",
            q=0.01,
        ),
    ]
    df = pd.DataFrame(rows)
    pd.DataFrame(
        {
            "period": "pre_post",
            "mag_id": df["mag_id"],
            "contig": CONTIG,
            "position": df["position"],
            "gene_id": df["gene_id"],
            "test_type": df["test_type"],
            "min_p_value": df["q"] / 2,
            "source_file": f"{MAG}_two_sample_paired.tsv.gz",
            "q_value": df["q"],
        }
    ).to_csv(path, sep="\t", index=False)


def _cli():
    """The console script from the SAME env as pytest, else ``python -m`` fallback.

    Prefer the entry point next to ``sys.executable`` so a stale copy earlier on
    PATH (another interpreter) can never be picked up by accident.
    """
    exe = Path(sys.executable).parent / "alleleflux-dnds-from-timepoints"
    if exe.exists():
        return [str(exe)]
    return [sys.executable, "-m", "alleleflux.scripts.evolution.dnds_from_timepoints"]


def _run_dnds(tmp_path, ref_base_override=None):
    """Build the fixture under ``tmp_path``, run the CLI, return (CompletedProcess, outdir)."""
    cache_dir = tmp_path / "allele_freq_cache"
    _write_cache(cache_dir, ref_base_override)
    prodigal = tmp_path / "genes.fna"
    _write_prodigal(prodigal)
    sig = tmp_path / "p_value_summary.tsv"
    _write_sig_sites(sig)
    outdir = tmp_path / "dnds_out"
    cmd = _cli() + [
        "--significant_sites", str(sig),
        "--mag_ids", MAG,
        "--p_value_column", "q_value",
        "--p_value_threshold", "0.05",
        "--test-type", TEST_TYPE,
        "--ancestral_sample_id", ANC_SAMPLE,
        "--derived_sample_id", DER_SAMPLE,
        "--cache_dir", str(cache_dir),
        "--ancestral_timepoint", ANC_TP,
        "--derived_timepoint", DER_TP,
        "--prodigal_fasta", str(prodigal),
        "--outdir", str(outdir),
        "--prefix", "dnds",
        "--cpus", "1",
    ]  # fmt: skip
    done = subprocess.run(cmd, capture_output=True, text=True)
    return done, outdir


def _table(outdir, suffix):
    """Read ``dnds_out/MAG_X/dnds_MAG_X_<test>_q_value_<suffix>`` (process_mag's naming)."""
    path = outdir / MAG / f"dnds_{MAG}_{TEST_TYPE}_q_value_{suffix}"
    assert (
        path.exists()
    ), f"missing output {path}; dir has {sorted(os.listdir(outdir / MAG))}"
    return pd.read_csv(path, sep="\t")


# --- Tests --------------------------------------------------------------------
class TestDndsCli:
    """The console script on the synthetic genome in the module docstring."""

    def test_planted_substitutions_are_classified_and_summarised(self, tmp_path):
        """Four planted k=1 events: S/NS per the table; both genes summarised.

        Worked example for the reverse-strand NS row: contig 23 is forward
        ``G`` ancestrally and ``T`` derived; gene index 26 - 23 = 3 is the
        first base of codon 2 of g2; complemented, C -> A turns CTG (L) into
        ATG (M): non-synonymous, frac_NS = 1, frac_S = 0.
        """
        done, outdir = _run_dnds(tmp_path)
        assert done.returncode == 0, done.stderr

        # --- codon events: one row per planted site, k=1, classified as planted
        events = _table(outdir, "codon_events_ng86.tsv")
        assert len(events) == 4, events
        by_pos = events.set_index("contig_position")
        by_pos.index = by_pos.index.astype(int)  # k=1 -> a single position per row
        assert sorted(by_pos.index) == sorted(PLANTED)
        for pos, cls in EXPECTED_CLASS.items():
            row = by_pos.loc[pos]
            assert row["mutation_type"] == cls, (pos, row.to_dict())
            assert row["k"] == 1 and bool(row["is_valid"])
            # fractional NG86 counts are whole for k=1: all on the planted side
            assert (row["frac_S"], row["frac_NS"]) == (
                (1.0, 0.0) if cls == "S" else (0.0, 1.0)
            )
            assert row["ancestral_allele"] == CONTIG_FWD[pos]  # forward-strand ref
            assert row["derived_allele"] == PLANTED[pos]
        # the strand reached the table, and the reverse gene's codons were read
        # from its coding sequence (not the forward contig)
        assert by_pos.loc[23, "strand"] == -1 and by_pos.loc[5, "strand"] == 1
        assert (by_pos.loc[23, "ancestral_codon"], by_pos.loc[23, "derived_codon"]) == (
            "CTG",
            "ATG",
        )
        assert (by_pos.loc[18, "ancestral_codon"], by_pos.loc[18, "derived_codon"]) == (
            "TTC",
            "TTT",
        )
        assert (by_pos.loc[5, "ancestral_codon"], by_pos.loc[5, "derived_codon"]) == (
            "AAA",
            "AAG",
        )
        assert (by_pos.loc[6, "ancestral_codon"], by_pos.loc[6, "derived_codon"]) == (
            "GGC",
            "AGC",
        )

        # --- gene summary: both genes, one S + one NS each, finite dN and dS
        genes = _table(outdir, "gene_summary_ng86.tsv").set_index("gene_id")
        assert sorted(genes.index) == sorted([G1, G2])
        for g in (G1, G2):
            assert genes.loc[g, "num_codon_events"] == 2
            assert (genes.loc[g, "total_frac_s"], genes.loc[g, "total_frac_ns"]) == (
                1.0,
                1.0,
            )
            assert genes.loc[g, "dN"] > 0 and genes.loc[g, "dS"] > 0
            assert genes.loc[g, "dN_dS_ratio"] > 0

        # --- MAG summary: one MAG, two genes, sums of the above
        mag = _table(outdir, "mag_summary_ng86.tsv")
        assert len(mag) == 1 and mag.loc[0, "mag_id"] == MAG
        assert mag.loc[0, "num_genes"] == 2
        assert (mag.loc[0, "total_frac_s"], mag.loc[0, "total_frac_ns"]) == (2.0, 2.0)

        # --- global summary exists and is non-empty
        assert len(_table(outdir, "global_dnds_ng86_summary.tsv")) > 0

        # --- ancestral ORFs: ancestral sample == reference, so the reconstructed
        # sequences must round-trip to the Prodigal coding sequences (this is
        # where a strand/complement slip in the reconstruction would show)
        ffn = (
            outdir / MAG / f"dnds_{MAG}_{TEST_TYPE}_q_value_ancestral_orfs.ffn"
        ).read_text()
        seqs, cur = {}, None
        for line in ffn.splitlines():
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = ""
            else:
                seqs[cur] += line.strip()
        assert seqs == {f"ancestral_{G1}": G1_CODING, f"ancestral_{G2}": G2_CODING}

        # no reference-base mismatch should be reported on a consistent fixture
        assert "Reference base mismatch" not in done.stderr, done.stderr

    def test_wrong_ref_base_in_cache_is_warned_not_fatal(self, tmp_path):
        """A cache ``ref_base`` that contradicts Prodigal logs the sanity warning, exit 0.

        Contig 6 is ``G`` in the reference; the cache claims ``C``.  Both
        sanity checks (ancestral reconstruction and derived-site grouping)
        compare the cache's base to Prodigal's forward-oriented base and must
        warn.  With 17-vs-3 counts there is no tie, so the wrong tie-breaker
        changes nothing: the NS call at 6 survives.
        """
        done, outdir = _run_dnds(tmp_path, ref_base_override={6: "C"})
        assert done.returncode == 0, done.stderr
        assert "Reference base mismatch" in done.stderr, done.stderr
        # the warning names the planted position and both sides of the disagreement
        mismatch_lines = [
            l for l in done.stderr.splitlines() if "Reference base mismatch" in l
        ]
        assert any(
            f"{CONTIG}:6" in l or " 6 in " in l for l in mismatch_lines
        ), mismatch_lines
        assert all("'G'" in l and "'C'" in l for l in mismatch_lines), mismatch_lines
        # only position 6 was corrupted, so nothing else should be flagged
        assert all("6" in l for l in mismatch_lines), mismatch_lines

        events = _table(outdir, "codon_events_ng86.tsv")
        assert len(events) == 4
        row6 = events[events["contig_position"].astype(int) == 6].iloc[0]
        assert row6["mutation_type"] == "NS" and row6["derived_allele"] == "A"
