#!/usr/bin/env python3
"""
Generate Synthetic AlleleFlux Test Data

Produces all input files needed to run the AlleleFlux pipeline WITHOUT real
sequencing data. Instead of starting from BAM files and running the profiling
step, this script writes the profile TSVs directly, simulating what the
profiler would produce.

What gets generated
-------------------
  reference/combined_mags.fasta    — reference genome (random sequences)
  reference/prodigal_genes.fna     — gene coordinates in Prodigal header format
  reference/mag_mapping.tsv        — maps every contig to its MAG
  reference/gtdbtk_taxonomy.tsv    — mock GTDB taxonomy rows
  metadata/sample_metadata.tsv     — 8-row sample sheet
  profiles/{sample}/               — one gzipped TSV per sample × MAG
  significant_sites/               — dummy significant-site table for visualisation
  config_generated.yml             — ready-to-run AlleleFlux config

Key design choice: treatment effect is injected directly into read counts.
This is impossible to do at the BAM level with a generic read simulator (wgsim
treats every sample identically). By writing profiles directly we can control
exactly which positions are "significant" in which samples.

Usage:
    # Regenerate the bundled longitudinal dataset (paired pre/post samples)
    python generate_synthetic_data.py \
        --output_dir example_data_longitudinal \
        --data_type longitudinal --seed 42

    # Regenerate the bundled single-timepoint dataset
    python generate_synthetic_data.py \
        --output_dir example_data_single \
        --data_type single --seed 42

    # Ad-hoc larger dataset for one-off experiments
    python generate_synthetic_data.py \
        --output_dir my_test_data \
        --num_mags 10 --num_samples 20
"""

import argparse
import gzip
import json
import logging
import random
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# The four DNA bases in a fixed order.
# The ORDER matters: the treatment effect shifts each position's dominant
# allele one step forward in this cycle (A→T→G→C→A).
BASES = ["A", "T", "G", "C"]

# ---- Planted-truth knobs ----------------------------------------------------
# Every fact the regression test asserts is planted with these numbers and
# written to planted_truth.json, so a test never has to guess them.
DEFAULT_COVERAGE_RANGE = (20, 80)  # per-position depth draw, uniform
N_STANDING = 12  # significant sites whose rising allele already exists at baseline
SHIFT_RISING_FRACTION = 0.48  # rising allele's read share in a carrier sample: large
#   enough for a paired t-test on 4 vs 4 mice, but BELOW one half so the
#   consensus base stays the reference and the shift is NOT a strain change
#   in the ANI step (conANI counts consensus differences only)
STANDING_FRACTION = 0.15  # rising allele's read share in the ONE control mouse that
#   carries each standing allele, at pre AND post alike: present at baseline
#   (far above min_freq 0.05 and the per-depth error bar), and unchanged over
#   time in that mouse, so the control group's pre-to-post change stays ~0 and
#   the paired tests keep their power (a pre-only carrier cost 11 of 12
#   standing sites their significance at 4 mice per group)
STANDING_MIN_COVERAGE = 40  # that mouse's depth floor at standing sites: >= 6 reads
BLOCK_SPACING_BP = 300  # three de novo alleles this far apart: one cluster at
#   the default 1000 bp window, so they must become one block_event
SWAP_FRACTION = 0.02  # share of positions whose CONSENSUS rotates in the one
#   strain-swapped sample: conANI ~0.98 against its own earlier sample, far
#   below the 0.999 same-strain line, while every other mouse stays above it
SPIKE_LENGTH = 40  # positions per depth spike
SPIKE_FACTOR = 10  # spike depth = this x the mean of the coverage range: the
#   site-depth filter's default ceiling is 3 x the sample mean, so every
#   spike position must be removed for that sample
MIN_PLANT_MARGIN = 50  # keep plants this far from contig ends


def generate_random_sequence(length: int) -> str:
    """Return a random DNA string of the given length."""
    return "".join(random.choices(BASES, k=length))


def generate_gene_sequence(length: int) -> str:
    """Return a plausible coding sequence: ATG … (non-stop codons) … stop.

    The sequence is padded/trimmed so its length is divisible by 3, starts
    with ATG (met), and ends with a stop codon. Stop codons are explicitly
    excluded from the internal codons so Prodigal-style parsers don't
    truncate the gene early.
    """
    length = (length // 3) * 3  # round down to codon boundary
    if length < 6:
        length = 6

    seq = "ATG"  # mandatory start codon

    # Standard genetic code stop codons (TAA/TAG/TGA), from Biopython's
    # translation table rather than a hardcoded literal.
    stop_codons = CodonTable.standard_dna_table.stop_codons
    stop_set = set(stop_codons)
    while len(seq) < length - 3:  # -3 leaves room for the stop codon
        codon = "".join(random.choices(BASES, k=3))
        if codon not in stop_set:
            seq += codon

    seq += random.choice(stop_codons)  # append stop codon
    return seq


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of an A/T/G/C string.

    Prodigal stores every gene in its own 5'→3' direction; for a strand=-1 gene
    that direction is the reverse strand, so the forward-strand bases occupying
    those coordinates on the contig are the reverse complement of the stored
    gene sequence.
    """
    return str(Seq(seq).reverse_complement())


def create_mag_structure(
    mag_id: str,
    num_contigs: int = 2,
    contig_length_range: Tuple[int, int] = (800, 2000),
    genes_per_contig: int = 3,
) -> Dict:
    """Build the in-memory representation of one MAG.

    Returns a dict:
      {
        "mag_id": "MAG_001",
        "contigs": [
          {
            "contig_id": "MAG_001.fa_contig1",
            "length": 1500,
            "sequence": "ATGC...",
            "genes": [
              {"gene_id": ..., "start": 100, "end": 280, "strand": 1,
               "sequence": "ATG...TAA"},
              ...
            ]
          },
          ...
        ]
      }

    Gene placement uses a rejection-sampling loop: draw a random start, check
    it doesn't overlap any previously placed gene, retry up to 50 times. Genes
    are placed at least 50 bp from either contig end to leave non-genic flanks.

    Contig/gene sequence consistency
    --------------------------------
    The contig sequence is built in TWO stages:
      1. Start with a random A/T/G/C background of the full contig length.
      2. After gene coordinates are chosen, overwrite contig[start:end] with
         the gene's coding sequence — reverse-complemented when strand == -1
         so the forward-strand bases on the contig really ARE the reverse
         complement of the gene sequence stored in prodigal_genes.fna.

    With this in place, contig.sequence, prodigal_genes.fna, and the profile
    ref_base column (derived from contig.sequence in generate_profile_data)
    all agree at gene coordinates — so dN/dS reconstruction's sanity check
    against the Prodigal reference doesn't fire warnings on every variable site.
    """
    mag = {"mag_id": mag_id, "contigs": []}

    for i in range(num_contigs):
        contig_id = f"{mag_id}.fa_contig{i+1}"
        length = random.randint(*contig_length_range)

        # Stage 1: random A/T/G/C background. Intergenic regions stay random;
        # gene regions get overwritten in stage 2 below.
        contig_chars = list(generate_random_sequence(length))

        genes = []
        used_regions = []  # list of (start, end) for overlap checking

        for j in range(genes_per_contig):
            gene_id = f"{contig_id}_gene{j+1}"
            gene_length = random.randint(150, min(400, length // 2))

            for _ in range(50):  # 50 placement attempts per gene
                start = random.randint(50, length - gene_length - 50)
                end = start + gene_length

                # Reject if this region overlaps any already-placed gene.
                # Condition: two intervals [s1,e1] and [s2,e2] do NOT overlap
                # iff e1 < s2 or s1 > e2 — so they DO overlap if neither holds.
                overlap = any(not (end < us or start > ue) for us, ue in used_regions)
                if not overlap:
                    strand = random.choice([1, -1])
                    gene_seq = generate_gene_sequence(gene_length)
                    # gene_seq may be 1–2 bp shorter than gene_length because
                    # generate_gene_sequence rounds down to a codon boundary.
                    # Trust the actual embedded length so end stays consistent
                    # with what's physically present in contig.sequence.
                    embed_seq = (
                        reverse_complement(gene_seq) if strand == -1 else gene_seq
                    )
                    embed_len = len(embed_seq)
                    end = start + embed_len

                    # Stage 2: weave the gene's bases into the contig at its
                    # coordinates, so contig.sequence[start:end] reads back as
                    # the (forward-strand) gene sequence.
                    contig_chars[start:end] = list(embed_seq)

                    used_regions.append((start, end))
                    genes.append(
                        {
                            "gene_id": gene_id,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "sequence": gene_seq,
                        }
                    )
                    break

        mag["contigs"].append(
            {
                "contig_id": contig_id,
                "length": length,
                "sequence": "".join(contig_chars),
                "genes": sorted(genes, key=lambda x: x["start"]),
            }
        )

    return mag


def write_reference_fasta(mags: List[Dict], output_path: Path) -> None:
    """Write a multi-FASTA with one entry per contig via Biopython's SeqIO.

    Sequences are wrapped at Biopython's default 60 bp/line; wrap width is
    cosmetic — every downstream consumer reads this file with SeqIO.parse.
    """
    records = [
        SeqRecord(Seq(contig["sequence"]), id=contig["contig_id"], description="")
        for mag in mags
        for contig in mag["contigs"]
    ]
    SeqIO.write(records, str(output_path), "fasta")
    logger.info(f"Wrote reference FASTA: {output_path}")


def write_prodigal_genes(mags: List[Dict], output_path: Path) -> None:
    """Write gene predictions in Prodigal FASTA format.

    Prodigal headers look like:
      >gene_id # start # end # strand # ID=...;partial=00;start_type=ATG;...

    AlleleFlux's profile_mags.py parses this format to build a position →
    gene_id lookup table. The gc_cont and other metadata fields are mocked
    but structurally correct.
    """
    records = []
    for mag in mags:
        for contig in mag["contigs"]:
            for gene in contig["genes"]:
                # Prodigal coordinates are 1-based INCLUSIVE; our gene dict
                # stores 0-based half-open [start, end). Convert before writing:
                #   1-based start          = start + 1
                #   1-based inclusive end  = end  (last base is 0-based end-1)
                # Both map_genes (start-1, end-1) and dnds_from_timepoints
                # (forward: pos-(start-1)) decode these as 1-based, so emitting
                # raw 0-based here misaligns forward-strand genes by one base.
                prodigal_start = gene["start"] + 1
                prodigal_end = gene["end"]
                # Biopython writes the header as ">{id} {description}", so the
                # coordinate block goes in description. On read-back,
                # record.description is the full line and splits on '#' into the
                # 5 fields map_genes expects (gene_id, start, end, strand, attrs).
                description = (
                    f"# {prodigal_start} # {prodigal_end} # {gene['strand']} # "
                    f"ID={gene['gene_id']};partial=00;start_type=ATG;"
                    f"rbs_motif=GGAG;gc_cont=0.55"
                )
                records.append(
                    SeqRecord(
                        Seq(gene["sequence"]),
                        id=gene["gene_id"],
                        description=description,
                    )
                )
    SeqIO.write(records, str(output_path), "fasta")
    logger.info(f"Wrote Prodigal genes: {output_path}")


def write_mag_mapping(mags: List[Dict], output_path: Path) -> None:
    """Write the contig → MAG mapping TSV (contig_id, mag_id).

    AlleleFlux uses this to group per-contig pileup data back into per-MAG
    profiles and to route each position to the correct output file.
    """
    with open(output_path, "w") as f:
        f.write("contig_id\tmag_id\n")
        for mag in mags:
            for contig in mag["contigs"]:
                f.write(f"{contig['contig_id']}\t{mag['mag_id']}\n")
    logger.info(f"Wrote MAG mapping: {output_path}")


def write_gtdb_taxonomy(mags: List[Dict], output_path: Path) -> None:
    """Write a GTDB taxonomy table with only the columns AlleleFlux consumes.

    Real GTDB-Tk classify summary.tsv has 19 columns, but the pipeline reads it
    via ``pd.read_csv(..., usecols=["user_genome", "classification"])`` (see
    utilities.read_gtdb) — it selects those two BY NAME and ignores everything
    else. So we emit just those two columns; the other 17 are dead weight here.

    The classification value must be a ';'-separated lineage with the standard
    d__/p__/c__/o__/f__/g__/s__ rank prefixes; parse_classification splits on
    those to build the per-rank columns used by taxon-level score aggregation.

    The lineages are deliberately arranged so they OVERLAP at lower ranks —
    otherwise every taxon-level group would have size 1 and the aggregation step
    (taxa_score_levels: phylum, genus) would have nothing to pool. Cycling order
    (used via taxonomies[i % len]):
      tax0, tax1  → same phylum (Bacteroidota) AND genus (Bacteroides),
                    differing only at species. Default num_mags=2 hits exactly
                    these two, giving an immediate 2-member phylum+genus group.
      tax2, tax3  → same phylum (Bacillota) but different genus, so at
                    num_mags=4 the phylum level has two 2-member groups while
                    genus has one 2-member group plus singletons.
    """
    taxonomies = [
        # --- Bacteroidota / Bacteroides : overlap down to genus --------------
        "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;"
        "f__Bacteroidaceae;g__Bacteroides;s__Bacteroides vulgatus",
        "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;"
        "f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis",
        # --- Bacillota : shared phylum, divergent genus ----------------------
        "d__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;"
        "f__Lachnospiraceae;g__Lachnoclostridium;s__Lachnoclostridium sp001",
        "d__Bacteria;p__Bacillota;c__Clostridia;o__Oscillospirales;"
        "f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium prausnitzii",
    ]
    with open(output_path, "w") as f:
        f.write("user_genome\tclassification\n")
        for i, mag in enumerate(mags):
            tax = taxonomies[i % len(taxonomies)]  # cycle through the 4 lineages
            f.write(f"{mag['mag_id']}\t{tax}\n")
    logger.info(f"Wrote GTDB taxonomy: {output_path}")


def generate_sample_structure(
    num_samples: int,
    groups: List[str] = ["control", "treatment"],
    timepoints: List[str] = ["pre", "post"],
    data_type: str = "longitudinal",
) -> List[Dict]:
    """Build the list of sample metadata dicts.

    For a longitudinal design with num_samples=8, groups=["control","treatment"],
    timepoints=["pre","post"], the function produces exactly 8 rows:

      control_subj1_pre   / subj1 / control / A / pre
      control_subj1_post  / subj1 / control / A / post
      control_subj2_pre   / subj2 / control / B / pre
      control_subj2_post  / subj2 / control / B / post
      treatment_subj3_pre / subj3 / treatment / A / pre
      treatment_subj3_post/ subj3 / treatment / A / post
      treatment_subj4_pre / subj4 / treatment / B / pre
      treatment_subj4_post/ subj4 / treatment / B / post

    subject_id increments globally (not per-group) so each subject has a unique
    ID across groups. The replicate letter is assigned per subject within a
    group (A for the first subject in each group, B for the second, etc.) and
    is used by the CMH and paired tests to match subjects across timepoints.
    """
    samples = []

    if data_type == "longitudinal":
        # How many subjects per group?
        # e.g. 8 samples / (2 groups × 2 timepoints) = 2 subjects/group
        samples_per_group = num_samples // (len(groups) * len(timepoints))
        if samples_per_group < 1:
            samples_per_group = 1

        subject_id = 1
        replicate_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        for group in groups:
            for rep_idx in range(samples_per_group):
                replicate = replicate_letters[rep_idx % len(replicate_letters)]
                # Both timepoints for this subject are emitted together so they
                # share the same subjectID, which is what the paired tests require.
                for timepoint in timepoints:
                    samples.append(
                        {
                            "sample_id": f"{group}_subj{subject_id}_{timepoint}",
                            "subjectID": f"subj{subject_id}",
                            "group": group,
                            "replicate": replicate,
                            "time": timepoint,
                        }
                    )
                subject_id += 1

    else:
        # Single-timepoint design: one row per sample, no pairing needed.
        samples_per_group = num_samples // len(groups)
        subject_id = 1
        replicate_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        for group in groups:
            for rep_idx in range(samples_per_group):
                replicate = replicate_letters[rep_idx % len(replicate_letters)]
                samples.append(
                    {
                        "sample_id": f"{group}_sample{subject_id}",
                        "subjectID": f"subj{subject_id}",
                        "group": group,
                        "replicate": replicate,
                        "time": timepoints[0] if timepoints else "t1",
                    }
                )
                subject_id += 1

    return samples


def collection_dates_for(timepoints: List[str]) -> Dict[str, str]:
    """One ISO calendar date per timepoint label, three months apart.

    The pipeline's cohort baseline scope orders samples by ``collection_date``
    (``YYYY-MM-DD``); without dates it falls back to label order with a
    warning, so the bundled dataset carries real dates.  Example:
    ``["pre", "post"]`` -> ``{"pre": "2024-01-15", "post": "2024-04-15"}``.
    """
    dates = {}
    for i, tp in enumerate(timepoints):
        month_index = 3 * i  # 0, 3, 6, ... months after January 2024
        dates[tp] = f"{2024 + month_index // 12}-{month_index % 12 + 1:02d}-15"
    return dates


def _ordered_unique(values):
    """Unique values in first-seen order (``dict.fromkeys`` idiom)."""
    return list(dict.fromkeys(values))


def _all_positions(mags: List[Dict]):
    """Every (mag_id, contig_id, position) at least MIN_PLANT_MARGIN from a contig end."""
    for mag in mags:
        for contig in mag["contigs"]:
            for pos in range(MIN_PLANT_MARGIN, contig["length"] - MIN_PLANT_MARGIN):
                yield mag["mag_id"], contig["contig_id"], pos


def plan_the_plants(
    mags: List[Dict],
    samples: List[Dict],
    data_type: str,
    *,
    n_standing: int = N_STANDING,
    coverage_range: Tuple[int, int] = DEFAULT_COVERAGE_RANGE,
) -> Dict:
    """Decide every planted fact up front and return the card that records them.

    The card is what ``generate_profile_data`` reads while writing counts and
    what ``generate_dataset`` saves as ``planted_truth.json``; the regression
    test reads the same card, so a fact exists in exactly one place.  Nothing
    in it depends on the profiles: it is pure bookkeeping over the MAG
    structures and the sample sheet, deterministic under the module's seed.

    What gets planted
    -----------------
    * ``significant_sites``: positions where the treatment group's LATER samples
      (the ``carrier_samples``) carry the ``rising_allele`` -- the next base in
      the A->T->G->C cycle after the reference -- at ``shift_rising_fraction``
      of their reads.  Longitudinal sites come in two kinds: ``standing``
      (``n_standing`` of them; the rising allele is also present at baseline,
      as a stable minority in ONE control mouse at every timepoint, see
      ``standing_mouse``) and ``de_novo`` (the
      rising allele has zero reads in every non-carrier sample).  The de novo
      sites are one ``lone`` allele and a ``block`` of three alleles
      ``block_spacing_bp`` apart on one contig; the lone one sits on a
      different contig so no cluster window can reach it.  A single-timepoint
      design has no baseline, so its sites are all kind ``shift``.
    * ``strain_swap``: the last treatment mouse's later sample rotates its
      CONSENSUS at ``swap_fraction`` of each contig's positions (never at a
      planted site).
    * ``depth_spikes``: two control samples each get one ``SPIKE_LENGTH``
      stretch at ``SPIKE_FACTOR`` x the mean depth, placed clear of every
      planted site.

    Example (longitudinal, 16 samples, two MAGs of two contigs): 12 standing
    sites, 1 lone + 3 block de novo sites, carrier_samples = the four
    ``treatment_subj*_post`` ids, standing_mouse = ``subj1`` (its pre and post
    samples, 15 % of reads), strain_swap on ``treatment_subj8_post``, spikes in
    ``control_subj1_pre`` and ``control_subj2_post``.
    """
    if data_type not in ("single", "longitudinal"):
        raise ValueError(f"data_type must be single or longitudinal, got {data_type!r}")
    groups = _ordered_unique(s["group"] for s in samples)
    if len(groups) != 2:
        raise ValueError(f"the plan needs exactly two groups, got {groups}")
    control, treatment = groups
    timepoints = _ordered_unique(s["time"] for s in samples)
    later = timepoints[-1]

    # Which samples carry the rising allele: the treatment group at the later
    # timepoint (longitudinal) or simply the treatment group (single).
    carriers = [s for s in samples if s["group"] == treatment and s["time"] == later]
    controls = [s for s in samples if s["group"] == control]
    if not carriers or len(controls) < 2:
        raise ValueError(
            "the plan needs treatment samples at the later timepoint and two control samples"
        )

    seq_of = {
        (mag["mag_id"], contig["contig_id"]): contig["sequence"]
        for mag in mags
        for contig in mag["contigs"]
    }
    len_of = {key: len(seq) for key, seq in seq_of.items()}

    def site(mag_id, contig_id, pos, kind, label, **extra):
        ref = seq_of[(mag_id, contig_id)][pos]
        return {
            "mag": mag_id,
            "contig": contig_id,
            "position": pos,
            "ref_base": ref,
            "rising_allele": BASES[(BASES.index(ref) + 1) % 4],
            "kind": kind,
            "label": label,
            **extra,
        }

    taken = set()  # (mag, contig, pos) already used by a plant
    sites = []

    if data_type == "longitudinal":
        # The block goes on the longest contig (it needs 2 x spacing of room);
        # the lone allele on the longest OTHER contig, so they can never share
        # a cluster key (MAG x contig).
        by_length = sorted(len_of, key=lambda k: len_of[k], reverse=True)
        if len(by_length) < 2:
            raise ValueError("the longitudinal plan needs at least two contigs")
        (block_mag, block_contig), (lone_mag, lone_contig) = by_length[0], by_length[1]
        block_len = len_of[(block_mag, block_contig)]
        room = block_len - 2 * MIN_PLANT_MARGIN - 2 * BLOCK_SPACING_BP
        if room < 1:
            raise ValueError(
                f"contig {block_contig} ({block_len} bp) is too short for a "
                f"block of three sites {BLOCK_SPACING_BP} bp apart"
            )
        start = random.randint(MIN_PLANT_MARGIN, MIN_PLANT_MARGIN + room - 1)
        for k in range(3):
            pos = start + k * BLOCK_SPACING_BP
            sites.append(
                site(block_mag, block_contig, pos, "de_novo", "block", block_index=k)
            )
            taken.add((block_mag, block_contig, pos))
        lone_len = len_of[(lone_mag, lone_contig)]
        lone_pos = random.randint(MIN_PLANT_MARGIN, lone_len - MIN_PLANT_MARGIN - 1)
        sites.append(site(lone_mag, lone_contig, lone_pos, "de_novo", "lone"))
        taken.add((lone_mag, lone_contig, lone_pos))
        standing_kind = "standing"
    else:
        standing_kind = "shift"

    # The remaining listed sites, spread over every contig.
    free = [p for p in _all_positions(mags) if p not in taken]
    for key in sorted(random.sample(free, n_standing)):
        sites.append(site(*key, standing_kind, standing_kind))
        taken.add(key)
    # Stable order for the card: by MAG, contig, position.
    sites.sort(key=lambda s: (s["mag"], s["contig"], s["position"]))

    # Depth spikes: two control samples, each one stretch clear of every plant.
    spike_samples = [controls[0]["sample_id"], controls[-1]["sample_id"]]
    spike_depth = SPIKE_FACTOR * (coverage_range[0] + coverage_range[1]) // 2
    spikes = []
    for sample_id, (mag_id, contig_id) in zip(spike_samples, sorted(len_of)):
        length = len_of[(mag_id, contig_id)]
        for _ in range(1000):  # placement attempts
            start = random.randint(
                MIN_PLANT_MARGIN, length - MIN_PLANT_MARGIN - SPIKE_LENGTH
            )
            stretch = range(start, start + SPIKE_LENGTH)
            if not any((mag_id, contig_id, p) in taken for p in stretch):
                break
        else:
            raise ValueError(
                f"could not place a depth spike on {contig_id} clear of the plants"
            )
        spikes.append(
            {
                "sample_id": sample_id,
                "mag": mag_id,
                "contig": contig_id,
                "start": start,
                "end": start + SPIKE_LENGTH,  # half-open, like gene coordinates
                "depth": spike_depth,
            }
        )

    # Strain swap: the last carrier's consensus rotates at SWAP_FRACTION of each
    # contig's positions, never at a planted site.
    swap_sample = carriers[-1]
    swap_positions: Dict[str, Dict[str, List[int]]] = {}
    n_swap = 0
    for (mag_id, contig_id), length in sorted(len_of.items()):
        candidates = [p for p in range(length) if (mag_id, contig_id, p) not in taken]
        n = round(SWAP_FRACTION * length)
        chosen = sorted(random.sample(candidates, n))
        swap_positions.setdefault(mag_id, {})[contig_id] = chosen
        n_swap += n

    truth = {
        "data_type": data_type,
        "groups": {"control": control, "treatment": treatment},
        "timepoints": timepoints,
        "collection_dates": collection_dates_for(timepoints),
        "rising_allele_rule": "the next base after ref_base in the A,T,G,C cycle",
        "shift_rising_fraction": SHIFT_RISING_FRACTION,
        "n_standing": n_standing,
        "block_spacing_bp": BLOCK_SPACING_BP,
        "carrier_samples": [s["sample_id"] for s in carriers],
        "significant_sites": sites,
        "strain_swap": {
            "sample_id": swap_sample["sample_id"],
            "subjectID": swap_sample["subjectID"],
            "timepoint": swap_sample["time"],
            "fraction": SWAP_FRACTION,
            "n_positions": n_swap,
            "positions": swap_positions,
        },
        "depth_spikes": spikes,
    }
    if data_type == "longitudinal":
        # The first control mouse, every one of its samples: the standing
        # allele is a stable minority there, not a pre-only visitor.
        subject = controls[0]["subjectID"]
        truth["standing_mouse"] = {
            "subjectID": subject,
            "sample_ids": [
                s["sample_id"] for s in controls if s["subjectID"] == subject
            ],
            "fraction": STANDING_FRACTION,
        }
    return truth


def write_metadata(
    samples: List[Dict],
    output_path: Path,
    *,
    path_resolver: Callable[[Dict], Path],
    collection_dates: Dict[str, str],
) -> None:
    """Write the sample metadata TSV that AlleleFlux reads at startup.

    The bam_path column holds the path AlleleFlux uses to locate each sample's
    data. Its meaning depends on the run mode:
      - "profile" mode: bam_path points at the sample's profile *directory*.
        When profiles_path is set in config.yml, AlleleFlux uses these as
        pointers to pre-generated profiles and skips profiling.
      - "BAM" mode: bam_path points at a sample's .bam *file* and AlleleFlux
        runs the profiling step.

    Rather than baking either convention into this function, the caller passes
    a `path_resolver(sample) -> Path` callable that yields the right path. Two
    typical call sites — both passing **relative** paths anchored at the
    dataset's example_dir (so they resolve against the pipeline's cwd, which
    by convention is set to example_dir):

        write_metadata(..., path_resolver=lambda s: Path("profiles") / s["sample_id"])
        write_metadata(..., path_resolver=lambda s: Path("bams") / f"{s['sample_id']}.bam")

    Paths are written AS-GIVEN (no ``.resolve()``).  The README + the
    regression test both ``cd`` into example_dir before running the
    pipeline, so relative-to-example_dir paths resolve correctly.  This makes
    the bundled dataset portable across machines (CI runners, contributor
    laptops, HPC clusters) without re-running the generator — which the
    earlier absolute-path design forbade.
    """
    # collection_date is one of the pipeline's optional metadata columns
    # (ISO YYYY-MM-DD, validated at alleleflux-metadata); every sample of a
    # timepoint shares that timepoint's date, see collection_dates_for.
    with open(output_path, "w") as f:
        f.write(
            "sample_id\tbam_path\tsubjectID\tgroup\treplicate\ttime\tcollection_date\n"
        )
        for sample in samples:
            bam_path = path_resolver(sample)
            f.write(
                f"{sample['sample_id']}\t{bam_path}\t{sample['subjectID']}\t"
                f"{sample['group']}\t{sample['replicate']}\t{sample['time']}\t"
                f"{collection_dates[sample['time']]}\n"
            )
    logger.info(f"Wrote metadata: {output_path}")


def _background_counts(
    dominant_base: str, coverage: int, noise_bases: List[str]
) -> Dict[str, int]:
    """Monomorphic-with-noise counts: the dominant base takes 85-98 % of the
    reads, the leftover is scattered over ``noise_bases`` at <= 3 reads each.

    random.uniform(0.85, 0.98) keeps every position below a variant caller's
    minor-allele threshold but leaves a sliver of background reads.  The
    3-read cap keeps each minor count below the usual polymorphism threshold
    so it reads as sequencer noise.  The realised depth can fall a few reads
    short of ``coverage`` when the minor rolls land on 0; total_coverage is
    recomputed from the counts downstream, so the row stays consistent.

    Example: dominant A, coverage 50, noise over T/G/C -> A=46 T=2 G=1 C=0.
    """
    counts = {"A": 0, "T": 0, "G": 0, "C": 0}
    dominant_count = int(coverage * random.uniform(0.85, 0.98))
    counts[dominant_base] = dominant_count
    remaining = coverage - dominant_count
    for b in noise_bases:
        if remaining > 0:
            c = random.randint(0, min(remaining, 3))
            counts[b] = c
            remaining -= c
    return counts


def generate_profile_data(
    sample: Dict,
    mags: List[Dict],
    coverage_range: Tuple[int, int] = DEFAULT_COVERAGE_RANGE,
    *,
    truth: Dict,
) -> Dict[str, List[str]]:
    """Simulate the per-position allele count table for one sample.

    This is the core function. It mimics what alleleflux-profile produces from
    a BAM file, but without any reads or alignment: for each genomic position
    it constructs the counts (A, C, G, T) a pileup would return.  What happens
    at a position is decided by the planted-truth card from ``plan_the_plants``
    -- nothing significant is left to a coin flip, so every fact the card
    records is really in the profiles and nothing else is.

    Reference base: read from the contig's forward-strand sequence built once
    in create_mag_structure, so every sample shares the same ref_base and the
    profiles agree with prodigal_genes.fna (no dN/dS mismatch warnings).

    Per position, in this order of precedence:

    1. **Planted significant site** (``truth["significant_sites"]``):
       - in a *carrier sample* (treatment at the later timepoint) the
         ``rising_allele`` takes ``shift_rising_fraction`` (0.48) of the reads
         and the reference the rest, no noise: ref A, coverage 50 ->
         A=26 T=24.  The consensus stays A, so the ANI step sees no strain
         change; the frequency change (about -0.4 against the ~0.9 baseline)
         is what the paired tests detect;
       - in the ``standing_mouse`` samples (one control mouse, pre and post),
         at a *standing* site, the rising allele gets ``fraction`` (0.15) of at
         least 40 reads and the reference the rest: A=34 T=6 -> present at
         baseline unambiguously, and unchanged in that mouse over time;
       - in every other sample the reference is dominant with the usual noise
         on the two bases that are NEITHER ref NOR rising, and the rising
         allele has exactly 0 reads: a de novo allele is absent everywhere it
         should be, and a standing one is present in exactly one baseline mouse.
    2. **Strain swap** (``truth["strain_swap"]``, that sample only): at the
       listed positions the dominant base is rotated one step in the BASES
       cycle (ref A -> dominant T) with the usual noise -- a consensus change,
       which is what conANI counts.
    3. **Depth spike** (``truth["depth_spikes"]``, that sample only): coverage
       is the planted ``depth`` instead of a draw from ``coverage_range``;
       counts are otherwise background.
    4. **Background** everywhere else: reference dominant, noise on the other
       three bases.

    Returns dict mapping mag_id -> list of tab-separated row strings.
    """
    sample_id = sample["sample_id"]
    is_carrier = sample_id in truth["carrier_samples"]
    standing = truth.get("standing_mouse") or {}
    is_standing_mouse = sample_id in standing.get("sample_ids", ())
    standing_fraction = standing.get("fraction", 0.0)

    # Lookups for this sample, built once per call (a few hundred entries).
    site_at = {(s["contig"], s["position"]): s for s in truth["significant_sites"]}
    swap = truth["strain_swap"]
    swap_at = set()
    if swap["sample_id"] == sample_id:
        swap_at = {
            (contig, pos)
            for contigs in swap["positions"].values()
            for contig, positions in contigs.items()
            for pos in positions
        }
    spike_depth_at = {}
    for spike in truth["depth_spikes"]:
        if spike["sample_id"] == sample_id:
            for pos in range(spike["start"], spike["end"]):
                spike_depth_at[(spike["contig"], pos)] = spike["depth"]

    profiles = {}

    for mag in mags:
        rows = []

        for contig in mag["contigs"]:
            contig_seq = contig["sequence"]
            contig_id = contig["contig_id"]
            # Pre-build a flat list of (start, end, gene_id) for O(n) gene lookup.
            # end is the half-open upper bound: gene occupies [start, end).
            gene_ranges = [
                (gene["start"], gene["end"], gene["gene_id"])
                for gene in contig["genes"]
            ]

            for pos in range(contig["length"]):
                key = (contig_id, pos)
                coverage = spike_depth_at.get(key) or random.randint(*coverage_range)
                ref_base = contig_seq[pos]
                site = site_at.get(key)

                if site is not None:
                    rising = site["rising_allele"]
                    if is_carrier:
                        # Rising allele just under half, reference the rest, no
                        # noise.  The cap keeps rising STRICTLY below ref: at
                        # even depths <= 24, round(0.48 x depth) equals half the
                        # depth and the consensus would be a tie.
                        rising_n = min(
                            round(coverage * truth["shift_rising_fraction"]),
                            (coverage - 1) // 2,
                        )
                        counts = {"A": 0, "T": 0, "G": 0, "C": 0}
                        counts[rising] = rising_n
                        counts[ref_base] = coverage - rising_n
                    elif is_standing_mouse and site["kind"] == "standing":
                        coverage = max(coverage, STANDING_MIN_COVERAGE)
                        standing_n = round(coverage * standing_fraction)
                        counts = {"A": 0, "T": 0, "G": 0, "C": 0}
                        counts[rising] = standing_n
                        counts[ref_base] = coverage - standing_n
                    else:
                        # Reference dominant; noise only on the bases that are
                        # neither ref nor rising, so the rising allele is 0.
                        noise_bases = [b for b in BASES if b not in (ref_base, rising)]
                        counts = _background_counts(ref_base, coverage, noise_bases)
                elif key in swap_at:
                    dominant_base = BASES[(BASES.index(ref_base) + 1) % 4]
                    counts = _background_counts(
                        dominant_base,
                        coverage,
                        [b for b in BASES if b != dominant_base],
                    )
                else:
                    counts = _background_counts(
                        ref_base, coverage, [b for b in BASES if b != ref_base]
                    )
                    if key in spike_depth_at:
                        # A spike must realise EXACTLY the planted depth (the
                        # test and the site filter both read it), so the reads
                        # the noise rounding dropped go back to the reference.
                        counts[ref_base] += coverage - sum(counts.values())

                total_coverage = counts["A"] + counts["C"] + counts["G"] + counts["T"]

                # Find which gene (if any) this position falls inside.
                # Half-open: a position equal to end is just past the gene's
                # last base, which lives at end - 1.
                gene_id = ""
                for start, end, gid in gene_ranges:
                    if start <= pos < end:
                        gene_id = gid
                        break

                rows.append(
                    f"{contig['contig_id']}\t{pos}\t{ref_base}\t{total_coverage}"
                    f"\t{counts['A']}\t{counts['C']}\t{counts['G']}\t{counts['T']}"
                    f"\t{gene_id}"
                )

        profiles[mag["mag_id"]] = rows

    return profiles


def write_profiles(
    sample: Dict, profiles: Dict[str, List[str]], output_dir: Path
) -> None:
    """Write one gzipped TSV per MAG for this sample.

    Output path: {output_dir}/{sample_id}/{sample_id}_{mag_id}_profiled.tsv.gz

    The column order matches what alleleflux-profile produces from a real BAM
    so downstream rules can read either source without modification.
    """
    sample_dir = output_dir / sample["sample_id"]
    sample_dir.mkdir(parents=True, exist_ok=True)

    for mag_id, rows in profiles.items():
        filepath = sample_dir / f"{sample['sample_id']}_{mag_id}_profiled.tsv.gz"
        with gzip.open(filepath, "wt") as f:
            f.write("contig\tposition\tref_base\ttotal_coverage\tA\tC\tG\tT\tgene_id\n")
            for row in rows:
                f.write(row + "\n")


def generate_significant_sites(
    truth: Dict, mags: List[Dict], output_path: Path
) -> None:
    """Write a mock significant-sites table for visualisation testing.

    One row per planted significant site from the card, with a fixed
    placeholder p-value.  The rows are NOT derived from running statistical
    tests -- they exist so that visualisation tools (e.g.
    alleleflux-terminal-nuc-analysis) have something to display without a
    full pipeline run first -- but they now name the positions that really
    ARE significant in the generated data, instead of random gene positions.
    """
    gene_at = {}
    for mag in mags:
        for contig in mag["contigs"]:
            for gene in contig["genes"]:
                for pos in range(gene["start"], gene["end"]):
                    gene_at[(contig["contig_id"], pos)] = gene["gene_id"]
    test_type = (
        "two_sample_paired_tTest"
        if truth["data_type"] == "longitudinal"
        else "two_sample_unpaired_tTest"
    )
    with open(output_path, "w") as f:
        f.write("mag_id\tcontig\tposition\tgene_id\ttest_type\tmin_p_value\tq_value\n")
        for site in truth["significant_sites"]:
            gene_id = gene_at.get((site["contig"], site["position"]), "")
            f.write(
                f"{site['mag']}\t{site['contig']}\t{site['position']}\t"
                f"{gene_id}\t{test_type}\t1.00e-06\t1.00e-05\n"
            )
    logger.info(f"Wrote significant sites: {output_path}")


def write_config(
    output_dir: Path,
    data_type: str,
    groups: List[str],
    timepoints: List[str],
    *,
    mode: str = "profiles",
) -> Path:
    """Write a ready-to-run AlleleFlux config pointing at the generated data.

    The two run modes differ in just five values; the rest of the YAML is
    identical, so they share one template here instead of duplicating the
    full ~70 lines twice.

    Parameters
    ----------
    mode : "profiles" or "bams"
        - "profiles" (default): config references the pre-generated profiles in
          <output_dir>/profiles/ via profiles_path. Profiling is skipped.
          Filename: config_generated_{data_type}.yml.
        - "bams": config omits profiles_path so the pipeline runs profiling
          from BAMs at <output_dir>/bams/<sample>.bam, with min_base_quality
          lowered to 0 (wgsim writes Phred ~23). metadata_path points at
          sample_metadata_bams.tsv. Output lands in output_bams/.
          Filename: config_with_bams_{data_type}.yml.

    dn_ds_test_type is chosen to match the experimental design:
      - longitudinal → "two_sample_paired_tTest" (pre/post pairing on the same
        subject is the natural unit of comparison).
      - single       → "two_sample_unpaired_tTest" (no pre/post pairing exists
        at a single timepoint). The pipeline silently skips dN/dS in single
        mode regardless (see dynamic_targets.smk), so this is informational,
        but the wrong value would confuse a reader inspecting the config.

    Filename reflects the data_type so single-mode and longitudinal-mode runs
    in adjacent dirs don't get confused, and a user inspecting a folder can
    tell at a glance which design the config targets.
    """
    if mode not in {"profiles", "bams"}:
        raise ValueError(f"mode must be 'profiles' or 'bams', got {mode!r}")

    if data_type == "single":
        dnds_test_type = "two_sample_unpaired_tTest"
    else:
        dnds_test_type = "two_sample_paired_tTest"

    # ---- data-type-specific feature block ------------------------------------
    # Every feature the data type supports is switched ON, so the regression
    # run exercises the whole surface and the planted facts have somewhere to
    # surface.  The four longitudinal-only features (dN/dS, regional contrast,
    # baseline presence, strain turnover) are written as explicit False on
    # single data: the workflow refuses an explicit True there at parse time,
    # and an explicit False reads clearer than a missing key.
    if data_type == "longitudinal":
        longitudinal_block = """  use_dnds: True
  # Regional contrast is experimental (not in the publication): stays off.
  use_regional_contrast: False
  # Pairwise conANI/popANI between each mouse's pre and post sample; the
  # strain-swapped mouse (planted_truth.json: strain_swap) is the one pair
  # whose conANI falls below con_threshold.
  use_pairwise_ani: True
  pairwise_ani:
    min_cov: 5
    min_freq: 0.05
    fdr: 1.0e-6
    pairs: "within_subject"
    store_snp_locations: "within_subject"
  use_strain_turnover: True
  strain_turnover:
    min_compared: 0.1
    pop_threshold: 0.99999
    con_threshold: 0.999
    # 4 mice per group: the real-cohort floor of 8 voters would leave every
    # verdict blank (too_few_voters), so the test data votes with 2.
    min_voters: 2
    vote_rule: "majority"
    tie: "not_replaced"
    replicate_rule: "average"
  # Standing vs de novo for every significant site of the paired t-test, judged
  # against every sample dated before the first carrier (cohort scope; the
  # metadata carries collection_date), with the carrier-level strain gate.
  use_baseline_presence: True
  baseline_presence:
    summary: "two_sample_paired"
    test_type: "two_sample_paired_tTest"
    threshold_column: "q_value"
    threshold: 0.05
    baseline_scope: "cohort"
    cluster_window_bp: 1000
    max_de_novo_neighbours: 0
    strain_gate: "carrier"
    # The planted de novo alleles have 4 carriers, but a lone allele carried by
    # fewer than 3 determined mice would be strain_undetermined at the house
    # default; 1 keeps every carried allele judgeable on this small design.
    min_determined_carriers: 1
  # Score heatmaps, boxplots and the HTML run report (images; excluded from
  # the regression comparison, but the rules run).
  figures:
    enabled: True"""
    else:
        longitudinal_block = """  # Longitudinal-only features: single data has no earlier timepoint.
  use_dnds: False
  use_regional_contrast: False
  use_baseline_presence: False
  use_strain_turnover: False
  # Pairwise conANI/popANI between EVERY pair of samples (the between-mouse
  # strain-sharing question; the only mode with a meaning on one timepoint).
  # The planted divergent sample (planted_truth.json: strain_swap) is the one
  # whose pairs all fall below con_threshold.
  use_pairwise_ani: True
  pairwise_ani:
    min_cov: 5
    min_freq: 0.05
    fdr: 1.0e-6
    pairs: "all"
    store_snp_locations: "none"
  # Off until the score-figure loader handles a single-timepoint CMH table
  # (it reports a score/summary MAG mismatch for cmh at t1).
  figures:
    enabled: False"""

    # ---- mode-specific overrides --------------------------------------------
    # Everything else in the YAML is shared. Tweak just these five values plus
    # the filename and header comment between the two modes.
    if mode == "profiles":
        header = (
            "# AlleleFlux Configuration (Generated)\n"
            "# =====================================\n"
            "# Paths are RELATIVE to this config's directory (example_data_<dt>/).\n"
            "# Run the pipeline from inside example_data_<dt>/ so they resolve correctly:\n"
            "#   cd docs/source/examples/example_data_" + data_type + "\n"
            "#   alleleflux run --config " + f"config_generated_{data_type}.yml\n"
        )
        run_name = "generated_test"
        metadata_filename = "sample_metadata.tsv"
        # profiles_path line is present in profile mode.
        profiles_line = (
            f"  # Use the pre-generated profiles directly — skips the profiling step.\n"
            f'  profiles_path: "profiles"'
        )
        output_subdir = "output"
        # 30 is the default upstream threshold; works for real sequencing data.
        min_base_quality = 30
        bq_comment = ""
        config_filename = f"config_generated_{data_type}.yml"
    else:  # mode == "bams"
        header = (
            "# AlleleFlux Configuration (Generated, BAM mode)\n"
            "# =================================================\n"
            "# Paths are RELATIVE to this config's directory (example_data_<dt>/).\n"
            "# This config runs the FULL pipeline starting from BAM files — the\n"
            "# profiling step is NOT skipped. Generate the BAMs first with:\n"
            "#\n"
            f"#   bash generate_bams.sh --data-dir example_data_{data_type}\n"
            "#\n"
            "# Then run the pipeline from inside example_data_<dt>/:\n"
            f"#   cd docs/source/examples/example_data_{data_type}\n"
            f"#   alleleflux run --config config_with_bams_{data_type}.yml\n"
        )
        run_name = "generated_test_bams"
        metadata_filename = "sample_metadata_bams.tsv"
        # profiles_path is intentionally absent — pipeline profiles from BAMs.
        profiles_line = "  # profiles_path intentionally omitted — pipeline profiles directly from BAMs."
        output_subdir = "output_bams"
        # wgsim writes Phred ~23 regardless of error rate; raise to ≥20 for real data.
        min_base_quality = 0
        bq_comment = (
            "  # wgsim writes Phred ~23 to every base; raise to ≥20 for real data.\n"
        )
        config_filename = f"config_with_bams_{data_type}.yml"

    config = f"""{header}
run_name: "{run_name}"

input:
  fasta_path: "reference/combined_mags.fasta"
  prodigal_path: "reference/prodigal_genes.fna"
  metadata_path: "metadata/{metadata_filename}"
  gtdb_path: "reference/gtdbtk_taxonomy.tsv"
  mag_mapping_path: "reference/mag_mapping.tsv"
{profiles_line}

output:
  root_dir: "./{output_subdir}"

log_level: "INFO"

analysis:
  data_type: "{data_type}"
  allele_analysis_only: False
  # Every feature the data type supports is ON (the bundled dataset is the
  # regression test's input; planted_truth.json says what each should find),
  # with one exception: LMM.  On this 4-mice-per-group design it took ~58
  # CPU-minutes for two 5 kb genomes (every other step together: ~10) and
  # called a third of all positions significant, so it is off here; its own
  # unit tests cover it.  The within-group across-time tests run
  # (run_within_group_tests default).
  use_lmm: False
  use_significance_tests: True
  use_cmh: True
  use_gene_scores: True
  use_outlier_detection: True
{longitudinal_block}

  taxa_score_levels:
    - phylum
    - genus

  timepoints_combinations:
    - timepoint: {timepoints}
      focus: "{timepoints[-1] if len(timepoints) > 1 else timepoints[0]}"

  groups_combinations:
    - treatment: "{groups[1]}"
      control: "{groups[0]}"

quality_control:
  min_sample_num: 2
  breadth_threshold: 0.1
  coverage_threshold: 1
  disable_zero_diff_filtering: False
  # MIDAS2's per-site depth rule.  The two planted depth spikes
  # (planted_truth.json: depth_spikes) sit at 10x the sample mean, so the
  # 3.0x ceiling must remove exactly those positions for those samples.
  site_depth_filter:
    enabled: True
    site_depth: 0
    site_ratio: 3.0
    min_allele_reads: 0
    keep_original_counts: True

profiling:
  ignore_orphans: True
{bq_comment}  min_base_quality: {min_base_quality}
  min_mapping_quality: 2
  ignore_overlaps: True

statistics:
  filter_type: "t-test"
  preprocess_between_groups: True
  preprocess_within_groups: True
  max_zero_count: 2
  p_value_threshold: 0.05
  fdr_group_by_mag_id: False
  min_positions_after_preprocess: 1

dnds:
  p_value_column: "q_value"
  # dN/dS is only run when data_type == "longitudinal"; ignored for single mode.
  dn_ds_test_type: "{dnds_test_type}"

regional_contrast:
  mode: "both"
  window_size: 500
  agg_method: "median"
  min_informative_sites: 2
  min_informative_fraction: 0.0
  use_fisher: True

resources:
  threads_per_job: 4
  mem_per_job: "4G"
  time: "01:00:00"
  retries: 1
  mem_step: "2G"
  time_step: "0:30:00"
"""
    config_path = output_dir / config_filename
    with open(config_path, "w") as f:
        f.write(config)
    logger.info(f"Wrote config ({mode} mode): {config_path}")
    return config_path


def generate_dataset(
    output_dir,
    *,
    data_type: str = "longitudinal",
    num_mags: int = 2,
    num_samples: int = 8,
    num_contigs: int = 2,
    genes_per_contig: int = 3,
    contig_length_range: Tuple[int, int] = (800, 2000),
    groups: Optional[List[str]] = None,
    timepoints: Optional[List[str]] = None,
    seed: int = 42,
    coverage_range: Tuple[int, int] = DEFAULT_COVERAGE_RANGE,
) -> Dict:
    """Generate one complete dataset into ``output_dir`` and return its card.

    This is ``main`` without the argument parsing, so tests and other scripts
    can build a dataset in a temp folder.  Steps, in order: seed both RNGs
    (every plant and every noise draw is reproducible from ``seed`` alone);
    build the MAG structures and reference files; build the sample sheet; plan
    the plants (the card); write both metadata files (profile mode and BAM
    mode, each with ``collection_date``); write every sample's profiles from
    the card; write the mock significant-sites table, the two configs, the
    README and ``planted_truth.json``.

    Example: ``generate_dataset(Path("ds"), data_type="longitudinal",
    num_samples=16, seed=11)`` writes 16 sample folders under ``ds/profiles/``
    and returns the same dict that ``ds/planted_truth.json`` holds.
    """
    if groups is None:
        groups = ["control", "treatment"]
    if timepoints is None:
        timepoints = ["t1"] if data_type == "single" else ["pre", "post"]

    random.seed(seed)
    np.random.seed(seed)

    # Resolve to absolute so the README's "cd" line points at a real location
    # whatever the caller's cwd was; configs and metadata stay RELATIVE.
    output_dir = Path(output_dir).resolve()
    for subdir in ["reference", "metadata", "profiles", "significant_sites"]:
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)

    logger.info(f"Generating synthetic data in: {output_dir}")
    logger.info(f"{num_mags} MAGs, {num_samples} samples, {data_type} design")

    # --- MAG structures -------------------------------------------------
    mags = []
    for i in range(num_mags):
        mag_id = f"MAG_{i+1:03d}"  # zero-padded: MAG_001, MAG_002, ...
        mags.append(
            create_mag_structure(
                mag_id,
                num_contigs=num_contigs,
                contig_length_range=contig_length_range,
                genes_per_contig=genes_per_contig,
            )
        )

    # --- Reference files ------------------------------------------------
    write_reference_fasta(mags, output_dir / "reference" / "combined_mags.fasta")
    write_prodigal_genes(mags, output_dir / "reference" / "prodigal_genes.fna")
    write_mag_mapping(mags, output_dir / "reference" / "mag_mapping.tsv")
    write_gtdb_taxonomy(mags, output_dir / "reference" / "gtdbtk_taxonomy.tsv")

    # --- Sample structure and the planted-truth card --------------------
    samples = generate_sample_structure(
        num_samples, groups=groups, timepoints=timepoints, data_type=data_type
    )
    truth = plan_the_plants(mags, samples, data_type, coverage_range=coverage_range)
    truth["seed"] = seed

    profiles_dir = output_dir / "profiles"
    # Profile-mode metadata: bam_path points at the per-sample profile directory,
    # written RELATIVE to example_dir so the dataset is portable.  Pipeline
    # convention is to cd into example_dir before invoking, so "profiles/<id>"
    # resolves correctly against cwd.
    write_metadata(
        samples,
        output_dir / "metadata" / "sample_metadata.tsv",
        path_resolver=lambda s: Path("profiles") / s["sample_id"],
        collection_dates=truth["collection_dates"],
    )
    # BAM-mode metadata: bam_path points at bams/<sample>.bam, relative to
    # example_dir (same convention as the profile-mode metadata above).  BAMs
    # don't exist yet (user runs generate_bams.sh after this); the file is
    # written eagerly so the BAM-mode config can reference it without a second
    # invocation.  The pipeline validates BAM existence at runtime.  NOTE: the
    # BAMs are simulated from the unmodified reference, so none of the planted
    # facts reach them -- BAM mode exercises profiling, not the plants.
    write_metadata(
        samples,
        output_dir / "metadata" / "sample_metadata_bams.tsv",
        path_resolver=lambda s: Path("bams") / f"{s['sample_id']}.bam",
        collection_dates=truth["collection_dates"],
    )

    # --- Profiles (the slow step for large datasets) --------------------
    logger.info(f"Generating profiles for {len(samples)} samples...")
    for i, sample in enumerate(samples):
        # Throttle progress logs: every 10th sample plus the final one, so a
        # 17-sample run logs "10/17" and then "17/17" instead of going silent.
        if (i + 1) % 10 == 0 or i == len(samples) - 1:
            logger.info(f"  {i+1}/{len(samples)}: {sample['sample_id']}")
        profiles = generate_profile_data(
            sample, mags, coverage_range=coverage_range, truth=truth
        )
        write_profiles(sample, profiles, profiles_dir)

    # --- Ancillary files ------------------------------------------------
    generate_significant_sites(
        truth, mags, output_dir / "significant_sites" / "significant_sites.tsv"
    )
    config_path = write_config(
        output_dir, data_type, groups, timepoints, mode="profiles"
    )
    config_bams_path = write_config(
        output_dir, data_type, groups, timepoints, mode="bams"
    )

    # The card: the only writer is this function, the readers are the tests.
    card_path = output_dir / "planted_truth.json"
    card_path.write_text(json.dumps(truth, indent=2, sort_keys=True) + "\n")
    logger.info(f"Wrote planted-truth card: {card_path}")

    # Write a minimal README next to the generated data showing BOTH run paths.
    # Configs use relative paths anchored at example_dir, so the pipeline must
    # be invoked from inside this directory (or with cwd set to it).
    config_basename = config_path.name
    config_bams_basename = config_bams_path.name
    readme = (
        f"# Generated AlleleFlux Test Data\n\n"
        f"Seed: {seed} | MAGs: {num_mags} | "
        f"Samples: {len(samples)} ({data_type})\n\n"
        f"All paths in the configs and metadata are RELATIVE to this directory.\n"
        f"`cd` into this directory before running the pipeline.\n\n"
        f"`planted_truth.json` lists every fact planted in the profiles\n"
        f"(significant sites, standing / de novo alleles, the strain-swapped\n"
        f"sample, the depth spikes); the regression test asserts them.\n\n"
        f"## Option A — pre-generated profiles (fastest; carries the planted facts)\n\n"
        f"```bash\n"
        f"cd {output_dir}\n"
        f"alleleflux run --config {config_basename}\n"
        f"```\n\n"
        f"## Option B — full pipeline from BAMs (requires wgsim + bwa)\n\n"
        f"The BAMs are simulated from the unmodified reference, so this path\n"
        f"exercises profiling but carries NONE of the planted facts.\n\n"
        f"```bash\n"
        f"# 1. Generate BAMs from the reference (one-time)\n"
        f"bash generate_bams.sh --data-dir {output_dir}\n\n"
        f"# 2. Run the pipeline (includes profiling step)\n"
        f"cd {output_dir}\n"
        f"alleleflux run --config {config_bams_basename}\n"
        f"```\n"
    )
    (output_dir / "README.md").write_text(readme)

    total_size = sum(f.stat().st_size for f in output_dir.rglob("*") if f.is_file())
    logger.info(f"Done. Total size: {total_size / 1024:.1f} KB")
    logger.info(
        f"Run with pre-generated profiles: "
        f"cd {output_dir} && alleleflux run --config {config_basename}"
    )
    return truth


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic AlleleFlux test data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        required=True,
        help="Output directory for generated data",
    )
    parser.add_argument(
        "--num_mags", "-m", type=int, default=2, help="Number of MAGs to generate"
    )
    parser.add_argument(
        "--num_samples",
        "-s",
        type=int,
        default=8,
        help="Number of samples (must be divisible by groups × timepoints)",
    )
    parser.add_argument(
        "--num_contigs_per_mag", type=int, default=2, help="Number of contigs per MAG"
    )
    parser.add_argument(
        "--genes_per_contig",
        type=int,
        default=3,
        help="Number of genes to place per contig",
    )
    parser.add_argument(
        "--contig_length_min", type=int, default=800, help="Minimum contig length (bp)"
    )
    parser.add_argument(
        "--contig_length_max", type=int, default=2000, help="Maximum contig length (bp)"
    )
    parser.add_argument(
        "--data_type",
        type=str,
        choices=["single", "longitudinal"],
        default="longitudinal",
        help="Experimental design: single timepoint or longitudinal",
    )
    parser.add_argument(
        "--groups",
        type=str,
        nargs="+",
        default=["control", "treatment"],
        help="Experimental group names",
    )
    parser.add_argument(
        "--timepoints",
        type=str,
        nargs="+",
        default=None,
        help="Timepoint labels. Auto-defaults based on --data_type if omitted: "
        "'pre post' for longitudinal, 't1' for single.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed — use the same seed to reproduce identical data",
    )
    parser.add_argument(
        "--coverage_min", type=int, default=20, help="Minimum per-position coverage"
    )
    parser.add_argument(
        "--coverage_max", type=int, default=80, help="Maximum per-position coverage"
    )

    args = parser.parse_args()

    # ---- Resolve & validate timepoints against data_type --------------------
    # Single-timepoint designs have no pre/post pairing, so they need exactly
    # one timepoint label. Longitudinal designs need ≥2. If --timepoints was
    # omitted we pick a sensible default based on --data_type; if it was given
    # explicitly we still validate the count matches the design.
    if args.timepoints is None:
        args.timepoints = ["t1"] if args.data_type == "single" else ["pre", "post"]
    if args.data_type == "single" and len(args.timepoints) != 1:
        parser.error(
            f"--data_type single requires exactly 1 timepoint, got "
            f"{len(args.timepoints)}: {args.timepoints}"
        )
    if args.data_type == "longitudinal" and len(args.timepoints) < 2:
        parser.error(
            f"--data_type longitudinal requires ≥2 timepoints, got "
            f"{len(args.timepoints)}: {args.timepoints}"
        )

    generate_dataset(
        args.output_dir,
        data_type=args.data_type,
        num_mags=args.num_mags,
        num_samples=args.num_samples,
        num_contigs=args.num_contigs_per_mag,
        genes_per_contig=args.genes_per_contig,
        contig_length_range=(args.contig_length_min, args.contig_length_max),
        groups=args.groups,
        timepoints=args.timepoints,
        seed=args.seed,
        coverage_range=(args.coverage_min, args.coverage_max),
    )


if __name__ == "__main__":
    main()
