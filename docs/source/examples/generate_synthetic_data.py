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
    python generate_synthetic_data.py --output_dir my_test_data
    python generate_synthetic_data.py --num_mags 10 --num_samples 20 --output_dir large_test
"""

import argparse
import gzip
import logging
import random
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

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
    length = (length // 3) * 3   # round down to codon boundary
    if length < 6:
        length = 6

    seq = "ATG"   # mandatory start codon

    stop_codons = {"TAA", "TAG", "TGA"}
    while len(seq) < length - 3:          # -3 leaves room for the stop codon
        codon = "".join(random.choices(BASES, k=3))
        if codon not in stop_codons:
            seq += codon

    seq += random.choice(["TAA", "TAG", "TGA"])   # append stop codon
    return seq


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
    """
    mag = {"mag_id": mag_id, "contigs": []}

    for i in range(num_contigs):
        contig_id = f"{mag_id}.fa_contig{i+1}"
        length = random.randint(*contig_length_range)
        sequence = generate_random_sequence(length)

        genes = []
        used_regions = []   # list of (start, end) for overlap checking

        for j in range(genes_per_contig):
            gene_id = f"{contig_id}_gene{j+1}"
            gene_length = random.randint(150, min(400, length // 2))

            for _ in range(50):   # 50 placement attempts per gene
                start = random.randint(50, length - gene_length - 50)
                end = start + gene_length

                # Reject if this region overlaps any already-placed gene.
                # Condition: two intervals [s1,e1] and [s2,e2] do NOT overlap
                # iff e1 < s2 or s1 > e2 — so they DO overlap if neither holds.
                overlap = any(
                    not (end < us or start > ue)
                    for us, ue in used_regions
                )
                if not overlap:
                    used_regions.append((start, end))
                    genes.append({
                        "gene_id": gene_id,
                        "start": start,
                        "end": end,
                        "strand": random.choice([1, -1]),
                        "sequence": generate_gene_sequence(gene_length),
                    })
                    break

        mag["contigs"].append({
            "contig_id": contig_id,
            "length": length,
            "sequence": sequence,
            "genes": sorted(genes, key=lambda x: x["start"]),
        })

    return mag


def write_reference_fasta(mags: List[Dict], output_path: Path) -> None:
    """Write a multi-FASTA with one entry per contig, lines wrapped at 80 bp."""
    with open(output_path, "w") as f:
        for mag in mags:
            for contig in mag["contigs"]:
                f.write(f">{contig['contig_id']}\n")
                seq = contig["sequence"]
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i + 80] + "\n")
    logger.info(f"Wrote reference FASTA: {output_path}")


def write_prodigal_genes(mags: List[Dict], output_path: Path) -> None:
    """Write gene predictions in Prodigal FASTA format.

    Prodigal headers look like:
      >gene_id # start # end # strand # ID=...;partial=00;start_type=ATG;...

    AlleleFlux's profile_mags.py parses this format to build a position →
    gene_id lookup table. The gc_cont and other metadata fields are mocked
    but structurally correct.
    """
    with open(output_path, "w") as f:
        for mag in mags:
            for contig in mag["contigs"]:
                for gene in contig["genes"]:
                    header = (
                        f">{gene['gene_id']} # {gene['start']} # {gene['end']} "
                        f"# {gene['strand']} # "
                        f"ID={gene['gene_id']};partial=00;start_type=ATG;"
                        f"rbs_motif=GGAG;gc_cont=0.55"
                    )
                    f.write(header + "\n")
                    f.write(gene["sequence"] + "\n")
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
    """Write a mock GTDB-Tk taxonomy table.

    The full GTDB-Tk output has ~20 columns; only 'user_genome' and
    'classification' are used downstream. The rest are filled with
    plausible-looking but arbitrary values. MAGs cycle through 4 different
    phyla so that the taxonomy-level score aggregation has something to group.
    """
    taxonomies = [
        "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;"
        "f__Bacteroidaceae;g__Bacteroides;s__Bacteroides vulgatus",
        "d__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;"
        "f__Lachnospiraceae;g__Lachnoclostridium;s__Lachnoclostridium sp001",
        "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;"
        "f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
        "d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Bifidobacteriales;"
        "f__Bifidobacteriaceae;g__Bifidobacterium;s__Bifidobacterium longum",
    ]
    header = (
        "user_genome\tclassification\tfastani_reference\tfastani_reference_radius\t"
        "fastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\t"
        "closest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\t"
        "closest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\t"
        "other_related_references\tmsa_percent\ttranslation_table\tred_value\twarnings"
    )
    with open(output_path, "w") as f:
        f.write(header + "\n")
        for i, mag in enumerate(mags):
            tax = taxonomies[i % len(taxonomies)]   # cycle through phyla
            f.write(
                f"{mag['mag_id']}\t{tax}\tGCF_000000001.1\t95.0\t{tax}\t98.5\t0.95\t"
                f"GCF_000000001.1\t95.0\t{tax}\t98.5\t0.95\t{tax}\tANI\tN/A\tN/A\t"
                f"85.0\t11\tN/A\tN/A\n"
            )
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
                    samples.append({
                        "sample_id": f"{group}_subj{subject_id}_{timepoint}",
                        "subjectID": f"subj{subject_id}",
                        "group": group,
                        "replicate": replicate,
                        "time": timepoint,
                    })
                subject_id += 1

    else:
        # Single-timepoint design: one row per sample, no pairing needed.
        samples_per_group = num_samples // len(groups)
        subject_id = 1
        replicate_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        for group in groups:
            for rep_idx in range(samples_per_group):
                replicate = replicate_letters[rep_idx % len(replicate_letters)]
                samples.append({
                    "sample_id": f"{group}_sample{subject_id}",
                    "subjectID": f"subj{subject_id}",
                    "group": group,
                    "replicate": replicate,
                    "time": timepoints[0] if timepoints else "t1",
                })
                subject_id += 1

    return samples


def write_metadata(samples: List[Dict], profiles_dir: Path, output_path: Path) -> None:
    """Write the sample metadata TSV that AlleleFlux reads at startup.

    The bam_path column stores the path to the profile *directory* for each
    sample (not a BAM file). When profiles_path is set in config.yml, AlleleFlux
    uses these paths as pointers to pre-generated profiles rather than running
    the profiling step.

    Paths are written relative to output_dir so the metadata file is portable —
    the same file works regardless of where the dataset is placed on disk.
    """
    output_dir = output_path.parent.parent   # metadata/ → output_dir
    with open(output_path, "w") as f:
        f.write("sample_id\tbam_path\tsubjectID\tgroup\treplicate\ttime\n")
        for sample in samples:
            profile_path = profiles_dir / sample["sample_id"]
            rel_path = profile_path.relative_to(output_dir)
            f.write(
                f"{sample['sample_id']}\t{rel_path}\t{sample['subjectID']}\t"
                f"{sample['group']}\t{sample['replicate']}\t{sample['time']}\n"
            )
    logger.info(f"Wrote metadata: {output_path}")


def generate_profile_data(
    sample: Dict,
    mags: List[Dict],
    coverage_range: Tuple[int, int] = (20, 80),
    treatment_effect_rate: float = 0.08,
) -> Dict[str, List[str]]:
    """Simulate the per-position allele count table for one sample.

    This is the core function. It mimics what alleleflux-profile produces from
    a BAM file, but without any reads or alignment. For each genomic position
    it directly constructs the counts (A, C, G, T) that a pileup would return.

    Reference base assignment
    -------------------------
    The reference base at each position cycles through BASES by index:
      pos 0 → A,  pos 1 → T,  pos 2 → G,  pos 3 → C,
      pos 4 → A,  pos 5 → T,  ...

    This is deterministic: the same position always gets the same ref_base
    across all samples, which is what AlleleFlux expects (the reference genome
    doesn't change between samples).

    Allele count model
    ------------------
    For each position, one allele is designated "dominant". It receives
    85–98% of the coverage; the remaining 2–15% reads are distributed randomly
    among the other three bases (capped at 3 reads each to simulate low-level
    sequencing noise rather than real polymorphism).

    Treatment effect — the allele shift
    ------------------------------------
    ONLY in samples where group == "treatment" AND time == "post" does the
    dominant allele have a chance of shifting.

    At each position, a Bernoulli draw with p = treatment_effect_rate (default
    0.08, i.e. 8%) decides whether the shift happens:

      if is_treatment and is_post and random.random() < 0.08:
          dominant_base = BASES[(pos % 4 + 1) % 4]   # rotate one step forward

    The rotation mapping is:
      ref A (pos%4 == 0) → dominant becomes T  (index 1)
      ref T (pos%4 == 1) → dominant becomes G  (index 2)
      ref G (pos%4 == 2) → dominant becomes C  (index 3)
      ref C (pos%4 == 3) → dominant becomes A  (index 0, wraps around)

    Concrete example — position 0, coverage = 50:
      All other samples:         A=44  T=2  G=1  C=0   (A frequency ≈ 0.88)
      treatment_subj3_post (shifted): T=46  A=2  G=1  C=0  (A frequency ≈ 0.04)

    AlleleFlux computes Δ = post_freq − pre_freq per subject, then tests
    whether those deltas differ between groups. At this position, the treatment
    group has Δ_A ≈ 0.04 − 0.88 = −0.84, while the control group has
    Δ_A ≈ 0.88 − 0.88 = 0 → a large, detectable difference.

    The shift is applied independently and randomly at each position, so ~8%
    of all positions in treatment post-timepoint samples carry this signal.

    Returns
    -------
    dict mapping mag_id → list of tab-separated row strings (one per position)
    """
    is_treatment = sample["group"] == "treatment"
    is_post = sample["time"] == "post"

    profiles = {}

    for mag in mags:
        rows = []

        for contig in mag["contigs"]:
            # Pre-build a flat list of (start, end, gene_id) for O(n) gene lookup.
            gene_ranges = [
                (gene["start"], gene["end"], gene["gene_id"])
                for gene in contig["genes"]
            ]

            for pos in range(contig["length"]):
                coverage = random.randint(*coverage_range)

                # Deterministic reference base: cycles A→T→G→C→A→...
                ref_base = BASES[pos % 4]
                dominant_base = ref_base

                # ---- Treatment effect ----------------------------------------
                # Independently for each position, flip a biased coin. If it
                # lands heads (prob = treatment_effect_rate), rotate the dominant
                # allele one step forward in the BASES cycle.
                # This only fires for treatment post-timepoint samples.
                if is_treatment and is_post and random.random() < treatment_effect_rate:
                    dominant_base = BASES[(pos % 4 + 1) % 4]
                # --------------------------------------------------------------

                counts = {"A": 0, "T": 0, "G": 0, "C": 0}

                # The dominant allele captures the vast majority of reads.
                dominant_count = int(coverage * random.uniform(0.85, 0.98))
                counts[dominant_base] = dominant_count
                remaining = coverage - dominant_count

                # Distribute leftover reads among the three non-dominant bases.
                # Cap at 3 each to keep minor allele frequencies low (noise, not signal).
                other_bases = [b for b in BASES if b != dominant_base]
                for b in other_bases:
                    if remaining > 0:
                        c = random.randint(0, min(remaining, 3))
                        counts[b] = c
                        remaining -= c

                total_coverage = counts["A"] + counts["C"] + counts["G"] + counts["T"]
                n_count = 0   # synthetic reads have no ambiguous bases

                # Find which gene (if any) this position falls inside.
                gene_id = ""
                for start, end, gid in gene_ranges:
                    if start <= pos <= end:
                        gene_id = gid
                        break

                rows.append(
                    f"{contig['contig_id']}\t{pos}\t{ref_base}\t{total_coverage}"
                    f"\t{counts['A']}\t{counts['C']}\t{counts['G']}\t{counts['T']}"
                    f"\t{n_count}\t{gene_id}"
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
            f.write("contig\tposition\tref_base\ttotal_coverage\tA\tC\tG\tT\tN\tgene_id\n")
            for row in rows:
                f.write(row + "\n")


def generate_significant_sites(
    mags: List[Dict], output_path: Path, num_sites_per_gene: int = 5
) -> None:
    """Write a mock significant-sites table for visualisation testing.

    These sites are randomly sampled from gene bodies with fabricated p-values.
    They are NOT derived from running statistical tests — they exist solely so
    that visualisation tools (e.g. alleleflux-terminal-nuc-analysis) have
    something to display without requiring a full pipeline run first.
    """
    with open(output_path, "w") as f:
        f.write("mag_id\tcontig\tposition\tgene_id\ttest_type\tmin_p_value\tq_value\n")
        for mag in mags:
            for contig in mag["contigs"]:
                for gene in contig["genes"]:
                    n = min(num_sites_per_gene, gene["end"] - gene["start"])
                    positions = sorted(random.sample(range(gene["start"], gene["end"]), n))
                    for pos in positions:
                        p_val = random.uniform(1e-7, 1e-4)
                        q_val = p_val * random.uniform(5, 20)
                        f.write(
                            f"{mag['mag_id']}\t{contig['contig_id']}\t{pos}\t"
                            f"{gene['gene_id']}\ttwo_sample_paired_tTest\t"
                            f"{p_val:.2e}\t{q_val:.2e}\n"
                        )
    logger.info(f"Wrote significant sites: {output_path}")


def write_config(
    output_dir: Path, data_type: str, groups: List[str], timepoints: List[str]
) -> None:
    """Write a ready-to-run AlleleFlux config pointing at the generated data.

    Note: LMM and CMH are enabled here (unlike the bundled example config which
    disables them for speed). Adjust as needed for your use case.
    """
    config = f"""# AlleleFlux Configuration (Generated)
# =====================================

run_name: "generated_test"

input:
  fasta_path: "{output_dir}/reference/combined_mags.fasta"
  prodigal_path: "{output_dir}/reference/prodigal_genes.fna"
  metadata_path: "{output_dir}/metadata/sample_metadata.tsv"
  gtdb_path: "{output_dir}/reference/gtdbtk_taxonomy.tsv"
  mag_mapping_path: "{output_dir}/reference/mag_mapping.tsv"
  # Use the pre-generated profiles directly — skips the profiling step.
  profiles_path: "{output_dir}/profiles"

output:
  root_dir: "{output_dir}/output"

log_level: "INFO"

analysis:
  data_type: "{data_type}"
  allele_analysis_only: False
  # LMM and CMH require more replicates than the small example provides.
  # Enable them when you generate larger datasets.
  use_lmm: False
  use_significance_tests: True
  use_cmh: False
  use_regional_contrast: True
  use_gene_scores: False
  use_outlier_detection: False

  taxa_score_levels:
    - phylum
    - genus

  timepoints_combinations:
    - timepoint: {timepoints}
      focus: "{timepoints[-1] if len(timepoints) > 1 else timepoints[0]}"

  groups_combinations:
    - {groups}

quality_control:
  min_sample_num: 2
  breadth_threshold: 0.1
  coverage_threshold: 1
  disable_zero_diff_filtering: False

profiling:
  ignore_orphans: True
  min_base_quality: 30
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
  dn_ds_test_type: "two_sample_paired_tTest"

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
    config_path = output_dir / "config_generated.yml"
    with open(config_path, "w") as f:
        f.write(config)
    logger.info(f"Wrote config: {config_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic AlleleFlux test data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--output_dir", "-o", type=str, required=True,
                        help="Output directory for generated data")
    parser.add_argument("--num_mags", "-m", type=int, default=2,
                        help="Number of MAGs to generate")
    parser.add_argument("--num_samples", "-s", type=int, default=8,
                        help="Number of samples (must be divisible by groups × timepoints)")
    parser.add_argument("--num_contigs_per_mag", type=int, default=2,
                        help="Number of contigs per MAG")
    parser.add_argument("--genes_per_contig", type=int, default=3,
                        help="Number of genes to place per contig")
    parser.add_argument("--contig_length_min", type=int, default=800,
                        help="Minimum contig length (bp)")
    parser.add_argument("--contig_length_max", type=int, default=2000,
                        help="Maximum contig length (bp)")
    parser.add_argument("--data_type", type=str, choices=["single", "longitudinal"],
                        default="longitudinal",
                        help="Experimental design: single timepoint or longitudinal")
    parser.add_argument("--groups", type=str, nargs="+",
                        default=["control", "treatment"],
                        help="Experimental group names")
    parser.add_argument("--timepoints", type=str, nargs="+",
                        default=["pre", "post"],
                        help="Timepoint labels (longitudinal only)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed — use the same seed to reproduce identical data")
    parser.add_argument("--coverage_min", type=int, default=20,
                        help="Minimum per-position coverage")
    parser.add_argument("--coverage_max", type=int, default=80,
                        help="Maximum per-position coverage")

    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    output_dir = Path(args.output_dir)
    for subdir in ["reference", "metadata", "profiles", "significant_sites"]:
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)

    logger.info(f"Generating synthetic data in: {output_dir}")
    logger.info(f"{args.num_mags} MAGs, {args.num_samples} samples, {args.data_type} design")

    # --- MAG structures -------------------------------------------------
    mags = []
    for i in range(args.num_mags):
        mag_id = f"MAG_{i+1:03d}"   # zero-padded: MAG_001, MAG_002, ...
        mags.append(create_mag_structure(
            mag_id,
            num_contigs=args.num_contigs_per_mag,
            contig_length_range=(args.contig_length_min, args.contig_length_max),
            genes_per_contig=args.genes_per_contig,
        ))

    # --- Reference files ------------------------------------------------
    write_reference_fasta(mags, output_dir / "reference" / "combined_mags.fasta")
    write_prodigal_genes(mags, output_dir / "reference" / "prodigal_genes.fna")
    write_mag_mapping(mags, output_dir / "reference" / "mag_mapping.tsv")
    write_gtdb_taxonomy(mags, output_dir / "reference" / "gtdbtk_taxonomy.tsv")

    # --- Sample structure -----------------------------------------------
    samples = generate_sample_structure(
        args.num_samples,
        groups=args.groups,
        timepoints=args.timepoints,
        data_type=args.data_type,
    )

    profiles_dir = output_dir / "profiles"
    write_metadata(samples, profiles_dir, output_dir / "metadata" / "sample_metadata.tsv")

    # --- Profiles (the slow step for large datasets) --------------------
    logger.info(f"Generating profiles for {len(samples)} samples...")
    for i, sample in enumerate(samples):
        if (i + 1) % 10 == 0 or i == len(samples) - 1:
            logger.info(f"  {i+1}/{len(samples)}: {sample['sample_id']}")
        profiles = generate_profile_data(
            sample, mags, coverage_range=(args.coverage_min, args.coverage_max)
        )
        write_profiles(sample, profiles, profiles_dir)

    # --- Ancillary files ------------------------------------------------
    generate_significant_sites(
        mags, output_dir / "significant_sites" / "significant_sites.tsv"
    )
    write_config(output_dir, args.data_type, args.groups, args.timepoints)

    # Write a minimal README next to the generated data
    readme = (
        f"# Generated AlleleFlux Test Data\n\n"
        f"Seed: {args.seed} | MAGs: {args.num_mags} | "
        f"Samples: {len(samples)} ({args.data_type})\n\n"
        f"```bash\nalleleflux run --config {output_dir}/config_generated.yml\n```\n"
    )
    (output_dir / "README.md").write_text(readme)

    total_size = sum(f.stat().st_size for f in output_dir.rglob("*") if f.is_file())
    logger.info(f"Done. Total size: {total_size / 1024:.1f} KB")
    logger.info(f"Run with: alleleflux run --config {output_dir}/config_generated.yml")


if __name__ == "__main__":
    main()
