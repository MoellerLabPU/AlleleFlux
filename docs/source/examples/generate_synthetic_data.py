#!/usr/bin/env python3
"""
Generate Synthetic AlleleFlux Test Data

This script generates synthetic metagenomic data for testing AlleleFlux.
It creates reference FASTAs, gene predictions, metadata, and profile files
that can be used to run the complete AlleleFlux pipeline.

Usage:
    python generate_synthetic_data.py --output_dir my_test_data
    python generate_synthetic_data.py --num_mags 10 --num_samples 20 --output_dir large_test

For visualization testing, use the pre-generated example data in:
    docs/source/examples/example_data/
"""

import argparse
import gzip
import logging
import os
import random
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# Constants
BASES = ["A", "T", "G", "C"]
CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def generate_random_sequence(length: int) -> str:
    """Generate a random DNA sequence."""
    return "".join(random.choices(BASES, k=length))


def generate_gene_sequence(length: int) -> str:
    """Generate a gene-like sequence (divisible by 3, starts with ATG)."""
    # Ensure length is divisible by 3
    length = (length // 3) * 3
    if length < 6:
        length = 6

    # Start with ATG (start codon)
    seq = "ATG"

    # Generate middle codons (avoiding stop codons)
    stop_codons = {"TAA", "TAG", "TGA"}
    while len(seq) < length - 3:
        codon = "".join(random.choices(BASES, k=3))
        if codon not in stop_codons:
            seq += codon

    # End with a stop codon
    seq += random.choice(["TAA", "TAG", "TGA"])

    return seq


def create_mag_structure(
    mag_id: str,
    num_contigs: int = 2,
    contig_length_range: Tuple[int, int] = (800, 2000),
    genes_per_contig: int = 3,
) -> Dict:
    """Create structure for a single MAG with contigs and genes."""
    mag = {"mag_id": mag_id, "contigs": []}

    for i in range(num_contigs):
        contig_id = f"{mag_id}.fa_contig{i+1}"
        length = random.randint(*contig_length_range)

        # Generate sequence
        sequence = generate_random_sequence(length)

        # Generate genes
        genes = []
        used_regions = []

        for j in range(genes_per_contig):
            gene_id = f"{contig_id}_gene{j+1}"
            gene_length = random.randint(150, min(400, length // 2))

            # Find non-overlapping position
            max_attempts = 50
            for _ in range(max_attempts):
                start = random.randint(50, length - gene_length - 50)
                end = start + gene_length

                # Check for overlap
                overlap = False
                for used_start, used_end in used_regions:
                    if not (end < used_start or start > used_end):
                        overlap = True
                        break

                if not overlap:
                    used_regions.append((start, end))
                    strand = random.choice([1, -1])
                    gene_seq = generate_gene_sequence(gene_length)

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
                "sequence": sequence,
                "genes": sorted(genes, key=lambda x: x["start"]),
            }
        )

    return mag


def write_reference_fasta(mags: List[Dict], output_path: Path) -> None:
    """Write combined reference FASTA file."""
    with open(output_path, "w") as f:
        for mag in mags:
            for contig in mag["contigs"]:
                f.write(f">{contig['contig_id']}\n")
                # Write sequence in 80-character lines
                seq = contig["sequence"]
                for i in range(0, len(seq), 80):
                    f.write(seq[i : i + 80] + "\n")
    logger.info(f"Wrote reference FASTA: {output_path}")


def write_prodigal_genes(mags: List[Dict], output_path: Path) -> None:
    """Write Prodigal-format gene predictions."""
    with open(output_path, "w") as f:
        for mag in mags:
            for contig in mag["contigs"]:
                for gene in contig["genes"]:
                    # Prodigal header format: >gene_id # start # end # strand # metadata
                    header = (
                        f">{gene['gene_id']} # {gene['start']} # {gene['end']} # {gene['strand']} # "
                        f"ID={gene['gene_id']};partial=00;start_type=ATG;rbs_motif=GGAG;gc_cont=0.55"
                    )
                    f.write(header + "\n")
                    f.write(gene["sequence"] + "\n")
    logger.info(f"Wrote Prodigal genes: {output_path}")


def write_mag_mapping(mags: List[Dict], output_path: Path) -> None:
    """Write contig-to-MAG mapping file."""
    with open(output_path, "w") as f:
        f.write("contig_name\tmag_id\n")
        for mag in mags:
            for contig in mag["contigs"]:
                f.write(f"{contig['contig_id']}\t{mag['mag_id']}\n")
    logger.info(f"Wrote MAG mapping: {output_path}")


def write_gtdb_taxonomy(mags: List[Dict], output_path: Path) -> None:
    """Write mock GTDB taxonomy file."""
    taxonomies = [
        "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides vulgatus",
        "d__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Lachnoclostridium;s__Lachnoclostridium sp001",
        "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
        "d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium;s__Bifidobacterium longum",
    ]

    header = "user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references\tmsa_percent\ttranslation_table\tred_value\twarnings"

    with open(output_path, "w") as f:
        f.write(header + "\n")
        for i, mag in enumerate(mags):
            tax = taxonomies[i % len(taxonomies)]
            f.write(
                f"{mag['mag_id']}\t{tax}\tGCF_000000001.1\t95.0\t{tax}\t98.5\t0.95\tGCF_000000001.1\t95.0\t{tax}\t98.5\t0.95\t{tax}\tANI\tN/A\tN/A\t85.0\t11\tN/A\tN/A\n"
            )
    logger.info(f"Wrote GTDB taxonomy: {output_path}")


def generate_sample_structure(
    num_samples: int,
    groups: List[str] = ["control", "treatment"],
    timepoints: List[str] = ["pre", "post"],
    data_type: str = "longitudinal",
) -> List[Dict]:
    """Generate sample metadata structure."""
    samples = []

    if data_type == "longitudinal":
        # Longitudinal: each subject has multiple timepoints
        samples_per_group = num_samples // (len(groups) * len(timepoints))
        if samples_per_group < 1:
            samples_per_group = 1

        subject_id = 1
        replicate_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        for group in groups:
            for rep_idx in range(samples_per_group):
                replicate = replicate_letters[rep_idx % len(replicate_letters)]
                for timepoint in timepoints:
                    sample_id = f"{group}_subj{subject_id}_{timepoint}"
                    samples.append(
                        {
                            "sample_id": sample_id,
                            "subjectID": f"subj{subject_id}",
                            "group": group,
                            "replicate": replicate,
                            "time": timepoint,
                        }
                    )
                subject_id += 1
    else:
        # Single timepoint
        samples_per_group = num_samples // len(groups)
        subject_id = 1
        replicate_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        for group in groups:
            for rep_idx in range(samples_per_group):
                replicate = replicate_letters[rep_idx % len(replicate_letters)]
                sample_id = f"{group}_sample{subject_id}"
                samples.append(
                    {
                        "sample_id": sample_id,
                        "subjectID": f"subj{subject_id}",
                        "group": group,
                        "replicate": replicate,
                        "time": timepoints[0] if timepoints else "t1",
                    }
                )
                subject_id += 1

    return samples


def write_metadata(samples: List[Dict], profiles_dir: Path, output_path: Path) -> None:
    """Write sample metadata file."""
    with open(output_path, "w") as f:
        f.write("sample_id\tbam_path\tsubjectID\tgroup\treplicate\ttime\n")
        for sample in samples:
            profile_path = profiles_dir / sample["sample_id"]
            f.write(
                f"{sample['sample_id']}\t{profile_path}\t{sample['subjectID']}\t{sample['group']}\t{sample['replicate']}\t{sample['time']}\n"
            )
    logger.info(f"Wrote metadata: {output_path}")


def generate_profile_data(
    sample: Dict,
    mags: List[Dict],
    coverage_range: Tuple[int, int] = (20, 80),
    treatment_effect_rate: float = 0.08,
) -> Dict[str, List[str]]:
    """Generate profile data for a sample."""
    is_treatment = sample["group"] == "treatment"
    is_post = sample["time"] == "post"

    profiles = {}

    for mag in mags:
        rows = []

        for contig in mag["contigs"]:
            # Build gene position lookup
            gene_ranges = []
            for gene in contig["genes"]:
                gene_ranges.append((gene["start"], gene["end"], gene["gene_id"]))

            for pos in range(contig["length"]):
                # Determine coverage
                coverage = random.randint(*coverage_range)

                # Reference base (use position modulo for consistency)
                ref_idx = pos % 4
                ref_base = BASES[ref_idx]

                # Allele counts
                counts = [0, 0, 0, 0]  # A, T, G, C
                dominant_idx = ref_idx

                # Treatment effect: some positions shift allele in treatment+post
                if is_treatment and is_post and random.random() < treatment_effect_rate:
                    dominant_idx = (ref_idx + 1) % 4

                # Dominant allele gets most reads
                counts[dominant_idx] = int(coverage * random.uniform(0.85, 0.98))
                remaining = coverage - counts[dominant_idx]

                # Distribute remaining among other alleles
                for i in range(4):
                    if i != dominant_idx and remaining > 0:
                        c = random.randint(0, min(remaining, 3))
                        counts[i] = c
                        remaining -= c

                # Find overlapping gene
                gene_id = ""
                for start, end, gid in gene_ranges:
                    if start <= pos <= end:
                        gene_id = gid
                        break

                row = f"{contig['contig_id']}\t{pos}\t{gene_id}\t{ref_base}\t{counts[0]}\t{counts[1]}\t{counts[2]}\t{counts[3]}"
                rows.append(row)

        profiles[mag["mag_id"]] = rows

    return profiles


def write_profiles(
    sample: Dict, profiles: Dict[str, List[str]], output_dir: Path
) -> None:
    """Write profile files for a sample."""
    sample_dir = output_dir / sample["sample_id"]
    sample_dir.mkdir(parents=True, exist_ok=True)

    for mag_id, rows in profiles.items():
        filename = f"{sample['sample_id']}_{mag_id}_profiled.tsv.gz"
        filepath = sample_dir / filename

        with gzip.open(filepath, "wt") as f:
            f.write("contig\tposition\tgene_id\tref_base\tA\tT\tG\tC\n")
            for row in rows:
                f.write(row + "\n")


def generate_significant_sites(
    mags: List[Dict], output_path: Path, num_sites_per_gene: int = 5
) -> None:
    """Generate example significant sites file for visualization testing."""
    with open(output_path, "w") as f:
        f.write("mag_id\tcontig\tposition\tgene_id\ttest_type\tmin_p_value\tq_value\n")

        for mag in mags:
            for contig in mag["contigs"]:
                for gene in contig["genes"]:
                    # Generate a few significant positions per gene
                    positions = sorted(
                        random.sample(
                            range(gene["start"], gene["end"]),
                            min(num_sites_per_gene, gene["end"] - gene["start"]),
                        )
                    )

                    for pos in positions:
                        p_val = random.uniform(1e-7, 1e-4)
                        q_val = p_val * random.uniform(5, 20)
                        f.write(
                            f"{mag['mag_id']}\t{contig['contig_id']}\t{pos}\t{gene['gene_id']}\ttwo_sample_paired_tTest\t{p_val:.2e}\t{q_val:.2e}\n"
                        )

    logger.info(f"Wrote significant sites: {output_path}")


def write_config(
    output_dir: Path, data_type: str, groups: List[str], timepoints: List[str]
) -> None:
    """Write example configuration file."""
    config = f"""# AlleleFlux Configuration (Generated)
# =====================================

run_name: "generated_test"

input:
  fasta_path: "{output_dir}/reference/combined_mags.fasta"
  prodigal_path: "{output_dir}/reference/prodigal_genes.fna"
  metadata_path: "{output_dir}/metadata/sample_metadata.tsv"
  gtdb_path: "{output_dir}/reference/gtdbtk_taxonomy.tsv"
  mag_mapping_path: "{output_dir}/reference/mag_mapping.tsv"

output:
  root_dir: "{output_dir}/output"

log_level: "INFO"

analysis:
  data_type: "{data_type}"
  allele_analysis_only: False
  use_lmm: True
  use_significance_tests: True
  use_cmh: True

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
  ignore_orphans: False
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

resources:
  threads_per_job: 4
  mem_per_job: "4G"
  time: "01:00:00"
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
        "--num_samples", "-s", type=int, default=8, help="Number of samples to generate"
    )
    parser.add_argument(
        "--num_contigs_per_mag", type=int, default=2, help="Number of contigs per MAG"
    )
    parser.add_argument(
        "--genes_per_contig", type=int, default=3, help="Number of genes per contig"
    )
    parser.add_argument(
        "--contig_length_min", type=int, default=800, help="Minimum contig length"
    )
    parser.add_argument(
        "--contig_length_max", type=int, default=2000, help="Maximum contig length"
    )
    parser.add_argument(
        "--data_type",
        type=str,
        choices=["single", "longitudinal"],
        default="longitudinal",
        help="Data type: single or longitudinal",
    )
    parser.add_argument(
        "--groups",
        type=str,
        nargs="+",
        default=["control", "treatment"],
        help="Experimental groups",
    )
    parser.add_argument(
        "--timepoints",
        type=str,
        nargs="+",
        default=["pre", "post"],
        help="Timepoints (for longitudinal data)",
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--coverage_min", type=int, default=20, help="Minimum coverage per position"
    )
    parser.add_argument(
        "--coverage_max", type=int, default=80, help="Maximum coverage per position"
    )

    args = parser.parse_args()

    # Set random seeds
    random.seed(args.seed)
    np.random.seed(args.seed)

    # Create output directory structure
    output_dir = Path(args.output_dir)
    (output_dir / "reference").mkdir(parents=True, exist_ok=True)
    (output_dir / "metadata").mkdir(parents=True, exist_ok=True)
    (output_dir / "profiles").mkdir(parents=True, exist_ok=True)
    (output_dir / "significant_sites").mkdir(parents=True, exist_ok=True)

    logger.info(f"Generating synthetic data in: {output_dir}")
    logger.info(
        f"Configuration: {args.num_mags} MAGs, {args.num_samples} samples, {args.data_type} data"
    )

    # Generate MAGs
    logger.info("Generating MAG structures...")
    mags = []
    for i in range(args.num_mags):
        mag_id = f"MAG_{i+1:03d}"
        mag = create_mag_structure(
            mag_id,
            num_contigs=args.num_contigs_per_mag,
            contig_length_range=(args.contig_length_min, args.contig_length_max),
            genes_per_contig=args.genes_per_contig,
        )
        mags.append(mag)

    # Write reference files
    write_reference_fasta(mags, output_dir / "reference" / "combined_mags.fasta")
    write_prodigal_genes(mags, output_dir / "reference" / "prodigal_genes.fna")
    write_mag_mapping(mags, output_dir / "reference" / "mag_mapping.tsv")
    write_gtdb_taxonomy(mags, output_dir / "reference" / "gtdbtk_taxonomy.tsv")

    # Generate samples
    logger.info("Generating sample structure...")
    samples = generate_sample_structure(
        args.num_samples,
        groups=args.groups,
        timepoints=args.timepoints,
        data_type=args.data_type,
    )

    # Write metadata
    profiles_dir = output_dir / "profiles"
    write_metadata(
        samples, profiles_dir, output_dir / "metadata" / "sample_metadata.tsv"
    )

    # Generate and write profiles
    logger.info(f"Generating profiles for {len(samples)} samples...")
    for i, sample in enumerate(samples):
        if (i + 1) % 10 == 0 or i == len(samples) - 1:
            logger.info(
                f"  Processing sample {i+1}/{len(samples)}: {sample['sample_id']}"
            )

        profiles = generate_profile_data(
            sample, mags, coverage_range=(args.coverage_min, args.coverage_max)
        )
        write_profiles(sample, profiles, profiles_dir)

    # Generate significant sites for visualization testing
    generate_significant_sites(
        mags, output_dir / "significant_sites" / "significant_sites.tsv"
    )

    # Write config file
    write_config(output_dir, args.data_type, args.groups, args.timepoints)

    # Write README
    readme_content = f"""# Generated AlleleFlux Test Data

Generated: {args.seed} (seed)

## Contents

- **{args.num_mags} MAGs** with {args.num_contigs_per_mag} contigs each
- **{len(samples)} samples** ({args.data_type} design)
- Groups: {', '.join(args.groups)}
- Timepoints: {', '.join(args.timepoints)}

## Usage

```bash
alleleflux run --config {output_dir}/config_generated.yml
```

## Files

- `reference/combined_mags.fasta` - Combined MAG reference
- `reference/prodigal_genes.fna` - Gene predictions
- `reference/mag_mapping.tsv` - Contig-to-MAG mapping
- `reference/gtdbtk_taxonomy.tsv` - Mock taxonomy
- `metadata/sample_metadata.tsv` - Sample metadata
- `profiles/` - Pre-generated profile files
- `significant_sites/significant_sites.tsv` - Example significant sites
- `config_generated.yml` - Configuration file
"""

    with open(output_dir / "README.md", "w") as f:
        f.write(readme_content)

    # Calculate and report size
    total_size = sum(f.stat().st_size for f in output_dir.rglob("*") if f.is_file())
    logger.info(f"Complete! Total size: {total_size / 1024:.1f} KB")
    logger.info(f"Run with: alleleflux run --config {output_dir}/config_generated.yml")


if __name__ == "__main__":
    main()
