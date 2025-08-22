#!/usr/bin/env python3
"""
dN/dS Analysis from Ancestral to Derived States
===============================================

This script performs dN/dS (non-synonymous to synonymous substitution ratio) analysis
by comparing ancestral and derived states at significant sites identified by AlleleFlux.

The script accepts the output from p_value_summary.py as input, allowing users to:
- Choose between raw p-values (min_p_value) or FDR-corrected values (q_value)
- Apply significance thresholds for site filtering
- Filter by specific test types if desired

Key Features:
- Reconstruction of ancestral sequences from profile data
- Calculation of potential synonymous and non-synonymous sites using Nei-Gojobori method
- Analysis of actual substitutions between ancestral and derived states
- Global dN/dS calculation avoiding double-counting of codons
- Gene-level and MAG-level summary statistics

Expected Input Format (from p_value_summary.py):
- TSV file with columns: mag_id, contig, position, gene_id, test_type, min_p_value, q_value
- Profile data from AlleleFlux profile_mags.py
- Prodigal-predicted gene sequences and coordinates

Usage:
    python dnds_from_timepoints.py --significant_sites p_value_summary_results.tsv \
                                  --mag_id MAG_001 \
                                  --p_value_column q_value \
                                  --p_value_threshold 0.05 \
                                  --ancestral_sample_id pre_sample \
                                  --derived_sample_id post_sample \
                                  --profile_dir /path/to/profiles \
                                  --prodigal_fasta genes.fasta \
                                  --outdir results/
"""

import argparse
import logging
import multiprocessing as mp
import os
import random
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from alleleflux.scripts.utilities.logging_config import setup_logging

# Set up logger for this script
logger = logging.getLogger(__name__)

# --- Global Constants ---
# Defines the complement for each DNA base. Used for reverse strand genes.
COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# --- dN/dS Pre-computation Functions (Nei-Gojobori Method) ---


def _calculate_codon_sites(
    codon: str, table: CodonTable.CodonTable
) -> tuple[float, float]:
    """
    Helper function to calculate S and N sites for a single codon.

    For a given codon, this function iterates through all 9 possible single-base
    changes and classifies each as synonymous or non-synonymous. It returns
    the fractional count of potential S and N sites for that codon.

    The Nei-Gojobori method calculates potential sites by considering all possible
    single nucleotide changes at each position and determining what fraction would
    be synonymous vs non-synonymous. Since there are 3 possible changes at each
    of the 3 positions (9 total), each position contributes 1/3 to the final count.

    Args:
        codon: A 3-base DNA string (e.g., "ATG").
        table: A Biopython CodonTable object for translation.

    Returns:
        A tuple containing the fractional number of synonymous sites and
        non-synonymous sites for the codon.
    """
    s_codon, n_codon = 0.0, 0.0

    # Explicitly handle stop codons and invalid codons for clarity.
    if codon in table.stop_codons:
        aa_orig = "*"  # Use '*' to represent a stop codon.
    else:
        # Use .get() to handle codons that might not be in the forward table.
        aa_orig = table.forward_table.get(codon)

    # If the codon is invalid (e.g., contains 'N' or '-'), we can't analyze it.
    if aa_orig is None:
        raise ValueError(f"Invalid codon: {codon}. Cannot analyze non-standard codons.")

    # Iterate through each of the 3 positions in the codon (0, 1, 2).
    for i in range(3):
        original_base = codon[i]
        # Iterate through the 3 possible alternative bases at this position.
        for new_base in "ATGC":
            if new_base == original_base:
                continue  # Skip the original base, as it's not a change.

            # Create the new codon by substituting the base.
            new_codon_list = list(codon)
            new_codon_list[i] = new_base
            new_codon_str = "".join(new_codon_list)

            # Check if the new codon is a stop codon.
            if new_codon_str in table.stop_codons:
                aa_new = "*"
            else:
                aa_new = table.forward_table.get(new_codon_str)

            # If the new codon is invalid, count it as non-synonymous.
            if aa_new is None:
                raise ValueError(
                    f"New codon is invalid: {new_codon_str}. Cannot analyze non-standard codons."
                )

            # Compare amino acids to classify the hypothetical mutation.
            if aa_new == aa_orig:
                s_codon += 1  # The change was synonymous.
            else:
                n_codon += 1  # The change was non-synonymous.

    # Each site contributes a fraction (s/3 or n/3) to the total.
    # This averages the potential outcomes over the 3 positions.
    # For example, if position 1 had 2 synonymous and 1 non-synonymous change,
    # position 2 had 1 synonymous and 2 non-synonymous, etc., we average them.
    return s_codon / 3.0, n_codon / 3.0


def _precompute_codon_sites_cache(table_id: int = 11) -> dict:
    """
    Pre-computes and caches the S and N site counts for all 64 possible codons.
    This avoids redundant calculations and significantly improves performance.

    Since there are only 64 possible codons (4^3), we can calculate the S/N sites
    for each one upfront and store them in a dictionary. This cache is used by
    calculate_potential_sites() to avoid recalculating the same codon multiple times.

    Args:
        table_id: The NCBI genetic code table ID (11 = bacterial/archaeal).

    Returns:
        A dictionary mapping each codon string to its (potential_S_sites, potential_N_sites) tuple.
        Example: {"ATG": (0.33, 2.67), "TTT": (0.67, 2.33), ...}
    """
    logger.info(
        f"Pre-computing potential S and N sites for all 64 codons using table {table_id}."
    )
    table = CodonTable.unambiguous_dna_by_id[table_id]
    bases = "ATGC"
    cache = {}
    # Generate all 64 possible codons programmatically.
    for a in bases:
        for b in bases:
            for c in bases:
                codon = a + b + c
                # Calculate S/N sites for each codon and store in the cache.
                # The key is the codon string, the value is the tuple of (S, N) sites.
                cache[codon] = _calculate_codon_sites(codon, table)
    return cache


# --- Pre-computed cache for all 64 codons, created once at script startup ---
_CODON_SITE_CACHE = _precompute_codon_sites_cache(table_id=11)


def calculate_potential_sites_for_gene(gene_seq_str: str) -> dict:
    """
    Calculates potential synonymous (S) and non-synonymous (N) sites for a gene
    using a pre-computed cache for maximum efficiency.

    Args:
        gene_seq_str: The string of the protein-coding DNA sequence.

    Returns:
        A dictionary containing the total counts for 'N' and 'S' sites.
        Example: {"S": 45.67, "N": 123.33}
    """
    potential_S_sites, potential_N_sites = 0.0, 0.0
    # Iterate over the sequence in 3-base steps (codons).
    for i in range(0, len(gene_seq_str), 3):
        # Slice out the current 3-base codon.
        codon = gene_seq_str[i : i + 3].upper()
        # Skip partial codons at the end of a sequence.
        if len(codon) != 3:
            logger.warning(
                f"Gene sequence '{gene_seq_str}' of length {len(gene_seq_str):,} is not a multiple of 3. "
                f"Skipping partial codon at position {i}."
            )
            continue

        # Look up the pre-calculated S and N values from the cache. This is much
        # faster than re-calculating for every codon in every gene.
        # 1. First, get the result from the cache into a single variable.
        potential_codon_sites = _CODON_SITE_CACHE.get(codon)

        # 2. Now, check if the result is None (meaning the codon was not found).
        if potential_codon_sites is None:
            # If the codon is invalid, log a warning and skip it.
            # This handles non-standard bases (like 'Y', 'N') gracefully.
            logger.warning(
                f"Invalid codon '{codon}' found in gene '{gene_seq_str}' at position {i}. Skipping."
            )
            continue

        # 3. If the result is valid, unpack it into the two variables.
        s_codon, n_codon = potential_codon_sites

        # Add the fractional counts for this codon to the gene's total.
        potential_S_sites += s_codon
        potential_N_sites += n_codon

    # Return the final summed totals for the gene.
    return {"S": potential_S_sites, "N": potential_N_sites}


# --- Core Data Loading and Sequence Reconstruction ---


def setup_and_load_data(
    significant_sites_path,
    mag_id,
    p_value_column,
    p_value_threshold,
    test_type,
    group_analyzed,
):
    """
    Loads and filters the initial significant sites data.
    This function reads the significant sites file (output from p_value_summary.py)
    and applies filtering to ensure we only analyze sites that meet our criteria.
    Args:
        significant_sites_path: Path to the significant sites file from p_value_summary.py
        mag_id: The MAG ID to filter for
        p_value_column: Either 'min_p_value' or 'q_value' for significance filtering
        p_value_threshold: Threshold value for significance filtering
        test_type: Optional test type to filter for.
        group_analyzed: Optional group name to filter for.
    Returns:
        A pandas DataFrame of sites to be processed, or None if no sites are found.
    """
    logger.info(f"Loading significant sites from {significant_sites_path}")
    sig_sites_df = pd.read_csv(significant_sites_path, sep="\t")

    # Validate required columns
    required_columns = ["mag_id", "contig", "position", "gene_id", p_value_column]
    missing_columns = [
        col for col in required_columns if col not in sig_sites_df.columns
    ]
    if missing_columns:
        logger.error(
            f"Required columns {missing_columns} not found in significant_sites file. Cannot continue."
        )
        return None

    logger.info(f"Filtering significant sites for mag_id '{mag_id}'")
    sig_sites_df = sig_sites_df[sig_sites_df["mag_id"] == mag_id].copy()
    logger.info(f"Found {len(sig_sites_df):,} sites for mag_id '{mag_id}'")

    if sig_sites_df.empty:
        logger.warning(f"No sites found for mag_id '{mag_id}'.")
        return None

    # Filter by test type if specified
    if test_type:
        if "test_type" in sig_sites_df.columns:
            logger.info(f"Filtering for test type: {test_type} ")
            sig_sites_df = sig_sites_df[sig_sites_df["test_type"] == test_type].copy()
            logger.info(
                f"After test type filtering, {len(sig_sites_df):,} sites remain for mag_id '{mag_id}'."
            )
        else:
            logger.warning(
                "`--test-type` provided but 'test_type' column not found. Skipping test type filtering."
            )

    # Filter by group_analyzed if specified
    if group_analyzed:
        if "group_analyzed" in sig_sites_df.columns:
            logger.info(
                f"Filtering for group_analyzed: {group_analyzed} for mag_id '{mag_id}'"
            )
            # Ensure the column doesn't have mixed types that would break comparison
            sig_sites_df = sig_sites_df[
                sig_sites_df["group_analyzed"] == group_analyzed
            ].copy()
            logger.info(
                f"After group filtering, {len(sig_sites_df)} sites remain for mag_id '{mag_id}'."
            )
        else:
            logger.warning(
                "`--group-analyzed` provided but 'group_analyzed' column not found. Skipping group filtering."
            )

    # Apply significance threshold filtering
    logger.info(
        f"Applying significance threshold: {p_value_column} <= {p_value_threshold}"
    )
    sig_sites_df = sig_sites_df[
        sig_sites_df[p_value_column] <= p_value_threshold
    ].copy()
    logger.info(
        f"After significance filtering ({p_value_column} <= {p_value_threshold}), "
        f"{len(sig_sites_df)} significant sites remain"
    )

    if sig_sites_df.empty:
        logger.warning(
            f"No sites meet the significance threshold ({p_value_column} <= {p_value_threshold}) "
            "and other filter criteria."
        )
        return None

    logger.info(
        "Pre-filtering sites to include only those with a valid, single gene_id."
    )
    sites_to_process = sig_sites_df.dropna(subset=["gene_id"]).copy()
    sites_to_process["gene_id"] = sites_to_process["gene_id"].astype(str).str.strip()
    sites_to_process = sites_to_process[~sites_to_process["gene_id"].str.contains(",")]

    logger.info(
        f"After pre-filtering, {len(sites_to_process)} sites remain with valid single gene_id for mag_id '{mag_id}'."
    )
    if sites_to_process.empty:
        logger.warning("No processable sites found after all filtering steps.")
        return None
    return sites_to_process


def get_major_allele(row, ref_base):
    """
    Determines the major allele from allele count data in a DataFrame row.

    If there's a tie in allele counts, it prioritizes the reference base.
    If the reference base is not part of the tie, it chooses one randomly.

    Args:
        row (pd.Series): A row from the profile DataFrame containing A, T, G, C counts.
        ref_base (str): The reference allele at this position.

    Returns:
        The single-character major allele (e.g., 'A').
    """
    # Create a dictionary of allele counts from the row.
    allele_counts = {
        "A": row.get("A", 0),
        "T": row.get("T", 0),
        "G": row.get("G", 0),
        "C": row.get("C", 0),
    }
    # Find the highest allele count.
    max_count = max(allele_counts.values())
    # If there's no coverage, default to the reference base.
    if max_count == 0:
        return ref_base
    # Find all alleles that have the maximum count.
    major_alleles = [
        allele for allele, count in allele_counts.items() if count == max_count
    ]
    # If there's only one major allele, return it.
    if len(major_alleles) == 1:
        return major_alleles[0]
    else:
        # In case of a tie, prefer the reference base if it's among the tied alleles.
        # Otherwise, break the tie randomly to avoid systematic bias.
        return ref_base if ref_base in major_alleles else random.choice(major_alleles)


def reconstruct_ancestral_sequences(unique_genes, prodigal_records, profile_by_gene):
    """
    Reconstructs ancestral sequences for a list of genes using profile data.

    For each gene, it starts with the reference sequence from Prodigal and then
    substitutes the major allele from the ancestral profile at variable sites.

    Args:
        unique_genes (list): A list of unique gene IDs to reconstruct.
        prodigal_records (dict): A dictionary containing parsed info for each gene.
        profile_by_gene (pd.DataFrameGroupBy): The ancestral profile data, grouped by gene_id.

    Returns:
        A tuple containing:
        - ancestral_orfs (dict): Maps gene_id to its reconstructed sequence string.
        - ancestral_major_alleles (dict): Maps (contig, position) to the major allele on the forward strand.
    """
    ancestral_orfs = {}
    ancestral_major_alleles = {}
    # Iterate through each unique gene that needs reconstruction.
    for gene_id in unique_genes:
        # Skip if the gene ID from the sites file isn't in the annotation file.
        if gene_id not in prodigal_records:
            logger.warning(
                f"Gene ID '{gene_id}' not found in prodigal FASTA, skipping reconstruction."
            )
            continue

        gene_info = prodigal_records[gene_id]
        # Start with the original reference sequence from the Prodigal FASTA file.
        ancestral_seq_list = list(gene_info["record"].seq)

        # Check if this gene has any sites in the ancestral profile data.
        if gene_id in profile_by_gene.groups:
            gene_profile = profile_by_gene.get_group(gene_id)
            # Iterate through each variable site for this gene.
            for _, row in gene_profile.iterrows():
                contig, contig_position = row["contig"], int(row["position"])
                # Get the position of this site within the gene's own coordinate system.
                _, pos_in_gene = get_codon_from_site(
                    contig_position, gene_info, ancestral_seq_list
                )

                if pos_in_gene is not None:
                    # This is the reference base from the original Prodigal FASTA file
                    ref_base = ancestral_seq_list[pos_in_gene]
                    # SANITY CHECK Make sure that the Reference base in profile is same as the one in prodigal
                    ref_base_from_profile = row["ref_base"].upper()
                    # The ref_base from prodigal needs to be oriented to the forward strand for comparison
                    ref_base_from_prodigal_fwd = (
                        COMPLEMENT_MAP.get(ref_base.upper())
                        if gene_info["strand"] == -1
                        else ref_base.upper()
                    )

                    if ref_base_from_prodigal_fwd != ref_base_from_profile:
                        logger.warning(
                            f"Reference base mismatch at {contig_position} in {gene_id}. "
                            f"Prodigal (fwd strand): '{ref_base_from_prodigal_fwd}', "
                            f"Profile: '{ref_base_from_profile}'. Using Prodigal base for reconstruction."
                        )
                    # Determine major allele, using the profile base as the tie-breaker
                    major_allele = get_major_allele(row, ref_base_from_profile)

                    # Store the forward-strand major allele for later comparison.
                    # This is simply the major allele from the profile, as it's always on the forward strand.
                    ancestral_major_allele_fwd = major_allele.upper()
                    ancestral_major_alleles[(contig, contig_position)] = (
                        ancestral_major_allele_fwd
                    )
                    # Determine the "effective allele" to insert into the gene sequence.
                    # For reverse strand genes, this must be the complement of the major allele.
                    effective_allele = (
                        COMPLEMENT_MAP.get(major_allele.upper())
                        if gene_info["strand"] == -1
                        else major_allele
                    )

                    if effective_allele:
                        ancestral_seq_list[pos_in_gene] = effective_allele
                    else:
                        logger.warning(
                            f"Could not determine effective allele for major allele '{major_allele}' in gene {gene_id}."
                        )

                else:
                    logger.warning(
                        f"Contig position {contig_position} for gene {gene_id} is out of bounds. Skipping."
                    )

        ancestral_orfs[gene_id] = "".join(ancestral_seq_list)
    return ancestral_orfs, ancestral_major_alleles


# --- Helper and Core dN/dS Analysis Functions ---
def get_codon_from_site(position, gene_info, sequence):
    """
    Helper to calculate a site's position within a gene and extract its codon.

    This function handles coordinate conversion based on whether the gene is on
    the forward (+) or reverse (-) strand.
    get_codon_from_site

    COORDINATE SYSTEM EXPLANATION:
    - Input position: 0-based position on the contig (from frequency file)
    - Prodigal coordinates: 1-based start/end positions on the contig
    - Gene sequence: Always 5'->3' coding sequence, regardless of strand

    For FORWARD genes: position in gene = SNV_pos - (gene_start - 1)
    For REVERSE genes: position in gene = (gene_end - 1) - SNV_pos

    STRAND HANDLING:
    - Forward strand: Use alleles directly from frequency file
    - Reverse strand: Complement alleles to match reverse-complement sequence

    EXAMPLE:
      - Contig is 10000 bp long.
      - SNV Position (0-based from freq file): 8025

      - FORWARD GENE EXAMPLE:
        - Prodigal Header: # 8000 # 8500 # 1
        - The gene sequence starts at contig position 8000.
        - Calculation: 8025 (SNV) - (8000 - 1) (gene start) = index 26
        - This is intuitive: the SNV is 26 bases into the gene.

      - REVERSE GENE EXAMPLE:
        - Prodigal Header: # 8000 # 8500 # -1
        - The provided gene sequence is the REVERSE COMPLEMENT. Its first base
          corresponds to contig position 8500.
        - Calculation: (8500 - 1) (gene end) - 8025 (SNV) = index 474
        - This is counter-intuitive but correct: the SNV is near the 'start'
          coordinate on the contig, which means it's at the *end* of the
          reverse-complemented sequence string.

    Args:
        position (int): The genomic position on the contig.
        gene_info (dict): Metadata for the gene (start, end, strand).
        sequence (list or str): The gene's sequence.

    Returns:
        A tuple of (codon_string, pos_in_gene) or (None, None) if the position is invalid.
    """
    strand = gene_info["strand"]
    # Calculate the 0-indexed position within the gene's sequence.
    if strand == 1:
        pos_in_gene = position - (gene_info["start"] - 1)
    elif strand == -1:
        pos_in_gene = (gene_info["end"] - 1) - position
    else:
        logger.warning(
            f"Unknown strand '{strand}' for gene {gene_info['record'].id}. Cannot get codon."
        )
        return None, None
    # Validate that the calculated position is within the bounds of the sequence.
    if not (0 <= pos_in_gene < len(sequence)):
        logger.warning(
            f"Calculated position {pos_in_gene} for gene {gene_info['record'].id} is out of sequence bounds (gene len: {len(sequence)})."
        )
        return None, None

    # Find the start of the 3-base codon containing this position.
    codon_start_index = (pos_in_gene // 3) * 3
    # Handle both string and list sequence types
    codon_slice = sequence[codon_start_index : codon_start_index + 3]
    codon = "".join(codon_slice) if isinstance(codon_slice, list) else codon_slice

    # Ensure a full codon was extracted.
    if len(codon) != 3:
        logger.warning(
            f"Site at contig pos {position} in gene {gene_info['record'].id} results in a partial codon ('{codon}')."
        )
        return None, None

    return codon, pos_in_gene


def analyze_mutation_effect(gene_info, ancestral_seq, position, allele_after):
    """
    Determines if a mutation is synonymous (S) or non-synonymous (NS).

    It compares the amino acid translated from the ancestral codon to the one
    translated from the derived (mutated) codon.

    Args:
        gene_info (dict): Metadata for the gene.
        ancestral_seq (str): The reconstructed ancestral gene sequence.
        position (int): The genomic position of the mutation.
        allele_after (str): The derived major allele (on the forward strand).

    Returns:
        A dictionary with details about the mutation, or an empty dict if analysis fails.
    """
    ancestral_codon_str, pos_in_gene = get_codon_from_site(
        position, gene_info, ancestral_seq
    )
    if ancestral_codon_str is None:
        return {}

    strand = gene_info["strand"]
    # Alleles from profiles are from the forward strand. They must be complemented for reverse-strand genes.
    effective_allele_after = (
        COMPLEMENT_MAP.get(allele_after.upper())
        if strand == -1
        else allele_after.upper()
    )

    if not effective_allele_after:
        logger.warning(
            f"Cannot complement derived allele for gene {gene_info['record'].id}: {allele_after}"
        )
        return {}

    # Find the position within the 3-base codon (0, 1, or 2).
    pos_in_codon = pos_in_gene % 3

    # The "after" codon is the ancestral codon with the derived allele substituted.
    mutable_codon_after = MutableSeq(ancestral_codon_str)
    mutable_codon_after[pos_in_codon] = effective_allele_after
    derived_codon_seq = Seq(mutable_codon_after)

    # Translate both codons to determine the effect.
    aa_before = Seq(ancestral_codon_str).translate(table=11, cds=False)
    aa_after = derived_codon_seq.translate(table=11, cds=False)
    mutation_type = "S" if aa_before == aa_after else "NS"

    return {
        "codon_before": ancestral_codon_str,
        "codon_after": str(derived_codon_seq),
        "aa_before": str(aa_before),
        "aa_after": str(aa_after),
        "mutation_type": mutation_type,
    }


def find_substitutions(
    sites_df,
    ancestral_seqs,
    prodigal_records,
    derived_profile_df,
    ancestral_major_alleles,
):
    """
    Identifies and analyzes substitutions by comparing ancestral and derived major alleles.

    Args:
        sites_df (pd.DataFrame): DataFrame of significant sites to check.
        ancestral_seqs (dict): Reconstructed ancestral sequences.
        prodigal_records (dict): Gene metadata.
        derived_profile_df (pd.DataFrame): Profile data for the derived sample.
        ancestral_major_alleles (dict): Pre-calculated major alleles for the ancestral sample.

    Returns:
        A list of dictionaries, where each dictionary details a confirmed substitution.
    """
    results = []
    # Indexing the derived profile allows for much faster lookups.
    derived_profile_indexed = derived_profile_df.set_index(
        ["contig", "position", "gene_id"]
    )

    for _, site_row in sites_df.iterrows():
        contig, position, gene_id = (
            site_row["contig"],
            site_row["position"],
            str(site_row["gene_id"]).strip(),
        )

        if gene_id not in prodigal_records:
            logger.warning(
                f"Gene ID '{gene_id}' from sites file at {contig}:{position} not found in records. Skipping."
            )
            continue

        gene_info = prodigal_records[gene_id]
        ancestral_seq = ancestral_seqs[gene_id]

        if ancestral_seq is None:
            logger.warning(
                f"Ancestral sequence for gene '{gene_id}' not found in derived profile results. Skipping site {contig}:{position}."
            )
            continue
        try:
            derived_site_row = derived_profile_indexed.loc[(contig, position, gene_id)]
        except KeyError:
            logger.warning(
                f"Site {contig}:{position} for gene {gene_id} not found in the derived profile. Skipping."
            )
            continue

        # 1. Get the ANCESTRAL major allele from the pre-calculated dictionary.
        ancestral_major_allele_fwd = ancestral_major_alleles.get((contig, position))
        if ancestral_major_allele_fwd is None:
            logger.warning(
                f"Ancestral site {contig}:{position} not found in ancestral profile results. Skipping."
            )
            continue

        # 2. Determine the DERIVED major allele using the derived_site_row as reference.
        profile_ref_base_fwd = derived_site_row["ref_base"].upper()
        # 2.1 SANITY CHECK Make sure that the Reference base in profile is same as the one in prodigal
        original_ref_seq = str(gene_info["record"].seq)
        # Get position within the gene's coordinate system (0-indexed)
        _, pos_in_gene = get_codon_from_site(position, gene_info, original_ref_seq)
        prodigal_base_on_strand = original_ref_seq[pos_in_gene]
        # The ref_base from prodigal needs to be oriented to the forward strand for comparison
        prodigal_ref_base_fwd = (
            COMPLEMENT_MAP.get(prodigal_base_on_strand.upper())
            if gene_info["strand"] == -1
            else prodigal_base_on_strand.upper()
        )

        if profile_ref_base_fwd != prodigal_ref_base_fwd:
            logger.warning(
                f"Reference base mismatch at {contig}:{position} in {gene_id}. "
                f"Prodigal says: '{prodigal_ref_base_fwd}', "
                f"Profile file says: '{profile_ref_base_fwd}'. "
                f"Using profile 'ref_base' for tie-breaking."
            )
        # 2.2 Determine the DERIVED major allele using the derived_site_row as reference.
        derived_major_allele_fwd = get_major_allele(
            derived_site_row, profile_ref_base_fwd
        )

        # 3. A substitution is a site where the major allele has CHANGED BETWEEN time points.
        if (
            ancestral_major_allele_fwd
            and derived_major_allele_fwd
            and ancestral_major_allele_fwd != derived_major_allele_fwd
        ):
            mutation_info = analyze_mutation_effect(
                gene_info,
                ancestral_seq,
                position,
                derived_major_allele_fwd,
            )
            if mutation_info:
                results.append(
                    {
                        "mag_id": site_row["mag_id"],
                        "contig": contig,
                        "position": position,
                        "gene_id": gene_id,
                        "pos_in_gene": pos_in_gene,
                        "codon_start_index": (pos_in_gene // 3) * 3,
                        "position_in_codon": pos_in_gene % 3,
                        "strand": gene_info["strand"],
                        "major_allele_ancestral": ancestral_major_allele_fwd,
                        "major_allele_derived": derived_major_allele_fwd,
                        "potential_N_sites": gene_info.get("potential_N_sites", np.nan),
                        "potential_S_sites": gene_info.get("potential_S_sites", np.nan),
                        **mutation_info,
                    }
                )

    return results


# --- Global dN/dS Calculation ---
def calculate_global_dnds_for_sites(
    substitution_results_df, ancestral_seqs, prodigal_records
):
    """
    Calculates a single, global dN/dS ratio for a specific set of sites.

    This method correctly avoids double-counting potential sites by tracking unique
    codons involved in substitutions.

    Args:
        substitution_results_df (pd.DataFrame): DataFrame of confirmed substitutions.
        ancestral_seqs (dict): Reconstructed ancestral sequences.
        prodigal_records (dict): Gene metadata.

    Returns:
        A tuple containing:
        - summary_df (pd.DataFrame): Global dN/dS calculation summary.
        - codon_specific_df (pd.DataFrame): Table with codon-specific potential S/N sites for each significant site.
    """
    logger.info(
        f"Calculating global dN/dS for {len(substitution_results_df):,} significant sites..."
    )
    if substitution_results_df.empty:
        logger.warning("No substitutions were found, cannot calculate global dN/dS.")
        return None, None

    total_potential_S = 0.0
    total_potential_N = 0.0

    # Use a set to track unique codons that have been processed to avoid double-counting.
    # The key will be a tuple: (gene_id, codon_start_index)
    processed_codons = set()

    # Store detailed information for each site to create the codon-specific table
    codon_specific_data = []

    # 1. Calculate the total potential S and N sites for the UNIQUE codons of interest
    # Only iterate over the sites where substitutions actually happened.
    for _, site_row in substitution_results_df.iterrows():
        gene_id = str(site_row["gene_id"]).strip()
        position = site_row["position"]

        # Ensure data is available for this gene.
        if gene_id not in prodigal_records:
            logger.warning(
                f"Gene ID '{gene_id}' at site {site_row['contig']}:{position} not found in records; skipping for global dN/dS."
            )
            continue

        gene_info = prodigal_records[gene_id]
        ancestral_seq = ancestral_seqs[gene_id]

        # get the ancestral codon
        codon, pos_in_gene = get_codon_from_site(position, gene_info, ancestral_seq)

        if codon:
            codon_start_index = (pos_in_gene // 3) * 3
            codon_identifier = (gene_id, codon_start_index)

            # Get codon-specific potential sites
            codon_sites = _CODON_SITE_CACHE.get(codon.upper())
            if codon_sites is None:
                logger.warning(
                    f"Codon '{codon}' at site {site_row['contig']}:{position} not in cache; skipping for global dN/dS."
                )
                continue

            s_codon, n_codon = codon_sites

            # Add detailed information for this site to the codon-specific table
            codon_specific_data.append(
                {
                    "mag_id": site_row["mag_id"],
                    "contig": site_row["contig"],
                    "position": position,
                    "gene_id": gene_id,
                    "pos_in_gene": pos_in_gene,
                    "codon_start_index": codon_start_index,
                    "position_in_codon": pos_in_gene % 3,
                    "strand": site_row["strand"],
                    "major_allele_ancestral": site_row["major_allele_ancestral"],
                    "major_allele_derived": site_row["major_allele_derived"],
                    "codon_before": site_row["codon_before"],
                    "codon_after": site_row["codon_after"],
                    "aa_before": site_row["aa_before"],
                    "aa_after": site_row["aa_after"],
                    "mutation_type": site_row["mutation_type"],
                    "potential_S_sites_codon": s_codon,
                    "potential_N_sites_codon": n_codon,
                }
            )

            # Only process this codon if we haven't seen it before for global totals
            # IMP: This is useful if there is more than one significant site on the same codon.
            # The (gene_id, codon_start_index) is unique for each codon in the gene.
            if codon_identifier not in processed_codons:
                total_potential_S += s_codon
                total_potential_N += n_codon

                # Add this codon to the set of processed codons to prevent re-counting
                processed_codons.add(codon_identifier)

    # 2. Count the observed S and NS substitutions from the analysis results
    observed_S = (substitution_results_df["mutation_type"] == "S").sum()
    observed_N = (substitution_results_df["mutation_type"] == "NS").sum()

    # 3. Calculate pN, pS, and the dN/dS ratio
    pN = observed_N / total_potential_N if total_potential_N > 0 else 0
    pS = observed_S / total_potential_S if total_potential_S > 0 else 0
    dNdS_ratio = pN / pS if pS > 0 else np.nan

    # 4. Create a summary DataFrame
    summary_data = {
        "Metric": [
            "Total Potential Non-Synonymous Sites (n)",
            "Total Potential Synonymous Sites (s)",
            "Observed Non-Synonymous Substitutions (NS)",
            "Observed Synonymous Substitutions (S)",
            "pN (NS / n)",
            "pS (S / s)",
            "Global dN/dS (pN/pS)",
        ],
        "Value": [
            f"{total_potential_N:.2f}",
            f"{total_potential_S:.2f}",
            observed_N,
            observed_S,
            f"{pN:.4f}",
            f"{pS:.4f}",
            f"{dNdS_ratio:.4f}",
        ],
    }
    summary_df = pd.DataFrame(summary_data)
    logger.info(
        "\n--- Global dN/dS Summary for Significant Sites ---\n"
        + summary_df.to_string(index=False)
    )

    # 5. Create the codon-specific DataFrame
    codon_specific_df = pd.DataFrame(codon_specific_data)

    return summary_df, codon_specific_df


# --- Summarization and Output Functions ---
def summarize_results(results_df):
    """Generates gene-level and MAG-level summary tables with dN/dS ratios."""
    if results_df.empty:
        logger.warning("Results dataframe is empty. Cannot generate summaries.")
        return {}

    # --- Gene-level Statistics ---
    logger.info("Summarizing results at the gene level...")
    gene_mutation_counts = (
        results_df.groupby(["mag_id", "gene_id"])
        .agg(
            s_count=("mutation_type", lambda s: (s == "S").sum()),
            ns_count=("mutation_type", lambda s: (s == "NS").sum()),
        )
        .reset_index()
    )

    # Get the potential S/N site counts for each gene (calculated earlier).
    gene_site_counts = results_df[
        ["mag_id", "gene_id", "potential_N_sites", "potential_S_sites"]
    ].drop_duplicates()
    gene_stats = pd.merge(
        gene_mutation_counts, gene_site_counts, on=["mag_id", "gene_id"]
    )

    # dN = (# non-synonymous substitutions) / (# potential non-synonymous sites)
    gene_stats["dN"] = gene_stats.apply(
        lambda r: (
            r["ns_count"] / r["potential_N_sites"] if r["potential_N_sites"] > 0 else 0
        ),
        axis=1,
    )
    # dS = (# synonymous substitutions) / (# potential synonymous sites)
    gene_stats["dS"] = gene_stats.apply(
        lambda r: (
            r["s_count"] / r["potential_S_sites"] if r["potential_S_sites"] > 0 else 0
        ),
        axis=1,
    )
    # dN/dS ratio. Return NaN if dS is 0 to avoid division by zero errors.
    gene_stats["dN_dS_ratio"] = gene_stats.apply(
        lambda r: r["dN"] / r["dS"] if r["dS"] > 0 else np.nan, axis=1
    )

    # --- MAG-level Statistics ---
    logger.info("Summarizing results at the MAG level...")
    mag_stats = (
        gene_stats.groupby("mag_id")
        .agg(
            total_ns_count=("ns_count", "sum"),
            total_s_count=("s_count", "sum"),
            potential_N_sites=("potential_N_sites", "sum"),
            potential_S_sites=("potential_S_sites", "sum"),
        )
        .reset_index()
    )

    mag_stats["dN"] = mag_stats.apply(
        lambda r: (
            r["total_ns_count"] / r["potential_N_sites"]
            if r["potential_N_sites"] > 0
            else 0
        ),
        axis=1,
    )
    mag_stats["dS"] = mag_stats.apply(
        lambda r: (
            r["total_s_count"] / r["potential_S_sites"]
            if r["potential_S_sites"] > 0
            else 0
        ),
        axis=1,
    )
    mag_stats["dN_dS_ratio"] = mag_stats.apply(
        lambda r: r["dN"] / r["dS"] if r["dS"] > 0 else np.nan, axis=1
    )

    return {"gene": gene_stats, "mag": mag_stats}


def write_output_files(
    summaries,
    results_df,
    args,
    ancestral_sequences,
    prodigal_records,
    global_dnds_summary_df,
    codon_specific_df,
):
    """
    Writes the main events table, summary tables, and ancestral sequences to disk.
    Skips writing files if the corresponding dataframe is empty or None.
    """
    args.outdir.mkdir(parents=True, exist_ok=True)

    # --- All substitutions ---
    path_all_events = args.outdir / f"{args.prefix}_all_substitutions.tsv"
    if results_df is not None and not results_df.empty:
        results_df.to_csv(path_all_events, sep="\t", index=False)
        logger.info(
            f"Successfully wrote {len(results_df)} substitution events to {path_all_events}"
        )
    else:
        logger.warning(f"No substitutions found. Skipping file: {path_all_events}")

    # --- Gene and MAG summaries ---
    if summaries:
        for summary_type in ["gene", "mag"]:
            path = args.outdir / f"{args.prefix}_{summary_type}_summary.tsv"
            df = summaries.get(summary_type)
            if df is not None and not df.empty:
                df.to_csv(path, sep="\t", index=False)
                logger.info(f"Wrote {summary_type} summary to {path}")
            else:
                logger.warning(
                    f"No {summary_type} summary was generated. Skipping file: {path}"
                )
    else:
        logger.warning(
            "No summaries were generated. Skipping gene and MAG summary files."
        )

    # --- Global dN/dS summary ---
    path_global_summary = args.outdir / f"{args.prefix}_global_dnds_summary.tsv"
    if global_dnds_summary_df is not None and not global_dnds_summary_df.empty:
        global_dnds_summary_df.to_csv(path_global_summary, sep="\t", index=False)
        logger.info(f"Wrote global dN/dS summary to {path_global_summary}")
    else:
        logger.warning(
            f"No global dN/dS summary was generated. Skipping file: {path_global_summary}"
        )

    # --- Codon-specific sites ---
    path_codon_specific = args.outdir / f"{args.prefix}_codon_specific_sites.tsv"
    if codon_specific_df is not None and not codon_specific_df.empty:
        codon_specific_df.to_csv(path_codon_specific, sep="\t", index=False)
        logger.info(
            f"Wrote codon-specific potential sites table with {len(codon_specific_df)} entries to {path_codon_specific}"
        )
    else:
        logger.warning(
            f"No codon-specific data was generated. Skipping file: {path_codon_specific}"
        )

    # --- Reconstructed ancestral sequences ---
    path_ancestral_fasta = args.outdir / f"{args.prefix}_ancestral_orfs.ffn"
    logger.info(f"Writing ancestral sequences to {path_ancestral_fasta}...")
    if ancestral_sequences:
        ancestral_records = []
        for gene_id, seq_str in ancestral_sequences.items():
            if gene_id in prodigal_records:
                original_description = prodigal_records[gene_id]["record"].description
                rest_of_description = original_description.replace(
                    gene_id, "", 1
                ).strip()
                record = SeqRecord(
                    Seq(seq_str),
                    id=f"ancestral_{gene_id}",
                    description=rest_of_description,
                )
                ancestral_records.append(record)
        with open(path_ancestral_fasta, "w") as output_handle:
            SeqIO.write(ancestral_records, output_handle, "fasta")
        logger.info(f"Successfully wrote {len(ancestral_records)} ancestral sequences.")
    else:
        logger.warning(
            f"No ancestral sequences to write. Skipping file: {path_ancestral_fasta}"
        )


def process_mag(mag_id, args, prodigal_records):
    """Process a single MAG and write output files."""
    # logger.info(f"--- Processing MAG: {mag_id} ---")

    # Create a copy of the args object to avoid modifying the original
    mag_args = argparse.Namespace(**vars(args))
    # Update the prefix and outdir for this specific MAG
    mag_args.prefix = f"{args.prefix}_{mag_id}_{args.test_type}_{args.p_value_column}"
    mag_args.outdir = args.outdir / mag_id

    # 1. Load and filter significant sites for this MAG
    sites_to_process = setup_and_load_data(
        args.significant_sites,
        mag_id,
        mag_args.p_value_column,
        mag_args.p_value_threshold,
        mag_args.test_type,
        mag_args.group_analyzed,
    )
    if sites_to_process is None or sites_to_process.empty:
        #     logger.warning(
        #         f"No processable sites found for MAG {mag_id} after filtering. Writing empty files and skipping."
        #     )
        #     write_output_files(None, None, mag_args, None, None, None, None)
        return

    # 2. Load profile data for this MAG
    ancestral_profile_path = (
        mag_args.profile_dir
        / f"{mag_args.ancestral_sample_id}/{mag_args.ancestral_sample_id}_{mag_id}_profiled.tsv.gz"
    )
    derived_profile_path = (
        mag_args.profile_dir
        / f"{mag_args.derived_sample_id}/{mag_args.derived_sample_id}_{mag_id}_profiled.tsv.gz"
    )

    if not ancestral_profile_path.exists() or not derived_profile_path.exists():
        logger.warning(
            f"Profile file(s) for MAG {mag_id} not found at {ancestral_profile_path} or {derived_profile_path}. Skipping."
        )
        return

    ancestral_profile_df = pd.read_csv(ancestral_profile_path, sep="\t")
    derived_profile_df = pd.read_csv(derived_profile_path, sep="\t")

    # 3. Reconstruct Ancestral Sequences
    logger.info(
        f"Reconstructing Ancestral Sequences for {mag_id} (genes: {sites_to_process['gene_id'].nunique()})"
    )
    ancestral_profile_df.dropna(subset=["gene_id"], inplace=True)
    ancestral_profile_df["gene_id"] = (
        ancestral_profile_df["gene_id"].astype(str).str.strip()
    )
    profile_by_gene = ancestral_profile_df.groupby("gene_id")
    unique_genes = sites_to_process["gene_id"].unique()
    ancestral_sequences, ancestral_major_alleles = reconstruct_ancestral_sequences(
        unique_genes, prodigal_records, profile_by_gene
    )

    if not ancestral_sequences:
        logger.warning(f"No ancestral sequences reconstructed for MAG {mag_id}.")

    # 4. Calculate potential S/N sites for each entire ancestral gene
    for gene_id, seq_str in ancestral_sequences.items():
        if gene_id in prodigal_records:
            potential_sites = calculate_potential_sites_for_gene(seq_str)
            prodigal_records[gene_id]["potential_S_sites"] = potential_sites["S"]
            prodigal_records[gene_id]["potential_N_sites"] = potential_sites["N"]

    # 5. Find substitutions
    derived_profile_df.dropna(subset=["gene_id"], inplace=True)
    derived_profile_df["gene_id"] = (
        derived_profile_df["gene_id"].astype(str).str.strip()
    )
    logger.info(f"Analyzing substitutions for {mag_id}...")
    substitution_results = find_substitutions(
        sites_to_process,
        ancestral_sequences,
        prodigal_records,
        derived_profile_df,
        ancestral_major_alleles,
    )
    results_df = (
        pd.DataFrame(substitution_results) if substitution_results else pd.DataFrame()
    )

    # 6. Summarize results
    summaries = summarize_results(results_df)

    # 7. Calculate global dN/dS
    global_dnds_summary_df, codon_specific_df = calculate_global_dnds_for_sites(
        results_df, ancestral_sequences, prodigal_records
    )

    # 8. Write outputs
    write_output_files(
        summaries,
        results_df,
        mag_args,
        ancestral_sequences,
        prodigal_records,
        global_dnds_summary_df,
        codon_specific_df,
    )


def main():
    """Main execution function to orchestrate the entire analysis."""
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Analyze dN/dS from a reconstructed ancestral state to a derived state.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--significant_sites",
        required=True,
        type=Path,
        help="Path to the significant sites dataframe (TSV) from p_value_summary.py. "
        "Expected columns: mag_id, contig, position, gene_id, test_type, min_p_value, q_value.",
    )
    parser.add_argument(
        "--mag_ids",
        nargs="+",
        help="One or more MAG IDs to process.",
    )
    parser.add_argument(
        "--p_value_column",
        choices=["min_p_value", "q_value"],
        default="q_value",
        help="Column to use for significance filtering: 'min_p_value' (raw p-values) or 'q_value' (FDR-corrected).",
    )
    parser.add_argument(
        "--p_value_threshold",
        type=float,
        default=0.05,
        help="Threshold for significance filtering using the selected p-value column.",
    )
    parser.add_argument(
        "--test-type",
        type=str,
        choices=[
            "two_sample_unpaired_tTest",
            "two_sample_unpaired_MannWhitney",
            "two_sample_paired_tTest",
            "two_sample_paired_Wilcoxon",
            "single_sample_tTest",
            "single_sample_Wilcoxon",
            "cmh",
            "lmm",
            "lmm_across_time",
            "cmh_across_time",
        ],
        default="two_sample_paired_tTest",
        help="A test type to filter the `significant_sites` dataframe.",
    )
    parser.add_argument(
        "--group-analyzed",
        type=str,
        help="Filter the `significant_sites` dataframe for a specific group, if the 'group_analyzed' column is present.",
    )
    parser.add_argument(
        "--ancestral_sample_id",
        required=True,
        help="The sample ID for the ANCESTRAL timepoint.",
    )
    parser.add_argument(
        "--derived_sample_id",
        required=True,
        help="The sample ID for the DERIVED timepoint for comparison.",
    )
    parser.add_argument(
        "--profile_dir",
        required=True,
        type=Path,
        help="Directory containing the profile files from profile_mags.py.",
    )
    parser.add_argument(
        "--prodigal_fasta",
        required=True,
        type=Path,
        help="Path to the Prodigal-predicted FASTA file for gene sequences and coordinates.",
    )
    parser.add_argument(
        "--outdir", required=True, type=Path, help="Directory to save the output files."
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="ancestral_dnds",
        help="Prefix for the output file names.",
    )
    parser.add_argument(
        "--cpus",
        type=int,
        default=mp.cpu_count(),
        help="Number of CPUs to use for multiprocessing. Default: use all available CPUs.",
    )
    args = parser.parse_args()

    # Load Prodigal records once for all MAGs
    logger.info("Loading gene annotation file...")
    prodigal_raw = SeqIO.to_dict(SeqIO.parse(args.prodigal_fasta, "fasta"))
    prodigal_records = {}
    for rec in prodigal_raw.values():
        header_parts = [p.strip() for p in rec.description.split("#")]
        prodigal_records[rec.id] = {
            "record": rec,
            "start": int(header_parts[1]),
            "end": int(header_parts[2]),
            "strand": int(header_parts[3]),
        }

    # Determine number of CPUs to use, ensuring it doesn't exceed the number of MAGs
    num_cpus = min(args.cpus, len(args.mag_ids))
    logger.info(f"Using {num_cpus} CPUs to process {len(args.mag_ids)} MAGs")

    # Process each MAG using multiprocessing
    sorted_mag_ids = sorted(args.mag_ids)

    # Create a partial function with fixed args and prodigal_records
    process_mag_partial = partial(
        process_mag, args=args, prodigal_records=prodigal_records
    )

    # Use multiprocessing with imap_unordered and tqdm
    with mp.Pool(processes=num_cpus) as pool:
        # Use imap_unordered to get results as they complete
        list(
            tqdm(
                pool.imap_unordered(process_mag_partial, sorted_mag_ids),
                total=len(sorted_mag_ids),
                desc="Processing MAGs",
            )
        )

    logger.info("Analysis completed successfully for all processed MAGs!")


if __name__ == "__main__":
    main()
