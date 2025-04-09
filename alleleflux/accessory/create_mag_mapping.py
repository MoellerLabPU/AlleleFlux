#!/usr/bin/env python3
"""
create_mag_mapping.py

A standalone preprocessing tool that should be run BEFORE the AlleleFlux workflow.

This script processes MAG FASTA files by:
1. Modifying contig headers to prepend the MAG filename
2. Concatenating all modified FASTAs into one output file
3. Generating a mapping file of scaffold IDs to MAG IDs

The output mapping file can be used directly in the AlleleFlux workflow
by specifying it in the config.yml file under input.mag_mapping_path.

Example usage:
    alleleflux-create-mag-mapping --dir /path/to/mags/ --output-fasta combined.fasta --output-mapping mag_mapping.tsv

After running this tool:
1. Set the path to the generated FASTA file in config.yml under input.fasta_path
2. Set the path to the generated mapping file in config.yml under input.mag_mapping_path
3. Run prodigal to generate gene predictions for the combined FASTA file
4. Run Bowtie2 to align the reads to the combined FASTA file
5. Run the AlleleFlux workflow with the updated config.yml
"""

import argparse
import logging
import os
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def find_fasta_files(directory, extension):
    """
    Find all FASTA files with the specified extension in a directory.

    Parameters:
        directory (str): Path to directory to search
        extension (str): File extension to match (without leading dot)

    Returns:
        list: List of Path objects for matching files
    """
    # Normalize extension
    ext = extension.lstrip(".")
    return sorted(Path(directory).glob(f"*.{ext}"))


def process_and_concatenate_fastas(input_files, output_fasta, mapping_file, extension):
    """
    Process each FASTA file, modify headers, and concatenate into one file.
    Also generates a mapping file.

    Parameters:
        input_files (list): List of input FASTA file paths
        output_fasta (str): Path to output concatenated FASTA file
        mapping_file (str): Path to output mapping file

    Returns:
        tuple: (total_seqs, total_mags) - Count of sequences and MAGs processed
    """

    # Create output directory if it doesn't exist
    output_fasta = Path(output_fasta)
    mapping_file = Path(mapping_file)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    total_seqs = 0
    total_mags = 0
    mag_ids = set()

    # Process each file and write to the output
    with open(output_fasta, "w") as out_fasta, open(mapping_file, "w") as out_map:
        # Write header to mapping file
        out_map.write("mag_id\tcontig_id\n")

        for mag_file in input_files:
            mag_basename = mag_file.name.replace(f".{extension}", "")
            # Check for duplicate MAG IDs
            if mag_basename in mag_ids:
                raise ValueError(
                    f"Duplicate MAG ID found: {mag_basename}. Skipping file."
                )
            mag_ids.add(mag_basename)

            logging.info(f"Processing {mag_file}")

            with mag_file.open("r") as handle:
                seqs_in_file = 0
                for record in SeqIO.parse(handle, "fasta"):
                    original_id = record.id
                    new_id = f"{mag_basename}_{original_id}"

                    new_record = SeqRecord(
                        record.seq,
                        id=new_id,
                        description="",
                    )
                    SeqIO.write(new_record, out_fasta, "fasta")
                    out_map.write(f"{mag_basename}\t{new_id}\n")

                    seqs_in_file += 1
                    total_seqs += 1

            total_mags += 1
            logging.info(f"Added {seqs_in_file} sequences from {mag_basename}")

    logging.info(f"Processing complete: {total_seqs} sequences from {total_mags} MAGs")
    return total_seqs, total_mags


def main():
    logging.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logging.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="""
        PREPROCESSING UTILITY: Create a combined FASTA file and MAG mapping file from individual MAG FASTA files.
        This tool should be run BEFORE starting the AlleleFlux workflow.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--dir",
        help="Path to directory containing MAG FASTA files",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--extension",
        help="Exact file extension to identify FASTA files (e.g., 'fasta', 'fa', 'fa.gz')",
        type=str,
        default="fa",
    )

    parser.add_argument(
        "--output-fasta",
        help="Path for the concatenated output FASTA file",
        type=str,
        default="megamag.fasta",
    )

    parser.add_argument(
        "--output-mapping",
        help="Path for the scaffold-to-bin mapping TSV file",
        type=str,
        default="megamag_mapping.tsv",
    )

    args = parser.parse_args()

    # Validate input directory
    if not os.path.isdir(args.dir):
        raise ValueError(f"Directory not found: {args.dir}")

    # Find FASTA files
    fasta_files = find_fasta_files(args.dir, args.extension)

    if not fasta_files:
        raise ValueError(
            f"No FASTA files with extension '{args.extension}' found in {args.dir}"
        )

    logging.info(f"Found {len(fasta_files)} FASTA files to process")

    # Process files
    total_seqs, total_mags = process_and_concatenate_fastas(
        fasta_files, args.output_fasta, args.output_mapping, args.extension
    )

    if total_seqs == 0:
        raise ValueError("No sequences were processed")

    logging.info(f"Created concatenated FASTA file: {args.output_fasta}")
    logging.info(f"Created mapping file: {args.output_mapping}")
    logging.info(f"\nNext steps:")
    logging.info(
        f"1. Set {os.path.abspath(args.output_fasta)} as fasta_path in config.yml"
    )
    logging.info(
        f"2. Set {os.path.abspath(args.output_mapping)} as mag_mapping_path in config.yml"
    )
    logging.info(
        f"3. Run prodigal to generate gene predictions for the combined FASTA file"
    )
    logging.info(f"4. Run Bowtie2 to align the reads to the combined FASTA file")
    logging.info(f"5. Run the AlleleFlux workflow with the updated config.yml")


if __name__ == "__main__":
    main()
