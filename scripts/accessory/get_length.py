import argparse
import logging

import pandas as pd
from Bio import SeqIO


def parse_fasta_to_table(fasta_file, output_csv):

    logging.info(f"Processing FASTA file: {fasta_file}")
    data = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_id = record.id  # Full contig ID
        sequence_length = len(record.seq)  # Length of the sequence
        mag_id = contig_id.split(".fa")[0]  # Extract MAG before ".fa"

        data.append({"contig": contig_id, "length": sequence_length, "MAG": mag_id})

    # Create a DataFrame
    df = pd.DataFrame(data, columns=["contig", "length", "MAG"])

    # Save the DataFrame to a CSV file
    df.to_csv(output_csv, index=False, sep="\t")
    logging.info(f"Table successfully saved to {output_csv}")


def main():
    logging.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logging.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Parse a FASTA file and create a table with contig, length, and MAG."
    )
    parser.add_argument(
        "--fasta", type=str, help="Path to the input FASTA file", required=True
    )
    parser.add_argument(
        "--out_fPath", type=str, help="Path to the output TSV file", required=True
    )

    # Parse the arguments
    args = parser.parse_args()
    parse_fasta_to_table(args.fasta, args.out_fPath)


if __name__ == "__main__":
    main()
