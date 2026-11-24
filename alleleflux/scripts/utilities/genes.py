"""Prodigal gene records and the one piece of codon-index arithmetic.

Two analyses read the same Prodigal gene file and map contig positions onto
coding sequences: dN/dS (``evolution/dnds_from_timepoints.py``) and the gene
de novo tables (``analysis/gene_de_novo.py``).  Before this module each kept
its own copy of the header parser and the strand arithmetic; a strand bug
could have been fixed in one and not the other.  Now both call here.

Coordinate systems, spelled out once
------------------------------------
* **Prodigal header** ``>id # start # end # strand # ...``: ``start`` and
  ``end`` are **1-based inclusive** contig coordinates; ``strand`` is ``1``
  or ``-1``.  ``parse_prodigal_records`` keeps these numbers as written.
* **Profile / cache positions** are **0-based** contig coordinates (samtools
  pileup minus one, as ``profile_mags`` writes them).  ``load_genes`` subtracts
  one from both ends so its rows live in this system, the same rule
  ``profile_mags.map_genes`` applies when it stamps ``gene_id`` on sites.
* **Coding sequence**: the FASTA body is always read 5' to 3' on the coding
  strand.  For a reverse-strand gene Prodigal has already
  reverse-complemented it, so its first letter sits at the gene's **end**
  coordinate on the contig and its letters are the complements of the
  contig's.  ``coding_index`` therefore counts from the end on strand -1, and
  callers complement contig letters with ``COMPLEMENT_MAP`` before comparing.
"""

from __future__ import annotations

import pandas as pd
from Bio import SeqIO

# Complement of each DNA base, for reverse-strand genes.  Only the four
# unambiguous bases: an N at a site is never classified by either consumer.
COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Column order of the frame ``load_genes`` returns (kept even when empty).
GENE_COLUMNS = ["gene_id", "contig", "start", "end", "strand", "length", "sequence"]


def _header_fields(description: str) -> tuple[str, int, int, int]:
    """Split a Prodigal description into ``(gene_id, start, end, strand)``.

    Prodigal separates its fields with ``#``: id, start, end, strand, then a
    ``;``-joined attribute string.  Fields are stripped because the file has
    a space on either side of each ``#``.  Raises ``ValueError`` unless there
    are exactly five ``#``-separated fields with numeric start / end / strand
    (every real Prodigal header has five), naming the offending header.

    Example: ``"c1_1 # 101 # 109 # 1 # ID=1_1;partial=00"`` ->
    ``("c1_1", 101, 109, 1)``.
    """
    parts = [p.strip() for p in description.split("#")]
    if len(parts) != 5:
        raise ValueError(f"unexpected Prodigal header {description!r}")
    try:
        return parts[0], int(parts[1]), int(parts[2]), int(parts[3])
    except ValueError as exc:
        raise ValueError(f"unexpected Prodigal header {description!r}") from exc


def parse_prodigal_records(prodigal_fasta: str) -> dict[str, dict]:
    """The gene file as ``{gene_id: {"record", "start", "end", "strand"}}``.

    This is dN/dS's view of the file: coordinates are left in Prodigal's
    **1-based inclusive** system and the Biopython ``SeqRecord`` is kept
    whole, because dN/dS reads its sequence and copies its description onto
    the ancestral sequences it writes out.  Nothing is validated beyond the
    header shape (dN/dS checks sequence bounds itself, per site), so the
    numbers downstream are exactly what Prodigal wrote.

    Later duplicates of an id replace earlier ones, which is what
    ``SeqIO.to_dict`` followed by a loop used to do; Prodigal ids are unique
    per file so this never matters in practice.

    Example: ``>c1_1 # 101 # 109 # -1 # ID=1_1`` over ``ATGGGATGG`` ->
    ``{"c1_1": {"record": <SeqRecord c1_1>, "start": 101, "end": 109,
    "strand": -1}}``.
    """
    records: dict[str, dict] = {}
    for rec in SeqIO.parse(prodigal_fasta, "fasta"):
        _, start, end, strand = _header_fields(rec.description)
        records[rec.id] = {"record": rec, "start": start, "end": end, "strand": strand}
    return records


def load_genes(prodigal_fasta: str) -> pd.DataFrame:
    """The gene file as one row per gene with **0-based inclusive** coordinates.

    Columns: ``gene_id, contig, start, end, strand, length, sequence``.  Both
    ends are Prodigal's minus one, so they are directly comparable with
    cache / profile positions; ``contig`` is the id up to its last ``_``
    (Prodigal appends ``_<n>`` to the contig name); ``sequence`` is
    upper-cased so it can be compared with reference letters.  A sequence
    whose length disagrees with its header span raises ``ValueError`` naming
    the gene, since every later index into it would be wrong; a header
    without the four ``#`` fields raises too.  An empty file gives an empty
    frame with the same columns.

    Example: ``>c1_1 # 101 # 109 # 1 # ...`` over ``ATGGGATGG`` -> gene
    ``c1_1`` on contig ``c1``, start 100, end 108, strand 1, length 9.
    """
    rows = []
    for rec in SeqIO.parse(prodigal_fasta, "fasta"):
        gene_id, start1, end1, strand = _header_fields(rec.description)
        # 1-based inclusive -> 0-based inclusive: subtract one from both ends.
        start, end = start1 - 1, end1 - 1
        seq = str(rec.seq).upper()
        if len(seq) != end - start + 1:
            raise ValueError(
                f"gene {gene_id}: header spans {end - start + 1} bases but the "
                f"sequence has {len(seq)}"
            )
        rows.append(
            {
                "gene_id": gene_id,
                "contig": "_".join(gene_id.split("_")[:-1]),
                "start": start,
                "end": end,
                "strand": strand,
                "length": len(seq),
                "sequence": seq,
            }
        )
    return pd.DataFrame(rows, columns=GENE_COLUMNS)


def coding_index(position: int, start0: int, end0: int, strand: int) -> int:
    """Index into the coding sequence of the base at contig ``position``.

    ``start0`` / ``end0`` are the gene's **0-based inclusive** ends and
    ``position`` a 0-based contig coordinate.  On strand ``1`` the coding
    sequence starts at ``start0``, so the index is ``position - start0``.  On
    strand ``-1`` the file holds the reverse complement, whose first letter is
    the contig base at ``end0``, so the index is ``end0 - position``.  The
    letter found there is the contig letter (strand 1) or its complement
    (strand -1).

    Raises ``ValueError`` when ``position`` lies outside ``[start0, end0]``
    (the site is not in this gene) or ``strand`` is not ``1`` / ``-1``.

    Example: gene at 100..108.  Strand 1: position 104 -> 4.  Strand -1:
    position 104 -> 108 - 104 = 4 as well, but position 100 -> 8 and
    108 -> 0: the strand flips which end is "first".
    """
    if not start0 <= position <= end0:
        raise ValueError(f"position {position} is outside the gene ({start0}..{end0})")
    if strand == 1:
        return position - start0
    if strand == -1:
        return end0 - position
    raise ValueError(f"unknown strand {strand!r}; expected 1 or -1")
