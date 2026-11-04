"""The shared Prodigal gene reader and the codon-index arithmetic.

Two consumers depend on this module agreeing with itself: dN/dS
(``parse_prodigal_records``, 1-based ends) and the gene de novo tables
(``load_genes``, 0-based ends).  ``coding_index`` is the one place that turns a
contig position into an index into the coding sequence, so a strand bug can
only exist once.
"""

import pytest
from Bio.Seq import Seq

from alleleflux.scripts.utilities.genes import (
    COMPLEMENT_MAP,
    coding_index,
    load_genes,
    parse_prodigal_records,
)

# A 9-base gene: ATG GGA TGG = M G W.
CODING = "ATGGGATGG"


def _write_fasta(tmp_path, text):
    fasta = tmp_path / "genes.fna"
    fasta.write_text(text)
    return str(fasta)


class TestComplementMap:
    def test_is_an_involution_over_the_four_bases(self):
        assert set(COMPLEMENT_MAP) == {"A", "C", "G", "T"}
        for base, comp in COMPLEMENT_MAP.items():
            assert COMPLEMENT_MAP[comp] == base


class TestParseProdigalRecords:
    def test_keeps_prodigal_one_based_ends_and_the_record(self, tmp_path):
        fasta = _write_fasta(
            tmp_path,
            ">c1_1 # 101 # 109 # 1 # ID=1_1;partial=00\n"
            + CODING
            + "\n>c1_2 # 201 # 209 # -1 # ID=1_2;partial=00\n"
            + CODING
            + "\n",
        )
        records = parse_prodigal_records(fasta)
        assert list(records) == ["c1_1", "c1_2"]
        fwd, rev = records["c1_1"], records["c1_2"]
        # Prodigal's own numbers, untouched: 1-based inclusive.
        assert (fwd["start"], fwd["end"], fwd["strand"]) == (101, 109, 1)
        assert (rev["start"], rev["end"], rev["strand"]) == (201, 209, -1)
        # The SeqRecord rides along: dN/dS writes its description into the
        # ancestral .ffn files and reads its sequence.
        assert fwd["record"].id == "c1_1"
        assert fwd["record"].description == "c1_1 # 101 # 109 # 1 # ID=1_1;partial=00"
        assert str(fwd["record"].seq) == CODING

    def test_header_parts_are_stripped_of_whitespace(self, tmp_path):
        fasta = _write_fasta(tmp_path, ">g1 #  5 #   13 #  -1  # x\n" + CODING + "\n")
        rec = parse_prodigal_records(fasta)["g1"]
        assert (rec["start"], rec["end"], rec["strand"]) == (5, 13, -1)


class TestLoadGenes:
    def test_rows_are_zero_based_inclusive_with_contig_and_sequence(self, tmp_path):
        fasta = _write_fasta(
            tmp_path,
            ">k1_7_1 # 101 # 109 # 1 # ID=1_1;partial=00\n"
            + CODING.lower()  # lower-case input is upper-cased
            + "\n>k1_7_2 # 201 # 209 # -1 # ID=1_2;partial=00\n"
            + CODING
            + "\n",
        )
        genes = load_genes(fasta)
        assert list(genes.columns) == [
            "gene_id",
            "contig",
            "start",
            "end",
            "strand",
            "length",
            "sequence",
        ]
        assert genes.to_dict("records") == [
            {
                "gene_id": "k1_7_1",
                "contig": "k1_7",  # id up to its LAST underscore
                "start": 100,
                "end": 108,
                "strand": 1,
                "length": 9,
                "sequence": CODING,
            },
            {
                "gene_id": "k1_7_2",
                "contig": "k1_7",
                "start": 200,
                "end": 208,
                "strand": -1,
                "length": 9,
                "sequence": CODING,
            },
        ]

    def test_empty_file_gives_an_empty_frame_with_the_columns(self, tmp_path):
        genes = load_genes(_write_fasta(tmp_path, ""))
        assert genes.empty
        assert "sequence" in genes.columns

    def test_sequence_length_disagreeing_with_the_header_raises_naming_the_gene(
        self, tmp_path
    ):
        fasta = _write_fasta(
            tmp_path, ">c1_1 # 101 # 112 # 1 # ID=1_1\n" + CODING + "\n"
        )
        with pytest.raises(ValueError, match="c1_1"):
            load_genes(fasta)

    def test_malformed_header_raises(self, tmp_path):
        # Only three '#'-separated fields: not a Prodigal header.
        fasta = _write_fasta(tmp_path, ">c1_1 # 101 # 109\n" + CODING + "\n")
        with pytest.raises(ValueError, match="header"):
            load_genes(fasta)


class TestCodingIndex:
    # Gene occupying 0-based contig positions 100..108 (9 bases).
    START0, END0 = 100, 108

    def test_forward_counts_from_the_start(self):
        assert coding_index(100, self.START0, self.END0, 1) == 0
        assert coding_index(104, self.START0, self.END0, 1) == 4
        assert coding_index(108, self.START0, self.END0, 1) == 8

    def test_reverse_counts_from_the_end(self):
        # The file holds the reverse complement, whose first letter sits at
        # the gene's END coordinate on the contig.
        assert coding_index(108, self.START0, self.END0, -1) == 0
        assert coding_index(103, self.START0, self.END0, -1) == 5
        assert coding_index(100, self.START0, self.END0, -1) == 8

    @pytest.mark.parametrize("position", [99, 109])
    @pytest.mark.parametrize("strand", [1, -1])
    def test_outside_the_gene_raises(self, position, strand):
        with pytest.raises(ValueError, match=str(position)):
            coding_index(position, self.START0, self.END0, strand)

    @pytest.mark.parametrize("strand", [0, 2, "+"])
    def test_unknown_strand_raises(self, strand):
        with pytest.raises(ValueError, match="strand"):
            coding_index(104, self.START0, self.END0, strand)


class TestBruteForceAgainstTheContig:
    def test_every_position_reads_the_complemented_contig_letter(self, tmp_path):
        """Build a contig, cut one forward and one reverse gene from it, and
        check that ``coding_index`` + ``COMPLEMENT_MAP`` recover the contig
        letter at every single position of both genes.  This is the property
        both dN/dS and the gene de novo classifier rely on."""
        contig = "TTGACCGATGCATTACGGCATAGCCTAGGTACCATGAGGTTCATCGGAATTCCGAT"
        # Forward gene: contig[5:17] copied as-is.
        fwd_start0, fwd_end0 = 5, 16
        fwd_seq = contig[fwd_start0 : fwd_end0 + 1]
        # Reverse gene: contig[25:40] reverse-complemented (what Prodigal writes).
        rev_start0, rev_end0 = 25, 39
        rev_seq = str(Seq(contig[rev_start0 : rev_end0 + 1]).reverse_complement())
        fasta = _write_fasta(
            tmp_path,
            f">c_1 # {fwd_start0 + 1} # {fwd_end0 + 1} # 1 # ID=1\n{fwd_seq}\n"
            f">c_2 # {rev_start0 + 1} # {rev_end0 + 1} # -1 # ID=2\n{rev_seq}\n",
        )
        genes = load_genes(fasta).set_index("gene_id")
        records = parse_prodigal_records(fasta)
        for gene_id in ("c_1", "c_2"):
            g = genes.loc[gene_id]
            r = records[gene_id]
            # The two readers describe the same span, one base apart.
            assert (r["start"] - 1, r["end"] - 1) == (g["start"], g["end"])
            seen = set()
            for position in range(g["start"], g["end"] + 1):
                idx = coding_index(position, g["start"], g["end"], g["strand"])
                seen.add(idx)
                contig_letter = contig[position]
                expected = (
                    contig_letter if g["strand"] == 1 else COMPLEMENT_MAP[contig_letter]
                )
                assert g["sequence"][idx] == expected, (gene_id, position, idx)
            # Every coding index hit exactly once: the map is a bijection.
            assert seen == set(range(g["length"]))
