#!/usr/bin/env python3
"""
Comprehensive unit tests for profile_mags.py.

This test suite covers:
- Individual function testing (extract_contigID, map_genes, create_intervalTree, etc.)
- Edge cases and error handling
- Integration-style tests with minimal mock data
- Multiprocessing components

Tests use mocking for pysam objects to avoid dependency on actual BAM/FASTA files.
"""

import os
import tempfile
import unittest
from collections import defaultdict
from io import StringIO
from pathlib import Path
from unittest import mock

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from intervaltree import IntervalTree

# Import functions from profile_mags
from alleleflux.scripts.analysis.profile_mags import (
    create_intervalTree,
    extract_contigID,
    map_genes,
    map_positions_to_genes,
    process_contig,
)


class TestExtractContigID(unittest.TestCase):
    """Tests for extract_contigID function."""

    def test_standard_gene_id(self):
        """Test extraction from standard Prodigal gene ID."""
        gene_id = "SLG1007_DASTool_bins_35.fa_k141_82760_1"
        expected = "SLG1007_DASTool_bins_35.fa_k141_82760"
        self.assertEqual(extract_contigID(gene_id), expected)

    def test_simple_gene_id(self):
        """Test extraction from simple gene ID."""
        gene_id = "contig_1_gene_5"
        expected = "contig_1_gene"
        self.assertEqual(extract_contigID(gene_id), expected)

    def test_single_underscore(self):
        """Test gene ID with single underscore."""
        gene_id = "contig_1"
        expected = "contig"
        self.assertEqual(extract_contigID(gene_id), expected)

    def test_no_underscore(self):
        """Test gene ID without underscore."""
        gene_id = "gene1"
        expected = ""
        self.assertEqual(extract_contigID(gene_id), expected)

    def test_multiple_underscores(self):
        """Test gene ID with many underscores."""
        gene_id = "mag_01_contig_123_scaffold_456_gene_789"
        expected = "mag_01_contig_123_scaffold_456_gene"
        self.assertEqual(extract_contigID(gene_id), expected)


class TestMapGenes(unittest.TestCase):
    """Tests for map_genes function."""

    def setUp(self):
        """Create temporary FASTA file for testing."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary files."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_valid_prodigal_fasta(self):
        """Test parsing of valid Prodigal FASTA file."""
        fasta_content = """>gene_1 # 100 # 300 # 1 # ID=1_1;partial=00;start_type=ATG
ATGCGTACGTAG
>gene_2 # 400 # 600 # -1 # ID=1_2;partial=00;start_type=ATG
ATGAAACCC
"""
        fasta_path = os.path.join(self.temp_dir, "genes.fasta")
        with open(fasta_path, "w") as f:
            f.write(fasta_content)

        genes_df = map_genes(fasta_path)

        # Check DataFrame structure
        self.assertIsInstance(genes_df, pd.DataFrame)
        self.assertEqual(len(genes_df), 2)
        self.assertIn("gene_id", genes_df.columns)
        self.assertIn("start", genes_df.columns)
        self.assertIn("end", genes_df.columns)
        self.assertIn("contig", genes_df.columns)

        # Check first gene (Prodigal uses 1-based, converted to 0-based)
        # Note: gene_id has trailing space due to how header is split by "#"
        gene1 = genes_df[genes_df["gene_id"].str.strip() == "gene_1"].iloc[0]
        self.assertEqual(gene1["start"], 99)  # 100-1
        self.assertEqual(gene1["end"], 299)  # 300-1
        self.assertEqual(gene1["contig"], "gene")

        # Check second gene
        gene2 = genes_df[genes_df["gene_id"].str.strip() == "gene_2"].iloc[0]
        self.assertEqual(gene2["start"], 399)  # 400-1
        self.assertEqual(gene2["end"], 599)  # 600-1

    def test_malformed_header(self):
        """Test that malformed header raises ValueError."""
        fasta_content = """>gene_1 # 100 # 300
ATGCGT
"""
        fasta_path = os.path.join(self.temp_dir, "malformed.fasta")
        with open(fasta_path, "w") as f:
            f.write(fasta_content)

        with self.assertRaises(ValueError) as context:
            map_genes(fasta_path)
        self.assertIn("Unexpected header format", str(context.exception))

    def test_coordinates_conversion(self):
        """Test that 1-based Prodigal coordinates are converted to 0-based."""
        fasta_content = """>test_gene # 1 # 3 # 1 # ID=1_1
ATG
"""
        fasta_path = os.path.join(self.temp_dir, "coords.fasta")
        with open(fasta_path, "w") as f:
            f.write(fasta_content)

        genes_df = map_genes(fasta_path)
        gene = genes_df.iloc[0]
        self.assertEqual(gene["start"], 0)  # 1-1
        self.assertEqual(gene["end"], 2)  # 3-1


class TestCreateIntervalTree(unittest.TestCase):
    """Tests for create_intervalTree function."""

    def test_single_contig_single_gene(self):
        """Test interval tree creation with single gene."""
        genes_df = pd.DataFrame(
            {
                "gene_id": ["gene_1"],
                "contig": ["contig_1"],
                "start": [100],
                "end": [200],
            }
        )

        contig_trees = create_intervalTree(genes_df)

        self.assertIn("contig_1", contig_trees)
        self.assertIsInstance(contig_trees["contig_1"], IntervalTree)
        self.assertEqual(len(contig_trees["contig_1"]), 1)

        # Check interval (IntervalTree uses half-open intervals, so end+1)
        intervals = list(contig_trees["contig_1"])
        self.assertEqual(intervals[0].begin, 100)
        self.assertEqual(intervals[0].end, 201)  # end+1
        self.assertEqual(intervals[0].data, "gene_1")

    def test_single_contig_multiple_genes(self):
        """Test interval tree with multiple genes on same contig."""
        genes_df = pd.DataFrame(
            {
                "gene_id": ["gene_1", "gene_2", "gene_3"],
                "contig": ["contig_1", "contig_1", "contig_1"],
                "start": [100, 300, 500],
                "end": [200, 400, 600],
            }
        )

        contig_trees = create_intervalTree(genes_df)

        self.assertEqual(len(contig_trees["contig_1"]), 3)

        # Test each gene is in the tree
        for gene_id, start, end in zip(
            ["gene_1", "gene_2", "gene_3"], [100, 300, 500], [200, 400, 600]
        ):
            overlaps = contig_trees["contig_1"].at(start)
            self.assertTrue(any(interval.data == gene_id for interval in overlaps))

    def test_multiple_contigs(self):
        """Test interval tree creation with multiple contigs."""
        genes_df = pd.DataFrame(
            {
                "gene_id": ["gene_1", "gene_2", "gene_3"],
                "contig": ["contig_1", "contig_2", "contig_1"],
                "start": [100, 100, 300],
                "end": [200, 200, 400],
            }
        )

        contig_trees = create_intervalTree(genes_df)

        self.assertEqual(len(contig_trees), 2)
        self.assertIn("contig_1", contig_trees)
        self.assertIn("contig_2", contig_trees)
        self.assertEqual(len(contig_trees["contig_1"]), 2)
        self.assertEqual(len(contig_trees["contig_2"]), 1)

    def test_empty_dataframe(self):
        """Test interval tree with empty DataFrame."""
        genes_df = pd.DataFrame(columns=["gene_id", "contig", "start", "end"])

        contig_trees = create_intervalTree(genes_df)

        self.assertIsInstance(contig_trees, defaultdict)
        self.assertEqual(len(contig_trees), 0)

    def test_overlapping_genes(self):
        """Test interval tree with overlapping gene coordinates."""
        genes_df = pd.DataFrame(
            {
                "gene_id": ["gene_1", "gene_2"],
                "contig": ["contig_1", "contig_1"],
                "start": [100, 150],
                "end": [200, 250],
            }
        )

        contig_trees = create_intervalTree(genes_df)

        # Position 175 should overlap both genes
        overlaps = contig_trees["contig_1"].at(175)
        gene_ids = {interval.data for interval in overlaps}
        self.assertEqual(gene_ids, {"gene_1", "gene_2"})


class TestMapPositionsToGenes(unittest.TestCase):
    """Tests for map_positions_to_genes function."""

    def setUp(self):
        """Create test data for position mapping."""
        self.genes_df = pd.DataFrame(
            {
                "gene_id": ["gene_1", "gene_2", "gene_3"],
                "contig": ["contig_1", "contig_1", "contig_2"],
                "start": [100, 300, 100],
                "end": [200, 400, 200],
            }
        )
        self.contig_trees = create_intervalTree(self.genes_df)

    def test_positions_within_genes(self):
        """Test mapping positions that fall within genes."""
        positions_df = pd.DataFrame(
            {
                "contig": ["contig_1", "contig_1", "contig_2"],
                "position": [150, 350, 150],
                "coverage": [10, 20, 30],
            }
        )

        mapped_df = map_positions_to_genes(positions_df, self.contig_trees)

        self.assertEqual(len(mapped_df), 3)
        self.assertIn("gene_id", mapped_df.columns)

        # Check mappings
        gene_ids = mapped_df["gene_id"].tolist()
        self.assertEqual(gene_ids[0], "gene_1")
        self.assertEqual(gene_ids[1], "gene_2")
        self.assertEqual(gene_ids[2], "gene_3")

    def test_positions_outside_genes(self):
        """Test mapping positions that don't overlap any genes."""
        positions_df = pd.DataFrame(
            {
                "contig": ["contig_1", "contig_1"],
                "position": [50, 250],
                "coverage": [10, 20],
            }
        )

        mapped_df = map_positions_to_genes(positions_df, self.contig_trees)

        # Positions should have None for gene_id
        self.assertTrue(mapped_df["gene_id"].isnull().all())

    def test_position_in_multiple_genes(self):
        """Test position that overlaps multiple genes."""
        # Create overlapping genes
        genes_df = pd.DataFrame(
            {
                "gene_id": ["gene_1", "gene_2"],
                "contig": ["contig_1", "contig_1"],
                "start": [100, 150],
                "end": [200, 250],
            }
        )
        trees = create_intervalTree(genes_df)

        positions_df = pd.DataFrame(
            {"contig": ["contig_1"], "position": [175], "coverage": [10]}
        )

        mapped_df = map_positions_to_genes(positions_df, trees)

        # Should contain both gene IDs, comma-separated
        gene_id = mapped_df["gene_id"].iloc[0]
        self.assertIn("gene_1", gene_id)
        self.assertIn("gene_2", gene_id)
        self.assertIn(",", gene_id)

    def test_unknown_contig(self):
        """Test mapping positions from contig not in tree."""
        positions_df = pd.DataFrame(
            {"contig": ["unknown_contig"], "position": [100], "coverage": [10]}
        )

        mapped_df = map_positions_to_genes(positions_df, self.contig_trees)

        # Should have None for gene_id
        self.assertTrue(mapped_df["gene_id"].isnull().all())

    def test_preserves_original_columns(self):
        """Test that original DataFrame columns are preserved."""
        positions_df = pd.DataFrame(
            {
                "contig": ["contig_1"],
                "position": [150],
                "coverage": [10],
                "ref_base": ["A"],
                "A": [5],
                "C": [3],
                "G": [2],
                "T": [0],
            }
        )

        mapped_df = map_positions_to_genes(positions_df, self.contig_trees)

        # All original columns should be present
        for col in positions_df.columns:
            self.assertIn(col, mapped_df.columns)


class TestProcessContig(unittest.TestCase):
    """Tests for process_contig function using mocked pysam objects."""

    def setUp(self):
        """Set up mocks for pysam objects."""
        self.mock_bamfile = mock.MagicMock()
        self.mock_reference_fasta = mock.MagicMock()

    @mock.patch("alleleflux.scripts.analysis.profile_mags.bamfile")
    @mock.patch("alleleflux.scripts.analysis.profile_mags.reference_fasta")
    def test_process_contig_basic(self, mock_ref_fasta, mock_bam):
        """Test basic contig processing with mocked data."""
        # Setup mock pileup column
        mock_pileup_column = mock.MagicMock()
        mock_pileup_column.reference_pos = 100
        mock_pileup_column.reference_name = "contig_1"

        # Setup mock pileup reads
        mock_pileup_read1 = mock.MagicMock()
        mock_pileup_read1.is_del = False
        mock_pileup_read1.is_refskip = False
        mock_pileup_read1.query_position = 0
        mock_pileup_read1.alignment.query_sequence = "ATGC"
        mock_pileup_read1.alignment.mapping_quality = 30

        mock_pileup_read2 = mock.MagicMock()
        mock_pileup_read2.is_del = False
        mock_pileup_read2.is_refskip = False
        mock_pileup_read2.query_position = 0
        mock_pileup_read2.alignment.query_sequence = "ATGC"
        mock_pileup_read2.alignment.mapping_quality = 35

        mock_pileup_column.pileups = [mock_pileup_read1, mock_pileup_read2]

        # Setup mock BAM file
        mock_bam.pileup.return_value = [mock_pileup_column]
        mock_bam.references = ["contig_1"]

        # Setup mock reference FASTA
        mock_ref_fasta.references = ["contig_1"]
        mock_ref_fasta.fetch.return_value = "A"

        result = process_contig("contig_1")

        # Check result structure
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)

        data = result[0]
        self.assertEqual(data["contig"], "contig_1")
        self.assertEqual(data["position"], 100)
        self.assertEqual(data["ref_base"], "A")
        self.assertEqual(data["total_coverage"], 2)
        self.assertEqual(data["A"], 2)
        self.assertEqual(data["C"], 0)
        self.assertEqual(data["G"], 0)
        self.assertEqual(data["T"], 0)
        self.assertIn("30", data["mapq_scores"])
        self.assertIn("35", data["mapq_scores"])

    @mock.patch("alleleflux.scripts.analysis.profile_mags.bamfile")
    @mock.patch("alleleflux.scripts.analysis.profile_mags.reference_fasta")
    def test_process_contig_zero_coverage(self, mock_ref_fasta, mock_bam):
        """Test that positions with zero coverage are skipped."""
        # Setup mock pileup column with deletions only
        mock_pileup_column = mock.MagicMock()
        mock_pileup_column.reference_pos = 100
        mock_pileup_column.reference_name = "contig_1"

        mock_pileup_read = mock.MagicMock()
        mock_pileup_read.is_del = True
        mock_pileup_read.is_refskip = False

        mock_pileup_column.pileups = [mock_pileup_read]

        mock_bam.pileup.return_value = [mock_pileup_column]
        mock_ref_fasta.references = ["contig_1"]
        mock_ref_fasta.fetch.return_value = "A"

        result = process_contig("contig_1")

        # Should return empty list (position skipped)
        self.assertEqual(len(result), 0)

    @mock.patch("alleleflux.scripts.analysis.profile_mags.bamfile")
    @mock.patch("alleleflux.scripts.analysis.profile_mags.reference_fasta")
    def test_process_contig_ambiguous_bases(self, mock_ref_fasta, mock_bam):
        """Test that ambiguous bases are converted to N."""
        mock_pileup_column = mock.MagicMock()
        mock_pileup_column.reference_pos = 100
        mock_pileup_column.reference_name = "contig_1"

        # Create reads with ambiguous bases
        mock_reads = []
        for base in ["R", "Y", "N", "W"]:
            mock_read = mock.MagicMock()
            mock_read.is_del = False
            mock_read.is_refskip = False
            mock_read.query_position = 0
            mock_read.alignment.query_sequence = base
            mock_read.alignment.mapping_quality = 30
            mock_reads.append(mock_read)

        mock_pileup_column.pileups = mock_reads

        mock_bam.pileup.return_value = [mock_pileup_column]
        mock_ref_fasta.references = ["contig_1"]
        mock_ref_fasta.fetch.return_value = "A"

        result = process_contig("contig_1")

        data = result[0]
        # All ambiguous bases should be counted as N
        self.assertEqual(data["N"], 4)
        self.assertEqual(data["A"], 0)
        self.assertEqual(data["C"], 0)
        self.assertEqual(data["G"], 0)
        self.assertEqual(data["T"], 0)

    @mock.patch("alleleflux.scripts.analysis.profile_mags.bamfile")
    @mock.patch("alleleflux.scripts.analysis.profile_mags.reference_fasta")
    def test_process_contig_missing_reference(self, mock_ref_fasta, mock_bam):
        """Test handling of contig not in reference FASTA."""
        mock_pileup_column = mock.MagicMock()
        mock_pileup_column.reference_pos = 100
        mock_pileup_column.reference_name = "contig_1"

        mock_pileup_read = mock.MagicMock()
        mock_pileup_read.is_del = False
        mock_pileup_read.is_refskip = False
        mock_pileup_read.query_position = 0
        mock_pileup_read.alignment.query_sequence = "A"
        mock_pileup_read.alignment.mapping_quality = 30

        mock_pileup_column.pileups = [mock_pileup_read]

        mock_bam.pileup.return_value = [mock_pileup_column]
        mock_ref_fasta.references = []  # Empty - contig not found

        result = process_contig("contig_1")

        # Should use 'X' as reference base
        self.assertEqual(result[0]["ref_base"], "X")

    @mock.patch("alleleflux.scripts.analysis.profile_mags.bamfile")
    @mock.patch("alleleflux.scripts.analysis.profile_mags.reference_fasta")
    def test_process_contig_name_mismatch(self, mock_ref_fasta, mock_bam):
        """Test that contig name mismatch raises ValueError."""
        mock_pileup_column = mock.MagicMock()
        mock_pileup_column.reference_pos = 100
        mock_pileup_column.reference_name = "wrong_contig"

        mock_bam.pileup.return_value = [mock_pileup_column]

        with self.assertRaises(ValueError) as context:
            process_contig("contig_1")
        self.assertIn("Contig name mismatch", str(context.exception))


class TestIntegrationScenarios(unittest.TestCase):
    """Integration-style tests combining multiple functions."""

    def setUp(self):
        """Create test data."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_full_pipeline_simple(self):
        """Test complete pipeline from FASTA to position mapping."""
        # Create Prodigal FASTA
        fasta_content = """>contig_1_gene_1 # 100 # 200 # 1 # ID=1_1
ATGCGT
>contig_1_gene_2 # 300 # 400 # 1 # ID=1_2
ATGAAA
"""
        fasta_path = os.path.join(self.temp_dir, "genes.fasta")
        with open(fasta_path, "w") as f:
            f.write(fasta_content)

        # Parse genes
        genes_df = map_genes(fasta_path)
        self.assertEqual(len(genes_df), 2)

        # Create interval tree
        contig_trees = create_intervalTree(genes_df)
        self.assertIn("contig_1_gene", contig_trees)

        # Create positions DataFrame
        positions_df = pd.DataFrame(
            {
                "contig": ["contig_1_gene", "contig_1_gene", "contig_1_gene"],
                "position": [150, 250, 350],  # 1st in gene1, 2nd outside, 3rd in gene2
                "coverage": [10, 20, 30],
            }
        )

        # Map positions
        mapped_df = map_positions_to_genes(positions_df, contig_trees)

        # Verify mappings (note: gene_id has trailing space from header parsing)
        self.assertEqual(mapped_df.iloc[0]["gene_id"].strip(), "contig_1_gene_1")
        # Second position (250) doesn't have any gene associated with it
        self.assertTrue(pd.isna(mapped_df.iloc[1]["gene_id"]))
        self.assertEqual(mapped_df.iloc[2]["gene_id"].strip(), "contig_1_gene_2")

    def test_multiple_contigs_workflow(self):
        """Test workflow with multiple contigs."""
        fasta_content = """>contig_1_gene_1 # 100 # 200 # 1 # ID=1_1
ATG
>contig_2_gene_1 # 100 # 200 # 1 # ID=2_1
GGG
"""
        fasta_path = os.path.join(self.temp_dir, "multi.fasta")
        with open(fasta_path, "w") as f:
            f.write(fasta_content)

        genes_df = map_genes(fasta_path)
        contig_trees = create_intervalTree(genes_df)

        # Multiple contigs
        self.assertEqual(len(contig_trees), 2)
        self.assertIn("contig_1_gene", contig_trees)
        self.assertIn("contig_2_gene", contig_trees)

        positions_df = pd.DataFrame(
            {
                "contig": ["contig_1_gene", "contig_2_gene"],
                "position": [150, 150],
                "coverage": [10, 20],
            }
        )

        mapped_df = map_positions_to_genes(positions_df, contig_trees)

        # Note: gene_id has trailing space from header parsing
        self.assertEqual(mapped_df.iloc[0]["gene_id"].strip(), "contig_1_gene_1")
        self.assertEqual(mapped_df.iloc[1]["gene_id"].strip(), "contig_2_gene_1")


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def test_extract_contig_empty_string(self):
        """Test extract_contigID with empty string."""
        self.assertEqual(extract_contigID(""), "")

    def test_map_genes_nonexistent_file(self):
        """Test map_genes with nonexistent file."""
        with self.assertRaises(FileNotFoundError):
            map_genes("/nonexistent/file.fasta")

    def test_create_intervaltree_with_invalid_coordinates(self):
        """Test interval tree with start > end.

        Note: IntervalTree raises ValueError for null intervals (where begin >= end).
        This happens at line 273 in create_intervalTree when calling addi().
        """
        genes_df = pd.DataFrame(
            {
                "gene_id": ["gene_1"],
                "contig": ["contig_1"],
                "start": [200],
                "end": [100],  # Invalid: end < start
            }
        )

        # IntervalTree raises ValueError for null intervals
        with self.assertRaises(ValueError) as context:
            contig_trees = create_intervalTree(genes_df)

        self.assertIn("Null Interval", str(context.exception))

    def test_map_positions_boundary_conditions(self):
        """Test position mapping at gene boundaries."""
        genes_df = pd.DataFrame(
            {
                "gene_id": ["gene_1"],
                "contig": ["contig_1"],
                "start": [100],
                "end": [200],
            }
        )
        contig_trees = create_intervalTree(genes_df)

        # Test positions at boundaries
        positions_df = pd.DataFrame(
            {
                "contig": ["contig_1", "contig_1", "contig_1", "contig_1"],
                "position": [99, 100, 200, 201],  # before, start, end, after
            }
        )

        mapped_df = map_positions_to_genes(positions_df, contig_trees)

        # Position 99 should be outside (before start)
        self.assertTrue(pd.isna(mapped_df.iloc[0]["gene_id"]))

        # Position 100 should be inside (at start, inclusive)
        self.assertEqual(mapped_df.iloc[1]["gene_id"], "gene_1")

        # Position 200 should be inside (at end, inclusive)
        self.assertEqual(mapped_df.iloc[2]["gene_id"], "gene_1")

        # Position 201 should be outside (after end, half-open interval)
        self.assertTrue(pd.isna(mapped_df.iloc[3]["gene_id"]))


class TestDataFrameConsistency(unittest.TestCase):
    """Test DataFrame column consistency and types."""

    def test_map_genes_output_columns(self):
        """Test that map_genes produces correct column structure."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">gene_1 # 100 # 200 # 1 # ID=1_1\nATG\n")
            fasta_path = f.name

        try:
            genes_df = map_genes(fasta_path)

            # Check required columns exist
            required_columns = ["gene_id", "start", "end", "contig"]
            for col in required_columns:
                self.assertIn(col, genes_df.columns)

            # Check data types
            self.assertEqual(genes_df["gene_id"].dtype, object)
            self.assertTrue(pd.api.types.is_integer_dtype(genes_df["start"]))
            self.assertTrue(pd.api.types.is_integer_dtype(genes_df["end"]))
            self.assertEqual(genes_df["contig"].dtype, object)
        finally:
            os.unlink(fasta_path)

    def test_map_positions_preserves_types(self):
        """Test that map_positions_to_genes preserves column types."""
        genes_df = pd.DataFrame(
            {"gene_id": ["g1"], "contig": ["c1"], "start": [100], "end": [200]}
        )
        trees = create_intervalTree(genes_df)

        positions_df = pd.DataFrame(
            {"contig": ["c1"], "position": [150], "value": [1.5]}  # float column
        )

        mapped_df = map_positions_to_genes(positions_df, trees)

        # Check that float column is preserved
        self.assertTrue(pd.api.types.is_float_dtype(mapped_df["value"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
