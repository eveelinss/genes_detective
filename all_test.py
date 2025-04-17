import os

import pytest
from Bio.Seq import Seq
from Bio.SeqIO import parse, write
from Bio.SeqRecord import SeqRecord

from genes_detective import (
    create_filtered_directory_and_file,
    process_bounds,
    filter_fastq,
)

TEST_FASTQ = "test_input.fastq"


@pytest.fixture
def setup_teardown():
    """Fixture for setup and teardown before/after tests"""
    records = [
        SeqRecord(
            Seq("ATGC"),
            id="seq1",
            description="",
            letter_annotations={"phred_quality": [40, 40, 40, 40]},
        ),
        SeqRecord(
            Seq("GCCG"),
            id="seq2",
            description="",
            letter_annotations={"phred_quality": [80, 80, 80, 80]},
        ),
        SeqRecord(
            Seq("TTTT"),
            id="seq3",
            description="",
            letter_annotations={"phred_quality": [10, 10, 10, 10]},
        ),
    ]

    with open(TEST_FASTQ, "w") as f:
        write(records, f, "fastq")

    yield

    if os.path.exists(TEST_FASTQ):
        os.remove(TEST_FASTQ)
    if os.path.exists("fastq_filter.log"):
        os.remove("fastq_filter.log")
    if os.path.exists("filtered"):
        for file in os.listdir("filtered"):
            os.remove(os.path.join("filtered", file))
        os.rmdir("filtered")


class TestFileOperations:
    """Tests for file operations"""

    def test_create_filtered_directory(self, setup_teardown):
        """Test creation of filtered directory"""
        output_file = "test_output.fastq"
        output_path = create_filtered_directory_and_file(output_file)

        assert os.path.exists("filtered")
        assert os.path.basename(output_path) == output_file

    def test_output_file_creation(self, setup_teardown):
        """Test output file creation"""
        filter_fastq(
            input_file=TEST_FASTQ,
            output_file="test_output.fastq",
            gc_bounds=(0, 100),
            length_bounds=(0, 100),
            quality_threshold=0,
        )

        output_path = os.path.join("filtered", "test_output.fastq")
        assert os.path.exists(output_path)

        with open(output_path) as f:
            records = list(parse(f, "fastq"))
            assert len(records) == 3


class TestFilteringFunctionality:
    """Tests for filtering functionality"""

    def test_all_sequences_pass(self, setup_teardown):
        """Test that all sequences pass with --gc 0 100 --length 0 100 -q 0"""
        filter_fastq(
            input_file=TEST_FASTQ,
            output_file="all_pass.fastq",
            gc_bounds=(0, 100),
            length_bounds=(0, 100),
            quality_threshold=0,
        )

        input_count = sum(1 for _ in parse(TEST_FASTQ, "fastq"))
        output_count = sum(1 for _ in parse("filtered/all_pass.fastq", "fastq"))
        assert input_count == output_count

    def test_length_filtering(self, setup_teardown):
        """Test sequence length filtering"""
        filter_fastq(
            input_file=TEST_FASTQ,
            output_file="length_filtered.fastq",
            gc_bounds=(0, 100),
            length_bounds=(4, 4),
            quality_threshold=0,
        )

        output_path = os.path.join("filtered", "length_filtered.fastq")
        with open(output_path) as f:
            records = list(parse(f, "fastq"))
            assert len(records) == 3
            assert all(len(rec.seq) == 4 for rec in records)

    def test_quality_filtering(self, setup_teardown):
        """Test quality filtering"""
        filter_fastq(
            input_file=TEST_FASTQ,
            output_file="quality_filtered.fastq",
            gc_bounds=(0, 100),
            length_bounds=(0, 100),
            quality_threshold=15,
        )

        output_path = os.path.join("filtered", "quality_filtered.fastq")
        with open(output_path) as f:
            records = list(parse(f, "fastq"))
            assert len(records) == 2
            assert all(rec.id in ["seq1", "seq2"] for rec in records)


class TestErrorHandling:
    """Tests for error handling"""

    def test_missing_input_file(self):
        """Test handling of missing input file"""
        with pytest.raises(FileNotFoundError):
            filter_fastq(
                input_file="nonexistent.fastq",
                output_file="output.fastq",
                gc_bounds=(0, 100),
                length_bounds=(0, 100),
                quality_threshold=0,
            )

    def test_invalid_gc_bounds(self):
        """Test handling of invalid GC bounds"""
        with pytest.raises(ValueError):
            process_bounds([10, 20, 30], (0, 100), float)


class TestLogging:
    """Tests for logging"""

    def test_error_logging(self, caplog):
        """Test error logging"""
        with pytest.raises(FileNotFoundError):
            filter_fastq(
                input_file="nonexistent.fastq",
                output_file="output.fastq",
                gc_bounds=(0, 100),
                length_bounds=(0, 100),
                quality_threshold=0,
            )

        assert "Input file not found" in caplog.text
