import argparse
import logging
import os
from abc import ABC, abstractmethod
from typing import Tuple, Union

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqUtils import gc_fraction as GC


class BiologicalSequence(ABC):
    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        if not self.check_alphabet():
            raise ValueError("Invalid alphabet")

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, index: int) -> str:
        return self.sequence[index]

    def __str__(self) -> str:
        return self.sequence

    @abstractmethod
    def check_alphabet(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __getattribute__(self, name: str):
        if name == "check_alphabet" and self.__class__ == NucleicAcidSequence:
            raise NotImplementedError("Method not supported for NucleicAcidSequence")
        return object.__getattribute__(self, name)

    def complement(self):
        complemented_seq = self.sequence.translate(str.maketrans(self.complement_map))
        return self.__class__(complemented_seq)

    def reverse(self):
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        return self.reverse().complement()

    def check_alphabet(self) -> bool:
        return all(nucleotide in self.possible_letters for nucleotide in self.sequence)


class DNASequence(NucleicAcidSequence):
    possible_letters = ("A", "T", "G", "C")
    complement_map = str.maketrans("ATGC", "TACG")

    def transcribe(self):
        return RNASequence(self.sequence.translate(str.maketrans("T", "U")))


class RNASequence(NucleicAcidSequence):
    possible_letters = ("A", "U", "G", "C")
    complement_map = str.maketrans("AUGC", "UACG")


class AminoAcidSequence(BiologicalSequence):
    possible_letters = (
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    )

    def check_alphabet(self) -> bool:
        return all(aa in self.possible_letters for aa in self.sequence)

    def find_motif(self, motif: str) -> int:
        """
        Finds the starting index of a first motif in the sequence.
        Returns:
            int: Index of the motif start, or -1 if not found.
        """
        motif = motif.upper()
        return self.sequence.find(motif)


def setup_logging():
    """Set up logging to file and console"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler("fastq_filter.log"), logging.StreamHandler()],
    )
    return logging.getLogger(__name__)


logger = setup_logging()


def create_filtered_directory_and_file(output_file: str) -> str:
    """
    Creates a 'filtered' directory and returns the path to the output file.

    Args:
        output_file: Name of the output file

    Returns:
        Full path to the output file in the filtered directory

    Raises:
        Exception: If directory creation fails
    """
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        filtered_dir = os.path.join(current_dir, "filtered")

        if not os.path.exists(filtered_dir):
            os.makedirs(filtered_dir)
            logger.info(f"Created results directory: {filtered_dir}")

        output_path = os.path.join(filtered_dir, output_file)
        return output_path

    except Exception as e:
        logger.error(f"Error creating directory: {str(e)}")
        raise


def parse_args() -> argparse.Namespace:
    """Parse command line arguments

    Returns:
        Namespace object containing parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Filter FastQ by GC content, length and quality"
    )

    parser.add_argument("-i", "--input", required=True, help="Input FastQ file")
    parser.add_argument("-o", "--output", required=True, help="Output file")
    parser.add_argument(
        "--gc",
        nargs="+",
        type=float,
        default=[0, 100],
        help="GC content bounds (0-100) as MIN MAX or MAX",
    )
    parser.add_argument(
        "--length",
        nargs="+",
        type=int,
        default=[0, 2**32],
        help="Length bounds as MIN MAX or MAX",
    )
    parser.add_argument(
        "-q", "--quality", type=int, default=0, help="Minimum average quality (0-40)"
    )

    return parser.parse_args()


def process_bounds(values: list, default: Tuple, type_cast) -> Tuple:
    """Process parameter bounds

    Args:
        values: List of input values
        default: Default bounds tuple
        type_cast: Type casting function

    Returns:
        Tuple of processed bounds (min, max)

    Raises:
        ValueError: If input values are invalid
    """
    try:
        if len(values) == 0:
            return default
        elif len(values) == 1:
            return (0, type_cast(values[0]))
        elif len(values) == 2:
            return (type_cast(values[0]), type_cast(values[1]))
        else:
            raise ValueError(f"Invalid number of values: {len(values)}")
    except ValueError as e:
        logger.error(f"Error processing bounds: {str(e)}")
        raise


def filter_fastq(
    input_file: str,
    output_file: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """
    Filter FastQ file by specified parameters

    Args:
        input_file: Path to input FastQ file
        output_file: Name of output file
        gc_bounds: GC content bounds (min, max)
        length_bounds: Length bounds (min, max)
        quality_threshold: Minimum average quality score

    Raises:
        FileNotFoundError: If input file doesn't exist
        Exception: If processing fails
    """
    try:
        output_path = create_filtered_directory_and_file(output_file)
        logger.info(f"Started processing file: {input_file}")

        total = 0
        passed = 0

        with open(input_file) as input_handle, open(output_path, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fastq"):
                total += 1
                gc_content = GC(record.seq)
                qualities = record.letter_annotations.get("phred_quality", [])
                avg_quality = sum(qualities) / len(qualities) if qualities else 0

                if (
                    gc_bounds[0] <= gc_content <= gc_bounds[1]
                    and length_bounds[0] <= len(record.seq) <= length_bounds[1]
                    and avg_quality >= quality_threshold
                ):
                    SeqIO.write(record, output_handle, "fastq")
                    passed += 1

        logger.info(
            f"Processing complete. Total sequences: {total}, "
            f"passed filter: {passed} ({passed/total*100:.1f}%)"
        )
        logger.info(f"Results saved to: {output_path}")

    except FileNotFoundError:
        logger.error(f"Input file not found: {input_file}")
        raise
    except Exception as e:
        logger.error(f"Error processing file: {str(e)}")
        raise


def main():
    """Main function to run the FastQ filtering"""
    try:
        args = parse_args()

        logger.info("Starting FastQ filtering with parameters:")
        logger.info(f"Input file: {args.input}")
        logger.info(f"Output file: {args.output}")
        logger.info(f"GC bounds: {args.gc}")
        logger.info(f"Length bounds: {args.length}")
        logger.info(f"Quality threshold: {args.quality}")

        gc_bounds = process_bounds(args.gc, (0, 100), float)
        length_bounds = process_bounds(args.length, (0, 2**32), int)

        filter_fastq(
            input_file=args.input,
            output_file=args.output,
            gc_bounds=gc_bounds,
            length_bounds=length_bounds,
            quality_threshold=args.quality,
        )

    except Exception as e:
        logger.error(f"Critical error: {str(e)}")
        raise


if __name__ == "__main__":
    main()
