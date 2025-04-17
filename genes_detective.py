import os
from abc import ABC, abstractmethod
from typing import Tuple, Union
from Bio import SeqIO
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
    "A", "C", "D", "E", "F", "G", "H",
    "I", "K", "L", "M", "N", "P", "Q",
    "R", "S", "T", "V", "W", "Y",
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


def create_filtered_directory_and_file(output_file: str) -> str:
    """
    Creates a 'filtered' directory and returns the path to the output file.

    Args:
        output_file (str): Name of the output file.

    Returns:
        str: Path to the output file.
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    filtered_dir = os.path.join(parent_dir, "filtered")
    if not os.path.exists(filtered_dir):
        os.makedirs(filtered_dir)
    return os.path.join(filtered_dir, output_file)


def filter_fastq(
    input_file: str,
    output_file: str,
    gc_bounds: Union[Tuple[float, float], float, int] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """
    Filters a FastQ file based on GC content, length, and average quality.

    Args:
        input_file (str): Path to the input FastQ file.
        output_file (str): Name of the output file.
        gc_bounds (Union[Tuple[float, float], float, int]): Bounds for GC content (default: (0, 100)).
        length_bounds (Union[Tuple[int, int], int]): Bounds for sequence length (default: (0, 2^32)).
        quality_threshold (int): Minimum average quality threshold (default: 0).
    """
    output_file_path = create_filtered_directory_and_file(output_file)
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    with open(output_file_path, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fastq"):
            gc_content = GC(record.seq) * 100
            qualities = record.letter_annotations["phred_quality"]
            avg_quality = sum(qualities) / len(qualities) if qualities else 0
            if (
                gc_bounds[0] <= gc_content <= gc_bounds[1]
                and length_bounds[0] <= len(record.seq) <= length_bounds[1]
                and avg_quality >= quality_threshold
            ):
                SeqIO.write(record, output_handle, "fastq")