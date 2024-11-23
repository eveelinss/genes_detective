"""
Hello! The function run_dna_rna_tool allows you to convert RNA and DNA.
The filter_fastq function is designed to help filter fastqc by the specified
parameters: gc_bounds, length_bounds and quality_threshold.
"""

import modules.tools_for_rna_dna
import modules.filter


def run_dna_rna_tools(*args: tuple) -> str | list[str | int]:
    """
    Function run_dna_rna_tools

    Args: *args: tuple
    The function accepts a number of arguments: a sequence of strings with RNA
    and/or DNA (at least one) and a string with the operation.
    The following operations are supported
    - reverse
    - complement
    - reverse_complement
    - transcribe
    - search_start_codon_in_rna
    - search_first_codon_in_rna.

    Returns: str or list
    A string or list of converted sequences is returned.
    """
    if len(args) == 0 or args[-1] not in modules.tools_for_rna_dna.fn_map:
        return "the data is entered incorrectly"

    seqs = args[:-1]
    operation = args[-1]
    result = []

    for seq in seqs:
        if modules.tools_for_rna_dna.check_na(seq):
            fn = modules.tools_for_rna_dna.fn_map[operation]
            result.append(fn(seq))
        else:
            result.append("not na")

    if len(result) == 1:
        return result[0]
    return result


def filter_fastq(
    input_file: str,
    output_file: str,
    gc_bounds: tuple | float | int = (0, 100),
    length_bounds: tuple | int = (0, 2**32),
    quality_threshold: int = 0,
):
    """
    Function filter_fastq

    Args: input_file: str, output_file: str,
    gc_bounds: tuple | float | int,
    length_bounds: tuple | int, quality_threshold: int.

    The function takes the input_file string containing the path
    to the source data file, filters the data and writes it to
    a new file in the filtered folder

    The value is a tuple of two strings: consistency and quality.
    The gc_bounds, length_bounds and quality_threshold parameters
    have default values.
    """
    output_file = modules.filter.create_filtered_directory_and_file(output_file)

    with open(input_file, "r") as file:
        while True:
            lines = [file.readline() for _ in range(4)]
            if not any(lines):
                break

            seq = lines[1].rstrip()
            quality = lines[3].rstrip()

            if modules.filter.filter_fastq(
                seq, quality, gc_bounds, length_bounds, quality_threshold
            ):
                modules.filter.write_file(output_file, lines)
