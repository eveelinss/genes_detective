"""
Hello! The function run_dna_rna_tool allows you to convert RNA and DNA.
The filter_fastq function is designed to help filter fastqc by the specified
parameters: gc_bounds, length_bounds and quality_threshold.
"""

import moduls.tools_for_rna_dna
import moduls.filter


def run_dna_rna_tools(*args: tuple) -> str | list:
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
    if len(args) == 0 or args[-1] not in moduls.tools_for_rna_dna.fn_map:
        return "the data is entered incorrectly"

    seqs = args[:-1]
    operation = args[-1]
    flags = moduls.tools_for_rna_dna.check_rna_dna(seqs)
    result = []

    for seq_number, (seqs, flags) in enumerate(zip(seqs, flags), start=1):
        fn = moduls.tools_for_rna_dna.fn_map[operation]
        result.append(fn(seqs, flags, seq_number))

    if len(result) == 1:
        return result[0]
    return result


def filter_fastq(
    seqs: dict,
    gc_bounds: tuple | float | int = (0, 100),
    length_bounds: tuple | int = (0, 2**32),
    quality_threshold: int = 0,
) -> dict:
    """
    Function filter_fastq

    Args: seqs: dict, gc_bounds: tuple | float | int,
    length_bounds: tuple | int, quality_threshold: int.
    A dictionary with the following structure is accepted:
    The key is a string, the name of the sequence.
    The value is a tuple of two strings: consistency and quality.
    The gc_bounds, length_bounds and quality_threshold parameters
    have default values.

    Returns: dict
    The function returns a dictionary with the same structure,
    but filtered data.
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    filtered_seqs = {}

    for name, (sequence, quality) in seqs.items():
        gc_content = moduls.filter.calculate_gc_content(sequence)
        if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
            continue
        seq_length = len(sequence)
        if not (length_bounds[0] <= seq_length <= length_bounds[1]):
            continue
        avg_quality = moduls.filter.calculate_average_quality(quality)
        if avg_quality < quality_threshold:
            continue
        filtered_seqs[name] = (sequence, quality)

    return filtered_seqs
