"""
Hi! The convert_multiline_fasta_to_oneline function connects broken sequences in a fasta format file.
The parse_blast_output function finds the best matches from the file for a certain sequence after the BLAST operation
"""

import modules.filter


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    """
    Function convert_multiline_fasta_to_oneline

    Args: input_fasta: str, output_fasta: str=None

    This function receives a fasta file as input, in which the sequence of nucleic acid or protein is
    broken down and broken down sequentially on different lines.
    The function overwrites the data in the correct format in a new file
    """
    if output_fasta is None:
        output_fasta = f"{input_fasta}_corrected"

    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        seq = []
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    outfile.write("".join(seq) + "\n")
                    seq = []
                outfile.write(line + "\n")
            else:
                seq.append(line)
        if seq:
            outfile.write("".join(seq) + "\n")


def parse_blast_output(input_file: str, output_file: str):
    """
    Function parse_blast_output

    Args: input_file: str, output_file: str

    The parse_blast_output function receives the input path to the file
    containing the result of the BLAST operation.
    The function records the best matches from a file for a specific
    sequence in a new file in alphabetical order.
    """
    with open(input_file, "r") as infile:
        for line in infile:
            if "Sequences producing significant alignments:" in line:
                line = infile.readline()
                number_chr_column = len(line.split("Scientific")[0])
                break

        lines = []
        for line in infile:
            if "Sequences producing significant alignments:" in line:
                line = infile.readline()
                line = infile.readline()
                lines.append(infile.readline()[0:number_chr_column] + "\n")

        sorted_lines = sorted(lines, key=str.lower)
        modules.filter.write_file(output_file, sorted_lines)
