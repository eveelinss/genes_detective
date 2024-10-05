"""
Contains auxiliary functions for filter_fastq
"""


def calculate_gc_content(sequence: str) -> float:
    gc_count = sequence.count("G") + sequence.count("C")
    if len(sequence) > 0:
       return (gc_count / len(sequence)) * 100
    return 0.0


def calculate_average_quality(quality: str) -> float:
    if len(quality) > 0:
        return sum(ord(char) - 33 for char in quality) / len(quality)
    return 0.0
