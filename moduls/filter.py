"""
Contains auxiliary functions for filter_fastq
"""


def calculate_gc_content(sequence: str) -> float:
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0.0


def calculate_average_quality(quality: str) -> float:
    return (
        sum(ord(char) - 33 for char in quality) / len(quality)
        if len(quality) > 0
        else 0.0
    )
