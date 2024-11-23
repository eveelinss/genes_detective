"""
Contains auxiliary functions for filter_fastq
"""

import os


def calculate_gc_content(sequence: str) -> float:
    gc_count = sequence.count("G") + sequence.count("C")
    if len(sequence) > 0:
        return (gc_count / len(sequence)) * 100
    return 0.0


def calculate_average_quality(quality: str) -> float:
    if len(quality) > 0:
        return sum(ord(char) - 33 for char in quality) / len(quality)
    return 0.0


def is_bounded(bounds: tuple[float | int], x: float) -> bool:
    return bounds[0] <= x <= bounds[1]


def create_filtered_directory_and_file(output_file: str) -> str:
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    filtered_dir = os.path.join(parent_dir, "filtered")
    if not os.path.exists(filtered_dir):
        os.makedirs(filtered_dir)
    output_file = os.path.join(filtered_dir, output_file)
    return output_file


def write_file(output_file: str, lines: list):
    with open(output_file, "a") as file:
        file.writelines(lines)


def filter_fastq(
    seq: str,
    quality: int | float,
    gc_bounds: tuple[int | float] | float,
    length_bounds: tuple[int | float] | float,
    quality_threshold: int,
):
    gc_content = calculate_gc_content(seq)
    avg_quality = calculate_average_quality(quality)

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    if (
        is_bounded(gc_bounds, gc_content)
        and is_bounded(length_bounds, len(seq))
        and avg_quality >= quality_threshold
    ):
        return True
