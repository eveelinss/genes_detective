## genes_detective - a tool for working with nucleic acid sequences 🔎

**genes_detective** allows you to make various transformations with DNA and RNA, as well as filter fastq sequences according to specified parameters.

![Alt text](https://a.d-cd.net/DqAAAgEQDeA-960.jpg)

## Content

* [Installation](#installation)
* [Main Features](#main-features)
* [Examples](#examples)
* [Acknowledgements](#acknowledgements)

## Installation

Run the following commands:

```sh
git clone git@github.com:eveelinss/genes_detective.git
cd genes_detective
```

You have gained access to the program `genes_detective`.
Profit!

## Main Features

### Biological Sequence Classes

The program provides an object-oriented approach to working with biological sequences, defining the following classes:

- `BiologicalSequence`: An abstract base class for different types of sequences.
- `NucleicAcidSequence`: A base class for DNA and RNA sequences.
- `DNASequence`: A class representing DNA sequences, supporting transcription to RNA.
- `RNASequence`: A class representing RNA sequences.
- `AminoAcidSequence`: A class representing protein sequences, with motif search functionality.

Each sequence class enforces strict alphabet validation and raises `ValueError` if an invalid character is found.

### Nucleic Acid Operations

- **Complementation**: Get the complementary strand.
- **Reverse and Reverse Complement**: Generate the reverse strand and its complement.
- **Transcription**: Convert DNA to RNA.

### Amino Acid Sequence Operations

- **Motif Search**: Find the first occurrence of a given motif in a protein sequence.

### FastQ Filtering

The `filter_fastq` function processes a FastQ file, filtering sequences based on:

- **GC Content (`gc_bounds`)**: Accepts a tuple `(min, max)`, a single value as max, or an int for strict filtering.
- **Sequence Length (`length_bounds`)**: Accepts a tuple `(min, max)`, a single value as max, or an int.
- **Quality Threshold (`quality_threshold`)**: Filters sequences based on average quality score.

The function saves the filtered sequences into a new directory `filtered/`.

## Examples

### Working with DNA Sequences
```python
dna = DNASequence("ATGCGT")
print(dna.complement())  # Outputs: TACGCA
print(dna.reverse_complement())  # Outputs: ACGCAT
print(dna.transcribe())  # Outputs: AUGCGU
```

### Working with RNA Sequences
```python
rna = RNASequence("AUGCGA")
print(rna.complement())  # Outputs: UACGCU
```

### Working with Protein Sequences
```python
protein = AminoAcidSequence("MKTW")
print(protein.find_motif("KT"))  # Outputs: 1
```

### Filtering FastQ Files
```python
filter_fastq("input.fastq", "output.fastq", gc_bounds=(40, 60), length_bounds=(50, 300), quality_threshold=30)
```

## Acknowledgements

I would like to express my gratitude to my nerve cells for their patience and work and to the following people:
- Nikita Vaulin (https://github.com/nvaulin) for understandable lectures,
- Anton Sidorin (https://github.com/SidorinAnton) for detailed seminars,
- Artyom Ershov (https://github.com/iam28th) for edits that allowed us to significantly improve the code. 
