"""
Contains auxiliary functions for run_dna_rna_tools
"""

complement_dict_for_dna: dict = {
    "A": "T",
    "G": "C",
    "C": "G",
    "T": "A",
    "a": "t",
    "g": "c",
    "c": "g",
    "t": "a",
}

complement_dict_for_rna: dict = {
    "A": "U",
    "G": "C",
    "C": "G",
    "U": "A",
    "a": "u",
    "g": "c",
    "c": "g",
    "u": "c",
}

transcribe_dict: dict = {
    "A": "A",
    "G": "G",
    "C": "C",
    "T": "U",
    "a": "a",
    "g": "g",
    "c": "c",
    "t": "u",
}


def check_na(seq: str) -> bool:
    seq = set(seq)
    if seq <= set("AGCTagct") or seq <= set("AGCUagcu"):
        return True
    return False


def check_rna(seq: str) -> bool:
    seq = set(seq)
    if seq <= set("AGCUagcu"):
        return True
    return False


def reverse(seq: str) -> str:
    return seq[::-1]


def complement(seq: str) -> str:
    if check_rna(seq):
        return "".join(complement_dict_for_rna[ch] for ch in seq)
    return "".join(complement_dict_for_dna[ch] for ch in seq)


def reverse_complement(seq: str) -> str:
    seq = reverse(seq)
    seq = complement(seq)
    return seq


def transcribe(seq: str) -> str:
    if not check_rna(seq):
        return "".join(transcribe_dict[ch] for ch in seq)
    return "the sequence is rna"


def search_start_codon_in_rna(seq: str) -> int | str:
    if check_rna(seq):
        if "aug" in seq.lower():
            start = seq.lower().find("aug")
            return start + 1
        return "no start codon"
    return "dna"


def search_first_stop_codon_in_rna(seq: str) -> int | str:
    idx = []
    if check_rna(seq):
        for cdn in ("uga", "uaa", "uag"):
            idx_cdn = seq.lower().find(cdn)
            if idx_cdn != -1:
                idx.append(idx_cdn)
        if len(idx) != 0:
            return min(idx) + 1
        return "no stop codon"
    return "dna"


fn_map: dict = {
    "reverse": reverse,
    "complement": complement,
    "reverse_complement": reverse_complement,
    "transcribe": transcribe,
    "search_start_codon_in_rna": search_start_codon_in_rna,
    "search_first_stop_codon_in_rna": search_first_stop_codon_in_rna,
}
