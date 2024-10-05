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


def check_rna_dna(seqs: tuple) -> list:
    flags = []
    for i in seqs:
        if set(i) <= set("AGCTagct"):
            flags.append("dna")
        elif set(i) <= set("AGCUagcu"):
            flags.append("rna")
        else:
            flags.append("not na")
    return flags


def reverse(seqs: str, flags: str, seq_number: int) -> str:
    if flags != "not na":
        return seqs[::-1]
    return f"the sequence №{seq_number} is not a nucleic acid"


def complement(seqs: str, flags: str, seq_number: int) -> str:
    if flags == "rna":
        return "".join(complement_dict_for_rna[ch] for ch in seqs)
    elif flags == "dna":
        return "".join(complement_dict_for_dna[ch] for ch in seqs)
    return f"the sequence №{seq_number} is not a nucleic acid"


def reverse_complement(seqs: str, flags: str, seq_number: int) -> str:
    if flags != "not na":
        seqs = reverse(seqs, flags, seq_number)
        seqs = complement(seqs, flags, seq_number)
        return seqs
    return f"the sequence №{seq_number} is not a nucleic acid"


def transcribe(seqs: str, flags: str, seq_number: int) -> str:
    if flags == "dna":
        return "".join(transcribe_dict[ch] for ch in seqs)
    elif flags == "rna":
        return f"the sequence №{seq_number} is rna"
    return f"the sequence №{seq_number} is not a nucleic acid"


def search_start_codon_in_rna(seqs: str, flags: str, seq_number: int) -> str:
    if flags == "rna":
        if "aug" in seqs.lower():
            start = seqs.lower().find("aug")
            return f"in sequense №{seq_number} the start codon starts with the {start + 1} character"
        return f"in sequense №{seq_number} no start codon"
    elif flags == "dna":
        return f"the sequense №{seq_number} is dna"
    return f"the sequence №{seq_number} is not a nucleic acid"


def search_first_stop_codon_in_rna(seqs: str, flags: str, seq_number: int) -> str:
    idx = []
    if flags == "rna":
        for cdn in ("uga", "uaa", "uag"):
            idx_cdn = seqs.lower().find(cdn)
            if idx_cdn != -1:
                idx.append(idx_cdn)
        if len(idx) != 0:
            return f"in sequense №{seq_number} the first stop codon starts with the {min(idx) + 1} character"
        return f"in sequense №{seq_number} no stop codon"
    elif flags == "dna":
        return f"the sequense №{seq_number} is dna"
    return f"the sequense №{seq_number} is not a nucleic acid"


fn_map: dict = {
    "reverse": reverse,
    "complement": complement,
    "reverse_complement": reverse_complement,
    "transcribe": transcribe,
    "search_start_codon_in_rna": search_start_codon_in_rna,
    "search_first_stop_codon_in_rna": search_first_stop_codon_in_rna,
}
