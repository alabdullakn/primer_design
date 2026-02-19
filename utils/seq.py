DNA = set("ACGT")

def clean_dna(s: str) -> str:
    return "".join([c for c in s.upper() if c in DNA or c == "|"])

def revcomp(seq: str) -> str:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp.get(b, "N") for b in seq.upper()[::-1])

def gc_pct(seq: str) -> float:
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = sum(1 for c in seq if c in ("G", "C"))
    return 100.0 * gc / len(seq)

