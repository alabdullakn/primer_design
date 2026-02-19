DNA = set("ACGT")

def clean_dna(s: str) -> str:
    # keep A/C/G/T and allow '|' for junction mode
    return "".join([c for c in s.upper() if c in DNA or c == "|"])

def gc_pct(seq: str) -> float:
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = sum(1 for c in seq if c in ("G", "C"))
    return 100.0 * gc / len(seq)

def tm_wallace(seq: str) -> float:
    seq = seq.upper()
    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")
    return 2.0 * (a + t) + 4.0 * (g + c)

def revcomp(seq: str) -> str:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp.get(b, "N") for b in seq.upper()[::-1])

def has_bad_runs(seq: str, max_run: int = 4) -> bool:
    seq = seq.upper()
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
            if run > max_run:
                return True
        else:
            run = 1
    return False
