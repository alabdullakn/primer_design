from utils.tm import tm_nn

DNA = set("ACGT")


def clean_dna(s: str) -> str:
    """Remove non-ACGT characters, allow | (junction marker used in some tools)."""
    return "".join(c for c in s.upper() if c in DNA or c == "|")


def clean_seq(s: str) -> str:
    """Remove non-ACGT characters strictly. Raises if result is empty."""
    s = "".join(c for c in (s or "").upper() if c in DNA)
    if not s:
        raise ValueError("No A/C/G/T bases found. Paste a DNA sequence.")
    return s


def gc_pct(seq: str) -> float:
    if not seq:
        return 0.0
    return 100.0 * sum(c in "GC" for c in seq.upper()) / len(seq)


# Alias used in some engines
gc_content = gc_pct


def revcomp(seq: str) -> str:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp.get(b, "N") for b in seq.upper()[::-1])


def max_run(seq: str) -> int:
    """Return the length of the longest homopolymer run."""
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def has_bad_runs(seq: str, max_run_len: int = 4) -> bool:
    return max_run(seq) > max_run_len


def has_3prime_gc_clamp(seq: str) -> bool:
    if len(seq) < 2:
        return False
    return seq[-1] in "GC" or seq[-2] in "GC"


def has_3prime_complementarity(p1: str, p2: str, k: int = 4) -> bool:
    if len(p1) < k or len(p2) < k:
        return False
    return revcomp(p1[-k:]) in p2


def dimer_risk_percent(p1: str, p2: str, max_k: int = 8) -> float:
    best = 0
    for k in range(3, max_k + 1):
        if has_3prime_complementarity(p1, p2, k=k):
            best = k
    if best == 0:
        return 0.0
    return min(100.0, (best - 2) / (max_k - 2) * 100.0)


def count_stripped(raw: str) -> int:
    """Count non-ACGT characters in a raw sequence string, ignoring FASTA headers and whitespace."""
    total = 0
    for line in raw.splitlines():
        if line.strip().startswith(">"):
            continue
        for c in line.strip().upper():
            if c not in DNA:
                total += 1
    return total


def primer_score(seq: str, tm_target: float) -> float:
    """Simple heuristic score — lower is better."""
    tm = tm_nn(seq)
    score = abs(tm - tm_target)
    if has_bad_runs(seq, max_run_len=4):
        score += 5.0
    gcp = gc_pct(seq)
    if gcp < 35 or gcp > 65:
        score += 3.0
    return score
