from .seq import gc_pct, has_bad_runs
from .tm import tm_nn

def primer_score(seq: str, tm_target: float) -> float:
    """
    Lower is better.
    Simple heuristic: closeness to target Tm,
    penalties for long homopolymers and extreme GC%.
    """
    tm = tm_nn(seq)
    score = abs(tm - tm_target)

    if has_bad_runs(seq, max_run_len=4):
        score += 5.0

    gcp = gc_pct(seq)
    if gcp < 35 or gcp > 65:
        score += 3.0

    return score
