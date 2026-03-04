# primer_engine.py
# Contains:
# - Exon primer helper (2 forward + 1 reverse)
# - Regular PCR primer pair design
# - qPCR junction primer design (legacy: design_qpcr_junction_primers)
# - qPCR SYBR vs TaqMan (design_qpcr_junction_pair)
# - Streamlit dimer reports

from dataclasses import dataclass
from typing import List, Tuple, Optional

from utils.tm import tm_nn, has_hairpin, has_homodimer
from utils.seq import (
    clean_seq, revcomp, gc_content, max_run,
    has_3prime_gc_clamp, has_3prime_complementarity, dimer_risk_percent,
)

DNA = set("ACGT")


# ---------------- Data structures ----------------

@dataclass
class PrimerHit:
    kind: str                      # "FWD" or "REV"
    exon_name: str                 # "Exon1", "Exon2", "Exon3", "Junction", etc
    start_0based: int
    length: int
    seq_5to3: str                  # primer you ORDER (5'->3')
    tm_c: float
    gc_pct: float
    score: float
    template_seq_5to3: Optional[str] = None  # only for REV: binding site on forward strand




# ---------------- Candidate scoring ----------------
def score_candidate(
    seq: str,
    tm_target: float,
    tm_tol: float,
    sodium_mM: float = 50.0,
    mg_mM: float = 1.5,
    dntp_mM: float = 0.2,
    gc_min: float = 35.0,
    gc_max: float = 65.0,
) -> Optional[Tuple[float, float, float]]:
    tm = tm_nn(seq, mv_conc=sodium_mM, dv_conc=mg_mM, dntp_conc=dntp_mM)
    gc = gc_content(seq)

    if abs(tm - tm_target) > tm_tol:
        return None
    if not (gc_min <= gc <= gc_max):
        return None
    if max_run(seq) >= 5:
        return None
    if has_hairpin(seq, mv_conc=sodium_mM, dv_conc=mg_mM, dntp_conc=dntp_mM):
        return None
    if has_homodimer(seq, mv_conc=sodium_mM, dv_conc=mg_mM, dntp_conc=dntp_mM):
        return None

    score = 0.0
    score += abs(tm - tm_target) * 3.0
    if not has_3prime_gc_clamp(seq):
        score += 2.0
    if max_run(seq) == 4:
        score += 2.0
    score += max(0, 20 - len(seq)) * 0.2  # prefer longer primers

    return score, tm, gc


# ---------------- Exon splicing design ----------------

def best_primer_from_exon(
    exon_seq: str,
    exon_name: str,
    kind: str,
    min_len: int,
    max_len: int,
    tm_target: float,
    tm_tol: float,
    sodium_mM: float = 50.0,
    mg_mM: float = 1.5,
    dntp_mM: float = 0.2,
    gc_min: float = 35.0,
    gc_max: float = 65.0,
) -> PrimerHit:
    exon_seq = clean_seq(exon_seq)
    best: Optional[PrimerHit] = None

    for L in range(min_len, max_len + 1):
        for i in range(0, len(exon_seq) - L + 1):
            window = exon_seq[i:i + L]
            primer = window if kind == "FWD" else revcomp(window)

            scored = score_candidate(
                primer,
                tm_target,
                tm_tol,
                sodium_mM=sodium_mM,
                mg_mM=mg_mM,
                dntp_mM=dntp_mM,
                gc_min=gc_min,
                gc_max=gc_max,
            )
            if scored is None:
                continue

            score, tm, gc = scored
            hit = PrimerHit(
                kind=kind,
                exon_name=exon_name,
                start_0based=i,
                length=L,
                seq_5to3=primer,
                tm_c=tm,
                gc_pct=gc,
                score=score,
            )

            if best is None or hit.score < best.score:
                best = hit

    if best is None:
        raise RuntimeError(f"No {kind} primer found for {exon_name} with your constraints.")
    return best


def best_reverse_avoiding_fwds(
    exon_seq_for_reverse: str,
    fwd1: PrimerHit,
    fwd2: PrimerHit,
    min_len: int,
    max_len: int,
    tm_target: float,
    tm_tol: float,
    dimer_k: int = 4,
    sodium_mM: float = 50.0,
    mg_mM: float = 1.5,
    dntp_mM: float = 0.2,
    gc_min: float = 35.0,
    gc_max: float = 65.0,
) -> PrimerHit:
    exon_seq = clean_seq(exon_seq_for_reverse)
    candidates: List[PrimerHit] = []

    for L in range(min_len, max_len + 1):
        for i in range(0, len(exon_seq) - L + 1):
            template_seq = exon_seq[i:i + L]
            rev_primer = revcomp(template_seq)

            scored = score_candidate(
                rev_primer,
                tm_target,
                tm_tol,
                sodium_mM=sodium_mM,
                mg_mM=mg_mM,
                dntp_mM=dntp_mM,
                gc_min=gc_min,
                gc_max=gc_max,
            )
            if scored is None:
                continue

            if has_3prime_complementarity(fwd1.seq_5to3, rev_primer, k=dimer_k):
                continue
            if has_3prime_complementarity(fwd2.seq_5to3, rev_primer, k=dimer_k):
                continue

            score, tm, gc = scored
            candidates.append(
                PrimerHit(
                    kind="REV",
                    exon_name="Exon3",
                    start_0based=i,
                    length=L,
                    seq_5to3=rev_primer,
                    tm_c=tm,
                    gc_pct=gc,
                    score=score,
                    template_seq_5to3=template_seq
                )
            )

    if not candidates:
        raise RuntimeError(
            "No reverse primer found that avoids 3' complementarity to both forward primers. "
            "Try dimer_k=3, widen tm_tol, or widen length range."
        )

    candidates.sort(key=lambda x: x.score)
    return candidates[0]


def design_exon_primers(
    exon1: str,
    exon2: str,
    exon3_for_reverse: str,
    min_len: int = 16,
    max_len: int = 25,
    tm_target: float = 60.0,
    tm_tol: float = 5.0,
    dimer_k: int = 4,
    sodium_mM: float = 50.0,
    mg_mM: float = 1.5,
    dntp_mM: float = 0.2,
    gc_min: float = 35.0,
    gc_max: float = 65.0,
) -> Tuple[PrimerHit, PrimerHit, PrimerHit]:
    fwd1 = best_primer_from_exon(
        exon1, "Exon1", "FWD", min_len, max_len, tm_target, tm_tol,
        sodium_mM=sodium_mM, mg_mM=mg_mM, dntp_mM=dntp_mM, gc_min=gc_min, gc_max=gc_max,
    )
    fwd2 = best_primer_from_exon(
        exon2, "Exon2", "FWD", min_len, max_len, tm_target, tm_tol,
        sodium_mM=sodium_mM, mg_mM=mg_mM, dntp_mM=dntp_mM, gc_min=gc_min, gc_max=gc_max,
    )
    rev3 = best_reverse_avoiding_fwds(
        exon3_for_reverse, fwd1, fwd2, min_len, max_len, tm_target, tm_tol,
        dimer_k=dimer_k, sodium_mM=sodium_mM, mg_mM=mg_mM, dntp_mM=dntp_mM,
        gc_min=gc_min, gc_max=gc_max,
    )
    return fwd1, fwd2, rev3


def print_dimer_report(f1: PrimerHit, f2: PrimerHit, r3: PrimerHit) -> None:
    try:
        import streamlit as st
    except Exception:
        return

    rows = []
    for f in [f1, f2]:
        risk = dimer_risk_percent(f.seq_5to3, r3.seq_5to3, max_k=8)
        rows.append({
            "Forward": f"{f.kind} {f.exon_name}",
            "Reverse": f"{r3.kind} {r3.exon_name}",
            "Heuristic dimer risk (%)": round(risk, 1)
        })

    st.table(rows)


# ---------------- Regular PCR design ----------------

def design_basic_pcr_primers(
    template_seq: str,
    min_len: int = 18,
    max_len: int = 25,
    tm_target: float = 60.0,
    tm_tol: float = 5.0,
    start_window: int = 300,
    end_window: int = 300,
    amplicon_min: int = 100,
    amplicon_max: int = 500,
    dimer_k: int = 4,
    sodium_mM: float = 50.0,
    mg_mM: float = 1.5,
    dntp_mM: float = 0.2,
    gc_min: float = 35.0,
    gc_max: float = 65.0,
    top_n: int = 5,
) -> List[Tuple[str, str, int, int, int]]:
    """
    Returns up to top_n (fwd, rev, amplicon_len, fwd_start, rev_bind_start) tuples,
    sorted best-first by combined score.

    Pre-scores all forward and reverse candidates independently (O(n) tm_nn calls),
    then pairs them with cheap range + dimer checks — avoids O(n²) thermodynamic calls.
    """

    seq = clean_seq(template_seq)
    n = len(seq)

    start_window = max(1, min(start_window, n))
    end_window = max(1, min(end_window, n))

    fwd_region_end = min(n, start_window)
    rev_region_start = max(0, n - end_window)

    if fwd_region_end < min_len:
        raise RuntimeError("Forward search window is too small for your min primer length.")
    if (n - rev_region_start) < min_len:
        raise RuntimeError("Reverse search window is too small for your min primer length.")
    if amplicon_min > amplicon_max:
        raise RuntimeError("Amplicon min cannot be greater than amplicon max.")

    # ── Step 1: score every forward candidate once ──────────────────────────
    # Each entry: (score, seq, fwd_start, fwd_3end)
    fwd_candidates: List[Tuple[float, str, int, int]] = []
    for Lf in range(min_len, max_len + 1):
        for i in range(0, fwd_region_end - Lf + 1):
            fwd = seq[i:i + Lf]
            scored = score_candidate(
                fwd, tm_target, tm_tol,
                sodium_mM=sodium_mM, mg_mM=mg_mM, dntp_mM=dntp_mM,
                gc_min=gc_min, gc_max=gc_max,
            )
            if scored is None:
                continue
            fwd_candidates.append((scored[0], fwd, i, i + Lf))

    # ── Step 2: score every reverse candidate once ───────────────────────────
    # Each entry: (score, rev_seq, bind_start, bind_end)
    rev_candidates: List[Tuple[float, str, int, int]] = []
    for Lr in range(min_len, max_len + 1):
        for j in range(rev_region_start, n - Lr + 1):
            bind_site = seq[j:j + Lr]
            rev = revcomp(bind_site)
            scored = score_candidate(
                rev, tm_target, tm_tol,
                sodium_mM=sodium_mM, mg_mM=mg_mM, dntp_mM=dntp_mM,
                gc_min=gc_min, gc_max=gc_max,
            )
            if scored is None:
                continue
            rev_candidates.append((scored[0], rev, j, j + Lr))

    if not fwd_candidates:
        raise RuntimeError(
            "No forward primer found with your constraints. "
            "Try increasing Tm tolerance, widening length range, or adjusting GC range."
        )
    if not rev_candidates:
        raise RuntimeError(
            "No reverse primer found with your constraints. "
            "Try increasing Tm tolerance, widening length range, or adjusting GC range."
        )

    # ── Step 3: pair candidates — cheap range + dimer check only ─────────────
    all_pairs = []
    for f_score, fwd, fwd_start, fwd_3 in fwd_candidates:
        best_rev_for_this_fwd = None
        for r_score, rev, j, j_end in rev_candidates:
            amp_len = j_end - fwd_start
            if amp_len < amplicon_min or amp_len > amplicon_max:
                continue
            if j < fwd_3:
                continue
            if has_3prime_complementarity(fwd, rev, k=dimer_k):
                continue
            if best_rev_for_this_fwd is None or r_score < best_rev_for_this_fwd[0]:
                best_rev_for_this_fwd = (r_score, rev, amp_len, j)

        if best_rev_for_this_fwd is None:
            continue
        r_score, rev, amp_len, j = best_rev_for_this_fwd
        all_pairs.append((f_score + r_score, fwd, rev, amp_len, fwd_start, j))

    if not all_pairs:
        raise RuntimeError(
            "No primer pair found with your constraints. "
            "Try increasing Tm tolerance, widening amplicon range, or widening the search windows."
        )

    all_pairs.sort(key=lambda x: x[0])
    return [
        (fwd, rev, amp_len, fwd_start, rev_bind_start)
        for _, fwd, rev, amp_len, fwd_start, rev_bind_start in all_pairs[:top_n]
    ]


def print_dimer_report_pair(fwd_seq: str, rev_seq: str) -> None:
    try:
        import streamlit as st
    except Exception:
        return

    risk_fwd_to_rev = dimer_risk_percent(fwd_seq, rev_seq, max_k=8)
    risk_rev_to_fwd = dimer_risk_percent(rev_seq, fwd_seq, max_k=8)

    st.table([{
        "Pair": "FWD vs REV",
        "Heuristic dimer risk (FWD->REV) %": round(risk_fwd_to_rev, 1),
        "Heuristic dimer risk (REV->FWD) %": round(risk_rev_to_fwd, 1),
    }])
