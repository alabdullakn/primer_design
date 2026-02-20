# primer_engine.py
# Contains:
# - Exon primer helper (2 forward + 1 reverse)  [your existing]
# - qPCR junction helper                        [your existing]
# - Regular PCR primer pair design              [added]
# - Streamlit dimer reports                     [added]

from dataclasses import dataclass
from typing import List, Tuple, Optional

DNA = set("ACGT")


# ---------------- Data structures ----------------

@dataclass
class BlastSummary:
    ok: bool
    hits_returned: Optional[int] = None
    top_hit_id: Optional[str] = None
    top_hit_def: Optional[str] = None
    top_evalue: Optional[float] = None
    top_identity_pct: Optional[float] = None
    note: Optional[str] = None


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
    template_seq_5to3: Optional[str] = None  # only for REV: binding site in exon (as given)
    blast: Optional[BlastSummary] = None


# ---------------- Basic helpers ----------------

def clean_seq(seq: str) -> str:
    seq = (seq or "").upper()
    seq = "".join(b for b in seq if b in DNA)
    if not seq:
        raise ValueError("No A/C/G/T bases found. Paste a DNA sequence.")
    return seq


def revcomp(seq: str) -> str:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp[b] for b in reversed(seq))


def gc_content(seq: str) -> float:
    return 100.0 * sum(b in "GC" for b in seq) / len(seq)


def tm_wallace(seq: str) -> float:
    return 2.0 * (seq.count("A") + seq.count("T")) + 4.0 * (seq.count("G") + seq.count("C"))


def max_run(seq: str) -> int:
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def has_3prime_gc_clamp(seq: str) -> bool:
    if len(seq) < 2:
        return False
    return (seq[-1] in "GC") or (seq[-2] in "GC")


def self_complementarity_flag(seq: str) -> bool:
    rc = revcomp(seq)
    for i in range(len(seq) - 3):
        if seq[i:i + 4] in rc:
            return True
    return False


def has_3prime_complementarity(p1: str, p2: str, k: int = 4) -> bool:
    if len(p1) < k:
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


# ---------------- Candidate scoring ----------------

def score_candidate(seq: str, tm_target: float, tm_tol: float) -> Optional[Tuple[float, float, float]]:
    tm = tm_wallace(seq)
    gc = gc_content(seq)

    if abs(tm - tm_target) > tm_tol:
        return None
    if not (35.0 <= gc <= 65.0):
        return None
    if max_run(seq) >= 5:
        return None
    if self_complementarity_flag(seq):
        return None

    score = 0.0
    score += abs(tm - tm_target) * 3.0
    if not has_3prime_gc_clamp(seq):
        score += 2.0
    if max_run(seq) == 4:
        score += 2.0

    return score, tm, gc


# ---------------- Exon splicing design (your existing) ----------------

def best_primer_from_exon(
    exon_seq: str,
    exon_name: str,
    kind: str,
    min_len: int,
    max_len: int,
    tm_target: float,
    tm_tol: float
) -> PrimerHit:
    exon_seq = clean_seq(exon_seq)
    best: Optional[PrimerHit] = None

    for L in range(min_len, max_len + 1):
        for i in range(0, len(exon_seq) - L + 1):
            window = exon_seq[i:i + L]
            primer = window if kind == "FWD" else revcomp(window)

            scored = score_candidate(primer, tm_target, tm_tol)
            if scored is None:
                continue

            score, tm, gc = scored
            hit = PrimerHit(kind=kind, exon_name=exon_name, start_0based=i,
                            length=L, seq_5to3=primer, tm_c=tm, gc_pct=gc, score=score)

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
    dimer_k: int = 4
) -> PrimerHit:
    exon_seq = clean_seq(exon_seq_for_reverse)
    candidates: List[PrimerHit] = []

    for L in range(min_len, max_len + 1):
        for i in range(0, len(exon_seq) - L + 1):
            template_seq = exon_seq[i:i + L]
            rev_primer = revcomp(template_seq)

            scored = score_candidate(rev_primer, tm_target, tm_tol)
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
    dimer_k: int = 4
) -> Tuple[PrimerHit, PrimerHit, PrimerHit]:
    fwd1 = best_primer_from_exon(exon1, "Exon1", "FWD", min_len, max_len, tm_target, tm_tol)
    fwd2 = best_primer_from_exon(exon2, "Exon2", "FWD", min_len, max_len, tm_target, tm_tol)
    rev3 = best_reverse_avoiding_fwds(exon3_for_reverse, fwd1, fwd2, min_len, max_len, tm_target, tm_tol, dimer_k=dimer_k)
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


# ---------------- Regular PCR design (ADDED) ----------------

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
    dimer_k: int = 4
) -> Tuple[str, str, int, int, int]:
    """
    Returns:
      (fwd_primer_5to3, rev_primer_5to3, amplicon_len, fwd_start, rev_bind_start)

    rev_bind_start is the 0-based start position of the reverse primer binding site
    on the forward strand (template).
    """

    seq = clean_seq(template_seq)
    n = len(seq)

    # clamp windows to template length
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

    best_pair = None  # (total_score, fwd_seq, rev_seq, amp_len, fwd_start, rev_bind_start)

    # Enumerate forward primer candidates in first start_window bases
    for Lf in range(min_len, max_len + 1):
        for i in range(0, fwd_region_end - Lf + 1):
            fwd = seq[i:i + Lf]
            f_scored = score_candidate(fwd, tm_target, tm_tol)
            if f_scored is None:
                continue
            f_score, f_tm, f_gc = f_scored
            fwd_3 = i + Lf

            # reverse primer binding site must be downstream to make an amplicon
            # reverse binding site start j gives amp_len = (j+Lr) - i
            search_j_min = fwd_3 + amplicon_min - 1
            search_j_max = fwd_3 + amplicon_max

            # constrain to last end_window region
            j_start = max(rev_region_start, search_j_min)
            j_end = min(n, search_j_max)

            if j_start >= j_end:
                continue

            best_rev_for_this_fwd = None  # (score, rev_seq, amp_len, j)
            for Lr in range(min_len, max_len + 1):
                for j in range(j_start, j_end - Lr + 1):
                    bind_site = seq[j:j + Lr]
                    rev = revcomp(bind_site)

                    r_scored = score_candidate(rev, tm_target, tm_tol)
                    if r_scored is None:
                        continue
                    r_score, r_tm, r_gc = r_scored

                    # quick 3' dimer filter between fwd and rev
                    if has_3prime_complementarity(fwd, rev, k=dimer_k):
                        continue

                    amp_len = (j + Lr) - i
                    if amp_len < amplicon_min or amp_len > amplicon_max:
                        continue

                    if best_rev_for_this_fwd is None or r_score < best_rev_for_this_fwd[0]:
                        best_rev_for_this_fwd = (r_score, rev, amp_len, j)

            if best_rev_for_this_fwd is None:
                continue

            r_score, rev, amp_len, j = best_rev_for_this_fwd
            total = f_score + r_score

            if best_pair is None or total < best_pair[0]:
                best_pair = (total, fwd, rev, amp_len, i, j)

    if best_pair is None:
        raise RuntimeError(
            "No primer pair found with your constraints. "
            "Try increasing Tm tolerance, widening length range, or widening amplicon range."
        )

    _, fwd, rev, amp_len, fwd_start, rev_bind_start = best_pair
    return fwd, rev, amp_len, fwd_start, rev_bind_start


def print_dimer_report_pair(fwd_seq: str, rev_seq: str) -> None:
    try:
        import streamlit as st
    except Exception:
        return

    risk = dimer_risk_percent(fwd_seq, rev_seq, max_k=8)
    st.table([{
        "FWD (5'->3')": fwd_seq,
        "REV (5'->3')": rev_seq,
        "Heuristic dimer risk (%)": round(risk, 1)
    }])


# ---------------- qPCR (your existing block, kept) ----------------

def parse_junction_marked_seq(seq_with_marker: str, marker: str = "^") -> Tuple[str, str]:
    s = (seq_with_marker or "")
    if marker not in s:
        raise ValueError("For qPCR, mark the exon-exon junction using '^' inside the sequence.")
    left_raw, right_raw = s.split(marker, 1)
    left = clean_seq(left_raw)
    right = clean_seq(right_raw)
    if len(left) < 10 or len(right) < 10:
        raise ValueError("Junction sides are too short. Paste more bases on each side of '^'.")
    return left, right


def build_junction_primer_forward(left: str, right: str, left_take: int, right_take: int) -> str:
    return left[-left_take:] + right[:right_take]


def build_junction_primer_reverse(left: str, right: str, left_take: int, right_take: int) -> str:
    template = left[-left_take:] + right[:right_take]
    return revcomp(template)


def design_qpcr_junction_primers(
    seq_with_junction_marker: str,
    span_primer: str = "FWD",
    min_len: int = 18,
    max_len: int = 25,
    tm_target: float = 60.0,
    tm_tol: float = 5.0,
    amplicon_min: int = 70,
    amplicon_max: int = 200,
    junction_min_overlap_each_side: int = 6
) -> Tuple[PrimerHit, PrimerHit]:

    left, right = parse_junction_marked_seq(seq_with_junction_marker, marker="^")
    template = left + right
    junction_index = len(left)

    best_pair = None

    for L in range(min_len, max_len + 1):
        for left_take in range(junction_min_overlap_each_side, L - junction_min_overlap_each_side + 1):
            right_take = L - left_take
            if right_take < junction_min_overlap_each_side:
                continue

            if span_primer.upper() == "FWD":
                jseq = build_junction_primer_forward(left, right, left_take, right_take)
                j_scored = score_candidate(jseq, tm_target, tm_tol)
                if j_scored is None:
                    continue
                j_score, j_tm, j_gc = j_scored
                fwd_hit = PrimerHit(
                    kind="FWD",
                    exon_name="Junction",
                    start_0based=junction_index - left_take,
                    length=L,
                    seq_5to3=jseq,
                    tm_c=j_tm,
                    gc_pct=j_gc,
                    score=j_score
                )

                fwd_3prime_pos = fwd_hit.start_0based + fwd_hit.length
                search_start = fwd_3prime_pos + amplicon_min - 1
                search_end = min(len(template), fwd_3prime_pos + amplicon_max)

                if search_start >= search_end:
                    continue

                best_rev = None
                for Lr in range(min_len, max_len + 1):
                    for i in range(search_start, search_end - Lr + 1):
                        bind_site = template[i:i + Lr]
                        rev_seq = revcomp(bind_site)
                        r_scored = score_candidate(rev_seq, tm_target, tm_tol)
                        if r_scored is None:
                            continue
                        r_score, r_tm, r_gc = r_scored
                        rev_hit = PrimerHit(
                            kind="REV",
                            exon_name="RightSide",
                            start_0based=i,
                            length=Lr,
                            seq_5to3=rev_seq,
                            tm_c=r_tm,
                            gc_pct=r_gc,
                            score=r_score,
                            template_seq_5to3=bind_site
                        )
                        if best_rev is None or rev_hit.score < best_rev.score:
                            best_rev = rev_hit

                if best_rev is None:
                    continue

                total = fwd_hit.score + best_rev.score
                if best_pair is None or total < best_pair[0]:
                    best_pair = (total, fwd_hit, best_rev)

            else:
                jseq = build_junction_primer_reverse(left, right, left_take, right_take)
                j_scored = score_candidate(jseq, tm_target, tm_tol)
                if j_scored is None:
                    continue
                j_score, j_tm, j_gc = j_scored
                rev_hit = PrimerHit(
                    kind="REV",
                    exon_name="Junction",
                    start_0based=junction_index - left_take,
                    length=L,
                    seq_5to3=jseq,
                    tm_c=j_tm,
                    gc_pct=j_gc,
                    score=j_score,
                    template_seq_5to3=(left[-left_take:] + right[:right_take])
                )

                rev_bind_start = rev_hit.start_0based
                search_end = max(0, rev_bind_start - amplicon_min + 1)
                search_start = max(0, rev_bind_start - amplicon_max)

                if search_start >= search_end:
                    continue

                best_fwd = None
                for Lf in range(min_len, max_len + 1):
                    for i in range(search_start, search_end - Lf + 1):
                        fseq = template[i:i + Lf]
                        f_scored = score_candidate(fseq, tm_target, tm_tol)
                        if f_scored is None:
                            continue
                        f_score, f_tm, f_gc = f_scored
                        fwd_hit = PrimerHit(
                            kind="FWD",
                            exon_name="LeftSide",
                            start_0based=i,
                            length=Lf,
                            seq_5to3=fseq,
                            tm_c=f_tm,
                            gc_pct=f_gc,
                            score=f_score
                        )
                        if best_fwd is None or fwd_hit.score < best_fwd.score:
                            best_fwd = fwd_hit

                if best_fwd is None:
                    continue

                total = best_fwd.score + rev_hit.score
                if best_pair is None or total < best_pair[0]:
                    best_pair = (total, best_fwd, rev_hit)

    if best_pair is None:
        raise RuntimeError("No qPCR junction primer pair found with your constraints. Try wider Tm tolerance or length range.")

    _, best_fwd, best_rev = best_pair
    return best_fwd, best_rev


def qpcr_amplicon_size(template_seq_with_marker: str, fwd: PrimerHit, rev: PrimerHit) -> int:
    left, right = parse_junction_marked_seq(template_seq_with_marker, marker="^")
    template = left + right
    rev_end = rev.start_0based + rev.length
    return max(0, rev_end - fwd.start_0based)
