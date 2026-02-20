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



import streamlit as st
import pandas as pd
import requests

from primer_engine import (
    design_qpcr_junction_primers,
    design_qpcr_basic_primers,
    qpcr_amplicon_size,
    print_dimer_report_pair,
)
from utils.blast import primer_blast_url_pair
from utils.primer_utils import gc_pct, tm_wallace, primer_score
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer


def _parse_fasta_text(text: str) -> str:
    if not text:
        return ""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    seq = []
    for l in lines:
        if l.startswith(">"):
            continue
        seq.append(l)
    raw = "".join(seq).upper()
    allowed = set("ACGT")
    return "".join([c for c in raw if c in allowed])


def _fetch_ncbi_fasta(accession: str) -> str:
    acc = (accession or "").strip()
    if not acc:
        return ""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={acc}&rettype=fasta&retmode=text"
    )
    r = requests.get(url, timeout=25)
    if r.status_code != 200:
        raise Exception("Failed to fetch accession from NCBI.")
    return r.text


def render():
    st.title("qPCR")
    st.write("Design primers for SYBR qPCR (short amplicon).")

    if "qpcr_template_seq" not in st.session_state:
        st.session_state["qpcr_template_seq"] = ""
    if "qpcr_junction_seq" not in st.session_state:
        st.session_state["qpcr_junction_seq"] = ""

    # -------------------------
    # Mode
    # -------------------------
    qpcr_mode = st.radio(
        "qPCR mode",
        ["Junction (cDNA, use ^)", "Basic (within exon)"],
        horizontal=True,
        key="qpcr_mode",
    )

    # -------------------------
    # Parameters
    # -------------------------
    with st.expander("Primer design parameters", expanded=False):
        c1, c2 = st.columns(2)

        with c1:
            min_len = st.number_input("Min primer length", 16, 60, 18, key="qp_min_len")
            max_len = st.number_input("Max primer length", 16, 80, 25, key="qp_max_len")
            tm_target = st.number_input("Target Tm (°C)", 45.0, 75.0, 60.0, key="qp_tm_target")

        with c2:
            tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 25.0, 5.0, key="qp_tm_tol")
            amp_min = st.number_input("Amplicon min (bp)", 50, 500, 70, key="qp_amp_min")
            amp_max = st.number_input("Amplicon max (bp)", 60, 800, 200, key="qp_amp_max")

        dimer_k = st.number_input("3' dimer check window (k)", 3, 8, 4, key="qp_dimer_k")

        if qpcr_mode == "Junction (cDNA, use ^)":
            st.caption("Put a ^ exactly at the exon exon junction. Example: ...CAGT^GTTAC...")
            span_primer = st.selectbox(
                "Which primer spans the junction?",
                ["FWD", "REV"],
                index=0,
                key="qp_span_primer",
            )
            overlap = st.number_input(
                "Min overlap on each side of junction (bp)",
                4, 10, 6,
                key="qp_overlap",
            )
        else:
            st.caption("Basic qPCR uses a single template region (no ^).")
            c3, c4 = st.columns(2)
            with c3:
                start_window = st.slider(
                    "Forward search window (bp from 5')",
                    min_value=50, max_value=20000, value=1000, step=50,
                    key="qp_start_window",
                )
            with c4:
                end_window = st.slider(
                    "Reverse search window (bp from 3')",
                    min_value=50, max_value=20000, value=1000, step=50,
                    key="qp_end_window",
                )

    st.markdown("---")

    # -------------------------
    # Input
    # -------------------------
    st.subheader("Sequence input")

    if qpcr_mode == "Junction (cDNA, use ^)":
        seq = st.text_area(
            "Paste sequence with ^ at junction (A/C/G/T plus ^)",
            height=220,
            key="qpcr_junction_text",
        )
        if seq:
            st.session_state["qpcr_junction_seq"] = seq

        if not st.session_state["qpcr_junction_seq"]:
            st.info("Paste a junction sequence containing ^.")
            add_footer()
            return

        run = st.button("Design qPCR primers", key="qpcr_run_junc")
        if not run:
            add_footer()
            return

        try:
            fwd_hit, rev_hit = design_qpcr_junction_primers(
                st.session_state["qpcr_junction_seq"],
                span_primer=st.session_state.get("qp_span_primer", "FWD"),
                min_len=int(min_len),
                max_len=int(max_len),
                tm_target=float(tm_target),
                tm_tol=float(tm_tol),
                amplicon_min=int(amp_min),
                amplicon_max=int(amp_max),
                junction_min_overlap_each_side=int(overlap),
            )

            amp_len = qpcr_amplicon_size(st.session_state["qpcr_junction_seq"], fwd_hit, rev_hit)

            st.success("qPCR junction primers designed successfully.")
            st.caption(SCORE_EXPLANATION)

            rows = [
                {
                    "Type": "FWD",
                    "Primer (5'→3')": fwd_hit.seq_5to3,
                    "Length": fwd_hit.length,
                    "Tm (°C)": round(tm_wallace(fwd_hit.seq_5to3), 1),
                    "GC (%)": round(gc_pct(fwd_hit.seq_5to3), 1),
                    "Score": round(primer_score(fwd_hit.seq_5to3, float(tm_target)), 2),
                },
                {
                    "Type": "REV",
                    "Primer (5'→3')": rev_hit.seq_5to3,
                    "Length": rev_hit.length,
                    "Tm (°C)": round(tm_wallace(rev_hit.seq_5to3), 1),
                    "GC (%)": round(gc_pct(rev_hit.seq_5to3), 1),
                    "Score": round(primer_score(rev_hit.seq_5to3, float(tm_target)), 2),
                },
            ]
            st.dataframe(pd.DataFrame(rows), use_container_width=True)
            st.write(f"Estimated qPCR amplicon length (spliced template): **{amp_len} bp**")

            st.subheader("Primer-BLAST link (NCBI)")
            url = primer_blast_url_pair(fwd_hit.seq_5to3, rev_hit.seq_5to3, "Homo sapiens")
            st.markdown(f"[Open in Primer-BLAST]({url})")
            st.info(BLAST_INSTRUCTIONS)

            st.subheader("Dimer check")
            print_dimer_report_pair(fwd_hit.seq_5to3, rev_hit.seq_5to3)

            add_footer()

        except Exception as e:
            st.error(str(e))
            add_footer()

    else:
        # Basic mode supports paste, FASTA, or NCBI accession like regular PCR
        mode = st.radio(
            "Choose input type",
            ["Paste sequence", "Upload FASTA", "NCBI accession"],
            horizontal=True,
            key="qpcr_input_mode",
        )

        if mode == "Paste sequence":
            text = st.text_area(
                "Paste template sequence (A/C/G/T only)",
                height=220,
                key="qpcr_seq_paste",
            )
            seq = _parse_fasta_text(">x\n" + (text or ""))
            if seq:
                st.session_state["qpcr_template_seq"] = seq
                st.success(f"Loaded sequence length: {len(seq)}")

        elif mode == "Upload FASTA":
            file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"], key="qpcr_fasta")
            if file:
                content = file.read().decode("utf-8", errors="ignore")
                seq = _parse_fasta_text(content)
                st.session_state["qpcr_template_seq"] = seq
                st.success(f"Loaded sequence length: {len(seq)}")

        elif mode == "NCBI accession":
            acc = st.text_input("Enter NCBI accession", key="qpcr_ncbi")
            if st.button("Fetch from NCBI", key="qpcr_fetch"):
                try:
                    fasta = _fetch_ncbi_fasta(acc)
                    seq = _parse_fasta_text(fasta)
                    if not seq:
                        st.error("Fetched FASTA but could not parse any A/C/G/T bases.")
                    else:
                        st.session_state["qpcr_template_seq"] = seq
                        st.success(f"Fetched {acc} | length: {len(seq)}")
                except Exception as e:
                    st.error(str(e))

        template = st.session_state.get("qpcr_template_seq", "")
        if not template:
            st.info("Paste a sequence, upload a FASTA file, or fetch an NCBI accession.")
            add_footer()
            return

        run = st.button("Design qPCR primers", key="qpcr_run_basic")
        if not run:
            add_footer()
            return

        try:
            n = len(template)
            min_window_needed = max(50, int(min_len) * 2)
            start_w = max(min_window_needed, min(int(start_window), n))
            end_w = max(min_window_needed, min(int(end_window), n))

            fwd, rev, amp_len, fwd_start, rev_bind_start = design_qpcr_basic_primers(
                template,
                min_len=int(min_len),
                max_len=int(max_len),
                tm_target=float(tm_target),
                tm_tol=float(tm_tol),
                amplicon_min=int(amp_min),
                amplicon_max=int(amp_max),
                start_window=int(start_w),
                end_window=int(end_w),
                dimer_k=int(dimer_k),
            )

            st.success("qPCR primers designed successfully.")
            st.caption(SCORE_EXPLANATION)

            rows = [
                {
                    "Type": "FWD",
                    "Primer (5'→3')": fwd,
                    "Length": len(fwd),
                    "Tm (°C)": round(tm_wallace(fwd), 1),
                    "GC (%)": round(gc_pct(fwd), 1),
                    "Score": round(primer_score(fwd, float(tm_target)), 2),
                    "Start (0-based)": fwd_start,
                },
                {
                    "Type": "REV",
                    "Primer (5'→3')": rev,
                    "Length": len(rev),
                    "Tm (°C)": round(tm_wallace(rev), 1),
                    "GC (%)": round(gc_pct(rev), 1),
                    "Score": round(primer_score(rev, float(tm_target)), 2),
                    "Bind start (0-based)": rev_bind_start,
                },
            ]
            st.dataframe(pd.DataFrame(rows), use_container_width=True)
            st.write(f"Estimated qPCR amplicon length: **{amp_len} bp**")

            st.subheader("Primer-BLAST link (NCBI)")
            url = primer_blast_url_pair(fwd, rev, "Homo sapiens")
            st.markdown(f"[Open in Primer-BLAST]({url})")
            st.info(BLAST_INSTRUCTIONS)

            st.subheader("Dimer check")
            print_dimer_report_pair(fwd, rev)

            add_footer()

        except Exception as e:
            st.error(str(e))
            add_footer()


