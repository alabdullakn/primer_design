import streamlit as st
import pandas as pd
import requests

from utils.blast import primer_blast_url_pair
from utils.primer_utils import gc_pct, tm_wallace, revcomp, primer_score
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer


# ============================
# Helpers: FASTA + cleaning
# ============================

def _parse_fasta_text(text: str) -> str:
    if not text:
        return ""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    seq_parts = []
    for l in lines:
        if l.startswith(">"):
            continue
        seq_parts.append(l)
    raw = "".join(seq_parts).upper()
    allowed = set("ACGT")
    return "".join([c for c in raw if c in allowed])


def _fetch_ncbi_fasta(accession: str) -> str:
    acc = accession.strip()
    if not acc:
        return ""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={acc}&rettype=fasta&retmode=text"
    )
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        raise RuntimeError("Failed to fetch accession from NCBI.")
    return r.text


# ============================
# Primer design (simple)
# ============================

def _max_run(seq: str) -> int:
    if not seq:
        return 0
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
    return best


def _has_3prime_gc_clamp(seq: str) -> bool:
    if len(seq) < 2:
        return False
    return (seq[-1] in "GC") or (seq[-2] in "GC")


def _has_3prime_complementarity(p1: str, p2: str, k: int = 4) -> bool:
    if len(p1) < k or len(p2) < k:
        return False
    return revcomp(p1[-k:]) in p2


def _score_primer(seq: str, tm_target: float, tm_tol: float) -> tuple[bool, float, float, float]:
    """
    Returns: (ok, score, tm, gc)
    Lower score is better.
    """
    tm = tm_wallace(seq)
    gc = gc_pct(seq)

    if abs(tm - tm_target) > tm_tol:
        return (False, 1e9, tm, gc)
    if not (35.0 <= gc <= 65.0):
        return (False, 1e9, tm, gc)
    if _max_run(seq) >= 5:
        return (False, 1e9, tm, gc)

    score = 0.0
    score += abs(tm - tm_target) * 3.0
    if not _has_3prime_gc_clamp(seq):
        score += 2.0
    if _max_run(seq) == 4:
        score += 2.0

    return (True, score, tm, gc)


def design_basic_pcr_pair(
    template: str,
    min_len: int = 18,
    max_len: int = 25,
    tm_target: float = 60.0,
    tm_tol: float = 5.0,
    amplicon_min: int = 100,
    amplicon_max: int = 500,
    start_window: int = 300,
    end_window: int = 300,
    dimer_k: int = 4,
) -> dict:
    """
    Returns dict with fwd, rev, amp_len, fwd_start, rev_bind_start, fwd_stats, rev_stats
    rev is returned as the primer you order (5'->3'), which is revcomp(binding_site).
    """
    seq = template.upper()
    seq = "".join([c for c in seq if c in "ACGT"])
    if len(seq) < max_len + amplicon_min:
        raise RuntimeError("Template too short for your primer/amplicon settings.")

    start_window = min(start_window, len(seq))
    end_window = min(end_window, len(seq))

    best = None  # (total_score, payload)

    # candidate forward sites in first start_window
    for Lf in range(min_len, max_len + 1):
        for fwd_start in range(0, start_window - Lf + 1):
            fwd = seq[fwd_start : fwd_start + Lf]
            ok_f, s_f, tm_f, gc_f = _score_primer(fwd, tm_target, tm_tol)
            if not ok_f:
                continue

            fwd_3prime = fwd_start + Lf

            # reverse binding site must start such that amplicon within [min, max]
            rev_bind_min = fwd_3prime + amplicon_min - 1
            rev_bind_max = fwd_3prime + amplicon_max - 1
            if rev_bind_min >= len(seq):
                continue
            rev_bind_max = min(rev_bind_max, len(seq) - min_len)

            # restrict reverse search to last end_window
            end_region_start = max(0, len(seq) - end_window)
            rev_search_start = max(rev_bind_min, end_region_start)
            rev_search_end = max(rev_search_start, rev_bind_max)

            for Lr in range(min_len, max_len + 1):
                last_start = rev_search_end
                if last_start + Lr > len(seq):
                    last_start = len(seq) - Lr
                for rev_bind_start in range(rev_search_start, last_start + 1):
                    bind_site = seq[rev_bind_start : rev_bind_start + Lr]
                    rev = revcomp(bind_site)

                    ok_r, s_r, tm_r, gc_r = _score_primer(rev, tm_target, tm_tol)
                    if not ok_r:
                        continue

                    # quick dimer screen
                    if _has_3prime_complementarity(fwd, rev, k=dimer_k):
                        continue

                    amp_len = (rev_bind_start + Lr) - fwd_start
                    if not (amplicon_min <= amp_len <= amplicon_max):
                        continue

                    total = s_f + s_r + abs(tm_f - tm_r)
                    payload = {
                        "fwd": fwd,
                        "rev": rev,
                        "amp_len": amp_len,
                        "fwd_start": fwd_start,
                        "rev_bind_start": rev_bind_start,
                        "fwd_tm": tm_f,
                        "fwd_gc": gc_f,
                        "fwd_score": s_f,
                        "rev_tm": tm_r,
                        "rev_gc": gc_r,
                        "rev_score": s_r,
                    }

                    if best is None or total < best[0]:
                        best = (total, payload)

    if best is None:
        raise RuntimeError(
            "No primer pair found with your constraints. "
            "Try increasing Tm tolerance, widening length range, or widening amplicon range."
        )

    return best[1]


def _dimer_risk_percent(p1: str, p2: str, max_k: int = 8) -> float:
    best_k = 0
    for k in range(3, max_k + 1):
        if _has_3prime_complementarity(p1, p2, k=k):
            best_k = k
    if best_k == 0:
        return 0.0
    return min(100.0, (best_k - 2) / (max_k - 2) * 100.0)


# ============================
# Streamlit page
# ============================

def render():
    st.title("Regular PCR")
    st.write("Design a standard forward and reverse primer pair from a single template sequence.")

    # keep template across reruns
    if "reg_template" not in st.session_state:
        st.session_state["reg_template"] = ""

    # ----------------------------
    # Parameters (collapsible)
    # ----------------------------
    with st.expander("Primer design parameters", expanded=False):
        c1, c2 = st.columns(2)

        with c1:
            min_len = st.number_input("Min primer length", 16, 40, 18, key="reg_min_len")
            max_len = st.number_input("Max primer length", 16, 60, 25, key="reg_max_len")
            dimer_k = st.number_input("3' dimer check window (k)", 3, 10, 4, key="reg_dimer_k")

        with c2:
            tm_target = st.number_input("Target Tm (°C)", 45.0, 75.0, 60.0, key="reg_tm_target")
            tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 20.0, 5.0, key="reg_tm_tol")

        c3, c4 = st.columns(2)
        with c3:
            amp_min = st.number_input("Amplicon min (bp)", 50, 5000, 100, key="reg_amp_min")
            start_window = st.number_input("Forward search window (bp from 5')", 50, 5000, 300, key="reg_start_win")
        with c4:
            amp_max = st.number_input("Amplicon max (bp)", 80, 8000, 500, key="reg_amp_max")
            end_window = st.number_input("Reverse search window (bp from 3')", 50, 5000, 300, key="reg_end_win")

    st.markdown("---")

    # ----------------------------
    # Input mode
    # ----------------------------
    st.subheader("Sequence input")

    mode = st.radio(
        "Choose input type",
        ["Paste sequence", "Upload FASTA", "NCBI accession"],
        horizontal=True,
        key="reg_input_mode",
    )

    if mode == "Paste sequence":
        text = st.text_area(
            "Paste template sequence (A/C/G/T only)",
            height=200,
            key="reg_seq_paste",
        )
        parsed = _parse_fasta_text(">x\n" + text)
        st.session_state["reg_template"] = parsed

        if parsed:
            st.success(f"Loaded sequence length: {len(parsed)}")

    elif mode == "Upload FASTA":
        file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"], key="reg_fasta")
        if file:
            content = file.read().decode("utf-8", errors="ignore")
            parsed = _parse_fasta_text(content)
            st.session_state["reg_template"] = parsed
            if parsed:
                st.success(f"Loaded sequence length: {len(parsed)}")

    else:
        acc = st.text_input("Enter NCBI accession", key="reg_ncbi")
        if st.button("Fetch from NCBI", key="reg_fetch"):
            try:
                fasta = _fetch_ncbi_fasta(acc)
                parsed = _parse_fasta_text(fasta)
                st.session_state["reg_template"] = parsed
                if parsed:
                    st.success(f"Fetched {acc} | length: {len(parsed)}")
                else:
                    st.error("Fetched data, but no A/C/G/T sequence was parsed.")
            except Exception as e:
                st.error(str(e))

        if st.session_state["reg_template"]:
            st.caption(f"Template loaded from NCBI. Length: {len(st.session_state['reg_template'])}")

    st.markdown("---")

    # ----------------------------
    # Run
    # ----------------------------
    template = st.session_state["reg_template"]
    if not template:
        st.info("Paste a sequence, upload a FASTA file, or fetch an NCBI accession.")
        add_footer()
        return

    run = st.button("Design primers", key="reg_run")
    if not run:
        add_footer()
        return

    try:
        result = design_basic_pcr_pair(
            template=template,
            min_len=int(min_len),
            max_len=int(max_len),
            tm_target=float(tm_target),
            tm_tol=float(tm_tol),
            amplicon_min=int(amp_min),
            amplicon_max=int(amp_max),
            start_window=int(start_window),
            end_window=int(end_window),
            dimer_k=int(dimer_k),
        )

        fwd = result["fwd"]
        rev = result["rev"]
        amp_len = result["amp_len"]

        st.success("Primers designed successfully.")
        st.caption(SCORE_EXPLANATION)

        rows = [
            {
                "Type": "FWD",
                "Primer (5'→3')": fwd,
                "Length": len(fwd),
                "Tm (°C)": round(tm_wallace(fwd), 1),
                "GC (%)": round(gc_pct(fwd), 1),
                "Score": round(primer_score(fwd, tm_target), 2),
            },
            {
                "Type": "REV",
                "Primer (5'→3')": rev,
                "Length": len(rev),
                "Tm (°C)": round(tm_wallace(rev), 1),
                "GC (%)": round(gc_pct(rev), 1),
                "Score": round(primer_score(rev, tm_target), 2),
            },
        ]
        st.dataframe(pd.DataFrame(rows), use_container_width=True)

        st.write(f"Estimated amplicon length: **{amp_len} bp**")

        st.subheader("Primer-BLAST link (NCBI)")
        url = primer_blast_url_pair(fwd, rev, "Homo sapiens")
        st.markdown(f"[Open in Primer-BLAST]({url})")
        st.info(BLAST_INSTRUCTIONS)

        st.subheader("Dimer screen")
        risk = _dimer_risk_percent(fwd, rev, max_k=8)
        st.table(
            [
                {
                    "Pair": "FWD vs REV",
                    "3' dimer window (k)": int(dimer_k),
                    "Heuristic dimer risk (%)": round(risk, 1),
                }
            ]
        )

        add_footer()

    except Exception as e:
        st.error(str(e))
        add_footer()
