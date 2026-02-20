# designers/regular_pcr.py

import streamlit as st
import pandas as pd
import requests

from utils.blast import primer_blast_url_pair, primer_blast_url_single
from utils.primer_utils import gc_pct, tm_wallace, revcomp, primer_score
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer


# ----------------------------
# FASTA + NCBI helpers
# ----------------------------

def _parse_fasta_text(text: str) -> str:
    """Accept raw sequence or FASTA text and return A/C/G/T only, uppercase."""
    if not text:
        return ""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    seq_parts = []
    for l in lines:
        if l.startswith(">"):
            continue
        seq_parts.append(l)
    raw = "".join(seq_parts).upper()
    return "".join([c for c in raw if c in set("ACGT")])


def _fetch_ncbi_fasta(accession: str) -> str:
    acc = (accession or "").strip()
    if not acc:
        return ""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={acc}&rettype=fasta&retmode=text"
    )
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        raise RuntimeError("Failed to fetch accession from NCBI.")
    return r.text


# ----------------------------
# Simple PCR engine (local)
# ----------------------------

def _max_run(seq: str) -> int:
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def _self_comp_flag(seq: str) -> bool:
    """
    Quick hairpin-ish heuristic: if any 4-mer appears in its own reverse complement.
    """
    rc = revcomp(seq)
    for i in range(0, len(seq) - 3):
        if seq[i : i + 4] in rc:
            return True
    return False


def _has_3p_complement(p1: str, p2: str, k: int) -> bool:
    """
    True if last k of p1 is complementary to any substring of p2.
    """
    if len(p1) < k:
        return False
    return revcomp(p1[-k:]) in p2


def _valid_primer(seq: str, tm_target: float, tm_tol: float) -> bool:
    tm = tm_wallace(seq)
    gc = gc_pct(seq)
    if abs(tm - tm_target) > tm_tol:
        return False
    if not (35.0 <= gc <= 65.0):
        return False
    if _max_run(seq) >= 5:
        return False
    if _self_comp_flag(seq):
        return False
    return True


def design_basic_pcr_primers(
    template: str,
    min_len: int = 18,
    max_len: int = 25,
    tm_target: float = 60.0,
    tm_tol: float = 5.0,
    start_window: int = 2000,
    end_window: int = 2000,
    amplicon_min: int = 100,
    amplicon_max: int = 500,
    dimer_k: int = 4,
    max_candidates: int = 300,
):
    """
    Returns:
      fwd_seq, rev_seq, amplicon_len, fwd_start_0based, rev_bind_start_0based

    Notes:
      - forward candidates scanned from 5' window
      - reverse candidates scanned from 3' window (binding site on template), primer is revcomp(binding_site)
      - picks best pair by combined primer_score with basic dimer screen
    """
    seq = _parse_fasta_text(template)
    if not seq:
        raise RuntimeError("No A/C/G/T bases found in the template.")

    n = len(seq)
    start_window = max(50, min(int(start_window), n))
    end_window = max(50, min(int(end_window), n))

    left_start = 0
    left_end = start_window

    right_start = max(0, n - end_window)
    right_end = n

    # Build forward candidates (5' region)
    fwd_cands = []
    for L in range(min_len, max_len + 1):
        for i in range(left_start, left_end - L + 1):
            p = seq[i : i + L]
            if not _valid_primer(p, tm_target, tm_tol):
                continue
            s = primer_score(p, tm_target)
            fwd_cands.append((s, i, p))

    # Build reverse candidates (3' region) by scanning binding sites on template
    rev_cands = []
    for L in range(min_len, max_len + 1):
        for i in range(right_start, right_end - L + 1):
            bind_site = seq[i : i + L]
            p = revcomp(bind_site)
            if not _valid_primer(p, tm_target, tm_tol):
                continue
            s = primer_score(p, tm_target)
            rev_cands.append((s, i, p))  # i = binding site start on template

    if not fwd_cands or not rev_cands:
        raise RuntimeError(
            "No primer candidates found with your constraints. "
            "Try increasing Tm tolerance, widening length range, or widening search windows."
        )

    fwd_cands.sort(key=lambda x: x[0])
    rev_cands.sort(key=lambda x: x[0])

    fwd_top = fwd_cands[: max_candidates]
    rev_top = rev_cands[: max_candidates]

    best = None  # (pair_score, fwd_seq, rev_seq, amp_len, fwd_start, rev_start)

    for f_score, f_i, f_seq in fwd_top:
        for r_score, r_i, r_seq in rev_top:
            amp_len = (r_i + len(r_seq)) - f_i
            if amp_len < amplicon_min or amp_len > amplicon_max:
                continue

            # Basic 3' dimer screen
            if _has_3p_complement(f_seq, r_seq, k=dimer_k):
                continue
            if _has_3p_complement(r_seq, f_seq, k=dimer_k):
                continue

            pair_score = f_score + r_score
            if best is None or pair_score < best[0]:
                best = (pair_score, f_seq, r_seq, amp_len, f_i, r_i)

    if best is None:
        raise RuntimeError(
            "No primer pair found with your constraints. "
            "Try increasing Tm tolerance, widening length range, widening amplicon range, "
            "or reducing the 3' dimer check window (k)."
        )

    _, fwd, rev, amp_len, fwd_start, rev_start = best
    return fwd, rev, amp_len, fwd_start, rev_start


# ----------------------------
# Streamlit UI
# ----------------------------

def render():
    st.title("Regular PCR")
    st.write("Design a standard forward and reverse primer pair from a single template sequence.")

    # Persist template across reruns (needed for NCBI fetch + Design primers click)
    if "reg_template_seq" not in st.session_state:
        st.session_state.reg_template_seq = ""
    if "reg_template_src" not in st.session_state:
        st.session_state.reg_template_src = ""

    # ============================
    # Primer parameters
    # ============================

    with st.expander("Primer design parameters", expanded=False):
        c1, c2 = st.columns(2)

        with c1:
            min_len = st.number_input("Min primer length", 16, 60, 18, key="reg_min_len")
            max_len = st.number_input("Max primer length", 16, 80, 25, key="reg_max_len")
            dimer_k = st.number_input("3' dimer check window (k)", 3, 10, 4, key="reg_dimer_k")

        with c2:
            tm_target = st.number_input("Target Tm (°C)", 45.0, 75.0, 60.0, key="reg_tm_target")
            tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 25.0, 5.0, key="reg_tm_tol")

        st.markdown("---")

        c3, c4 = st.columns(2)
        with c3:
            amp_min = st.number_input("Amplicon min (bp)", 50, 50000, 100, key="reg_amp_min")
            amp_max = st.number_input("Amplicon max (bp)", 80, 200000, 500, key="reg_amp_max")

        with c4:
            search_mode = st.selectbox(
                "Search mode",
                ["Fast (recommended)", "Thorough"],
                key="reg_search_mode",
            )

            # Sliders feel better than huge number inputs
            win_cap_ui = 20000
            start_window = st.slider(
                "Forward search window (bp from 5')",
                min_value=50,
                max_value=win_cap_ui,
                value=2000,
                step=100,
                key="reg_start_window",
            )
            end_window = st.slider(
                "Reverse search window (bp from 3')",
                min_value=50,
                max_value=win_cap_ui,
                value=2000,
                step=100,
                key="reg_end_window",
            )

    st.markdown("---")

    # ============================
    # Sequence input
    # ============================

    st.subheader("Sequence input")

    mode = st.radio(
        "Choose input type",
        ["Paste sequence", "Upload FASTA", "NCBI accession"],
        horizontal=True,
        key="reg_input_mode",
    )

    if mode == "Paste sequence":
        raw = st.text_area(
            "Paste template sequence (A/C/G/T only)",
            height=200,
            key="reg_seq_paste",
            placeholder="Paste DNA sequence here...",
        )
        seq = _parse_fasta_text(">x\n" + (raw or ""))
        if seq:
            st.session_state.reg_template_seq = seq
            st.session_state.reg_template_src = f"Paste sequence | length: {len(seq)}"

    elif mode == "Upload FASTA":
        file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"], key="reg_fasta")
        if file:
            content = file.read().decode("utf-8", errors="ignore")
            seq = _parse_fasta_text(content)
            if seq:
                st.session_state.reg_template_seq = seq
                st.session_state.reg_template_src = f"Upload FASTA | length: {len(seq)}"
                st.success(f"Loaded sequence length: {len(seq)}")

    else:
        acc = st.text_input("Enter NCBI accession", key="reg_ncbi")
        if st.button("Fetch from NCBI", key="reg_fetch"):
            try:
                fasta = _fetch_ncbi_fasta(acc)
                seq = _parse_fasta_text(fasta)
                if not seq:
                    st.error("Fetched data but could not parse A/C/G/T from it.")
                else:
                    st.session_state.reg_template_seq = seq
                    st.session_state.reg_template_src = f"Fetched {acc} | length: {len(seq)}"
                    st.success(st.session_state.reg_template_src)
            except Exception as e:
                st.error(str(e))

        # No extra template box in NCBI mode (you asked to remove it)
        if st.session_state.reg_template_seq:
            st.success(st.session_state.reg_template_src)
            with st.expander("Preview fetched sequence", expanded=False):
                st.text_area(
                    "Fetched template (read-only preview)",
                    st.session_state.reg_template_seq[:2000] + ("..." if len(st.session_state.reg_template_seq) > 2000 else ""),
                    height=200,
                    disabled=True,
                )

    seq_final = st.session_state.reg_template_seq

    if not seq_final:
        st.info("Paste a sequence, upload a FASTA file, or fetch an NCBI accession.")
        add_footer()
        return

    st.markdown("---")

    run = st.button("Design primers", key="reg_run")
    if not run:
        add_footer()
        return

    # ============================
    # Run design
    # ============================

    try:
        # Clamp windows to actual sequence length
        seq_len = len(seq_final)
        start_window_eff = min(int(start_window), seq_len)
        end_window_eff = min(int(end_window), seq_len)

        max_candidates = 300 if search_mode == "Fast (recommended)" else 1200

        fwd, rev, amp_len, fwd_start, rev_start = design_basic_pcr_primers(
            seq_final,
            min_len=int(min_len),
            max_len=int(max_len),
            tm_target=float(tm_target),
            tm_tol=float(tm_tol),
            start_window=int(start_window_eff),
            end_window=int(end_window_eff),
            amplicon_min=int(amp_min),
            amplicon_max=int(amp_max),
            dimer_k=int(dimer_k),
            max_candidates=int(max_candidates),
        )

        st.success("Primers designed successfully.")
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
                "Bind start (0-based)": rev_start,
            },
        ]

        st.dataframe(pd.DataFrame(rows), use_container_width=True)

        st.write(f"Estimated amplicon length: **{amp_len} bp**")

        st.subheader("Primer-BLAST links (NCBI)")
        org = "Homo sapiens"
        st.markdown(f"[Open in Primer-BLAST]({primer_blast_url_pair(fwd, rev, org)})")
        with st.expander("Single-primer Primer-BLAST links", expanded=False):
            st.markdown(f"**FWD**: [Primer-BLAST]({primer_blast_url_single(fwd, org)})")
            st.markdown(f"**REV**: [Primer-BLAST]({primer_blast_url_single(rev, org)})")
        st.info(BLAST_INSTRUCTIONS)

        add_footer()

    except Exception as e:
        st.error(str(e))
        add_footer()
