import streamlit as st
import pandas as pd
import requests

import primer_engine  # must contain design_basic_pcr_primers + print_dimer_report_pair

from utils.blast import primer_blast_url_pair
from utils.primer_utils import gc_pct, tm_wallace, primer_score
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer


# ----------------------------
# FASTA + NCBI helpers
# ----------------------------

def _parse_fasta_text(text: str) -> str:
    """Parses FASTA or raw sequence text and returns A/C/G/T/N only (uppercase)."""
    if not text:
        return ""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    seq_parts = []
    for l in lines:
        if l.startswith(">"):
            continue
        seq_parts.append(l)
    raw = "".join(seq_parts).upper()
    allowed = set("ACGTN")
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


# ----------------------------
# Main Render
# ----------------------------

def render():
    st.title("Regular PCR")
    st.write("Design a standard forward and reverse primer pair from a single template sequence.")

    # session state to persist fetched/uploaded sequence across reruns
    if "reg_sequence" not in st.session_state:
        st.session_state["reg_sequence"] = ""

    # ============================
    # Primer parameters
    # ============================
    with st.expander("Primer design parameters", expanded=False):
        c1, c2 = st.columns(2)

        with c1:
            min_len = st.number_input("Min primer length", 16, 40, 18, key="reg_min_len")
            max_len = st.number_input("Max primer length", 16, 60, 25, key="reg_max_len")
            tm_target = st.number_input("Target Tm (°C)", 45.0, 75.0, 60.0, key="reg_tm_target")

        with c2:
            tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 20.0, 5.0, key="reg_tm_tol")
            amp_min = st.number_input("Amplicon min (bp)", 50, 5000, 100, key="reg_amp_min")
            amp_max = st.number_input("Amplicon max (bp)", 80, 8000, 500, key="reg_amp_max")

        dimer_k = st.number_input("3' dimer check window (k)", 3, 10, 4, key="reg_dimer_k")

    st.markdown("---")

    # ============================
    # Sequence input
    # ============================
    st.header("Sequence input")

    mode = st.radio(
        "Choose input type",
        ["Paste sequence", "Upload FASTA", "NCBI accession"],
        horizontal=True,
        key="reg_input_mode",
    )

    if mode == "Paste sequence":
        pasted = st.text_area(
            "Paste template sequence (A/C/G/T/N only)",
            height=180,
            key="reg_seq_paste",
        )
        seq = _parse_fasta_text(">x\n" + pasted)
        st.session_state["reg_sequence"] = seq

    elif mode == "Upload FASTA":
        file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"], key="reg_fasta")
        if file is not None:
            content = file.read().decode("utf-8", errors="ignore")
            seq = _parse_fasta_text(content)
            st.session_state["reg_sequence"] = seq
            if seq:
                st.success(f"Loaded sequence length: {len(seq)}")

    else:  # NCBI accession
        acc = st.text_input("Enter NCBI accession", key="reg_ncbi")
        if st.button("Fetch from NCBI", key="reg_fetch"):
            try:
                fasta = _fetch_ncbi_fasta(acc)
                seq = _parse_fasta_text(fasta)
                st.session_state["reg_sequence"] = seq
                if seq:
                    st.success(f"Fetched {acc} | length: {len(seq)}")
                else:
                    st.error("Fetched data, but no sequence was parsed.")
            except Exception as e:
                st.error(str(e))

    sequence = st.session_state["reg_sequence"]

    if not sequence:
        st.info("Add a template sequence to continue.")
        add_footer()
        return

    st.caption(f"Template length: {len(sequence)} bp")
    st.markdown("---")

    # ============================
    # Run design
    # ============================
    run = st.button("Design primers", key="reg_run")
    if not run:
        add_footer()
        return

    try:
        fwd, rev, amp_len, fwd_i, rev_i = primer_engine.design_basic_pcr_primers(
            sequence,
            min_len=int(min_len),
            max_len=int(max_len),
            tm_target=float(tm_target),
            tm_tol=float(tm_tol),
            start_window=300,
            end_window=300,
            amplicon_min=int(amp_min),
            amplicon_max=int(amp_max),
            dimer_k=int(dimer_k),
        )

        st.success("Primers designed successfully.")
        st.caption(SCORE_EXPLANATION)

        rows = [
            {
                "Type": "FWD",
                "Primer (5'->3')": fwd,
                "Length": len(fwd),
                "Tm (°C)": round(tm_wallace(fwd), 1),
                "GC (%)": round(gc_pct(fwd), 1),
                "Score": round(primer_score(fwd, tm_target), 2),
                "Start (0-based)": int(fwd_i),
            },
            {
                "Type": "REV",
                "Primer (5'->3')": rev,
                "Length": len(rev),
                "Tm (°C)": round(tm_wallace(rev), 1),
                "GC (%)": round(gc_pct(rev), 1),
                "Score": round(primer_score(rev, tm_target), 2),
                "Bind start (0-based)": int(rev_i),
            },
        ]

        st.dataframe(pd.DataFrame(rows), use_container_width=True)
        st.write(f"Estimated amplicon length: **{amp_len} bp**")

        st.subheader("Primer-BLAST link (NCBI)")
        url = primer_blast_url_pair(fwd, rev, "Homo sapiens")
        st.markdown(f"[Open in Primer-BLAST]({url})")
        st.info(BLAST_INSTRUCTIONS)

        st.subheader("Dimer check")
        primer_engine.print_dimer_report_pair(fwd, rev)

        add_footer()

    except Exception as e:
        st.error(str(e))
        add_footer()
