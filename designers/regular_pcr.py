mport streamlit as st
import pandas as pd
import requests

from primer_engine import design_basic_pcr_primers, print_dimer_report_pair
from utils.blast import primer_blast_url_pair
from utils.primer_utils import gc_pct, tm_wallace, primer_score
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer


# ============================
# FASTA + NCBI helpers
# ============================

def _parse_fasta_or_raw(text: str) -> str:
    if not text:
        return ""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    seq_parts = []
    for line in lines:
        if line.startswith(">"):
            continue
        seq_parts.append(line)
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
        raise RuntimeError("Failed to fetch accession from NCBI (bad accession or NCBI blocked).")
    return r.text


# ============================
# Main
# ============================

def render():
    st.title("Regular PCR")
    st.write("Paste a DNA template, upload a FASTA, or fetch from an NCBI accession, then design one primer pair.")

    # init session state
    if "reg_template" not in st.session_state:
        st.session_state.reg_template = ""

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

        start_window = st.number_input("Search window at 5' end (bp)", 50, 5000, 300, key="reg_start_win")
        end_window = st.number_input("Search window at 3' end (bp)", 50, 5000, 300, key="reg_end_win")
        dimer_k = st.number_input("3' dimer block window (k)", 3, 10, 4, key="reg_dimer_k")

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

    template = ""

    if mode == "Paste sequence":
        raw = st.text_area(
            "Paste template sequence (A/C/G/T only)",
            height=220,
            key="reg_seq_paste",
        )
        template = _parse_fasta_or_raw(raw)

    elif mode == "Upload FASTA":
        file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"], key="reg_fasta")
        if file:
            content = file.read().decode("utf-8", errors="ignore")
            template = _parse_fasta_or_raw(content)
            st.session_state.reg_template = template

        if st.session_state.reg_template:
            st.success(f"Loaded sequence length: {len(st.session_state.reg_template)}")
            template = st.text_area(
                "Template (you can edit before designing)",
                value=st.session_state.reg_template,
                height=180,
                key="reg_template_edit_upload",
            )
            template = _parse_fasta_or_raw(template)

    else:
        acc = st.text_input("Enter NCBI accession", key="reg_ncbi")
        fetch = st.button("Fetch from NCBI", key="reg_fetch")

        if fetch:
            try:
                fasta = _fetch_ncbi_fasta(acc)
                st.session_state.reg_template = _parse_fasta_or_raw(fasta)
                st.success(f"Fetched {acc} | length: {len(st.session_state.reg_template)}")
            except Exception as e:
                st.error(str(e))

        if st.session_state.reg_template:
            template = st.text_area(
                "Template (auto-filled after fetch, you can edit before designing)",
                value=st.session_state.reg_template,
                height=180,
                key="reg_template_edit_ncbi",
            )
            template = _parse_fasta_or_raw(template)

    if not template:
        add_footer()
        return

    st.caption(f"Template length: {len(template)} bp")
    st.markdown("---")

    run = st.button("Design primers", key="reg_run")
    if not run:
        add_footer()
        return

    # ============================
    # Run design
    # ============================

    try:
        fwd_hit, rev_hit, amp_len, _, _ = design_basic_pcr_primers(
            template,
            min_len=int(min_len),
            max_len=int(max_len),
            tm_target=float(tm_target),
            tm_tol=float(tm_tol),
            start_window=int(start_window),
            end_window=int(end_window),
            amplicon_min=int(amp_min),
            amplicon_max=int(amp_max),
            dimer_k=int(dimer_k),
        )

        fwd = fwd_hit.seq_5to3
        rev = rev_hit.seq_5to3

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

        st.subheader("Heterodimer check")
        print_dimer_report_pair(fwd, rev)

        add_footer()

    except Exception as e:
        st.error(str(e))
        add_footer()
