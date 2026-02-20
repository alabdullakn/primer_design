import streamlit as st
import pandas as pd
import requests

from primer_engine import design_basic_pcr_primers, print_dimer_report_pair
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
    st.title("Regular PCR")
    st.write("Design a standard forward and reverse primer pair from a single template sequence.")

    # persistent template storage across reruns
    if "reg_template_seq" not in st.session_state:
        st.session_state["reg_template_seq"] = ""

    # ============================
    # Primer parameters
    # ============================

    with st.expander("Primer design parameters", expanded=False):
        c1, c2 = st.columns(2)

        with c1:
            min_len = st.number_input("Min primer length", 16, 60, 18, key="reg_min_len")
            max_len = st.number_input("Max primer length", 16, 80, 25, key="reg_max_len")
            tm_target = st.number_input("Target Tm (°C)", 45.0, 75.0, 60.0, key="reg_tm_target")

        with c2:
            tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 25.0, 5.0, key="reg_tm_tol")
            amp_min = st.number_input("Amplicon min (bp)", 50, 5000, 100, key="reg_amp_min")
            amp_max = st.number_input("Amplicon max (bp)", 80, 8000, 500, key="reg_amp_max")

        st.markdown("### Search windows (Option A)")
        st.caption("Smaller values are allowed, but windows will be clamped safely at runtime.")

        # Option A: allow smaller inputs (min 50), then clamp before calling engine
        c3, c4 = st.columns(2)

        with c3:
            start_window_ui = st.slider(
                "Forward search window (bp from 5')",
                min_value=50,
                max_value=20000,
                value=300,
                step=50,
                key="reg_start_window_ui",
            )

        with c4:
            end_window_ui = st.slider(
                "Reverse search window (bp from 3')",
                min_value=50,
                max_value=20000,
                value=300,
                step=50,
                key="reg_end_window_ui",
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
        text = st.text_area(
            "Paste template sequence (A/C/G/T only)",
            height=220,
            key="reg_seq_paste",
        )
        seq = _parse_fasta_text(">x\n" + (text or ""))
        if seq:
            st.session_state["reg_template_seq"] = seq
            st.success(f"Loaded sequence length: {len(seq)}")

    elif mode == "Upload FASTA":
        file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"], key="reg_fasta")
        if file:
            content = file.read().decode("utf-8", errors="ignore")
            seq = _parse_fasta_text(content)
            st.session_state["reg_template_seq"] = seq
            st.success(f"Loaded sequence length: {len(seq)}")

    elif mode == "NCBI accession":
        acc = st.text_input("Enter NCBI accession", key="reg_ncbi")
        if st.button("Fetch from NCBI", key="reg_fetch"):
            try:
                fasta = _fetch_ncbi_fasta(acc)
                seq = _parse_fasta_text(fasta)
                if not seq:
                    st.error("Fetched FASTA but could not parse any A/C/G/T bases.")
                else:
                    st.session_state["reg_template_seq"] = seq
                    st.success(f"Fetched {acc} | length: {len(seq)}")
            except Exception as e:
                st.error(str(e))

    template = st.session_state.get("reg_template_seq", "")

    if not template:
        st.info("Paste a sequence, upload a FASTA file, or fetch an NCBI accession.")
        add_footer()
        return

    st.markdown("---")

    run = st.button("Design primers", key="reg_run")
    if not run:
        add_footer()
        return

    # ============================
    # Option A clamp logic
    # ============================

    n = len(template)

    # clamp windows to safe values
    # must be at least 2x min_len, and cannot exceed template length
    min_window_needed = max(50, int(min_len) * 2)

    start_window = max(min_window_needed, min(int(start_window_ui), n))
    end_window = max(min_window_needed, min(int(end_window_ui), n))

    # also clamp amplicon constraints to sane ordering
    if int(amp_min) > int(amp_max):
        st.error("Amplicon min cannot be greater than amplicon max.")
        add_footer()
        return

    # ============================
    # Run design
    # ============================

    try:
        fwd, rev, amp_len, fwd_start, rev_bind_start = design_basic_pcr_primers(
            template,
            min_len=int(min_len),
            max_len=int(max_len),
            tm_target=float(tm_target),
            tm_tol=float(tm_tol),
            start_window=int(start_window),
            end_window=int(end_window),
            amplicon_min=int(amp_min),
            amplicon_max=int(amp_max),
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
                "Bind start (0-based)": rev_bind_start,
            },
        ]
        st.dataframe(pd.DataFrame(rows), use_container_width=True)

        st.write(f"Estimated amplicon length: **{amp_len} bp**")

        st.subheader("Primer-BLAST link (NCBI)")
        url = primer_blast_url_pair(fwd, rev, "Homo sapiens")
        st.markdown(f"[Open in Primer-BLAST]({url})")
        st.info(BLAST_INSTRUCTIONS)

        st.subheader("Dimer check")
        print_dimer_report_pair(fwd, rev)

        add_footer()

    except Exception as e:
        st.error(str(e))
        st.caption(f"Debug: template length={n}, start_window={start_window}, end_window={end_window}")
        add_footer()
