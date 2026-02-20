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
