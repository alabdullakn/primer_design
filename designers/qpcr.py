import streamlit as st
import pandas as pd

from primer_engine import (
    design_qpcr_junction_pair,
    qpcr_amplicon_size_from_hits,
    print_dimer_report_pair,
)
from utils.blast import primer_blast_url_pair
from ui.text import BLAST_INSTRUCTIONS
from ui.footer import add_footer


def render():
    st.title("qPCR primers")
    st.write("Design qPCR primer pairs for cDNA. Supports SYBR Green and TaqMan (hydrolysis probe).")

    # ============================
    # Chemistry choice
    # ============================

    chemistry = st.selectbox(
        "Detection chemistry",
        ["SYBR Green", "TaqMan (hydrolysis probe)"],
        key="qpcr_chemistry",
    )
    chem_key = "SYBR" if chemistry.startswith("SYBR") else "TAQMAN"

    st.markdown("---")

    # ============================
    # Parameters
    # ============================

    with st.expander("qPCR design parameters", expanded=False):
        c1, c2 = st.columns(2)

        with c1:
            min_len = st.number_input("Min primer length", 16, 30, 18, key="qpcr_min_len")
            max_len = st.number_input("Max primer length", 16, 40, 24, key="qpcr_max_len")
            primer_tm_target = st.number_input("Primer target Tm (°C)", 50.0, 70.0, 60.0, key="qpcr_tm_target")
            primer_tm_tol = st.number_input("Primer Tm tolerance (± °C)", 0.5, 10.0, 2.0, key="qpcr_tm_tol")

        with c2:
            amplicon_min = st.number_input("Amplicon min (bp)", 50, 400, 70, key="qpcr_amp_min")
            amplicon_max = st.number_input("Amplicon max (bp)", 60, 600, 200, key="qpcr_amp_max")
            max_tm_diff_pair = st.number_input("Max Tm difference between primers (°C)", 0.5, 5.0, 1.0, key="qpcr_tm_diff")
            dimer_k = st.number_input("3' dimer screen window (k)", 3, 10, 4, key="qpcr_dimer_k")

        st.caption("Typical qPCR: 70–200 bp amplicon, 18–24 nt primers, primer Tm around 60°C, primer Tm difference ≤ 1°C.")

        strict = st.checkbox("Strict mode (recommended)", value=True, key="qpcr_strict")
        if strict:
            primer_gc_min = 40.0
            primer_gc_max = 60.0
            max_homopolymer = 3
        else:
            primer_gc_min = 35.0
            primer_gc_max = 65.0
            max_homopolymer = 4

        junction_min_overlap = st.slider(
            "Minimum bases on each side of junction for junction primer",
            min_value=4,
            max_value=10,
            value=6,
            key="qpcr_overlap",
        )

        junction_primer = st.selectbox(
            "Which primer spans the junction",
            ["AUTO (try both)", "FWD spans junction", "REV spans junction"],
            key="qpcr_junction_choice",
        )
        if junction_primer.startswith("AUTO"):
            junction_key = "AUTO"
        elif junction_primer.startswith("FWD"):
            junction_key = "FWD"
        else:
            junction_key = "REV"

    st.markdown("---")

    # ============================
    # Input
    # ============================

    st.subheader("Template input (cDNA junction marked)")
    st.write("Paste sequence with '^' marking the exon-exon junction. Example: `...ACCTG^GTTCA...`")
    seq_with_marker = st.text_area(
        "Junction-marked sequence",
        height=200,
        key="qpcr_seq_marker",
        placeholder="Paste here with ^ at the junction",
    )

    st.markdown("---")

    run = st.button("Design qPCR primers", key="qpcr_run")
    if not run:
        add_footer()
        return

    try:
        fwd, rev, probe = design_qpcr_junction_pair(
            seq_with_junction_marker=seq_with_marker,
            chemistry=chem_key,
            junction_primer=junction_key,
            min_len=int(min_len),
            max_len=int(max_len),
            primer_tm_target=float(primer_tm_target),
            primer_tm_tol=float(primer_tm_tol),
            primer_gc_min=float(primer_gc_min),
            primer_gc_max=float(primer_gc_max),
            max_homopolymer=int(max_homopolymer),
            amplicon_min=int(amplicon_min),
            amplicon_max=int(amplicon_max),
            junction_min_overlap_each_side=int(junction_min_overlap),
            max_tm_diff_pair=float(max_tm_diff_pair),
            dimer_k=int(dimer_k),
        )

        amp_len = qpcr_amplicon_size_from_hits(fwd, rev)

        st.success("qPCR primers designed.")
        st.write(f"Amplicon size (estimated on spliced template): **{amp_len} bp**")

        rows = [
            {
                "Type": "FWD",
                "Role": fwd.exon_name,
                "Primer (5'→3')": fwd.seq_5to3,
                "Length": fwd.length,
                "Tm (°C)": round(fwd.tm_c, 1),
                "GC (%)": round(fwd.gc_pct, 1),
                "Score": round(fwd.score, 2),
            },
            {
                "Type": "REV",
                "Role": rev.exon_name,
                "Primer (5'→3')": rev.seq_5to3,
                "Length": rev.length,
                "Tm (°C)": round(rev.tm_c, 1),
                "GC (%)": round(rev.gc_pct, 1),
                "Score": round(rev.score, 2),
            },
        ]

        st.dataframe(pd.DataFrame(rows), use_container_width=True)

        if probe is not None:
            st.subheader("TaqMan probe")
            st.write("Probe is reported 5'→3' and sits between primers.")
            st.dataframe(
                pd.DataFrame([{
                    "Probe (5'→3')": probe.seq_5to3,
                    "Length": probe.length,
                    "Tm (°C)": round(probe.tm_c, 1),
                    "GC (%)": round(probe.gc_pct, 1),
                    "Score": round(probe.score, 2),
                }]),
                use_container_width=True
            )

        st.subheader("Primer-BLAST link (NCBI)")
        url = primer_blast_url_pair(fwd.seq_5to3, rev.seq_5to3, "Homo sapiens")
        st.markdown(f"[Open in Primer-BLAST]({url})")
        st.info(BLAST_INSTRUCTIONS)

        st.subheader("Dimer check")
        print_dimer_report_pair(fwd.seq_5to3, rev.seq_5to3)

        add_footer()

    except Exception as e:
        st.error(str(e))
        add_footer()
