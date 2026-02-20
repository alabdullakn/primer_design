# pages/qpcr_primers.py
import streamlit as st
import pandas as pd

from primer_engine import (
    design_qpcr_junction_pair,
    qpcr_amplicon_size_from_hits,
    print_dimer_report_pair,
)
from utils.blast import primer_blast_url_pair, primer_blast_url_single
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer


def render():
    st.title("qPCR primers")
    st.write("Design qPCR primers with strict rules and an exon-exon junction option.")

    st.markdown(
        """
**Input format**
- Paste the spliced cDNA sequence and mark the exon-exon junction with `^`
- Example: `...CTGAC^GTTCCA...`

**What this does**
- You can force the junction primer to be **Forward**, **Reverse**, or let it try **both** and pick the best
- Rules are stricter than regular PCR (short amplicon, tighter GC, tighter Tm matching)
- Optional: generate a TaqMan probe between primers
        """
    )

    st.markdown("---")

    # ============================
    # Chemistry + junction options
    # ============================

    c0, c1 = st.columns(2)
    with c0:
        chemistry = st.selectbox(
            "Chemistry",
            ["SYBR", "TAQMAN"],
            index=0,
            help="SYBR uses only primers. TAQMAN also designs a probe between primers.",
            key="qpcr_chem",
        )
    with c1:
        junction_mode = st.selectbox(
            "Junction primer",
            ["AUTO", "FWD", "REV"],
            index=0,
            help="AUTO tries both FWD-junction and REV-junction and returns the best pair.",
            key="qpcr_jmode",
        )

    st.markdown("---")

    # ============================
    # Strict defaults (editable)
    # ============================

    with st.expander("qPCR rules (strict defaults)", expanded=False):
        cA, cB, cC = st.columns(3)

        with cA:
            min_len = st.number_input("Min primer length", 16, 30, 18, key="qpcr_min_len")
            max_len = st.number_input("Max primer length", 18, 40, 24, key="qpcr_max_len")
            max_homopolymer = st.number_input("Max homopolymer run", 2, 6, 3, key="qpcr_max_hpoly")
            dimer_k = st.number_input("3' dimer screen k", 3, 10, 4, key="qpcr_dimer_k")

        with cB:
            primer_tm_target = st.number_input("Primer target Tm (°C)", 50.0, 70.0, 60.0, key="qpcr_tm_tgt")
            primer_tm_tol = st.number_input("Primer Tm tolerance (±°C)", 0.5, 8.0, 2.0, key="qpcr_tm_tol")
            max_tm_diff_pair = st.number_input("Max Tm diff (FWD vs REV)", 0.0, 5.0, 1.0, key="qpcr_tm_pairdiff")
            junction_min_overlap = st.number_input("Min overlap each side of junction", 4, 12, 6, key="qpcr_jover")

        with cC:
            primer_gc_min = st.number_input("Primer GC min (%)", 20.0, 60.0, 40.0, key="qpcr_gcmin")
            primer_gc_max = st.number_input("Primer GC max (%)", 40.0, 80.0, 60.0, key="qpcr_gcmax")
            amplicon_min = st.number_input("Amplicon min (bp)", 40, 250, 70, key="qpcr_ampmin")
            amplicon_max = st.number_input("Amplicon max (bp)", 60, 400, 200, key="qpcr_ampmax")

    st.markdown("---")

    # ============================
    # Sequence input
    # ============================

    seq = st.text_area(
        "Spliced cDNA sequence with one junction marker ^",
        height=220,
        placeholder="Paste spliced cDNA and put ^ at the exon-exon junction\nExample:\n...AAGGACCTGATGCTGAC^GTTCCAGGAGTCTGACT...",
        key="qpcr_seq",
    )

    st.markdown("---")

    run = st.button("Design qPCR primers", key="qpcr_run")
    if not run:
        add_footer()
        return

    if not seq or "^" not in seq:
        st.error("Paste the spliced sequence and include exactly one '^' at the exon-exon junction.")
        add_footer()
        return

    if seq.count("^") != 1:
        st.error("Use exactly ONE '^' marker. Pick one junction only.")
        add_footer()
        return

    # ============================
    # Run design
    # ============================

    try:
        fwd, rev, probe = design_qpcr_junction_pair(
            seq_with_junction_marker=seq,
            chemistry=chemistry,
            junction_primer=junction_mode,
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

        st.success("qPCR design complete.")
        st.caption(SCORE_EXPLANATION)

        # ============================
        # Output table
        # ============================

        rows = [
            {
                "Type": "FWD",
                "Primer (5'→3')": fwd.seq_5to3,
                "Start (0-based)": fwd.start_0based,
                "End (0-based, exclusive)": fwd.start_0based + fwd.length,
                "Length": fwd.length,
                "Tm (°C)": round(fwd.tm_c, 1),
                "GC (%)": round(fwd.gc_pct, 1),
                "Score": round(fwd.score, 2),
                "Note": fwd.exon_name,
            },
            {
                "Type": "REV",
                "Primer (5'→3')": rev.seq_5to3,
                "Start (0-based)": rev.start_0based,
                "End (0-based, exclusive)": rev.start_0based + rev.length,
                "Length": rev.length,
                "Tm (°C)": round(rev.tm_c, 1),
                "GC (%)": round(rev.gc_pct, 1),
                "Score": round(rev.score, 2),
                "Note": rev.exon_name,
            },
        ]

        st.dataframe(pd.DataFrame(rows), use_container_width=True)
        st.write(f"Estimated amplicon length (spliced template): **{amp_len} bp**")

        if chemistry == "TAQMAN" and probe is not None:
            st.subheader("TaqMan probe (between primers)")
            probe_rows = [{
                "Type": "PROBE",
                "Sequence (5'→3')": probe.seq_5to3,
                "Start (0-based)": probe.start_0based,
                "End (0-based, exclusive)": probe.start_0based + probe.length,
                "Length": probe.length,
                "Tm (°C)": round(probe.tm_c, 1),
                "GC (%)": round(probe.gc_pct, 1),
                "Score": round(probe.score, 2),
            }]
            st.dataframe(pd.DataFrame(probe_rows), use_container_width=True)

        st.markdown("---")

        # ============================
        # BLAST links
        # ============================

        org = "Homo sapiens"
        st.subheader("Primer-BLAST links (NCBI)")

        st.markdown(f"[Open pair in Primer-BLAST]({primer_blast_url_pair(fwd.seq_5to3, rev.seq_5to3, org)})")
        st.info(BLAST_INSTRUCTIONS)

        with st.expander("Single-primer Primer-BLAST links"):
            st.markdown(f"**FWD**: [Primer-BLAST]({primer_blast_url_single(fwd.seq_5to3, org)})")
            st.markdown(f"**REV**: [Primer-BLAST]({primer_blast_url_single(rev.seq_5to3, org)})")

        st.markdown("---")

        # ============================
        # Dimer check
        # ============================

        st.subheader("Heterodimer check")
        print_dimer_report_pair(fwd.seq_5to3, rev.seq_5to3)

        add_footer()

    except Exception as e:
        st.error(str(e))
        add_footer()
