# designers/qpcr_primers.py
import streamlit as st
import pandas as pd

import primer_engine

from utils.blast import primer_blast_url_pair
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer

    
def render():
    st.title("qPCR primers")

    st.write(
        "Paste spliced cDNA and mark the exon exon junction using ^. "
        "Example: ...CTGAC^GTTCCA..."
    )

    st.markdown("---")

    c1, c2 = st.columns(2)
    with c1:
        chemistry = st.selectbox("Chemistry", ["SYBR", "TAQMAN"], index=0, key="qpcr_chemistry")
    with c2:
        junction_primer = st.selectbox(
            "Junction primer",
            ["AUTO", "FWD", "REV"],
            index=0,
            help="AUTO tries both and picks the best valid pair.",
            key="qpcr_junction_primer",
        )

    with st.expander("Strict qPCR parameters", expanded=False):
        a, b, c = st.columns(3)

        with a:
            min_len = st.number_input("Min primer length", 16, 30, 18, key="qpcr_min_len")
            max_len = st.number_input("Max primer length", 18, 40, 24, key="qpcr_max_len")
            max_homopolymer = st.number_input("Max homopolymer run", 2, 6, 3, key="qpcr_max_hpoly")
            dimer_k = st.number_input("3' dimer screen k", 3, 10, 4, key="qpcr_dimer_k")

        with b:
            primer_tm_target = st.number_input("Primer target Tm (C)", 50.0, 70.0, 60.0, key="qpcr_tm_tgt")
            primer_tm_tol = st.number_input("Primer Tm tolerance (plus minus)", 0.5, 8.0, 2.0, key="qpcr_tm_tol")
            max_tm_diff_pair = st.number_input("Max Tm diff between primers", 0.0, 5.0, 1.0, key="qpcr_tm_pairdiff")
            junction_min_overlap = st.number_input("Min overlap each side of junction", 4, 12, 6, key="qpcr_jover")

        with c:
            primer_gc_min = st.number_input("Primer GC min percent", 20.0, 60.0, 40.0, key="qpcr_gcmin")
            primer_gc_max = st.number_input("Primer GC max percent", 40.0, 80.0, 60.0, key="qpcr_gcmax")
            amplicon_min = st.number_input("Amplicon min bp", 40, 250, 70, key="qpcr_ampmin")
            amplicon_max = st.number_input("Amplicon max bp", 60, 400, 200, key="qpcr_ampmax")

        if chemistry == "TAQMAN":
            st.markdown("### Probe parameters (TaqMan)")
            p1, p2, p3 = st.columns(3)
            with p1:
                probe_min_len = st.number_input("Probe min length", 15, 40, 18, key="probe_min_len")
                probe_max_len = st.number_input("Probe max length", 15, 60, 30, key="probe_max_len")
            with p2:
                probe_tm_target = st.number_input("Probe target Tm (C)", 55.0, 80.0, 68.0, key="probe_tm_tgt")
                probe_tm_tol = st.number_input("Probe Tm tolerance", 1.0, 10.0, 5.0, key="probe_tm_tol")
            with p3:
                probe_gc_min = st.number_input("Probe GC min percent", 20.0, 80.0, 30.0, key="probe_gc_min")
                probe_gc_max = st.number_input("Probe GC max percent", 20.0, 90.0, 80.0, key="probe_gc_max")
        else:
            probe_min_len = 18
            probe_max_len = 30
            probe_tm_target = 68.0
            probe_tm_tol = 5.0
            probe_gc_min = 30.0
            probe_gc_max = 80.0

    st.markdown("---")

    seq = st.text_area(
        "Spliced cDNA sequence with one ^ marker",
        height=220,
        key="qpcr_seq",
        placeholder="...CTGAC^GTTCCA...",
    )

    st.markdown("---")

    if not st.button("Design qPCR primers", key="qpcr_run"):
        add_footer()
        return

    if not seq or "^" not in seq:
        st.error("You must include exactly one ^ marker at the exon exon junction.")
        add_footer()
        return

    if seq.count("^") != 1:
        st.error("Use exactly one ^ marker.")
        add_footer()
        return

    try:
        fwd, rev, probe = primer_engine.design_qpcr_junction_pair(
            seq_with_junction_marker=seq,
            chemistry=chemistry,
            junction_primer=junction_primer,
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
            probe_min_len=int(probe_min_len),
            probe_max_len=int(probe_max_len),
            probe_tm_target=float(probe_tm_target),
            probe_tm_tol=float(probe_tm_tol),
            probe_gc_min=float(probe_gc_min),
            probe_gc_max=float(probe_gc_max),
        )

        amp_len = primer_engine.qpcr_amplicon_size_from_hits(fwd, rev)

        st.success("qPCR primers designed.")
        st.caption(SCORE_EXPLANATION)

        rows = [
            {
                "Type": "FWD",
                "Sequence (5 to 3)": fwd.seq_5to3,
                "Length": fwd.length,
                "Tm": round(fwd.tm_c, 1),
                "GC percent": round(fwd.gc_pct, 1),
                "Start": fwd.start_0based,
                "Score": round(fwd.score, 2),
                "Note": fwd.exon_name,
            },
            {
                "Type": "REV",
                "Sequence (5 to 3)": rev.seq_5to3,
                "Length": rev.length,
                "Tm": round(rev.tm_c, 1),
                "GC percent": round(rev.gc_pct, 1),
                "Start": rev.start_0based,
                "Score": round(rev.score, 2),
                "Note": rev.exon_name,
            },
        ]

        st.dataframe(pd.DataFrame(rows), use_container_width=True)
        st.write(f"Amplicon length (spliced): **{amp_len} bp**")

        if chemistry == "TAQMAN" and probe is not None:
            st.subheader("Probe (TaqMan)")
            st.dataframe(
                pd.DataFrame([{
                    "Type": "PROBE",
                    "Sequence (5 to 3)": probe.seq_5to3,
                    "Length": probe.length,
                    "Tm": round(probe.tm_c, 1),
                    "GC percent": round(probe.gc_pct, 1),
                    "Start": probe.start_0based,
                    "Score": round(probe.score, 2),
                }]),
                use_container_width=True
            )

        st.markdown("---")
        st.subheader("Primer BLAST link")
        url = primer_blast_url_pair(fwd.seq_5to3, rev.seq_5to3, "Homo sapiens")
        st.markdown(f"[Open in Primer BLAST]({url})")
        st.info(BLAST_INSTRUCTIONS)

        st.markdown("---")
        st.subheader("Heterodimer check")
        primer_engine.print_dimer_report_pair(fwd.seq_5to3, rev.seq_5to3)

        add_footer()

    except Exception as e:
        st.error(str(e))
        add_footer()
