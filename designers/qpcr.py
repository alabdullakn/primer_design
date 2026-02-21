# designers/qpcr.py

import streamlit as st
import primer_engine


def render():
    st.header("qPCR primers")

    st.write(
        "For junction qPCR, paste sequence with a junction marker like this:\n\n"
        "`...EXON1_END^EXON2_START...`\n\n"
        "One primer will span the junction. Use cDNA sequence context."
    )

    colA, colB = st.columns(2)
    with colA:
        chemistry = st.radio("Chemistry", ["SYBR", "TAQMAN"], horizontal=True)
    with colB:
        junction_primer = st.radio("Junction primer", ["AUTO", "FWD", "REV"], horizontal=True)

    seq_with_marker = st.text_area(
        "Junction sequence (use ^ at exon-exon junction)",
        height=140,
        placeholder="Example: ATG...TTG^GAA...CCT"
    )

    with st.expander("qPCR strict rules (recommended defaults)", expanded=True):
        c1, c2, c3 = st.columns(3)
        with c1:
            min_len = st.number_input("Min primer length", 15, 30, 18, 1)
            max_len = st.number_input("Max primer length", 15, 35, 24, 1)
            max_homopolymer = st.number_input("Max homopolymer run", 2, 6, 3, 1)
        with c2:
            primer_tm_target = st.number_input("Primer Tm target (C)", 50.0, 70.0, 60.0, 0.5)
            primer_tm_tol = st.number_input("Primer Tm tolerance (C)", 0.5, 10.0, 2.0, 0.5)
            max_tm_diff_pair = st.number_input("Max Tm diff between primers (C)", 0.0, 10.0, 1.0, 0.5)
        with c3:
            primer_gc_min = st.number_input("Primer GC min (%)", 0.0, 100.0, 40.0, 1.0)
            primer_gc_max = st.number_input("Primer GC max (%)", 0.0, 100.0, 60.0, 1.0)
            dimer_k = st.number_input("3' dimer check window (k)", 3, 8, 4, 1)

        c4, c5 = st.columns(2)
        with c4:
            amplicon_min = st.number_input("Amplicon min (bp)", 30, 400, 70, 5)
        with c5:
            amplicon_max = st.number_input("Amplicon max (bp)", 30, 500, 200, 5)

        junction_min_overlap_each_side = st.number_input("Junction overlap each side (bp)", 3, 12, 6, 1)

    design = st.button("Design qPCR primers")

    if not design:
        return

    try:
        fwd, rev, probe = primer_engine.design_qpcr_junction_pair(
            seq_with_junction_marker=seq_with_marker,
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
            junction_min_overlap_each_side=int(junction_min_overlap_each_side),
            max_tm_diff_pair=float(max_tm_diff_pair),
            dimer_k=int(dimer_k),
        )

        amp_len = primer_engine.qpcr_amplicon_size(seq_with_marker, fwd, rev)

        st.subheader("Result")
        st.write(f"Amplicon size (bp): **{amp_len}**")

        st.table([
            {
                "Type": "Forward",
                "Seq (5'->3')": fwd.seq_5to3,
                "Tm (C)": round(fwd.tm_c, 2),
                "GC (%)": round(fwd.gc_pct, 2),
                "Start (0-based)": fwd.start_0based,
                "Length": fwd.length,
            },
            {
                "Type": "Reverse",
                "Seq (5'->3')": rev.seq_5to3,
                "Tm (C)": round(rev.tm_c, 2),
                "GC (%)": round(rev.gc_pct, 2),
                "Start (0-based)": rev.start_0based,
                "Length": rev.length,
            },
        ])

        primer_engine.print_dimer_report_pair(fwd.seq_5to3, rev.seq_5to3)

        if chemistry == "TAQMAN" and probe is not None:
            st.subheader("TaqMan probe")
            st.table([{
                "Probe (5'->3')": probe.seq_5to3,
                "Tm (C)": round(probe.tm_c, 2),
                "GC (%)": round(probe.gc_pct, 2),
                "Start (0-based)": probe.start_0based,
                "Length": probe.length,
            }])

        st.info(
            "Practical qPCR notes:\n"
            "1) Amplicon 70 to 200 bp is ideal.\n"
            "2) Keep primer Tm close (usually within 1 C).\n"
            "3) Try AUTO first because sometimes only REV spanning works.\n"
            "4) For SYBR, specificity matters a lot, so junction spanning helps reduce gDNA signal."
        )

    except Exception as e:
        st.error(str(e))
