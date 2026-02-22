# designers/qpcr.py
import streamlit as st
import pandas as pd

from qpcr_engine import (
    amplicon_size_from_hits,
    design_qpcr_junction_pair,
)
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from utils.blast import primer_blast_url_pair, primer_blast_url_single

def render():
    st.title("qPCR primers")
    
    with st.expander("Primer parameters", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            chemistry_label = st.selectbox(
                "Chemistry",
                ["SYBR Green", "TaqMan probe"],
                index=0,
                help="SYBR Green designs primer pairs only. TaqMan also designs an internal hydrolysis probe.",
                key="qpcr_chemistry",
            )
        with c2:
            junction_primer = st.selectbox(
                "Junction spanning primer",
                ["AUTO", "FWD", "REV"],
                index=0,
                help="Force the forward or reverse primer to span the junction, or let the algorithm choose.",
                key="qpcr_junction_primer",
            )
        with c3:
            dimer_k = st.number_input("3' dimer check window (k)", 3, 10, 4, key="qpcr_dimer_k")

        st.subheader("Primer constraints")
        r1, r2 = st.columns(2)
        with r1:
            min_len = st.number_input("Min length", 16, 28, 18, key="qpcr_min_len")
            tm_target = st.number_input("Tm target (°C)", 50.0, 70.0, 60.0, key="qpcr_tm_target")
            gc_min = st.number_input("GC min percent", 20.0, 80.0, 40.0, key="qpcr_gc_min")
            amp_min = st.number_input("Amplicon min bp", 40, 300, 70, key="qpcr_amp_min")
            j_ov = st.number_input("Min overlap each side of junction", 3, 12, 6, key="qpcr_j_ov")

        with r2:
            max_len = st.number_input("Max length", 16, 32, 24, key="qpcr_max_len")
            tm_tol = st.number_input("Tm tolerance (± °C)", 0.5, 10.0, 2.0, key="qpcr_tm_tol")
            gc_max = st.number_input("GC max percent", 20.0, 80.0, 60.0, key="qpcr_gc_max")
            amp_max = st.number_input("Amplicon max bp", 60, 500, 200, key="qpcr_amp_max")
            junction_3p_max_distance = st.number_input(
                "Max distance of junction from 3' end (strict)",
                1,
                10,
                8,
                help="Smaller values are stricter and force the junction to be closer to the primer 3' end.",
                key="qpcr_junction_3p_max_distance",
            )

        c4, _ = st.columns([1, 1])
        with c4:
            max_hpoly = st.number_input("Max homopolymer run", 2, 6, 3, key="qpcr_max_hpoly")
            max_tm_diff = st.number_input("Max Tm difference pair", 0.0, 10.0, 1.0, key="qpcr_max_tm_diff")

        chemistry = "TAQMAN" if chemistry_label == "TaqMan probe" else "SYBR"

        probe_tm_target = 69.0
        probe_tm_tol = 3.0
        probe_min_len = 18
        probe_max_len = 30

        if chemistry == "TAQMAN":
            st.subheader("TaqMan probe constraints")
            p1, p2 = st.columns(2)
            with p1:
                probe_tm_target = st.number_input("Probe Tm target", 60.0, 78.0, 69.0, key="qpcr_probe_tm_target")
                probe_min_len = st.number_input("Probe min length", 16, 28, 18, key="qpcr_probe_min_len")
            with p2:
                probe_tm_tol = st.number_input("Probe Tm tolerance", 0.5, 8.0, 3.0, key="qpcr_probe_tm_tol")
                probe_max_len = st.number_input("Probe max length", 18, 35, 30, key="qpcr_probe_max_len")

    st.write("Paste your sequence and mark the exon-exon junction with ^")
    st.write("Example: ...ACCTG^GTTCA...")

    seq = st.text_area(
        "Sequence with junction marker ^",
        height=180,
        placeholder="Paste DNA sequence here with ^ at the junction",
    )

    if junction_3p_max_distance < j_ov:
        st.warning("Current settings are contradictory: max 3' distance must be at least min overlap each side.")

    st.divider()

    if st.button("Design qPCR primers"):
        try:
            if not seq or "^" not in seq:
                st.error("Your sequence must include the junction marker ^")
                return

            fwd, rev, probe = design_qpcr_junction_pair(
                seq_with_junction_marker=seq,
                chemistry=chemistry,
                junction_primer=junction_primer,
                min_len=int(min_len),
                max_len=int(max_len),
                primer_tm_target=float(tm_target),
                primer_tm_tol=float(tm_tol),
                primer_gc_min=float(gc_min),
                primer_gc_max=float(gc_max),
                max_homopolymer=int(max_hpoly),
                amplicon_min=int(amp_min),
                amplicon_max=int(amp_max),
                junction_min_overlap_each_side=int(j_ov),
                junction_max_3prime_distance=int(junction_3p_max_distance),
                max_tm_diff_pair=float(max_tm_diff),
                dimer_k=int(dimer_k),
                probe_tm_target=float(probe_tm_target),
                probe_tm_tol=float(probe_tm_tol),
                probe_min_len=int(probe_min_len),
                probe_max_len=int(probe_max_len),
            )

            amp_size = amplicon_size_from_hits(fwd, rev)

            st.success("Designed qPCR primers")
            st.write(f"Chemistry: {chemistry_label}")
            st.caption(SCORE_EXPLANATION)

            rows = [
                {
                    "Type": "FWD",
                    "Primer (5'→3')": fwd.seq_5to3,
                    "Length": fwd.length,
                    "Tm (°C)": round(fwd.tm_c, 1),
                    "GC (%)": round(fwd.gc_pct, 1),
                    "Score": round(fwd.score, 2),
                    "Start (0-based)": fwd.start_0based,
                    "Role": fwd.role,
                },
                {
                    "Type": "REV",
                    "Primer (5'→3')": rev.seq_5to3,
                    "Length": rev.length,
                    "Tm (°C)": round(rev.tm_c, 1),
                    "GC (%)": round(rev.gc_pct, 1),
                    "Score": round(rev.score, 2),
                    "Start (0-based)": rev.start_0based,
                    "Role": rev.role,
                },
            ]
            
            if probe is not None:
                rows.append(
                    {
                        "Type": "PROBE",
                        "Primer (5'→3')": probe.seq_5to3,
                        "Length": probe.length,
                        "Tm (°C)": round(probe.tm_c, 1),
                        "GC (%)": round(probe.gc_pct, 1),
                        "Score": round(probe.score, 2),
                        "Start (0-based)": probe.start_0based,
                        "Role": "Internal probe",
                    }
                )

            st.dataframe(pd.DataFrame(rows), use_container_width=True)
            st.write(f"Estimated amplicon length: **{amp_size} bp**")

            st.subheader("Primer-BLAST links (NCBI)")
            organism = "Homo sapiens"
            pair_url = primer_blast_url_pair(fwd.seq_5to3, rev.seq_5to3, organism)
            st.markdown(f"**qPCR pair**: [Open in Primer-BLAST]({pair_url})")

            with st.expander("Single-primer Primer-BLAST links"):
                st.markdown(
                    f"**Forward primer**: [Primer-BLAST]({primer_blast_url_single(fwd.seq_5to3, organism)})"
                )
                st.markdown(
                    f"**Reverse primer**: [Primer-BLAST]({primer_blast_url_single(rev.seq_5to3, organism)})"
                )
                if probe is not None:
                    st.markdown(
                        f"**Probe**: [Primer-BLAST]({primer_blast_url_single(probe.seq_5to3, organism)})"
                    )

            st.info(BLAST_INSTRUCTIONS)
        except Exception as e:
            st.error(f"qPCR design failed: {e}")
            st.info("Try widening Tm tolerance/GC range, relaxing strict junction distance, or increasing amplicon max")
