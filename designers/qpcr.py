# designers/qpcr.py
import streamlit as st
import primer_engine as st

from primer_engine import (
    clean_seq,
    design_qpcr_junction_pair,
    qpcr_amplicon_size_from_hits,
)

def render():
    # Always show something so "blank page" cannot happen
    st.title("qPCR primers")

    st.caption(
        "Input format: paste the spliced sequence and mark the exon-exon junction with '^'. "
        "Example: ...ACCTG^GTTACA..."
    )

    # Debug line: if you see this, render() is being executed
    st.info("qPCR tab loaded successfully.")

    with st.expander("qPCR rules this tool enforces", expanded=True):
        st.markdown(
            """
- Amplicon usually **70 to 200 bp**
- Primer length usually **18 to 24 bp**
- Primer Tm usually **~60 C** with **tight tolerance** (default 2 C)
- GC% usually **40 to 60%**
- Avoid long homopolymers (default max run 3)
- Avoid obvious self complementarity
- Primer pair Tm difference should be small (default 1 C)
- Junction spanning primer:
  - You can force **Forward** spans junction, **Reverse** spans junction, or **AUTO** (tries both)
- Chemistry:
  - **SYBR**: returns primer pair
  - **TAQMAN**: returns primer pair + probe between primers (if possible)
"""
        )

    # Inputs
    seq_with_marker = st.text_area(
        "Spliced template with junction marker '^'",
        height=160,
        placeholder="Paste sequence with ^ at the exon-exon junction",
    )

    col1, col2, col3 = st.columns(3)
    with col1:
        chemistry = st.selectbox("Chemistry", ["SYBR", "TAQMAN"], index=0)
    with col2:
        junction_primer = st.selectbox("Junction primer", ["AUTO", "FWD", "REV"], index=0)
    with col3:
        dimer_k = st.number_input("3' dimer check window (k)", min_value=3, max_value=8, value=4, step=1)

    st.subheader("Primer constraints")

    c1, c2, c3, c4 = st.columns(4)
    with c1:
        min_len = st.number_input("Min primer length", min_value=14, max_value=40, value=18, step=1)
    with c2:
        max_len = st.number_input("Max primer length", min_value=14, max_value=40, value=24, step=1)
    with c3:
        primer_tm_target = st.number_input("Target primer Tm (C)", min_value=45.0, max_value=75.0, value=60.0, step=0.5)
    with c4:
        primer_tm_tol = st.number_input("Primer Tm tolerance (C)", min_value=0.5, max_value=10.0, value=2.0, step=0.5)

    c5, c6, c7, c8 = st.columns(4)
    with c5:
        primer_gc_min = st.number_input("Primer GC min (%)", min_value=0.0, max_value=100.0, value=40.0, step=1.0)
    with c6:
        primer_gc_max = st.number_input("Primer GC max (%)", min_value=0.0, max_value=100.0, value=60.0, step=1.0)
    with c7:
        max_homopolymer = st.number_input("Max homopolymer run", min_value=2, max_value=8, value=3, step=1)
    with c8:
        max_tm_diff_pair = st.number_input("Max Tm diff in pair (C)", min_value=0.0, max_value=10.0, value=1.0, step=0.5)

    st.subheader("Amplicon constraints")
    a1, a2, a3 = st.columns(3)
    with a1:
        amplicon_min = st.number_input("Amplicon min (bp)", min_value=30, max_value=1000, value=70, step=5)
    with a2:
        amplicon_max = st.number_input("Amplicon max (bp)", min_value=30, max_value=1000, value=200, step=5)
    with a3:
        junction_min_overlap_each_side = st.number_input("Min overlap on each side of junction", min_value=3, max_value=12, value=6, step=1)

    st.divider()

    if st.button("Design qPCR primers"):
        try:
            if not seq_with_marker.strip():
                st.error("Paste a sequence first.")
                return

            # quick sanity check: should contain ^
            if "^" not in seq_with_marker:
                st.error("You must include '^' at the exon-exon junction in the sequence.")
                return

            # run design
            fwd, rev, probe = design_qpcr_junction_pair(
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

            amp = qpcr_amplicon_size_from_hits(fwd, rev)

            st.success("Primer design complete.")

            st.subheader("Results")
            st.write(f"Predicted amplicon size: **{amp} bp**")

            st.markdown("**Forward primer (5'->3')**")
            st.code(fwd.seq_5to3)
            st.write(f"Tm: {fwd.tm_c:.1f} C, GC: {fwd.gc_pct:.1f}%, start: {fwd.start_0based}, len: {fwd.length}, role: {fwd.exon_name}")

            st.markdown("**Reverse primer (5'->3')**")
            st.code(rev.seq_5to3)
            st.write(f"Tm: {rev.tm_c:.1f} C, GC: {rev.gc_pct:.1f}%, start: {rev.start_0based}, len: {rev.length}, role: {rev.exon_name}")

            if chemistry == "TAQMAN":
                st.subheader("Probe (TaqMan)")
                if probe is None:
                    st.warning("No probe returned. Try widening probe rules or increasing amplicon size.")
                else:
                    st.code(probe.seq_5to3)
                    st.write(f"Tm: {probe.tm_c:.1f} C, GC: {probe.gc_pct:.1f}%, start: {probe.start_0based}, len: {probe.length}")

        except Exception as e:
            st.error("qPCR designer crashed. The error is below.")
            st.exception(e)
