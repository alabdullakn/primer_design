# designers/qpcr.py
import streamlit as st

from primer_engine import (
    clean_seq,
    design_qpcr_junction_pair,
    qpcr_amplicon_size,
)


def _render_rules():
    with st.expander("qPCR rules (strict, practical)", expanded=True):
        st.markdown(
            """
**Core rules**
- Amplicon: **70–200 bp** (typical best range: 80–150 bp)
- Primer length: **18–24 nt**
- Primer Tm: **~60 C** (keep pair within **1 C** if possible)
- GC: **40–60%**
- Avoid homopolymers: **no runs > 3**
- Avoid strong self-complementarity and 3' complementarity (primer dimers)

**Junction rule**
- For cDNA qPCR, at least one primer should **span the exon-exon junction** (this reduces gDNA amplification)
- Sometimes the forward spanning primer fails, so this tool supports:
  - **FWD spans junction**
  - **REV spans junction**
  - **AUTO** (tries both and picks best)

**Chemistry**
- **SYBR Green**: just primers, specificity comes from primer design + melt curve
- **TaqMan probe**: primers + a probe inside the amplicon
  - Probe Tm should be **higher than primers** (about **68–70 C**)
  - Avoid probe starting with **G** (common practical rule)
"""
        )


def render():
    st.header("qPCR primers")
    st.caption("Paste sequence with '^' at the exon-exon junction. Example: ...AACTG^GTTCA...")

    _render_rules()

    col1, col2 = st.columns([2, 1], vertical_alignment="top")

    with col1:
        seq_in = st.text_area(
            "Spliced template with junction marker '^'",
            height=160,
            placeholder="Paste exon1_end...^...exon2_start (include at least 30 bases on each side if you can)",
        )

    with col2:
        chemistry = st.selectbox("Chemistry", ["SYBR", "TAQMAN"], index=0)
        junction_mode = st.selectbox("Which primer spans the junction?", ["AUTO", "FWD", "REV"], index=0)

        st.markdown("### Primer constraints")
        min_len, max_len = st.slider("Primer length (nt)", 16, 30, (18, 24), step=1)
        primer_tm_target = st.number_input("Primer target Tm (C)", value=60.0, min_value=45.0, max_value=75.0, step=0.5)
        primer_tm_tol = st.number_input("Primer Tm tolerance (C)", value=2.0, min_value=0.5, max_value=10.0, step=0.5)
        gc_min, gc_max = st.slider("GC percent", 20, 80, (40, 60), step=1)
        max_hpoly = st.selectbox("Max homopolymer run", [3, 4, 5], index=0)
        max_tm_diff = st.number_input("Max Tm difference (pair) (C)", value=1.0, min_value=0.0, max_value=10.0, step=0.5)

        st.markdown("### Amplicon constraints")
        amp_min, amp_max = st.slider("Amplicon size (bp)", 40, 400, (70, 200), step=5)

        st.markdown("### Dimer filter")
        dimer_k = st.selectbox("3' dimer check window (k)", [3, 4, 5, 6], index=1)

        st.markdown("### Junction overlap")
        j_overlap = st.slider("Min overlap on each side of junction (nt)", 4, 10, 6, step=1)

        st.markdown("### Probe (TaqMan only)")
        probe_tm_target = st.number_input("Probe target Tm (C)", value=69.0, min_value=55.0, max_value=80.0, step=0.5)
        probe_tm_tol = st.number_input("Probe Tm tolerance (C)", value=3.0, min_value=0.5, max_value=10.0, step=0.5)
        probe_min_len, probe_max_len = st.slider("Probe length (nt)", 12, 40, (18, 30), step=1)

    st.divider()

    if st.button("Design qPCR primers", type="primary"):
        try:
            if not seq_in or "^" not in seq_in:
                raise ValueError("You must include '^' to mark the exon-exon junction in the input sequence.")

            # This also cleans A/C/G/T and will raise if empty
            _ = clean_seq(seq_in.replace("^", ""))

            fwd, rev, probe = design_qpcr_junction_pair(
                seq_with_junction_marker=seq_in,
                chemistry=chemistry,
                junction_primer=junction_mode,
                min_len=min_len,
                max_len=max_len,
                primer_tm_target=primer_tm_target,
                primer_tm_tol=primer_tm_tol,
                primer_gc_min=float(gc_min),
                primer_gc_max=float(gc_max),
                max_homopolymer=int(max_hpoly),
                amplicon_min=int(amp_min),
                amplicon_max=int(amp_max),
                junction_min_overlap_each_side=int(j_overlap),
                max_tm_diff_pair=float(max_tm_diff),
                dimer_k=int(dimer_k),
                probe_min_len=int(probe_min_len),
                probe_max_len=int(probe_max_len),
                probe_tm_target=float(probe_tm_target),
                probe_tm_tol=float(probe_tm_tol),
            )

            amp_size = qpcr_amplicon_size(seq_in, fwd, rev)

            st.success("Designed qPCR set")
            c1, c2 = st.columns(2)

            with c1:
                st.subheader("Forward primer (5'->3')")
                st.code(fwd.seq_5to3)
                st.write({"Tm (C)": round(fwd.tm_c, 2), "GC %": round(fwd.gc_pct, 1), "Length": fwd.length, "Start (0-based)": fwd.start_0based})

            with c2:
                st.subheader("Reverse primer (5'->3')")
                st.code(rev.seq_5to3)
                st.write({"Tm (C)": round(rev.tm_c, 2), "GC %": round(rev.gc_pct, 1), "Length": rev.length, "Bind start (0-based)": rev.start_0based})

            st.subheader("Amplicon")
            st.write({"Amplicon size (bp)": int(amp_size)})

            if chemistry == "TAQMAN":
                if probe is None:
                    st.warning("TaqMan selected but no probe returned.")
                else:
                    st.subheader("Probe (5'->3')")
                    st.code(probe.seq_5to3)
                    st.write({"Tm (C)": round(probe.tm_c, 2), "GC %": round(probe.gc_pct, 1), "Length": probe.length, "Start (0-based)": probe.start_0based})

        except Exception as e:
            st.error(str(e))
            st.info("Try: widen primer Tm tolerance, widen amplicon range, or reduce junction overlap requirements.")
