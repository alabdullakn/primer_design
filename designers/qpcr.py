import streamlit as st
import primer_engine

def render():
    st.title("qPCR primers")

    try:
        st.markdown("### Input")
        st.info("Paste sequence with the exon junction marked using ^  (example: AAAACCC^GGGTTT).")

        seq = st.text_area("Junction-marked template", height=160, placeholder="...AAAA^TTTT...")

        col1, col2, col3 = st.columns(3)
        with col1:
            chemistry = st.selectbox("Chemistry", ["SYBR", "TAQMAN"], index=0)
        with col2:
            junction_primer = st.selectbox("Junction primer", ["AUTO", "FWD", "REV"], index=0)
        with col3:
            amplicon_min = st.number_input("Amplicon min (bp)", min_value=40, max_value=300, value=70, step=1)

        amplicon_max = st.number_input("Amplicon max (bp)", min_value=60, max_value=400, value=200, step=1)

        st.markdown("### Primer rules")
        c1, c2, c3 = st.columns(3)
        with c1:
            min_len = st.number_input("Min primer length", min_value=16, max_value=30, value=18, step=1)
            max_len = st.number_input("Max primer length", min_value=18, max_value=35, value=24, step=1)
        with c2:
            tm_target = st.number_input("Primer Tm target (C)", min_value=50.0, max_value=75.0, value=60.0, step=0.5)
            tm_tol = st.number_input("Primer Tm tolerance (C)", min_value=0.5, max_value=10.0, value=2.0, step=0.5)
        with c3:
            dimer_k = st.number_input("3' dimer check window (k)", min_value=3, max_value=8, value=4, step=1)
            max_tm_diff = st.number_input("Max Tm diff pair (C)", min_value=0.0, max_value=5.0, value=1.0, step=0.5)

        st.markdown("### GC and homopolymer")
        g1, g2, g3 = st.columns(3)
        with g1:
            gc_min = st.number_input("Primer GC min (%)", min_value=20.0, max_value=70.0, value=40.0, step=1.0)
        with g2:
            gc_max = st.number_input("Primer GC max (%)", min_value=30.0, max_value=80.0, value=60.0, step=1.0)
        with g3:
            max_hpoly = st.number_input("Max homopolymer run", min_value=2, max_value=6, value=3, step=1)

        st.markdown("### Junction overlap")
        j1, j2 = st.columns(2)
        with j1:
            j_overlap = st.number_input("Min overlap each side (bp)", min_value=4, max_value=12, value=6, step=1)

        if st.button("Design qPCR primers", type="primary", disabled=(not seq.strip())):
            fwd, rev, probe = primer_engine.design_qpcr_junction_pair(
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
                amplicon_min=int(amplicon_min),
                amplicon_max=int(amplicon_max),
                junction_min_overlap_each_side=int(j_overlap),
                max_tm_diff_pair=float(max_tm_diff),
                dimer_k=int(dimer_k),
            )

            st.success("Found a qPCR primer set")

            amp = primer_engine.qpcr_amplicon_size(seq, fwd, rev)

            st.markdown("### Output")
            st.write({
                "FWD (5'->3')": fwd.seq_5to3,
                "FWD Tm": round(fwd.tm_c, 2),
                "FWD GC%": round(fwd.gc_pct, 1),
                "REV (5'->3')": rev.seq_5to3,
                "REV Tm": round(rev.tm_c, 2),
                "REV GC%": round(rev.gc_pct, 1),
                "Amplicon (bp)": amp,
                "Junction spans": junction_primer
            })

            if chemistry == "TAQMAN":
                st.markdown("### Probe")
                if probe is None:
                    st.warning("No probe returned.")
                else:
                    st.write({
                        "Probe (5'->3')": probe.seq_5to3,
                        "Probe Tm": round(probe.tm_c, 2),
                        "Probe GC%": round(probe.gc_pct, 1),
                        "Probe start (0-based)": probe.start_0based
                    })

    except Exception as e:
        st.error("qPCR tab crashed")
        st.exception(e)
