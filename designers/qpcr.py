# designers/qpcr.py
import streamlit as st
import primer_engine

def _clean_keep_marker(seq: str) -> str:
    """
    Keep A/C/G/T and '^' only (so user can paste junction-marked sequence).
    """
    seq = (seq or "").upper()
    out = []
    for ch in seq:
        if ch in "ACGT^":
            out.append(ch)
    return "".join(out)

def render():
    # Always show something so you can confirm the tab is working
    st.markdown("## qPCR primers")
    st.caption("Loaded designers/qpcr.py")

    st.info(
        "How to use: paste a cDNA junction sequence and mark the exon-exon junction with '^'.\n"
        "Example:  ...EXON1_LAST_BASES^EXON2_FIRST_BASES...\n\n"
        "qPCR is stricter than regular PCR: short amplicon (70-200 bp), tight Tm match, low homopolymers."
    )

    with st.expander("qPCR strict rules (what this tool enforces)", expanded=False):
        st.markdown(
            "- Amplicon: typically **70 to 200 bp**\n"
            "- Primer length: typically **18 to 24/25 nt**\n"
            "- Primer Tm: typically **~60 C** with **tight tolerance** (often 1 to 2 C)\n"
            "- Pair Tm difference: typically **<= 1 C**\n"
            "- GC: typically **40 to 60%**\n"
            "- Homopolymers: ideally **<= 3**\n"
            "- One primer should span the junction for isoform specificity. Some junctions work better with FWD spanning, others with REV spanning, so AUTO tries both."
        )

    # Chemistry choice
    colA, colB = st.columns(2)
    with colA:
        chemistry = st.selectbox(
            "Chemistry",
            ["SYBR", "TAQMAN"],
            index=0,
            help="SYBR: primers only. TAQMAN: primers + an internal probe between them."
        )
    with colB:
        junction_mode = st.selectbox(
            "Which primer spans the junction",
            ["AUTO", "FWD", "REV"],
            index=0,
            help="AUTO tries both spanning options and picks best. Use FWD or REV if you need to force it."
        )

    st.markdown("### Junction-marked input (cDNA)")
    seq_in = st.text_area(
        "Paste sequence with '^' at the exon-exon junction",
        height=180,
        placeholder="Example: ACTG...TTGCA^GGTAA...CTGA",
    )
    seq_in = _clean_keep_marker(seq_in)

    if seq_in and "^" not in seq_in:
        st.error("You must include '^' to mark the exon-exon junction.")
        return

    with st.expander("Parameters (strict defaults)", expanded=True):
        c1, c2, c3 = st.columns(3)
        with c1:
            min_len = st.number_input("Min primer length", 15, 30, 18, 1)
            max_len = st.number_input("Max primer length", 15, 35, 24, 1)
            max_hpoly = st.number_input("Max homopolymer run", 2, 6, 3, 1)
        with c2:
            tm_target = st.number_input("Primer Tm target (C)", 50.0, 70.0, 60.0, 0.5)
            tm_tol = st.number_input("Primer Tm tolerance (+/- C)", 0.5, 10.0, 2.0, 0.5)
            max_tm_diff = st.number_input("Max Tm difference in pair (C)", 0.0, 10.0, 1.0, 0.5)
        with c3:
            gc_min = st.number_input("Primer GC min (%)", 0.0, 100.0, 40.0, 1.0)
            gc_max = st.number_input("Primer GC max (%)", 0.0, 100.0, 60.0, 1.0)
            dimer_k = st.number_input("3' dimer check window (k)", 3, 8, 4, 1)

        d1, d2, d3 = st.columns(3)
        with d1:
            amp_min = st.number_input("Amplicon min (bp)", 40, 400, 70, 5)
        with d2:
            amp_max = st.number_input("Amplicon max (bp)", 60, 600, 200, 5)
        with d3:
            j_overlap = st.number_input("Min overlap each side of junction (bp)", 4, 12, 6, 1)

        # Probe params only if TAQMAN
        probe_min_len = probe_max_len = probe_tm_target = probe_tm_tol = probe_gc_min = probe_gc_max = None
        if chemistry == "TAQMAN":
            st.markdown("#### TaqMan probe settings")
            p1, p2, p3 = st.columns(3)
            with p1:
                probe_min_len = st.number_input("Probe min length", 12, 40, 18, 1)
                probe_max_len = st.number_input("Probe max length", 12, 60, 30, 1)
            with p2:
                probe_tm_target = st.number_input("Probe Tm target (C)", 55.0, 80.0, 69.0, 0.5)
                probe_tm_tol = st.number_input("Probe Tm tolerance (+/- C)", 0.5, 15.0, 3.0, 0.5)
            with p3:
                probe_gc_min = st.number_input("Probe GC min (%)", 0.0, 100.0, 30.0, 1.0)
                probe_gc_max = st.number_input("Probe GC max (%)", 0.0, 100.0, 80.0, 1.0)

    if st.button("Design qPCR primers", type="primary", disabled=not bool(seq_in)):
        try:
            if chemistry == "TAQMAN":
                fwd, rev, probe = primer_engine.design_qpcr_junction_pair(
                    seq_with_junction_marker=seq_in,
                    chemistry="TAQMAN",
                    junction_primer=junction_mode,
                    min_len=int(min_len),
                    max_len=int(max_len),
                    primer_tm_target=float(tm_target),
                    primer_tm_tol=float(tm_tol),
                    primer_gc_min=float(gc_min),
                    primer_gc_max=float(gc_max),
                    max_homopolymer=int(max_hpoly),
                    amplicon_min=int(amp_min),
                    amplicon_max=int(amp_max),
                    junction_min_overlap_each_side=int(j_overlap),
                    max_tm_diff_pair=float(max_tm_diff),
                    dimer_k=int(dimer_k),
                )

                # If your primer_engine supports probe tuning through args, keep this simple for now.
                # You can expose those later by wiring the probe params into primer_engine.

            else:
                fwd, rev, probe = primer_engine.design_qpcr_junction_pair(
                    seq_with_junction_marker=seq_in,
                    chemistry="SYBR",
                    junction_primer=junction_mode,
                    min_len=int(min_len),
                    max_len=int(max_len),
                    primer_tm_target=float(tm_target),
                    primer_tm_tol=float(tm_tol),
                    primer_gc_min=float(gc_min),
                    primer_gc_max=float(gc_max),
                    max_homopolymer=int(max_hpoly),
                    amplicon_min=int(amp_min),
                    amplicon_max=int(amp_max),
                    junction_min_overlap_each_side=int(j_overlap),
                    max_tm_diff_pair=float(max_tm_diff),
                    dimer_k=int(dimer_k),
                )

            amp_len = (rev.start_0based + rev.length) - fwd.start_0based

            st.success("Designed qPCR set")

            out_rows = [
                {
                    "Type": "FWD",
                    "Spans junction?": "YES" if fwd.exon_name == "Junction" else "NO",
                    "Seq (5'->3')": fwd.seq_5to3,
                    "Len": fwd.length,
                    "Tm (C)": round(fwd.tm_c, 2),
                    "GC %": round(fwd.gc_pct, 1),
                    "Start (0-based)": fwd.start_0based,
                },
                {
                    "Type": "REV",
                    "Spans junction?": "YES" if rev.exon_name == "Junction" else "NO",
                    "Seq (5'->3')": rev.seq_5to3,
                    "Len": rev.length,
                    "Tm (C)": round(rev.tm_c, 2),
                    "GC %": round(rev.gc_pct, 1),
                    "Bind start (0-based)": rev.start_0based,
                },
            ]

            st.write(f"Amplicon length (approx): **{amp_len} bp**")
            st.table(out_rows)

            st.markdown("### Dimer check (heuristic)")
            primer_engine.print_dimer_report_pair(fwd.seq_5to3, rev.seq_5to3)

            if chemistry == "TAQMAN":
                if probe is None:
                    st.warning("No probe returned (probe placement failed). Try increasing amplicon max, or relax probe rules in primer_engine.")
                else:
                    st.markdown("### TaqMan probe")
                    st.table([{
                        "Type": "PROBE",
                        "Seq (5'->3')": probe.seq_5to3,
                        "Len": probe.length,
                        "Tm (C)": round(probe.tm_c, 2),
                        "GC %": round(probe.gc_pct, 1),
                        "Start (0-based)": probe.start_0based,
                    }])

        except Exception as e:
            st.error(str(e))
            st.exception(e)
