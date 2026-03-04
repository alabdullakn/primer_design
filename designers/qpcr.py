# designers/qpcr.py
import streamlit as st
import pandas as pd

from qpcr_engine import (
    amplicon_size_from_hits,
    design_qpcr_junction_pair,
    design_intron_flanking_pair,
)
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from utils.blast import primer_blast_url_pair, primer_blast_url_single
from utils.qblast import run_qblast, specificity_label
from ui.footer import add_footer
from ui.blocks import render_diagram, download_primers_csv

def render():
    col_title, col_btn = st.columns([5, 1])
    with col_title:
        st.title("qPCR primers")
    with col_btn:
        st.write("")  # vertical alignment spacer
        if st.button("How to use", key="qpcr_help_btn"):
            st.session_state["tutorial_topic"] = "qpcr"
            st.session_state["requested_tab"] = "How to Use"
            st.rerun()
    
    # ============================
    # Strategy selector
    # ============================
    strategy = st.radio(
        "qPCR strategy",
        ["Junction-spanning", "Intron-flanking"],
        horizontal=True,
        key="qpcr_strategy",
        help=(
            "Junction-spanning: one primer straddles the exon-exon junction — "
            "mark it with ^ in your sequence.\n\n"
            "Intron-flanking: primers sit in two different exons; "
            "the intron between them makes the genomic amplicon too large to amplify."
        ),
    )
    is_junction = strategy == "Junction-spanning"

    # ============================
    # Shared parameters expander
    # ============================
    with st.expander("Primer parameters", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            chemistry_label = st.selectbox(
                "Chemistry",
                ["SYBR Green", "TaqMan probe"] if is_junction else ["SYBR Green"],
                index=0,
                help=(
                    "SYBR Green: primer pair only. "
                    "TaqMan: adds an internal probe (junction-spanning mode only)."
                ),
                key="qpcr_chemistry",
            )
            organism = st.selectbox(
                "Organism (for BLAST)",
                ["Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Drosophila melanogaster"],
                index=0,
                key="qpcr_organism",
            )
        with c2:
            if is_junction:
                junction_primer = st.selectbox(
                    "Junction spanning primer",
                    ["AUTO", "FWD", "REV"],
                    index=0,
                    help="Force FWD or REV to span the junction, or let the algorithm choose.",
                    key="qpcr_junction_primer",
                )
            else:
                junction_primer = "AUTO"
        with c3:
            dimer_k = st.number_input("3' dimer check window (k)", 3, 10, 4, key="qpcr_dimer_k")

        st.subheader("Primer constraints")
        r1, r2 = st.columns(2)
        with r1:
            min_len = st.number_input("Min length", 16, 28, 18, key="qpcr_min_len")
            tm_target = st.number_input("Tm target (°C)", 50.0, 70.0, 60.0, key="qpcr_tm_target")
            gc_min = st.number_input("GC min percent", 20.0, 80.0, 40.0, key="qpcr_gc_min")
            amp_min = st.number_input("Amplicon min bp", 40, 300, 70, key="qpcr_amp_min")
            if is_junction:
                j_ov = st.number_input("Min overlap each side of junction", 3, 12, 6, key="qpcr_j_ov")
            else:
                j_ov = 6

        with r2:
            max_len = st.number_input("Max length", 16, 32, 24, key="qpcr_max_len")
            tm_tol = st.number_input("Tm tolerance (± °C)", 0.5, 10.0, 2.0, key="qpcr_tm_tol")
            gc_max = st.number_input("GC max percent", 20.0, 80.0, 60.0, key="qpcr_gc_max")
            amp_max = st.number_input("Amplicon max bp", 60, 500, 200, key="qpcr_amp_max")
            if is_junction:
                junction_3p_max_distance = st.number_input(
                    "Max distance of junction from 3' end",
                    1, 10, 8,
                    help="Smaller = stricter. Primers with junction at position -1 or -2 are scored best.",
                    key="qpcr_junction_3p_max_distance",
                )
            else:
                junction_3p_max_distance = 8

        c4, c5 = st.columns([1, 1])
        with c4:
            max_hpoly = st.number_input("Max homopolymer run", 2, 6, 3, key="qpcr_max_hpoly")
            max_tm_diff = st.number_input("Max Tm difference pair", 0.0, 10.0, 1.0, key="qpcr_max_tm_diff")

        st.markdown("---")
        st.caption("Reaction buffer (used for Tm and structure calculations)")
        edit_buffer = st.checkbox("Edit reaction buffer", value=False, key="qpcr_edit_buffer",
                                  help="Enable to customize salt, Mg²⁺, and dNTP concentrations.")
        b1, b2, b3 = st.columns(3)
        with b1:
            mv_conc = st.number_input("Monovalent salt [Na⁺] mM", 1.0, 500.0, 50.0,
                                      key="qpcr_na_mM", disabled=not edit_buffer)
        with b2:
            dv_conc = st.number_input("Mg²⁺ concentration mM", 0.0, 20.0, 1.5,
                                      key="qpcr_mg_mM", disabled=not edit_buffer)
        with b3:
            dntp_conc = st.number_input("Total dNTP mM", 0.0, 10.0, 0.2,
                                        key="qpcr_dntp_mM", disabled=not edit_buffer)

        chemistry = "TAQMAN" if chemistry_label == "TaqMan probe" else "SYBR"

        probe_tm_target = 69.0
        probe_tm_tol = 3.0
        probe_min_len = 18
        probe_max_len = 30

        if is_junction and chemistry == "TAQMAN":
            st.subheader("TaqMan probe constraints")
            p1, p2 = st.columns(2)
            with p1:
                probe_tm_target = st.number_input("Probe Tm target", 60.0, 78.0, 69.0, key="qpcr_probe_tm_target")
                probe_min_len = st.number_input("Probe min length", 16, 28, 18, key="qpcr_probe_min_len")
            with p2:
                probe_tm_tol = st.number_input("Probe Tm tolerance", 0.5, 8.0, 3.0, key="qpcr_probe_tm_tol")
                probe_max_len = st.number_input("Probe max length", 18, 35, 30, key="qpcr_probe_max_len")

    # ============================
    # Sequence input (mode-specific)
    # ============================
    if is_junction:
        st.write("Mark the exon-exon junction with `^` in your cDNA/mRNA sequence.")
        st.caption("Example: `...AAGCTAGTCCAG^GTTAGCCTAAGG...`  — provide ≥10 bases on each side.")
        seq = st.text_area(
            "Sequence with junction marker ^",
            height=160,
            placeholder="Paste cDNA sequence here with ^ at the junction",
            key="qpcr_seq_junction",
        )
        exon_a = exon_b = ""

        if is_junction and junction_3p_max_distance < j_ov:
            st.warning("Contradictory settings: max 3' distance must be ≥ min overlap each side.")
    else:
        st.write("Paste the two exon sequences separately. The FWD primer will bind Exon A, REV will bind Exon B.")
        st.caption(
            "Tip: provide the last ~100 bp of Exon A and the first ~100 bp of Exon B. "
            "The cDNA amplicon spans the junction; genomic DNA won't amplify (large intron between them)."
        )
        col_a, col_b = st.columns(2)
        with col_a:
            exon_a = st.text_area("Exon A sequence (FWD primer here)", height=160, key="qpcr_exon_a")
        with col_b:
            exon_b = st.text_area("Exon B sequence (REV primer here)", height=160, key="qpcr_exon_b")
        seq = ""

    if int(max_len) < int(min_len):
        st.error("Max primer length must be ≥ min primer length.")
    if int(amp_max) < int(amp_min):
        st.error("Amplicon max must be ≥ amplicon min.")

    st.divider()

    # ============================
    # Design button
    # ============================
    if st.button("Design qPCR primers"):
        if is_junction:
            if not seq or "^" not in seq:
                st.error("Your sequence must include the junction marker ^")
            else:
                try:
                    pairs = design_qpcr_junction_pair(
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
                        mv_conc=float(mv_conc),
                        dv_conc=float(dv_conc),
                        dntp_conc=float(dntp_conc),
                        top_n=5,
                    )
                    # Derive junction position for diagram
                    from qpcr_engine import parse_junction_marked_seq
                    _left, _ = parse_junction_marked_seq(seq)
                    st.session_state["qpcr_last_result"] = {
                        "mode": "junction",
                        "pairs": pairs,
                        "junction_pos": len(_left),
                        "template_len": len(_left) + len(seq.split("^", 1)[1].replace(" ", "")),
                        "chemistry_label": chemistry_label,
                        "organism": organism,
                    }
                    st.session_state.pop("qpcr_blast_results", None)
                except Exception as e:
                    st.error(f"Design failed: {e}")
                    st.info("Try widening Tm tolerance/GC range, relaxing junction distance, or increasing amplicon max.")
        else:
            if not exon_a.strip() or not exon_b.strip():
                st.error("Provide both Exon A and Exon B sequences.")
            else:
                try:
                    from utils.seq import clean_seq as _cs
                    _ea = _cs(exon_a.strip())
                    pairs = design_intron_flanking_pair(
                        exon_a=exon_a.strip(),
                        exon_b=exon_b.strip(),
                        min_len=int(min_len),
                        max_len=int(max_len),
                        tm_target=float(tm_target),
                        tm_tol=float(tm_tol),
                        gc_min=float(gc_min),
                        gc_max=float(gc_max),
                        max_homopolymer=int(max_hpoly),
                        amplicon_min=int(amp_min),
                        amplicon_max=int(amp_max),
                        dimer_k=int(dimer_k),
                        max_tm_diff_pair=float(max_tm_diff),
                        mv_conc=float(mv_conc),
                        dv_conc=float(dv_conc),
                        dntp_conc=float(dntp_conc),
                        top_n=5,
                    )
                    # Wrap as (fwd, rev, probe=None) triples for uniform handling
                    pairs = [(f, r, None) for f, r in pairs]
                    from utils.seq import clean_seq as _cs2
                    _eb = _cs2(exon_b.strip())
                    st.session_state["qpcr_last_result"] = {
                        "mode": "intron_flanking",
                        "pairs": pairs,
                        "junction_pos": len(_ea),
                        "template_len": len(_ea) + len(_eb),
                        "chemistry_label": "SYBR Green",
                        "organism": organism,
                    }
                    st.session_state.pop("qpcr_blast_results", None)
                except Exception as e:
                    st.error(f"Design failed: {e}")
                    st.info("Try widening Tm tolerance/GC range, amplicon range, or provide longer exon sequences.")

    # ============================
    # Show results (persists across reruns)
    # ============================
    result = st.session_state.get("qpcr_last_result")
    if result:
        pairs = result["pairs"]
        _chemistry_label = result["chemistry_label"]
        _organism = result["organism"]
        _mode = result.get("mode", "junction")
        _template_len = result.get("template_len", 0)
        _junction_pos = result.get("junction_pos")

        st.success(f"Found {len(pairs)} primer pair(s) — showing top {len(pairs)}.")
        st.write(f"Chemistry: {_chemistry_label}")
        st.caption(SCORE_EXPLANATION)

        # Pair selector
        if len(pairs) > 1:
            pair_labels = [f"Pair {i + 1}" + (" — best" if i == 0 else "") for i in range(len(pairs))]
            pair_label = st.radio("Select primer pair", pair_labels, horizontal=True, key="qpcr_pair_sel")
            sel = int(pair_label.split()[1]) - 1
        else:
            sel = 0

        fwd, rev, probe = pairs[sel]
        amp_size = amplicon_size_from_hits(fwd, rev)

        if _mode == "intron_flanking":
            st.info(
                f"cDNA amplicon: **{amp_size} bp**. "
                "The genomic amplicon will be larger by the intron size — "
                "run on a gel to confirm mRNA specificity."
            )

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
            rows.append({
                "Type": "PROBE",
                "Primer (5'→3')": probe.seq_5to3,
                "Length": probe.length,
                "Tm (°C)": round(probe.tm_c, 1),
                "GC (%)": round(probe.gc_pct, 1),
                "Score": round(probe.score, 2),
                "Start (0-based)": probe.start_0based,
                "Role": "Internal probe",
            })

        st.dataframe(pd.DataFrame(rows), use_container_width=True)
        download_primers_csv(rows, filename="qpcr_primers.csv")
        if _mode == "junction":
            st.write(f"Estimated amplicon length: **{amp_size} bp**")

        st.subheader("Primer positions")
        if _template_len > 0:
            render_diagram(
                _template_len,
                fwd.start_0based, fwd.length,
                rev.start_0based, rev.length,
                amp_size,
                junction_pos=_junction_pos,
            )

        st.subheader("Primer-BLAST links (NCBI)")
        pair_url = primer_blast_url_pair(fwd.seq_5to3, rev.seq_5to3, _organism)
        st.markdown(f"**qPCR pair**: [Open in Primer-BLAST]({pair_url})")

        with st.expander("Single-primer Primer-BLAST links"):
            st.markdown(f"**Forward primer**: [Primer-BLAST]({primer_blast_url_single(fwd.seq_5to3, _organism)})")
            st.markdown(f"**Reverse primer**: [Primer-BLAST]({primer_blast_url_single(rev.seq_5to3, _organism)})")
            if probe is not None:
                st.markdown(f"**Probe**: [Primer-BLAST]({primer_blast_url_single(probe.seq_5to3, _organism)})")

        st.info(BLAST_INSTRUCTIONS)

        # ============================
        # In-app specificity (QBLAST)
        # ============================
        st.subheader("In-app Specificity Check (QBLAST)")
        st.caption("Queries NCBI BLAST directly — takes ~30 s per primer. Results persist until you redesign.")

        blast_seqs = {"FWD": fwd.seq_5to3, "REV": rev.seq_5to3}
        if probe is not None:
            blast_seqs["PROBE"] = probe.seq_5to3

        if st.button("Check Specificity with BLAST", key="qpcr_qblast_run"):
            blast_results = {}
            for label, pseq in blast_seqs.items():
                with st.spinner(f"BLASTing {label}…"):
                    blast_results[label] = run_qblast(pseq, organism=_organism)
            st.session_state["qpcr_blast_results"] = blast_results

        blast_results = st.session_state.get("qpcr_blast_results")
        if blast_results:
            for label, hits in blast_results.items():
                st.markdown(f"**{label}** — {specificity_label(hits)}")
                if hits:
                    hit_rows = [
                        {
                            "Title": h.title,
                            "Accession": h.accession,
                            "Identity (%)": h.identity_pct,
                            "Query Coverage (%)": h.query_coverage_pct,
                            "E-value": h.evalue,
                        }
                        for h in hits
                    ]
                    st.dataframe(pd.DataFrame(hit_rows), use_container_width=True)
                else:
                    st.info("No hits returned.")

    add_footer()
