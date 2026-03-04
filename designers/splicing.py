import streamlit as st
import pandas as pd
import base64
from pathlib import Path

from primer_engine import design_exon_primers, print_dimer_report
from utils.blast import primer_blast_url_pair, primer_blast_url_single
from utils.qblast import run_qblast, specificity_label
from utils.seq import count_stripped

from ui.blocks import render_diagram, download_primers_csv
from utils.examples import EXAMPLE_SPLICING_FWD_A, EXAMPLE_SPLICING_FWD_B, EXAMPLE_SPLICING_REV_A

ASSETS_DIR = Path(__file__).resolve().parents[1] / "assets"
SPLICING_IMG = ASSETS_DIR / "splicing_examples.png"


def render():
    st.title("Splicing primers")
    st.write("Design primers for exon skipping, intron retention, or alternative splicing.")
    if st.button("How to use", key="splicing_help_btn"):
        st.session_state["tutorial_topic"] = "splicing"
        st.session_state["requested_tab"] = "How to Use"
        st.rerun()

    st.subheader("Examples")
    if SPLICING_IMG.exists():
        image_b64 = base64.b64encode(SPLICING_IMG.read_bytes()).decode("utf-8")
        st.markdown(
            f"""
            <div style=\"max-width:700px; margin:0 auto; padding:10px; border:1px solid rgba(255, 255, 255, 0.35); border-radius:10px; background:rgba(255, 255, 255, 0.08);\">
                <img src=\"data:image/png;base64,{image_b64}\" style=\"width:100%; display:block; border-radius:6px;\" />
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.info("Image not found: assets/splicing_examples.png")
        
    st.markdown('<div style="height: 0.75rem;"></div>', unsafe_allow_html=True)

    with st.expander("Primer design parameters", expanded=False):
        c1, c2 = st.columns(2)

        with c1:
            min_len = st.number_input(
                "Min primer length", 16, 40, 18, key="splicing_min_len"
            )
            max_len = st.number_input(
                "Max primer length", 16, 60, 25, key="splicing_max_len"
            )
            dimer_k = st.number_input(
                "3' dimer check window (k)", 3, 10, 4, key="splicing_dimer_k"
            )
            organism = st.selectbox(
                "Organism (for BLAST)",
                ["Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Drosophila melanogaster"],
                index=0,
                key="splicing_organism",
            )

        with c2:
            tm_target = st.number_input(
                "Target Tm (°C)", 45.0, 75.0, 60.0, key="splicing_tm_target"
            )
            tm_tol = st.number_input(
                "Tm tolerance (± °C)", 1.0, 20.0, 5.0, key="splicing_tm_tol"
            )
            gc_min = st.number_input("GC min (%)", 20.0, 80.0, 35.0, key="splicing_gc_min")
            gc_max = st.number_input("GC max (%)", 20.0, 80.0, 65.0, key="splicing_gc_max")

        st.markdown("---")
        st.caption("Reaction buffer (used for Tm correction)")
        edit_reaction_buffer = st.checkbox(
            "Edit reaction buffer",
            value=False,
            key="splicing_edit_reaction_buffer",
            help="Enable to customize salt, Mg²⁺, and dNTP concentrations.",
        )
        c3, c4, c5 = st.columns(3)
        with c3:
            sodium_mM = st.number_input(
                "Monovalent salt [Na+] mM",
                1.0,
                500.0,
                50.0,
                key="splicing_na_mM",
                disabled=not edit_reaction_buffer,
            )
        with c4:
            mg_mM = st.number_input(
                "Mg²⁺ concentration mM",
                0.0,
                20.0,
                1.5,
                key="splicing_mg_mM",
                disabled=not edit_reaction_buffer,
            )
        with c5:
            dntp_mM = st.number_input(
                "Total dNTP mM",
                0.0,
                10.0,
                0.2,
                key="splicing_dntp_mM",
                disabled=not edit_reaction_buffer,
            )
    st.subheader("Choose which condition to generate")
    st.caption("Look at the image above for reference.")
    
    st.markdown(
        """
        <style>
        .st-key-splicing_condition [role="radiogroup"] {
            gap: 0.75rem;
            flex-wrap: wrap;
        }

        .st-key-splicing_condition [role="radiogroup"] > label {
            flex: 1 1 320px;
            min-width: 0;
            min-height: 3.4rem;
            align-items: center;
            justify-content: center;
            padding: 0.4rem 0.7rem;
        }
        .st-key-splicing_condition [role="radiogroup"] > label p {
            white-space: normal;
            overflow-wrap: anywhere;
            text-align: center;
            line-height: 1.15;
            font-size: clamp(0.78rem, 0.45vw + 0.55rem, 1.02rem);
        }

        @media (max-width: 640px) {
            .st-key-splicing_condition [role="radiogroup"] > label {
                flex-basis: 100%;
            }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

    condition = st.radio(
        "Choose condition",
        options=[
            "Condition 1: FWD A + FWD B + REV A",
            "Condition 2: FWD A + REV A + REV B",
        ],
        index=0,
        key="splicing_condition",
        label_visibility="collapsed",
    )

    if condition.startswith("Condition 1"):
        selected_pairs = ["FWD A + REV A", "FWD B + REV A"]
    else:
        selected_pairs = ["FWD A + REV A", "FWD A + REV B"]

    if st.button("Load example (human CD44 splicing)", key="splicing_load_example"):
        st.session_state["splicing_condition"] = "Condition 1: FWD A + FWD B + REV A"
        st.session_state["splicing_fwd_a"] = EXAMPLE_SPLICING_FWD_A
        st.session_state["splicing_fwd_b"] = EXAMPLE_SPLICING_FWD_B
        st.session_state["splicing_rev_a"] = EXAMPLE_SPLICING_REV_A
        st.rerun()

    needs_fwd_a = any(p.startswith("FWD A") for p in selected_pairs)
    needs_fwd_b = any(p.startswith("FWD B") for p in selected_pairs)
    needs_rev_a = any(p.endswith("REV A") for p in selected_pairs)
    needs_rev_b = any(p.endswith("REV B") for p in selected_pairs)

    st.subheader("Paste sequences (A/C/G/T only)")

    exon_fwd_a = ""
    exon_fwd_b = ""
    exon_rev_a = ""
    exon_rev_b = ""

    if needs_fwd_a:
        exon_fwd_a = st.text_area("Forward A sequence", height=140, key="splicing_fwd_a")
        if exon_fwd_a and count_stripped(exon_fwd_a):
            st.warning(f"Forward A: {count_stripped(exon_fwd_a)} non-ACGT character(s) ignored.")
    if needs_fwd_b:
        exon_fwd_b = st.text_area("Forward B sequence", height=140, key="splicing_fwd_b")
        if exon_fwd_b and count_stripped(exon_fwd_b):
            st.warning(f"Forward B: {count_stripped(exon_fwd_b)} non-ACGT character(s) ignored.")
    if needs_rev_a:
        exon_rev_a = st.text_area("Reverse A sequence", height=140, key="splicing_rev_a")
        if exon_rev_a and count_stripped(exon_rev_a):
            st.warning(f"Reverse A: {count_stripped(exon_rev_a)} non-ACGT character(s) ignored.")
    if needs_rev_b:
        exon_rev_b = st.text_area("Reverse B sequence", height=140, key="splicing_rev_b")
        if exon_rev_b and count_stripped(exon_rev_b):
            st.warning(f"Reverse B: {count_stripped(exon_rev_b)} non-ACGT character(s) ignored.")

    missing = []
    if needs_fwd_a and not exon_fwd_a.strip():
        missing.append("Forward A")
    if needs_fwd_b and not exon_fwd_b.strip():
        missing.append("Forward B")
    if needs_rev_a and not exon_rev_a.strip():
        missing.append("Reverse A")
    if needs_rev_b and not exon_rev_b.strip():
        missing.append("Reverse B")

    if missing:
        st.error("Missing: " + ", ".join(missing))
        return
    st.markdown("---")

    if int(max_len) < int(min_len):
        st.error("Max primer length must be ≥ min primer length.")
        return

    run = st.button("Design primers", key="splicing_design_btn")

    if run:
        try:
            exon1_safe = exon_fwd_a.strip() if exon_fwd_a.strip() else exon_fwd_b.strip()
            exon2_safe = exon_fwd_b.strip() if exon_fwd_b.strip() else exon1_safe

            res_A = None
            res_B = None

            if needs_rev_a:
                p1A, p2A, p3A = design_exon_primers(
                    exon1_safe, exon2_safe, exon_rev_a.strip(),
                    min_len=int(min_len), max_len=int(max_len),
                    tm_target=float(tm_target), tm_tol=float(tm_tol),
                    dimer_k=int(dimer_k), sodium_mM=float(sodium_mM),
                    mg_mM=float(mg_mM), dntp_mM=float(dntp_mM),
                    gc_min=float(gc_min), gc_max=float(gc_max),
                )
                res_A = (p1A, p2A, p3A)

            if needs_rev_b:
                p1B, p2B, p3B = design_exon_primers(
                    exon1_safe, exon2_safe, exon_rev_b.strip(),
                    min_len=int(min_len), max_len=int(max_len),
                    tm_target=float(tm_target), tm_tol=float(tm_tol),
                    dimer_k=int(dimer_k), sodium_mM=float(sodium_mM),
                    mg_mM=float(mg_mM), dntp_mM=float(dntp_mM),
                    gc_min=float(gc_min), gc_max=float(gc_max),
                )
                res_B = (p1B, p2B, p3B)

            from utils.seq import clean_seq as _cs
            st.session_state["splicing_last_result"] = {
                "res_A": res_A, "res_B": res_B,
                "selected_pairs": selected_pairs, "organism": organism,
                "len_exon1": len(_cs(exon1_safe)),
                "len_exon_rev_a": len(_cs(exon_rev_a.strip())) if exon_rev_a.strip() else 0,
                "len_exon_rev_b": len(_cs(exon_rev_b.strip())) if exon_rev_b.strip() else 0,
            }
            st.session_state.pop("splicing_blast_results", None)
        except Exception as e:
            st.error(str(e))

    # ============================
    # Show results (persists across reruns)
    # ============================
    result = st.session_state.get("splicing_last_result")
    if result:
        res_A = result["res_A"]
        res_B = result["res_B"]
        _selected_pairs = result["selected_pairs"]
        _organism = result["organism"]

        st.success("Primers designed successfully.")

        out_rows = []
        blast_links = []

        def add_pair(pair_name: str, fwd_obj, rev_obj):
            out_rows.append({
                "Pair": pair_name, "Type": "FWD",
                "Primer (5'→3')": fwd_obj.seq_5to3, "Length": fwd_obj.length,
                "Tm (°C)": round(fwd_obj.tm_c, 1), "GC (%)": round(fwd_obj.gc_pct, 1),
                "Score": round(fwd_obj.score, 2),
            })
            out_rows.append({
                "Pair": pair_name, "Type": "REV",
                "Primer (5'→3')": rev_obj.seq_5to3, "Length": rev_obj.length,
                "Tm (°C)": round(rev_obj.tm_c, 1), "GC (%)": round(rev_obj.gc_pct, 1),
                "Score": round(rev_obj.score, 2),
            })
            blast_links.append(
                (pair_name, primer_blast_url_pair(fwd_obj.seq_5to3, rev_obj.seq_5to3, _organism))
            )

        for p in _selected_pairs:
            if p == "FWD A + REV A" and res_A is not None:
                p1A, p2A, p3A = res_A
                add_pair("FWD A + REV A", p1A, p3A)
            if p == "FWD B + REV A" and res_A is not None:
                p1A, p2A, p3A = res_A
                add_pair("FWD B + REV A", p2A, p3A)
            if p == "FWD A + REV B" and res_B is not None:
                p1B, p2B, p3B = res_B
                add_pair("FWD A + REV B", p1B, p3B)

        if out_rows:
            st.dataframe(pd.DataFrame(out_rows), use_container_width=True)
            download_primers_csv(out_rows, filename="splicing_primers.csv")

        if res_A is not None:
            p1A, p2A, p3A = res_A
            _len_a = result.get("len_exon1", 0) + result.get("len_exon_rev_a", 0)
            if _len_a > 0:
                st.subheader("Primer positions — Reverse A set")
                st.caption("FWD A (Exon 1)")
                render_diagram(_len_a, p1A.start_0based, p1A.length, p3A.start_0based, p3A.length,
                               max(0, (p3A.start_0based + p3A.length) - p1A.start_0based))
                st.caption("FWD B (Exon 2)")
                render_diagram(_len_a, p2A.start_0based, p2A.length, p3A.start_0based, p3A.length,
                               max(0, (p3A.start_0based + p3A.length) - p2A.start_0based))
            st.subheader("Dimer check (Reverse A set)")
            print_dimer_report(p1A, p2A, p3A)
        if res_B is not None:
            p1B, p2B, p3B = res_B
            _len_b = result.get("len_exon1", 0) + result.get("len_exon_rev_b", 0)
            if _len_b > 0:
                st.subheader("Primer positions — Reverse B set")
                render_diagram(_len_b, p1B.start_0based, p1B.length, p3B.start_0based, p3B.length,
                               max(0, (p3B.start_0based + p3B.length) - p1B.start_0based))
            st.subheader("Dimer check (Reverse B set)")
            print_dimer_report(p1B, p2B, p3B)

        # ============================
        # In-app specificity (QBLAST)
        # ============================
        st.subheader("Specificity Check (QBLAST)")

        blast_seqs = {}
        if res_A is not None:
            p1A, p2A, p3A = res_A
            blast_seqs["FWD A"] = p1A.seq_5to3
            blast_seqs["FWD B"] = p2A.seq_5to3
            blast_seqs["REV A"] = p3A.seq_5to3
        if res_B is not None:
            p1B, p2B, p3B = res_B
            blast_seqs["REV B"] = p3B.seq_5to3

        if st.button("Check Specificity with BLAST", key="splicing_qblast_run"):
            blast_results = {}
            try:
                for label, pseq in blast_seqs.items():
                    with st.spinner(f"BLASTing {label} (~30 s)…"):
                        blast_results[label] = run_qblast(pseq, organism=_organism)
                st.session_state["splicing_blast_results"] = blast_results
            except RuntimeError as e:
                st.error(str(e))

        blast_results = st.session_state.get("splicing_blast_results")
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

        st.subheader("Primer-BLAST (NCBI)")
        for name, url in blast_links:
            st.markdown(f"**{name}**: [Open in Primer-BLAST ↗]({url})")
