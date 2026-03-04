import streamlit as st
import pandas as pd
import requests

from primer_engine import design_basic_pcr_primers, print_dimer_report_pair
from utils.blast import primer_blast_url_pair
from utils.seq import gc_pct, primer_score, count_stripped
from utils.tm import tm_nn
from utils.qblast import run_qblast, specificity_label
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer
from ui.blocks import render_diagram, download_primers_csv


def _parse_fasta_text(text: str) -> str:
    if not text:
        return ""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    seq = []
    for l in lines:
        if l.startswith(">"):
            continue
        seq.append(l)
    raw = "".join(seq).upper()
    allowed = set("ACGT")
    return "".join([c for c in raw if c in allowed])


def _fetch_ncbi_fasta(accession: str) -> str:
    acc = (accession or "").strip()
    if not acc:
        return ""
    r = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params={"db": "nuccore", "id": acc, "rettype": "fasta", "retmode": "text"},
        timeout=25,
    )
    if r.status_code != 200:
        raise Exception("Failed to fetch accession from NCBI.")
    return r.text


def render():
    st.title("Regular PCR")
    st.write("Design a standard forward and reverse primer pair from a single template sequence.")
    if st.button("How to use", key="reg_help_btn"):
        st.session_state["tutorial_topic"] = "regular_pcr"
        st.session_state["requested_tab"] = "How to Use"
        st.rerun()

    # persistent template storage across reruns
    if "reg_template_seq" not in st.session_state:
        st.session_state["reg_template_seq"] = ""

    # ============================
    # Primer parameters
    # ============================

    with st.expander("Primer design parameters", expanded=False):
        c1, c2 = st.columns(2)

        with c1:
            min_len = st.number_input("Min primer length", 16, 60, 18, key="reg_min_len")
            max_len = st.number_input("Max primer length", 16, 80, 25, key="reg_max_len")
            tm_target = st.number_input("Target Tm (°C)", 45.0, 75.0, 60.0, key="reg_tm_target")
            organism = st.selectbox(
                "Organism (for BLAST)",
                ["Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Drosophila melanogaster"],
                index=0,
                key="reg_organism",
            )

        with c2:
            tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 25.0, 5.0, key="reg_tm_tol")
            amp_min = st.number_input("Amplicon min (bp)", 50, 5000, 100, key="reg_amp_min")
            amp_max = st.number_input("Amplicon max (bp)", 80, 8000, 500, key="reg_amp_max")
            gc_min = st.number_input("GC min (%)", 20.0, 80.0, 35.0, key="reg_gc_min")
            gc_max = st.number_input("GC max (%)", 20.0, 80.0, 65.0, key="reg_gc_max")

        
        st.caption("Smaller values are allowed, but windows will be clamped safely at runtime.")

        # Option A: allow smaller inputs (min 50), then clamp before calling engine
        c3, c4 = st.columns(2)

        with c3:
            start_window_ui = st.slider(
                "Forward search window (bp from 5')",
                min_value=50,
                max_value=20000,
                value=300,
                step=50,
                key="reg_start_window_ui",
            )

        with c4:
            end_window_ui = st.slider(
                "Reverse search window (bp from 3')",
                min_value=50,
                max_value=20000,
                value=1000,
                step=50,
                key="reg_end_window_ui",
            )
        st.markdown("---")
        st.caption("Reaction buffer (used for Tm correction)")
        edit_reaction_buffer = st.checkbox(
            "Edit reaction buffer",
            value=False,
            key="reg_edit_reaction_buffer",
            help="Enable to customize salt, Mg²⁺, and dNTP concentrations.",
        )
        c5, c6, c7 = st.columns(3)
        with c5:
            sodium_mM = st.number_input(
                "Monovalent salt [Na+] mM",
                1.0,
                500.0,
                50.0,
                key="reg_na_mM",
                disabled=not edit_reaction_buffer,
            )
        with c6:
            mg_mM = st.number_input(
                "Mg²⁺ concentration mM",
                0.0,
                20.0,
                1.5,
                key="reg_mg_mM",
                disabled=not edit_reaction_buffer,
            )
        with c7:
            dntp_mM = st.number_input(
                "Total dNTP mM",
                0.0,
                10.0,
                0.2,
                key="reg_dntp_mM",
                disabled=not edit_reaction_buffer,
            )


    

    # ============================
    # Sequence input
    # ============================

    st.subheader("Sequence input")

    mode = st.radio(
        "Choose input type",
        ["Paste sequence", "Upload FASTA", "NCBI accession"],
        horizontal=True,
        key="reg_input_mode",
    )

    if mode == "Paste sequence":
        text = st.text_area(
            "Paste template sequence (A/C/G/T only)",
            height=220,
            key="reg_seq_paste",
        )
        seq = _parse_fasta_text(">x\n" + (text or ""))
        if seq:
            st.session_state["reg_template_seq"] = seq
            n_stripped = count_stripped(text or "")
            if n_stripped:
                st.warning(f"{n_stripped} non-ACGT character(s) were ignored (ambiguity codes, spaces, etc.).")
            st.success(f"Loaded sequence length: {len(seq)} bp")

    elif mode == "Upload FASTA":
        file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"], key="reg_fasta")
        if file:
            content = file.read().decode("utf-8", errors="ignore")
            seq = _parse_fasta_text(content)
            n_stripped = count_stripped(content)
            if n_stripped:
                st.warning(f"{n_stripped} non-ACGT character(s) were ignored.")
            st.session_state["reg_template_seq"] = seq
            st.success(f"Loaded sequence length: {len(seq)} bp")

    elif mode == "NCBI accession":
        acc = st.text_input("Enter NCBI accession", key="reg_ncbi")
        if st.button("Fetch from NCBI", key="reg_fetch"):
            try:
                fasta = _fetch_ncbi_fasta(acc)
                seq = _parse_fasta_text(fasta)
                if not seq:
                    st.error("Fetched FASTA but could not parse any A/C/G/T bases.")
                else:
                    st.session_state["reg_template_seq"] = seq
                    st.success(f"Fetched {acc} | length: {len(seq)}")
            except Exception as e:
                st.error(str(e))

    template = st.session_state.get("reg_template_seq", "")

    if not template:
        st.info("Paste a sequence, upload a FASTA file, or fetch an NCBI accession.")
        add_footer()
        return


    if int(max_len) < int(min_len):
        st.error("Max primer length must be ≥ min primer length.")
        add_footer()
        return

    run = st.button("Design primers", key="reg_run")

    if run:
        # ============================
        # Option A clamp logic
        # ============================
        n = len(template)
        min_window_needed = max(50, int(min_len) * 2)
        start_window = max(min_window_needed, min(int(start_window_ui), n))
        end_window = max(min_window_needed, min(int(end_window_ui), n))

        if int(amp_min) > int(amp_max):
            st.error("Amplicon min cannot be greater than amplicon max.")
        else:
            try:
                pairs = design_basic_pcr_primers(
                    template,
                    min_len=int(min_len),
                    max_len=int(max_len),
                    tm_target=float(tm_target),
                    tm_tol=float(tm_tol),
                    start_window=int(start_window),
                    end_window=int(end_window),
                    amplicon_min=int(amp_min),
                    amplicon_max=int(amp_max),
                    sodium_mM=float(sodium_mM),
                    mg_mM=float(mg_mM),
                    dntp_mM=float(dntp_mM),
                    gc_min=float(gc_min),
                    gc_max=float(gc_max),
                    top_n=5,
                )
                st.session_state["reg_last_result"] = {
                    "pairs": pairs,
                    "template_len": n,
                    "sodium_mM": float(sodium_mM), "mg_mM": float(mg_mM),
                    "dntp_mM": float(dntp_mM), "tm_target": float(tm_target),
                    "organism": organism,
                }
                st.session_state.pop("reg_blast_results", None)
            except Exception as e:
                st.error(str(e))
                st.caption(f"Debug: template length={n}, start_window={start_window}, end_window={end_window}")

    # ============================
    # Show results (persists across reruns)
    # ============================
    result = st.session_state.get("reg_last_result")
    if result:
        pairs = result["pairs"]
        template_len = result["template_len"]
        _na = result["sodium_mM"]
        _mg = result["mg_mM"]
        _dntp = result["dntp_mM"]
        _tm_target = result["tm_target"]
        _organism = result["organism"]

        st.success(f"Found {len(pairs)} primer pair(s) — showing top {len(pairs)}.")
        st.caption(SCORE_EXPLANATION)

        # Pair selector
        if len(pairs) > 1:
            pair_labels = [f"Pair {i + 1}" + (" — best" if i == 0 else "") for i in range(len(pairs))]
            pair_label = st.radio("Select primer pair", pair_labels, horizontal=True, key="reg_pair_sel")
            sel = int(pair_label.split()[1]) - 1
        else:
            sel = 0

        fwd, rev, amp_len, fwd_start, rev_bind_start = pairs[sel]

        rows = [
            {
                "Type": "FWD",
                "Primer (5'→3')": fwd,
                "Length": len(fwd),
                "Tm (°C)": round(tm_nn(fwd, mv_conc=_na, dv_conc=_mg, dntp_conc=_dntp), 1),
                "GC (%)": round(gc_pct(fwd), 1),
                "Score": round(primer_score(fwd, _tm_target), 2),
                "Start (0-based)": fwd_start,
            },
            {
                "Type": "REV",
                "Primer (5'→3')": rev,
                "Length": len(rev),
                "Tm (°C)": round(tm_nn(rev, mv_conc=_na, dv_conc=_mg, dntp_conc=_dntp), 1),
                "GC (%)": round(gc_pct(rev), 1),
                "Score": round(primer_score(rev, _tm_target), 2),
                "Bind start (0-based)": rev_bind_start,
            },
        ]
        st.dataframe(pd.DataFrame(rows), use_container_width=True)
        download_primers_csv(rows, filename="primers.csv")
        st.write(f"Estimated amplicon length: **{amp_len} bp**")

        st.subheader("Primer positions")
        render_diagram(template_len, fwd_start, len(fwd), rev_bind_start, len(rev), amp_len)

        st.subheader("Primer-BLAST link (NCBI)")
        url = primer_blast_url_pair(fwd, rev, _organism)
        st.markdown(f"[Open in Primer-BLAST]({url})")
        st.info(BLAST_INSTRUCTIONS)

        st.subheader("Dimer check")
        print_dimer_report_pair(fwd, rev)

        # ============================
        # In-app specificity (QBLAST)
        # ============================
        st.subheader("In-app Specificity Check (QBLAST)")
        st.caption("Queries NCBI BLAST directly — takes ~30 s per primer. Results persist until you redesign.")
        if st.button("Check Specificity with BLAST", key="reg_qblast_run"):
            blast_results = {}
            with st.spinner("BLASTing forward primer…"):
                blast_results["FWD"] = run_qblast(fwd, organism=_organism)
            with st.spinner("BLASTing reverse primer…"):
                blast_results["REV"] = run_qblast(rev, organism=_organism)
            st.session_state["reg_blast_results"] = blast_results

        blast_results = st.session_state.get("reg_blast_results")
        if blast_results:
            for label, hits in blast_results.items():
                st.markdown(f"**{label} primer** — {specificity_label(hits)}")
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
