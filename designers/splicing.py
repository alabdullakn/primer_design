treamlit as st
import pandas as pd
from pathlib import Path

from primer_engine import design_exon_primers, print_dimer_report
from primer_engine import design_exon_primers, print_dimer_report, estimate_tm
from utils.blast import primer_blast_url_pair, primer_blast_url_single
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer

ASSETS_DIR = Path(__file__).resolve().parents[1] / "assets"
SPLICING_IMG = ASSETS_DIR / "splicing_examples.png"


def render():
    st.title("Splicing primers")
    st.write("Design primers for exon skipping, intron retention, or alternative splicing.")

    st.subheader("Examples")
    if SPLICING_IMG.exists():
        st.image(str(SPLICING_IMG), use_container_width=True)
    else:
        st.info("Image not found: assets/splicing_examples.png")

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

        with c2:
            tm_target = st.number_input(
                "Target Tm (°C)", 45.0, 75.0, 60.0, key="splicing_tm_target"
            )
            tm_tol = st.number_input(
                "Tm tolerance (± °C)", 1.0, 20.0, 5.0, key="splicing_tm_tol"
            )


        st.markdown("**Reaction and dimer constraints**")
        c3, c4 = st.columns(2)
        with c3:
            tm_model = st.selectbox(
                "Tm model",
                ["Wallace (quick)", "Salt-adjusted"],
                index=0,
                key="splicing_tm_model",
            )
            salt_mM = st.number_input("Buffer/salt concentration (mM)", 1.0, 500.0, 50.0, key="splicing_salt_mM")
            primer_nM = st.number_input("Primer concentration (nM)", 10.0, 2000.0, 250.0, key="splicing_primer_nM")

        with c4:
            max_self_comp_k = st.slider("Max allowed self-complementary stretch (bp)", 4, 12, 8, key="splicing_max_self_comp_k")
            max_heterodimer_risk = st.slider("Max allowed heterodimer risk (%)", 0, 100, 60, key="splicing_max_heterodimer_risk")

    st.subheader("Choose which primer pairs to generate (max 2)")

    colA, colB = st.columns(2)
    with colA:
        pair_aa = st.checkbox("FWD A + REV A", value=True, key="splicing_pair_aa")
        pair_ab = st.checkbox("FWD A + REV B", key="splicing_pair_ab")
    with colB:
        pair_ba = st.checkbox("FWD B + REV A", key="splicing_pair_ba")

    selected_pairs = []
    if pair_aa:
        selected_pairs.append("FWD A + REV A")
    if pair_ab:
        selected_pairs.append("FWD A + REV B")
    if pair_ba:
        selected_pairs.append("FWD B + REV A")

    if len(selected_pairs) == 0:
        st.error("Please select at least one primer pair.")
        return

    if len(selected_pairs) > 2:
        st.error("Maximum 2 primer pairs allowed.")
        return

        return
    st.markdown("---")

    run = st.button("Design primers", key="splicing_design_btn")
    if not run:
        add_footer()
        return

    try:
        exon1_safe = exon_fwd_a.strip() if exon_fwd_a.strip() else exon_fwd_b.strip()
        exon2_safe = exon_fwd_b.strip() if exon_fwd_b.strip() else exon1_safe

        res_A = None
        res_B = None

        if needs_rev_a:
            p1A, p2A, p3A = design_exon_primers(
                exon1_safe,
                exon2_safe,
                exon_rev_a.strip(),
                min_len=int(min_len),
                max_len=int(max_len),
                tm_target=float(tm_target),
                tm_tol=float(tm_tol),
                dimer_k=int(dimer_k),
                tm_model="wallace" if tm_model.startswith("Wallace") else "salt_adjusted",
                salt_mM=float(salt_mM),
                primer_nM=float(primer_nM),
                max_self_comp_k=int(max_self_comp_k),
                max_heterodimer_risk=float(max_heterodimer_risk),
            )
            res_A = (p1A, p2A, p3A)

        if needs_rev_b:
            p1B, p2B, p3B = design_exon_primers(
                exon1_safe,
                exon2_safe,
                exon_rev_b.strip(),
                min_len=int(min_len),
                max_len=int(max_len),
                tm_target=float(tm_target),
                tm_tol=float(tm_tol),
                dimer_k=int(dimer_k),
                tm_model="wallace" if tm_model.startswith("Wallace") else "salt_adjusted",
                salt_mM=float(salt_mM),
                primer_nM=float(primer_nM),
                max_self_comp_k=int(max_self_comp_k),
                max_heterodimer_risk=float(max_heterodimer_risk),
            )
            res_B = (p1B, p2B, p3B)

        st.success("Primers designed successfully.")
        st.caption(SCORE_EXPLANATION)

        out_rows = []
        blast_links = []
        org = "Homo sapiens"

        def add_pair(pair_name: str, fwd_obj, rev_obj):
            out_rows.append(
                {
                    "Pair": pair_name,
                    "Type": "FWD",
                    "Primer (5'→3')": fwd_obj.seq_5to3,
                    "Length": fwd_obj.length,
                    "Tm (°C)": round(fwd_obj.tm_c, 1),
                    "Tm (°C)": round(estimate_tm(fwd_obj.seq_5to3, tm_model="wallace" if tm_model.startswith("Wallace") else "salt_adjusted", salt_mM=float(salt_mM), primer_nM=float(primer_nM)), 1),
                    "GC (%)": round(fwd_obj.gc_pct, 1),
                    "Score": round(fwd_obj.score, 2),
                }
            )
            out_rows.append(
                {
                    "Pair": pair_name,
                    "Type": "REV",
                    "Primer (5'→3')": rev_obj.seq_5to3,
                    "Length": rev_obj.length,
                    "Tm (°C)": round(rev_obj.tm_c, 1),
                    "Tm (°C)": round(estimate_tm(rev_obj.seq_5to3, tm_model="wallace" if tm_model.startswith("Wallace") else "salt_adjusted", salt_mM=float(salt_mM), primer_nM=float(primer_nM)), 1),
                    "GC (%)": round(rev_obj.gc_pct, 1),
                    "Score": round(rev_obj.score, 2),
                }
            )
            blast_links.append(
                (pair_name, primer_blast_url_pair(fwd_obj.seq_5to3, rev_obj.seq_5to3, org))
            )

        for p in selected_pairs:
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

        st.subheader("Primer-BLAST links (NCBI)")
