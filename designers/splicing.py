import streamlit as st
import pandas as pd

from utils.primer_utils import clean_dna, gc_pct, tm_wallace, revcomp, primer_score
from utils.blast import primer_blast_url_pair, primer_blast_url_single

from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer
from pathlib import Path
import streamlit as st

ASSETS_DIR = Path(__file__).resolve().parents[1] / "assets"
SPLICING_IMG = ASSETS_DIR / "splicing_examples.png"


def render():
    st.title("Splicing primers")
    st.write("Design primers for exon skipping, intron retention, or alternative splicing.")

    st.subheader("Examples")
    st.image(str(SPLICING_IMG), use_container_width=True)

    st.markdown("---")
def render():
    st.title("Splicing primers")
    st.write("Design primers for exon skipping, intron retention, or alternative splicing.")

    st.subheader("Examples")
    st.image(str(SPLICING_IMG), use_container_width=True)

    st.markdown("---")

    # ============================
    # Step 2: Pair selection (max 2)
    # ============================

    st.subheader("Choose which primer pairs to generate (max 2)")

    c1, c2 = st.columns(2)

    with c1:
        pair_aa = st.checkbox("FWD A + REV A", value=True, key="splicing_pair_aa")
        pair_ab = st.checkbox("FWD A + REV B", key="splicing_pair_ab")

    with c2:
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
        st.stop()

    if len(selected_pairs) > 2:
        st.error("Maximum 2 primer pairs allowed.")
        st.stop()

    # ============================
    # Determine required sequences
    # ============================

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
        exon_fwd_a = st.text_area(
            "Forward A sequence",
            height=140,
            key="splicing_fwd_a",
        )

    if needs_fwd_b:
        exon_fwd_b = st.text_area(
            "Forward B sequence",
            height=140,
            key="splicing_fwd_b",
        )

    if needs_rev_a:
        exon_rev_a = st.text_area(
            "Reverse A sequence",
            height=140,
            key="splicing_rev_a",
        )

    if needs_rev_b:
        exon_rev_b = st.text_area(
            "Reverse B sequence",
            height=140,
            key="splicing_rev_b",
        )

    # ============================
    # Validation
    # ============================

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
        st.stop()

    st.markdown("---")
    st.button("Design primers", key="splicing_design_btn")




    try:
        # If user only selected pairs that need one forward, we still need to pass
        # something for exon1/exon2 into the engine. Use a safe fallback.
        exon1_safe = exon_fwd_a.strip() if exon_fwd_a.strip() else exon_fwd_b.strip()
        exon2_safe = exon_fwd_b.strip() if exon_fwd_b.strip() else exon1_safe

        res_A = None  # results when using Reverse A sequence
        res_B = None  # results when using Reverse B sequence

        # Run engine for Reverse A set if any selected pair ends with REV A
        if needs_rev_a:
            p1A, p2A, p3A = design_exon_primers(
                exon1_safe,
                exon2_safe,
                exon_rev_a.strip(),
                min_len=min_len,
                max_len=max_len,
                tm_target=tm_target,
                tm_tol=tm_tol,
                dimer_k=dimer_k,
            )
            res_A = (p1A, p2A, p3A)

        # Run engine for Reverse B set if any selected pair ends with REV B
        if needs_rev_b:
            p1B, p2B, p3B = design_exon_primers(
                exon1_safe,
                exon2_safe,
                exon_rev_b.strip(),
                min_len=min_len,
                max_len=max_len,
                tm_target=tm_target,
                tm_tol=tm_tol,
                dimer_k=dimer_k,
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
                    "GC (%)": round(rev_obj.gc_pct, 1),
                    "Score": round(rev_obj.score, 2),
                }
            )
            blast_links.append(
                (pair_name, primer_blast_url_pair(fwd_obj.seq_5to3, rev_obj.seq_5to3, org))
            )

        # Build outputs for only the selected pairs
        for p in selected_pairs:
            if p == "FWD A + REV A":
                if res_A is None:
                    continue
                p1A, p2A, p3A = res_A
                add_pair("FWD A + REV A", p1A, p3A)

            elif p == "FWD B + REV A":
                if res_A is None:
                    continue
                p1A, p2A, p3A = res_A
                add_pair("FWD B + REV A", p2A, p3A)

            elif p == "FWD A + REV B":
                if res_B is None:
                    continue
                p1B, p2B, p3B = res_B
                add_pair("FWD A + REV B", p1B, p3B)

        if out_rows:
            st.dataframe(pd.DataFrame(out_rows), use_container_width=True)

        st.subheader("Primer-BLAST links (NCBI)")
        for name, url in blast_links:
            st.markdown(f"**{name}**: [Open in Primer-BLAST]({url})")
        st.info(BLAST_INSTRUCTIONS)

        # Dimer reports for each reverse set (only if generated)
        if res_A is not None:
            st.subheader("Dimer check (Reverse A set)")
            p1A, p2A, p3A = res_A
            print_dimer_report(p1A, p2A, p3A)

        if res_B is not None:
            st.subheader("Dimer check (Reverse B set)")
            p1B, p2B, p3B = res_B
            print_dimer_report(p1B, p2B, p3B)

        add_footer()

    except Exception as e:
        st.error(str(e))
