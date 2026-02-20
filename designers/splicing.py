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
