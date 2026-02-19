import streamlit as st

BLAST_INSTRUCTIONS = (
    "How to use Primer-BLAST:\n"
    "1) Click the **Open in Primer-BLAST** link.\n"
    "2) On the NCBI page, do NOT change any settings.\n"
    "3) Just click **Get Primers**.\n"
    "4) Check the top hit matches your intended gene/transcript.\n"
    "This is a quick specificity check."
)

SCORE_EXPLANATION = (
    "Score (lower is better): internal ranking only.\n"
    "It is based on:\n"
    "• |Tm − target Tm|\n"
    "• +5 penalty for long homopolymer runs (e.g. AAAAA)\n"
    "• +3 penalty for extreme GC% (<35% or >65%)\n"
    "This is not a BLAST score and not experimental validation."
)

FOOTER_TEXT = (
    "Free and open to everyone. Built by Khalid Alabdulla. "
    "If you find it useful, please share it. "
    "Feedback or bugs: alabdulla8932@gmail.com"
)

def add_footer():
    st.markdown("---")
    st.caption(FOOTER_TEXT)
