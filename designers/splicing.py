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



import streamlit as st

def render():
    st.title("Splicing primers")
    st.write("Design primers for exon skipping, intron retention, or alternative splicing.")
    st.subheader("Examples")
    st.image(str(SPLICING_IMG), use_container_width=True)
    # UI + logic here

