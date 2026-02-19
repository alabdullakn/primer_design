import streamlit as st
import pandas as pd

from utils import clean_dna, gc_pct, tm_wallace, revcomp, primer_score
from utils.checks import primer_blast_url_pair, primer_blast_url_single
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer


import streamlit as st

def render():
    st.title("Regular PCR")
    st.write("Paste a sequence and design primers (with heterodimer checks).")
    # UI + logic here

