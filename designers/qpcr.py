import streamlit as st
import pandas as pd

from utils import clean_dna, gc_pct, tm_wallace, revcomp, primer_score
from utils.blast import primer_blast_url_pair
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION


import streamlit as st

def render():
    st.title("qPCR primers")
    st.write("Design qPCR primers and check specificity with Primer-BLAST.")
    # UI + logic here

