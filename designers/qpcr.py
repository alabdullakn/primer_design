import streamlit as st
import pandas as pd

from utils import clean_dna, revcomp, gc_pct, tm_wallace, has_bad_runs
from utils.blast import primer_blast_url_pair, primer_blast_url_single
from ui.blocks import BLAST_INSTRUCTIONS, SCORE_EXPLANATION, add_footer


import streamlit as st

def render():
    st.title("qPCR primers")
    st.write("Design qPCR primers and check specificity with Primer-BLAST.")
    # UI + logic here

