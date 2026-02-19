import streamlit as st
import pandas as pd

from ..utils.primer_utils import clean_dna, gc_pct, tm_wallace, revcomp, primer_score
from ..utils.blast import primer_blast_url_pair
from ..ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ..ui.footer import add_footer


import streamlit as st

def render():
    st.title("Splicing primers")
    st.write("Design primers for exon skipping, intron retention, or alternative splicing.")
    # UI + logic here

