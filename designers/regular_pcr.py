import streamlit as st
import pandas as pd

from utils import clean_dna, revcomp, gc_pct, tm_wallace, has_bad_runs
from utils.blast import primer_blast_url_pair, primer_blast_url_single
from ui.blocks import BLAST_INSTRUCTIONS, SCORE_EXPLANATION, add_footer

import streamlit as st

def render():
    st.title("Regular PCR")
    st.write("Paste a sequence and design primers (with heterodimer checks).")
    # UI + logic here

