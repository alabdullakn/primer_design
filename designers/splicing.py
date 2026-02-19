import streamlit as st
import pandas as pd

from utils import clean_dna, revcomp, gc_pct, tm_wallace, has_bad_runs
from utils.blast import primer_blast_url_pair, primer_blast_url_single
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION

import streamlit as st

def render():
    st.title("Splicing primers")
    st.write("Design primers for exon skipping, intron retention, or alternative splicing.")
    # UI + logic here

