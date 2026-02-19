import streamlit as st
import pandas as pd

from primer_engine import design_exon_primers, print_dimer_report
from utils.checks import primer_blast_url_pair, primer_blast_url_single
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer


import streamlit as st

def render():
    st.title("Splicing primers")
    st.write("Design primers for exon skipping, intron retention, or alternative splicing.")
    # UI + logic here

