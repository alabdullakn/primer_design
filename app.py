from ui.styles import apply_styles

apply_styles()

import streamlit as st
from ui.styles import apply_styles
from designers.regular_pcr import render as render_regular
from designers.splicing import render as render_splicing
from designers.qpcr import render as render_qpcr

st.set_page_config(page_title="Primer Designer", layout="wide")
apply_styles()

tabs = st.tabs(["Regular PCR", "Splicing primers", "qPCR primers"])
with tabs[0]:
    render_regular()
with tabs[1]:
    render_splicing()
with tabs[2]:
    render_qpcr()
