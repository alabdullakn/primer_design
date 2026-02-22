from ui.styles import apply_styles

apply_styles()

import streamlit as st
from designers.regular_pcr import render as render_regular
from designers.splicing import render as render_splicing
from designers.home import render as render_home

st.set_page_config(page_title="Qprimer", layout="wide")
apply_styles()

tabs = st.tabs(["Home", "Regular PCR", "Splicing primers", "qPCR primers"])

with tabs[0]:
    render_home()

with tabs[1]:
    render_regular()

with tabs[2]:
    render_splicing()

with tabs[3]:
    try:
        from designers.qpcr import render as render_qpcr
        render_qpcr()
    except Exception as e:
        st.error("qPCR tab crashed during import or rendering.")
        st.exception(e)
designers/home.py
