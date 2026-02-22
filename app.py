from ui.styles import apply_styles

apply_styles()

import streamlit as st
from designers.regular_pcr import render as render_regular
from designers.splicing import render as render_splicing
from designers.home import render as render_home

st.set_page_config(page_title="Qprimer", layout="wide")
apply_styles()

TAB_OPTIONS = ["Home", "Regular PCR", "Splicing primers", "qPCR primers"]
if "active_tab" not in st.session_state:
    st.session_state["active_tab"] = "Home"

if "requested_tab" in st.session_state:
    st.session_state["active_tab"] = st.session_state.pop("requested_tab")

st.radio(
    "Navigation",
    TAB_OPTIONS,
    key="active_tab",
    horizontal=True,
    label_visibility="collapsed",
)

if st.session_state["active_tab"] == "Home":
    render_home()
elif st.session_state["active_tab"] == "Regular PCR":
    render_regular()
elif st.session_state["active_tab"] == "Splicing primers":
    render_splicing()
else:
    try:
        from designers.qpcr import render as render_qpcr

        render_qpcr()
    except Exception as e:
        st.error("qPCR tab crashed during import or rendering.")
        st.exception(e)
