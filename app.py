from ui.styles import apply_styles

apply_styles()

import streamlit as st

st.set_page_config(
    page_title="PrimerQ",
    page_icon="🧬",
)
import streamlit as st
from designers.regular_pcr import render as render_regular
from designers.splicing import render as render_splicing
from designers.home import render as render_home
from designers.about import render as render_about
from ui.footer import add_footer


st.set_page_config(page_title="Qprimer", layout="wide")
apply_styles()

st.markdown(
    """
    <style>
    /* Equal-size top navigation tabs */
    div[role="radiogroup"] {
        display: grid !important;
        grid-template-columns: repeat(5, minmax(0, 1fr));
        gap: 0.55rem;
        width: 100%;
    }


    div[role="radiogroup"] > label {
        width: 100%;
        margin: 0;
        display: flex;
        justify-content: center;
        background: rgba(15, 23, 42, 0.82);
        border: 1px solid rgba(148, 163, 184, 0.35);
        border-radius: 10px;
        padding: 0.3rem 0.35rem;
    }


    div[role="radiogroup"] > label > div:first-child {
        display: none;
    }

    div[role="radiogroup"] > label p {
        color: #e2e8f0;
        font-weight: 700;
        font-size: 1.02rem;
        text-align: center;
        margin: 0;
        white-space: nowrap;
    }

    div[role="radiogroup"] > label:has(input:checked) {
        border-color: #7dd3fc;
        background: linear-gradient(180deg, rgba(37, 99, 235, 0.96), rgba(30, 64, 175, 0.96));
        box-shadow: 0 8px 20px rgba(37, 99, 235, 0.35);
    }
    </style>
    """,
    unsafe_allow_html=True,
)

TAB_OPTIONS = ["Home", "Regular PCR", "Splicing primers", "qPCR primers", "About us"]
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
elif st.session_state["active_tab"] == "qPCR primers":
    try:
        from designers.qpcr import render as render_qpcr

        render_qpcr()
    except Exception as e:
        st.error("qPCR tab crashed during import or rendering.")
        st.exception(e)
else:
    render_about()

add_footer()
