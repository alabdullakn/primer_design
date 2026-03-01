import streamlit as st

# ---------------- PAGE CONFIG (MUST BE FIRST STREAMLIT CALL) ----------------
st.set_page_config(
    page_title="PrimerQ | qPCR Primer Design and Splicing Primer Tool",
    page_icon="🧬",
    layout="wide",
)

# ---------------- STYLES ----------------
from ui.styles import apply_styles
apply_styles()

# ---------------- SEO TEXT (GOOGLE READS THIS) ----------------
st.markdown("""
# PrimerQ

PrimerQ is a free online primer design platform for molecular biology research.

Design primers for:
- qPCR primer design
- exon junction primers
- alternative splicing detection
- SYBR Green qPCR primers
- TaqMan probe primers
- standard PCR primers

Built for gene expression analysis, RNA splicing studies, and cancer research workflows.
""")

# Hidden keywords to help search engines understand the site topic
st.markdown("""
<div style="display:none;">
qPCR primer design tool
splicing primer designer
exon junction primer design
PCR primer design online
SYBR Green primer design
TaqMan probe design tool
gene expression primer design software
primer design web app
</div>
""", unsafe_allow_html=True)

# ---------------- IMPORT APP MODULES ----------------
from designers.regular_pcr import render as render_regular
from designers.splicing import render as render_splicing
from designers.home import render as render_home
from designers.about import render as render_about
from ui.footer import add_footer

# ---------------- YOUR EXISTING CSS ----------------
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

# ---------------- NAVIGATION ----------------
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

# ---------------- ROUTING ----------------
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

# ---------------- FOOTER ----------------
add_footer()
