import base64
from pathlib import Path

import streamlit as st


def _background_style() -> str:
    image_path = Path("assets/dna_background.png")
    if image_path.exists():
        encoded = base64.b64encode(image_path.read_bytes()).decode("utf-8")
        return f"""
        .stApp {{
            background-image: linear-gradient(rgba(2, 6, 23, 0.68), rgba(2, 6, 23, 0.68)),
                              url('data:image/png;base64,{encoded}');
            background-size: cover;
            background-position: center;
            background-attachment: fixed;
        }}
        """

    return """
    .stApp {
        background: radial-gradient(circle at 20% 10%, #334155 0%, #0f172a 55%, #020617 100%);
    }
    """


def render():
    st.markdown(
        f"""
        <style>
        {_background_style()}

        .stApp, .stApp * {{
            font-family: "Inter", "Segoe UI", "Helvetica Neue", Arial, sans-serif;
        }}

        .hero-title {{
            text-align: center;
            font-size: 4.3rem;
            font-weight: 900;
            color: #f8fafc;
            letter-spacing: 0.02em;
            margin-top: 0.6rem;
            margin-bottom: 0.3rem;
            text-shadow: 0 6px 30px rgba(0, 0, 0, 0.45);
        }}

        .hero-sub {{
            text-align: center;
            color: #dbeafe;
            font-size: 1.14rem;
            margin-bottom: 1.25rem;
        }}

        .workflow-note {{
            text-align: center;
            color: #dbeafe;
            font-size: 1.02rem;
            line-height: 1.35;
            margin-top: 0.65rem;
        }}
        /* Keep top navigation tabs visually identical to the other pages. */
        div[role="radiogroup"] > label {{
            background: rgba(15, 23, 42, 0.82) !important;
            border: 1px solid rgba(148, 163, 184, 0.35) !important;
            border-radius: 10px !important;
            padding: 0.3rem 0.35rem !important;
        }}

        div[role="radiogroup"] > label p {{
            color: #e2e8f0 !important;
            font-weight: 700 !important;
            font-size: 1.02rem !important;
        }}

        div[role="radiogroup"] > label[data-selected="true"] {{
            border-color: #7dd3fc !important;
            background: linear-gradient(180deg, rgba(37, 99, 235, 0.96), rgba(30, 64, 175, 0.96)) !important;
            box-shadow: 0 8px 20px rgba(37, 99, 235, 0.35) !important;
        }}


        div[data-testid="stButton"] > button {{
            min-height: 150px;
            border-radius: 14px;
            font-size: 2rem;
            font-weight: 900;
            border: 1px solid rgba(148, 163, 184, 0.42);
            background: linear-gradient(180deg, rgba(30, 41, 80, 0.98) 0%, rgba(23, 32, 63, 0.98) 100%);
            color: #f8fafc;
            box-shadow: 0 10px 24px rgba(15, 23, 42, 0.42);
            line-height: 1.15;
            white-space: nowrap;
            text-wrap: balance;
        }}

        div[data-testid="stButton"] > button p {{
            margin: 0;
            line-height: 1.15;
            font-weight: 900;
            font-size: 1.75rem;
            letter-spacing: 0.01em;
        }}

        div[data-testid="stButton"] > button:hover {{
            border-color: #60a5fa;
            background: linear-gradient(180deg, rgba(37, 99, 235, 0.92) 0%, rgba(30, 64, 175, 0.92) 100%);
            color: #ffffff;
        }}
        </style>
        """,
        unsafe_allow_html=True,
    )

    st.markdown('<div class="hero-title">QPrimer</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="hero-sub">Choose one workflow below and jump directly into design.</div>',
        unsafe_allow_html=True,
    )

    c1, c2, c3 = st.columns(3)

    with c1:
        go_regular = st.button("Regular PCR", key="home_regular", width="stretch")
        st.markdown('<div class="workflow-note">Design a standard forward/reverse primer pair.</div>', unsafe_allow_html=True)

    with c2:
        go_splicing = st.button("Splicing primers", key="home_splicing", width="stretch")
        st.markdown('<div class="workflow-note">Build overlap/splicing primer designs.</div>', unsafe_allow_html=True)

    with c3:
        go_qpcr = st.button("qPCR primers", key="home_qpcr", width="stretch")
        st.markdown('<div class="workflow-note">Generate qPCR-friendly primer candidates.</div>', unsafe_allow_html=True)

    if go_regular:
        st.session_state["requested_tab"] = "Regular PCR"
        st.rerun()
    if go_splicing:
        st.session_state["requested_tab"] = "Splicing primers"
        st.rerun()
    if go_qpcr:
        st.session_state["requested_tab"] = "qPCR primers"
        st.rerun()
