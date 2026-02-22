import base64
from pathlib import Path

import streamlit as st


def _background_style() -> str:
    image_path = Path("assets/dna_background.png")
    if image_path.exists():
        encoded = base64.b64encode(image_path.read_bytes()).decode("utf-8")
        return f"""
        .stApp {{
            background-image: linear-gradient(rgba(2, 6, 23, 0.78), rgba(2, 6, 23, 0.78)),
                              url('data:image/png;base64,{encoded}');
            background-size: cover;
            background-position: center;
            background-attachment: fixed;
        }}
        """

    return """
    .stApp {
        background: radial-gradient(circle at 20% 10%, #1e293b 0%, #020617 55%, #01030a 100%);
    }
    """


def render():
    st.markdown(
        f"""
        <style>
        {_background_style()}

        .hero-title {{
            text-align: center;
            font-size: 4.2rem;
            font-weight: 900;
            color: #f8fafc;
            letter-spacing: 0.04em;
            margin-top: 0.6rem;
            margin-bottom: 0.3rem;
            text-shadow: 0 6px 30px rgba(0, 0, 0, 0.5);
        }}

        .hero-sub {{
            text-align: center;
            color: #cbd5e1;
            font-size: 1.1rem;
            margin-bottom: 1.25rem;
        }}

        .workflow-label {{
            text-align: center;
            font-weight: 700;
            color: #f1f5f9;
            margin-bottom: 0.3rem;
            font-size: 1.9rem;
        }}

        .workflow-note {{
            text-align: center;
            color: #cbd5e1;
            font-size: 0.98rem;
            min-height: 2.8rem;
            margin-bottom: 0.65rem;
        }}

        .workflow-card {{
            background: linear-gradient(180deg, rgba(15, 23, 42, 0.75) 0%, rgba(15, 23, 42, 0.55) 100%);
            border: 1px solid rgba(148, 163, 184, 0.35);
            border-radius: 16px;
            padding: 1.1rem 0.8rem 0.9rem 0.8rem;
            margin-bottom: 0.8rem;
            backdrop-filter: blur(3px);
        }}

        div[data-testid="stButton"] > button {{
            height: 68px;
            border-radius: 12px;
            font-size: 1.15rem;
            font-weight: 700;
            border: 1px solid rgba(96, 165, 250, 0.9);
            background: linear-gradient(180deg, #1d4ed8 0%, #1e40af 100%);
            color: #eff6ff;
            box-shadow: 0 10px 24px rgba(30, 64, 175, 0.45);
        }}

        div[data-testid="stButton"] > button:hover {{
            border-color: #93c5fd;
            background: linear-gradient(180deg, #2563eb 0%, #1d4ed8 100%);
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
        st.markdown('<div class="workflow-card">', unsafe_allow_html=True)
        st.markdown('<div class="workflow-label">Regular PCR</div>', unsafe_allow_html=True)
        st.markdown('<div class="workflow-note">Design a standard forward/reverse primer pair.</div>', unsafe_allow_html=True)
        go_regular = st.button("Open Regular PCR", key="home_regular", width="stretch")
        st.markdown('</div>', unsafe_allow_html=True)

    with c2:
        st.markdown('<div class="workflow-card">', unsafe_allow_html=True)
        st.markdown('<div class="workflow-label">Splicing primers</div>', unsafe_allow_html=True)
        st.markdown('<div class="workflow-note">Build overlap/splicing primer designs.</div>', unsafe_allow_html=True)
        go_splicing = st.button("Open Splicing Primers", key="home_splicing", width="stretch")
        st.markdown('</div>', unsafe_allow_html=True)

    with c3:
        st.markdown('<div class="workflow-card">', unsafe_allow_html=True)
        st.markdown('<div class="workflow-label">qPCR primers</div>', unsafe_allow_html=True)
        st.markdown('<div class="workflow-note">Generate qPCR-friendly primer candidates.</div>', unsafe_allow_html=True)
        go_qpcr = st.button("Open qPCR Primers", key="home_qpcr", width="stretch")
        st.markdown('</div>', unsafe_allow_html=True)

    if go_regular:
        st.session_state["requested_tab"] = "Regular PCR"
        st.rerun()
    if go_splicing:
        st.session_state["requested_tab"] = "Splicing primers"
        st.rerun()
    if go_qpcr:
        st.session_state["requested_tab"] = "qPCR primers"
        st.rerun()
