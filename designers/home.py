import base64
from pathlib import Path

import streamlit as st


def _background_style() -> str:
    image_path = Path("assets/dna_background.png")
    if not image_path.exists():
        return ""

    encoded = base64.b64encode(image_path.read_bytes()).decode("utf-8")
    return f"""
    .stApp {{
        background-image: linear-gradient(rgba(255,255,255,0.90), rgba(255,255,255,0.90)),
                          url('data:image/png;base64,{encoded}');
        background-size: cover;
        background-position: center;
        background-attachment: fixed;
    }}
    """


def render():
    st.markdown(
        f"""
        <style>
        {_background_style()}

        .hero-title {{
            text-align: center;
            font-size: 4rem;
            font-weight: 900;
            color: #0f172a;
            letter-spacing: 0.03em;
            margin-top: 0.4rem;
            margin-bottom: 0.35rem;
        }}

        .hero-sub {{
            text-align: center;
            color: #334155;
            font-size: 1.05rem;
            margin-bottom: 1.25rem;
        }}

        .workflow-label {{
            text-align: center;
            font-weight: 700;
            color: #0f172a;
            margin-bottom: 0.35rem;
        }}

        .workflow-note {{
            text-align: center;
            color: #475569;
            font-size: 0.9rem;
            min-height: 2.8rem;
            margin-bottom: 0.55rem;
        }}

        div[data-testid="stButton"] > button {{
            height: 120px;
            border-radius: 14px;
            font-size: 1.1rem;
            font-weight: 700;
            border: 1px solid #93c5fd;
            background: linear-gradient(180deg, #f8fbff 0%, #eaf2ff 100%);
            color: #0f172a;
        }}

        div[data-testid="stButton"] > button:hover {{
            border-color: #2563eb;
            color: #1d4ed8;
        }}
        </style>
        """,
        unsafe_allow_html=True,
    )

    st.markdown('<div class="hero-title">QPrimer</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="hero-sub">Choose one workflow below. Each rectangle contains the option to start.</div>',
        unsafe_allow_html=True,
    )

    c1, c2, c3 = st.columns(3)

    with c1:
        st.markdown('<div class="workflow-label">Regular PCR</div>', unsafe_allow_html=True)
        st.markdown('<div class="workflow-note">Design a standard forward/reverse primer pair.</div>', unsafe_allow_html=True)
        go_regular = st.button("Open Regular PCR", key="home_regular", width="stretch")

    with c2:
        st.markdown('<div class="workflow-label">Splicing primers</div>', unsafe_allow_html=True)
        st.markdown('<div class="workflow-note">Build overlap/splicing primer designs.</div>', unsafe_allow_html=True)
        go_splicing = st.button("Open Splicing Primers", key="home_splicing", width="stretch")

    with c3:
        st.markdown('<div class="workflow-label">qPCR primers</div>', unsafe_allow_html=True)
        st.markdown('<div class="workflow-note">Generate qPCR-friendly primer candidates.</div>', unsafe_allow_html=True)
        go_qpcr = st.button("Open qPCR Primers", key="home_qpcr", width="stretch")

    if go_regular:
        st.success("Go to the **Regular PCR** tab above to continue.")
    if go_splicing:
        st.success("Go to the **Splicing primers** tab above to continue.")
    if go_qpcr:
        st.success("Go to the **qPCR primers** tab above to continue.")
