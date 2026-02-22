import streamlit as st
def add_footer():
    st.markdown(
        """
        <style>
        .stApp {
            padding-bottom: 92px;
        }

        .qp-footer {
            position: fixed;
            left: 0;
            right: 0;
            bottom: 0;
            z-index: 9999;
            padding: 0.75rem 1rem;
            border-top: 1px solid rgba(148, 163, 184, 0.25);
            background: linear-gradient(180deg, rgba(2, 6, 23, 0.9), rgba(2, 6, 23, 0.98));
            color: rgba(226, 232, 240, 0.92);
            font-size: 0.96rem;
            line-height: 1.45;
        }

        .qp-footer a {
            color: #60a5fa;
            text-decoration: underline;
        }
        </style>

        <div class="qp-footer">
            This tool is free and open to everyone. It was designed and built by Khalid Alabdulla.
            If you find it useful, please consider sharing it. For feedback or bug reports, contact
            <a href="mailto:alabdulla8932@gmail.com">alabdulla8932@gmail.com</a>. All images were created through
            <a href="https://www.biorender.com/" target="_blank" rel="noopener noreferrer">BioRender</a>.
        </div>
        """,
        unsafe_allow_html=True,
    )
