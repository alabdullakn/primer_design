import streamlit as st

def apply_styles():
    st.markdown(
        """
        <style>
        /* Wider content */
        .block-container { padding-top: 2rem; padding-bottom: 2rem; }

        /* Slightly nicer headers */
        h1, h2, h3 { letter-spacing: 0.2px; }

        /* Optional: subtle “DNA vibe” background */
        .stApp {
            background: radial-gradient(circle at 20% 10%, rgba(120, 70, 255, 0.10), transparent 35%),
                        radial-gradient(circle at 80% 30%, rgba(0, 200, 255, 0.08), transparent 40%),
                        linear-gradient(180deg, rgba(10,10,12,1) 0%, rgba(6,6,8,1) 100%);
        }

        /* Make expanders a bit cleaner */
        details { border-radius: 12px; }
        </style>
        """,
        unsafe_allow_html=True
    )
