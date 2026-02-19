import streamlit as st

FOOTER_TEXT = (
    "This tool is free and open to everyone. "
    "It was designed and built by Khalid Alabdulla. "
    "If you find it useful, please consider sharing it. "
    "For feedback or bug reports, contact "
    "[alabdulla8932@gmail.com](mailto:alabdulla8932@gmail.com)."
)

def add_footer():
    st.markdown("---")
    st.caption(FOOTER_TEXT)
