import streamlit as st


def render():
    st.markdown(
        """
        <style>
        .about-card {
            background: rgba(15, 23, 42, 0.78);
            border: 1px solid rgba(148, 163, 184, 0.35);
            border-radius: 14px;
            padding: 1.35rem 1.2rem;
            color: #e2e8f0;
            box-shadow: 0 10px 24px rgba(2, 6, 23, 0.35);
            line-height: 1.55;
        }

        .about-card h2 {
            margin-top: 0;
            color: #f8fafc;
            font-size: 1.8rem;
            margin-bottom: 0.75rem;
        }

        .about-card h3 {
            color: #cbd5e1;
            font-size: 1.15rem;
            margin-top: 1.1rem;
            margin-bottom: 0.4rem;
        }

        .about-card p {
            margin: 0.35rem 0;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div class="about-card">
            <h2>About Us</h2>
            <p>
                My name is <strong>Khalid Alabdulla</strong> from <strong>Carnegie Mellon University Qatar</strong>.
                I created this website because I found it difficult to find easy-to-use websites for designing
                splicing and qPCR primers.
            </p>
            <p>
                Hopefully this website is helpful for students, researchers, and anyone working with primer design.
                I would also like to thank my professors for their guidance and support.
            </p>

            References:
            Images used were created with BioRender.com

            Citation request:
            We request, but do not require, that use of this software be cited in publications as:
            Alabdulla, K. (n.d.). PrimerQ [Software].
            
        </div>
        """,
        unsafe_allow_html=True,
    )
