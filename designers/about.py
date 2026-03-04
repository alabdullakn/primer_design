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
            <h2>About PrimerQ</h2>
            <p>
                My name is <strong>Khalid Alabdulla</strong> from <strong>Carnegie Mellon University Qatar</strong>.
                I created this website because I found it difficult to find easy-to-use tools for designing
                splicing and qPCR primers.
            </p>
            <p>
                PrimerQ supports <strong>Regular PCR</strong>, <strong>Splicing primers</strong> (exon skipping,
                intron retention, alternative splice sites), and <strong>qPCR primers</strong> with both
                junction-spanning and intron-flanking strategies. It includes <strong>TaqMan probe design</strong>,
                <strong>in-app NCBI QBLAST specificity checking</strong>, nearest-neighbour Tm thermodynamics,
                primer-dimer screening, and CSV export.
            </p>
            <p>
                Hopefully this tool is helpful for students, researchers, and anyone working with primer design.
                I would also like to thank my professors for their guidance and support.
            </p>
            <h3>References</h3>
            <p>Images were created with <a href="https://www.biorender.com/" target="_blank" rel="noopener noreferrer">BioRender.com</a>.</p>
            <h3>Citation</h3>
            <p>
                We request, but do not require, that use of this software be cited in publications as:<br>
                <em>Alabdulla, K. (2025). PrimerQ [Software]. Carnegie Mellon University Qatar.</em>
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )
