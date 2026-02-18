import streamlit as st
import pandas as pd
import urllib.parse

from primer_engine import (
    design_exon_primers,
    print_dimer_report
)

# ---------------- Primer-BLAST URL helpers ----------------

def primer_blast_url_pair(fwd_seq: str, rev_seq: str, organism: str = "Homo sapiens") -> str:
    params = {
        "PRIMER_LEFT_INPUT": fwd_seq,
        "PRIMER_RIGHT_INPUT": rev_seq,
        "ORGANISM": organism
    }
    return "https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?" + urllib.parse.urlencode(params)

def primer_blast_url_single(primer_seq: str, organism: str = "Homo sapiens") -> str:
    params = {
        "PRIMER_LEFT_INPUT": primer_seq,
        "ORGANISM": organism
    }
    return "https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?" + urllib.parse.urlencode(params)

# ---------------- Streamlit UI ----------------

st.set_page_config(page_title="Exon Primer Designer", layout="wide")

st.title("Exon Primer Design Tool")
st.write(
    "Design exon-specific primers with Tm optimization, dimer checks, "
    "and direct Primer-BLAST specificity links."
)

with st.sidebar:
    st.header("Primer settings")

    min_len = st.number_input("Min primer length", 16, 30, 16)
    max_len = st.number_input("Max primer length", 16, 40, 25)
    tm_target = st.number_input("Target Tm (°C)", 50.0, 70.0, 60.0)
    tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 10.0, 5.0)
    dimer_k = st.number_input("3' dimer check k", 3, 8, 4)

st.subheader("Paste exon sequences (A/C/G/T only)")

exon1 = st.text_area("Exon 1 sequence", height=150)
exon2 = st.text_area("Exon 2 sequence", height=150)
exon3 = st.text_area("Exon 3 sequence (reverse primer)", height=150)

run = st.button("Design primers")

if run:
    if not exon1 or not exon2 or not exon3:
        st.error("Please paste all three exon sequences.")
    else:
        try:
            p1, p2, p3 = design_exon_primers(
                exon1,
                exon2,
                exon3,
                min_len=min_len,
                max_len=max_len,
                tm_target=tm_target,
                tm_tol=tm_tol,
                dimer_k=dimer_k
            )

            st.success("Primers designed successfully.")

            # Results table
            rows = [
                {
                    "Type": "FWD",
                    "Exon": "Exon 1",
                    "Primer (5'→3')": p1.seq_5to3,
                    "Length": p1.length,
                    "Tm (°C)": round(p1.tm_c, 1),
                    "GC (%)": round(p1.gc_pct, 1),
                    "Score": round(p1.score, 2)
                },
                {
                    "Type": "FWD",
                    "Exon": "Exon 2",
                    "Primer (5'→3')": p2.seq_5to3,
                    "Length": p2.length,
                    "Tm (°C)": round(p2.tm_c, 1),
                    "GC (%)": round(p2.gc_pct, 1),
                    "Score": round(p2.score, 2)
                },
                {
                    "Type": "REV",
                    "Exon": "Exon 3",
                    "Primer (5'→3')": p3.seq_5to3,
                    "Length": p3.length,
                    "Tm (°C)": round(p3.tm_c, 1),
                    "GC (%)": round(p3.gc_pct, 1),
                    "Score": round(p3.score, 2)
                }
            ]

            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True)

            # Primer-BLAST links
            st.subheader("Primer-BLAST links (NCBI)")

            organism_name = "Homo sapiens"
            url_exon1_pair = primer_blast_url_pair(p1.seq_5to3, p3.seq_5to3, organism_name)
            url_exon2_pair = primer_blast_url_pair(p2.seq_5to3, p3.seq_5to3, organism_name)

            st.markdown(f"**Exon 1 (FWD) + Exon 3 (REV)**: [Open in Primer-BLAST]({url_exon1_pair})")
            st.markdown(f"**Exon 2 (FWD) + Exon 3 (REV)**: [Open in Primer-BLAST]({url_exon2_pair})")

            with st.expander("Single-primer Primer-BLAST links"):
                st.markdown(f"**Exon 1 (FWD)**: [Primer-BLAST]({primer_blast_url_single(p1.seq_5to3, organism_name)})")
                st.markdown(f"**Exon 2 (FWD)**: [Primer-BLAST]({primer_blast_url_single(p2.seq_5to3, organism_name)})")
                st.markdown(f"**Exon 3 (REV)**: [Primer-BLAST]({primer_blast_url_single(p3.seq_5to3, organism_name)})")

            # Dimer report
            st.subheader("Dimer check (forward vs reverse)")
            print_dimer_report(p1, p2, p3)

        except Exception as e:
            st.error(str(e))
