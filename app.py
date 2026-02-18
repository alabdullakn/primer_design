import streamlit as st
import pandas as pd
import urllib.parse

from primer_engine import (
    design_exon_primers,
    print_dimer_report
)

def primer_blast_url(primer_seq: str, organism: str = "Homo sapiens") -> str:
    params = {
        "PRIMER_LEFT_INPUT": primer_seq,
        "ORGANISM": organism
    }
    return "https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?" + urllib.parse.urlencode(params)


st.set_page_config(page_title="Exon Primer Designer", layout="wide")

st.title("Exon Primer Design Tool")
st.write("Design exon-specific primers with Tm optimization and dimer checks. Includes Primer-BLAST links for specificity.")

with st.sidebar:
    st.header("Primer settings")
    min_len = st.number_input("Min primer length", 16, 30, 16)
    max_len = st.number_input("Max primer length", 16, 40, 25)
    tm_target = st.number_input("Target Tm (°C)", 50.0, 70.0, 60.0)
    tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 10.0, 5.0)
    dimer_k = st.number_input("3' dimer check k", 3, 8, 4)

    organism_name = st.selectbox(
        "Organism for Primer-BLAST links",
        ["Homo sapiens", "Mus musculus", "Rattus norvegicus"],
        index=0
    )

st.subheader("Paste exon sequences (A/C/G/T only)")

exon1 = st.text_area("Exon 1 sequence", height=150)
exon2 = st.text_area("Exon 2 sequence", height=150)
exon3 = st.text_area("Exon for reverse primer (Exon 3)", height=150)

run = st.button("Design primers")

if run:
    if not exon1 or not exon2 or not exon3:
        st.error("Please paste all three exon sequences.")
    else:
        try:
            f1, f2, r3 = design_exon_primers(
                exon1,
                exon2,
                exon3,
                min_len=min_len,
                max_len=max_len,
                tm_target=tm_target,
                tm_tol=tm_tol,
                dimer_k=dimer_k
            )

            primers = [f1, f2, r3]
            st.success("Primers designed successfully.")

            rows = []
            for p in primers:
                rows.append({
                    "Type": p.kind,
                    "Exon": p.exon_name,
                    "Primer (5'→3')": p.seq_5to3,
                    "Length": p.length,
                    "Tm (°C)": round(p.tm_c, 1),
                    "GC (%)": round(p.gc_pct, 1),
                    "Score": round(p.score, 2)
                })

            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True)

            st.subheader("Primer-BLAST links (NCBI)")
            for p in primers:
                url = primer_blast_url(p.seq_5to3, organism=organism_name)
                st.markdown(f"**{p.kind} ({p.exon_name})**: [Open in Primer-BLAST]({url})")

            st.subheader("Dimer check (forward vs reverse)")
            print_dimer_report(f1, f2, r3)

        except Exception as e:
            st.error(str(e))
