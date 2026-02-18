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

# ---------------- Simple qPCR junction design helpers ----------------

DNA = set("ACGT")

def clean_dna(s: str) -> str:
    return "".join([c for c in s.upper() if c in DNA or c == "|"])

def gc_pct(seq: str) -> float:
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = sum(1 for c in seq if c in ("G", "C"))
    return 100.0 * gc / len(seq)

def tm_wallace(seq: str) -> float:
    seq = seq.upper()
    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")
    return 2.0 * (a + t) + 4.0 * (g + c)

def revcomp(seq: str) -> str:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp.get(b, "N") for b in seq.upper()[::-1])

def has_bad_runs(seq: str, max_run: int = 4) -> bool:
    seq = seq.upper()
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
            if run > max_run:
                return True
        else:
            run = 1
    return False

def primer_score(seq: str, tm_target: float) -> float:
    tm = tm_wallace(seq)
    score = abs(tm - tm_target)
    if has_bad_runs(seq, max_run=4):
        score += 5.0
    gcp = gc_pct(seq)
    if gcp < 35 or gcp > 65:
        score += 3.0
    return score

def design_qpcr_junction_primers(full_with_bar: str,
                                min_len: int,
                                max_len: int,
                                tm_target: float,
                                tm_tol: float,
                                min_junction_overlap: int,
                                amplicon_min: int,
                                amplicon_max: int,
                                downstream_window: int):
    s = clean_dna(full_with_bar)
    if s.count("|") != 1:
        raise ValueError("For qPCR: paste ONE full sequence and mark the junction with exactly one '|' like EXON1|EXON2.")

    j = s.index("|")
    left = s[:j]
    right = s[j + 1:]
    if len(left) < min_junction_overlap or len(right) < min_junction_overlap:
        raise ValueError("Not enough bases on one side of the junction for the overlap you chose.")

    full = left + right
    junction_pos = len(left)

    # Forward primer must span junction: take x bases from left tail + y bases from right head
    fwd_candidates = []
    for L in range(min_len, max_len + 1):
        for x in range(min_junction_overlap, min(L - min_junction_overlap, len(left)) + 1):
            y = L - x
            if y < min_junction_overlap or y > len(right):
                continue
            seq = left[-x:] + right[:y]
            tm = tm_wallace(seq)
            if abs(tm - tm_target) <= tm_tol:
                fwd_candidates.append(seq)

    if not fwd_candidates:
        raise ValueError("No junction-spanning forward primer found. Try increasing Tm tolerance or length range.")

    fwd_best = min(fwd_candidates, key=lambda q: primer_score(q, tm_target))

    # Reverse primer: choose in right side downstream, then reverse-complement it
    search_start = junction_pos + max(amplicon_min - len(fwd_best), 10)
    search_end = min(len(full), junction_pos + downstream_window)

    if search_start >= search_end - min_len:
        raise ValueError("Downstream search window is too small for a reverse primer.")

    rev_candidates = []
    for L in range(min_len, max_len + 1):
        for start in range(search_start, search_end - L + 1):
            site = full[start:start + L]  # plus strand site
            primer = revcomp(site)
            tm = tm_wallace(primer)
            if abs(tm - tm_target) <= tm_tol:
                amplicon_len = (start + L) - (junction_pos - (len(fwd_best) // 2))
                if amplicon_min <= amplicon_len <= amplicon_max:
                    rev_candidates.append((primer, amplicon_len))

    if not rev_candidates:
        raise ValueError("No reverse primer found that matches amplicon constraints. Try wider amplicon range or downstream window.")

    rev_best, best_amp = min(rev_candidates, key=lambda x: primer_score(x[0], tm_target))

    return fwd_best, rev_best, best_amp

# ---------------- Streamlit UI ----------------

st.set_page_config(page_title="Primer Designer", layout="wide")

tabs = st.tabs(["Alternative splicing (Exon primers)", "qPCR (junction primers)"])

with st.sidebar:
    st.header("Primer settings")
    min_len = st.number_input("Min primer length", 16, 30, 16)
    max_len = st.number_input("Max primer length", 16, 40, 25)
    tm_target = st.number_input("Target Tm (°C)", 50.0, 70.0, 60.0)
    tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 15.0, 5.0)

# -------- Tab 1: Alternative splicing (your current exon tool) --------

with tabs[0]:
    st.title("Exon Primer Design Tool")
    st.write("Design exon-specific primers with Tm optimization, dimer checks, and direct Primer-BLAST links.")

    with st.sidebar:
        st.subheader("Alt splicing options")
        dimer_k = st.number_input("3' dimer check k", 3, 8, 4)

    st.subheader("Paste exon sequences (A/C/G/T only)")
    exon1 = st.text_area("Exon 1 sequence", height=150, key="as_ex1")
    exon2 = st.text_area("Exon 2 sequence", height=150, key="as_ex2")
    exon3 = st.text_area("Exon 3 sequence (reverse primer)", height=150, key="as_ex3")

    run = st.button("Design primers", key="as_run")

    if run:
        if not exon1 or not exon2 or not exon3:
            st.error("Please paste all three exon sequences.")
        else:
            try:
                p1, p2, p3 = design_exon_primers(
                    exon1, exon2, exon3,
                    min_len=min_len,
                    max_len=max_len,
                    tm_target=tm_target,
                    tm_tol=tm_tol,
                    dimer_k=dimer_k
                )

                st.success("Primers designed successfully.")

                rows = [
                    {"Type": "FWD", "Exon": "Exon 1", "Primer (5'→3')": p1.seq_5to3, "Length": p1.length, "Tm (°C)": round(p1.tm_c, 1), "GC (%)": round(p1.gc_pct, 1), "Score": round(p1.score, 2)},
                    {"Type": "FWD", "Exon": "Exon 2", "Primer (5'→3')": p2.seq_5to3, "Length": p2.length, "Tm (°C)": round(p2.tm_c, 1), "GC (%)": round(p2.gc_pct, 1), "Score": round(p2.score, 2)},
                    {"Type": "REV", "Exon": "Exon 3", "Primer (5'→3')": p3.seq_5to3, "Length": p3.length, "Tm (°C)": round(p3.tm_c, 1), "GC (%)": round(p3.gc_pct, 1), "Score": round(p3.score, 2)},
                ]
                st.dataframe(pd.DataFrame(rows), use_container_width=True)

                st.subheader("Primer-BLAST links (NCBI)")
                org = "Homo sapiens"
                st.markdown(f"**Exon 1 (FWD) + Exon 3 (REV)**: [Open in Primer-BLAST]({primer_blast_url_pair(p1.seq_5to3, p3.seq_5to3, org)})")
                st.markdown(f"**Exon 2 (FWD) + Exon 3 (REV)**: [Open in Primer-BLAST]({primer_blast_url_pair(p2.seq_5to3, p3.seq_5to3, org)})")

                with st.expander("Single-primer Primer-BLAST links"):
                    st.markdown(f"**Exon 1 (FWD)**: [Primer-BLAST]({primer_blast_url_single(p1.seq_5to3, org)})")
                    st.markdown(f"**Exon 2 (FWD)**: [Primer-BLAST]({primer_blast_url_single(p2.seq_5to3, org)})")
                    st.markdown(f"**Exon 3 (REV)**: [Primer-BLAST]({primer_blast_url_single(p3.seq_5to3, org)})")

                st.subheader("Dimer check (forward vs reverse)")
                print_dimer_report(p1, p2, p3)

            except Exception as e:
                st.error(str(e))

# -------- Tab 2: qPCR junction primers --------

with tabs[1]:
    st.title("qPCR Junction Primer Tool")
    st.write("Design qPCR primers where the forward primer spans an exon–exon junction.")

    st.markdown(
        """
**How to get a transcript sequence with exon boundaries**
- Use AceView (NCBI): https://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/index.html?human  
- Search your gene and open a transcript to view the spliced mRNA/cDNA sequence.
- Copy the sequence and identify the exon boundary you want to target.

**How to use this tab**
1) Paste the full spliced sequence (cDNA) into the box below.  
2) Mark the exon junction with a single `|` character (exactly one).  
3) Example junction formatting:
- `...AAGGACCTGATGCTGAC|GTTCCAGGAGTCTGACT...`
- Left of `|` = upstream exon end
- Right of `|` = downstream exon start

**What the tool does**
- Designs a junction-spanning **forward** primer (must include bases from both sides of `|`).  
- Designs a **reverse** primer downstream to hit your amplicon size range.  
- Provides a paired Primer-BLAST link (FWD + REV) so you can check specificity.
        """
    )

    st.info(
        "Tip: If no primer is found, try increasing Tm tolerance, widening length range, or increasing the downstream window."
    )

    with st.sidebar:
        st.subheader("qPCR options")
        min_junc = st.number_input("Min bases overlapping each side of junction", 4, 12, 6)
        amp_min = st.number_input("Amplicon min (bp)", 50, 300, 70)
        amp_max = st.number_input("Amplicon max (bp)", 80, 600, 200)
        down_win = st.number_input("Downstream search window (bp)", 150, 2000, 600)

    seq_full = st.text_area(
        "Full sequence with junction marker |",
        height=220,
        key="qpcr_full",
        placeholder="Paste here, e.g.\n...AAGGACCTGATGCTGAC|GTTCCAGGAGTCTGACT..."
    )

    run_q = st.button("Design junction primers", key="qpcr_run")

    if run_q:
        if not seq_full:
            st.error("Please paste the full sequence and include a single '|' at the junction.")
        else:
            try:
                fwd, rev, amp_len = design_qpcr_junction_primers(
                    seq_full,
                    min_len=min_len,
                    max_len=max_len,
                    tm_target=tm_target,
                    tm_tol=tm_tol,
                    min_junction_overlap=min_junc,
                    amplicon_min=amp_min,
                    amplicon_max=amp_max,
                    downstream_window=down_win
                )

                st.success("qPCR primers designed successfully.")

                rows = [
                    {"Type": "FWD (junction)", "Primer (5'→3')": fwd, "Length": len(fwd), "Tm (°C)": round(tm_wallace(fwd), 1), "GC (%)": round(gc_pct(fwd), 1), "Score": round(primer_score(fwd, tm_target), 2)},
                    {"Type": "REV", "Primer (5'→3')": rev, "Length": len(rev), "Tm (°C)": round(tm_wallace(rev), 1), "GC (%)": round(gc_pct(rev), 1), "Score": round(primer_score(rev, tm_target), 2)},
                ]
                st.dataframe(pd.DataFrame(rows), use_container_width=True)

                st.write(f"Estimated amplicon length (approx): **{amp_len} bp**")

                st.subheader("Primer-BLAST link (pair)")
                org = "Homo sapiens"
                st.markdown(f"[Open in Primer-BLAST]({primer_blast_url_pair(fwd, rev, org)})")

            except Exception as e:
                st.error(str(e))
