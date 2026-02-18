import streamlit as st
import pandas as pd
import urllib.parse

from primer_engine import (
    design_exon_primers,
    print_dimer_report
)

# ===========================
# Text blocks: instructions + score + footer
# ===========================

BLAST_INSTRUCTIONS = (
    "How to use Primer-BLAST:\n"
    "1) Click the **Open in Primer-BLAST** link.\n"
    "2) On the NCBI page, do NOT change any settings.\n"
    "3) Just click **Get Primers**.\n"
    "4) Check the top hit matches your intended gene/transcript.\n"
    "This is a quick specificity check."
)

SCORE_EXPLANATION = (
    "Score (lower is better): internal ranking only.\n"
    "It is based on:\n"
    "• |Tm − target Tm|\n"
    "• +5 penalty for long homopolymer runs (e.g. AAAAA)\n"
    "• +3 penalty for extreme GC% (<35% or >65%)\n"
    "This is not a BLAST score and not experimental validation."
)

FOOTER_TEXT = (
    "This tool is free and open to everyone. "
    "It was designed and built by Khalid Alabdulla. "
    "If you find it useful, please consider sharing it."
    ""
    "If you have feedback or find a bug, email alabdulla8932@gmail.com"
)

def add_footer():
    st.markdown("---")
    st.caption(FOOTER_TEXT)

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

# ---------------- Simple primer helpers ----------------

DNA = set("ACGT")

def clean_dna(s: str) -> str:
    # keep A/C/G/T and also allow '|' for junction tab
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
    """
    Lower is better.
    This is a simple heuristic: closeness to target Tm,
    penalties for long homopolymers and extreme GC%.
    """
    tm = tm_wallace(seq)
    score = abs(tm - tm_target)

    if has_bad_runs(seq, max_run=4):
        score += 5.0

    gcp = gc_pct(seq)
    if gcp < 35 or gcp > 65:
        score += 3.0

    return score

# ---------------- qPCR junction design ----------------

def design_qpcr_junction_primers(
    full_with_bar: str,
    min_len: int,
    max_len: int,
    tm_target: float,
    tm_tol: float,
    min_junction_overlap: int,
    amplicon_min: int,
    amplicon_max: int,
    downstream_window: int
):
    s = clean_dna(full_with_bar)
    if s.count("|") != 1:
        raise ValueError("Paste ONE sequence and mark the junction with exactly one '|' like EXON1|EXON2.")

    j = s.index("|")
    left = s[:j]
    right = s[j + 1:]
    if len(left) < min_junction_overlap or len(right) < min_junction_overlap:
        raise ValueError("Not enough bases on one side of the junction for the overlap you chose.")

    full = left + right
    junction_pos = len(left)

    # Forward primer spans junction
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

    # Reverse primer downstream
    search_start = junction_pos + max(amplicon_min - len(fwd_best), 10)
    search_end = min(len(full), junction_pos + downstream_window)

    if search_start >= search_end - min_len:
        raise ValueError("Downstream window too small for a reverse primer. Increase downstream window.")

    rev_candidates = []
    for L in range(min_len, max_len + 1):
        for start in range(search_start, search_end - L + 1):
            site = full[start:start + L]
            primer = revcomp(site)
            tm = tm_wallace(primer)
            if abs(tm - tm_target) <= tm_tol:
                amplicon_len = (start + L) - (junction_pos - (len(fwd_best) // 2))
                if amplicon_min <= amplicon_len <= amplicon_max:
                    rev_candidates.append((primer, amplicon_len))

    if not rev_candidates:
        raise ValueError("No reverse primer found. Try wider amplicon range or downstream window.")

    rev_best, best_amp = min(rev_candidates, key=lambda x: primer_score(x[0], tm_target))
    return fwd_best, rev_best, best_amp

# ---------------- Basic PCR design (one sequence -> FWD + REV) ----------------

def design_basic_pcr_primers(
    template: str,
    min_len: int,
    max_len: int,
    tm_target: float,
    tm_tol: float,
    start_window: int,
    end_window: int,
    amplicon_min: int,
    amplicon_max: int
):
    s = clean_dna(template).replace("|", "")
    if len(s) < max(amplicon_min, 80):
        raise ValueError("Sequence is too short. Paste a longer template sequence.")

    start_window = min(start_window, len(s))
    fwd_candidates = []
    for L in range(min_len, max_len + 1):
        for i in range(0, max(0, start_window - L + 1)):
            seq = s[i:i + L]
            tm = tm_wallace(seq)
            if abs(tm - tm_target) <= tm_tol:
                fwd_candidates.append((seq, i))

    if not fwd_candidates:
        raise ValueError("No forward primer found. Increase Tm tolerance or start window or length range.")

    end_window = min(end_window, len(s))
    rev_candidates = []
    end_start = max(0, len(s) - end_window)
    for L in range(min_len, max_len + 1):
        for i in range(end_start, len(s) - L + 1):
            site = s[i:i + L]
            primer = revcomp(site)
            tm = tm_wallace(primer)
            if abs(tm - tm_target) <= tm_tol:
                rev_candidates.append((primer, i, L))

    if not rev_candidates:
        raise ValueError("No reverse primer found. Increase Tm tolerance or end window or length range.")

    best = None
    best_key = None
    for fwd, f_i in fwd_candidates:
        for rev, r_i, r_L in rev_candidates:
            amp_len = (r_i + r_L) - f_i
            if amplicon_min <= amp_len <= amplicon_max:
                key = (
                    primer_score(fwd, tm_target) + primer_score(rev, tm_target),
                    abs(amp_len - ((amplicon_min + amplicon_max) / 2))
                )
                if best is None or key < best_key:
                    best = (fwd, rev, amp_len, f_i, r_i)
                    best_key = key

    if best is None:
        raise ValueError("Could not find a primer pair that matches your amplicon size. Widen amplicon range or windows.")

    return best  # fwd, rev, amp_len, fwd_start, rev_site_start

# ---------------- Streamlit UI ----------------

st.set_page_config(page_title="Primer Designer", layout="wide")

tabs = st.tabs([
    "Basic PCR (two primers from one sequence)",
    "Alternative splicing (Exon primers)",
    "qPCR (junction primers)"
])

# -------- Tab 1: Basic PCR --------

with tabs[0]:
    st.title("Basic PCR Primer Tool")
    st.write("Paste one DNA sequence and get a forward and reverse primer (no junction required).")

    with st.expander("Edit primer conditions", expanded=False):
        c1, c2 = st.columns(2)
        with c1:
            min_len = st.number_input("Min primer length", 16, 30, 16, key="basic_min_len")
            tm_target = st.number_input("Target Tm (°C)", 50.0, 70.0, 60.0, key="basic_tm_target")
        with c2:
            max_len = st.number_input("Max primer length", 16, 40, 25, key="basic_max_len")
            tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 15.0, 5.0, key="basic_tm_tol")

        st.subheader("Basic PCR options")
        basic_start_win = st.number_input("Forward search window (bp from start)", 50, 3000, 300, key="basic_start_win")
        basic_end_win = st.number_input("Reverse search window (bp from end)", 50, 3000, 300, key="basic_end_win")
        basic_amp_min = st.number_input("Amplicon min (bp)", 50, 5000, 100, key="basic_amp_min")
        basic_amp_max = st.number_input("Amplicon max (bp)", 80, 8000, 400, key="basic_amp_max")

    st.markdown(
        """
**How to use**
1) Paste your template sequence (A/C/G/T only).  
2) The tool searches for a forward primer near the beginning and a reverse primer near the end.  
3) You control the amplicon size range.  
        """
    )

    template = st.text_area(
        "Template sequence (A/C/G/T only)",
        height=240,
        key="basic_template",
        placeholder="Paste a gene region, plasmid region, or any DNA template here..."
    )

    run_b = st.button("Design basic PCR primers", key="basic_run")

    if run_b:
        if not template:
            st.error("Please paste a template sequence.")
        else:
            try:
                fwd, rev, amp_len, f_i, r_i = design_basic_pcr_primers(
                    template,
                    min_len=min_len,
                    max_len=max_len,
                    tm_target=tm_target,
                    tm_tol=tm_tol,
                    start_window=int(basic_start_win),
                    end_window=int(basic_end_win),
                    amplicon_min=int(basic_amp_min),
                    amplicon_max=int(basic_amp_max),
                )

                st.success("Basic PCR primers designed successfully.")
                st.caption(SCORE_EXPLANATION)

                rows = [
                    {"Type": "FWD", "Primer (5'→3')": fwd, "Length": len(fwd), "Tm (°C)": round(tm_wallace(fwd), 1), "GC (%)": round(gc_pct(fwd), 1), "Score": round(primer_score(fwd, tm_target), 2)},
                    {"Type": "REV", "Primer (5'→3')": rev, "Length": len(rev), "Tm (°C)": round(tm_wallace(rev), 1), "GC (%)": round(gc_pct(rev), 1), "Score": round(primer_score(rev, tm_target), 2)},
                ]
                st.dataframe(pd.DataFrame(rows), use_container_width=True)

                st.write(f"Estimated amplicon length (approx): **{amp_len} bp**")

                st.subheader("Primer-BLAST link (pair)")
                org = "Homo sapiens"
                st.markdown(f"[Open in Primer-BLAST]({primer_blast_url_pair(fwd, rev, org)})")
                st.info(BLAST_INSTRUCTIONS)

            except Exception as e:
                st.error(str(e))

    add_footer()

# -------- Tab 2: Alternative splicing --------

with tabs[1]:
    st.title("Exon Primer Design Tool")
    st.write("Design exon-specific primers with Tm optimization, dimer checks, and direct Primer-BLAST links.")

    with st.expander("Edit primer conditions", expanded=False):
        c1, c2 = st.columns(2)
        with c1:
            min_len = st.number_input("Min primer length", 16, 30, 16, key="as_min_len")
            tm_target = st.number_input("Target Tm (°C)", 50.0, 70.0, 60.0, key="as_tm_target")
            dimer_k = st.number_input("3' dimer check k", 3, 8, 4, key="as_dimer_k")
        with c2:
            max_len = st.number_input("Max primer length", 16, 40, 25, key="as_max_len")
            tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 15.0, 5.0, key="as_tm_tol")

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
                st.caption(SCORE_EXPLANATION)

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
                st.info(BLAST_INSTRUCTIONS)

                with st.expander("Single-primer Primer-BLAST links"):
                    st.markdown(f"**Exon 1 (FWD)**: [Primer-BLAST]({primer_blast_url_single(p1.seq_5to3, org)})")
                    st.markdown(f"**Exon 2 (FWD)**: [Primer-BLAST]({primer_blast_url_single(p2.seq_5to3, org)})")
                    st.markdown(f"**Exon 3 (REV)**: [Primer-BLAST]({primer_blast_url_single(p3.seq_5to3, org)})")

                st.subheader("Dimer check (forward vs reverse)")
                print_dimer_report(p1, p2, p3)

            except Exception as e:
                st.error(str(e))

    add_footer()

# -------- Tab 3: qPCR junction primers --------

with tabs[2]:
    st.title("qPCR Junction Primer Tool")
    st.write("Design qPCR primers where the forward primer spans an exon–exon junction.")

    st.markdown(
        """
**How to get a transcript sequence with exon boundaries**
- Use AceView (NCBI): https://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/index.html?human  
- Search your gene and open a transcript to view the spliced sequence  
- Copy the spliced cDNA sequence and choose the junction you want  

**How to use this tab**
1) Paste the full spliced sequence (cDNA) below  
2) Mark the exon junction with exactly one `|`  
3) Example:
- `...AAGGACCTGATGCTGAC|GTTCCAGGAGTCTGACT...`

**What you get**
- Junction-spanning forward primer  
- Reverse primer downstream (amplicon size controlled)  
- Paired Primer-BLAST link  
        """
    )

    with st.expander("Edit primer conditions", expanded=False):
        c1, c2 = st.columns(2)
        with c1:
            min_len = st.number_input("Min primer length", 16, 30, 16, key="qpcr_min_len")
            tm_target = st.number_input("Target Tm (°C)", 50.0, 70.0, 60.0, key="qpcr_tm_target")
        with c2:
            max_len = st.number_input("Max primer length", 16, 40, 25, key="qpcr_max_len")
            tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 15.0, 5.0, key="qpcr_tm_tol")

        st.subheader("qPCR options")
        min_junc = st.number_input("Min bases overlapping each side of junction", 4, 12, 6, key="qpcr_min_junc")
        amp_min = st.number_input("Amplicon min (bp)", 50, 300, 70, key="qpcr_amp_min")
        amp_max = st.number_input("Amplicon max (bp)", 80, 600, 200, key="qpcr_amp_max")
        down_win = st.number_input("Downstream search window (bp)", 150, 2000, 600, key="qpcr_down_win")

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
                st.caption(SCORE_EXPLANATION)

                rows = [
                    {"Type": "FWD (junction)", "Primer (5'→3')": fwd, "Length": len(fwd), "Tm (°C)": round(tm_wallace(fwd), 1), "GC (%)": round(gc_pct(fwd), 1), "Score": round(primer_score(fwd, tm_target), 2)},
                    {"Type": "REV", "Primer (5'→3')": rev, "Length": len(rev), "Tm (°C)": round(tm_wallace(rev), 1), "GC (%)": round(gc_pct(rev), 1), "Score": round(primer_score(rev, tm_target), 2)},
                ]
                st.dataframe(pd.DataFrame(rows), use_container_width=True)

                st.write(f"Estimated amplicon length (approx): **{amp_len} bp**")

                st.subheader("Primer-BLAST link (pair)")
                org = "Homo sapiens"
                st.markdown(f"[Open in Primer-BLAST]({primer_blast_url_pair(fwd, rev, org)})")
                st.info(BLAST_INSTRUCTIONS)

            except Exception as e:
                st.error(str(e))

    add_footer()
