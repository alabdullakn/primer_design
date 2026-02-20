import streamlit as st
import pandas as pd
import requests

from utils.primer_utils import gc_pct, tm_wallace, revcomp, primer_score
from utils.blast import primer_blast_url_pair
from ui.text import BLAST_INSTRUCTIONS, SCORE_EXPLANATION
from ui.footer import add_footer


# ----------------------------
# FASTA + NCBI helpers
# ----------------------------

def _parse_fasta_text(text: str) -> str:
    if not text:
        return ""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    seq = []
    for l in lines:
        if l.startswith(">"):
            continue
        seq.append(l)
    raw = "".join(seq).upper()
    allowed = set("ACGTN")
    return "".join([c for c in raw if c in allowed])


def _fetch_ncbi_fasta(accession: str) -> str:
    acc = (accession or "").strip()
    if not acc:
        return ""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={acc}&rettype=fasta&retmode=text"
    )
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        raise Exception("Failed to fetch accession from NCBI.")
    return r.text


# ----------------------------
# Simple regular PCR designer
# ----------------------------

def _max_run(seq: str) -> int:
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def _is_valid_primer(seq: str, tm_target: float, tm_tol: float) -> bool:
    if not seq:
        return False
    tm = tm_wallace(seq)
    if abs(tm - tm_target) > tm_tol:
        return False
    gc = gc_pct(seq)
    if gc < 35 or gc > 65:
        return False
    if _max_run(seq) >= 5:
        return False
    return True


def design_basic_pcr_primers(
    template: str,
    min_len: int = 18,
    max_len: int = 25,
    tm_target: float = 60.0,
    tm_tol: float = 5.0,
    start_window: int = 300,
    end_window: int = 300,
    amplicon_min: int = 100,
    amplicon_max: int = 500,
):
    """
    Returns: (fwd_5to3, rev_5to3, amplicon_len)
    Heuristic:
      - choose FWD from first start_window
      - choose REV from last end_window
      - enforce amplicon size between min/max
    """
    seq = _parse_fasta_text(template)
    if len(seq) < (min_len + min_len + amplicon_min):
        raise Exception("Template too short for your constraints.")

    left = seq[: min(len(seq), start_window)]
    right = seq[max(0, len(seq) - end_window) :]
    right_start = len(seq) - len(right)

    best = None  # (score, fwd, rev, amp_len)

    for Lf in range(min_len, max_len + 1):
        for i in range(0, len(left) - Lf + 1):
            fwd = left[i : i + Lf]
            if not _is_valid_primer(fwd, tm_target, tm_tol):
                continue
            fwd_tm = tm_wallace(fwd)

            for Lr in range(min_len, max_len + 1):
                for j in range(0, len(right) - Lr + 1):
                    bind_site = right[j : j + Lr]
                    rev = revcomp(bind_site)
                    if not _is_valid_primer(rev, tm_target, tm_tol):
                        continue
                    rev_tm = tm_wallace(rev)

                    amp_len = (right_start + j + Lr) - i
                    if amp_len < amplicon_min or amp_len > amplicon_max:
                        continue

                    score = abs(fwd_tm - tm_target) + abs(rev_tm - tm_target)
                    score += abs(fwd_tm - rev_tm) * 0.5
                    score += (amp_len / 1000.0)

                    if best is None or score < best[0]:
                        best = (score, fwd, rev, amp_len)

    if best is None:
        raise Exception("No primer pair found. Try widening Tm tolerance, length range, or amplicon range.")

    _, fwd, rev, amp_len = best
    return fwd, rev, amp_len


def _quick_dimer_note(fwd: str, rev: str) -> str:
    k = 4
    if len(fwd) >= k and revcomp(fwd[-k:]) in rev:
        return "Warning: possible 3' complementarity (k=4)."
    if len(rev) >= k and revcomp(rev[-k:]) in fwd:
        return "Warning: possible 3' complementarity (k=4)."
    return "No obvious 3' complementarity (k=4) detected."


# ----------------------------
# Streamlit render
# ----------------------------

def render():
    st.title("Regular PCR")
    st.write("Design one forward and one reverse primer from a single template sequence.")

    # Session state
    if "reg_seq" not in st.session_state:
        st.session_state["reg_seq"] = ""
    if "reg_status" not in st.session_state:
        st.session_state["reg_status"] = ""

    # Params
    with st.expander("Primer design parameters", expanded=False):
        c1, c2 = st.columns(2)
        with c1:
            min_len = st.number_input("Min primer length", 16, 40, 18, key="reg_min_len")
            max_len = st.number_input("Max primer length", 16, 60, 25, key="reg_max_len")
            tm_target = st.number_input("Target Tm (°C)", 45.0, 75.0, 60.0, key="reg_tm_target")
        with c2:
            tm_tol = st.number_input("Tm tolerance (± °C)", 1.0, 20.0, 5.0, key="reg_tm_tol")
            amp_min = st.number_input("Amplicon min (bp)", 50, 5000, 100, key="reg_amp_min")
            amp_max = st.number_input("Amplicon max (bp)", 80, 8000, 500, key="reg_amp_max")

    st.markdown("---")

    st.subheader("Sequence input")

    mode = st.radio(
        "Choose input type",
        ["Paste sequence", "Upload FASTA", "NCBI accession"],
        horizontal=True,
        key="reg_input_mode",
    )

    if mode == "Paste sequence":
        raw = st.text_area(
            "Paste template sequence (A/C/G/T only)",
            height=220,
            key="reg_paste_box",
            placeholder="Paste sequence here...",
        )
        seq = _parse_fasta_text(">x\n" + raw)
        if seq:
            st.session_state["reg_seq"] = seq
            st.session_state["reg_status"] = f"Using pasted sequence | length: {len(seq)}"

    elif mode == "Upload FASTA":
        file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"], key="reg_fasta_file")
        if file is not None:
            content = file.read().decode("utf-8", errors="ignore")
            seq = _parse_fasta_text(content)
            if seq:
                st.session_state["reg_seq"] = seq
                st.session_state["reg_status"] = f"Loaded FASTA | length: {len(seq)}"

    else:  # NCBI accession
        acc = st.text_input("Enter NCBI accession", key="reg_ncbi_acc")
        if st.button("Fetch from NCBI", key="reg_fetch_btn"):
            try:
                fasta = _fetch_ncbi_fasta(acc)
                seq = _parse_fasta_text(fasta)
                if not seq:
                    raise Exception("NCBI returned empty sequence.")
                st.session_state["reg_seq"] = seq
                st.session_state["reg_status"] = f"Fetched {acc} | length: {len(seq)}"
                st.rerun()
            except Exception as e:
                st.session_state["reg_status"] = ""
                st.error(str(e))

        # No big template box here
        if st.session_state["reg_seq"]:
            with st.expander("Preview fetched sequence", expanded=False):
                st.text_area("Fetched sequence (read-only)", st.session_state["reg_seq"], height=180, disabled=True)

    if st.session_state["reg_status"]:
        st.success(st.session_state["reg_status"])

    seq_used = _parse_fasta_text(st.session_state["reg_seq"])
    if not seq_used:
        st.info("Paste a sequence, upload a FASTA file, or fetch an NCBI accession.")
        add_footer()
        return
st.markdown("---")

    run = st.button("Design primers", key="reg_run_btn")
    if not run:
        add_footer()
        return

    try:
        fwd, rev, amp_len = design_basic_pcr_primers(
            seq_used,
            min_len=int(min_len),
            max_len=int(max_len),
            tm_target=float(tm_target),
            tm_tol=float(tm_tol),
            start_window=300,
            end_window=300,
            amplicon_min=int(amp_min),
            amplicon_max=int(amp_max),
        )

        st.success("Primers designed successfully.")
        st.caption(SCORE_EXPLANATION)

        rows = [
            {
                "Type": "FWD",
                "Primer (5'→3')": fwd,
                "Length": len(fwd),
                "Tm (°C)": round(tm_wallace(fwd), 1),
                "GC (%)": round(gc_pct(fwd), 1),
                "Score": round(primer_score(fwd, tm_target), 2),
            },
            {
                "Type": "REV",
                "Primer (5'→3')": rev,
        "Length": len(rev),
                "Tm (°C)": round(tm_wallace(rev), 1),
                "GC (%)": round(gc_pct(rev), 1),
                "Score": round(primer_score(rev, tm_target), 2),
            },
        ]
        st.dataframe(pd.DataFrame(rows), use_container_width=True)

        st.write(f"Estimated amplicon length: **{amp_len} bp**")

        st.subheader("Primer-BLAST link (pair)")
        url = primer_blast_url_pair(fwd, rev, "Homo sapiens")
        st.markdown(f"[Open in Primer-BLAST]({url})")
        st.info(BLAST_INSTRUCTIONS)

        st.subheader("Quick dimer screen")
        st.info(_quick_dimer_note(fwd, rev))

        add_footer()

    except Exception as e:
        st.error(str(e))
        add_footer()
