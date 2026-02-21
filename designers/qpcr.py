# designers/qpcr.py
import streamlit as st

from qpcr_engine import (
   design_qpcr_junction_pair,
   amplicon_size_from_hits,
)

def render():
   st.title("qPCR primers")
   st.info("qPCR tab is running.")

   st.caption("Paste spliced sequence and mark exon-exon junction with '^'.")

   seq = st.text_area(
       "Sequence with '^' at junction",
       height=160,
       placeholder="...ACCTG^GTTACA...",
   )

   c1, c2, c3 = st.columns(3)
   with c1:
       chemistry = st.selectbox("Chemistry", ["SYBR", "TAQMAN"], index=0)
   with c2:
       junction_primer = st.selectbox("Junction primer", ["AUTO", "FWD", "REV"], index=0)
   with c3:
       dimer_k = st.number_input("3' dimer check (k)", min_value=3, max_value=8, value=4, step=1)

   st.subheader("Strict qPCR constraints")

   r1, r2, r3, r4 = st.columns(4)
   with r1:
       min_len = st.number_input("Min len", 14, 40, 18, 1)
   with r2:
       max_len = st.number_input("Max len", 14, 40, 24, 1)
   with r3:
       tm_target = st.number_input("Tm target", 45.0, 75.0, 60.0, 0.5)
   with r4:
       tm_tol = st.number_input("Tm tol", 0.5, 10.0, 2.0, 0.5)

   r5, r6, r7, r8 = st.columns(4)
   with r5:
       gc_min = st.number_input("GC min %", 0.0, 100.0, 40.0, 1.0)
   with r6:
       gc_max = st.number_input("GC max %", 0.0, 100.0, 60.0, 1.0)
   with r7:
       max_homo = st.number_input("Max homopolymer run", 2, 8, 3, 1)
   with r8:
       max_tm_diff = st.number_input("Max Tm diff", 0.0, 10.0, 1.0, 0.5)

   a1, a2, a3 = st.columns(3)
   with a1:
       amp_min = st.number_input("Amplicon min (bp)", 30, 1000, 70, 5)
   with a2:
       amp_max = st.number_input("Amplicon max (bp)", 30, 1000, 200, 5)
   with a3:
       overlap = st.number_input("Min overlap each side", 3, 12, 6, 1)

   if st.button("Design qPCR primers"):
       try:
           if "^" not in (seq or ""):
               st.error("You must include '^' in the sequence.")
               return

           fwd, rev, probe = design_qpcr_junction_pair(
               seq_with_junction_marker=seq,
               chemistry=chemistry,
               junction_primer=junction_primer,
               min_len=int(min_len),
               max_len=int(max_len),
               primer_tm_target=float(tm_target),
               primer_tm_tol=float(tm_tol),
               primer_gc_min=float(gc_min),
               primer_gc_max=float(gc_max),
               max_homopolymer=int(max_homo),
               amplicon_min=int(amp_min),
               amplicon_max=int(amp_max),
               junction_min_overlap_each_side=int(overlap),
               max_tm_diff_pair=float(max_tm_diff),
               dimer_k=int(dimer_k),
           )

           amp = amplicon_size_from_hits(fwd, rev)

           st.success(f"Designed primers. Amplicon: {amp} bp")

           st.markdown("### Forward primer (5'->3')")
           st.code(fwd.seq_5to3)
           st.write(f"Tm {fwd.tm_c:.1f} C, GC {fwd.gc_pct:.1f}%, role {fwd.role}")

           st.markdown("### Reverse primer (5'->3')")
           st.code(rev.seq_5to3)
           st.write(f"Tm {rev.tm_c:.1f} C, GC {rev.gc_pct:.1f}%, role {rev.role}")

           if chemistry == "TAQMAN":
               st.markdown("### Probe (5'->3')")
               if probe is None:
                   st.warning("No probe found.")
               else:
                   st.code(probe.seq_5to3)
                   st.write(f"Tm {probe.tm_c:.1f} C, GC {probe.gc_pct:.1f}%")

       except Exception as e:
           st.error("qPCR failed. See full error below.")
           st.exception(e)
