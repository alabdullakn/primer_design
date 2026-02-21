# designers/qpcr.py
import streamlit as st

from qpcr_engine import (
   amplicon_size_from_hits,
)

def render():
   st.header("qPCR primer designer")

   st.write("Paste your spliced sequence and mark the exon exon junction with ^")
   st.write("Example: ...ACCTG^GTTCA...")

   seq = st.text_area(
       "Spliced sequence with junction marker ^",
       height=180,
       placeholder="Paste DNA sequence here with ^ at the junction",
   )

   c1, c2, c3 = st.columns(3)
   with c1:
       chemistry = st.selectbox("Chemistry", ["SYBR", "TAQMAN"], index=0)
   with c2:
       junction_primer = st.selectbox("Junction spanning primer", ["AUTO", "FWD", "REV"], index=0)
   with c3:
       dimer_k = st.slider("3 prime complementarity check (k)", 3, 8, 4, 1)

   st.subheader("Primer constraints")
   r1, r2, r3, r4 = st.columns(4)
   with r1:
       min_len = st.slider("Min length", 16, 28, 18, 1)
   with r2:
       max_len = st.slider("Max length", 16, 32, 24, 1)
   with r3:
       tm_target = st.slider("Tm target", 50.0, 70.0, 60.0, 0.5)
   with r4:
       tm_tol = st.slider("Tm tolerance", 0.5, 10.0, 2.0, 0.5)

   g1, g2, g3 = st.columns(3)
   with g1:
       gc_min = st.slider("GC min percent", 20.0, 80.0, 40.0, 1.0)
   with g2:
       gc_max = st.slider("GC max percent", 20.0, 80.0, 60.0, 1.0)
   with g3:
       max_hpoly = st.slider("Max homopolymer run", 2, 6, 3, 1)

   st.subheader("Amplicon window")
   a1, a2, a3 = st.columns(3)
   with a1:
       amp_min = st.slider("Amplicon min bp", 40, 300, 70, 5)
   with a2:
       amp_max = st.slider("Amplicon max bp", 60, 500, 200, 5)
   with a3:
       max_tm_diff = st.slider("Max Tm difference pair", 0.0, 10.0, 1.0, 0.5)

   j_ov = st.slider("Min overlap each side of junction", 3, 12, 6, 1)

   st.divider()

   if st.button("Design qPCR primers"):
       try:
           if not seq or "^" not in seq:
               st.error("Your sequence must include the junction marker ^")
               return

           fwd, rev, probe = design_qpcr_junction_pair(
               seq_with_junction_marker=seq,
               chemistry=chemistry,
               junction_primer=junction_primer,
               min_len=min_len,
               max_len=max_len,
               primer_tm_target=tm_target,
               primer_tm_tol=tm_tol,
               primer_gc_min=gc_min,
               primer_gc_max=gc_max,
               max_homopolymer=max_hpoly,
               amplicon_min=amp_min,
               amplicon_max=amp_max,
               junction_min_overlap_each_side=j_ov,
               max_tm_diff_pair=max_tm_diff,
               dimer_k=dimer_k,
           )

           amp_size = amplicon_size_from_hits(fwd, rev)

           st.success("Designed primers")
           st.write(f"Amplicon size bp: {amp_size}")

           st.markdown("Forward primer 5 to 3")
           st.code(fwd.seq_5to3)
           st.write(f"Tm {fwd.tm_c:.1f} C, GC {fwd.gc_pct:.1f} percent")

           st.markdown("Reverse primer 5 to 3")
           st.code(rev.seq_5to3)
           st.write(f"Tm {rev.tm_c:.1f} C, GC {rev.gc_pct:.1f} percent")

           if probe is not None:
               st.markdown("Probe 5 to 3")
               st.code(probe.seq_5to3)
               st.write(f"Tm {probe.tm_c:.1f} C, GC {probe.gc_pct:.1f} percent")

       except Exception as e:
           st.error(f"qPCR design failed: {e}")
           st.info("Try widening Tm tolerance, widening GC range, or increasing amplicon max")
