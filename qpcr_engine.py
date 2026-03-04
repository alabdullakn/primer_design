from dataclasses import dataclass
from typing import Optional, Tuple

from utils.tm import tm_nn, has_hairpin, has_homodimer
from utils.seq import clean_seq, revcomp, gc_content, max_run, has_3prime_complementarity

DNA = set("ACGT")


@dataclass
class PrimerHit:
   kind: str            # "FWD" or "REV"
   role: str            # "Junction" or "Partner"
   start_0based: int
   length: int
   seq_5to3: str
   tm_c: float
   gc_pct: float
   score: float
   template_bind_site_5to3: Optional[str] = None  # only meaningful for REV


@dataclass
class ProbeHit:
   start_0based: int
   length: int
   seq_5to3: str
   tm_c: float
   gc_pct: float
   score: float


def parse_junction_marked_seq(seq_with_marker: str, marker: str = "^") -> Tuple[str, str]:
   s = seq_with_marker or ""
   if marker not in s:
       raise ValueError("Missing '^' junction marker in the sequence.")
   left_raw, right_raw = s.split(marker, 1)
   left = clean_seq(left_raw)
   right = clean_seq(right_raw)
   if len(left) < 10 or len(right) < 10:
       raise ValueError("Add more bases on each side of '^' (at least 10 each side).")
   return left, right


def _qpcr_score_candidate(
   seq: str,
   tm_target: float,
   tm_tol: float,
   gc_min: float,
   gc_max: float,
   max_homopolymer: int,
   mv_conc: float = 50.0,
   dv_conc: float = 1.5,
   dntp_conc: float = 0.2,
) -> Optional[Tuple[float, float, float]]:
   tm = tm_nn(seq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc)
   gc = gc_content(seq)

   if abs(tm - tm_target) > tm_tol:
       return None
   if not (gc_min <= gc <= gc_max):
       return None
   if max_run(seq) > max_homopolymer:
       return None
   if has_hairpin(seq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc):
       return None
   if has_homodimer(seq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc):
       return None

   score = 0.0
   score += abs(tm - tm_target) * 4.0
   if seq[-1] not in "GC" and seq[-2] not in "GC":
       score += 1.5  # penalise lack of 3' GC clamp
   score += max(0, 20 - len(seq)) * 0.2  # prefer longer primers
   return score, tm, gc


def _pair_ok_tm_diff(fwd_tm: float, rev_tm: float, max_diff: float) -> bool:
   return abs(fwd_tm - rev_tm) <= max_diff


def amplicon_size_from_hits(fwd: PrimerHit, rev: PrimerHit) -> int:
   return max(0, (rev.start_0based + rev.length) - fwd.start_0based)


def _find_partner_reverse(
   template: str,
   fwd: PrimerHit,
   min_len: int,
   max_len: int,
   tm_target: float,
   tm_tol: float,
   gc_min: float,
   gc_max: float,
   max_homopolymer: int,
   amplicon_min: int,
   amplicon_max: int,
   dimer_k: int,
   mv_conc: float = 50.0,
   dv_conc: float = 1.5,
   dntp_conc: float = 0.2,
) -> Optional[PrimerHit]:

   best = None
   for L in range(min_len, max_len + 1):
       # amplicon_size = (rev_start + rev_len) - fwd_start
       i_min = fwd.start_0based + amplicon_min - L
       i_max = fwd.start_0based + amplicon_max - L
       search_start = max(0, i_min)
       search_end = min(len(template) - L, i_max)
       if search_start > search_end:
           continue

       for i in range(search_start, search_end + 1):
           bind_site = template[i:i + L]
           rev_seq = revcomp(bind_site)

           scored = _qpcr_score_candidate(
               rev_seq, tm_target, tm_tol, gc_min, gc_max, max_homopolymer,
               mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc,
           )
           if scored is None:
               continue

           if has_3prime_complementarity(fwd.seq_5to3, rev_seq, k=dimer_k):
               continue
           if has_3prime_complementarity(rev_seq, fwd.seq_5to3, k=dimer_k):
               continue

           score, tm, gc = scored
           hit = PrimerHit(
               kind="REV",
               role="Partner",
               start_0based=i,
               length=L,
               seq_5to3=rev_seq,
               tm_c=tm,
               gc_pct=gc,
               score=score,
               template_bind_site_5to3=bind_site,
           )
           if best is None or hit.score < best.score:
               best = hit
   return best


def _find_partner_forward(
   template: str,
   rev: PrimerHit,
   min_len: int,
   max_len: int,
   tm_target: float,
   tm_tol: float,
   gc_min: float,
   gc_max: float,
   max_homopolymer: int,
   amplicon_min: int,
   amplicon_max: int,
   dimer_k: int,
   mv_conc: float = 50.0,
   dv_conc: float = 1.5,
   dntp_conc: float = 0.2,
) -> Optional[PrimerHit]:
   best = None
   for L in range(min_len, max_len + 1):
       # amplicon_size = (rev_start + rev_len) - fwd_start
       i_min = rev.start_0based + rev.length - amplicon_max
       i_max = rev.start_0based + rev.length - amplicon_min
       search_start = max(0, i_min)
       search_end = min(len(template) - L, i_max)
       if search_start > search_end:
           continue

       for i in range(search_start, search_end + 1):
           fwd_seq = template[i:i + L]

           scored = _qpcr_score_candidate(
               fwd_seq, tm_target, tm_tol, gc_min, gc_max, max_homopolymer,
               mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc,
           )
           if scored is None:
               continue

           if has_3prime_complementarity(fwd_seq, rev.seq_5to3, k=dimer_k):
               continue
           if has_3prime_complementarity(rev.seq_5to3, fwd_seq, k=dimer_k):
               continue

           score, tm, gc = scored
           hit = PrimerHit(
               kind="FWD",
               role="Partner",
               start_0based=i,
               length=L,
               seq_5to3=fwd_seq,
               tm_c=tm,
               gc_pct=gc,
               score=score,
           )
           if best is None or hit.score < best.score:
               best = hit
   return best


def design_taqman_probe_between(
   template_spliced: str,
   fwd: PrimerHit,
   rev: PrimerHit,
   probe_tm_target: float = 69.0,
   probe_tm_tol: float = 3.0,
   probe_min_len: int = 18,
   probe_max_len: int = 30,
   gc_min: float = 30.0,
   gc_max: float = 80.0,
   mv_conc: float = 50.0,
   dv_conc: float = 1.5,
   dntp_conc: float = 0.2,
) -> ProbeHit:
   fwd_end = fwd.start_0based + fwd.length
   rev_start = rev.start_0based

   if rev_start - fwd_end < probe_min_len:
       raise RuntimeError("Amplicon too short to place a probe between primers.")

   region = template_spliced[fwd_end:rev_start]
   best = None

   for L in range(probe_min_len, probe_max_len + 1):
       for i in range(0, len(region) - L + 1):
           pseq = region[i:i + L]

           if pseq.startswith("G"):
               continue
           if max_run(pseq) > 3:
               continue
           if has_hairpin(pseq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc):
               continue
           if has_homodimer(pseq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc):
               continue

           tm = tm_nn(pseq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc)
           gc = gc_content(pseq)
           if abs(tm - probe_tm_target) > probe_tm_tol:
               continue
           if not (gc_min <= gc <= gc_max):
               continue

           score = abs(tm - probe_tm_target) * 4.0

           hit = ProbeHit(
               start_0based=fwd_end + i,
               length=L,
               seq_5to3=pseq,
               tm_c=tm,
               gc_pct=gc,
               score=score,
           )
           if best is None or hit.score < best.score:
               best = hit

   if best is None:
       raise RuntimeError("No valid probe found between primers.")
   return best


def design_qpcr_junction_pair(
   seq_with_junction_marker: str,
   chemistry: str = "SYBR",        # "SYBR" or "TAQMAN"
   junction_primer: str = "AUTO",  # "AUTO", "FWD", "REV"
   min_len: int = 18,
   max_len: int = 24,
   primer_tm_target: float = 60.0,
   primer_tm_tol: float = 2.0,
   primer_gc_min: float = 40.0,
   primer_gc_max: float = 60.0,
   max_homopolymer: int = 3,
   amplicon_min: int = 70,
   amplicon_max: int = 200,
   junction_min_overlap_each_side: int = 6,
   junction_max_3prime_distance: int = 8,
   max_tm_diff_pair: float = 1.0,
   dimer_k: int = 4,
   probe_tm_target: float = 69.0,
   probe_tm_tol: float = 3.0,
   probe_min_len: int = 18,
   probe_max_len: int = 30,
   mv_conc: float = 50.0,
   dv_conc: float = 1.5,
   dntp_conc: float = 0.2,
   top_n: int = 5,
) -> List[Tuple[PrimerHit, PrimerHit, Optional[ProbeHit]]]:

   left, right = parse_junction_marked_seq(seq_with_junction_marker, marker="^")
   template = left + right
   junction_index = len(left)

   def try_span(span: str):
       collected = []  # list of (total_score, fwd, rev)
       for L in range(min_len, max_len + 1):
           for left_take in range(junction_min_overlap_each_side, L - junction_min_overlap_each_side + 1):
               right_take = L - left_take
               if right_take < junction_min_overlap_each_side:
                   continue

               if span == "FWD":
                   if right_take > junction_max_3prime_distance:
                       continue

                   jseq = left[-left_take:] + right[:right_take]
                   scored = _qpcr_score_candidate(
                       jseq, primer_tm_target, primer_tm_tol,
                       primer_gc_min, primer_gc_max, max_homopolymer,
                       mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc,
                   )
                   if scored is None:
                       continue
                   j_score, j_tm, j_gc = scored
                   # Reward junction closer to 3' end (right_take=1 is ideal)
                   j_score += (right_take - 1) * 0.5
                   fwd = PrimerHit(
                       kind="FWD",
                       role="Junction",
                       start_0based=junction_index - left_take,
                       length=L,
                       seq_5to3=jseq,
                       tm_c=j_tm,
                       gc_pct=j_gc,
                       score=j_score,
                   )

                   rev = _find_partner_reverse(
                       template, fwd, min_len, max_len,
                       primer_tm_target, primer_tm_tol,
                       primer_gc_min, primer_gc_max,
                       max_homopolymer, amplicon_min, amplicon_max, dimer_k,
                       mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc,
                   )
                   if rev is None:
                       continue
                   if not _pair_ok_tm_diff(fwd.tm_c, rev.tm_c, max_tm_diff_pair):
                       continue

                   total = fwd.score + rev.score
                   collected.append((total, fwd, rev))

               else:
                   if left_take > junction_max_3prime_distance:
                       continue

                   template_span = left[-left_take:] + right[:right_take]
                   rev_span = revcomp(template_span)
                   scored = _qpcr_score_candidate(
                       rev_span, primer_tm_target, primer_tm_tol,
                       primer_gc_min, primer_gc_max, max_homopolymer,
                       mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc,
                   )
                   if scored is None:
                       continue
                   j_score, j_tm, j_gc = scored
                   # For REV, left_take is the distance from junction to 3' end
                   j_score += (left_take - 1) * 0.5
                   rev = PrimerHit(
                       kind="REV",
                       role="Junction",
                       start_0based=junction_index - left_take,
                       length=L,
                       seq_5to3=rev_span,
                       tm_c=j_tm,
                       gc_pct=j_gc,
                       score=j_score,
                       template_bind_site_5to3=template_span,
                   )

                   fwd = _find_partner_forward(
                       template, rev, min_len, max_len,
                       primer_tm_target, primer_tm_tol,
                       primer_gc_min, primer_gc_max,
                       max_homopolymer, amplicon_min, amplicon_max, dimer_k,
                       mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc,
                   )
                   if fwd is None:
                       continue
                   if not _pair_ok_tm_diff(fwd.tm_c, rev.tm_c, max_tm_diff_pair):
                       continue

                   total = fwd.score + rev.score
                   collected.append((total, fwd, rev))

       return collected

   mode = (junction_primer or "AUTO").upper().strip()
   all_candidates = []

   if mode == "FWD":
       all_candidates.extend(try_span("FWD"))
   elif mode == "REV":
       all_candidates.extend(try_span("REV"))
   else:
       all_candidates.extend(try_span("FWD"))
       all_candidates.extend(try_span("REV"))

   if not all_candidates:
       raise RuntimeError("No qPCR junction primer pair found. Relax constraints.")

   all_candidates.sort(key=lambda x: x[0])
   chemistry_mode = (chemistry or "SYBR").upper().replace(" ", "")

   results = []
   for _, fwd, rev in all_candidates[:top_n]:
       probe = None
       if chemistry_mode in {"TAQMAN", "TAQMANPROBE"}:
           try:
               probe = design_taqman_probe_between(
                   template, fwd, rev,
                   probe_tm_target=probe_tm_target,
                   probe_tm_tol=probe_tm_tol,
                   probe_min_len=probe_min_len,
                   probe_max_len=probe_max_len,
                   gc_min=primer_gc_min,
                   gc_max=primer_gc_max,
                   mv_conc=mv_conc,
                   dv_conc=dv_conc,
                   dntp_conc=dntp_conc,
               )
           except RuntimeError:
               probe = None
       results.append((fwd, rev, probe))

   return results


def design_intron_flanking_pair(
   exon_a: str,
   exon_b: str,
   min_len: int = 18,
   max_len: int = 24,
   tm_target: float = 60.0,
   tm_tol: float = 2.0,
   gc_min: float = 40.0,
   gc_max: float = 60.0,
   max_homopolymer: int = 3,
   amplicon_min: int = 70,
   amplicon_max: int = 200,
   dimer_k: int = 4,
   max_tm_diff_pair: float = 1.0,
   mv_conc: float = 50.0,
   dv_conc: float = 1.5,
   dntp_conc: float = 0.2,
   top_n: int = 5,
) -> List[Tuple[PrimerHit, PrimerHit]]:
   """
   Design a qPCR primer pair where FWD is entirely within exon_a and REV
   is entirely within exon_b.  The cDNA amplicon spans the exon junction;
   the genomic amplicon is much larger (intron between the two exons).
   Returns (fwd_hit, rev_hit).
   """
   exon_a = clean_seq(exon_a)
   exon_b = clean_seq(exon_b)
   if len(exon_a) < min_len:
       raise RuntimeError("Exon A sequence is shorter than the minimum primer length.")
   if len(exon_b) < min_len:
       raise RuntimeError("Exon B sequence is shorter than the minimum primer length.")

   junction = len(exon_a)
   template = exon_a + exon_b
   n = len(template)

   all_pairs = []  # (total_score, fwd_hit, rev_hit)

   for Lf in range(min_len, max_len + 1):
       # FWD primer must sit entirely within exon_a
       for i in range(0, junction - Lf + 1):
           fwd_seq = template[i:i + Lf]
           f_scored = _qpcr_score_candidate(
               fwd_seq, tm_target, tm_tol, gc_min, gc_max, max_homopolymer,
               mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc,
           )
           if f_scored is None:
               continue
           f_score, f_tm, f_gc = f_scored
           fwd_3 = i + Lf  # 3' end of FWD on the template

           for Lr in range(min_len, max_len + 1):
               # REV primer must sit entirely within exon_b: j >= junction
               j_min = max(junction, fwd_3 + amplicon_min - Lr)
               j_max = min(n - Lr, fwd_3 + amplicon_max - Lr)
               if j_min > j_max:
                   continue

               for j in range(j_min, j_max + 1):
                   bind_site = template[j:j + Lr]
                   rev_seq = revcomp(bind_site)

                   r_scored = _qpcr_score_candidate(
                       rev_seq, tm_target, tm_tol, gc_min, gc_max, max_homopolymer,
                       mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc,
                   )
                   if r_scored is None:
                       continue
                   r_score, r_tm, r_gc = r_scored

                   if not _pair_ok_tm_diff(f_tm, r_tm, max_tm_diff_pair):
                       continue
                   if has_3prime_complementarity(fwd_seq, rev_seq, k=dimer_k):
                       continue
                   if has_3prime_complementarity(rev_seq, fwd_seq, k=dimer_k):
                       continue

                   amp_len = (j + Lr) - i
                   if amp_len < amplicon_min or amp_len > amplicon_max:
                       continue

                   total = f_score + r_score
                   fwd_hit = PrimerHit(
                       kind="FWD", role="Exon A",
                       start_0based=i, length=Lf, seq_5to3=fwd_seq,
                       tm_c=f_tm, gc_pct=f_gc, score=f_score,
                   )
                   rev_hit = PrimerHit(
                       kind="REV", role="Exon B",
                       start_0based=j, length=Lr, seq_5to3=rev_seq,
                       tm_c=r_tm, gc_pct=r_gc, score=r_score,
                       template_bind_site_5to3=bind_site,
                   )
                   all_pairs.append((total, fwd_hit, rev_hit))

   if not all_pairs:
       raise RuntimeError(
           "No intron-flanking primer pair found. "
           "Try widening Tm tolerance, GC range, amplicon range, or provide longer exon sequences."
       )

   all_pairs.sort(key=lambda x: x[0])
   return [(fwd, rev) for _, fwd, rev in all_pairs[:top_n]]
