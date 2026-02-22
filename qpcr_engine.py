from dataclasses import dataclass
from typing import Optional, Tuple

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


def clean_seq(seq: str) -> str:
   s = (seq or "").upper()
   s = "".join(b for b in s if b in DNA)
   if not s:
       raise ValueError("No A/C/G/T bases found.")
   return s


def revcomp(seq: str) -> str:
   comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
   return "".join(comp[b] for b in reversed(seq))


def gc_content(seq: str) -> float:
   return 100.0 * sum(b in "GC" for b in seq) / len(seq)


def tm_wallace(seq: str) -> float:
   return 2.0 * (seq.count("A") + seq.count("T")) + 4.0 * (seq.count("G") + seq.count("C"))


def max_run(seq: str) -> int:
   best = 1
   cur = 1
   for i in range(1, len(seq)):
       if seq[i] == seq[i - 1]:
           cur += 1
           if cur > best:
               best = cur
       else:
           cur = 1
   return best


def self_complementarity_flag(seq: str) -> bool:
   rc = revcomp(seq)
   for i in range(len(seq) - 3):
       if seq[i:i + 4] in rc:
           return True
   return False


def has_3prime_complementarity(p1: str, p2: str, k: int = 4) -> bool:
   if len(p1) < k or len(p2) < k:
       return False
   return revcomp(p1[-k:]) in p2


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
) -> Optional[Tuple[float, float, float]]:
   tm = tm_wallace(seq)
   gc = gc_content(seq)

   if abs(tm - tm_target) > tm_tol:
       return None
   if not (gc_min <= gc <= gc_max):
       return None
   if max_run(seq) > max_homopolymer:
       return None
   if self_complementarity_flag(seq):
       return None

   score = 0.0
   score += abs(tm - tm_target) * 4.0
   if seq[-1] not in "GC":
       score += 1.5
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
) -> Optional[PrimerHit]:
   fwd_3p = fwd.start_0based + fwd.length
   search_start = fwd_3p + amplicon_min - 1
   search_end = min(len(template), fwd_3p + amplicon_max)
   if search_start >= search_end:
       return None

   best = None
   for L in range(min_len, max_len + 1):
       for i in range(search_start, search_end - L + 1):
           bind_site = template[i:i + L]
           rev_seq = revcomp(bind_site)

           scored = _qpcr_score_candidate(
               rev_seq, tm_target, tm_tol, gc_min, gc_max, max_homopolymer
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
) -> Optional[PrimerHit]:
   rev_bind_start = rev.start_0based
   search_end = max(0, rev_bind_start - amplicon_min + 1)
   search_start = max(0, rev_bind_start - amplicon_max)
   if search_start >= search_end:
       return None

   best = None
   for L in range(min_len, max_len + 1):
       for i in range(search_start, search_end - L + 1):
           fwd_seq = template[i:i + L]

           scored = _qpcr_score_candidate(
               fwd_seq, tm_target, tm_tol, gc_min, gc_max, max_homopolymer
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
           if self_complementarity_flag(pseq):
               continue

           tm = tm_wallace(pseq)
           gc = gc_content(pseq)
           if abs(tm - probe_tm_target) > probe_tm_tol:
               continue
           if not (30.0 <= gc <= 80.0):
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
   junction_max_3prime_distance: int = 4,
   max_tm_diff_pair: float = 1.0,
   dimer_k: int = 4,
   probe_tm_target: float = 69.0,
   probe_tm_tol: float = 3.0,
   probe_min_len: int = 18,
   probe_max_len: int = 30,
) -> Tuple[PrimerHit, PrimerHit, Optional[ProbeHit]]:

   left, right = parse_junction_marked_seq(seq_with_junction_marker, marker="^")
   template = left + right
   junction_index = len(left)

   def try_span(span: str):
       best = None  # (total_score, fwd, rev)
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
                       primer_gc_min, primer_gc_max, max_homopolymer
                   )
                   if scored is None:
                       continue
                   j_score, j_tm, j_gc = scored
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
                       max_homopolymer, amplicon_min, amplicon_max, dimer_k
                   )
                   if rev is None:
                       continue
                   if not _pair_ok_tm_diff(fwd.tm_c, rev.tm_c, max_tm_diff_pair):
                       continue

                   total = fwd.score + rev.score
                   if best is None or total < best[0]:
                       best = (total, fwd, rev)

               else:
                   if left_take > junction_max_3prime_distance:
                       continue

                   template_span = left[-left_take:] + right[:right_take]
                   rev_span = revcomp(template_span)
                   scored = _qpcr_score_candidate(
                       rev_span, primer_tm_target, primer_tm_tol,
                       primer_gc_min, primer_gc_max, max_homopolymer
                   )
                   if scored is None:
                       continue
                   j_score, j_tm, j_gc = scored
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
                       max_homopolymer, amplicon_min, amplicon_max, dimer_k
                   )
                   if fwd is None:
                       continue
                   if not _pair_ok_tm_diff(fwd.tm_c, rev.tm_c, max_tm_diff_pair):
                       continue

                   total = fwd.score + rev.score
                   if best is None or total < best[0]:
                       best = (total, fwd, rev)

       return best

   mode = (junction_primer or "AUTO").upper().strip()
   candidates = []

   if mode == "FWD":
       r = try_span("FWD")
       if r:
           candidates.append(r)
   elif mode == "REV":
       r = try_span("REV")
       if r:
           candidates.append(r)
   else:
       r1 = try_span("FWD")
       r2 = try_span("REV")
       if r1:
           candidates.append(r1)
       if r2:
           candidates.append(r2)

   if not candidates:
       raise RuntimeError("No qPCR junction primer pair found. Relax constraints.")

   candidates.sort(key=lambda x: x[0])
   _, best_fwd, best_rev = candidates[0]

   probe = None
   chemistry_mode = (chemistry or "SYBR").upper().replace(" ", "")
   if chemistry_mode in {"TAQMAN", "TAQMANPROBE"}:
       probe = design_taqman_probe_between(
           template,
           best_fwd,
           best_rev,
           probe_tm_target=probe_tm_target,
           probe_tm_tol=probe_tm_tol,
           probe_min_len=probe_min_len,
           probe_max_len=probe_max_len,
       )

   return best_fwd, best_rev, probe
