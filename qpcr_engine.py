 (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a/qpcr_engine.py b/qpcr_engine.py
index 851554e2e489d7c55f832eeaad825e0bebd3091e..9085496e12bbf5088ffb5a77b6c37fb3aabbd239 100644
--- a/qpcr_engine.py
+++ b/qpcr_engine.py
@@ -274,102 +274,113 @@ def design_taqman_probe_between(
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
+   junction_max_3prime_distance: int = 4,
    max_tm_diff_pair: float = 1.0,
    dimer_k: int = 4,
+   probe_tm_target: float = 69.0,
+   probe_tm_tol: float = 3.0,
+   probe_min_len: int = 18,
+   probe_max_len: int = 30,
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
+                   if right_take > junction_max_3prime_distance:
+                       continue
+
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
+                   if left_take > junction_max_3prime_distance:
+                       continue
+
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
@@ -390,29 +401,38 @@ def design_qpcr_junction_pair(
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
-   if (chemistry or "SYBR").upper().strip() == "TAQMAN":
-       probe = design_taqman_probe_between(template, best_fwd, best_rev)
+   chemistry_mode = (chemistry or "SYBR").upper().replace(" ", "")
+   if chemistry_mode in {"TAQMAN", "TAQMANPROBE"}:
+       probe = design_taqman_probe_between(
+           template,
+           best_fwd,
+           best_rev,
+           probe_tm_target=probe_tm_target,
+           probe_tm_tol=probe_tm_tol,
+           probe_min_len=probe_min_len,
+           probe_max_len=probe_max_len,
+       )
 
    return best_fwd, best_rev, probe
 
EOF
)
