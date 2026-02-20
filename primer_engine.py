# primer_engine.py
# Exon primer helper (2 forward + 1 reverse) with:
# - length range (default 16..25)
# - Tm filter (Wallace by default)
# - basic quality filters
# - 3' complementarity (dimer) checks between forward and reverse
# - optional NCBI BLAST summary (Biopython)

from dataclasses import dataclass
from typing import List, Tuple, Optional

DNA = set("ACGT")


# ---------------- Data structures ----------------

@dataclass
class BlastSummary:
    ok: bool
    hits_returned: Optional[int] = None
    top_hit_id: Optional[str] = None
    top_hit_def: Optional[str] = None
    top_evalue: Optional[float] = None
    top_identity_pct: Optional[float] = None
    note: Optional[str] = None


@dataclass
class PrimerHit:
    kind: str                      # "FWD" or "REV"
    exon_name: str                 # "Exon1", "Exon2", "Exon3"
    start_0based: int
    length: int
    seq_5to3: str                  # primer you ORDER (5'->3')
    tm_c: float
    gc_pct: float
    score: float
    template_seq_5to3: Optional[str] = None  # only for REV: binding site in exon (as given)
    blast: Optional[BlastSummary] = None


# ---------------- Basic helpers ----------------

def clean_seq(seq: str) -> str:
    """
    Keep only A/C/G/T and uppercase.
    This lets users paste FASTA headers or spaces without crashing.
    """
    seq = seq.upper()
    seq = "".join(b for b in seq if b in DNA)
    if not seq:
        raise ValueError("No A/C/G/T bases found. Paste a DNA exon sequence.")
    return seq


def revcomp(seq: str) -> str:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp[b] for b in reversed(seq))


def gc_content(seq: str) -> float:
    return 100.0 * sum(b in "GC" for b in seq) / len(seq)


def tm_wallace(seq: str) -> float:
    # Simple and fine for quick ranking; OneTaq will tolerate small mismatch anyway.
    return 2.0 * (seq.count("A") + seq.count("T")) + 4.0 * (seq.count("G") + seq.count("C"))


def max_run(seq: str) -> int:
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def has_3prime_gc_clamp(seq: str) -> bool:
    # soft preference only
    if len(seq) < 2:
        return False
    return (seq[-1] in "GC") or (seq[-2] in "GC")


def self_complementarity_flag(seq: str) -> bool:
    """
    Quick heuristic: flags if any 4-mer appears in its own reverse complement.
    Not a full thermodynamic hairpin check, but catches obvious cases.
    """
    rc = revcomp(seq)
    for i in range(len(seq) - 3):
        if seq[i:i + 4] in rc:
            return True
    return False


def has_3prime_complementarity(p1: str, p2: str, k: int = 4) -> bool:
    """
    True if last k bases of p1 (3' end) are complementary to any part of p2.
    This is a simple dimer screen.
    """
    if len(p1) < k:
        return False
    return revcomp(p1[-k:]) in p2


def dimer_risk_percent(p1: str, p2: str, max_k: int = 8) -> float:
    """
    Heuristic "risk %" based on the longest 3' complementarity length found.
    Returns 0..100. This is not a thermodynamic probability.
    """
    best = 0
    for k in range(3, max_k + 1):
        if has_3prime_complementarity(p1, p2, k=k):
            best = k
    # map k=3..8 to 20..100 roughly
    if best == 0:
        return 0.0
    return min(100.0, (best - 2) / (max_k - 2) * 100.0)


# ---------------- Candidate scoring ----------------

def score_candidate(seq: str, tm_target: float, tm_tol: float) -> Optional[Tuple[float, float, float]]:
    tm = tm_wallace(seq)
    gc = gc_content(seq)

    # hard filters
    if abs(tm - tm_target) > tm_tol:
        return None
    if not (35.0 <= gc <= 65.0):
        return None
    if max_run(seq) >= 5:
        return None
    if self_complementarity_flag(seq):
        return None

    # soft score (lower is better)
    score = 0.0
    score += abs(tm - tm_target) * 3.0
    if not has_3prime_gc_clamp(seq):
        score += 2.0
    if max_run(seq) == 4:
        score += 2.0

    return score, tm, gc


def best_primer_from_exon(
    exon_seq: str,
    exon_name: str,
    kind: str,
    min_len: int,
    max_len: int,
    tm_target: float,
    tm_tol: float
) -> PrimerHit:
    exon_seq = clean_seq(exon_seq)
    best: Optional[PrimerHit] = None

    for L in range(min_len, max_len + 1):
        for i in range(0, len(exon_seq) - L + 1):
            window = exon_seq[i:i + L]
            primer = window if kind == "FWD" else revcomp(window)

            scored = score_candidate(primer, tm_target, tm_tol)
            if scored is None:
                continue

            score, tm, gc = scored
            hit = PrimerHit(kind=kind, exon_name=exon_name, start_0based=i,
                            length=L, seq_5to3=primer, tm_c=tm, gc_pct=gc, score=score)

            if best is None or hit.score < best.score:
                best = hit

    if best is None:
        raise RuntimeError(f"No {kind} primer found for {exon_name} with your constraints.")
    return best


def best_reverse_avoiding_fwds(
    exon_seq_for_reverse: str,
    fwd1: PrimerHit,
    fwd2: PrimerHit,
    min_len: int,
    max_len: int,
    tm_target: float,
    tm_tol: float,
    dimer_k: int = 4
) -> PrimerHit:
    exon_seq = clean_seq(exon_seq_for_reverse)
    candidates: List[PrimerHit] = []

    for L in range(min_len, max_len + 1):
        for i in range(0, len(exon_seq) - L + 1):
            template_seq = exon_seq[i:i + L]     # as provided (5'->3' exon text)
            rev_primer = revcomp(template_seq)   # primer you order (5'->3')

            scored = score_candidate(rev_primer, tm_target, tm_tol)
            if scored is None:
                continue

            # avoid 3' complementarity to either forward primer
            if has_3prime_complementarity(fwd1.seq_5to3, rev_primer, k=dimer_k):
                continue
            if has_3prime_complementarity(fwd2.seq_5to3, rev_primer, k=dimer_k):
                continue

            score, tm, gc = scored
            candidates.append(
                PrimerHit(
                    kind="REV",
                    exon_name="Exon3",
                    start_0based=i,
                    length=L,
                    seq_5to3=rev_primer,
                    tm_c=tm,
                    gc_pct=gc,
                    score=score,
                    template_seq_5to3=template_seq
                )
            )

    if not candidates:
        raise RuntimeError(
            "No reverse primer found that avoids 3' complementarity to both forward primers. "
            "Try dimer_k=3, widen tm_tol, or widen length range."
        )

    candidates.sort(key=lambda x: x.score)
    return candidates[0]


def design_exon_primers(
    exon1: str,
    exon2: str,
    exon3_for_reverse: str,
    min_len: int = 16,
    max_len: int = 25,
    tm_target: float = 60.0,
    tm_tol: float = 5.0,
    dimer_k: int = 4
) -> Tuple[PrimerHit, PrimerHit, PrimerHit]:
    """
    Returns: (fwd_exon1, fwd_exon2, rev_exon3)
    """
    fwd1 = best_primer_from_exon(exon1, "Exon1", "FWD", min_len, max_len, tm_target, tm_tol)
    fwd2 = best_primer_from_exon(exon2, "Exon2", "FWD", min_len, max_len, tm_target, tm_tol)
    rev3 = best_reverse_avoiding_fwds(exon3_for_reverse, fwd1, fwd2, min_len, max_len, tm_target, tm_tol, dimer_k=dimer_k)
    return fwd1, fwd2, rev3


# ---------------- BLAST (optional) ----------------

def blast_primer_summary_ncbi(
    primer_seq: str,
    database: str = "refseq_rna",
    organism_filter: str = "Homo sapiens[Organism]",
    hitlist_size: int = 10
) -> BlastSummary:
    """
    Remote BLAST via NCBI using Biopython qblast.
    This can fail on hosted platforms (rate limits / blocked outbound).
    """
    try:
        from Bio.Blast import NCBIWWW, NCBIXML
    except Exception as e:
        return BlastSummary(ok=False, note=f"Biopython not available: {e}")

    try:
        handle = NCBIWWW.qblast(
            program="blastn",
            database=database,
            sequence=primer_seq,
            entrez_query=organism_filter,
            hitlist_size=hitlist_size,
            short_query=True,
            format_type="XML"
        )
        record = NCBIXML.read(handle)
        handle.close()

        hits = len(record.alignments)
        if hits == 0:
            return BlastSummary(ok=True, hits_returned=0, note="No hits returned (with current filters).")

        top = record.alignments[0]
        top_hsp = top.hsps[0]
        ident = 100.0 * top_hsp.identities / max(1, top_hsp.align_length)

        return BlastSummary(
            ok=True,
            hits_returned=hits,
            top_hit_id=top.hit_id,
            top_hit_def=top.hit_def,
            top_evalue=float(top_hsp.expect),
            top_identity_pct=float(ident),
            note=f"db={database}, filter={organism_filter}"
        )
    except Exception as e:
        return BlastSummary(ok=False, note=f"BLAST failed: {e}")


def attach_blast(
    primers: List[PrimerHit],
    do_blast: bool = False,
    database: str = "refseq_rna",
    organism_filter: str = "Homo sapiens[Organism]",
    hitlist_size: int = 10
) -> None:
    if not do_blast:
        return
    for p in primers:
        p.blast = blast_primer_summary_ncbi(
            p.seq_5to3,
            database=database,
            organism_filter=organism_filter,
            hitlist_size=hitlist_size
        )


# ---------------- Streamlit dimer report ----------------

def print_dimer_report(f1: PrimerHit, f2: PrimerHit, r3: PrimerHit) -> None:
    """
    Uses Streamlit directly so app.py can call it and render nicely.
    """
    try:
        import streamlit as st
    except Exception:
        return

    rows = []
    for f in [f1, f2]:
        risk = dimer_risk_percent(f.seq_5to3, r3.seq_5to3, max_k=8)
        rows.append({
            "Forward": f"{f.kind} {f.exon_name}",
            "Reverse": f"{r3.kind} {r3.exon_name}",
            "3' complementarity (k)": ">=4 blocked by design",
            "Heuristic dimer risk (%)": round(risk, 1)
        })

    st.table(rows)





### qpcr 

from dataclasses import dataclass
from typing import Optional, Tuple, List

# Reuse your existing helpers if they already exist in primer_engine.py
# If you already have these functions, do NOT duplicate them.
# clean_seq, revcomp, gc_content, tm_wallace, max_run, has_3prime_gc_clamp, self_complementarity_flag, score_candidate
# and PrimerHit dataclass

def parse_junction_marked_seq(seq_with_marker: str, marker: str = "^") -> Tuple[str, str]:
    s = clean_seq(seq_with_marker.replace(marker, ""))
    if marker not in seq_with_marker:
        raise ValueError("For qPCR, mark the exon-exon junction using '^' inside the sequence.")
    left_raw, right_raw = seq_with_marker.split(marker, 1)
    left = clean_seq(left_raw)
    right = clean_seq(right_raw)
    if len(left) < 10 or len(right) < 10:
        raise ValueError("Junction sides are too short. Paste more bases on each side of '^'.")
    return left, right

def build_junction_primer_forward(left: str, right: str, left_take: int, right_take: int) -> str:
    # Primer is written 5'->3' as ordered, crossing junction
    return left[-left_take:] + right[:right_take]

def build_junction_primer_reverse(left: str, right: str, left_take: int, right_take: int) -> str:
    # Reverse primer spans junction on the minus strand
    template = left[-left_take:] + right[:right_take]
    return revcomp(template)

def design_qpcr_junction_primers(
    seq_with_junction_marker: str,
    span_primer: str = "FWD",           # "FWD" or "REV" spans the junction
    min_len: int = 18,
    max_len: int = 25,
    tm_target: float = 60.0,
    tm_tol: float = 5.0,
    amplicon_min: int = 70,
    amplicon_max: int = 200,
    junction_min_overlap_each_side: int = 6
) -> Tuple["PrimerHit", "PrimerHit"]:
    """
    Returns (forward, reverse) primers.
    One primer is forced to span the junction.
    The partner primer is chosen to yield a short amplicon on the spliced template.
    """

    left, right = parse_junction_marked_seq(seq_with_junction_marker, marker="^")

    # Build full spliced template used for partner search
    template = left + right
    junction_index = len(left)  # 0-based index in template where right starts

    best_pair = None  # (score_sum, fwd_hit, rev_hit)

    # Enumerate junction spanning primer lengths by choosing how many bases from each side
    for L in range(min_len, max_len + 1):
        # need at least junction_min_overlap_each_side from each side
        for left_take in range(junction_min_overlap_each_side, L - junction_min_overlap_each_side + 1):
            right_take = L - left_take
            if right_take < junction_min_overlap_each_side:
                continue

            if span_primer.upper() == "FWD":
                jseq = build_junction_primer_forward(left, right, left_take, right_take)
                j_scored = score_candidate(jseq, tm_target, tm_tol)
                if j_scored is None:
                    continue
                j_score, j_tm, j_gc = j_scored
                fwd_hit = PrimerHit(
                    kind="FWD",
                    exon_name="Junction",
                    start_0based=junction_index - left_take,
                    length=L,
                    seq_5to3=jseq,
                    tm_c=j_tm,
                    gc_pct=j_gc,
                    score=j_score
                )

                # Partner reverse primer should bind on the right side to make small amplicon
                # We search reverse primer binding sites on the template within amplicon window
                # Reverse primer binds to template region downstream of forward 3' end
                fwd_3prime_pos = fwd_hit.start_0based + fwd_hit.length  # first base after primer on template
                search_start = fwd_3prime_pos + amplicon_min - 1
                search_end = min(len(template), fwd_3prime_pos + amplicon_max)

                if search_start >= search_end:
                    continue

                best_rev = None
                for Lr in range(min_len, max_len + 1):
                    for i in range(search_start, search_end - Lr + 1):
                        bind_site = template[i:i + Lr]
                        rev_seq = revcomp(bind_site)
                        r_scored = score_candidate(rev_seq, tm_target, tm_tol)
                        if r_scored is None:
                            continue
                        r_score, r_tm, r_gc = r_scored
                        rev_hit = PrimerHit(
                            kind="REV",
                            exon_name="RightSide",
                            start_0based=i,
                            length=Lr,
                            seq_5to3=rev_seq,
                            tm_c=r_tm,
                            gc_pct=r_gc,
                            score=r_score,
                            template_seq_5to3=bind_site
                        )
                        if best_rev is None or rev_hit.score < best_rev.score:
                            best_rev = rev_hit

                if best_rev is None:
                    continue

                total = fwd_hit.score + best_rev.score
                if best_pair is None or total < best_pair[0]:
                    best_pair = (total, fwd_hit, best_rev)

            else:
                # span_primer == "REV": reverse primer spans junction, forward primer upstream on left side
                jseq = build_junction_primer_reverse(left, right, left_take, right_take)
                j_scored = score_candidate(jseq, tm_target, tm_tol)
                if j_scored is None:
                    continue
                j_score, j_tm, j_gc = j_scored
                rev_hit = PrimerHit(
                    kind="REV",
                    exon_name="Junction",
                    start_0based=junction_index - left_take,
                    length=L,
                    seq_5to3=jseq,
                    tm_c=j_tm,
                    gc_pct=j_gc,
                    score=j_score,
                    template_seq_5to3=(left[-left_take:] + right[:right_take])
                )

                # Partner forward primer should be on the left side to make small amplicon
                # Forward primer binds upstream of reverse binding site.
                rev_bind_start = rev_hit.start_0based  # reverse binds here (binding site start on template)
                search_end = max(0, rev_bind_start - amplicon_min + 1)
                search_start = max(0, rev_bind_start - amplicon_max)

                if search_start >= search_end:
                    continue

                best_fwd = None
                for Lf in range(min_len, max_len + 1):
                    for i in range(search_start, search_end - Lf + 1):
                        fseq = template[i:i + Lf]
                        f_scored = score_candidate(fseq, tm_target, tm_tol)
                        if f_scored is None:
                            continue
                        f_score, f_tm, f_gc = f_scored
                        fwd_hit = PrimerHit(
                            kind="FWD",
                            exon_name="LeftSide",
                            start_0based=i,
                            length=Lf,
                            seq_5to3=fseq,
                            tm_c=f_tm,
                            gc_pct=f_gc,
                            score=f_score
                        )
                        if best_fwd is None or fwd_hit.score < best_fwd.score:
                            best_fwd = fwd_hit

                if best_fwd is None:
                    continue

                total = best_fwd.score + rev_hit.score
                if best_pair is None or total < best_pair[0]:
                    best_pair = (total, best_fwd, rev_hit)

    if best_pair is None:
        raise RuntimeError("No qPCR junction primer pair found with your constraints. Try wider Tm tolerance or length range.")

    _, best_fwd, best_rev = best_pair
    return best_fwd, best_rev

def qpcr_amplicon_size(template_seq_with_marker: str, fwd: "PrimerHit", rev: "PrimerHit") -> int:
    # Amplicon size on spliced template = from fwd start to end of rev binding site
    left, right = parse_junction_marked_seq(template_seq_with_marker, marker="^")
    template = left + right
    rev_end = rev.start_0based + rev.length
    return max(0, rev_end - fwd.start_0based)




# primer_engine.py
# Exon primer helper (2 forward + 1 reverse) with:
# - length range (default 16..25)
# - Tm filter (Wallace by default)
# - basic quality filters
# - 3' complementarity (dimer) checks between forward and reverse
# - optional NCBI BLAST summary (Biopython)

from dataclasses import dataclass
from typing import List, Tuple, Optional

DNA = set("ACGT")


# ---------------- Data structures ----------------

@dataclass
class BlastSummary:
    ok: bool
    hits_returned: Optional[int] = None
    top_hit_id: Optional[str] = None
    top_hit_def: Optional[str] = None
    top_evalue: Optional[float] = None
    top_identity_pct: Optional[float] = None
    note: Optional[str] = None


@dataclass
class PrimerHit:
    kind: str                      # "FWD" or "REV"
    exon_name: str                 # "Exon1", "Exon2", "Exon3"
    start_0based: int
    length: int
    seq_5to3: str                  # primer you ORDER (5'->3')
    tm_c: float
    gc_pct: float
    score: float
    template_seq_5to3: Optional[str] = None  # only for REV: binding site in exon (as given)
    blast: Optional[BlastSummary] = None


# ---------------- Basic helpers ----------------

def clean_seq(seq: str) -> str:
    """
    Keep only A/C/G/T and uppercase.
    This lets users paste FASTA headers or spaces without crashing.
    """
    seq = seq.upper()
    seq = "".join(b for b in seq if b in DNA)
    if not seq:
        raise ValueError("No A/C/G/T bases found. Paste a DNA exon sequence.")
    return seq


def revcomp(seq: str) -> str:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp[b] for b in reversed(seq))


def gc_content(seq: str) -> float:
    return 100.0 * sum(b in "GC" for b in seq) / len(seq)


def tm_wallace(seq: str) -> float:
    # Simple and fine for quick ranking; OneTaq will tolerate small mismatch anyway.
    return 2.0 * (seq.count("A") + seq.count("T")) + 4.0 * (seq.count("G") + seq.count("C"))


def max_run(seq: str) -> int:
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def has_3prime_gc_clamp(seq: str) -> bool:
    # soft preference only
    if len(seq) < 2:
        return False
    return (seq[-1] in "GC") or (seq[-2] in "GC")


def self_complementarity_flag(seq: str) -> bool:
    """
    Quick heuristic: flags if any 4-mer appears in its own reverse complement.
    Not a full thermodynamic hairpin check, but catches obvious cases.
    """
    rc = revcomp(seq)
    for i in range(len(seq) - 3):
        if seq[i:i + 4] in rc:
            return True
    return False


def has_3prime_complementarity(p1: str, p2: str, k: int = 4) -> bool:
    """
    True if last k bases of p1 (3' end) are complementary to any part of p2.
    This is a simple dimer screen.
    """
    if len(p1) < k:
        return False
    return revcomp(p1[-k:]) in p2


def dimer_risk_percent(p1: str, p2: str, max_k: int = 8) -> float:
    """
    Heuristic "risk %" based on the longest 3' complementarity length found.
    Returns 0..100. This is not a thermodynamic probability.
    """
    best = 0
    for k in range(3, max_k + 1):
        if has_3prime_complementarity(p1, p2, k=k):
            best = k
    # map k=3..8 to 20..100 roughly
    if best == 0:
        return 0.0
    return min(100.0, (best - 2) / (max_k - 2) * 100.0)


# ---------------- Candidate scoring ----------------

def score_candidate(seq: str, tm_target: float, tm_tol: float) -> Optional[Tuple[float, float, float]]:
    tm = tm_wallace(seq)
    gc = gc_content(seq)

    # hard filters
    if abs(tm - tm_target) > tm_tol:
        return None
    if not (35.0 <= gc <= 65.0):
        return None
    if max_run(seq) >= 5:
        return None
    if self_complementarity_flag(seq):
        return None

    # soft score (lower is better)
    score = 0.0
    score += abs(tm - tm_target) * 3.0
    if not has_3prime_gc_clamp(seq):
        score += 2.0
    if max_run(seq) == 4:
        score += 2.0

    return score, tm, gc


def best_primer_from_exon(
    exon_seq: str,
    exon_name: str,
    kind: str,
    min_len: int,
    max_len: int,
    tm_target: float,
    tm_tol: float
) -> PrimerHit:
    exon_seq = clean_seq(exon_seq)
    best: Optional[PrimerHit] = None

    for L in range(min_len, max_len + 1):
        for i in range(0, len(exon_seq) - L + 1):
            window = exon_seq[i:i + L]
            primer = window if kind == "FWD" else revcomp(window)

            scored = score_candidate(primer, tm_target, tm_tol)
            if scored is None:
                continue

            score, tm, gc = scored
            hit = PrimerHit(kind=kind, exon_name=exon_name, start_0based=i,
                            length=L, seq_5to3=primer, tm_c=tm, gc_pct=gc, score=score)

            if best is None or hit.score < best.score:
                best = hit

    if best is None:
        raise RuntimeError(f"No {kind} primer found for {exon_name} with your constraints.")
    return best


def best_reverse_avoiding_fwds(
    exon_seq_for_reverse: str,
    fwd1: PrimerHit,
    fwd2: PrimerHit,
    min_len: int,
    max_len: int,
    tm_target: float,
    tm_tol: float,
    dimer_k: int = 4
) -> PrimerHit:
    exon_seq = clean_seq(exon_seq_for_reverse)
    candidates: List[PrimerHit] = []

    for L in range(min_len, max_len + 1):
        for i in range(0, len(exon_seq) - L + 1):
            template_seq = exon_seq[i:i + L]     # as provided (5'->3' exon text)
            rev_primer = revcomp(template_seq)   # primer you order (5'->3')

            scored = score_candidate(rev_primer, tm_target, tm_tol)
            if scored is None:
                continue

            # avoid 3' complementarity to either forward primer
            if has_3prime_complementarity(fwd1.seq_5to3, rev_primer, k=dimer_k):
                continue
            if has_3prime_complementarity(fwd2.seq_5to3, rev_primer, k=dimer_k):
                continue

            score, tm, gc = scored
            candidates.append(
                PrimerHit(
                    kind="REV",
                    exon_name="Exon3",
                    start_0based=i,
                    length=L,
                    seq_5to3=rev_primer,
                    tm_c=tm,
                    gc_pct=gc,
                    score=score,
                    template_seq_5to3=template_seq
                )
            )

    if not candidates:
        raise RuntimeError(
            "No reverse primer found that avoids 3' complementarity to both forward primers. "
            "Try dimer_k=3, widen tm_tol, or widen length range."
        )

    candidates.sort(key=lambda x: x.score)
    return candidates[0]


def design_exon_primers(
    exon1: str,
    exon2: str,
    exon3_for_reverse: str,
    min_len: int = 16,
    max_len: int = 25,
    tm_target: float = 60.0,
    tm_tol: float = 5.0,
    dimer_k: int = 4
) -> Tuple[PrimerHit, PrimerHit, PrimerHit]:
    """
    Returns: (fwd_exon1, fwd_exon2, rev_exon3)
    """
    fwd1 = best_primer_from_exon(exon1, "Exon1", "FWD", min_len, max_len, tm_target, tm_tol)
    fwd2 = best_primer_from_exon(exon2, "Exon2", "FWD", min_len, max_len, tm_target, tm_tol)
    rev3 = best_reverse_avoiding_fwds(exon3_for_reverse, fwd1, fwd2, min_len, max_len, tm_target, tm_tol, dimer_k=dimer_k)
    return fwd1, fwd2, rev3


# ---------------- BLAST (optional) ----------------

def blast_primer_summary_ncbi(
    primer_seq: str,
    database: str = "refseq_rna",
    organism_filter: str = "Homo sapiens[Organism]",
    hitlist_size: int = 10
) -> BlastSummary:
    """
    Remote BLAST via NCBI using Biopython qblast.
    This can fail on hosted platforms (rate limits / blocked outbound).
    """
    try:
        from Bio.Blast import NCBIWWW, NCBIXML
    except Exception as e:
        return BlastSummary(ok=False, note=f"Biopython not available: {e}")

    try:
        handle = NCBIWWW.qblast(
            program="blastn",
            database=database,
            sequence=primer_seq,
            entrez_query=organism_filter,
            hitlist_size=hitlist_size,
            short_query=True,
            format_type="XML"
        )
        record = NCBIXML.read(handle)
        handle.close()

        hits = len(record.alignments)
        if hits == 0:
            return BlastSummary(ok=True, hits_returned=0, note="No hits returned (with current filters).")

        top = record.alignments[0]
        top_hsp = top.hsps[0]
        ident = 100.0 * top_hsp.identities / max(1, top_hsp.align_length)

        return BlastSummary(
            ok=True,
            hits_returned=hits,
            top_hit_id=top.hit_id,
            top_hit_def=top.hit_def,
            top_evalue=float(top_hsp.expect),
            top_identity_pct=float(ident),
            note=f"db={database}, filter={organism_filter}"
        )
    except Exception as e:
        return BlastSummary(ok=False, note=f"BLAST failed: {e}")


def attach_blast(
    primers: List[PrimerHit],
    do_blast: bool = False,
    database: str = "refseq_rna",
    organism_filter: str = "Homo sapiens[Organism]",
    hitlist_size: int = 10
) -> None:
    if not do_blast:
        return
    for p in primers:
        p.blast = blast_primer_summary_ncbi(
            p.seq_5to3,
            database=database,
            organism_filter=organism_filter,
            hitlist_size=hitlist_size
        )


# ---------------- Streamlit dimer report ----------------

def print_dimer_report(f1: PrimerHit, f2: PrimerHit, r3: PrimerHit) -> None:
    """
    Uses Streamlit directly so app.py can call it and render nicely.
    """
    try:
        import streamlit as st
    except Exception:
        return

    rows = []
    for f in [f1, f2]:
        risk = dimer_risk_percent(f.seq_5to3, r3.seq_5to3, max_k=8)
        rows.append({
            "Forward": f"{f.kind} {f.exon_name}",
            "Reverse": f"{r3.kind} {r3.exon_name}",
            "3' complementarity (k)": ">=4 blocked by design",
            "Heuristic dimer risk (%)": round(risk, 1)
        })

    st.table(rows)





### qpcr 

from dataclasses import dataclass
from typing import Optional, Tuple, List

# Reuse your existing helpers if they already exist in primer_engine.py
# If you already have these functions, do NOT duplicate them.
# clean_seq, revcomp, gc_content, tm_wallace, max_run, has_3prime_gc_clamp, self_complementarity_flag, score_candidate
# and PrimerHit dataclass

def parse_junction_marked_seq(seq_with_marker: str, marker: str = "^") -> Tuple[str, str]:
    s = clean_seq(seq_with_marker.replace(marker, ""))
    if marker not in seq_with_marker:
        raise ValueError("For qPCR, mark the exon-exon junction using '^' inside the sequence.")
    left_raw, right_raw = seq_with_marker.split(marker, 1)
    left = clean_seq(left_raw)
    right = clean_seq(right_raw)
    if len(left) < 10 or len(right) < 10:
        raise ValueError("Junction sides are too short. Paste more bases on each side of '^'.")
    return left, right

def build_junction_primer_forward(left: str, right: str, left_take: int, right_take: int) -> str:
    # Primer is written 5'->3' as ordered, crossing junction
    return left[-left_take:] + right[:right_take]

def build_junction_primer_reverse(left: str, right: str, left_take: int, right_take: int) -> str:
    # Reverse primer spans junction on the minus strand
    template = left[-left_take:] + right[:right_take]
    return revcomp(template)

def design_qpcr_junction_primers(
    seq_with_junction_marker: str,
    span_primer: str = "FWD",           # "FWD" or "REV" spans the junction
    min_len: int = 18,
    max_len: int = 25,
    tm_target: float = 60.0,
    tm_tol: float = 5.0,
    amplicon_min: int = 70,
    amplicon_max: int = 200,
    junction_min_overlap_each_side: int = 6
) -> Tuple["PrimerHit", "PrimerHit"]:
    """
    Returns (forward, reverse) primers.
    One primer is forced to span the junction.
    The partner primer is chosen to yield a short amplicon on the spliced template.
    """

    left, right = parse_junction_marked_seq(seq_with_junction_marker, marker="^")

    # Build full spliced template used for partner search
    template = left + right
    junction_index = len(left)  # 0-based index in template where right starts

    best_pair = None  # (score_sum, fwd_hit, rev_hit)

    # Enumerate junction spanning primer lengths by choosing how many bases from each side
    for L in range(min_len, max_len + 1):
        # need at least junction_min_overlap_each_side from each side
        for left_take in range(junction_min_overlap_each_side, L - junction_min_overlap_each_side + 1):
            right_take = L - left_take
            if right_take < junction_min_overlap_each_side:
                continue

            if span_primer.upper() == "FWD":
                jseq = build_junction_primer_forward(left, right, left_take, right_take)
                j_scored = score_candidate(jseq, tm_target, tm_tol)
                if j_scored is None:
                    continue
                j_score, j_tm, j_gc = j_scored
                fwd_hit = PrimerHit(
                    kind="FWD",
                    exon_name="Junction",
                    start_0based=junction_index - left_take,
                    length=L,
                    seq_5to3=jseq,
                    tm_c=j_tm,
                    gc_pct=j_gc,
                    score=j_score
                )

                # Partner reverse primer should bind on the right side to make small amplicon
                # We search reverse primer binding sites on the template within amplicon window
                # Reverse primer binds to template region downstream of forward 3' end
                fwd_3prime_pos = fwd_hit.start_0based + fwd_hit.length  # first base after primer on template
                search_start = fwd_3prime_pos + amplicon_min - 1
                search_end = min(len(template), fwd_3prime_pos + amplicon_max)

                if search_start >= search_end:
                    continue

                best_rev = None
                for Lr in range(min_len, max_len + 1):
                    for i in range(search_start, search_end - Lr + 1):
                        bind_site = template[i:i + Lr]
                        rev_seq = revcomp(bind_site)
                        r_scored = score_candidate(rev_seq, tm_target, tm_tol)
                        if r_scored is None:
                            continue
                        r_score, r_tm, r_gc = r_scored
                        rev_hit = PrimerHit(
                            kind="REV",
                            exon_name="RightSide",
                            start_0based=i,
                            length=Lr,
                            seq_5to3=rev_seq,
                            tm_c=r_tm,
                            gc_pct=r_gc,
                            score=r_score,
                            template_seq_5to3=bind_site
                        )
                        if best_rev is None or rev_hit.score < best_rev.score:
                            best_rev = rev_hit

                if best_rev is None:
                    continue

                total = fwd_hit.score + best_rev.score
                if best_pair is None or total < best_pair[0]:
                    best_pair = (total, fwd_hit, best_rev)

            else:
                # span_primer == "REV": reverse primer spans junction, forward primer upstream on left side
                jseq = build_junction_primer_reverse(left, right, left_take, right_take)
                j_scored = score_candidate(jseq, tm_target, tm_tol)
                if j_scored is None:
                    continue
                j_score, j_tm, j_gc = j_scored
                rev_hit = PrimerHit(
                    kind="REV",
                    exon_name="Junction",
                    start_0based=junction_index - left_take,
                    length=L,
                    seq_5to3=jseq,
                    tm_c=j_tm,
                    gc_pct=j_gc,
                    score=j_score,
                    template_seq_5to3=(left[-left_take:] + right[:right_take])
                )

                # Partner forward primer should be on the left side to make small amplicon
                # Forward primer binds upstream of reverse binding site.
                rev_bind_start = rev_hit.start_0based  # reverse binds here (binding site start on template)
                search_end = max(0, rev_bind_start - amplicon_min + 1)
                search_start = max(0, rev_bind_start - amplicon_max)

                if search_start >= search_end:
                    continue

                best_fwd = None
                for Lf in range(min_len, max_len + 1):
                    for i in range(search_start, search_end - Lf + 1):
                        fseq = template[i:i + Lf]
                        f_scored = score_candidate(fseq, tm_target, tm_tol)
                        if f_scored is None:
                            continue
                        f_score, f_tm, f_gc = f_scored
                        fwd_hit = PrimerHit(
                            kind="FWD",
                            exon_name="LeftSide",
                            start_0based=i,
                            length=Lf,
                            seq_5to3=fseq,
                            tm_c=f_tm,
                            gc_pct=f_gc,
                            score=f_score
                        )
                        if best_fwd is None or fwd_hit.score < best_fwd.score:
                            best_fwd = fwd_hit

                if best_fwd is None:
                    continue

                total = best_fwd.score + rev_hit.score
                if best_pair is None or total < best_pair[0]:
                    best_pair = (total, best_fwd, rev_hit)

    if best_pair is None:
        raise RuntimeError("No qPCR junction primer pair found with your constraints. Try wider Tm tolerance or length range.")

    _, best_fwd, best_rev = best_pair
    return best_fwd, best_rev

def qpcr_amplicon_size(template_seq_with_marker: str, fwd: "PrimerHit", rev: "PrimerHit") -> int:
    # Amplicon size on spliced template = from fwd start to end of rev binding site
    left, right = parse_junction_marked_seq(template_seq_with_marker, marker="^")
    template = left + right
    rev_end = rev.start_0based + rev.length
    return max(0, rev_end - fwd.start_0based)

# ============================
# Regular PCR (single template)
# ============================

def _collect_forward_candidates(template: str, start_window: int, min_len: int, max_len: int,
                                tm_target: float, tm_tol: float) -> List[PrimerHit]:
    template = clean_seq(template)
    window_end = min(len(template), max(0, start_window))
    candidates: List[PrimerHit] = []

    for L in range(min_len, max_len + 1):
        for i in range(0, max(0, window_end - L + 1)):
            seq = template[i:i + L]  # forward primer is same as template window
            scored = score_candidate(seq, tm_target, tm_tol)
            if scored is None:
                continue
            score, tm, gc = scored
            candidates.append(
                PrimerHit(
                    kind="FWD",
                    exon_name="Template",
                    start_0based=i,
                    length=L,
                    seq_5to3=seq,
                    tm_c=tm,
                    gc_pct=gc,
                    score=score
                )
            )

    if not candidates:
        raise RuntimeError("No forward primer found in the start window. Try wider Tm tolerance or length range.")
    candidates.sort(key=lambda x: x.score)
    return candidates[:200]


def _collect_reverse_candidates(template: str, end_window: int, min_len: int, max_len: int,
                                tm_target: float, tm_tol: float) -> List[PrimerHit]:
    template = clean_seq(template)
    n = len(template)
    win = min(n, max(0, end_window))
    window_start = max(0, n - win)
    candidates: List[PrimerHit] = []

    for L in range(min_len, max_len + 1):
        for i in range(window_start, n - L + 1):
            bind_site = template[i:i + L]      # binding site on template (5'->3')
            rev_seq = revcomp(bind_site)       # primer you order (5'->3')
            scored = score_candidate(rev_seq, tm_target, tm_tol)
            if scored is None:
                continue
            score, tm, gc = scored
            candidates.append(
                PrimerHit(
                    kind="REV",
                    exon_name="Template",
                    start_0based=i,
                    length=L,
                    seq_5to3=rev_seq,
                    tm_c=tm,
                    gc_pct=gc,
                    score=score,
                    template_seq_5to3=bind_site
                )
            )

    if not candidates:
        raise RuntimeError("No reverse primer found in the end window. Try wider Tm tolerance or length range.")
    candidates.sort(key=lambda x: x.score)
    return candidates[:200]


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
    dimer_k: int = 4
) -> Tuple[PrimerHit, PrimerHit, int, int, int]:
    """
    Returns (fwd_hit, rev_hit, amplicon_len_bp, fwd_start, rev_start)

    Notes:
    - forward primer searched in first start_window bases
    - reverse primer searched in last end_window bases
    - amplicon length computed on the template as:
        (rev_binding_start + rev_len) - fwd_start
    """
    template = clean_seq(template)
    n = len(template)

    if n < (amplicon_min + min_len * 2):
        raise RuntimeError("Template too short for your amplicon/primer constraints.")

    fwds = _collect_forward_candidates(template, start_window, min_len, max_len, tm_target, tm_tol)
    revs = _collect_reverse_candidates(template, end_window, min_len, max_len, tm_target, tm_tol)

    best = None  # (total_score, fwd_hit, rev_hit, amp_len)

    mid = (amplicon_min + amplicon_max) / 2.0
    span = max(1.0, (amplicon_max - amplicon_min))

    for f in fwds:
        for r in revs:
            # reverse binding site must be downstream of forward start
            amp_len = (r.start_0based + r.length) - f.start_0based
            if amp_len < amplicon_min or amp_len > amplicon_max:
                continue

            # quick heterodimer screen (3' complementarity)
            if has_3prime_complementarity(f.seq_5to3, r.seq_5to3, k=dimer_k):
                continue
            if has_3prime_complementarity(r.seq_5to3, f.seq_5to3, k=dimer_k):
                continue

            # score: primer quality + small penalty for amplicon far from midpoint
            length_penalty = abs(amp_len - mid) / span
            total = f.score + r.score + (length_penalty * 2.0)

            if best is None or total < best[0]:
                best = (total, f, r, amp_len)

    if best is None:
        raise RuntimeError(
            "No primer pair found with your constraints. "
            "Try widening amplicon range, widening Tm tolerance, or reducing dimer_k."
        )

    _, f_best, r_best, amp_best = best
    return f_best, r_best, int(amp_best), f_best.start_0based, r_best.start_0based


def print_dimer_report_pair(fwd_seq: str, rev_seq: str) -> None:
    """
    Streamlit table for heterodimer risk between a single forward and reverse primer sequence.
    """
    try:
        import streamlit as st
    except Exception:
        return

    f = clean_seq(fwd_seq)
    r = clean_seq(rev_seq)

    risk_fr = dimer_risk_percent(f, r, max_k=8)
    risk_rf = dimer_risk_percent(r, f, max_k=8)

    st.table([
        {
            "Pair": "FWD vs REV",
            "Heuristic dimer risk (%)": round(max(risk_fr, risk_rf), 1),
            "Max direction": "FWD->REV" if risk_fr >= risk_rf else "REV->FWD"
        }
    ])
