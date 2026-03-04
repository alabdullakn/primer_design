from dataclasses import dataclass
from typing import List

from Bio.Blast import NCBIWWW, NCBIXML


@dataclass
class BlastHit:
    title: str
    accession: str
    identity_pct: float
    query_coverage_pct: float
    evalue: float
    score: float


def run_qblast(
    seq: str,
    organism: str = "Homo sapiens",
    max_hits: int = 10,
) -> List[BlastHit]:
    """
    Run NCBI QBLAST for a primer sequence and return the top hits.

    Uses word_size=7 (required for short sequences like primers) and a high
    e-value cutoff so imperfect off-target matches are still returned.
    """
    handle = NCBIWWW.qblast(
        "blastn",
        "nt",
        seq,
        word_size=7,
        expect=1000,
        perc_ident=70,
        hitlist_size=max_hits,
        entrez_query=f"{organism}[Organism]",
    )

    records = list(NCBIXML.parse(handle))
    if not records:
        return []

    record = records[0]
    query_len = record.query_length
    hits = []

    for alignment in record.alignments:
        hsp = alignment.hsps[0]
        identity_pct = 100.0 * hsp.identities / hsp.align_length
        query_coverage_pct = 100.0 * hsp.align_length / query_len
        hits.append(BlastHit(
            title=alignment.title[:100],
            accession=alignment.accession,
            identity_pct=round(identity_pct, 1),
            query_coverage_pct=round(query_coverage_pct, 1),
            evalue=hsp.expect,
            score=float(hsp.score),
        ))

    return hits


def specificity_label(hits: List[BlastHit]) -> str:
    """
    Return a simple specificity interpretation based on hit pattern.
    """
    if not hits:
        return "No hits — check organism setting"
    perfect = [h for h in hits if h.identity_pct == 100.0 and h.query_coverage_pct >= 90.0]
    if len(perfect) == 1:
        return "Specific (single perfect match)"
    if len(perfect) == 0:
        return "No perfect match found"
    return f"Potentially non-specific ({len(perfect)} perfect matches)"
