# utils/blast.py
import urllib.parse

def primer_blast_url_pair(fwd_seq: str, rev_seq: str, organism: str = "Homo sapiens") -> str:
    params = {
        "PRIMER_LEFT_INPUT": fwd_seq,
        "PRIMER_RIGHT_INPUT": rev_seq,
        "ORGANISM": organism,
    }
    return "https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?" + urllib.parse.urlencode(params)

def primer_blast_url_single(primer_seq: str, organism: str = "Homo sapiens") -> str:
    params = {
        "PRIMER_LEFT_INPUT": primer_seq,
        "ORGANISM": organism,
    }
    return "https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?" + urllib.parse.urlencode(params)
