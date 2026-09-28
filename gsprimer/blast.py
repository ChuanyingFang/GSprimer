"""
GSprimer — Primer specificity.

Two layers:

  1. `blast_primer`      : NCBI blastn-short against a nucleotide DB
                           (same invocation as SpacerFinder primer_blast.py)
  2. `insilico_pcr`      : PAIR-level specificity. A primer that has many
                           hits is only a problem when a second primer sits
                           in the opposite orientation within amplifiable
                           distance. This is what actually predicts spurious
                           bands, so it drives the reported risk level.

Fallback: when BLAST is unavailable, `seed_search` does a fast seed-and-extend
scan over a FASTA dictionary (cDNA set or a small genome) in pure Python.
3'-end integrity is what matters for extension, so the seed is taken from the
3' end of the primer.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional

from .thermo import revcomp

MAX_SPURIOUS_AMPLICON = 5000   # bp; anything longer will not amplify in a
                               # standard 2-3 min extension
THREE_PRIME_WINDOW = 5         # nt at the 3' end that must match to extend


# ============================================================
# BLAST
# ============================================================

def blast_available(blastn: str = "blastn", db: str = "") -> bool:
    exe = shutil.which(blastn) or (blastn if os.path.exists(blastn) else "")
    if not exe:
        return False
    if db:
        # a BLAST DB is a set of sibling files sharing the given prefix
        d = os.path.dirname(db) or "."
        base = os.path.basename(db)
        if not os.path.isdir(d):
            return False
        if not any(f.startswith(base) for f in os.listdir(d)):
            return False
    return True


def blast_primer(seq: str, blastn: str = "blastn", db: str = "",
                 label: str = "primer", max_hits: int = 100) -> Dict:
    """blastn-short for a single primer. Returns {'hits': [...]} or {'error':}."""
    if not blast_available(blastn, db):
        return {"error": "blastn or BLAST database unavailable", "hits": []}
    qlen = len(seq)
    fd, fa = tempfile.mkstemp(suffix=".fa", prefix="gsp_")
    with os.fdopen(fd, "w") as f:
        f.write(f">{label}\n{seq}\n")
    fd2, out = tempfile.mkstemp(suffix=".tsv", prefix="gsp_")
    os.close(fd2)
    try:
        subprocess.run(
            [blastn, "-task", "blastn-short", "-db", db,
             "-query", fa, "-out", out,
             "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                        "qstart qend sstart send evalue bitscore",
             "-evalue", "1000", "-word_size", "7", "-dust", "no",
             "-max_target_seqs", str(max_hits), "-num_alignments", str(max_hits)],
            capture_output=True, timeout=120)
        hits = []
        with open(out) as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) < 12:
                    continue
                ss, se = int(c[8]), int(c[9])
                qs, qe = int(c[6]), int(c[7])
                strand = "+" if ss <= se else "-"
                hits.append({
                    "subject": c[1], "pident": float(c[2]),
                    "aln_len": int(c[3]), "mismatch": int(c[4]),
                    "gaps": int(c[5]), "q_start": qs, "q_end": qe,
                    "s_start": ss, "s_end": se, "strand": strand,
                    "evalue": float(c[10]),
                    "covers_3p": qe == qlen,
                    # genomic position of the primer 3' base
                    "three_prime_pos": se,
                    # total mismatches counting unaligned query tails
                    "total_mm": int(c[4]) + (qlen - (qe - qs + 1)),
                })
        return {"hits": hits, "qlen": qlen}
    except subprocess.TimeoutExpired:
        return {"error": "BLAST timed out", "hits": []}
    except Exception as e:
        return {"error": str(e), "hits": []}
    finally:
        for p in (fa, out):
            try:
                os.unlink(p)
            except OSError:
                pass


# ============================================================
# Pure-Python fallback
# ============================================================

def seed_search(seq: str, db: Dict[str, str], seed_len: int = 12,
                max_mm: int = 4, max_hits: int = 500) -> List[Dict]:
    """Seed-and-extend scan. Seed = 3'-most `seed_len` nt (extension-critical).

    Returns hit dicts compatible with `blast_primer` output.
    """
    s = seq.upper()
    n = len(s)
    seed_len = min(seed_len, n)
    seed_f = s[-seed_len:]                 # primer 3' end, plus orientation
    rc = revcomp(s)
    seed_r = rc[:seed_len]                 # same 3' end, minus orientation
    hits: List[Dict] = []

    for name, subject in db.items():
        # ---- plus strand ----
        pos = subject.find(seed_f)
        while pos != -1 and len(hits) < max_hits:
            start = pos + seed_len - n     # where the primer 5' end would sit
            if start >= 0 and start + n <= len(subject):
                window = subject[start:start + n]
                mm = sum(1 for a, b in zip(window, s) if a != b)
                if mm <= max_mm:
                    hits.append(_mk_hit(name, start + 1, start + n, "+",
                                        mm, n))
            pos = subject.find(seed_f, pos + 1)
        # ---- minus strand ----
        pos = subject.find(seed_r)
        while pos != -1 and len(hits) < max_hits:
            if pos + n <= len(subject):
                window = subject[pos:pos + n]
                mm = sum(1 for a, b in zip(window, rc) if a != b)
                if mm <= max_mm:
                    # primer 3' base maps to the LEFT edge on the minus strand
                    hits.append(_mk_hit(name, pos + n, pos + 1, "-", mm, n))
            pos = subject.find(seed_r, pos + 1)
    return hits


def _mk_hit(subject: str, s_start: int, s_end: int, strand: str,
            mm: int, qlen: int) -> Dict:
    return {"subject": subject, "pident": 100.0 * (qlen - mm) / qlen,
            "aln_len": qlen, "mismatch": mm, "gaps": 0,
            "q_start": 1, "q_end": qlen, "s_start": s_start, "s_end": s_end,
            "strand": strand, "evalue": 0.0, "covers_3p": True,
            "three_prime_pos": s_end, "total_mm": mm}


# ============================================================
# Hit classification
# ============================================================

def _gene_key(subject: str) -> str:
    """Collapse transcript isoforms to their parent gene.

    LOC_Os01g01010.1 / LOC_Os01g01010.2 -> LOC_Os01g01010.
    Chromosome / scaffold names (chr1, chrD2) are returned unchanged.
    """
    if "." in subject:
        head, tail = subject.rsplit(".", 1)
        if tail.isdigit() and head:
            return head
    return subject


def classify(hits: List[Dict], qlen: int) -> Dict:
    """Bucket hits by how likely they are to prime."""
    perfect, strong, weak = [], [], []
    for h in hits:
        full = h["aln_len"] >= qlen - 1 and h["gaps"] == 0
        if h["total_mm"] == 0 and full:
            perfect.append(h)
        elif h["covers_3p"] and h["total_mm"] <= 3:
            strong.append(h)
        else:
            weak.append(h)
    perfect_genes = {_gene_key(h["subject"]) for h in perfect}
    return {"perfect": perfect, "strong": strong, "weak": weak,
            "n_perfect": len(perfect), "n_perfect_genes": len(perfect_genes),
            "n_strong": len(strong), "n_weak": len(weak),
            "n_total": len(hits)}


def primer_risk(cls: Dict) -> str:
    """HIGH / MEDIUM / LOW for a single primer.

    HIGH only when the primer perfectly matches MORE THAN ONE distinct gene
    (isoforms of the same gene do not count — amplifying every isoform of the
    target locus is expected and desirable for CDS cloning).
    """
    if cls.get("n_perfect_genes", 0) > 1:
        return "HIGH"
    if cls["n_perfect"] == 1 and cls["n_strong"] >= 3:
        return "MEDIUM"
    if cls["n_perfect"] == 1 and cls["n_strong"] >= 1:
        return "LOW-MED"
    if cls["n_perfect"] == 0:
        return "UNKNOWN"
    return "LOW"


# ============================================================
# In-silico PCR
# ============================================================

def insilico_pcr(f_hits: List[Dict], r_hits: List[Dict],
                 expect_len: int = 0,
                 max_len: int = MAX_SPURIOUS_AMPLICON) -> Dict:
    """Predict every amplicon a primer pair could generate.

    An amplicon requires a plus-oriented 3' end upstream of a minus-oriented
    3' end on the same subject, within `max_len`.
    Also checks F+F and R+R self-amplification.
    """
    products = []
    tagged = ([dict(h, primer="F") for h in f_hits]
              + [dict(h, primer="R") for h in r_hits])
    by_subject: Dict[str, List[Dict]] = {}
    for h in tagged:
        if not h["covers_3p"] or h["total_mm"] > 3:
            continue      # cannot extend efficiently
        by_subject.setdefault(h["subject"], []).append(h)

    for subj, hs in by_subject.items():
        plus = [h for h in hs if h["strand"] == "+"]
        minus = [h for h in hs if h["strand"] == "-"]
        for p in plus:
            for m in minus:
                size = m["three_prime_pos"] - p["three_prime_pos"] + 1
                if 0 < size <= max_len:
                    products.append({
                        "subject": subj, "size": size,
                        "left": p["primer"], "right": m["primer"],
                        "left_mm": p["total_mm"], "right_mm": m["total_mm"],
                        "left_pos": p["three_prime_pos"],
                        "right_pos": m["three_prime_pos"],
                        "is_target": (p["primer"] == "F" and m["primer"] == "R"
                                      and p["total_mm"] == 0
                                      and m["total_mm"] == 0
                                      and (expect_len == 0
                                           or abs(size - expect_len) <= 5)),
                    })
    products.sort(key=lambda x: (not x["is_target"], x["left_mm"] + x["right_mm"],
                                 x["size"]))
    on_target = [p for p in products if p["is_target"]]
    off_target = [p for p in products if not p["is_target"]]
    perfect_off = [p for p in off_target
                   if p["left_mm"] == 0 and p["right_mm"] == 0]

    if perfect_off:
        level = "HIGH"
    elif len(on_target) != 1 and expect_len:
        level = "HIGH" if not on_target else "MEDIUM"
    elif off_target:
        level = "MEDIUM" if len(off_target) > 2 else "LOW-MED"
    else:
        level = "LOW"

    return {"products": products, "on_target": on_target,
            "off_target": off_target, "n_off": len(off_target),
            "n_perfect_off": len(perfect_off), "level": level}


# ============================================================
# Orchestration
# ============================================================

def check_pair(f_seq: str, r_seq: str, *, blastn: str = "blastn",
               db: str = "", fallback_db: Optional[Dict[str, str]] = None,
               expect_len: int = 0) -> Dict:
    """Full specificity workup for one primer pair."""
    method = "blast"
    if blast_available(blastn, db):
        fr = blast_primer(f_seq, blastn, db, "F")
        rr = blast_primer(r_seq, blastn, db, "R")
        if fr.get("error") or rr.get("error"):
            method = "error"
            f_hits, r_hits = fr.get("hits", []), rr.get("hits", [])
            err = fr.get("error") or rr.get("error")
        else:
            f_hits, r_hits = fr["hits"], rr["hits"]
            err = ""
    elif fallback_db:
        method = "seed-scan"
        f_hits = seed_search(f_seq, fallback_db)
        r_hits = seed_search(r_seq, fallback_db)
        err = ""
    else:
        return {"method": "none", "available": False,
                "error": "no BLAST database configured and no local FASTA "
                         "provided; specificity cannot be assessed",
                "level": "UNKNOWN"}

    f_cls = classify(f_hits, len(f_seq))
    r_cls = classify(r_hits, len(r_seq))
    pcr = insilico_pcr(f_hits, r_hits, expect_len=expect_len)

    levels = {"LOW": 0, "LOW-MED": 1, "MEDIUM": 2, "HIGH": 3, "UNKNOWN": 2}
    worst = max(primer_risk(f_cls), primer_risk(r_cls), pcr["level"],
                key=lambda x: levels.get(x, 0))

    return {"method": method, "available": True, "error": err,
            "f": f_cls, "r": r_cls,
            "f_risk": primer_risk(f_cls), "r_risk": primer_risk(r_cls),
            "pcr": pcr, "level": worst}
