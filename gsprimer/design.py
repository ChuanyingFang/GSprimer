"""
GSprimer — CDS amplification primer design for GS (Golden Gate / BsaI) vectors.

Design rules
------------
Thermodynamic rules are inherited verbatim from SpacerFinder v2.4.5
(see thermo.py).  The ONE deliberate deviation is the product-size window:
SpacerFinder enforces 300-2000 bp, GSprimer does not constrain amplicon size
because it is dictated by the ORF.

GS-specific anchoring
---------------------
  Forward primer : 5' end anchored at the ATG.
                   Small shifts are allowed but the offset MUST be a multiple
                   of 3 so that the reading frame relative to the vector's
                   CGAG overhang is preserved:
                       ATG_index - F_start_index  ==  0 (mod 3)
                   offset < 0  -> primer starts upstream in the 5'UTR
                   offset > 0  -> primer starts inside the CDS (N-terminal
                                  residues, including Met, are LOST)

  Reverse primer : 3' end anchored on the last base of the stop codon.
                   Only the 5' end floats (length 19-24 nt).
                   Downstream extension into the 3'UTR is possible via
                   `r_ext` but is OFF by default because it would insert
                   extra bases between the stop codon and the vector's GGAT
                   overhang.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from . import thermo as T
from .template import Transcript

# relaxed band used to classify "borderline" primers (Tier B)
GC_SOFT = (0.35, 0.65)
TM_SOFT = (52.0, 63.0)


@dataclass
class PrimerPair:
    # forward
    f_seq: str = ""
    f_start: int = 0            # 0-based index into mrna
    f_len: int = 0
    f_qc: Dict = field(default_factory=dict)
    offset: int = 0             # f_start - cds_start  (multiple of 3)
    # reverse
    r_seq: str = ""             # actual primer (reverse complement strand)
    r_template: str = ""        # sense-strand footprint
    r_end: int = 0              # 0-based exclusive end index into mrna
    r_len: int = 0
    r_qc: Dict = field(default_factory=dict)
    r_ext: int = 0              # nt of 3'UTR appended past the stop codon
    # pair-level
    tm_diff: float = 0.0
    hetero_dimer_dg: float = 0.0
    amplicon: str = ""
    amplicon_len: int = 0
    score: float = 0.0
    tier: str = "C"             # A = fully compliant, B = borderline, C = forced
    is_anchor: bool = False     # strict ATG..stop pair
    flags: List[str] = field(default_factory=list)
    # filled in later by blast / risk modules
    blast: Dict = field(default_factory=dict)
    risk: Dict = field(default_factory=dict)

    @property
    def label(self) -> str:
        tag = "ATG..STOP" if self.is_anchor else f"offset{self.offset:+d}"
        return f"{tag}/F{self.f_len}-R{self.r_len}"


# ============================================================
# Candidate enumeration
# ============================================================

def _forward_candidates(tx: Transcript, offsets: List[int],
                        len_range: Tuple[int, int]) -> List[Dict]:
    out = []
    n = len(tx.mrna)
    for off in offsets:
        if off % 3 != 0:
            continue
        start = tx.cds_start + off
        if start < 0 or start >= tx.cds_end:
            continue
        for plen in range(len_range[0], len_range[1] + 1):
            if start + plen > n:
                continue
            seq = tx.mrna[start:start + plen]
            if "N" in seq:
                continue
            qc = T.evaluate_primer(seq)
            out.append({"seq": seq, "start": start, "len": plen,
                        "offset": off, "qc": qc})
    return out


def _reverse_candidates(tx: Transcript, exts: List[int],
                        len_range: Tuple[int, int]) -> List[Dict]:
    out = []
    n = len(tx.mrna)
    for ext in exts:
        end = tx.cds_end + ext
        if end > n or end <= tx.cds_start:
            continue
        for plen in range(len_range[0], len_range[1] + 1):
            start = end - plen
            if start < 0:
                continue
            foot = tx.mrna[start:end]
            if "N" in foot:
                continue
            seq = T.revcomp(foot)
            qc = T.evaluate_primer(seq)
            out.append({"seq": seq, "template": foot, "end": end,
                        "len": plen, "ext": ext, "qc": qc})
    return out


def _tier_of(f_qc: Dict, r_qc: Dict, tm_diff: float, hetero_ok: bool,
             hetero_dg: float, structural_ok: bool = True) -> str:
    """A = every inherited SpacerFinder filter passes; B = soft band; C = worse.

    `structural_ok` is False when the pair breaks the GS anchoring contract
    (N-terminus truncated, or extra bases past the stop codon). Such pairs can
    never be Tier A/B no matter how good the thermodynamics are.
    """
    if not structural_ok:
        return "C"
    if f_qc["pass"] and r_qc["pass"] and tm_diff < T.TM_DIFF_MAX and hetero_ok:
        return "A"
    for qc in (f_qc, r_qc):
        if not (GC_SOFT[0] <= qc["gc"] <= GC_SOFT[1]):
            return "C"
        if not (TM_SOFT[0] <= qc["tm"] <= TM_SOFT[1]):
            return "C"
        if qc["homopolymer"] >= 6:
            return "C"
        if qc["hairpin_dg"] < -6.0 or qc["self_dimer_dg"] < -9.0:
            return "C"
    if tm_diff >= 4.0 or hetero_dg < -9.0:
        return "C"
    return "B"


def _score(f: Dict, r: Dict, tm_diff: float) -> float:
    """SpacerFinder scoring + GS-specific anchoring penalties."""
    f_qc, r_qc = f["qc"], r["qc"]
    s = (tm_diff * 10.0
         + (2.0 - f_qc["end_bonus"] - r_qc["end_bonus"])
         + abs(f_qc["gc"] - r_qc["gc"]) * 5.0
         + abs(f["len"] - 21) * 0.5
         + abs(r["len"] - 21) * 0.5)
    # GS anchoring penalties: staying on the ATG is worth a lot
    s += abs(f["offset"]) / 3.0 * 2.0
    if f["offset"] > 0:
        s += 40.0                     # loses the native N-terminus + Met
    if r["ext"] > 0:
        s += 40.0                     # inserts extra bases before the overhang
    # over-length primers are a deliberate relaxation, mild penalty
    for cand in (f, r):
        if cand["len"] > T.LEN_MAX:
            s += (cand["len"] - T.LEN_MAX) * 0.4
    # soft-band penalties so Tier B never outranks Tier A
    for qc in (f_qc, r_qc):
        if not qc["pass"]:
            s += 3.0 * len(qc["flags"])
    return s


def _build_pair(tx: Transcript, f: Dict, r: Dict) -> PrimerPair:
    tm_diff = abs(f["qc"]["tm"] - r["qc"]["tm"])
    hetero_dg = T.worst_dimer_dg(f["seq"], r["seq"])
    hetero_ok = not T.has_dimer(f["seq"], r["seq"])
    amp = tx.mrna[f["start"]:r["end"]]
    structural_ok = (f["offset"] <= 0 and r["ext"] == 0)
    pair = PrimerPair(
        f_seq=f["seq"], f_start=f["start"], f_len=f["len"], f_qc=f["qc"],
        offset=f["offset"],
        r_seq=r["seq"], r_template=r["template"], r_end=r["end"],
        r_len=r["len"], r_qc=r["qc"], r_ext=r["ext"],
        tm_diff=tm_diff, hetero_dimer_dg=hetero_dg,
        amplicon=amp, amplicon_len=len(amp),
        score=_score(f, r, tm_diff),
        tier=_tier_of(f["qc"], r["qc"], tm_diff, hetero_ok, hetero_dg,
                      structural_ok),
        is_anchor=(f["offset"] == 0 and r["ext"] == 0),
    )
    flags = []
    flags += [f"F:{x}" for x in f["qc"]["flags"]]
    flags += [f"R:{x}" for x in r["qc"]["flags"]]
    if tm_diff >= T.TM_DIFF_MAX:
        flags.append(f"dTm={tm_diff:.1f}")
    if not hetero_ok:
        flags.append(f"HeteroDimer={hetero_dg:.1f}")
    if f["offset"] < 0:
        flags.append(f"5'UTR+{-f['offset']}nt")
        utr = tx.mrna[f["start"]:tx.cds_start]
        stops = [i for i in range(0, len(utr) - 2, 3)
                 if utr[i:i + 3] in ("TAA", "TAG", "TGA")]
        if stops:
            flags.append("in-frame stop in UTR")
    if f["offset"] > 0:
        flags.append(f"N-term missing {f['offset'] // 3} aa")
    if r["ext"] > 0:
        flags.append(f"+{r['ext']} nt past stop")
    pair.flags = flags
    return pair


# ============================================================
# Main entry
# ============================================================

def design(tx: Transcript, *, max_shift: int = 12, r_ext_max: int = 0,
           len_range: Tuple[int, int] = (T.LEN_MIN, T.LEN_MAX),
           top_n: int = 10, allow_downstream: bool = True) -> List[PrimerPair]:
    """Enumerate and rank GS cloning primer pairs.

    Parameters
    ----------
    max_shift        : max |offset| for the forward primer (nt, rounded to 3)
    r_ext_max        : max 3'UTR extension for the reverse primer (0 = strict)
    allow_downstream : permit offset > 0 (truncates the N-terminus)
    top_n            : number of pairs returned (the strict ATG..stop anchor
                       pair is ALWAYS included, on top of / within top_n)

    Returns pairs sorted by (tier, score). Element 0 is guaranteed to be a
    strict ATG..stop pair when one can be built at all.
    """
    max_shift = (max_shift // 3) * 3
    upstream = [-o for o in range(3, max_shift + 1, 3)
                if o <= tx.utr5_len]
    downstream = [o for o in range(3, max_shift + 1, 3)] if allow_downstream else []
    offsets = [0] + upstream + downstream

    exts = [e for e in range(0, r_ext_max + 1, 3) if e <= tx.utr3_len]
    if not exts:
        exts = [0]

    f_cands = _forward_candidates(tx, offsets, len_range)
    r_cands = _reverse_candidates(tx, exts, len_range)
    if not f_cands or not r_cands:
        return []

    pairs = [_build_pair(tx, f, r) for f in f_cands for r in r_cands]

    tier_rank = {"A": 0, "B": 1, "C": 2}
    pairs.sort(key=lambda p: (tier_rank[p.tier], p.score))

    # --- requirement: the candidate set MUST contain a strict ATG..stop pair
    anchors = [p for p in pairs if p.is_anchor]
    best_anchor = anchors[0] if anchors else None

    selected: List[PrimerPair] = []
    if best_anchor is not None:
        best_anchor.flags = ["\u2605 standard anchor pair"] + best_anchor.flags
        selected.append(best_anchor)
    for p in pairs:
        if len(selected) >= top_n:
            break
        if p is best_anchor:
            continue
        selected.append(p)
    return selected


def summarize_isoforms(txs: List[Transcript]) -> List[Dict]:
    """Table rows describing every isoform of the queried gene."""
    rows = []
    for t in sorted(txs, key=lambda x: (-x.cds_len, x.tx_id)):
        rows.append({
            "transcript": t.tx_id,
            "mRNA_len": len(t.mrna),
            "CDS_len": t.cds_len,
            "protein_aa": t.protein_len,
            "utr5": t.utr5_len,
            "utr3": t.utr3_len,
            "exons": len(t.exons),
            "issues": "; ".join(t.validate()) or "OK",
        })
    return rows
