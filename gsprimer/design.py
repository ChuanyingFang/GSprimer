"""
GSprimer — CDS amplification primer design for GS (Golden Gate / BsaI) vectors.

Design rules
------------
Thermodynamic rules are inherited verbatim from SpacerFinder v2.4.5
(see thermo.py).  The ONE deliberate deviation is the product-size window:
SpacerFinder enforces 300-2000 bp, GSprimer does not constrain amplicon size
for CDS because it is dictated by the ORF; promoter mode (v2.0) instead targets
a 1200–1700 bp window upstream of the TSS.

GS-specific anchoring (CDS mode)
--------------------------------
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

Promoter-mode anchoring (v2.0)
------------------------------
  Reverse primer : 5' end (the annealing/distal end) is anchored one nucleotide
                   BEFORE the 5'-UTR start (the TSS) when ideal (offset == 0);
                   if not ideal it may (a) extend upstream into the promoter
                   region (offset < 0, up to 100 nt) or (b) fall back into the
                   transcript 5'-UTR (offset > 0, up to 200 nt; two-file
                   strategy). The forward primer is then placed in the promoter
                   region so the product lands in the 1200–1700 bp window.
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
    offset: int = 0             # f_start - cds_start (CDS mode) /
                               # R 5' end shift vs TSS (promoter mode)
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
    is_anchor: bool = False     # strict ATG..stop (CDS) / strict TSS anchor (promoter)
    flags: List[str] = field(default_factory=list)
    # filled in later by blast / risk modules
    blast: Dict = field(default_factory=dict)
    risk: Dict = field(default_factory=dict)
    # promoter-mode fields
    mode: str = "cds"           # "cds" | "promoter"
    utr5_idx: int = 0           # 0-based TSS position in the amplicon template
    r_in_utr: bool = False      # R 5' end fell back into the 5'-UTR (2 kb exhausted)

    @property
    def label(self) -> str:
        tag = "ATG..STOP" if self.is_anchor else f"offset{self.offset:+d}"
        return f"{tag}/F{self.f_len}-R{self.r_len}"


# ============================================================
# Candidate enumeration (CDS mode)
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
# Main entry (CDS mode)
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


# ============================================================
# Promoter-mode design (v2.0)
# ============================================================

def _promoter_tier(f_qc: Dict, r_qc: Dict, tm_diff: float, hetero_ok: bool,
                  hetero_dg: float, structural_ok: bool) -> str:
    """A = R precisely anchored at TSS and all thermodynamics pass; B = within
    tolerance but soft band; C = otherwise."""
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


def _promoter_score(f_qc: Dict, r_len: int, r_shift: int, tm_diff: float,
                    amplicon_len: int, mid: float = 1450.0) -> float:
    """Preference: R precisely anchored at TSS, F/R Tm close, product near the
    window centre, GC near 50%."""
    s = (abs(r_shift) / 10.0 * 3.0
         + tm_diff * 8.0
         + abs(amplicon_len - mid) * 0.02
         + abs(f_qc["gc"] - 0.5) * 10.0
         + abs(f_qc["len"] - 21) * 0.5
         + abs(r_len - 21) * 0.5)
    if not f_qc["pass"]:
        s += 3.0 * len(f_qc["flags"])
    return s


def _build_promoter_pair(s: str, fc: Dict, r_seq: str, r_qc: Dict,
                        r_start: int, r_5end: int, r_shift: int,
                        utr5: int) -> PrimerPair:
    f_seq = fc["seq"]; f_start = fc["start"]; f_qc = fc["qc"]; f_len = fc["len"]
    r_len = r_5end - r_start
    tm_diff = abs(f_qc["tm"] - r_qc["tm"])
    hetero_dg = T.worst_dimer_dg(f_seq, r_seq)
    hetero_ok = not T.has_dimer(f_seq, r_seq)
    amp = s[f_start:r_5end]
    amplicon_len = len(amp)
    structural_ok = (r_shift == 0)
    tier = _promoter_tier(f_qc, r_qc, tm_diff, hetero_ok, hetero_dg, structural_ok)
    score = _promoter_score(f_qc, r_len, r_shift, tm_diff, amplicon_len)
    p = PrimerPair(
        f_seq=f_seq, f_start=f_start, f_len=f_len, f_qc=f_qc,
        offset=r_shift,                     # reused: R 5' end shift vs TSS
        r_seq=r_seq, r_template=s[r_start:r_5end], r_end=r_5end, r_len=r_len,
        r_qc=r_qc, r_ext=max(0, r_shift),  # r_shift>0 means falling back into 5'-UTR
        tm_diff=tm_diff, hetero_dimer_dg=hetero_dg,
        amplicon=amp, amplicon_len=amplicon_len,
        score=score, tier=tier,
        is_anchor=(r_shift == 0),
        mode="promoter", utr5_idx=utr5,
        r_in_utr=(r_shift > 0),
    )
    flags = [f"F:{x}" for x in f_qc["flags"]] + [f"R:{x}" for x in r_qc["flags"]]
    if tm_diff >= T.TM_DIFF_MAX:
        flags.append(f"dTm={tm_diff:.1f}")
    if not hetero_ok:
        flags.append(f"HeteroDimer={hetero_dg:.1f}")
    if r_shift != 0:
        flags.append(f"R@TSS{'%+d' % r_shift}nt")
    p.flags = flags
    return p


def design_promoter(tx: Transcript, *,
                    product_min: int = 1200, product_max: int = 1700,
                    r_up_max: int = 100, r_utr_max: int = 200,
                    len_range: Tuple[int, int] = (T.LEN_MIN, T.LEN_MAX),
                    top_n: int = 10) -> List[PrimerPair]:
    """Promoter-amplification primer design (GSprimer v2.0).

    The input is a promoter-mode `Transcript`:
      - ``tx.mrna`` = the upstream promoter region (2 kb) concatenated before the
        5'-UTR
      - ``tx.utr5_start`` = the TSS (5'-UTR start) 0-based index into
        ``tx.mrna``, i.e. the length of the upstream region; the sequence after
        it is the transcript UTR (used only for the fallback).

    Design rules (per user spec):
      1) The reverse primer's 5' end (the annealing/distal end) is anchored
         precisely one nucleotide before the 5'-UTR start (r_shift == 0); if not
         ideal:
           - upstream extension into the promoter region is allowed (r_shift<0,
             up to 100 nt);
           - if no ideal R is found within the 2 kb, fall back to the transcript
             UTR, extending into the 5'-UTR (r_shift>0, up to 200 nt; two-file
             strategy).
      2) PCR product defaults to 1200~1700 bp; the forward primer is placed in
         the promoter region accordingly (F always stays upstream of the TSS).
      3) Adapter structure is identical to CDS mode (controlled by run_design's
         adapter arguments).
    Promoter is non-coding, so there is no frame / synonymous-mutation constraint.

    Two-file strategy: design R in the 2 kb upstream region first (r_shift <= 0);
    only when no Tier A/B R exists within the 2 kb are the r_shift>0 (transcript
    UTR fallback) candidates admitted to the candidate set.
    """
    s = tx.mrna.upper()
    n = len(s)
    utr5 = tx.utr5_start               # TSS index (upstream-region length)
    if not (0 <= utr5 <= n):
        return []

    mid = (product_min + product_max) / 2.0

    # ---- pre-scan forward-primer candidates (every position potentially covered
    #      by the R 5' end - product window) ----
    # R 5' end max reaches utr5 + r_utr_max, so F start min = utr5 + r_utr_max - product_max
    # R 5' end min is utr5 - r_up_max, so F start max = utr5 - r_up_max - product_min
    f_lo = max(0, (utr5 - r_utr_max) - product_max)
    f_hi = max(0, (utr5 - r_up_max) - product_min)
    f_cands: List[Dict] = []
    for start in range(f_lo, f_hi + 1):
        best = None
        for plen in range(len_range[0], len_range[1] + 1):
            foot = s[start:start + plen]
            if len(foot) < plen or "N" in foot:
                continue
            qc = T.evaluate_primer(foot)
            if not qc["pass"]:
                continue
            # pick the length that brings Tm closest to 58C and GC to 50%
            sc = abs(qc["tm"] - 58.0) + abs(qc["gc"] - 0.5) * 20
            if best is None or sc < best[0]:
                best = (sc, foot, plen, qc)
        if best:
            f_cands.append({"start": start, "seq": best[1],
                            "len": best[2], "qc": best[3]})

    # ---- enumerate reverse-primer candidates (two passes: 2 kb r_shift<=0,
    #      then UTR fallback r_shift>0) ----
    # Promoter mode has NO frame constraint, so scan R 5' end nucleotide by
    # nucleotide (including the offset=0 exactly-anchored-to-TSS pair).
    pairs: List[PrimerPair] = []
    for r_shift in range(-r_up_max, r_utr_max + 1):
        r_5end = utr5 + r_shift
        if r_5end <= 0 or r_5end > n:
            continue
        for r_len in range(len_range[0], len_range[1] + 1):
            r_start = r_5end - r_len
            if r_start < 0:
                continue
            foot = s[r_start:r_5end]
            if "N" in foot:
                continue
            r_seq = T.revcomp(foot)
            r_qc = T.evaluate_primer(r_seq)
            if not r_qc["pass"]:
                continue
            f_lo_w = r_5end - product_max
            f_hi_w = r_5end - product_min
            cands = [fc for fc in f_cands if f_lo_w <= fc["start"] <= f_hi_w]
            if not cands:
                continue
            cands.sort(key=lambda fc: (abs(fc["start"] - (r_5end - mid))
                                      + abs(fc["qc"]["tm"] - r_qc["tm"]) * 2))
            for fc in cands[:6]:
                pairs.append(_build_promoter_pair(
                    s, fc, r_seq, r_qc, r_start, r_5end, r_shift, utr5))

    if not pairs:
        return []

    # ---- two-file strategy resolution ----
    pass1 = [p for p in pairs if p.offset <= 0]          # within 2 kb (incl. exact anchor)
    pass2 = [p for p in pairs if p.offset > 0]           # fallback to transcript UTR
    good1 = [p for p in pass1 if p.tier in ("A", "B")]
    if good1:
        base = pass1                                    # 2 kb already has an ideal R; no UTR fallback
    else:
        base = pass1 + pass2                            # 2 kb not ideal; admit UTR fallback candidates

    # ---- rank and keep the strict-anchor pair ----
    tier_rank = {"A": 0, "B": 1, "C": 2}
    base.sort(key=lambda p: (tier_rank[p.tier], p.score))
    anchors = [p for p in base if p.is_anchor]
    selected: List[PrimerPair] = []
    if anchors:
        anchors[0].flags = ["\u2605 standard anchor pair (R precisely anchored at TSS)"] \
                          + anchors[0].flags
        selected.append(anchors[0])
    for p in base:
        if len(selected) >= top_n:
            break
        if anchors and p is anchors[0]:
            continue
        selected.append(p)
    return selected[:top_n]


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
