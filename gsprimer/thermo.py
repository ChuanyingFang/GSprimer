"""
GSprimer — Primer thermodynamics.

Cloned verbatim (parameters + thresholds) from SpacerFinder v2.4.5
`spacerfinder/primers.py` so that GSprimer primers are directly comparable
with SpacerFinder genotyping primers.

Inherited rules
---------------
  - Salt + Mg2+ adjusted Tm (Wallace rule, Owczarzy correction)
  - NN dG (SantaLucia 1998 unified parameters)
  - Self-dimer / hetero-dimer / hairpin detection
  - Length 19-24 nt, GC 40-60%
  - Tm 55-60 C, |Tm_F - Tm_R| < 2 C
  - 3' scoring: last base G/C +0.5; exactly 3 G/C in last 5 bp and not a
    3-run  ->  +0.8
  - dimer dG cutoff -6.0 kcal/mol (-5.0 when the duplex reaches a 3' end)
  - hairpin dG cutoff -4.0 kcal/mol

NOT inherited: product size window (SpacerFinder uses 300-2000 bp; a GS CDS
amplicon is whatever the ORF happens to be).
"""

import math
from typing import Dict, List

COMP_TABLE = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
              "a": "t", "t": "a", "g": "c", "c": "g", "n": "n"}

# ---- inherited thresholds -------------------------------------------------
LEN_MIN, LEN_MAX = 19, 24
GC_MIN, GC_MAX = 0.40, 0.60
TM_MIN, TM_MAX = 55.0, 60.0
TM_DIFF_MAX = 2.0
HAIRPIN_DG = -4.0
DIMER_DG = -6.0
DIMER_DG_3P = -5.0

# ============================================================
# NN parameters (SantaLucia 1998)
# ============================================================
_NN_P = {
    "AA/TT": (-7.9, -22.2), "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -21.3), "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4), "GG/CC": (-8.0, -19.9),
}
_NN_INIT = (0.2, -5.7)
_NN_TAT = (2.2, 6.9)


def revcomp(seq: str) -> str:
    return "".join(COMP_TABLE.get(b, "N") for b in reversed(seq))


def gc_fraction(seq: str) -> float:
    s = seq.upper()
    if not s:
        return 0.0
    return (s.count("G") + s.count("C")) / len(s)


def wallace_tm(seq: str, na_mM: float = 50.0, mg_mM: float = 1.5) -> float:
    """Salt + Mg2+ adjusted Tm (Wallace rule, Owczarzy correction)."""
    s = seq.upper()
    n = len(s)
    if n == 0:
        return 0.0
    gc_pct = (s.count("G") + s.count("C")) / n * 100
    na_eq = na_mM + 4.0 * math.sqrt(max(mg_mM, 0.0))
    return 81.5 + 16.6 * math.log10(max(na_eq, 0.001) / 1000.0) \
        + 0.41 * gc_pct - 600.0 / n


def nn_delta_g(seq1: str, seq2: str) -> float:
    """NN dG (kcal/mol) for a complementary duplex at 37 C."""
    s1 = seq1.upper()
    s2 = seq2.upper()
    n = len(s1)
    if n < 2:
        return 0.0
    dH, dS = _NN_INIT[0], _NN_INIT[1]
    if s1[0] in "AT" or s2[-1] in "AT":
        dH += _NN_TAT[0]
        dS += _NN_TAT[1]
    if s1[-1] in "AT" or s2[0] in "AT":
        dH += _NN_TAT[0]
        dS += _NN_TAT[1]
    for i in range(n - 1):
        p = s1[i:i + 2] + "/" + s2[i:i + 2]
        dh, ds = _NN_P.get(p, _NN_P.get(p[::-1], (0, 0)))
        dH += dh
        dS += ds
    return dH - 310.15 * dS / 1000.0


def find_dimers(seq1: str, seq2: str, min_bp: int = 3) -> List[Dict]:
    """All dimer alignments sorted by dG (most negative first)."""
    s1 = seq1.upper()
    s2 = seq2.upper()
    rc = revcomp(s2)
    n1, n2 = len(s1), len(s2)
    dimers = []
    for sh in range(-n1 + 1, n2):
        if sh < 0:
            r1s, r1e = -sh, n1
            r2s, r2e = 0, min(n1 + sh, n2)
        else:
            r1s, r1e = 0, min(n1, n2 - sh)
            r2s, r2e = sh, min(sh + n1, n2)
        if r1e - r1s < min_bp:
            continue
        sub1 = s1[r1s:r1e]
        sub2 = rc[r2s:r2e]
        bp = sum(1 for a, b in zip(sub1, sub2) if a == b)
        if bp < min_bp:
            continue
        dimers.append({"dG": nn_delta_g(sub1, sub2), "bp": bp,
                       "aligned1": sub1, "aligned2": sub2})
    dimers.sort(key=lambda x: x["dG"])
    return dimers


def worst_dimer_dg(seq1: str, seq2: str, min_bp: int = 3) -> float:
    """Most negative dimer dG, or 0.0 when nothing pairs."""
    d = find_dimers(seq1, seq2, min_bp)
    return d[0]["dG"] if d else 0.0


def has_dimer(seq1: str, seq2: str, dg_cutoff: float = DIMER_DG) -> bool:
    """NN dG-based dimer check (SpacerFinder logic, unchanged)."""
    s1 = seq1.upper()
    s2 = seq2.upper()
    rc = revcomp(s2)
    n1, n2 = len(s1), len(s2)
    for sh in range(-n1 + 1, n2):
        if sh < 0:
            r1s, r1e = -sh, n1
            r2s, r2e = 0, min(n1 + sh, n2)
        else:
            r1s, r1e = 0, min(n1, n2 - sh)
            r2s, r2e = sh, min(sh + n1, n2)
        if r1e - r1s < 3:
            continue
        sub1 = s1[r1s:r1e]
        sub2 = rc[r2s:r2e]
        if sum(1 for a, b in zip(sub1, sub2) if a == b) < 3:
            continue
        dg = nn_delta_g(sub1, sub2)
        at_3prime = (r1e == n1 or r2e == n2)
        cutoff = DIMER_DG_3P if at_3prime else dg_cutoff
        if dg < cutoff:
            return True
    return False


def hairpin_dg(seq: str, min_stem: int = 3) -> float:
    """Most negative hairpin dG found (0.0 if none)."""
    s = seq.upper()
    n = len(s)
    best = 0.0
    for ls in range(min_stem, max(min_stem + 1, n - min_stem - 3)):
        for le in range(ls + 3, n - min_stem + 1):
            stem_len = 0
            left = ls - 1
            right = le
            dH, dS = _NN_INIT[0], _NN_INIT[1]
            while left >= 0 and right < n and \
                    s[left] == COMP_TABLE.get(s[right], ''):
                stem_len += 1
                left -= 1
                right += 1
            if stem_len >= min_stem:
                dg = dH - 310.15 * dS / 1000.0
                if dg < best:
                    best = dg
    return best


def has_hairpin(seq: str, min_stem: int = 3,
                dg_cutoff: float = HAIRPIN_DG) -> bool:
    """NN hairpin check (SpacerFinder logic, unchanged)."""
    s = seq.upper()
    n = len(s)
    for ls in range(min_stem, max(min_stem + 1, n - min_stem - 3)):
        for le in range(ls + 3, n - min_stem + 1):
            stem_len = 0
            left = ls - 1
            right = le
            dH, dS = _NN_INIT[0], _NN_INIT[1]
            while left >= 0 and right < n and \
                    s[left] == COMP_TABLE.get(s[right], ''):
                stem_len += 1
                left -= 1
                right += 1
            if stem_len >= min_stem:
                dg = dH - 310.15 * dS / 1000.0
                if dg < dg_cutoff:
                    return True
    return False


def end_bonus(seq: str) -> float:
    """3' scoring, inherited from SpacerFinder.

    last base G/C -> +0.5
    exactly 3 G/C within the last 5 bp AND no run of 3 consecutive G/C -> +0.8
    """
    s = seq.upper()
    if not s:
        return 0.0
    bonus = 0.5 if s[-1] in ("G", "C") else 0.0
    tail5 = s[-5:] if len(s) >= 5 else s
    gc_count = sum(1 for b in tail5 if b in ("G", "C"))
    if gc_count == 3:
        has_run3 = any(
            tail5[i] in ("G", "C") and tail5[i + 1] in ("G", "C")
            and tail5[i + 2] in ("G", "C")
            for i in range(len(tail5) - 2)
        )
        if not has_run3:
            bonus += 0.8
    return bonus


def max_homopolymer(seq: str) -> int:
    """Longest run of a single base."""
    s = seq.upper()
    if not s:
        return 0
    best = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    return best


def evaluate_primer(seq: str) -> Dict:
    """Full QC record for one primer, with inherited pass/fail flags."""
    s = seq.upper()
    n = len(s)
    gc = gc_fraction(s)
    tm = wallace_tm(s)
    hp = hairpin_dg(s)
    sd = worst_dimer_dg(s, s)
    homo = max_homopolymer(s)
    flags = []
    if not (LEN_MIN <= n <= LEN_MAX):
        flags.append(f"Len={n}")
    if gc < GC_MIN or gc > GC_MAX:
        flags.append(f"GC={gc * 100:.0f}%")
    if tm < TM_MIN or tm > TM_MAX:
        flags.append(f"Tm={tm:.1f}")
    if has_hairpin(s):
        flags.append(f"Hairpin={hp:.1f}")
    if has_dimer(s, s):
        flags.append(f"SelfDimer={sd:.1f}")
    if homo >= 5:
        flags.append(f"PolyN={homo}")
    return {
        "seq": s, "len": n, "gc": gc, "tm": tm,
        "hairpin_dg": hp, "self_dimer_dg": sd,
        "homopolymer": homo, "end_bonus": end_bonus(s),
        "flags": flags, "pass": len(flags) == 0,
    }
