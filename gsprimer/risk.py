"""
GSprimer — Risk assessment for GS cloning amplicons.

Three mandatory risk blocks (per skill spec):

  1. primer specificity      -> from blast.py (BLAST + in-silico PCR)
  2. amplicon GC content     -> global + sliding-window + primer-proximal
  3. internal type-IIS sites -> BOTH strands, scanned separately for the
                               seamless (BsaI/GGTCTC) and Golden Gate
                               (PaqCI/CACCTGC) schemes, with synonymous
                               mutation rescue when a site falls inside the CDS

Plus two practical extras that bite in real labs:
  4. template choice (cDNA vs gDNA) — the template here is a transcript, so
     genomic DNA would carry introns
  5. amplification difficulty (length, homopolymer, repeats)
"""

import re
from typing import Dict, List, Optional

from .thermo import gc_fraction, max_homopolymer, revcomp, TM_DIFF_MAX

BSAI = "GGTCTC"
BSAI_RC = "GAGACC"

GC_GLOBAL_WARN = (0.35, 0.65)
GC_WINDOW = 50
GC_WINDOW_STEP = 10
GC_WINDOW_WARN = 0.75
GC_WINDOW_LOW = 0.25

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(seq: str) -> str:
    s = seq.upper()
    return "".join(CODON_TABLE.get(s[i:i + 3], "X")
                   for i in range(0, len(s) - 2, 3))


# ============================================================
# 1. Internal type-IIS sites
# ============================================================

def scan_sites(seq: str, site: str = BSAI) -> List[Dict]:
    """Find `site` on BOTH strands of `seq`. Positions are 0-based on the
    sense strand; strand '-' means the site reads 5'->3' on the antisense."""
    s = seq.upper()
    rc_site = revcomp(site)
    out = []
    for m in re.finditer(f"(?={re.escape(site)})", s):
        out.append({"pos": m.start(), "strand": "+", "match": site})
    for m in re.finditer(f"(?={re.escape(rc_site)})", s):
        out.append({"pos": m.start(), "strand": "-", "match": rc_site})
    out.sort(key=lambda x: x["pos"])
    return out


def _synonymous_fixes(amplicon: str, site_pos: int, site_len: int,
                      cds_offset: int, cds_len: int,
                      site: str = BSAI, rc_site: str = BSAI_RC) -> List[Dict]:
    """Propose single-base synonymous changes that destroy the site."""
    fixes: List[Dict] = []
    s = list(amplicon.upper())
    cds_stop = cds_offset + cds_len
    for i in range(site_pos, site_pos + site_len):
        if i < cds_offset or i >= cds_stop:
            continue
        codon_start = cds_offset + ((i - cds_offset) // 3) * 3
        if codon_start + 3 > len(s):
            continue
        codon = "".join(s[codon_start:codon_start + 3])
        aa = CODON_TABLE.get(codon)
        if aa is None:
            continue
        p_in_codon = i - codon_start
        for alt in "ACGT":
            if alt == codon[p_in_codon]:
                continue
            new_codon = codon[:p_in_codon] + alt + codon[p_in_codon + 1:]
            if CODON_TABLE.get(new_codon) != aa:
                continue
            trial = amplicon[:i] + alt + amplicon[i + 1:]
            window = trial[max(0, site_pos - 6):site_pos + site_len + 6]
            if site in window.upper() or rc_site in window.upper():
                if trial[site_pos:site_pos + site_len].upper() in (site, rc_site):
                    continue
            fixes.append({
                "amplicon_pos": i,                       # 0-based
                "codon_no": (i - cds_offset) // 3 + 1,   # 1-based aa number
                "codon": codon, "new_codon": new_codon,
                "aa": aa, "from": codon[p_in_codon], "to": alt,
                "wobble": p_in_codon == 2,
            })
    fixes.sort(key=lambda x: (not x["wobble"], x["codon_no"]))
    return fixes


def enzyme_risk(amplicon: str, *, cds_offset: int = 0, cds_len: int = 0,
                site: str = BSAI, enzyme: str = "BsaI") -> Dict:
    """Internal type-IIS sites inside the amplicon (adapters excluded)."""
    sites = scan_sites(amplicon, site)
    rc_site = revcomp(site)
    cds_len = cds_len or (len(amplicon) - cds_offset)
    detail = []
    for st in sites:
        fixes = _synonymous_fixes(amplicon, st["pos"], len(site),
                                  cds_offset, cds_len, site, rc_site)
        in_cds = cds_offset <= st["pos"] < cds_offset + cds_len
        detail.append({
            **st,
            "pos_1based": st["pos"] + 1,
            "in_cds": in_cds,
            "codon_no": ((st["pos"] - cds_offset) // 3 + 1) if in_cds else 0,
            "context": amplicon[max(0, st["pos"] - 6):
                                st["pos"] + len(site) + 6].upper(),
            "fixes": fixes[:4],
            "rescuable": bool(fixes),
        })
    n_plus = sum(1 for d in detail if d["strand"] == "+")
    n_minus = sum(1 for d in detail if d["strand"] == "-")
    n = len(detail)
    if n == 0:
        level = "LOW"
        msg = f"no internal {enzyme} site in the amplicon; safe for Golden Gate"
    else:
        unrescuable = [d for d in detail if not d["rescuable"]]
        level = "HIGH"
        msg = (f"{n} internal {enzyme} site(s) found in the amplicon "
               f"(sense {n_plus} + antisense {n_minus}); Golden Gate would "
               f"fragment the insert")
        if unrescuable:
            msg += (f"; {len(unrescuable)} of them cannot be removed by a "
                    f"single-base synonymous change")
    return {"enzyme": enzyme, "site": site, "n_total": n,
            "n_plus": n_plus, "n_minus": n_minus, "sites": detail,
            "level": level, "message": msg}


def construct_site_audit(full_product: str, site: str = BSAI,
                         enzyme: str = "BsaI") -> Dict:
    """Sanity check on the complete PCR product: exactly 2 sites expected."""
    sites = scan_sites(full_product, site)
    return {"n_total": len(sites), "expected": 2, "sites": sites,
            "ok": len(sites) == 2,
            "message": (f"complete PCR product contains {len(sites)} "
                        f"{enzyme} site(s) (expected 2, from the two adapters)")}


# ============================================================
# 2. GC content
# ============================================================

def gc_risk(amplicon: str, window: int = GC_WINDOW,
            step: int = GC_WINDOW_STEP) -> Dict:
    s = amplicon.upper()
    n = len(s)
    overall = gc_fraction(s)
    windows = []
    if n >= window:
        for i in range(0, n - window + 1, step):
            windows.append((i, gc_fraction(s[i:i + window])))
    else:
        windows.append((0, overall))
    gcs = [g for _, g in windows]
    wmax = max(windows, key=lambda x: x[1])
    wmin = min(windows, key=lambda x: x[1])

    head = gc_fraction(s[:min(100, n)])
    tail = gc_fraction(s[max(0, n - 100):])

    issues = []
    level = "LOW"
    if overall > GC_GLOBAL_WARN[1]:
        issues.append(f"Overall GC {overall * 100:.1f}% is high — add DMSO "
                      "(3–5%) or betaine (1 M) and raise denaturation "
                      "temperature")
        level = "MEDIUM"
    elif overall < GC_GLOBAL_WARN[0]:
        issues.append(f"Overall GC {overall * 100:.1f}% is low — lower the "
                      "annealing temperature and shorten denaturation time")
        level = "MEDIUM"
    if wmax[1] >= GC_WINDOW_WARN:
        issues.append(f"Local GC-rich region: {wmax[0] + 1}-"
                      f"{wmax[0] + window} nt, GC {wmax[1] * 100:.0f}% — a "
                      f"common PCR stall point")
        level = "HIGH" if wmax[1] >= 0.80 else "MEDIUM"
    if wmin[1] <= GC_WINDOW_LOW:
        issues.append(f"Local GC-poor region: {wmin[0] + 1}-"
                      f"{wmin[0] + window} nt, GC {wmin[1] * 100:.0f}% — "
                      f"prone to non-specific annealing")
        level = max(level, "MEDIUM", key=lambda x: ["LOW", "MEDIUM", "HIGH"].index(x))
    if not issues:
        issues.append(f"GC {overall * 100:.1f}%, no extreme regions across "
                      f"the fragment")
    return {"overall": overall, "head100": head, "tail100": tail,
            "window": window, "windows": windows,
            "max_window": {"start": wmax[0] + 1, "gc": wmax[1]},
            "min_window": {"start": wmin[0] + 1, "gc": wmin[1]},
            "mean_window": sum(gcs) / len(gcs),
            "level": level, "issues": issues}


# ============================================================
# 3. Amplification difficulty / template choice
# ============================================================

def amplification_risk(tx, pair, *, no_intron: Optional[bool] = None) -> Dict:
    """Template / intron status and length related warnings.

    `has_introns`
      True  : gene has introns -> cDNA/RT-PCR mandatory
      False : single-exon gene -> gDNA is a valid, simpler alternative
      None  : unknown (CDS/cDNA/direct input, no genome annotation)
    Pass `no_intron=True` to force the single-exon (gDNA-friendly) branch
    for FASTA / direct inputs where exon structure is not in the data.
    """
    issues = []
    level = "LOW"
    n = pair.amplicon_len

    # ---- exon / intron status ----
    if no_intron is not None:
        has_introns = (not no_intron)
    elif tx.exons:
        has_introns = len(tx.exons) > 1
    else:
        has_introns = None
    n_exon = len(tx.exons)

    # ---- length feasibility ----
    if n > 3000:
        issues.append(f"Amplicon {n} bp — needs a long-fragment high-fidelity "
                      "polymerase (e.g. KOD FX / PrimeSTAR GXL) with "
                      "~1 min/kb extension")
        level = "MEDIUM"
    elif n > 2000:
        issues.append(f"Amplicon {n} bp — use a high-fidelity polymerase and "
                      "extend the extension time")
    if n < 150:
        issues.append(f"Amplicon is only {n} bp — keep it distinct from primer "
                      "dimers")
        level = "MEDIUM"

    # ---- template / intron choice ----
    if has_introns is True:
        genomic_span = 0
        if tx.exons:
            coords = [c for iv in tx.exons for c in iv]
            genomic_span = max(coords) - min(coords) + 1
        issues.append(
            f"Template is a transcript ({n_exon} exons, genomic span "
            f"~{genomic_span} bp): a cDNA / RT-PCR product must be used as "
            f"template; genomic DNA would pull in introns and the product "
            f"would no longer be the target CDS")
        level = "MEDIUM"
        template_advice = {
            "recommendation": "cDNA (RT-PCR product) required",
            "reason": f"This gene has {n_exon} exons (genomic span "
                     f"~{genomic_span} bp); gDNA would co-amplify the introns "
                     f"and the product would not match the CDS",
            "has_introns": True,
        }
    elif has_introns is False:
        issues.append("This gene has no introns: cDNA and genomic DNA give "
                      "identical products, either is fine")
        template_advice = {
            "recommendation": "gDNA (genomic DNA) recommended",
            "reason": "With no introns, gDNA and cDNA products are identical; "
                      "gDNA is simpler — it skips reverse transcription, avoids "
                      "RT-enzyme bias and RNA degradation, and can be PCR'd "
                      "directly",
            "has_introns": False,
        }
    else:
        issues.append("Input provides no genome annotation; cannot determine "
                      "whether the gene has introns")
        template_advice = {
            "recommendation": "Template type pending (cDNA or gDNA)",
            "reason": "Current input (CDS/cDNA/direct sequence) cannot resolve "
                      "introns; if the gene is known to be intronless, use gDNA "
                      "directly to skip RT",
            "has_introns": None,
        }

    # ---- homopolymer in primers ----
    homo_f = max_homopolymer(pair.f_seq)
    homo_r = max_homopolymer(pair.r_seq)
    if max(homo_f, homo_r) >= 5:
        issues.append(f"Primer has a {max(homo_f, homo_r)}-nt homopolymer run, "
                      "which may cause slippage")
        level = "MEDIUM"

    # ---- homopolymer inside amplicon ----
    amp_homo = max_homopolymer(pair.amplicon)
    if amp_homo >= 10:
        issues.append(f"A {amp_homo}-nt homopolymer run inside the amplicon will "
                      "add read-length noise in sequencing")

    return {"level": level, "issues": issues, "amplicon_len": n,
            "n_exons": n_exon, "has_introns": has_introns,
            "template_advice": template_advice}


# ============================================================
# 4. PCR protocol advice (primer + amplicon characteristics)
# ============================================================

GC_HIGH = 0.65          # overall GC considered "high"
GC_WIN_VERYHIGH = 0.80  # local window GC that stalls PCR
GC_LOW = 0.35
LEN_LONG = 3000
LEN_MED = 2000
LEN_SHORT = 150


def pcr_advice(*, gc_overall: float, max_win_gc: float, min_win_gc: float,
               amplicon_len: int, tm_f: float, tm_r: float, tm_diff: float,
               amp_homopolymer: int, has_introns: Optional[bool]) -> Dict:
    """Concrete, actionable PCR protocol tips derived from the primer pair
    and amplicon.  Driven by: GC (global + local), fragment length, primer
    Tm mismatch, and intra-fragment homopolymer runs.

    The template choice (intron-aware) is handled separately in
    `amplification_risk` so it is not duplicated here.
    """
    tips: List[str] = []
    level = "LOW"

    def bump(lv: str) -> None:
        nonlocal level
        level = _worst(level, lv)

    # ---- high / low GC ----
    if max_win_gc >= GC_WIN_VERYHIGH:
        tips.append(
            f"Local GC up to {max_win_gc * 100:.0f}% — a classic high-GC stall "
            f"zone: use a GC-enhanced high-fidelity polymerase "
            f"(Q5® / Phusion GC / KAPA HiFi GC), add DMSO (3–5%) or betaine "
            f"(1 M), and apply touchdown PCR (65→55°C) to suppress "
            f"non-specific products and improve yield")
        bump("HIGH")
    if gc_overall >= GC_HIGH:
        tips.append(
            f"Overall GC {gc_overall * 100:.0f}% is high — use a high-fidelity "
            f"polymerase with 3% DMSO, set the annealing temperature at the "
            f"lower-Tm primer, and consider a GC-rich buffer system")
        bump("MEDIUM")
    if gc_overall <= GC_LOW or min_win_gc <= 0.25:
        tips.append(
            f"AT-rich / low GC (overall {gc_overall * 100:.0f}%): moderately "
            f"lower the annealing temperature and shorten denaturation; watch "
            f"for primer mis-priming and non-specific amplification")
        bump("MEDIUM")

    # ---- fragment length ----
    if amplicon_len > LEN_LONG:
        tips.append(
            f"Amplicon {amplicon_len} bp is long — use a long-fragment "
            f"high-fidelity polymerase (KOD FX neo / PrimeSTAR GXL), "
            f"45–60 s/kb extension, 25–30 cycles")
        bump("MEDIUM")
    elif amplicon_len > LEN_MED:
        tips.append(
            f"Fragment {amplicon_len} bp: high-fidelity polymerase with "
            f"moderately extended extension (~30 s/kb)")
    elif amplicon_len < LEN_SHORT:
        tips.append(
            f"Fragment is only {amplicon_len} bp — keep it distinct from primer "
            f"dimers; use 2–3% agarose gel for resolution")

    # ---- primer Tm mismatch ----
    if tm_diff >= 3.0:
        tips.append(
            f"Forward/reverse primer Tm gap {tm_diff:.1f}°C is large — set "
            f"annealing at the lower Tm, or use a 3-step protocol "
            f"(denature-anneal-extend) with extended annealing for the "
            f"high-Tm primer")
    elif tm_diff >= TM_DIFF_MAX:
        tips.append(
            f"ΔTm {tm_diff:.1f}°C is near the cap — set annealing at the mean "
            f"of the two primer Tms and validate with a gradient")

    # ---- homopolymer inside amplicon ----
    if amp_homopolymer >= 10:
        tips.append(
            f"A {amp_homopolymer}-nt homopolymer run inside the fragment can "
            f"slip during amplification/sequencing: design primers to flank it, "
            f"or switch to a slip-resistant polymerase")
        bump("MEDIUM")

    if not tips:
        tips.append("Fragment and primer traits are routine; standard "
                    "high-fidelity PCR (Phusion / Q5) with standard annealing "
                    "temperature suffices")
    return {"level": level, "tips": tips}


# ============================================================
# Aggregation
# ============================================================

_LEVEL_ORDER = ["LOW", "LOW-MED", "MEDIUM", "HIGH", "UNKNOWN"]


def _worst(*levels: str) -> str:
    known = [l for l in levels if l in _LEVEL_ORDER]
    if not known:
        return "UNKNOWN"
    if "UNKNOWN" in known and len(set(known)) == 1:
        return "UNKNOWN"
    ranked = sorted((l for l in known if l != "UNKNOWN"),
                    key=_LEVEL_ORDER.index)
    return ranked[-1] if ranked else "UNKNOWN"


def assess(tx, pair, *, adapter_f: str = "", adapter_r: str = "",
           full_product: str = "", enzyme: str = "BsaI",
           site: str = BSAI, no_intron: Optional[bool] = None) -> Dict:
    """Build the complete risk record attached to a primer pair."""
    cds_offset = max(0, -pair.offset)
    cds_len = min(tx.cds_len, pair.amplicon_len - cds_offset)

    enz = enzyme_risk(pair.amplicon, cds_offset=cds_offset,
                      cds_len=cds_len, site=site, enzyme=enzyme)
    gc = gc_risk(pair.amplicon)
    amp = amplification_risk(tx, pair, no_intron=no_intron)
    pcr = pcr_advice(gc_overall=gc["overall"],
                     max_win_gc=gc["max_window"]["gc"],
                     min_win_gc=gc["min_window"]["gc"],
                     amplicon_len=pair.amplicon_len,
                     tm_f=pair.f_qc["tm"], tm_r=pair.r_qc["tm"],
                     tm_diff=pair.tm_diff,
                     amp_homopolymer=max_homopolymer(pair.amplicon),
                     has_introns=amp["has_introns"])
    audit = (construct_site_audit(full_product, site, enzyme)
             if full_product else {})

    spec = pair.blast or {}
    spec_level = spec.get("level", "UNKNOWN")

    # frame / N-terminus integrity
    frame_issues = []
    frame_level = "LOW"
    if pair.offset % 3 != 0:
        frame_issues.append("Fatal: forward primer shift is not a multiple of 3; "
                            "reading frame is broken")
        frame_level = "HIGH"
    if pair.offset > 0:
        frame_issues.append(f"Forward primer lies {pair.offset} nt downstream of "
                            f"ATG, losing the N-terminal "
                            f"{pair.offset // 3} aa (including the start Met); "
                            f"only suitable for signal-peptide-truncated / "
                            f"truncation constructs")
        frame_level = "HIGH"
    if pair.offset < 0:
        utr = pair.amplicon[:cds_offset]
        stops = [i // 3 + 1 for i in range(0, len(utr) - 2, 3)
                 if utr[i:i + 3] in ("TAA", "TAG", "TGA")]
        frame_issues.append(f"Forward primer lies {-pair.offset} nt upstream of "
                            f"ATG; the insert will carry {-pair.offset} nt of "
                            f"5'UTR")
        if stops:
            frame_issues.append("Warning: an in-frame stop codon exists in the "
                                "carried 5'UTR; an N-terminal fusion tag would "
                                "not be translated through")
            frame_level = "HIGH"
        else:
            frame_level = "MEDIUM"
    if pair.r_ext > 0:
        frame_issues.append(f"Reverse primer extends {pair.r_ext} nt past the "
                            f"stop codon; the insert end will carry extra 3'UTR "
                            f"bases, breaking a C-terminal fusion")
        frame_level = "HIGH"
    if not frame_issues:
        frame_issues.append("Strict ATG..STOP anchor; frame matches the vector "
                            "overhang")

    overall = _worst(enz["level"], gc["level"], amp["level"], pcr["level"],
                     spec_level, frame_level)
    return {
        "enzyme": enz, "gc": gc, "amplification": amp, "pcr": pcr,
        "construct_audit": audit,
        "specificity": {"level": spec_level, "detail": spec},
        "frame": {"level": frame_level, "issues": frame_issues},
        "overall": overall,
    }
