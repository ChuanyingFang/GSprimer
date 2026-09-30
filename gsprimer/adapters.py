"""
GSprimer — GS vector adapter handling (GSprimer v2.0, GS series).

Default adapters (user-supplied, GS series):

    Golden Gate (PaqCI / AarI) : agCACCTGCagtcattc / agCACCTGCagtgctc
    Seamless    (BsaI)        : gataagcttGGTCTCTattc / CATggatccGGTCTCAgctc

Anatomy (BsaI = GGTCTC(1/5), leaves a 4-nt 5' overhang):

    gtgatatc  A GGTCTC T cgag  |  <forward primer, starts at ATG>
    ~~~~~~~~    ~~~~~~ ^ ~~~~
    spacer      BsaI   N overhang
    (protects the site from
     end-fraying during digest)

    gccgcggg  T GGTCTC A atcc  |  <reverse primer, ends on the stop codon>

On the sense strand the released insert therefore reads:

    attc [ promoter / ORF ..... ] gagc
    ^^^^                          ^^^^
    left overhang (= F overhang)  right overhang (= revcomp of `gctc`)

Because the adapters supply the ONLY two BsaI sites the assembly needs, any
additional BsaI site inside the amplicon will fragment the insert — that is
why risk.py scans both strands.
"""

import re
from dataclasses import dataclass
from typing import Dict, List, Optional

from .thermo import revcomp

ADAPTER_F = "gataagcttGGTCTCTattc"
ADAPTER_R = "CATggatccGGTCTCAgctc"
BSAI_SITE = "GGTCTC"
BSAI_SPACER = 1          # N nucleotides between site and cut
BSAI_OVERHANG = 4

# Golden Gate adapters (PaqCI / AarI isoschizomer, recognition CACCTGC(4/8)).
# Same 4-nt overhangs (attc / gagc) as the seamless set, so both schemes are
# compatible with the same GS destination vector; only the enzyme differs.
ADAPTER_F_GG = "agCACCTGCagtcattc"
ADAPTER_R_GG = "agCACCTGCagtgctc"
GG_ENZYME = "PaqCI"
GG_SITE = "CACCTGC"


@dataclass
class AdapterInfo:
    seq: str
    site_start: int          # 0-based index of GGTCTC
    spacer: str
    overhang: str
    tail: str                # everything 5' of the BsaI site
    valid: bool
    problems: List[str]


# Type IIS cut parameters per enzyme. `top_cut` = number of nucleotides the
# top strand is cleaved downstream of the recognition site; the 5' overhang is
# the `overhang` nucleotides immediately inside that cut:
#   BsaI  GGTCTC(1/5)  -> overhang = [site_end+1 : site_end+5]
#   PaqCI CACCTGC(4/8) -> overhang = [site_end+4 : site_end+8]  (last 4 nt of
#                          an 8-nt tail; this is what keeps the GG overhangs
#                          attc/gagc identical between the two cloning schemes)
ENZYME_PARAMS = {
    "BsaI":  {"site": BSAI_SITE, "top_cut": 1, "overhang": BSAI_OVERHANG},
    "PaqCI": {"site": GG_SITE,   "top_cut": 4, "overhang": 4},
}


def _detect_enzyme(adapter: str) -> str:
    up = adapter.upper()
    if GG_SITE in up:
        return "PaqCI"
    if BSAI_SITE in up:
        return "BsaI"
    return ""


def parse_adapter(adapter: str, site: str = "", top_cut: int = -1,
                  overhang_len: int = BSAI_OVERHANG) -> AdapterInfo:
    """Locate the type-IIS site and derive the overhang it will generate.

    When `site`/`top_cut` are omitted the enzyme is auto-detected from the
    recognition sequence embedded in the adapter (PaqCI for CACCTGC, BsaI for
    GGTCTC), so the seamless and Golden Gate schemes share one parser and
    always emit the correct 4-nt sticky end.
    """
    s = adapter
    up = s.upper()
    problems: List[str] = []

    if not site:
        enz = _detect_enzyme(s)
        if not enz:
            return AdapterInfo(s, -1, "", "", s, False,
                               ["no BsaI/PaqCI recognition site found in adapter"])
        params = ENZYME_PARAMS[enz]
        site = params["site"]; top_cut = params["top_cut"]; overhang_len = params["overhang"]
    elif top_cut < 0:
        enz = _detect_enzyme(s)
        if enz and ENZYME_PARAMS[enz]["site"] == site:
            top_cut = ENZYME_PARAMS[enz]["top_cut"]
        else:
            top_cut = BSAI_SPACER

    idx = up.find(site)
    if idx < 0:
        return AdapterInfo(s, -1, "", "", s, False,
                           [f"no {site} site found in adapter"])
    if up.count(site) > 1 or up.count(revcomp(site)) > 0:
        problems.append("adapter contains multiple restriction sites")
    site_end = idx + len(site)
    # Everything 3' of the recognition site in the adapter.
    tail3 = s[site_end:]
    if len(tail3) < overhang_len:
        problems.append("adapter 3' end too short to yield a 4-nt overhang")
        return AdapterInfo(s, idx, "", "", s[:idx], False, problems)
    # The overhang is the LAST `overhang_len` nucleotides of the tail after the
    # recognition site. This is correct for both BsaI(1/5) and PaqCI(4/8)
    # adapters (including the GG_R with only a 7-nt tail), and naturally yields
    # the shared GS-vector 4-nt overhang (attc/gagc).
    overhang = tail3[-overhang_len:]
    spacer = tail3[:-overhang_len] if len(tail3) > overhang_len else ""
    if idx < 6:
        problems.append(f"only {idx} nt of 5' protection upstream of the site "
                        f"(recommended >= 6 nt)")
    return AdapterInfo(
        seq=s, site_start=idx,
        spacer=spacer,
        overhang=overhang,
        tail=s[:idx], valid=not problems, problems=problems)


def build_ordered_primer(adapter: str, primer: str) -> str:
    """Adapter (as given, mixed case is meaningful) + primer (upper case)."""
    return adapter + primer.upper()


@dataclass
class Construct:
    f_full: str
    r_full: str
    f_adapter: AdapterInfo
    r_adapter: AdapterInfo
    amplicon: str            # sense strand, no adapters
    full_product: str        # sense strand, with adapters
    insert: str              # sense strand after BsaI digestion
    left_overhang: str
    right_overhang: str
    atg_in_insert: int       # 1-based position of A of ATG (0 = absent)
    frame_ok: bool
    notes: List[str]


def build_construct(amplicon: str, *, adapter_f: str = ADAPTER_F,
                    adapter_r: str = ADAPTER_R,
                    atg_in_amplicon: int = 0) -> Construct:
    """Assemble ordered primers, full PCR product, and the digested insert.

    `atg_in_amplicon` is the 0-based index of the A of ATG within `amplicon`
    (equals 0 for a strict ATG-anchored forward primer, or the length of the
    5'UTR carried along when the forward primer was shifted upstream).
    """
    fa = parse_adapter(adapter_f)
    ra = parse_adapter(adapter_r)
    amp = amplicon.upper()

    f_full = build_ordered_primer(adapter_f, amp[:1] and amp or "")
    # forward ordered primer is adapter + the forward primer, filled by caller
    r_tail_sense = revcomp(adapter_r)          # 5'->3' on the sense strand
    full_product = adapter_f + amp + r_tail_sense

    left_ov = fa.overhang.upper()
    right_ov = revcomp(ra.overhang).upper()
    insert = left_ov + amp + right_ov

    notes: List[str] = []
    notes += [f"F adapter: {p}" for p in fa.problems]
    notes += [f"R adapter: {p}" for p in ra.problems]

    if atg_in_amplicon >= 0 and amp[atg_in_amplicon:atg_in_amplicon + 3] == "ATG":
        atg_pos = len(left_ov) + atg_in_amplicon + 1     # 1-based in insert
    else:
        atg_pos = 0
        notes.append("ATG not located in the insert (forward primer shifted into "
                     "the CDS interior?)")

    frame_ok = (atg_in_amplicon % 3 == 0) if atg_pos else False
    return Construct(f_full=f_full, r_full="", f_adapter=fa, r_adapter=ra,
                     amplicon=amp, full_product=full_product, insert=insert,
                     left_overhang=left_ov, right_overhang=right_ov,
                     atg_in_insert=atg_pos, frame_ok=frame_ok, notes=notes)


def order_sheet(pair, adapter_f: str = ADAPTER_F,
                adapter_r: str = ADAPTER_R) -> Dict:
    """Everything needed to place a synthesis order for one primer pair."""
    f_full = build_ordered_primer(adapter_f, pair.f_seq)
    r_full = build_ordered_primer(adapter_r, pair.r_seq)
    atg_in_amp = max(0, -pair.offset) if pair.offset <= 0 else -1
    con = build_construct(pair.amplicon, adapter_f=adapter_f,
                          adapter_r=adapter_r, atg_in_amplicon=atg_in_amp)
    con.f_full = f_full
    con.r_full = r_full
    return {
        "F_name": "GS-F", "R_name": "GS-R",
        "F_full": f_full, "R_full": r_full,
        "F_len": len(f_full), "R_len": len(r_full),
        "F_core": pair.f_seq, "R_core": pair.r_seq,
        "adapter_f": adapter_f, "adapter_r": adapter_r,
        "left_overhang": con.left_overhang,
        "right_overhang": con.right_overhang,
        "insert_len": len(con.insert),
        "product_len": len(con.full_product),
        "atg_in_insert": con.atg_in_insert,
        "frame_ok": con.frame_ok,
        "construct": con,
    }


# ------------------------------------------------------------
# Two cloning schemes share the same GS destination vector:
#   seamless      -> BsaI,     overhangs attc/gagc
#   golden_gate   -> PaqCI,    overhangs attc/gagc (identical 4-nt ends)
# Both are offered to the user at finalize time so they can pick.
# ------------------------------------------------------------
ADAPTER_SCHEMES = {
    "seamless": {
        "key": "seamless",
        "name": "Scarless cloning adapter (BsaI)",
        "enzyme": "BsaI", "site": BSAI_SITE,
        "f": ADAPTER_F, "r": ADAPTER_R,
    },
    "golden_gate": {
        "key": "golden_gate",
        "name": "Golden Gate adapter (PaqCI / AarI)",
        "enzyme": GG_ENZYME, "site": GG_SITE,
        "f": ADAPTER_F_GG, "r": ADAPTER_R_GG,
    },
}
SCHEME_ORDER = ["seamless", "golden_gate"]


def order_sheets_for(pair) -> Dict[str, Dict]:
    """Build an order sheet for every registered adapter scheme."""
    return {k: order_sheet(pair, s["f"], s["r"])
            for k, s in ADAPTER_SCHEMES.items()}
