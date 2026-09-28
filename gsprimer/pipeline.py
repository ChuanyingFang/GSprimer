"""
GSprimer — Orchestration.

Two-stage design so the tool works both as a terminal program and as a
callable step inside an AI agent conversation:

    stage 1  `run_design`   -> enumerate + BLAST + risk-score every candidate,
                               persist state, print the candidate table
    [user picks a pair]
    stage 2  `run_finalize` -> attach GS adapters to the chosen pair(s),
                               rebuild the report with the order sheet

State is cached as JSON next to the report so stage 2 never recomputes BLAST.
"""

import json
import os
import sys
import time
from dataclasses import asdict
from typing import Dict, List, Optional, Tuple

from . import blast as B
from . import risk as R
from .adapters import (ADAPTER_F, ADAPTER_R, order_sheet, ADAPTER_SCHEMES,
                       GG_SITE, GG_ENZYME)
from .config import get
from .design import PrimerPair, design, summarize_isoforms
from .report import generate_report, write_tsv
from .template import (Transcript, load_transcripts, pick_longest,
                       read_fasta)

STATE_NAME = "gsprimer_state.json"


def _log(msg: str) -> None:
    print(msg, file=sys.stderr)


# ============================================================
# Stage 1 — design
# ============================================================

def run_design(query: str, *, genome_fa: str = "", gff3: str = "",
               cdna_fa: str = "", cds_fa: str = "", sequence: str = "",
               transcript: str = "", outdir: str = ".",
               prefix: str = "", max_shift: int = 12, r_ext_max: int = 0,
               top_n: int = 10, allow_downstream: bool = True,
               blastn: str = "", blast_db: str = "",
               fallback_fasta: str = "", do_blast: bool = True,
               adapter_f: str = ADAPTER_F, adapter_r: str = ADAPTER_R,
               enzyme: str = "BsaI", site: str = "GGTCTC",
               no_intron: Optional[bool] = None,
               longest_by: str = "cds") -> Dict:
    """Enumerate candidates, run specificity + risk, write report/TSV/state."""
    t0 = time.time()
    os.makedirs(outdir, exist_ok=True)
    prefix = prefix or (query or "gsprimer").replace(":", "_")

    # ---- template ----
    txs, mode = load_transcripts(
        query, genome_fa=genome_fa, gff3=gff3, cdna_fa=cdna_fa,
        cds_fa=cds_fa, sequence=sequence)
    if not txs:
        return {"ok": False,
                "error": f"Could not obtain a transcript sequence for {query}. "
                         f"Check the gene/transcript ID, GFF3/FASTA paths, or "
                         f"supply the sequence directly via --sequence."}

    if transcript:
        hit = [t for t in txs if t.tx_id == transcript
               or t.tx_id.split()[0] == transcript]
        if not hit:
            return {"ok": False,
                    "error": f"Specified transcript {transcript} not found; "
                             f"available: " + ", ".join(t.tx_id for t in txs)}
        tx = hit[0]
    else:
        tx = pick_longest(txs, by=longest_by)
    _log(f"[Template] {tx.tx_id}  mRNA={len(tx.mrna)}  CDS={tx.cds_len}  "
         f"protein={tx.protein_len}aa  ({mode})")
    issues = tx.validate()
    for m in issues:
        _log(f"  [!] {m}")

    # ---- design ----
    pairs = design(tx, max_shift=max_shift, r_ext_max=r_ext_max,
                   top_n=top_n, allow_downstream=allow_downstream)
    if not pairs:
        return {"ok": False,
                "error": "Could not build any primer pair within the length "
                         "range; try relaxing --max-shift or check the CDS length"}
    _log(f"[Design] {len(pairs)} candidate pairs "
         f"(Tier A={sum(1 for p in pairs if p.tier == 'A')}, "
         f"B={sum(1 for p in pairs if p.tier == 'B')}, "
         f"C={sum(1 for p in pairs if p.tier == 'C')})")

    # ---- specificity for EVERY candidate (spec requirement) ----
    blastn = blastn or get("blastn")
    blast_db = blast_db or get("blast_db")
    fb: Optional[Dict[str, str]] = None
    if fallback_fasta and os.path.exists(fallback_fasta):
        _log(f"[Spec] loaded local FASTA as specificity search DB: "
             f"{fallback_fasta}")
        fb = read_fasta(fallback_fasta)
    if do_blast:
        usable = B.blast_available(blastn, blast_db) or bool(fb)
        if not usable:
            _log("[Spec] no usable BLAST DB found; skipping specificity "
                 "search (report will mark it UNVERIFIED)")
        cache: Dict[Tuple[str, str], Dict] = {}
        for i, p in enumerate(pairs, 1):
            key = (p.f_seq, p.r_seq)
            if key in cache:
                p.blast = cache[key]
                continue
            p.blast = B.check_pair(p.f_seq, p.r_seq, blastn=blastn,
                                   db=blast_db, fallback_db=fb,
                                   expect_len=p.amplicon_len
                                   - (p.f_len - 1) - (p.r_len - 1))
            cache[key] = p.blast
            _log(f"  [{i}/{len(pairs)}] {p.label} -> "
                 f"{p.blast.get('level', 'NA')}")

    # ---- risk ----
    for p in pairs:
        sheet = order_sheet(p, adapter_f, adapter_r)
        p.risk = R.assess(tx, p, adapter_f=adapter_f, adapter_r=adapter_r,
                          full_product=sheet["construct"].full_product,
                          enzyme=enzyme, site=site, no_intron=no_intron)

    isoforms = summarize_isoforms(txs)
    intron_status = (pairs[0].risk.get("amplification", {}).get("has_introns")
                     if pairs else None)
    intron_cn = {True: "has introns (cDNA required)",
                 False: "intronless (gDNA usable)",
                 None: "unknown (no genome annotation in input)"}.get(
                     intron_status, "unknown")
    params = {
        "Query": query, "Template mode": mode, "Chosen transcript": tx.tx_id,
        "Transcript selection rule": f"longest {longest_by.upper()}",
        "Primer length": "19-24 nt (inherited from SpacerFinder)",
        "GC range": "40-60% (inherited)",
        "Tm range": "55-60 °C (Wallace + salt correction, inherited)",
        "ΔTm cap": "< 2 °C (inherited)",
        "Hairpin/dimer ΔG": "≥ -4.0 / ≥ -6.0 kcal·mol⁻¹ (3' end -5.0, inherited)",
        "Amplicon size limit": "none (fixed by full CDS length; the one "
                               "deviation from SpacerFinder)",
        "Forward anchor": f"ATG start, ±{max_shift} nt shift allowed, must be "
                          f"a multiple of 3",
        "Reverse anchor": f"stop codon end, 3' UTR extension cap {r_ext_max} nt",
        "F adapter": adapter_f, "R adapter": adapter_r,
        "Assembly enzyme": f"{enzyme} ({site})",
        "Specificity search": (pairs[0].blast.get("method", "not performed")
                               if pairs and pairs[0].blast else "not performed"),
        "BLAST DB": blast_db or fallback_fasta or "not configured",
        "Intron status": intron_cn,
    }

    html_path = os.path.join(outdir, f"{prefix}_GSprimer.html")
    tsv_path = os.path.join(outdir, f"{prefix}_GSprimer.tsv")
    write_tsv(pairs, tsv_path)
    generate_report(tx=tx, isoforms=isoforms, pairs=pairs, selected=[],
                    params=params, out_path=html_path, gene_query=query,
                    adapter_f=adapter_f, adapter_r=adapter_r)

    state = {
        "query": query, "mode": mode, "tx": _tx_to_dict(tx),
        "isoforms": isoforms, "params": params,
        "pairs": [_pair_to_dict(p) for p in pairs],
        "adapter_f": adapter_f, "adapter_r": adapter_r,
        "enzyme": enzyme, "site": site, "prefix": prefix,
        "html": html_path, "tsv": tsv_path,
    }
    state_path = os.path.join(outdir, STATE_NAME)
    with open(state_path, "w", encoding="utf-8") as f:
        json.dump(state, f, ensure_ascii=False, indent=1)

    _log(f"[Done] {time.time() - t0:.1f}s  report: {html_path}")
    return {"ok": True, "tx": tx, "pairs": pairs, "isoforms": isoforms,
            "params": params, "html": html_path, "tsv": tsv_path,
            "state": state_path, "template_issues": issues}


# ============================================================
# Stage 2 — finalize
# ============================================================

def run_finalize(selection: List[int], *, outdir: str = ".",
                 state_path: str = "") -> Dict:
    """Attach GS adapters to the user-selected pair(s) and rebuild the report."""
    state_path = state_path or os.path.join(outdir, STATE_NAME)
    if not os.path.exists(state_path):
        return {"ok": False,                 "error": f"design state file {state_path} not found; "
                         f"run the design stage first"}
    with open(state_path, encoding="utf-8") as f:
        st = json.load(f)

    tx = _tx_from_dict(st["tx"])
    pairs = [_pair_from_dict(d) for d in st["pairs"]]
    bad = [i for i in selection if i < 1 or i > len(pairs)]
    if bad:
        return {"ok": False,
                "error": f"number {bad} out of range (1-{len(pairs)})"}
    selected = [pairs[i - 1] for i in selection]

    sheets = []
    scheme_sheets: Dict[int, Dict[str, Dict]] = {}
    for p in selected:
        pid = pairs.index(p) + 1
        sheets.append(order_sheet(p, st["adapter_f"], st["adapter_r"]))
        # every registered scheme (seamless + golden_gate) gets its own sheet
        scheme_sheets[pid] = {k: order_sheet(p, s["f"], s["r"])
                              for k, s in ADAPTER_SCHEMES.items()}
        # Golden Gate (PaqCI) internal-site risk on the bare amplicon, so the
        # report can warn per-scheme independently of the seamless (BsaI) check
        gg = R.enzyme_risk(p.amplicon, site=GG_SITE, enzyme=GG_ENZYME)
        p.risk = dict(p.risk or {})
        p.risk["enzyme_gg"] = gg

    html_path = os.path.join(outdir, f"{st['prefix']}_GSprimer_final.html")
    generate_report(tx=tx, isoforms=st["isoforms"], pairs=pairs,
                    selected=selected, params=st["params"],
                    out_path=html_path, gene_query=st["query"],
                    adapter_f=st["adapter_f"], adapter_r=st["adapter_r"],
                    scheme_sheets=scheme_sheets)

    order_path = os.path.join(outdir, f"{st['prefix']}_order.tsv")
    with open(order_path, "w", encoding="utf-8") as f:
        f.write("Name\tScheme\tEnzyme\tSequence_5to3\tLength\tTm_core\t"
                "GC_core\tNote\n")
        for i, p in enumerate(selected, 1):
            ss = scheme_sheets[pairs.index(p) + 1]
            for k, s in ss.items():
                sc = ADAPTER_SCHEMES[k]
                f.write(f"{st['prefix']}-{i}-{sc['key']}-F\t{sc['key']}\t"
                        f"{sc['enzyme']}\t{s['F_full']}\t{s['F_len']}\t"
                        f"{p.f_qc['tm']:.1f}\t{p.f_qc['gc'] * 100:.1f}\t"
                        f"{sc['name']}+"
                        f"{'ATG start' if p.offset == 0 else 'offset%+d' % p.offset}\n")
                f.write(f"{st['prefix']}-{i}-{sc['key']}-R\t{sc['key']}\t"
                        f"{sc['enzyme']}\t{s['R_full']}\t{s['R_len']}\t"
                        f"{p.r_qc['tm']:.1f}\t{p.r_qc['gc'] * 100:.1f}\t"
                        f"{sc['name']}+"
                        f"{'stop codon end' if p.r_ext == 0 else 'stop+%dnt' % p.r_ext}\n")

    _log(f"[Final] report: {html_path}\n[Final] order sheet: {order_path}")
    return {"ok": True, "html": html_path, "order": order_path,
            "selected": selected, "sheets": sheets,
            "scheme_sheets": scheme_sheets, "tx": tx}


# ============================================================
# (de)serialisation
# ============================================================

def _tx_to_dict(tx: Transcript) -> Dict:
    d = asdict(tx)
    d["exons"] = [list(e) for e in tx.exons]
    return d


def _tx_from_dict(d: Dict) -> Transcript:
    t = Transcript(**{k: v for k, v in d.items() if k != "exons"})
    t.exons = [tuple(e) for e in d.get("exons", [])]
    return t


def _pair_to_dict(p: PrimerPair) -> Dict:
    return asdict(p)


def _pair_from_dict(d: Dict) -> PrimerPair:
    return PrimerPair(**d)


# ============================================================
# Console table (shared by CLI and interactive mode)
# ============================================================

def candidate_table(pairs: List[PrimerPair]) -> str:
    hdr = (f"{'#':>3} {'★':1} {'T':1} {'off':>4} {'F primer':<25}{'nt':>3}"
           f"{'Tm':>6}{'GC':>5}  {'R primer':<25}{'nt':>3}{'Tm':>6}{'GC':>5}"
           f"{'dTm':>5}{'bp':>7}{'BsaI':>5}  {'risk':<7} flags")
    lines = [hdr, "-" * len(hdr)]
    for i, p in enumerate(pairs, 1):
        r = p.risk or {}
        enz = r.get("enzyme", {})
        lines.append(
            f"{i:>3} {'*' if p.is_anchor else ' '} {p.tier} {p.offset:>+4} "
            f"{p.f_seq:<25}{p.f_len:>3}{p.f_qc['tm']:>6.1f}"
            f"{p.f_qc['gc'] * 100:>5.0f}  "
            f"{p.r_seq:<25}{p.r_len:>3}{p.r_qc['tm']:>6.1f}"
            f"{p.r_qc['gc'] * 100:>5.0f}"
            f"{p.tm_diff:>5.1f}{p.amplicon_len:>7}"
            f"{enz.get('n_total', '-'):>5}  "
            f"{r.get('overall', '-'):<7} {'; '.join(p.flags)[:46]}")
    return "\n".join(lines)
