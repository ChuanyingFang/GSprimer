"""
GSprimer — Report generation (HTML + TSV).

Colour scheme follows the lab's Fv1.4 palette:
  OE red  #B2182B  (high risk / emphasis)
  CR blue #2166AC  (good / primary)
  grey    #999999  (neutral / unknown)
"""

import html
import os
import time
from typing import Dict, List, Optional

from .adapters import order_sheet, ADAPTER_SCHEMES
from .risk import _worst
from .thermo import revcomp

C_RED = "#B2182B"
C_BLUE = "#2166AC"
C_GREY = "#999999"
C_ORANGE = "#E08214"
C_LBLUE = "#92C5DE"

LEVEL_COLOR = {"LOW": C_BLUE, "LOW-MED": C_LBLUE, "MEDIUM": C_ORANGE,
               "HIGH": C_RED, "UNKNOWN": C_GREY}
LEVEL_EN = {"LOW": "Low", "LOW-MED": "Low-Med", "MEDIUM": "Medium",
            "HIGH": "High", "UNKNOWN": "Unknown"}


def _e(x) -> str:
    return html.escape(str(x))


def _badge(level: str) -> str:
    c = LEVEL_COLOR.get(level, C_GREY)
    return (f'<span class="badge" style="background:{c}">'
            f'{LEVEL_EN.get(level, level)}</span>')


# ============================================================
# TSV
# ============================================================

TSV_COLS = ["Rank", "Anchor", "Tier", "Score", "Offset", "F_seq", "F_len",
            "F_Tm", "F_GC", "R_seq", "R_len", "R_Tm", "R_GC", "dTm",
            "HeteroDimer_dG", "Amplicon_bp", "BsaI_internal", "BsaI_plus",
            "BsaI_minus", "PaqCI_internal", "PaqCI_plus", "PaqCI_minus",
            "Amplicon_GC", "MaxWin_GC", "Spec_level",
            "Spec_offtarget_products", "Template_advice", "PCR_advice",
            "Overall_risk", "Flags"]


def write_tsv(pairs: List, path: str) -> str:
    lines = ["\t".join(TSV_COLS)]
    for i, p in enumerate(pairs, 1):
        r = p.risk or {}
        enz = r.get("enzyme", {}) or {}
        enz_gg = r.get("enzyme_gg", {}) or {}
        gc = r.get("gc", {})
        spec = r.get("specificity", {})
        pcr = (spec.get("detail") or {}).get("pcr", {})
        amp = r.get("amplification", {})
        tpl_rec = (amp.get("template_advice", {}) or {}).get("recommendation", "NA")
        pcr_tips = " | ".join((r.get("pcr", {}) or {}).get("tips", [])) or "NA"
        lines.append("\t".join(str(x) for x in [
            i, "YES" if p.is_anchor else "no", p.tier, f"{p.score:.2f}",
            p.offset, p.f_seq, p.f_len, f"{p.f_qc['tm']:.1f}",
            f"{p.f_qc['gc'] * 100:.1f}", p.r_seq, p.r_len,
            f"{p.r_qc['tm']:.1f}", f"{p.r_qc['gc'] * 100:.1f}",
            f"{p.tm_diff:.2f}", f"{p.hetero_dimer_dg:.1f}", p.amplicon_len,
            enz.get("n_total", "NA"), enz.get("n_plus", "NA"),
            enz.get("n_minus", "NA"),
            enz_gg.get("n_total", "NA"), enz_gg.get("n_plus", "NA"),
            enz_gg.get("n_minus", "NA"),
            f"{gc.get('overall', 0) * 100:.1f}" if gc else "NA",
            f"{gc.get('max_window', {}).get('gc', 0) * 100:.0f}" if gc else "NA",
            spec.get("level", "NA"), pcr.get("n_off", "NA"),
            tpl_rec, pcr_tips,
            r.get("overall", "NA"), ";".join(p.flags) or "OK",
        ]))
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return path


# ============================================================
# HTML fragments
# ============================================================

def _gc_svg(gc_block: Dict, width: int = 900, height: int = 170) -> str:
    wins = gc_block.get("windows", [])
    if not wins:
        return ""
    xs = [w[0] for w in wins]
    ys = [w[1] for w in wins]
    x0, x1 = min(xs), max(xs) or 1
    pad_l, pad_b, pad_t = 46, 26, 12
    pw = width - pad_l - 14
    ph = height - pad_b - pad_t

    def px(x):
        return pad_l + (x - x0) / max(1, (x1 - x0)) * pw

    def py(y):
        return pad_t + (1 - y) * ph

    pts = " ".join(f"{px(x):.1f},{py(y):.1f}" for x, y in wins)
    grid = []
    for lvl, lab in [(0.25, "25%"), (0.5, "50%"), (0.75, "75%")]:
        yy = py(lvl)
        col = C_RED if lvl == 0.75 else C_GREY
        dash = "4,3" if lvl == 0.75 else "2,4"
        grid.append(f'<line x1="{pad_l}" y1="{yy:.1f}" x2="{width - 14}" '
                    f'y2="{yy:.1f}" stroke="{col}" stroke-width="1" '
                    f'stroke-dasharray="{dash}" opacity="0.7"/>')
        grid.append(f'<text x="{pad_l - 6}" y="{yy + 4:.1f}" '
                    f'text-anchor="end" font-size="10" fill="#555">{lab}</text>')
    hi = [f'<circle cx="{px(x):.1f}" cy="{py(y):.1f}" r="3" fill="{C_RED}"/>'
          for x, y in wins if y >= 0.75]
    return f'''<svg viewBox="0 0 {width} {height}" class="gcplot">
  <rect x="{pad_l}" y="{pad_t}" width="{pw}" height="{ph}" fill="#fafafa" stroke="#e3e3e3"/>
  {"".join(grid)}
  <polyline points="{pts}" fill="none" stroke="{C_BLUE}" stroke-width="1.8"/>
  {"".join(hi)}
  <text x="{pad_l}" y="{height - 8}" font-size="10" fill="#555">1</text>
  <text x="{width - 14}" y="{height - 8}" text-anchor="end" font-size="10"
        fill="#555">{x1 + gc_block.get("window", 50)} nt</text>
  <text x="10" y="{pad_t + ph / 2}" font-size="10" fill="#555"
        transform="rotate(-90 10 {pad_t + ph / 2})" text-anchor="middle">GC</text>
</svg>'''


def _construct_diagram(sheet: Dict, pair) -> str:
    con = sheet["construct"]
    fa, ra = con.f_adapter, con.r_adapter
    amp_show = 18
    amp = con.amplicon
    amp_head = amp[:amp_show]
    amp_tail = amp[-amp_show:]
    return f'''<div class="construct">
  <div class="crow">
    <span class="seg tail" title="protection bases">{_e(fa.tail)}</span><span
      class="seg site" title="type-IIS recognition site">{_e(fa.seq[fa.site_start:fa.site_start + 6])}</span><span
      class="seg spacer" title="spacer N1">{_e(fa.spacer)}</span><span
      class="seg ovh" title="left 4-nt overhang">{_e(fa.overhang)}</span><span
      class="seg core" title="forward primer core (template sequence)">{_e(amp_head)}</span><span
      class="seg dots">… {pair.amplicon_len} bp …</span><span
      class="seg core">{_e(amp_tail)}</span><span
      class="seg ovh" title="right 4-nt overhang">{_e(con.right_overhang.lower())}</span><span
      class="seg spacer">{_e(revcomp(ra.spacer))}</span><span
      class="seg site">{_e(revcomp(ra.seq[ra.site_start:ra.site_start + 6]))}</span><span
      class="seg tail">{_e(revcomp(ra.tail))}</span>
  </div>
  <div class="legend">
    <span><i style="background:#eeeeee"></i>protection</span>
    <span><i style="background:#f6c6c9"></i>type-IIS site</span>
    <span><i style="background:#fde7c3"></i>N1 spacer</span>
    <span><i style="background:#cfe3f3"></i>4-nt overhang</span>
    <span><i style="background:#dff0d8"></i>template</span>
  </div>
  <div class="note">Insert after enzyme digestion: <code>{_e(con.left_overhang.lower())}</code>
    + template {pair.amplicon_len} bp + <code>{_e(con.right_overhang.lower())}</code>
    = <b>{sheet['insert_len']} bp</b>; ATG at position
    <b>{sheet['atg_in_insert'] or 'NA'}</b> of the insert; reading frame
    {'<b style="color:%s">consistent</b>' % C_BLUE if sheet['frame_ok'] else '<b style="color:%s">inconsistent</b>' % C_RED}
  </div>
</div>'''


def _amp_body(amp: Dict, risk: Dict) -> str:
    """Render risk block 5: feasibility issues + intron-aware template choice
    + concrete PCR protocol advice."""
    amp_issues = amp.get("issues", [])
    tpl = amp.get("template_advice", {})
    pcr = (risk.get("pcr", {}) or {})
    pcr_tips = pcr.get("tips", [])

    tpl_html = ""
    if tpl:
        rec_color = (C_BLUE if tpl.get("has_introns") is False
                     else C_RED if tpl.get("has_introns") is True else "#333")
        tpl_html = (f'<div class="sub2"><b>Template choice:</b>'
                    f'<span style="color:{rec_color};font-weight:600">'
                    f'{_e(tpl.get("recommendation", ""))}</span>'
                    f'<span style="color:#666"> — '
                    f'{_e(tpl.get("reason", ""))}</span></div>')

    pcr_html = ("<ul>" + "".join(f"<li>{_e(t)}</li>" for t in pcr_tips)
                + "</ul>") if pcr_tips else ""

    return ("<ul>" + "".join(f"<li>{_e(i)}</li>" for i in amp_issues) + "</ul>"
            + tpl_html
            + '<div style="margin-top:10px"><b>PCR recipe advice</b></div>'
            + (pcr_html or "<ul><li>Fragment and primer traits are routine; "
                          "standard high-fidelity PCR suffices.</li></ul>"))


def _enz_section(enz: Dict) -> str:
    """Render the internal type-IIS site table for one enzyme."""
    if not enz:
        return ""
    rows = ""
    for s in enz.get("sites", []):
        fixes = "；".join(
            f"{f['aa']}{f['codon_no']} {f['codon']}→{f['new_codon']} "
            f"({f['from']}{f['amplicon_pos'] + 1}{f['to']})"
            for f in s.get("fixes", [])[:2]) or \
            '<span style="color:%s">no single-base synonymous fix</span>' % C_RED
        rows += (f"<tr><td>{s['pos_1based']}</td><td>{s['strand']} strand</td>"
                 f"<td><code>{_e(s['context'])}</code></td>"
                 f"<td>{'in CDS, codon %d' % s['codon_no'] if s['in_cds'] else 'outside CDS'}</td>"
                 f"<td>{fixes}</td></tr>")
    body = f'<p>{_e(enz.get("message", ""))}</p>'
    if rows:
        body += (f'<table class="mini"><tr><th>Position (nt)</th><th>Strand</th>'
                 f'<th>Context</th><th>Locus attribute</th>'
                 f'<th>Synonymous-mutation rescue</th></tr>{rows}</table>')
    audit = enz.get("construct_audit") or {}
    if audit:
        col = C_BLUE if audit.get("ok") else C_RED
        body += (f'<p style="color:{col}">{_e(audit.get("message", ""))}'
                 f'{"" if audit.get("ok") else " — assembly will fail; must be fixed"}</p>')
    return body


def _risk_block(pair) -> str:
    r = pair.risk or {}
    if not r:
        return ""
    enz = r.get("enzyme", {}) or {}
    enz_gg = r.get("enzyme_gg", {}) or {}
    gc = r.get("gc", {})
    spec = r.get("specificity", {})
    det = spec.get("detail") or {}
    pcr = det.get("pcr", {})
    amp = r.get("amplification", {})
    frame = r.get("frame", {})

    # --- specificity ---
    if det.get("available"):
        f_cls, r_cls = det.get("f", {}), det.get("r", {})
        spec_body = f'''<table class="mini">
      <tr><th></th><th>Exact match</th><th>3'-extendable near-match</th><th>Weak hit</th><th>Verdict</th></tr>
      <tr><td>Forward primer</td><td>{f_cls.get('n_perfect', '-')}</td>
          <td>{f_cls.get('n_strong', '-')}</td><td>{f_cls.get('n_weak', '-')}</td>
          <td>{_badge(det.get('f_risk', 'UNKNOWN'))}</td></tr>
      <tr><td>Reverse primer</td><td>{r_cls.get('n_perfect', '-')}</td>
          <td>{r_cls.get('n_strong', '-')}</td><td>{r_cls.get('n_weak', '-')}</td>
          <td>{_badge(det.get('r_risk', 'UNKNOWN'))}</td></tr>
    </table>
    <p>In-silico PCR: <b>{len(pcr.get('on_target', []))}</b> on-target product(s),
       <b style="color:{C_RED if pcr.get('n_off') else C_BLUE}">
       {pcr.get('n_off', 0)}</b> off-target product(s)
       (of which {pcr.get('n_perfect_off', 0)} with both ends exact match);
       search method: {_e(det.get('method', '-'))}</p>'''
        offs = pcr.get("off_target", [])[:5]
        if offs:
            rows = "".join(
                f"<tr><td>{_e(o['subject'])}</td><td>{o['size']}</td>"
                f"<td>{o['left']}/{o['right']}</td>"
                f"<td>{o['left_mm']}+{o['right_mm']}</td>"
                f"<td>{o['left_pos']}-{o['right_pos']}</td></tr>"
                for o in offs)
            spec_body += (f'<table class="mini"><tr><th>Sequence</th><th>Product bp</th>'
                          f'<th>Primer combo</th><th>Mismatches</th><th>Location</th></tr>'
                          f'{rows}</table>')
    else:
        spec_body = (f'<p class="warn">{_e(det.get("error", "Specificity search not performed"))}'
                     f'. <b>BLAST verification is required before the experiment.</b></p>')

    # --- enzyme sites (both BsaI for the seamless scheme and PaqCI for the
    #     Golden Gate scheme; each scanned on both strands) ---
    enz_body = (f'<h4 style="margin:6px 0 2px">Seamless cloning — '
                f'{_e(enz.get("enzyme", "BsaI"))} internal sites</h4>'
                + _enz_section(enz))
    if enz_gg:
        enz_body += (f'<h4 style="margin:12px 0 2px">Golden Gate — '
                     f'{_e(enz_gg.get("enzyme", "PaqCI"))} internal sites</h4>'
                     + _enz_section(enz_gg))

    # --- GC ---
    gc_body = (f'<p>Overall GC <b>{gc.get("overall", 0) * 100:.1f}%</b>; '
               f'window mean {gc.get("mean_window", 0) * 100:.1f}%; '
               f'max {gc.get("max_window", {}).get("gc", 0) * 100:.0f}% '
               f'@{gc.get("max_window", {}).get("start", 0)} nt; '
               f'min {gc.get("min_window", {}).get("gc", 0) * 100:.0f}% '
               f'@{gc.get("min_window", {}).get("start", 0)} nt</p>'
               + _gc_svg(gc)
               + "<ul>" + "".join(f"<li>{_e(i)}</li>"
                                  for i in gc.get("issues", [])) + "</ul>")

    def sec(title, level, body):
        return (f'<div class="risk"><h4>{title} {_badge(level)}</h4>{body}</div>')

    enz_level = _worst(enz.get("level", "LOW"), enz_gg.get("level", "LOW"))
    n_both = enz.get("n_total", 0) + enz_gg.get("n_total", 0)
    return (sec("Risk 1 · Primer specificity", spec.get("level", "UNKNOWN"), spec_body)
            + sec("Risk 2 · Amplicon GC content", gc.get("level", "UNKNOWN"), gc_body)
            + sec(f"Risk 3 · Internal type-IIS sites "
                  f"({n_both} total across BsaI + PaqCI, both strands)",
                  enz_level, enz_body)
            + sec("Risk 4 · Reading frame & terminal integrity",
                  frame.get("level", "UNKNOWN"),
                  "<ul>" + "".join(f"<li>{_e(i)}</li>"
                                   for i in frame.get("issues", [])) + "</ul>")
            + sec("Risk 5 · Amplification feasibility",
                  _worst(amp.get("level", "LOW"),
                         (r.get("pcr", {}) or {}).get("level", "LOW")),
                  _amp_body(amp, r)))


# ============================================================
# Main report
# ============================================================

def _pair_card(p, adapter_f: str, adapter_r: str, schemes,
              title: str = "Selected primer pair") -> str:
    """One selected/anchor pair card.

    When `schemes` (scheme_key -> order sheet) is supplied, the card lists the
    FULL ordered primer sequences for BOTH cloning schemes so the user can pick
    at order time; otherwise it falls back to the single seamless sheet.
    """
    r = p.risk or {}
    if schemes:
        blocks = []
        for k, s in schemes.items():
            sc = ADAPTER_SCHEMES[k]
            blocks.append(
                '<div class="seqrow">'
                f'<div class="scname">{_e(sc["name"])} '
                f'<span class="scenz">{_e(sc["enzyme"])}</span> '
                f'<span class="scover">sticky end {_e(s["left_overhang"].lower())} '
                f'/ {_e(s["right_overhang"].lower())}</span></div>'
                f'<div class="seq">F ({s["F_len"]} nt): '
                f'<span class="ad">{_e(s["adapter_f"])}</span>'
                f'<b>{_e(s["F_core"])}</b></div>'
                f'<div class="seq">R ({s["R_len"]} nt): '
                f'<span class="ad">{_e(s["adapter_r"])}</span>'
                f'<b>{_e(s["R_core"])}</b></div>')
        seq_blocks = "".join(blocks)
        sheet = schemes.get("seamless") or next(iter(schemes.values()))
    else:
        sheet = order_sheet(p, adapter_f, adapter_r)
        seq_blocks = (
            f'<div class="seq">GS-F ({sheet["F_len"]} nt): '
            f'<span class="ad">{_e(sheet["adapter_f"])}</span>'
            f'<b>{_e(sheet["F_core"])}</b></div>'
            f'<div class="seq">GS-R ({sheet["R_len"]} nt): '
            f'<span class="ad">{_e(sheet["adapter_r"])}</span>'
            f'<b>{_e(sheet["R_core"])}</b></div>')

    # per-scheme enzyme-site summary (BsaI for seamless, PaqCI for golden_gate)
    enz_summary = ""
    if schemes:
        for k in schemes:
            sc = ADAPTER_SCHEMES[k]
            er = (r.get("enzyme_gg") if k == "golden_gate"
                  else r.get("enzyme", {}))
            if not er:
                continue
            n = er.get("n_total", 0)
            col = LEVEL_COLOR.get(er.get("level"), C_GREY)
            msg = (f' — {_e(er.get("message", ""))}' if n else ", assembly safe")
            enz_summary += (
                f'<div class="sub2"><b>{_e(sc["name"])}:</b> '
                f'internal <b style="color:{col}">{_e(sc["enzyme"])} {n} site(s)</b>'
                f'(sense {er.get("n_plus", 0)} / antisense {er.get("n_minus", 0)})'
                f'{msg}</div>')

    return f'''<div class="card">
  <h3>{_e(title)} · {_e(p.label)} {_badge(r.get("overall", "UNKNOWN"))}</h3>
  <div class="order">
    <div><b>Order sequences (adapters attached, ready to order — two schemes to choose)</b></div>
    <div style="margin-top:6px">{seq_blocks}</div>
    <div style="margin-top:8px;font-size:12px;color:#666">
      Grey = fixed adapter; bold = template-specific region; PCR product {sheet['product_len']} bp,
      insert after digestion {sheet['insert_len']} bp.
      Tm / GC below refer to the <b>core primer (adapter excluded)</b>; the full ordered primer's
      Tm is higher and is not used for annealing planning.
    </div>
  </div>
  {_construct_diagram(sheet, p)}
  {enz_summary}
  <div class="kv" style="margin-top:14px">
    <b>Forward primer</b><span class="mono">{_e(p.f_seq)} | {p.f_len} nt |
      Tm {p.f_qc['tm']:.1f}℃ | GC {p.f_qc['gc'] * 100:.1f}% |
      hairpin ΔG {p.f_qc['hairpin_dg']:.1f} | self-dimer ΔG {p.f_qc['self_dimer_dg']:.1f}</span>
    <b>Reverse primer</b><span class="mono">{_e(p.r_seq)} | {p.r_len} nt |
      Tm {p.r_qc['tm']:.1f}℃ | GC {p.r_qc['gc'] * 100:.1f}% |
      hairpin ΔG {p.f_qc['hairpin_dg']:.1f} | self-dimer ΔG {p.r_qc['self_dimer_dg']:.1f}</span>
    <b>Pair</b><span class="mono">ΔTm {p.tm_diff:.2f}℃ |
      hetero-dimer ΔG {p.hetero_dimer_dg:.1f} kcal/mol |
      amplicon {p.amplicon_len} bp</span>
    <b>Anchoring</b><span>forward offset {p.offset:+d} nt (multiple of 3
      {'✓' if p.offset % 3 == 0 else '✗'}); reverse
      {'strict stop-codon end' if p.r_ext == 0 else '+%d nt past stop' % p.r_ext}</span>
  </div>
  {_risk_block(p)}
</div>'''


CSS = """
body{font-family:"Segoe UI",system-ui,sans-serif;
 margin:0;padding:28px 34px;background:#ffffff;color:#232323;line-height:1.6;font-size:14px}
h1{font-size:22px;margin:0 0 4px;color:#1a1a1a;font-weight:650}
h2{font-size:17px;margin:30px 0 10px;padding-bottom:6px;
 border-bottom:2px solid #2166AC;color:#2166AC}
h3{font-size:15px;margin:20px 0 8px;color:#1a1a1a}
h4{font-size:13.5px;margin:0 0 8px;color:#333}
.sub{color:#777;font-size:12.5px;margin-bottom:18px}
.sub2{font-size:12.5px;margin:8px 0;color:#333}
.sub2 b{color:#1a1a1a}
table{border-collapse:collapse;width:100%;margin:10px 0;font-size:12.5px}
th,td{border:1px solid #e2e2e2;padding:6px 9px;text-align:left;vertical-align:top}
th{background:#f4f7fa;font-weight:600;color:#2166AC;white-space:nowrap}
tr.anchor td{background:#eef5fb}
tr.sel td{background:#fff6e8}
code,.mono{font-family:"Cascadia Mono",Consolas,monospace;font-size:12px}
.badge{display:inline-block;color:#fff;border-radius:9px;padding:1px 9px;
 font-size:11px;font-weight:600;margin-left:6px;vertical-align:middle}
.card{border:1px solid #e2e2e2;border-radius:7px;padding:16px 18px;margin:14px 0;
 background:#fff}
.risk{border-left:3px solid #d8d8d8;padding:2px 0 2px 14px;margin:16px 0}
.mini{font-size:12px}
.mini th{background:#fafafa;color:#555}
.warn{color:#B2182B;font-weight:600}
.kv{display:grid;grid-template-columns:130px 1fr;gap:3px 12px;font-size:13px}
.kv b{color:#555;font-weight:500}
.construct{margin:12px 0}
.crow{font-family:"Cascadia Mono",Consolas,monospace;font-size:12.5px;
 word-break:break-all;line-height:2.1}
.seg{padding:3px 2px;border-radius:2px}
.seg.tail{background:#eeeeee}
.seg.site{background:#f6c6c9;font-weight:700}
.seg.spacer{background:#fde7c3}
.seg.ovh{background:#cfe3f3;font-weight:700}
.seg.core{background:#dff0d8}
.seg.dots{color:#999;padding:0 6px}
.legend{margin-top:8px;font-size:11.5px;color:#666}
.legend span{margin-right:14px}
.legend i{display:inline-block;width:11px;height:11px;border-radius:2px;
 margin-right:4px;vertical-align:-1px}
.note{margin-top:10px;font-size:12.5px;color:#444;background:#f8f9fa;
 padding:8px 11px;border-radius:5px}
.order{background:#f8f9fa;border:1px dashed #bbb;border-radius:6px;
 padding:12px 14px;margin:12px 0}
.order .seq{font-family:"Cascadia Mono",Consolas,monospace;font-size:13px;
 word-break:break-all}
.seqrow{margin:8px 0;padding:8px 10px;border:1px solid #e2e2e2;
 border-radius:5px;background:#fff}
.seqrow .scname{font-weight:600;font-size:12.5px;color:#1a1a1a}
.seqrow .scenz{font-size:11px;color:#2166AC;font-weight:600}
.seqrow .scover{font-size:11px;color:#666}
.seqrow .seq{font-family:"Cascadia Mono",Consolas,monospace;font-size:12.5px;
 margin-top:3px;word-break:break-all}
.seqrow .ad{color:#999}
.gcplot{width:100%;max-width:900px;height:auto;margin:6px 0}
ul{margin:6px 0 6px 20px;padding:0}
li{margin:2px 0}
footer{margin-top:36px;padding-top:12px;border-top:1px solid #e5e5e5;
 color:#999;font-size:11.5px}
"""


def generate_report(*, tx, isoforms: List[Dict], pairs: List,
                    selected: Optional[List] = None, params: Dict,
                    out_path: str, gene_query: str = "",
                    adapter_f: str = "", adapter_r: str = "",
                    scheme_sheets: Optional[Dict] = None) -> str:
    selected = selected or []
    sel_ids = {id(p) for p in selected}

    iso_rows = "".join(
        f'<tr class="{"anchor" if r["transcript"] == tx.tx_id else ""}">'
        f'<td>{"► " if r["transcript"] == tx.tx_id else ""}{_e(r["transcript"])}</td>'
        f'<td>{r["mRNA_len"]}</td><td>{r["CDS_len"]}</td>'
        f'<td>{r["protein_aa"]}</td><td>{r["utr5"]}</td><td>{r["utr3"]}</td>'
        f'<td>{r["exons"]}</td><td>{_e(r["issues"])}</td></tr>'
        for r in isoforms)

    cand_rows = ""
    for i, p in enumerate(pairs, 1):
        cls = "sel" if id(p) in sel_ids else ("anchor" if p.is_anchor else "")
        r = p.risk or {}
        enz = r.get("enzyme", {}) or {}
        enz_gg = r.get("enzyme_gg", {}) or {}
        cand_rows += (
            f'<tr class="{cls}"><td>{i}{" ★" if p.is_anchor else ""}</td>'
            f'<td>{p.tier}</td><td>{p.offset:+d}</td>'
            f'<td class="mono">{_e(p.f_seq)}</td>'
            f'<td>{p.f_len}</td><td>{p.f_qc["tm"]:.1f}</td>'
            f'<td>{p.f_qc["gc"] * 100:.0f}</td>'
            f'<td class="mono">{_e(p.r_seq)}</td>'
            f'<td>{p.r_len}</td><td>{p.r_qc["tm"]:.1f}</td>'
            f'<td>{p.r_qc["gc"] * 100:.0f}</td>'
            f'<td>{p.tm_diff:.1f}</td><td>{p.amplicon_len}</td>'
            f'<td>{enz.get("n_total", "-")}</td>'
            f'<td>{enz_gg.get("n_total", "-")}</td>'
            f'<td>{_badge(r.get("overall", "UNKNOWN")) if r else "-"}</td>'
            f'<td>{_e("; ".join(p.flags) or "OK")}</td></tr>')

    sel_html = ""
    # Stage 1 (no selection): surface the full risk work-up for the best-ranked
    # candidate (always an ATG..stop anchor) so the report itself flags
    # specificity, GC, internal type-IIS (both BsaI and PaqCI strands), frame and
    # amplification risks without duplicating a card per candidate.
    risk_pairs = selected if selected else (pairs[:1] if pairs else [])
    for p in risk_pairs:
        pid = pairs.index(p) + 1
        sheets = (scheme_sheets or {}).get(pid)
        sel_html += _pair_card(
            p, adapter_f or params.get("adapter_f", ""),
            adapter_r or params.get("adapter_r", ""), sheets,
            title=('Selected primer pair' if selected else 'Standard anchor pair (risk detail)'))

    p_ = params
    param_rows = "".join(f"<tr><td>{_e(k)}</td><td>{_e(v)}</td></tr>"
                         for k, v in p_.items() if not k.startswith("_"))

    doc = f'''<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>GSprimer report · {_e(gene_query or tx.tx_id)}</title>
<style>{CSS}</style></head><body>
<h1>GSprimer · CDS-amplification primer design for GS / Golden Gate vectors</h1>
<div class="sub">Target: <b>{_e(gene_query or tx.tx_id)}</b> ·
  template transcript <b>{_e(tx.tx_id)}</b> ({_e(tx.source)} mode) ·
  generated {time.strftime('%Y-%m-%d %H:%M:%S')}</div>

<h2>1. Template transcript</h2>
<div class="kv">
  <b>Transcript</b><span>{_e(tx.tx_id)} (gene {_e(tx.gene_id or '-')},
    {_e(tx.chrom or '-')}{' ' + tx.strand if tx.chrom else ''})</span>
  <b>mRNA length</b><span>{len(tx.mrna)} nt (5'UTR {tx.utr5_len} /
    CDS {tx.cds_len} / 3'UTR {tx.utr3_len})</span>
  <b>Protein</b><span>{tx.protein_len} aa</span>
  <b>Structure check</b><span>{_e('; '.join(tx.validate()) or 'pass: ATG start, stop-codon end, length multiple of 3, no internal stop')}</span>
  <b>Notes</b><span>{_e('; '.join(tx.notes) or '-')}</span>
</div>
<table><tr><th>Transcript</th><th>mRNA</th><th>CDS</th><th>Protein(aa)</th>
<th>5'UTR</th><th>3'UTR</th><th>Exons</th><th>Structure issues</th></tr>
{iso_rows}</table>
<p style="font-size:12.5px;color:#666">► marks the chosen transcript (default rule: longest CDS;
tie-break by longer mRNA).</p>

<h2>2. Candidate primer pairs</h2>
<p style="font-size:12.5px;color:#666">★ = standard anchor pair with strict ATG start / stop-codon
end (mandated by the design spec to be present in candidates). Tier A = all SpacerFinder-derived
filters pass; B = boundaries acceptable; C = needs trade-off. offset is the forward primer 5' shift
relative to ATG, always a multiple of 3. <b>BsaI</b> / <b>PaqCI</b> = number of internal
type-IIS sites found inside the PCR product for the seamless (BsaI) and Golden Gate (PaqCI) schemes
respectively — 0 in both is required for clean assembly.</p>
<table><tr><th>#</th><th>Tier</th><th>offset</th>
<th>Forward 5'→3'</th><th>nt</th><th>Tm (core)</th><th>GC%</th>
<th>Reverse 5'→3'</th><th>nt</th><th>Tm (core)</th><th>GC%</th>
<th>ΔTm</th><th>Product bp</th><th>BsaI</th><th>PaqCI</th><th>Risk</th><th>Flags</th></tr>
{cand_rows}</table>
<p style="font-size:12.5px;color:#666"><b>Tm (core)</b>: melting temperature of the <b>primer core only</b> (template-specific region), computed <b>without</b> the 5' adapter / type-IIS module (≈17 nt). Ordered primers include the adapter, but its length and Tm are <b>excluded</b> from this column — the full ordered primer's Tm is ~6–9 °C higher and is <b>not</b> used for annealing-temperature planning. ΔTm is computed between the two core Tms.</p>

{'<h2>3. Standard anchor pair &amp; candidate risk assessment</h2>' + sel_html if sel_html else
 '<h2>3. Selected primer pairs</h2><p class="warn">No primer pair selected yet '
 '— please select in the interactive step and regenerate the report.</p>'}

<h2>{'4' if sel_html else '4'}. Design parameters</h2>
<table><tr><th>Parameter</th><th>Value</th></tr>{param_rows}</table>

<footer>GSprimer v1.2.0 · thermodynamics &amp; filters inherited from SpacerFinder v2.4.5
(amplicon-size limit excluded) · Tm shown is the core-primer Tm (adapter excluded) · palette Fv1.4</footer>
</body></html>'''

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(doc)
    return out_path
