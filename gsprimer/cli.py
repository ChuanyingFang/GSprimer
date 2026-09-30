#!/usr/bin/env python3
"""
GSprimer — command line interface.

Examples
--------
  # design from genome + GFF3 (rice MSU7.0)
  python -m gsprimer.cli --gene LOC_Os01g01010 \
      --genome MSU_dna.fa --gff3 MSU.gff3 --blast-db MSU_dna_db \
      --outdir out

  # design from a CDS FASTA, then pick the pair interactively and confirm
  python -m gsprimer.cli --gene LOC_Os01g01010 --cds MSU_cds.fa \
      --outdir out --interactive

  # promoter mode (v2.0): amplify the region upstream of the TSS
  #   upstream_2k.fa = ~2 kb upstream region (sequence end = TSS)
  #   utr5.fa         = transcript 5'-UTR fallback (only used when no ideal R in 2k)
  python -m gsprimer.cli --promoter upstream_2k.fa --utr5 utr5.fa \
      --outdir out --product-min 1200 --product-max 1700

  # finalize later, reusing the cached design state (lists the chosen primers
  # and asks for confirmation BEFORE any adapter is attached)
  python -m gsprimer.cli --finalize --select 1,3 --outdir out

  # scripted pipelines: skip the confirmation prompt with --yes
  python -m gsprimer.cli --gene X --sequence-file cds.fa --outdir out \
      --select 1 --yes
"""

import argparse
import os
import sys

from .adapters import ADAPTER_F, ADAPTER_R, ADAPTER_SCHEMES
from .config import get, set_config
from .pipeline import (candidate_table, load_state, pairs_from_state,
                       run_design, run_finalize)
from .template import read_fasta


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="gsprimer",
        description="CDS / promoter amplification primer design for GS / "
                    "Golden Gate (BsaI/PaqCI) vectors "
                    "(primer thermodynamics inherited from SpacerFinder)",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    g = p.add_argument_group("Target")
    g.add_argument("--gene", "-g", default="", help="gene / transcript ID")
    g.add_argument("--transcript", "-t", default="",
                   help="specific transcript (default: longest)")
    g.add_argument("--longest-by", choices=["cds", "mrna"], default="cds",
                   help="basis for longest-transcript selection (default cds)")

    # ---- promoter / functional-element amplification mode (v2.0) ----
    pr = p.add_argument_group("Promoter mode (v2.0, amplify the region upstream of the TSS)")
    pr.add_argument("--promoter", default="",
                    help="upstream promoter region sequence (5'->3' sense; "
                         "sequence end = TSS), or a FASTA / text file holding it. "
                         "Mutually exclusive with --gene/--sequence")
    pr.add_argument("--utr5", default="",
                    help="optional: transcript 5'-UTR sequence (used as fallback when "
                         "no ideal R is found within the 2 kb upstream region; "
                         "extend at most 200 nt into the 5'-UTR)")
    pr.add_argument("--promoter-file", default="",
                    help="upstream promoter region file (same as --promoter, pick one)")
    pr.add_argument("--utr5-file", default="",
                    help="5'-UTR sequence file (same as --utr5, pick one)")
    pr.add_argument("--product-min", type=int, default=1200,
                    help="minimum PCR product length (default 1200)")
    pr.add_argument("--product-max", type=int, default=1700,
                    help="maximum PCR product length (default 1700)")
    pr.add_argument("--r-up-max", type=int, default=100,
                    help="max nt the reverse primer 5' end may extend upstream of the "
                         "TSS (within the promoter region; default 100)")
    pr.add_argument("--r-utr-max", type=int, default=200,
                    help="max nt the reverse primer 5' end may extend into the 5'-UTR "
                         "(fallback mode only; default 200)")

    s = p.add_argument_group("Template source (pick one)")
    s.add_argument("--genome", default="", help="genome FASTA")
    s.add_argument("--gff3", default="", help="GFF3 annotation")
    s.add_argument("--cdna", default="", help="cDNA/mRNA multi-FASTA")
    s.add_argument("--cds", default="", help="CDS multi-FASTA")
    s.add_argument("--sequence", default="", help="paste raw sequence string")
    s.add_argument("--sequence-file", default="",
                   help="single-sequence FASTA / plain-text file")

    d = p.add_argument_group("Design parameters")
    d.add_argument("--max-shift", type=int, default=12,
                   help="max forward-primer shift in nt (auto-snapped to a "
                        "multiple of 3; default 12)")
    d.add_argument("--r-ext", type=int, default=0,
                   help="max reverse-primer extension past the stop codon in nt "
                        "(default 0 = strict)")
    d.add_argument("--top", type=int, default=10, help="number of candidate pairs")
    d.add_argument("--no-downstream", action="store_true",
                   help="forbid forward primer shifting into the CDS (avoids "
                        "N-terminal truncation)")

    b = p.add_argument_group("Specificity")
    b.add_argument("--blast-db", default="", help="nucleotide BLAST DB prefix")
    b.add_argument("--blastn", default="", help="path to the blastn executable")
    b.add_argument("--fallback-fasta", default="",
                   help="FASTA used for seed-and-extend when BLAST is absent "
                        "(e.g. full cDNA set)")
    b.add_argument("--no-blast", action="store_true", help="skip specificity search")

    v = p.add_argument_group("Vector adapters")
    v.add_argument("--adapter-f", default=ADAPTER_F, help="forward adapter")
    v.add_argument("--adapter-r", default=ADAPTER_R, help="reverse adapter")
    v.add_argument("--enzyme", default="BsaI", help="assembly enzyme name")
    v.add_argument("--site", default="GGTCTC", help="enzyme recognition sequence")
    v.add_argument("--scheme", default="all",
                   choices=["all", "seamless", "golden_gate"],
                   help="terminal print filter: which adapter scheme to show; "
                        "default 'all' = both the seamless (BsaI) and Golden "
                        "Gate (PaqCI) sets, for the user to choose at order time")

    o = p.add_argument_group("Output & workflow")
    o.add_argument("--outdir", "-o", default="gsprimer_out", help="output directory")
    o.add_argument("--prefix", default="", help="output file prefix")
    o.add_argument("--select", default="",
                   help="selected candidate number(s), e.g. 1 or 1,3")
    o.add_argument("--interactive", "-i", action="store_true",
                   help="enter interactive selection after design")
    o.add_argument("--finalize", action="store_true",
                   help="skip design and finalize directly from cached state")
    o.add_argument("--yes", action="store_true",
                   help="skip the interactive primer-sequence confirmation "
                        "prompt (attach adapters without asking; for scripts)")

    n = p.add_argument_group("Structure hints (optional, affects template advice)")
    n.add_argument("--no-intron", action="store_true",
                   help="declare the gene is intronless (CDS/cDNA/direct mode): "
                        "report will recommend gDNA as template, skipping RT")
    n.add_argument("--has-introns", action="store_true",
                   help="declare the gene has introns: force cDNA template even "
                        "outside GFF3 mode")
    return p


def _read_seq_file(path: str) -> str:
    if not path:
        return ""
    with open(path, encoding="utf-8", errors="replace") as f:
        txt = f.read()
    if txt.lstrip().startswith(">"):
        d = read_fasta(path)
        return next(iter(d.values())) if d else ""
    return txt


def _resolve_seq(arg: str) -> str:
    """Resolve a CLI argument that may be a raw sequence, an inline FASTA string,
    or a path to a FASTA / plain-text file into the bare sequence (uppercased,
    whitespace stripped)."""
    if not arg:
        return ""
    if arg.lstrip().startswith(">"):
        # inline FASTA pasted on the command line
        return "".join(l.strip() for l in arg.splitlines() if not l.startswith(">"))
    if os.path.exists(arg):
        return _read_seq_file(arg)
    # otherwise treat as a raw sequence string
    return "".join(arg.split()).upper()


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.finalize:
        st = load_state(args.outdir)
        if not st:
            print(f"Error: design state not found in {args.outdir}; "
                  f"run the design stage first", file=sys.stderr)
            return 2
        pairs = pairs_from_state(st)
        if args.select:
            sel = [int(x) for x in args.select.replace(" ", "").split(",") if x]
        else:
            # No explicit selection: present the list and ask the user to pick.
            print(candidate_table(pairs))
            from .interactive import ask_selection
            sel = ask_selection(pairs)
        bad = [i for i in sel if i < 1 or i > len(pairs)]
        if bad:
            print(f"Error: number {bad} out of range (1-{len(pairs)})",
                  file=sys.stderr)
            return 2
        if not _confirm_primers(pairs, sel, args.yes):
            print("Aborted: GS adapters were NOT attached; no order sheet "
                  "written. Re-run with your confirmed selection when ready.")
            return 0
        res = run_finalize(sel, outdir=args.outdir)
        if not res["ok"]:
            print(f"Error: {res['error']}", file=sys.stderr)
            return 1
        _print_final(res, args.scheme)
        return 0

    sequence = _resolve_seq(args.sequence) or _read_seq_file(args.sequence_file)
    promoter_seq = _resolve_seq(args.promoter) or _read_seq_file(args.promoter_file)
    utr5_seq = _resolve_seq(args.utr5) or _read_seq_file(args.utr5_file)

    if promoter_seq:
        if args.gene or sequence:
            print("Error: promoter mode (--promoter) is mutually exclusive "
                  "with --gene/--sequence", file=sys.stderr)
            return 2
        if not (1200 <= args.product_max):
            print("Error: --product-max must be >= 1200", file=sys.stderr)
            return 2
        res = run_design(
            args.gene or "promoter",
            promoter_seq=promoter_seq, utr5_seq=utr5_seq, outdir=args.outdir,
            prefix=args.prefix, top_n=args.top,
            blastn=args.blastn, blast_db=args.blast_db,
            fallback_fasta=args.fallback_fasta, do_blast=not args.no_blast,
            adapter_f=args.adapter_f, adapter_r=args.adapter_r,
            enzyme=args.enzyme, site=args.site,
            product_min=args.product_min, product_max=args.product_max,
            r_up_max=args.r_up_max, r_utr_max=args.r_utr_max)
    else:
        if not (args.gene or sequence):
            print("Error: --gene or --sequence/--sequence-file is required",
                  file=sys.stderr)
            return 2

        if args.no_intron and args.has_introns:
            print("Error: --no-intron and --has-introns are mutually exclusive",
                  file=sys.stderr)
            return 2
        no_intron = True if args.no_intron else (False if args.has_introns else None)

        res = run_design(
            args.gene or "user_seq",
            genome_fa=args.genome or get("genome"),
            gff3=args.gff3 or get("gff3"),
            cdna_fa=args.cdna or get("cdna"),
            cds_fa=args.cds or get("cds"),
            sequence=sequence,
            transcript=args.transcript, outdir=args.outdir, prefix=args.prefix,
            max_shift=args.max_shift, r_ext_max=args.r_ext, top_n=args.top,
            allow_downstream=not args.no_downstream,
            blastn=args.blastn, blast_db=args.blast_db,
            fallback_fasta=args.fallback_fasta, do_blast=not args.no_blast,
            adapter_f=args.adapter_f, adapter_r=args.adapter_r,
            enzyme=args.enzyme, site=args.site, no_intron=no_intron,
            longest_by=args.longest_by)

    if not res["ok"]:
        print(f"Error: {res['error']}", file=sys.stderr)
        return 1

    print()
    print(candidate_table(res["pairs"]))
    print()
    print(f"Candidate table: {res['tsv']}")
    print(f"Report         : {res['html']}")

    sel = []
    if args.select:
        sel = [int(x) for x in args.select.replace(" ", "").split(",") if x]
    elif args.interactive:
        from .interactive import ask_selection
        sel = ask_selection(res["pairs"])

    if sel:
        if not _confirm_primers(res["pairs"], sel, args.yes):
            print("Aborted: GS adapters were NOT attached. The design report "
                  "above is still valid; finalize later once the primer "
                  "sequences are confirmed.")
            return 0
        fin = run_finalize(sel, outdir=args.outdir)
        if not fin["ok"]:
            print(f"Error: {fin['error']}", file=sys.stderr)
            return 1
        _print_final(fin, args.scheme)
    else:
        print("\n>>> Next step: after the user picks a primer-pair number, run")
        print(f"    python -m gsprimer.cli --finalize --select <N> "
              f"--outdir {args.outdir}")
    return 0


def _confirm_primers(pairs, selection, yes: bool = False) -> bool:
    """Show the selected primer CORE sequences and require explicit confirmation
    before GS adapters are attached (the user must sign off on the actual
    sequences, not just a candidate number)."""
    print("\n" + "=" * 72)
    print("PRIMER CONFIRMATION — review the sequences below")
    print("=" * 72)
    print("Adapters are NOT attached yet. Confirm these exact primer sequences:")
    for i in selection:
        p = pairs[i - 1]
        r = p.risk or {}
        enz = (r.get("enzyme", {}) or {})
        enz_gg = (r.get("enzyme_gg", {}) or {})
        is_promo = getattr(p, "mode", "cds") == "promoter"
        print(f"\n[{i}] {p.label}   product {p.amplicon_len} bp   "
              f"risk {r.get('overall', 'NA')}")
        print(f"  F ({p.f_len} nt): {p.f_seq}")
        print(f"  R ({p.r_len} nt): {p.r_seq}")
        if is_promo:
            if p.offset == 0:
                r_txt = "precisely anchored at TSS (5'-UTR start)"
            elif p.offset < 0:
                r_txt = ("extends %d nt upstream of the 5'-UTR "
                         "(within the promoter region)" % -p.offset)
            else:
                r_txt = ("extends %d nt into the 5'-UTR "
                         "(fallback to transcript UTR)" % p.offset)
            print(f"  anchor  R 5' end: {r_txt}")
        else:
            print(f"  anchor  F: {'ATG start' if p.offset == 0 else '%+d nt' % p.offset}"
                  f" | R: {'stop-codon end' if p.r_ext == 0 else '+%d nt past stop' % p.r_ext}")
        print(f"  internal BsaI  sites: {enz.get('n_total', '?')} "
              f"(+{enz.get('n_plus', '?')} / -{enz.get('n_minus', '?')})")
        print(f"  internal PaqCI sites: {enz_gg.get('n_total', '?')} "
              f"(+{enz_gg.get('n_plus', '?')} / -{enz_gg.get('n_minus', '?')})")
    if yes:
        print("\n  [--yes] confirmation skipped (adapters will be attached)")
        return True
    try:
        ans = input("\nAttach GS adapters to the primer pair(s) above and write "
                    "the order sheet? [y/N] ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        print("\n[cancelled]")
        return False
    return ans in ("y", "yes")


def _print_final(res: dict, scheme: str = "all") -> None:
    print("\n" + "=" * 72)
    print("Selected primers (GS adapters attached, ready to order)")
    print("=" * 72)
    s_items = list(res["scheme_sheets"].items())
    for j, (pid, ss) in enumerate(s_items):
        p = res["selected"][j]
        print(f"\n[{pid}] {p.label}   product {p.amplicon_len} bp   "
              f"risk {(p.risk or {}).get('overall', 'NA')}")
        # filter by the user-chosen scheme (all = both sets)
        want = [k for k in ADAPTER_SCHEMES if scheme in ("all", k)]
        for k in want:
            s = ss[k]
            sc = ADAPTER_SCHEMES[k]
            print(f"  ── {sc['name']} / {sc['enzyme']} ──")
            print(f"    F ({s['F_len']} nt): {s['F_full']}")
            print(f"    R ({s['R_len']} nt): {s['R_full']}")
        enz = (p.risk or {}).get("enzyme", {})
        enz_gg = (p.risk or {}).get("enzyme_gg", {})
        print(f"  Internal BsaI sites : {enz.get('n_total', '?')} "
              f"(+{enz.get('n_plus', '?')} / -{enz.get('n_minus', '?')})")
        print(f"  Internal PaqCI sites: {enz_gg.get('n_total', '?')} "
              f"(+{enz_gg.get('n_plus', '?')} / -{enz_gg.get('n_minus', '?')})")
    print(f"\nFinal report: {res['html']}")
    print(f"Order sheet : {res['order']}")


if __name__ == "__main__":
    sys.exit(main())
