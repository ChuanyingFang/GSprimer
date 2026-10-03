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
        description="CDS-amplification primer design for GS / Golden Gate (BsaI/PaqCI) vectors "
                    "(primer thermodynamics inherited from SpacerFinder)",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    g = p.add_argument_group("Target")
    g.add_argument("--gene", "-g", default="", help="gene / transcript ID")
    g.add_argument("--transcript", "-t", default="",
                   help="specific transcript (default: longest)")
    g.add_argument("--longest-by", choices=["cds", "mrna"], default="cds",
                   help="basis for longest-transcript selection (default cds)")

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
    b.add_argument("--blast-db-cdna", default="",
                   help="transcriptome/cDNA BLAST DB prefix; preferred over "
                        "--blast-db for in-silico PCR when the gene is spliced "
                        "(a genomic DB cannot yield the cDNA-sized product)")
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

    sequence = args.sequence or _read_seq_file(args.sequence_file)
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
        blast_db_cdna=args.blast_db_cdna,
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
        print(f"\n[{i}] {p.label}   product {p.amplicon_len} bp   "
              f"risk {r.get('overall', 'NA')}")
        print(f"  F ({p.f_len} nt): {p.f_seq}")
        print(f"  R ({p.r_len} nt): {p.r_seq}")
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
