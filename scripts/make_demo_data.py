#!/usr/bin/env python3
"""
Build a self-contained demo dataset for GSprimer smoke tests.

Produces, in the target directory:
  demo_cds.fa      two CDS entries (one clean, one carrying internal BsaI sites)
  demo_cdna.fa     the same genes with 5'/3' UTRs
  demo_genome.fa   a 2-"chromosome" mini genome
  demo.gff3        gene / mRNA / exon / CDS records, incl. a 2-isoform gene
"""

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gsprimer.risk import CODON_TABLE  # noqa: E402
from gsprimer.thermo import evaluate_primer, revcomp  # noqa: E402

SENSE_CODONS = [c for c, a in CODON_TABLE.items() if a != "*" and c != "ATG"]
STOPS = ("TAA", "TAG", "TGA")


def make_cds(n_codons: int, seed: int, gc_target: float = 0.50,
             inserts=()) -> str:
    """Random in-frame CDS with an optional list of (codon_index, motif)."""
    rng = random.Random(seed)
    body = []
    for _ in range(n_codons):
        pool = sorted(SENSE_CODONS,
                      key=lambda c: abs((c.count("G") + c.count("C")) / 3
                                        - gc_target) + rng.random() * 0.4)
        body.append(pool[0])
    seq = "ATG" + "".join(body)
    for ci, motif in inserts:
        i = 3 * ci
        seq = seq[:i] + motif + seq[i + len(motif):]
    seq = seq + "TGA"
    # guarantee no premature stop
    fixed = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    for k in range(1, len(fixed) - 1):
        while fixed[k] in STOPS:
            fixed[k] = random.Random(seed + k).choice(SENSE_CODONS)
    return "".join(fixed)


def gc_block(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("GCGCGCGCAT") for _ in range(n))


def _no_site(seq: str) -> bool:
    return "GGTCTC" not in seq and "GAGACC" not in seq


def _find_anchor_head(rng: random.Random, n_cand: int = 400):
    """Find a 24-nt segment after ATG so that ATG+head[:k] passes every
    inherited SpacerFinder filter with Tm in [55,60] C for some k."""
    heads = []
    while len(heads) < n_cand:
        seg = "".join(rng.choice("ACGT") for _ in range(24))
        if not _no_site(seg):
            continue
        for k in range(16, 22):              # total primer len 19..24
            pr = evaluate_primer("ATG" + seg[:k])
            if pr["pass"] and 55.0 <= pr["tm"] <= 60.0:
                heads.append((seg, pr["tm"]))
                break
    return heads


def _find_anchor_tail(rng: random.Random, n_cand: int = 400):
    """Find a 27-nt segment before the stop codon so that the reverse primer
    (revcomp of tail[:k]+TGA) passes every inherited filter with Tm in [55,60]."""
    tails = []
    while len(tails) < n_cand:
        seg = "".join(rng.choice("ACGT") for _ in range(27))
        if not _no_site(seg):
            continue
        for k in range(19, 25):              # foot length 19..24 + 3 nt stop
            pr = evaluate_primer(revcomp(seg[:k] + "TGA"))
            if pr["pass"] and 55.0 <= pr["tm"] <= 60.0:
                tails.append((seg, pr["tm"]))
                break
    return tails


def build_clean_cds(n_codons: int, seed: int) -> str:
    """Assemble a CDS whose ATG-head and stop-tail each yield a Tier-A anchor
    primer (dTm < 2 C), with an in-frame random body in between."""
    rng = random.Random(seed)
    heads = _find_anchor_head(rng)
    tails = _find_anchor_tail(rng)
    chosen = None
    for h in heads:
        for t in tails:
            if abs(h[1] - t[1]) < 2.0:
                chosen = (h[0], t[0])
                break
        if chosen:
            break
    if chosen is None:                  # extremely unlikely, relax dTm
        chosen = (heads[0][0], tails[0][0])
    head_seg, tail_seg = chosen
    head_codons = len(head_seg) // 3
    tail_codons = len(tail_seg) // 3
    mid = make_cds(n_codons - 1 - head_codons - tail_codons, seed=seed + 1,
                   gc_target=0.50)
    # mid already starts with ATG; drop it and the stop
    mid_body = mid[3:-3]
    seq = "ATG" + head_seg + mid_body + tail_seg + "TGA"
    # fix any premature stop, then scrub BsaI sites
    fixed = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    for k in range(1, len(fixed) - 1):
        while fixed[k] in STOPS:
            fixed[k] = rng.choice(SENSE_CODONS)
    seq = "".join(fixed)
    while "GGTCTC" in seq or "GAGACC" in seq:
        seq = seq.replace("GGTCTC", "GGTCTG").replace("GAGACC", "GAGACG")
    return seq


def utr(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_fasta(path, records):
    with open(path, "w", encoding="utf-8") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


def main(outdir="demo_data"):
    os.makedirs(outdir, exist_ok=True)

    # --- GENE1: clean CDS, no internal BsaI, with a Tier-A anchor pair ---
    cds1 = build_clean_cds(298, seed=11)

    # --- GENE2: CDS carrying one BsaI site per strand + a GC-rich island ---
    cds2 = make_cds(240, seed=23, gc_target=0.52,
                    inserts=[(40, "GGTCTC"),       # sense-strand BsaI
                             (120, "GAGACC")])     # antisense-strand BsaI
    cds2 = cds2[:3 * 180] + gc_block(90, 5) + cds2[3 * 180 + 90:]
    cds2 = cds2[:len(cds2) - len(cds2) % 3]
    if cds2[-3:] not in STOPS:
        cds2 = cds2[:-3] + "TGA"

    write_fasta(os.path.join(outdir, "demo_cds.fa"),
                [("DEMO_G1.1", cds1), ("DEMO_G2.1", cds2)])

    # --- cDNA (with UTRs) ---
    u5a, u3a = utr(120, 1), utr(210, 2)
    u5b, u3b = utr(75, 3), utr(160, 4)
    cdna1 = u5a + cds1 + u3a
    cdna2 = u5b + cds2 + u3b
    # a shorter second isoform of GENE1 (tests longest-transcript selection)
    cds1_short = cds1[:3 * 150] + "TGA"
    cdna1b = u5a + cds1_short + utr(90, 6)
    write_fasta(os.path.join(outdir, "demo_cdna.fa"),
                [("DEMO_G1.1", cdna1), ("DEMO_G1.2", cdna1b),
                 ("DEMO_G2.1", cdna2)])

    # --- mini genome: gene1 has 2 exons, gene2 is single-exon ---
    intron = ("GT" + utr(196, 7) + "AG")
    pad = utr(500, 8)
    # gene1 layout: pad | ex1(u5a+cds1[:450]) | intron | ex2(cds1[450:]+u3a) | pad
    ex1 = u5a + cds1[:450]
    ex2 = cds1[450:] + u3a
    chr1 = pad + ex1 + intron + ex2 + pad
    g1_start = len(pad) + 1
    ex1_s, ex1_e = g1_start, g1_start + len(ex1) - 1
    ex2_s = ex1_e + len(intron) + 1
    ex2_e = ex2_s + len(ex2) - 1
    cds1_s = g1_start + len(u5a)
    cds1_e_ex1 = ex1_e
    cds1_s_ex2 = ex2_s
    cds1_e = ex2_s + (len(cds1) - 450) - 1

    # gene2 on the minus strand, single exon
    pad2 = utr(400, 9)
    g2_body = revcomp(cdna2)
    chr2 = pad2 + g2_body + pad2
    g2_start = len(pad2) + 1
    g2_end = g2_start + len(g2_body) - 1
    # minus strand: CDS start (ATG) is at the RIGHT end
    c2_hi = g2_end - len(u5b)
    c2_lo = c2_hi - len(cds2) + 1

    write_fasta(os.path.join(outdir, "demo_genome.fa"),
                [("chrD1", chr1), ("chrD2", chr2)])

    gff = ["##gff-version 3"]

    def row(chrom, ftype, s, e, strand, attrs, phase="."):
        gff.append(f"{chrom}\tdemo\t{ftype}\t{s}\t{e}\t.\t{strand}\t{phase}\t{attrs}")

    row("chrD1", "gene", g1_start, ex2_e, "+", "ID=DEMO_G1;Name=DEMO_G1")
    row("chrD1", "mRNA", g1_start, ex2_e, "+", "ID=DEMO_G1.1;Parent=DEMO_G1")
    row("chrD1", "exon", ex1_s, ex1_e, "+", "ID=DEMO_G1.1:e1;Parent=DEMO_G1.1")
    row("chrD1", "exon", ex2_s, ex2_e, "+", "ID=DEMO_G1.1:e2;Parent=DEMO_G1.1")
    row("chrD1", "CDS", cds1_s, cds1_e_ex1, "+",
        "ID=DEMO_G1.1:c1;Parent=DEMO_G1.1", "0")
    row("chrD1", "CDS", cds1_s_ex2, cds1_e, "+",
        "ID=DEMO_G1.1:c2;Parent=DEMO_G1.1", "0")
    # short isoform: exon1 only
    row("chrD1", "mRNA", g1_start, ex1_e, "+", "ID=DEMO_G1.2;Parent=DEMO_G1")
    row("chrD1", "exon", ex1_s, ex1_e, "+", "ID=DEMO_G1.2:e1;Parent=DEMO_G1.2")
    row("chrD1", "CDS", cds1_s, cds1_s + 3 * 150 + 2, "+",
        "ID=DEMO_G1.2:c1;Parent=DEMO_G1.2", "0")

    row("chrD2", "gene", g2_start, g2_end, "-", "ID=DEMO_G2;Name=DEMO_G2")
    row("chrD2", "mRNA", g2_start, g2_end, "-", "ID=DEMO_G2.1;Parent=DEMO_G2")
    row("chrD2", "exon", g2_start, g2_end, "-", "ID=DEMO_G2.1:e1;Parent=DEMO_G2.1")
    row("chrD2", "CDS", c2_lo, c2_hi, "-", "ID=DEMO_G2.1:c1;Parent=DEMO_G2.1", "0")

    with open(os.path.join(outdir, "demo.gff3"), "w", encoding="utf-8") as f:
        f.write("\n".join(gff) + "\n")

    print(f"demo data -> {outdir}")
    print(f"  DEMO_G1.1  CDS {len(cds1)} bp  ({len(cds1)//3-1} aa), 2 exons, + strand")
    print(f"  DEMO_G1.2  CDS {len(cds1_short)} bp (short isoform)")
    print(f"  DEMO_G2.1  CDS {len(cds2)} bp  ({len(cds2)//3-1} aa), 1 exon, - strand, "
          f"internal BsaI sites")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "demo_data")
