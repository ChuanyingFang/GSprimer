# GSprimer

**CDS / promoter amplification primer design for GS / Golden Gate (BsaI & PaqCI) vector cloning.**

GSprimer designs primers that amplify a full coding sequence (CDS) **or the
region upstream of a transcription start site (promoter / functional element)**
so it can be cloned into a **Golden Gate / GS-type vector**.

- **CDS mode** — the forward primer is anchored at the **ATG** and the reverse
  primer at the **stop codon**, with the reading frame held on a 3-nt grid, so the
  translated protein comes out complete.
- **Promoter mode (v2.0)** — the reverse primer is anchored at the **TSS**
  (1 nt before the 5′-UTR start); the forward primer is placed so the PCR product
  lands in a 1200–1700 bp window, with a two-file fallback into the transcript
  5′-UTR when no ideal R exists in the upstream region.

At finalize time GSprimer always emits **both** adapter schemes — a scarless
**BsaI** scheme and a **Golden Gate (PaqCI / AarI)** scheme — so you can pick the
one that matches your destination vector.

> Primer thermodynamics and all quality filters are inherited verbatim from
> **SpacerFinder**; the only deliberate deviation is that GSprimer imposes **no
> amplicon-size window** (it is fixed by the ORF).

---

## Why GSprimer

| Need | GSprimer |
|---|---|
| Clone a full CDS into a GS / Golden Gate vector | ✅ anchored ATG…stop, frame preserved |
| Choose between BsaI and PaqCI assembly | ✅ both adapter sets delivered, you decide |
| Audit a CDS amplicon for internal type-IIS sites | ✅ BsaI **and** PaqCI, both strands |
| Decide gDNA vs cDNA template | ✅ intron-aware template advice |
| Know the right PCR recipe (high GC, long fragment…) | ✅ auto-generated, concrete advice |
| No private dependency / no BLAST install | ✅ pure-Python; BLAST optional |

---

## Installation

```bash
git clone https://github.com/ChuanyingFang/GSprimer.git
cd GSprimer
pip install -e .
```

This installs the `gsprimer` command. There are **no required third-party
dependencies**. [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
is optional — used only for specificity checks. Without it, GSprimer falls back
to a seed-and-extend heuristic over a user-supplied cDNA FASTA, or marks
specificity as `UNKNOWN`.

---

## Quick start

```bash
# Stage 1 — design + specificity + risk report
python -m gsprimer.cli --gene LOC_Os01g01010 \
    --gff3 annotation.gff3 --genome genome.fa \
    --blast-db genome_db --outdir out

# (review the candidate table, pick pair N)

# Stage 2 — attach GS adapters, emit order sheet + final report
#   GSprimer lists the chosen primer CORE sequences and asks you to confirm
#   BEFORE any adapter is attached. Answer 'y' to proceed.
python -m gsprimer.cli --finalize --select 1 --outdir out
#   scripted pipelines: skip the prompt with --yes
python -m gsprimer.cli --finalize --select 1 --outdir out --yes
```

> **Confirmation is mandatory.** GSprimer never attaches adapters until you
> have reviewed the exact primer sequences it is about to order. At Stage 2 it
> prints the selected F/R core sequences (with anchoring and the internal
> BsaI / PaqCI site counts) and waits for an explicit `y` confirmation. Use
> `--yes` only in automated pipelines where the selection is already vetted.

### Three template input modes

1. **GFF3 + genome FASTA** (default for annotated genomes). Multi-isoform genes
   → longest-CDS transcript. Handles +/− strands and spliced mRNA.
2. **cDNA / CDS multi-FASTA** — looked up by ID (ORF found by longest ATG..stop
   for cDNA).
3. **Direct sequence** — paste a raw sequence; ORF inferred when no full
   ATG..stop is present.

### Minimal example with the bundled demo data (no BLAST)

```bash
python -m gsprimer.cli --gene DEMO_G1.1 \
    --cds tests/demo_data/demo_cds.fa --outdir out_demo --no-blast
# --yes bypasses the interactive primer confirmation (use it in this headless demo)
python -m gsprimer.cli --finalize --select 1 --outdir out_demo --yes
```

---

## Design principles

Thermodynamic rules are cloned verbatim from SpacerFinder:

| Parameter | Value |
|---|---|
| Length | 19–24 nt |
| GC% | 40–60% |
| Tm (Wallace + salt/Mg²⁺ correction) | 55–60 ℃ |
| \|ΔTm\| between F and R | < 2 ℃ |
| Self-dimer ΔG | ≥ −6.0 kcal·mol⁻¹ (−5.0 at the 3′ end) |
| Hairpin ΔG | ≥ −4.0 kcal·mol⁻¹ |

**GS-specific anchoring**

- **Forward primer** — 5′ end anchored at the **ATG**. Small shifts are allowed
  but `ATG_index − F_start_index` **must be a multiple of 3** (frame preserved).
  A shift into the CDS (N-terminal loss) is heavily penalised and never Tier A/B.
- **Reverse primer** — 3′ end anchored on the **last base of the stop codon**.
  Only the 5′ end floats (length 19–24 nt).

The candidate set **always includes the strict ATG..stop pair** (flagged
`★ standard anchor pair`).

### Promoter mode (v2.0)

Amplify the region **upstream of the TSS** (promoter / enhancer / other functional
element) instead of a CDS. Invoke with `--promoter` (a ~2 kb upstream sequence
whose 3′ end is the TSS) and, optionally, `--utr5` (the transcript 5′-UTR, used
only as a fallback):

```bash
python -m gsprimer.cli --promoter upstream_2k.fa --utr5 utr5.fa \
    --product-min 1200 --product-max 1700 --outdir out
```

**Anchoring rules**

- **Reverse primer** — 5′ end anchored at the **TSS** (1 nt before the 5′-UTR
  start). If no ideal R exists there, it may:
  - extend **≤ 100 nt upstream** of the 5′-UTR (into the promoter region), or
  - fall back **≤ 200 nt into the 5′-UTR** (using the `--utr5` transcript
    sequence) when the 2 kb upstream region yields no Tier A/B R.
- **Forward primer** — placed so the PCR product lands in the **1200–1700 bp**
  window; there is no reading frame in non-coding sequence.
- The candidate set **always includes the TSS-anchored pair** (flagged
  `★` in the report).

**Two-file strategy:** GSprimer first scans the upstream region for R; it only
falls back to the transcript 5′-UTR when no Tier A/B R is found there, keeping the
product as close to the native promoter as possible.

---

## Two compatible cloning schemes

Both schemes share the **same GS destination vector**: they generate identical
4-nt sticky ends (`attc` left / `gagc` right), differing only in the type-IIS
enzyme. At finalize time **both** sets are attached and listed so you choose
which to order.

**Scarless / BsaI (GGTCTC, 1/5 cut)**

- Forward: `gataagcttGGTCTCTattc`
- Reverse: `CATggatccGGTCTCAgctc`

**Golden Gate / PaqCI–AarI (CACCTGC, 4/8 cut)**

- Forward: `agCACCTGCagtcattc`
- Reverse: `agCACCTGCagtgctc`

> Both adapters carry ample 5′ protection, so the type-IIS sites are digested
> cleanly. The seamless and Golden Gate sets emit the **same** 4-nt sticky ends
> (`attc` / `gagc`), so either order drops into the same destination vector.

---

## Risk report (the 5 required blocks)

The HTML/TSV report always flags:

1. **Specificity** — BLAST / in-silico-PCR level (LOW / MEDIUM / HIGH / UNKNOWN).
2. **GC content** — overall + 50-nt sliding window + first/last 100 nt.
3. **Internal type-IIS sites** — counted on **both strands** for **both schemes**:
   BsaI (`GGTCTC` / `GAGACC`) and PaqCI/AarI (`CACCTGC` / `GCAGGTG`). Both
   produce the identical 4-nt sticky ends, so either adapter set drops into the
   same destination vector. Each scheme carries positions, codons, and
   **synonymous wobble rescue** options.
4. **Frame integrity** — offset, N-terminal truncation, premature in-frame stops.
5. **Amplification feasibility** — intron-aware template choice (cDNA vs gDNA)
   **plus concrete PCR recipe advice** (high/low GC handling, long-fragment
   extension, primer Tm mismatch, homopolymers).

For the final construct it also audits that the full ordered oligo has exactly
the expected type-IIS sites of the chosen scheme (the two adapters) — internal
sites would be cut during Golden Gate.

---

## Intron-aware template choice & PCR recipe advice

- **Gene has introns** (multi-exon) → must use **cDNA / RT-PCR**.
- **Gene has no introns** (single-exon, or `--no-intron`) → **gDNA recommended**:
  cDNA and gDNA give identical products, and gDNA skips reverse transcription.
- **Unknown** (CDS/cDNA/direct input) → reported as pending; use `--no-intron`
  when you know the gene is intronless.

PCR recipe advice is auto-generated from fragment + primer traits: high GC →
GC-enhanced polymerase + DMSO/betaine + touchdown; long fragment → long-range
polymerase; Tm mismatch → anneal at the lower Tm or 3-step protocol; homopolymer
runs → design around them.

---

## CLI reference

| Flag | Meaning |
|---|---|
| `--gene / --transcript` | gene/transcript ID (or `--longest-by cds\|mrna`) |
| `--gff3 --genome` | annotation + genome FASTA |
| `--cdna / --cds` | cDNA / CDS multi-FASTA |
| `--sequence / --sequence-file` | raw sequence |
| `--promoter / --promoter-file` | upstream promoter region (3′ end = TSS), mutually exclusive with `--gene/--sequence` |
| `--utr5 / --utr5-file` | optional transcript 5′-UTR fallback for promoter mode |
| `--product-min / --product-max` | PCR product window for promoter mode (default 1200 / 1700) |
| `--r-up-max` | max nt the R 5′ end may extend upstream of the TSS in promoter mode (default 100) |
| `--r-utr-max` | max nt R may extend into the 5′-UTR fallback (default 200) |
| `--max-shift` | max \|offset\| for the forward primer (multiple of 3) |
| `--top` | number of candidate pairs (default 10) |
| `--blast-db / --blastn` | BLAST database + executable |
| `--fallback-fasta` | cDNA/transcriptome FASTA for seed-and-extend when BLAST absent |
| `--no-blast` | skip specificity entirely (marks UNKNOWN) |
| `--no-intron` | declare the gene is intronless → report recommends gDNA |
| `--has-introns` | declare the gene has introns → force cDNA template |
| `--adapter-f / --adapter-r` | override GS adapters |
| `--select / --finalize` | choose pair(s); Stage 2 then lists the primers and asks for confirmation before attaching adapters |
| `--yes` | skip the interactive confirmation prompt (adapters attached without asking; scripting only) |
| `--scheme` | terminal print filter: `all` (default, both sets), `seamless`, or `golden_gate` |

---

## Output files (per run)

- `<gene>_GSprimer.html` — interactive report.
- `<gene>_GSprimer.tsv` — candidate table.
- `gsprimer_state.json` — cached design (lets `--finalize` skip recomputation).
- `<gene>_GSprimer_final.html` + `<gene>_order.tsv` — after `--finalize`.
  The `order.tsv` lists **4 rows per pair** (F/R × seamless + Golden Gate) so
  both schemes are ready to order; the HTML cards show both sets of full ordered
  sequences for the user to choose.

---

## Programmatic use

```python
from gsprimer import run_design, run_finalize

res = run_design("LOC_Os01g01010", cds_fa="MSU_cds.fa", outdir="out")
fin = run_finalize([1], outdir="out")
print(fin["html"])   # final report path
print(fin["order"])  # order sheet path
```

```python
# Promoter mode: amplify the 2 kb upstream region (TSS at its 3' end),
# with the transcript 5'-UTR as an optional fallback for the reverse primer.
res = run_design("GeneX", promoter_seq=upstream_2k, utr5_seq=utr5, outdir="out")
fin = run_finalize([1], outdir="out")
```

---

## Module map

```
GSprimer/
  gsprimer/
    config.py     config search (env > ~/.gsprimer > ~/.spacerfinder > defaults)
    thermo.py     SpacerFinder-primer thermodynamics (cloned verbatim)
    template.py   Transcript model + GFF3/cDNA/CDS/direct loaders
    design.py     anchored candidate enumeration + A/B/C tiering + scoring
    blast.py      BLAST + in-silico PCR + seed-and-extend fallback
    risk.py       specificity/GC/type-IIS(BsaI+PaqCI)/frame/amplification + rescue + PCR recipe
    adapters.py   dual-scheme GS adapter parsing (BsaI / PaqCI) + order sheet
    report.py     HTML + TSV report
    pipeline.py   run_design / run_finalize orchestration
    cli.py        command line
    interactive.py selection prompt
  tests/demo_data/  self-contained smoke-test dataset
  scripts/         make_demo_data.py
```

---

## License

MIT — see [LICENSE](LICENSE).
