---
name: GSprimer
description: |
  CDS-amplification primer design skill for GS / Golden Gate vector cloning
  (BsaI & PaqCI dual schemes). For a full-length CDS, the forward primer is
  anchored at the ATG and the reverse primer at the stop codon, keeping the
  reading frame on a 3-nt grid, emitting order-ready primer sequences (with
  adapters) plus a five-module risk report (specificity / GC / internal
  type-IIS sites / frame integrity / amplification feasibility).

  Trigger phrase (must match exactly):
  - GSprimer (F skills)

  Features:
  - Three template sources: GFF3 + genome FASTA / cDNA or CDS multi-FASTA / raw
    sequence input
  - Two-stage workflow: design (candidate enumeration + specificity + risk
    assessment) -> user confirmation -> finalize (attach adapters + order sheet)
  - Both cloning schemes delivered at once: scarless (BsaI) and Golden Gate
    (PaqCI/AarI), sharing the same GS destination vector
  - Dual-strand internal type-IIS site scan + synonymous rescue suggestions
  - Intron-aware template advice (cDNA vs gDNA) + concrete PCR recipe advice
  - Primer thermodynamics and all QC thresholds inherited from SpacerFinder
    v2.4.5 (the only deviation: no amplicon-size limit)

  Dependencies: pure Python (>=3.8), no required third-party libraries; BLAST+
  is optional (used only for specificity checks)

agent_created: true
version: 1.2.0
---

# GSprimer — CDS-amplification primer design for GS / Golden Gate vectors

## Overview

GSprimer designs primers that amplify a **full coding sequence (CDS)** so it can
be cloned into a **Golden Gate / GS-type vector**.

- The **forward primer** is 5′-anchored at the **ATG** (small shifts allowed, but
  `ATG_index − F_start_index` must be a multiple of 3 to preserve the frame).
- The **reverse primer** is 3′-anchored on the **last base of the stop codon**
  (only the 5′ end floats, length 19–24 nt).
- At finalize time GSprimer attaches **both** adapter schemes at once: the
  scarless **BsaI** scheme and the **Golden Gate (PaqCI/AarI)** scheme, which
  share the same GS destination vector (identical 4-nt sticky ends `cgag`/`ggat`).
  You choose which to order based on your destination vector.

> Primer thermodynamics and all quality-control filters are inherited verbatim
> from **SpacerFinder v2.4.5**; the only deliberate deviation is that GSprimer
> imposes **no amplicon-size window** (it is fixed by the ORF).

## Use cases

| Need | GSprimer |
|---|---|
| Clone a full CDS into a GS / Golden Gate vector | ✅ ATG…stop anchor, frame preserved |
| Choose between BsaI and PaqCI assembly | ✅ both adapter sets delivered |
| Audit a CDS amplicon for internal type-IIS sites | ✅ BsaI **and** PaqCI, both strands |
| Decide gDNA vs cDNA template | ✅ intron-aware template advice |
| Need a concrete PCR recipe (high GC, long fragment…) | ✅ auto-generated executable advice |
| No private dependency / no BLAST install | ✅ pure Python; BLAST optional |

## Environment & running

Skill root: `~/.workbuddy/skills/GSprimer/`

**No installation required** — invoke it as a module from the skill root:

```bash
cd ~/.workbuddy/skills/GSprimer
python -m gsprimer.cli --help
```

On Windows, to pin the interpreter, use the managed Python:

```
C:\Users\Administrator\.workbuddy\binaries\python\versions\3.13.12\python.exe -m gsprimer.cli --help
```

To get a runnable `gsprimer` command, install into a **virtual environment**
(do not pollute the global environment):

```bash
cd ~/.workbuddy/skills/GSprimer
python -m venv .venv && .venv/Scripts/pip install -e .
```

**Config search order** (`config.py`, only known keys are recognised):
file from env var `GSPRIMER_CONFIG` → `~/.gsprimer/config.yaml` →
`~/.config/gsprimer/config.yaml` → `gsprimer.yaml` in the current directory →
`~/.spacerfinder/config.yaml` (SpacerFinder compatibility — reuse an existing lab
config). Path-style env vars: `GSPRIMER_GENOME / GFF3 / CDNA / CDS / BLAST_DB /
BLASTN` (plus `SPACERFINDER_*` alias keys).

## Two-stage workflow (must be followed)

```bash
# Stage 1 — design + specificity + risk report
# For intron-bearing genes: --blast-db gives the genome DB for primer
#   specificity scanning; --blast-db-cdna gives the transcriptome DB for
#   in-silico PCR (avoids mis-flagging the target as an off-target).
python -m gsprimer.cli --gene LOC_Os01g01010 \
    --gff3 annotation.gff3 --genome genome.fa \
    --blast-db genome_db --blast-db-cdna cdna_db --outdir out

# (review the candidate table, pick pair N)

# Stage 2 — attach GS adapters, emit the order sheet + final report
python -m gsprimer.cli --finalize --select 1 --outdir out
```

**⚠️ Confirmation is mandatory.** At Stage 2 GSprimer first prints the selected
F/R **core sequences** (with anchoring positions and internal BsaI / PaqCI site
counts), then waits for an explicit `y` before attaching adapters. Use `--yes`
only in automated pipelines that have already been vetted by a human.

**Agent invocation rules:**

1. Run Stage 1 first and **present** the candidate table (number, Tier, anchor,
   Tm, GC, ΔTm, amplicon length, internal BsaI/PaqCI site counts, risk level) to
   the user.
2. Let the **user** decide which pair — do not pick on their behalf.
3. After receiving the number, run `--finalize --select N`; **omit `--yes` by
   default**, adding it only when the user explicitly asks to skip confirmation.
4. Before ordering, echo the full primer sequences from `order.tsv` / the report
   back to the user for a second check.

## Three template input modes

1. **GFF3 + genome FASTA** (default for annotated genomes). Multi-isoform genes
   → longest-CDS transcript; supports +/− strands and spliced mRNA.
2. **cDNA / CDS multi-FASTA** — looked up by ID (cDNA uses longest ATG..stop ORF
   inference).
3. **Raw sequence** — paste the raw sequence; ORF is inferred when no full
   ATG..stop is present.

Bundled demo data (no BLAST needed):

```bash
python -m gsprimer.cli --gene DEMO_G1.1 \
    --cds tests/demo_data/demo_cds.fa --outdir out_demo --no-blast
python -m gsprimer.cli --finalize --select 1 --outdir out_demo --yes
```

## Design principles (inherited thresholds)

| Parameter | Value |
|---|---|
| Length | 19–24 nt |
| GC% | 40–60% |
| Tm (Wallace + salt/Mg²⁺ correction) | 55–60 ℃ |
| \|ΔTm\| (F vs R) | < 2 ℃ |
| Self-dimer ΔG | ≥ −6.0 kcal·mol⁻¹ (−5.0 at the 3′ end) |
| Hairpin ΔG | ≥ −4.0 kcal·mol⁻¹ |

**GS-specific anchoring**

- Forward: 5′ anchored at ATG; shift must be a multiple of 3; shifting into the
  CDS (N-terminal loss) is heavily penalised and never reaches Tier A/B.
- Reverse: **5′ end** anchored on the stop codon's last base; the 3′ end is the
  extension end (length 19–24 nt determines where it lands inside the CDS);
  `--r-ext` allows the 5′ end to extend into the 3′UTR.
- The candidate set **always includes the strict ATG..stop pair** (flagged
  `★ standard anchor pair`).
- Tier: A = all inherited filters pass; B = soft boundary; C = needs trade-offs
  (structural defects that break the anchoring contract are always C).

## Two compatible cloning schemes

Same GS destination vector, identical 4-nt sticky ends (`cgag` left / `ggat`
right), differing only in the type-IIS enzyme; at finalize time **both** sets are
attached and listed.

| Scheme | Enzyme | Recognition site | Adapter |
|---|---|---|---|
| Scarless | BsaI (GGTCTC, 1/5 cut) | `GGTCTC` | F: `gtgatatcAGGTCTCTcgag` / R: `gccgcgggTGGTCTCAatcc` |
| Golden Gate | PaqCI / AarI (CACCTGC, 4/8 cut) | `CACCTGC` | F: `agCACCTGCagtccgag` / R: `agCACCTGCagtcatcc` |

> The Golden Gate adapter leaves only 2 nt of 5′ protection upstream of
> `CACCTGC` (recommended ≥ 6 nt); the report flags this. The sticky ends are
> correct and compatible.

## Risk report (the 5 required modules)

1. **Specificity** — BLAST / in-silico-PCR level (LOW / MEDIUM / HIGH / UNKNOWN).
2. **GC content** — overall + 50-nt sliding window + first/last 100 nt.
3. **Internal type-IIS sites** — **dual-strand** scan for **both schemes**: BsaI
   (`GGTCTC`/`GAGACC`) and PaqCI/AarI (`CACCTGC`/`GCAGGTG`); gives position,
   codon, and **synonymous wobble rescue** options.
4. **Frame integrity** — offset, N-terminal truncation, premature in-frame stops.
5. **Amplification feasibility** — intron-aware template choice (cDNA vs gDNA)
   **plus** concrete PCR recipe advice (high/low GC handling, long-fragment
   extension, primer Tm mismatch, homopolymers).

For the final construct it also audits that the full oligo contains exactly the
two type-IIS sites of the chosen scheme (internal sites would be cut during
Golden Gate).

## CLI reference

| Flag | Meaning |
|---|---|
| `--gene / --transcript` | gene/transcript ID (or `--longest-by cds\|mrna`) |
| `--gff3 --genome` | annotation + genome FASTA |
| `--cdna / --cds` | cDNA / CDS multi-FASTA |
| `--sequence / --sequence-file` | raw sequence |
| `--max-shift` | max \|offset\| for the forward primer (multiple of 3, default 12) |
| `--r-ext` | max nt the reverse primer may cross past the stop codon (default 0 = strict) |
| `--top` | number of candidate pairs (default 10) |
| `--no-downstream` | forbid forward primer shifting into the CDS (avoids N-terminal truncation) |
| `--blast-db / --blastn` | BLAST database + executable |
| `--fallback-fasta` | cDNA/transcriptome FASTA for seed-and-extend when BLAST is absent |
| `--no-blast` | skip specificity entirely (marks UNKNOWN) |
| `--no-intron` | declare the gene intronless → report recommends gDNA |
| `--has-introns` | declare the gene has introns → force cDNA template |
| `--adapter-f / --adapter-r` | override GS adapters |
| `--select / --finalize` | choose candidate pair(s); Stage 2 lists primers and waits for confirmation |
| `--yes` | skip the interactive confirmation prompt (**scripted / already-vetted flows only**) |
| `--scheme` | terminal print filter: `all` (default, both sets) / `seamless` / `golden_gate` |

## Output files (per run)

| File | Content |
|---|---|
| `<gene>_GSprimer.html` | interactive report (Fv1.4 palette: OE red `#B2182B` / CR blue `#2166AC` / neutral grey `#999999`) |
| `<gene>_GSprimer.tsv` | candidate table |
| `gsprimer_state.json` | cached design state (lets `--finalize` skip recomputation) |
| `<gene>_GSprimer_final.html` + `<gene>_order.tsv` | after `--finalize`; order.tsv lists 4 rows per pair (F/R × seamless + Golden Gate) |

## Programmatic use

```python
from gsprimer import run_design, run_finalize

res = run_design("LOC_Os01g01010", cds_fa="MSU_cds.fa", outdir="out")
fin = run_finalize([1], outdir="out")
print(fin["html"])   # final report path
print(fin["order"])  # order sheet path
```

## Module map

```
GSprimer/
  SKILL.md        ← this file
  gsprimer/
    config.py     config search (env > ~/.gsprimer > ~/.spacerfinder > defaults)
    thermo.py     primer thermodynamics inherited from SpacerFinder
    template.py   Transcript model + GFF3/cDNA/CDS/direct loaders
    design.py     anchored candidate enumeration + A/B/C tiering + scoring
    blast.py      BLAST + in-silico PCR + seed-and-extend fallback
    risk.py       specificity/GC/type-IIS/frame/amplification + rescue + PCR recipe
    adapters.py   dual-scheme GS adapter parsing (BsaI / PaqCI) + order sheet
    report.py     HTML + TSV report
    pipeline.py   run_design / run_finalize orchestration
    cli.py        command line
    interactive.py interactive selection
  tests/demo_data/  self-contained smoke-test dataset
  scripts/make_demo_data.py
  docs/SECURITY_AUDIT.md  security audit report for installation
```

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError: gsprimer` | current dir is not the skill root (where the package dir lives) | `cd ~/.workbuddy/skills/GSprimer` then `python -m gsprimer.cli`; or use `pip install -e .` inside a venv |
| **adapter case got mangled** | in `adapters.py`, case is **meaningful**: lowercase = protection bases, uppercase = enzyme site / sticky end | never `.upper()` the `ADAPTER_F/ADAPTER_R` constants in any preprocessing; override with `--adapter-f/--adapter-r` verbatim |
| specificity shows `UNKNOWN` in the report | no usable BLAST DB and no `--fallback-fasta` | provide `--blast-db` (after `makeblastdb`) or `--fallback-fasta` (seed-and-extend fallback); if neither, honestly mark UNKNOWN — do not treat it as LOW |
| `blastn: command not found` | BLAST+ not on PATH | pass `--blastn /path/to/blastn` explicitly, or `--no-blast` (marks UNKNOWN) |
| forward primer cannot move up into the 5′UTR | CDS-FASTA mode has no UTR; GFF3 missing `exon` records falls back to CDS as exon | to capture the 5′UTR use GFF3 + genome, or cDNA multi-FASTA |
| no Tier A in candidates | strict ATG..stop anchor + 40–60% GC + 55–60 ℃ Tm is inherently strict (especially when CDS 5′ GC is extreme) | look at Tier B/C and read the flags; widen with `--max-shift` into the 5′UTR if needed, or check for synonymous wobble space |
| Golden Gate assembly fails | internal type-IIS site remains in the amplicon | read module 3's `rescuable` + synonymous mutation suggestions in the report, edit the CDS before ordering; check both schemes separately |
| `--finalize` reports state missing | `--outdir` differs from the design stage | keep `--outdir` identical; state file is `gsprimer_state.json` |
| internal site count disagrees with expectation | only one strand scanned | this skill scans both strands (`scan_sites` scans each `site` and its reverse complement); keep this if you modify it |
| version number mismatch | `pyproject.toml` says 1.1.0, `gsprimer/__init__.py` says 1.2.0 | `__init__.__version__` is authoritative; unify when releasing upstream |
| **all candidates flag specificity HIGH but `Spec_offtarget_products=0`** (contradiction) | **BLAST+ ≥ 2.16 parameter conflict**: `-max_target_seqs` and `-num_alignments` cannot both be passed; `blastn` returns rc=1 + empty output, and old code silently treated empty output as "0 hits" | fixed: only `-max_target_seqs` is passed and `returncode` is checked; on failure it degrades to UNKNOWN rather than HIGH |
| in-silico PCR lists the **target itself** as a perfect off-target → specificity HIGH | with a **genome** BLAST DB, the gDNA product contains introns (e.g. 564 bp CDS → 2.3 kb), never == `expect_len` → `is_target=False` | add `--blast-db-cdna` (after `makeblastdb -in MSU_cdna.fa -dbtype nucl`); `pipeline.run_design` auto-switches when it detects a multi-exon gene and this DB is supplied |
| an isoform of the same gene counted as off-target → still HIGH | old `insilico_pcr()` only recognised `is_target` (exact length match); same-gene different-length products fell into off-target | fixed: added `target_gene` parameter and `is_isoform` flag; same-gene products (isoforms, intron-containing gDNA products) all count as on-target |
| output dir cannot be deleted (`[safe-delete][SAFE_DELETE_FAIL_CLOSED] trash-failed`) | local safe-delete protection blocks it | do not delete the dir — just re-run in place to overwrite (same `--outdir` overwrites html/tsv/state) |

## Security statement

This skill is **fully local computation**: no network requests, no `eval/exec`,
no remote script download/execution, no credential reads. The only subprocess
call is the local `blastn` (list-form arguments, no `shell=True`, 120 s timeout),
optional and disableable.
