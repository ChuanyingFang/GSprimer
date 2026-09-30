"""
GSprimer — Template acquisition.

Builds a unified `Transcript` model from any of three input modes:

  1. GFF3 mode   : genome FASTA + GFF3  -> full mRNA (exon-spliced) + CDS bounds
                   Multi-isoform genes default to the LONGEST transcript.
  2. FASTA mode  : cDNA and/or CDS multi-FASTA, looked up by ID
  3. Direct mode : a raw sequence string / single FASTA file

A fourth construction, `from_promoter`, builds a promoter-mode Transcript used
by the v2.0 promoter-amplification workflow (upstream 2 kb + optional 5'-UTR
fallback). Promoter mode carries no ORF; TSS is marked by `utr5_start`.

Coordinate conventions
----------------------
  * `mrna` is the mature transcript, already in 5'->3' sense orientation.
  * `cds_start` / `cds_end` are 0-based, half-open indices INTO `mrna`,
     i.e. mrna[cds_start:cds_end] is the ORF including the stop codon.
  * `exons` are 1-based inclusive genome intervals in BIOLOGICAL order
     (first transcribed exon first), used to map mRNA index -> genome coord.
"""

import re
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .thermo import revcomp

STOP_CODONS = {"TAA", "TAG", "TGA"}


# ============================================================
# Genome
# ============================================================

class Genome:
    """Multi-FASTA loader with optional chromosome filtering."""

    def __init__(self, fasta_path: str = "", only: Optional[set] = None):
        self.seq: Dict[str, str] = {}
        self.sizes: Dict[str, int] = {}
        if fasta_path:
            self.load(fasta_path, only)

    def load(self, fasta_path: str, only: Optional[set] = None) -> None:
        t0 = time.time()
        chrom = None
        keep = True
        buf: List[str] = []
        with open(fasta_path, encoding="utf-8", errors="replace") as f:
            for raw in f:
                if not raw:
                    continue
                if raw[0] == ">":
                    if chrom is not None and keep:
                        self.seq[chrom] = "".join(buf)
                        self.sizes[chrom] = len(self.seq[chrom])
                    chrom = raw[1:].strip().split()[0]
                    keep = (only is None) or (chrom in only)
                    buf = []
                elif keep:
                    buf.append(raw.strip().upper())
        if chrom is not None and keep:
            self.seq[chrom] = "".join(buf)
            self.sizes[chrom] = len(self.seq[chrom])
        total = sum(self.sizes.values())
        print(f"[Genome] {len(self.seq)} seq, {total:,} bp "
              f"({time.time() - t0:.1f}s)", file=sys.stderr)

    def subseq(self, chrom: str, start: int, end: int) -> str:
        """1-based inclusive slice."""
        s = self.seq.get(chrom, "")
        start = max(1, start)
        end = min(end, len(s))
        if end < start:
            return ""
        return s[start - 1:end]


# ============================================================
# Transcript model
# ============================================================

@dataclass
class Transcript:
    tx_id: str
    gene_id: str = ""
    chrom: str = ""
    strand: str = "+"
    mrna: str = ""              # mature transcript, sense orientation
    cds_start: int = 0          # 0-based index into mrna (A of ATG)
    cds_end: int = 0            # 0-based exclusive, past the stop codon
    exons: List[Tuple[int, int]] = field(default_factory=list)  # biological order
    source: str = ""            # gff3 | cdna | cds | direct | promoter
    mode: str = "cds"           # "cds" | "promoter"
    utr5_start: int = 0         # 0-based TSS position within mrna (promoter mode)
    notes: List[str] = field(default_factory=list)

    # ---------- derived ----------
    @property
    def cds(self) -> str:
        return self.mrna[self.cds_start:self.cds_end]

    @property
    def cds_len(self) -> int:
        return self.cds_end - self.cds_start

    @property
    def utr5_len(self) -> int:
        return self.cds_start

    @property
    def utr3_len(self) -> int:
        return len(self.mrna) - self.cds_end

    @property
    def protein_len(self) -> int:
        return max(0, self.cds_len // 3 - 1)

    def validate(self) -> List[str]:
        """Structural sanity checks; branches on CDS vs promoter mode."""
        if self.mode == "promoter":
            return self._validate_promoter()
        msgs = []
        cds = self.cds
        if not cds:
            msgs.append("CDS is empty")
            return msgs
        if not cds.startswith("ATG"):
            msgs.append(f"CDS does not start with ATG (found {cds[:3]})")
        if cds[-3:] not in STOP_CODONS:
            msgs.append(f"CDS does not end with a stop codon (found {cds[-3:]})")
        if len(cds) % 3 != 0:
            msgs.append(f"CDS length {len(cds)} is not a multiple of 3")
        internal = [i for i in range(0, len(cds) - 3, 3)
                    if cds[i:i + 3] in STOP_CODONS]
        if internal:
            msgs.append(f"CDS contains {len(internal)} internal in-frame stop "
                        f"codon(s) (first at codon {internal[0] // 3 + 1})")
        if "N" in cds:
            msgs.append(f"CDS contains {cds.count('N')} ambiguous base(s)")
        return msgs

    def _validate_promoter(self) -> List[str]:
        """Light structural checks for promoter mode (non-coding, ORF-free)."""
        msgs = []
        if not (0 <= self.utr5_start <= len(self.mrna)):
            msgs.append(f"TSS position {self.utr5_start} is outside the sequence "
                        f"(length {len(self.mrna)})")
        if self.utr5_start < 200:
            msgs.append("Upstream promoter sequence is short (<200 nt); the "
                        "1200–1700 bp product window may not be satisfiable")
        return msgs

    def mrna_to_genome(self, idx: int) -> int:
        """Map 0-based mRNA index -> 1-based genome coordinate (0 if unknown)."""
        if not self.exons:
            return 0
        remaining = idx
        for (s, e) in self.exons:
            length = e - s + 1
            if remaining < length:
                return (s + remaining) if self.strand == "+" else (e - remaining)
            remaining -= length
        return 0


# ============================================================
# GFF3 parsing
# ============================================================

_ATTR_RE = re.compile(r"([^=;]+)=([^;]*)")


def _parse_attrs(s: str) -> Dict[str, str]:
    return {m.group(1).strip(): m.group(2).strip()
            for m in _ATTR_RE.finditer(s)}


class Annotation:
    """Minimal GFF3 model: gene -> mRNA -> exon/CDS."""

    def __init__(self):
        self.tx_info: Dict[str, Dict] = {}          # tx_id -> meta
        self.tx_exons: Dict[str, List] = defaultdict(list)
        self.tx_cds: Dict[str, List] = defaultdict(list)
        self.gene_tx: Dict[str, List[str]] = defaultdict(list)
        self.gene_alias: Dict[str, str] = {}        # Name -> gene ID

    def parse(self, gff3_path: str, gene_filter: Optional[set] = None) -> None:
        """Parse GFF3. gene_filter restricts to specific gene/transcript IDs
        (substring match on the attribute block) for speed on large files."""
        t0 = time.time()
        n_line = 0
        with open(gff3_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line or line[0] == "#":
                    continue
                if gene_filter and not any(g in line for g in gene_filter):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                n_line += 1
                chrom, _src, ftype, start, end, _sc, strand, phase, attrs = cols[:9]
                if ftype not in ("mRNA", "transcript", "exon", "CDS", "gene"):
                    continue
                a = _parse_attrs(attrs)
                start, end = int(start), int(end)
                if ftype == "gene":
                    gid = a.get("ID", "")
                    name = a.get("Name", "")
                    if gid and name:
                        self.gene_alias[name] = gid
                elif ftype in ("mRNA", "transcript"):
                    tid = a.get("ID", "")
                    gid = a.get("Parent", "").split(",")[0]
                    if not tid:
                        continue
                    self.tx_info[tid] = {"gene": gid, "chrom": chrom,
                                         "strand": strand,
                                         "start": start, "end": end}
                    self.gene_tx[gid].append(tid)
                    if name := a.get("Name", ""):
                        self.gene_alias.setdefault(name, tid)
                elif ftype == "exon":
                    for p in a.get("Parent", "").split(","):
                        if p:
                            self.tx_exons[p].append((start, end))
                elif ftype == "CDS":
                    ph = 0 if phase in (".", "") else int(phase)
                    for p in a.get("Parent", "").split(","):
                        if p:
                            self.tx_cds[p].append((start, end, ph))
        print(f"[GFF3] {len(self.tx_info)} transcripts, "
              f"{len(self.gene_tx)} genes ({time.time() - t0:.1f}s)",
              file=sys.stderr)

    # ---------------------------------------------------------
    def transcripts_of(self, query: str) -> List[str]:
        """Resolve a gene ID / transcript ID / Name to transcript IDs."""
        if query in self.tx_info:
            return [query]
        if query in self.gene_tx:
            return list(self.gene_tx[query])
        if query in self.gene_alias:
            resolved = self.gene_alias[query]
            if resolved in self.gene_tx:
                return list(self.gene_tx[resolved])
            if resolved in self.tx_info:
                return [resolved]
        # loose match: LOC_Os01g01010 -> LOC_Os01g01010.1
        hits = [t for t in self.tx_info if t.split(".")[0] == query]
        if hits:
            return sorted(hits)
        hits = [g for g in self.gene_tx if g.split(".")[0] == query]
        if hits:
            out = []
            for g in hits:
                out.extend(self.gene_tx[g])
            return sorted(out)
        return []

    def build(self, genome: Genome, tx_id: str) -> Optional[Transcript]:
        """Assemble a Transcript from parsed features."""
        meta = self.tx_info.get(tx_id)
        cds_iv = sorted(self.tx_cds.get(tx_id, []), key=lambda x: x[0])
        if not cds_iv:
            return None
        if meta is None:  # CDS-only annotation without an mRNA line
            meta = {"gene": tx_id.split(".")[0], "chrom": "", "strand": "+"}
        chrom = meta["chrom"]
        strand = meta["strand"]

        exon_iv = sorted(self.tx_exons.get(tx_id, []), key=lambda x: x[0])
        notes = []
        if not exon_iv:
            exon_iv = [(s, e) for s, e, _ in cds_iv]
            notes.append("GFF3 has no exon records; CDS used instead (no UTR, "
                         "forward primer cannot shift upstream)")

        bio_exons = exon_iv if strand == "+" else exon_iv[::-1]

        # mature mRNA
        parts = []
        for (s, e) in exon_iv:
            parts.append(genome.subseq(chrom, s, e))
        mrna = "".join(parts)
        if strand == "-":
            mrna = revcomp(mrna)

        # ATG genome coordinate (biological first CDS base)
        if strand == "+":
            atg_g = min(s for s, _e, _p in cds_iv)
            stop_g = max(e for _s, e, _p in cds_iv)
            first_phase = [p for s, _e, p in cds_iv if s == atg_g][0]
        else:
            atg_g = max(e for _s, e, _p in cds_iv)
            stop_g = min(s for s, _e, _p in cds_iv)
            first_phase = [p for _s, e, p in cds_iv if e == atg_g][0]
        if first_phase != 0:
            notes.append(f"first CDS phase={first_phase} (annotation may be "
                         f"5'-incomplete)")

        cds_start = _genome_to_mrna(bio_exons, strand, atg_g)
        cds_total = sum(e - s + 1 for s, e, _p in cds_iv)
        if cds_start < 0:
            return None
        cds_end = cds_start + cds_total

        # some annotations exclude the stop codon from CDS -> try to extend
        tx = Transcript(tx_id=tx_id, gene_id=meta.get("gene", ""),
                        chrom=chrom, strand=strand, mrna=mrna,
                        cds_start=cds_start, cds_end=cds_end,
                        exons=bio_exons, source="gff3", notes=notes)
        if tx.cds[-3:] not in STOP_CODONS and cds_end + 3 <= len(mrna):
            if mrna[cds_end:cds_end + 3] in STOP_CODONS:
                tx.cds_end = cds_end + 3
                tx.notes.append("CDS annotation lacks a stop codon; "
                                "auto-extended 3 nt downstream")
        _ = stop_g
        return tx


def _genome_to_mrna(bio_exons: List[Tuple[int, int]], strand: str,
                    gpos: int) -> int:
    """1-based genome coord -> 0-based mRNA index (-1 if not in an exon)."""
    off = 0
    for (s, e) in bio_exons:
        if s <= gpos <= e:
            return off + ((gpos - s) if strand == "+" else (e - gpos))
        off += e - s + 1
    return -1


# ============================================================
# FASTA-based modes
# ============================================================

def read_fasta(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    name = None
    buf: List[str] = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line:
                continue
            if line[0] == ">":
                if name:
                    out[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip().upper())
    if name:
        out[name] = "".join(buf)
    return out


def _lookup(db: Dict[str, str], query: str) -> List[Tuple[str, str]]:
    if query in db:
        return [(query, db[query])]
    hits = [(k, v) for k, v in db.items() if k.split()[0] == query]
    if hits:
        return hits
    hits = [(k, v) for k, v in db.items() if k.split(".")[0] == query]
    return sorted(hits)


def from_cds_fasta(path: str, query: str) -> List[Transcript]:
    """Build transcripts from a CDS multi-FASTA (no UTR available)."""
    db = read_fasta(path)
    out = []
    for name, seq in _lookup(db, query):
        out.append(Transcript(
            tx_id=name, gene_id=name.split(".")[0], mrna=seq,
            cds_start=0, cds_end=len(seq), source="cds",
            notes=["template from CDS FASTA: no UTR sequence, forward primer can "
                   "only shift downstream"]))
    return out


def from_cdna_fasta(path: str, query: str) -> List[Transcript]:
    """Build transcripts from a cDNA multi-FASTA; ORF found by longest ATG..stop."""
    db = read_fasta(path)
    out = []
    for name, seq in _lookup(db, query):
        s, e = longest_orf(seq)
        if s < 0:
            continue
        out.append(Transcript(
            tx_id=name, gene_id=name.split(".")[0], mrna=seq,
            cds_start=s, cds_end=e, source="cdna",
            notes=["template from cDNA FASTA: ORF inferred from longest "
                   "ATG..stop"]))
    return out


def utr5_from_transcript(seq: str) -> str:
    """Extract the 5'-UTR (sequence before the CDS) from a cDNA/mRNA string.

    Used by the promoter-mode UTR fallback: when no ideal R primer is found
    within the 2 kb upstream region, up to 200 nt are taken from here to extend
    the R primer's 5' end into the 5'-UTR.
    """
    s = "".join(seq.split()).upper().replace("U", "T")
    cs, _ce = longest_orf(s)
    if cs < 0:
        return ""
    return s[:cs]


def from_promoter(upstream_seq: str, utr5_seq: str = "") -> Transcript:
    """Build a promoter-mode Transcript.

    `upstream_seq` is the promoter region upstream of the TSS (typically 2 kb,
    5'->3' on the sense strand, its tail being the TSS); `utr5_seq` is the
    5'-UTR from the transcript (optional, used only for the UTR fallback when no
    ideal R primer sits within the 2 kb). The two are concatenated into `mrna`,
    with the TSS at the join.
    """
    up = "".join(upstream_seq.split()).upper().replace("U", "T")
    utr = "".join(utr5_seq.split()).upper().replace("U", "T")
    mrna = up + utr
    notes = [f"promoter mode: upstream region {len(up)} bp, UTR fallback "
             f"{len(utr)} bp"]
    if not utr:
        notes.append("no UTR sequence supplied; R primer cannot fall back into "
                     "the 5'-UTR (only upstream extension within the promoter "
                     "region is allowed)")
    return Transcript(tx_id="promoter", gene_id="promoter", chrom="", strand="+",
                     mrna=mrna, cds_start=len(up), cds_end=len(mrna),
                     source="promoter", mode="promoter",
                     utr5_start=len(up), notes=notes)


def longest_orf(seq: str) -> Tuple[int, int]:
    """Return (start, end) of the longest ATG..stop ORF; (-1,-1) if none."""
    s = seq.upper()
    best = (-1, -1)
    best_len = 0
    for frame in range(3):
        i = frame
        while i < len(s) - 2:
            if s[i:i + 3] == "ATG":
                j = i
                while j < len(s) - 2:
                    if s[j:j + 3] in STOP_CODONS:
                        if (j + 3 - i) > best_len:
                            best_len = j + 3 - i
                            best = (i, j + 3)
                        break
                    j += 3
                i = j + 3 if j > i else i + 3
            else:
                i += 3
    return best


def from_sequence(seq: str, name: str = "user_seq") -> Transcript:
    """Build a transcript from a raw sequence pasted by the user."""
    s = "".join(seq.split()).upper()
    s = re.sub(r"[^ACGTNU]", "", s).replace("U", "T")
    notes = []
    if s.startswith("ATG") and s[-3:] in STOP_CODONS and len(s) % 3 == 0:
        cs, ce = 0, len(s)
        notes.append("input sequence recognised as a complete CDS")
    else:
        cs, ce = longest_orf(s)
        if cs < 0:
            cs, ce = 0, len(s)
            notes.append("no complete ORF found; treating the whole sequence as CDS")
        else:
            notes.append(f"ORF detected in the input sequence: {cs + 1}..{ce} "
                         f"(1-based)")
    return Transcript(tx_id=name, gene_id=name, mrna=s,
                      cds_start=cs, cds_end=ce, source="direct", notes=notes)


# ============================================================
# Isoform selection
# ============================================================

def pick_longest(txs: List[Transcript], by: str = "cds") -> Transcript:
    """Default isoform rule: longest CDS, tie-broken by longest mRNA, then ID."""
    if by == "mrna":
        key = lambda t: (len(t.mrna), t.cds_len, -_id_rank(t.tx_id))
    else:
        key = lambda t: (t.cds_len, len(t.mrna), -_id_rank(t.tx_id))
    return sorted(txs, key=key, reverse=True)[0]


def _id_rank(tx_id: str) -> int:
    m = re.search(r"\.(\d+)$", tx_id)
    return int(m.group(1)) if m else 999


def load_transcripts(query: str, *, genome_fa: str = "", gff3: str = "",
                     cdna_fa: str = "", cds_fa: str = "",
                     sequence: str = "") -> Tuple[List[Transcript], str]:
    """Unified entry point. Returns (transcripts, mode)."""
    if sequence:
        return [from_sequence(sequence, query or "user_seq")], "direct"

    if gff3 and genome_fa:
        ann = Annotation()
        ann.parse(gff3, gene_filter={query.split('.')[0]})
        tx_ids = ann.transcripts_of(query)
        if not tx_ids:
            ann2 = Annotation()
            ann2.parse(gff3)
            ann = ann2
            tx_ids = ann.transcripts_of(query)
        if tx_ids:
            chroms = {ann.tx_info[t]["chrom"] for t in tx_ids
                      if t in ann.tx_info}
            genome = Genome(genome_fa, only=chroms or None)
            txs = [t for t in (ann.build(genome, tid) for tid in tx_ids) if t]
            if txs:
                return txs, "gff3"

    if cdna_fa:
        txs = from_cdna_fasta(cdna_fa, query)
        if txs:
            return txs, "cdna"

    if cds_fa:
        txs = from_cds_fasta(cds_fa, query)
        if txs:
            return txs, "cds"

    return [], "none"
