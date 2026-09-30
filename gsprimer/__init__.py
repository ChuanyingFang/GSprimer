"""
GSprimer — CDS / promoter amplification primer design for GS (Golden Gate /
BsaI / PaqCI) vectors.

Primer thermodynamics and all quality filters are inherited from
SpacerFinder v2.4.5; the only intentional deviation is that amplicon size is
unconstrained for CDS (it is fixed by the ORF), while promoter mode targets a
1200–1700 bp window upstream of the TSS.

Quick start
-----------
    from gsprimer import run_design, run_finalize
    # CDS mode
    res = run_design("LOC_Os01g01010", cds_fa="MSU_cds.fa", outdir="out")
    # Promoter mode (two-file: 2 kb upstream + transcript 5'-UTR fallback)
    res = run_design("GeneX", promoter_seq=upstream_2k, utr5_seq=utr5,
                     outdir="out")
    fin = run_finalize([1], outdir="out")
"""

__version__ = "2.0.0"

from .adapters import ADAPTER_F, ADAPTER_R, order_sheet, parse_adapter
from .design import PrimerPair, design, design_promoter
from .pipeline import candidate_table, run_design, run_finalize
from .risk import assess, enzyme_risk, gc_risk, scan_sites
from .template import (Transcript, load_transcripts, pick_longest,
                       read_fasta, from_promoter, utr5_from_transcript)

__all__ = [
    "__version__",
    "run_design", "run_finalize", "candidate_table",
    "design", "design_promoter", "PrimerPair",
    "Transcript", "load_transcripts", "pick_longest", "read_fasta",
    "from_promoter", "utr5_from_transcript",
    "assess", "enzyme_risk", "gc_risk", "scan_sites",
    "order_sheet", "parse_adapter", "ADAPTER_F", "ADAPTER_R",
]
