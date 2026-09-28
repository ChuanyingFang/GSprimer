"""
GSprimer — CDS amplification primer design for GS (Golden Gate / BsaI) vectors.

Primer thermodynamics and all quality filters are inherited from
SpacerFinder v2.4.5; the only intentional deviation is that amplicon size is
unconstrained (it is fixed by the ORF).

Quick start
-----------
    from gsprimer import run_design, run_finalize
    res = run_design("LOC_Os01g01010", cds_fa="MSU_cds.fa", outdir="out")
    fin = run_finalize([1], outdir="out")
"""

__version__ = "1.2.0"

from .adapters import ADAPTER_F, ADAPTER_R, order_sheet, parse_adapter
from .design import PrimerPair, design
from .pipeline import candidate_table, run_design, run_finalize
from .risk import assess, enzyme_risk, gc_risk, scan_sites
from .template import (Transcript, load_transcripts, pick_longest,
                       read_fasta)

__all__ = [
    "__version__",
    "run_design", "run_finalize", "candidate_table",
    "design", "PrimerPair",
    "Transcript", "load_transcripts", "pick_longest", "read_fasta",
    "assess", "enzyme_risk", "gc_risk", "scan_sites",
    "order_sheet", "parse_adapter", "ADAPTER_F", "ADAPTER_R",
]
