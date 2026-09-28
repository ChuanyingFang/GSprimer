"""
GSprimer — Configuration management.

Paths can be set via config file, environment variables, or programmatically.
Search order:  env var  >  ~/.gsprimer/config.yaml  >  ~/.spacerfinder/config.yaml  >  defaults

GSprimer intentionally falls back to the SpacerFinder config so that a lab that
already runs SpacerFinder does not need to duplicate genome/BLAST paths.
"""

import os
import sys
from pathlib import Path
from typing import Optional

_DEFAULTS = {
    # --- sequence sources ---
    "genome": "",       # genome FASTA (for GFF3 mode + BLAST context)
    "gff3": "",         # GFF3 annotation
    "cdna": "",         # mRNA/cDNA multi-FASTA (fast mode)
    "cds": "",          # CDS multi-FASTA (fast mode)
    # --- specificity ---
    "blast_db": "",     # nucleotide BLAST DB (genome)
    "blast_db_cdna": "",  # optional transcriptome DB
    "blastn": "blastn",  # assume in PATH
    # --- GS vector adapters (overridable for other vector series) ---
    "adapter_f": "gtgatatcAGGTCTCTcgag",
    "adapter_r": "gccgcgggTGGTCTCAatcc",
    "enzyme": "BsaI",
    "enzyme_site": "GGTCTC",
}


def _find_config() -> Optional[Path]:
    """Find config file in standard locations."""
    candidates = [
        Path.home() / ".gsprimer" / "config.yaml",
        Path.home() / ".config" / "gsprimer" / "config.yaml",
        Path("gsprimer.yaml"),
        Path.home() / ".spacerfinder" / "config.yaml",  # inherit SpacerFinder
    ]
    if env := os.environ.get("GSPRIMER_CONFIG"):
        candidates.insert(0, Path(env))
    for p in candidates:
        if p.exists():
            return p
    return None


def _load_yaml(path: Path) -> dict:
    """Load YAML config (pure-Python fallback if PyYAML is absent)."""
    try:
        import yaml
        with open(path, encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    except ImportError:
        config = {}
        with open(path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if ":" in line:
                    k, v = line.split(":", 1)
                    v = v.split("#")[0] if v.strip().startswith(("/", "~")) else v
                    config[k.strip()] = v.strip().strip('"').strip("'")
        return config
    except Exception:
        return {}


def _get_config() -> dict:
    cfg = dict(_DEFAULTS)
    cfg_path = _find_config()
    if cfg_path:
        try:
            loaded = _load_yaml(cfg_path)
            # only adopt keys we know about
            for k, v in loaded.items():
                if k in _DEFAULTS and v:
                    cfg[k] = v
            cfg["_config_file"] = str(cfg_path)
        except Exception:
            pass

    env_map = {
        "GSPRIMER_GENOME": "genome",
        "GSPRIMER_GFF3": "gff3",
        "GSPRIMER_CDNA": "cdna",
        "GSPRIMER_CDS": "cds",
        "GSPRIMER_BLAST_DB": "blast_db",
        "GSPRIMER_BLASTN": "blastn",
        # SpacerFinder compatibility
        "SPACERFINDER_GENOME": "genome",
        "SPACERFINDER_GFF3": "gff3",
        "SPACERFINDER_BLAST_DB": "blast_db",
        "SPACERFINDER_BLASTN": "blastn",
    }
    for env_var, key in env_map.items():
        if val := os.environ.get(env_var):
            cfg[key] = val
    return cfg


_CONFIG = None


def get(key: str) -> str:
    global _CONFIG
    if _CONFIG is None:
        _CONFIG = _get_config()
    return _CONFIG.get(key, _DEFAULTS.get(key, ""))


def set_config(key: str, value: str) -> None:
    global _CONFIG
    if _CONFIG is None:
        _CONFIG = _get_config()
    if value:
        _CONFIG[key] = value


def all_config() -> dict:
    global _CONFIG
    if _CONFIG is None:
        _CONFIG = _get_config()
    return dict(_CONFIG)


def require(key: str) -> str:
    val = get(key)
    if not val:
        print(f"[GSprimer] Error: required config '{key}' is not set.",
              file=sys.stderr)
        print("  Set it in ~/.gsprimer/config.yaml, via --{} CLI flag, "
              "or an environment variable.".format(key.replace('_', '-')),
              file=sys.stderr)
        sys.exit(1)
    return val
