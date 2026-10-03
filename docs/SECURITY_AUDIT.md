# 🔍 Security Audit Report — GSprimer

## 📊 Executive summary

- **Audited object**: GSprimer (source `C:\Users\Administrator\Desktop\GSprimer-github.zip`, upstream `github.com/ChuanyingFang/GSprimer`)
- **Install path**: `~/.workbuddy/skills/GSprimer/`
- **Audit date**: 2026-10-03
- **Audit method**: pure static text analysis (read-only, no execution interaction)
- **Total findings**: 0
  - 🔴 P0 blocking: 0
  - ⚠️ P1 noteworthy: 0
  - 📝 Informational notes: 2 (non-risk items)
- **Security score**: 98 / 100

---

## 🔴 P0 blocking findings

✅ No P0 risk found

---

## ⚠️ P1 noteworthy findings

✅ No P1 risk found

---

## 📝 Informational notes (non-risk items)

1. **Dependency-install instructions do not mention a virtual environment** (informational)
   - **Location**: `README.md:38`, `.github/workflows/ci.yml:23`
   - **Snippet**: `pip install -e .`
   - **Explanation**: Appears in docs and the CI workflow as **manual user actions**; the skill
     never runs the install automatically. CI installs run on an isolated GitHub Actions runner.
     Under the audit boundary ("examples only, must be manually copied by the user"), this is not
     a poisoning risk.
   - **Suggestion**: already changed in `SKILL.md` to recommend a venv or direct
     `python -m gsprimer.cli`, avoiding pollution of the global environment.

2. **Version number mismatch** (informational)
   - **Location**: `pyproject.toml:7` = `1.1.0` vs `gsprimer/__init__.py:15` = `1.2.0`
   - **Explanation**: code-consistency issue, not a security risk.
   - **Suggestion**: unify when releasing upstream; `__init__.__version__` is authoritative.

---

## 📋 Detailed check results

### Command execution & privilege checks

- Occurrences: 1 (legitimate use)
- Details:
  - `gsprimer/blast.py:65` — `subprocess.run([blastn, "-task", "blastn-short", ...], capture_output=True, timeout=120)`
    - ✅ arguments passed as a **list** (no shell-injection surface), **no `shell=True`**,
      **120 s timeout**, temp files cleaned up in `finally`
    - Purpose: optional specificity check (local `blastn`), fully disableable with `--no-blast`
- Not found: `eval` / `exec` / `os.system` / `os.popen` / `pickle` / `__import__` / `curl|bash` / `nohup` / `disown`

### File operations & sensitive-path checks

- Occurrences: 1 deletion call (safe)
- Details:
  - `gsprimer/blast.py:102` — `os.unlink(p)`, inside a `finally` block, where `p` is a
    `gsp_*.fa` / `gsp_*.tsv` temp file created by `tempfile.mkstemp()`
  - All write targets are user-specified `--outdir` (report HTML/TSV, `gsprimer_state.json`, `order.tsv`)
  - All read targets are user-explicitly-passed GFF3 / FASTA / config paths
- Not found: `~/.ssh`, `~/.gnupg`, `/etc/passwd`, `~/.aws`, `.env`, `credentials`,
  `secrets`, `api_key`, `private_key`, `shutil.rmtree`, `rm -rf` (except CI cleanup of `tests/out_ci`)
- Sensitive-path hits: none

### Network-request checks

- URLs found (all project links in docs/metadata, no executable network calls):
  - `https://github.com/ChuanyingFang/GSprimer` (`pyproject.toml` URL metadata ×3, `README.md:36` clone instructions)
  - `https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/` (`README.md:42`, BLAST+ download instructions)
- Network calls in code: **0** — no `requests` / `urllib` / `httpx` / `aiohttp` / `socket` / `fetch` / `curl` / `wget`
- Base64 detection: no suspicious encoded strings
- Report HTML (`report.py`) contains no external CDN links / `<script src=...>`; entirely inline SVG and styles

### Remote-script deep analysis

Not applicable — this skill never auto-downloads or executes any remote script (no `curl | bash`, no
download-then-`eval`), so step 4 is not triggered.

### Dependency-install risk checks

- **Global-install detection**: no skill-auto-executed global install commands found
- **Virtual-environment check**: `SKILL.md` already provides a venv-isolated option; the core
  pipeline needs no install (`python -m gsprimer.cli` works directly)
- **Dependency-source check**: `requirements.txt` declares **zero required core dependencies**; the
  only optional one is `pyyaml` (falls back to a pure-Python parser when absent, `config.py:57-67`).
  No `--index-url`, no `git+https://` install

---

## 💡 Overall recommendations

1. This skill is a local-computation tool (primer design + report generation) with a minimal attack
   surface; safe to use directly.
2. The only optional external process is `blastn`; if your lab environment does not install BLAST,
   always pass `--no-blast` or `--fallback-fasta` explicitly, to avoid the report silently marking
   specificity as `UNKNOWN` and it being misread as safe.
3. Recommended to unify the version numbers in `pyproject.toml` and `__init__.py` upstream.

---

## ✅ Audit conclusion

**Risk level**: **P2 — Safe (no poisoning risk)**

**Usage recommendation**:
- ✅ **P2 - safe to use**: no network requests, no hidden execution, no sensitive-path access, no
  credential exfiltration, no automatic global install; the only subprocess call is a parameterized,
  timeout-guarded local `blastn`.

---

**📌 Audit-scope reminder**: reports only the skill's own supply-chain poisoning risk; does not
report teaching-code quality issues.
