"""
GSprimer — terminal interaction.

The selection step is MANDATORY by design: GSprimer never silently picks a
primer pair for the user. When driven by an AI agent the same decision point
is served by `--select`, with the agent asking the user instead.
"""

import sys
from typing import List

from .pipeline import candidate_table


def _input(prompt: str, default: str = "") -> str:
    try:
        v = input(prompt).strip()
        return v or default
    except (EOFError, KeyboardInterrupt):
        print("\n[cancelled]", file=sys.stderr)
        sys.exit(0)


def ask_selection(pairs: List, max_pick: int = 3) -> List[int]:
    """Prompt the user to choose candidate pair(s). Returns 1-based indices."""
    print("\n" + "=" * 72)
    print("Select the primer pair(s) to attach adapters and finalize "
          "(this step is your decision)")
    print("=" * 72)
    print("  ★ = strict ATG..STOP standard anchor pair")
    print("  Tier A passes all SpacerFinder-derived filters; B = boundaries; "
          "C = needs trade-off")
    print("  The Risk column takes the worst of the five risk blocks\n")

    anchor = next((i for i, p in enumerate(pairs, 1) if p.is_anchor), None)
    hint = f" (Enter = {anchor}, the standard anchor pair)" if anchor else ""

    while True:
        raw = _input(f"  Enter number(s), comma-separated, up to {max_pick} "
                     f"{hint}: ", str(anchor) if anchor else "")
        if not raw:
            print("  [!] please select at least one")
            continue
        try:
            sel = [int(x) for x in raw.replace(" ", "").split(",") if x]
        except ValueError:
            print("  [!] enter only digits and commas")
            continue
        bad = [i for i in sel if i < 1 or i > len(pairs)]
        if bad:
            print(f"  [!] number(s) {bad} out of range 1-{len(pairs)}")
            continue
        if len(sel) > max_pick:
            print(f"  [!] at most {max_pick} selections")
            continue

        for i in sel:
            p = pairs[i - 1]
            r = p.risk or {}
            if r.get("overall") == "HIGH":
                reasons = []
                for k, cn in [("enzyme", "internal type-IIS sites"),
                              ("gc", "extreme GC region"),
                              ("specificity", "specificity"),
                              ("frame", "reading frame / termini"),
                              ("amplification", "amplification difficulty")]:
                    if (r.get(k) or {}).get("level") == "HIGH":
                        reasons.append(cn)
                print(f"  [!] candidate {i} is HIGH risk: "
                      f"{'; '.join(reasons)}")
        return sel


def show_table(pairs: List) -> str:
    s = candidate_table(pairs)
    print(s)
    return s
