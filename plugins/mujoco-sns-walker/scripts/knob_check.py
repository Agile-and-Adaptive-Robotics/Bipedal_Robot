"""Four-place wiring audit for a new params.G knob (the recurring trap: a knob
added in one place but missing in the others silently misbehaves - the v7
renshaw omission cost a silent 0.14 kine delta).

Checks, READ-ONLY (wiring edits belong to Ben):
  1. params.py      default present in the G dict
  2. runner.py      consumed by the --best JSON loader (the JSON RULE tuple
                    or a dedicated `if key in best` branch)
  3. _curriculum.py appears in the stage KEYS tuple
  4. _curriculum.py has a suggest_* call AND a seed-dict entry (missing seed
                    keys get SAMPLED - the stage-2 lesson)

Usage:  knob_check.py contra_swing pm_add NEW_KNOB
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import wcommon  # noqa: E402


def find_lines(src: str, pattern: str, ctx: int = 0) -> list[tuple[int, str]]:
    out = []
    lines = src.splitlines()
    for i, line in enumerate(lines):
        if re.search(pattern, line):
            for j in range(max(0, i - ctx), min(len(lines), i + ctx + 1)):
                out.append((j + 1, lines[j].strip()[:130]))
            if ctx == 0:
                continue
    return out


def check_knob(spinal: Path, knob: str) -> bool:
    print(f"== knob '{knob}' ==")
    ok = True

    # 1. params.py G default
    src_p = (spinal / "params.py").read_text(encoding="utf-8", errors="replace")
    hits = find_lines(src_p, rf"""["']?{re.escape(knob)}["']?\s*[:=]\s*[-0-9.]""")
    print(f"  params.py default        : {'PRESENT' if hits else 'MISSING'} {hits[:1]}")
    ok &= bool(hits)

    # 2. runner.py loader consumption
    src_r = (spinal / "runner.py").read_text(encoding="utf-8", errors="replace")
    in_tuple = find_lines(src_r, rf"""["']{re.escape(knob)}["']\s*,""")
    in_branch = find_lines(src_r, rf"""["']{re.escape(knob)}["']\s+in\s+best""")
    in_assign = find_lines(src_r, rf"""G\[["']{re.escape(knob)}["']\]""")
    loader_ok = bool(in_tuple or in_branch)
    print(f"  runner.py --best loader  : {'PRESENT' if loader_ok else 'MISSING'}"
          f"{' (tuple)' if in_tuple else ''}{' (branch)' if in_branch else ''} "
          f"{(in_tuple or in_branch or in_assign)[:1]}")
    ok &= loader_ok

    # 3. curriculum KEYS
    src_c = (spinal / "_curriculum.py").read_text(encoding="utf-8", errors="replace")
    key_hits = find_lines(src_c, rf"""["']{re.escape(knob)}["']""")
    in_keys = any("KEYS" in h or "STAGE_KEYS" in h for _, h in key_hits) or bool(key_hits)
    print(f"  _curriculum KEYS         : {'PRESENT' if in_keys else 'MISSING'} {key_hits[:2]}")
    ok &= bool(key_hits)

    # 4. suggest + seed
    sug = find_lines(src_c, rf"""suggest_(float|int|categorical)\([^)]*["']{re.escape(knob)}["']""")
    print(f"  _curriculum suggest_*    : {'PRESENT' if sug else 'MISSING'} {sug[:1]}")
    ok &= bool(sug)
    # seed dict: dict-literal key OR seed["knob"] = ... assignment
    seed_pat = (r"""^\s*["']?""" + re.escape(knob) + r"""["']?\s*:"""
                r"""|\[['"]""" + re.escape(knob) + r"""['"]\]\s*=""")
    seed = find_lines(src_c, seed_pat)
    print(f"  _curriculum seed dict    : {'PRESENT' if seed else 'MISSING (optuna will SAMPLE it)'} {seed[:1]}")

    if not ok:
        print("  -> at least one place MISSING: fix wiring BEFORE the next study "
              "(edits to wiring are Ben's call)")
    return ok


def main() -> int:
    wcommon.utf8_stdio()
    knobs = sys.argv[1:]
    if not knobs:
        raise SystemExit(__doc__)
    spinal = wcommon.resolve_spinal()
    wcommon.banner(f"spinal={spinal}")
    results = {k: check_knob(spinal, k) for k in knobs}
    print("== summary ==")
    for k, v in results.items():
        print(f"  {k:24s} {'WIRED' if v else 'INCOMPLETE'}")
    return 0 if all(results.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
