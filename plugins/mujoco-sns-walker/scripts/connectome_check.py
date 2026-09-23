"""Validate a connectome spec JSON against the runner's REAL consumption
contract (runner.py ~lines 858-881) BEFORE a run silently applies 0 rules.

Catches the known failure classes:
  - block-editor export (nodes/edges/synapses) saved as connectome_gains.json
    -> parses fine, applies 0 rules, only trace is a startup print  [ERROR]
  - gain_key not in params.G -> rule silently skipped               [ERROR]
  - two enabled rules sharing one gain_key -> last wins             [WARN]
  - disabled rules are SKIPPED, not zeroed (runner `continue`s)     [INFO]
  - any enabled rule with hops>=1 sets the single global full_rules [INFO]
  - sign of `gain` is discarded (runner takes abs)                  [INFO]

Usage:
  connectome_check.py                 # checks spinal/connectome_gains.json
  connectome_check.py --file X.json
Exit code: 0 ok (warnings allowed), 1 errors.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import wcommon  # noqa: E402

KNOWN_RULE_IDS = {"ia_homo", "ia_recip", "ii_exc", "ii_inh", "ib_auto", "ib_rev",
                  "heel", "toe", "ib_load", "rc", "cross"}


def params_g_keys(spinal: Path) -> tuple[set[str], str]:
    """Gain-key allowlist from params.G. Try importing (exact); fall back to
    a regex scan of the G block (works under any python)."""
    try:
        sys.path.insert(0, str(spinal))
        for m in [m for m in list(sys.modules) if m == "params"]:
            del sys.modules[m]
        import params  # noqa: E402
        return set(params.G.keys()), "imported params.G"
    except Exception:
        src = (spinal / "params.py").read_text(encoding="utf-8", errors="replace")
        g0 = src.find("G = {")
        g1 = src.find("\n}", g0)
        block = src[g0:g1] if g0 >= 0 and g1 > g0 else src
        keys = set(re.findall(r"^\s*\"?([A-Za-z_][A-Za-z0-9_]*)\"?\s*[:=]\s*[-0-9.]", block, re.M))
        return keys, "regex-scanned params.py G block"


def main() -> int:
    wcommon.utf8_stdio()
    ap = argparse.ArgumentParser()
    ap.add_argument("--file")
    ap.add_argument("--spinal")
    args = ap.parse_args()

    spinal = Path(args.spinal) if args.spinal else wcommon.resolve_spinal()
    path = Path(args.file) if args.file else spinal / "connectome_gains.json"
    wcommon.banner(f"checking {path}")

    if not path.exists():
        print(f"FILE ABSENT: {path}")
        print("  (no connectome spec applied - runner uses params defaults as-is; that is valid)")
        return 0
    try:
        spec = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        print(f"ERROR: JSON parse failed: {e}")
        return 1

    errors, warnings = [], []

    top = set(spec.keys())
    if "nodes" in top and ("edges" in top or "synapses" in top):
        errors.append(
            "This is a BLOCK-EDITOR export (nodes/edges/synapses). runner.py only reads "
            "top-level 'rules' - it would parse this file and silently apply 0 rules. "
            "Re-export from connectome_editor.html (the RULES editor), or keep this "
            "file under a different name as a design document.")
    if "rules" not in spec:
        errors.append("no top-level 'rules' key - runner applies nothing from this file")
        _report(spec, errors, warnings)
        return 1

    g_keys, how = params_g_keys(spinal)
    print(f"params.G allowlist: {len(g_keys)} keys ({how})")

    rules = spec["rules"]
    if not isinstance(rules, dict):
        errors.append("'rules' must be an object {rule_id: {...}}")
        _report(spec, errors, warnings)
        return 1

    gain_to_rules: dict[str, list[str]] = {}
    n_enabled = 0
    full_rules_effect = False
    applied: list[str] = []

    for rid, rr in rules.items():
        if not isinstance(rr, dict):
            errors.append(f"rule '{rid}': value must be an object")
            continue
        if rid not in KNOWN_RULE_IDS:
            warnings.append(f"rule '{rid}': not one of the 11 editor rule ids "
                            f"({sorted(KNOWN_RULE_IDS)}) - runner does NOT validate ids, "
                            "so this only applies if gain_key is real")
        enabled = rr.get("enabled", False)
        gk = rr.get("gain_key")
        gain = rr.get("gain", 0.0)
        hops = rr.get("hops", 0)
        if not isinstance(enabled, bool):
            errors.append(f"rule '{rid}': 'enabled' must be true/false")
        if gk is None or not isinstance(gk, str):
            errors.append(f"rule '{rid}': missing string 'gain_key'")
        if not isinstance(gain, (int, float)):
            errors.append(f"rule '{rid}': 'gain' must be numeric")
        if not isinstance(hops, int):
            warnings.append(f"rule '{rid}': 'hops' should be an int (0=direct, >=1=interneuron)")
        if enabled and isinstance(gk, str):
            n_enabled += 1
            if gk not in g_keys:
                errors.append(f"rule '{rid}': gain_key '{gk}' NOT in params.G "
                              "-> runner silently skips this rule")
            else:
                applied.append(f"{rid}: G['{gk}'] = {abs(float(gain)):g}")
                gain_to_rules.setdefault(gk, []).append(rid)
            if isinstance(hops, int) and hops >= 1:
                full_rules_effect = True
        if enabled is False:
            warnings.append(f"rule '{rid}': disabled -> runner SKIPS it (does not zero the gain; "
                            "the params default stays)")

    for gk, rids in gain_to_rules.items():
        if len(rids) > 1:
            warnings.append(f"gain_key '{gk}' set by multiple enabled rules {rids} - "
                            "dict order decides (last wins)")

    print(f"runner would apply {n_enabled} rule(s):")
    for a in applied or ["  (none enabled)"]:
        print(f"  {a}")
    if full_rules_effect:
        print("full_rules -> 1.0 (any enabled rule with hops>=1; single global toggle, "
              "no per-rule mixing)")

    _report(spec, errors, warnings)
    return 1 if errors else 0


def _report(spec, errors, warnings) -> None:
    unknown = set(spec.keys()) - {"rules", "_comment"}
    if unknown:
        warnings.append(f"unknown top-level keys ignored by runner: {sorted(unknown)}")
    for w in warnings:
        print(f"WARN: {w}")
    for e in errors:
        print(f"ERROR: {e}")
    print("VERDICT: " + ("FAIL - fix before running" if errors else "OK (warnings above, if any)"))


if __name__ == "__main__":
    raise SystemExit(main())
