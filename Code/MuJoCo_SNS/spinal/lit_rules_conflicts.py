"""Literature rules conflict audit (Ben 2026-09-30): cross-compare the
four rule sets (rybak / shevtsova / shinohara replication JSONs + Ben's
own 09-24 rules) for SIGN CONFLICTS on homologous connections and
gain disagreements > 2x. Output: campaigns/20260930/lit_rules_conflict_audit.md

Homology: labels normalized (lowercase, side suffix stripped, family
synonyms mapped: IniF~IN-F, ia_recip~IaIN, ib_inh~IBIN, etc.).
"""
import json
import os
import re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = {
    "rybak": "replication/rybak_rules.json",
    "shevtsova": "replication/shevtsova_rules.json",
    "shinohara": "replication/shinohara_rules.json",
    "ben": "ben_rules_20260924.json",
}

SYN = {
    "in-f": "in-f", "inf": "in-f", "inflex": "in-f",
    "in-e": "in-e", "ine": "in-e",
    "ia_recip": "iain", "ia_a": "iain", "iain": "iain",
    "ib_inh": "ibin", "ibin": "ibin", "ib_a": "ibin",
    "ii_inh": "iiin", "iiin": "iiin",
    "renshaw": "rc", "rc": "rc",
    "mn": "mn", "motor": "mn",
    "rg-e": "rg-e", "rge": "rg-e", "extensor hc": "rg-e",
    "rg-f": "rg-f", "rgf": "rg-f", "flexor hc": "rg-f",
}


def norm(label):
    s = label.lower().strip()
    s = re.sub(r"[_\s](l|r)$", "", s)          # side suffix
    s = re.sub(r"[^a-z0-9 ]", " ", s).strip()
    # family extraction: keep tokens that name a cell class
    toks = s.split()
    for t in toks:
        if t in SYN:
            return SYN[t]
    if "rg" in s and ("ext" in s or "e" == s[-1]):
        return "rg-e"
    if "rg" in s and "fl" in s:
        return "rg-f"
    if "mn" in s or "motor" in s:
        # agonist/antagonist distinction matters for Ia circuits
        return "mn-ant" if "antagon" in s else "mn"
    if "ia" in toks:
        return "ia"
    if "ii" in toks:
        return "ii"
    if "ib" in toks:
        return "ib"
    if "vest" in s:
        return "vest"
    if "v0d" in s:
        return "v0d"
    if "v0v" in s:
        return "v0v"
    if "v2a" in s:
        return "v2a"
    if "v3" in s:
        return "v3"
    if "v1" in s:
        return "v1"
    if "v2b" in s:
        return "v2b"
    if "di6" in s or "dI6" in s.lower():
        return "di6"
    if "mlr" in s or "brainstem" in s or "supraspinal" in s:
        return "supra"
    if "drive" in s or "descend" in s:
        return "supra"
    return s


def load_edges():
    all_edges = defaultdict(list)     # (src, dst) -> [(file, sign, gain, tag)]
    for name, fn in FILES.items():
        d = json.load(open(os.path.join(HERE, fn), encoding="utf-8"))
        labels = {n["id"]: n["label"] for n in d["nodes"]}
        for e in d.get("edges") or d.get("synapses") or []:
            src = norm(labels.get(e["from"], e["from"]))
            dst = norm(labels.get(e["to"], e["to"]))
            sign = e.get("sign", "?")
            gain = float(e.get("gain", 0) or 0)
            all_edges[(src, dst)].append(
                (name, sign, gain, e.get("tag", "")))
    return all_edges


def main():
    edges = load_edges()
    sign_conflicts, gain_conflicts, multi = [], [], []
    for (src, dst), recs in sorted(edges.items()):
        if len(recs) < 2:
            continue
        files = {r[0] for r in recs}
        signs = {r[1] for r in recs}
        if len(signs) > 1:
            sign_conflicts.append((src, dst, recs))
        elif len(files) > 1:
            gains = [r[2] for r in recs if r[2] > 0]
            if gains and max(gains) > 2 * max(min(gains), 1e-9) and \
                    min(gains) > 0:
                gain_conflicts.append((src, dst, recs))
            else:
                multi.append((src, dst, recs))
    out = ["# Literature rules conflict audit (2026-09-30)", "",
           "Cross-comparison of rybak_rules.json, shevtsova_rules.json,",
           "shinohara_rules.json vs ben_rules_20260924.json (the four",
           "connectome rule sets in spinal\\). Normalized homolog labels;",
           f"{len(edges)} distinct (src,dst) families across the sets.", ""]
    out += [f"## SIGN CONFLICTS ({len(sign_conflicts)})", "",
            "| src | dst | per-file |", "|---|---|---|"]
    for src, dst, recs in sign_conflicts:
        cell = "; ".join(f"{f}:{s}(g={g}) [{t}]" for f, s, g, t in recs)
        out.append(f"| {src} | {dst} | {cell} |")
    out += ["", f"## GAIN DISAGREEMENTS >2x, same sign ({len(gain_conflicts)})",
            "", "| src | dst | per-file |", "|---|---|---|"]
    for src, dst, recs in gain_conflicts:
        cell = "; ".join(f"{f}:{s}(g={g})" for f, s, g, t in recs)
        out.append(f"| {src} | {dst} | {cell} |")
    out += ["", f"## Agreements (same sign, compatible gains; {len(multi)})",
            ""]
    for src, dst, recs in multi:
        cell = "; ".join(f"{f}:{s}(g={g})" for f, s, g, t in recs)
        out.append(f"- {src} -> {dst}: {cell}")
    p = os.path.join(HERE, "campaigns", "20260930",
                     "lit_rules_conflict_audit.md")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w", encoding="utf-8").write("\n".join(out) + "\n")
    print("\n".join(out[:60]))
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
