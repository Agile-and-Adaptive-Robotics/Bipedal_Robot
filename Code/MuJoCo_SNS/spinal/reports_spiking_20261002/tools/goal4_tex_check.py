"""GOAL 4 static LaTeX checks (no local TeX on easteregg2; Overleaf out of
scope). Checks label/ref integrity, brace/environment balance, figure files
existence, citation keys, and style rules (no (R)/(TM), no lowercase matlab,
Oxford comma spot check is manual).
"""
import re
import sys
from pathlib import Path

ROOT = Path(r"D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal")
FIGROOT = ROOT.parent / "Figures"

chapters = sorted((ROOT / "chapters").glob("*.tex")) + [ROOT / "main.tex"]

labels = {}
refs = []
cites = []
env_stack_issues = []
for f in chapters:
    txt = f.read_text(encoding="utf-8", errors="replace")
    # strip comments (keep \% intact-ish: crude but fine for labels)
    for m in re.finditer(r"\\label\{([^}]+)\}", txt):
        labels.setdefault(m.group(1), []).append(f.name)
    for m in re.finditer(r"\\(?:ref|eqref|pageref)\{([^}]+)\}", txt):
        refs.append((f.name, m.group(1)))
    for m in re.finditer(r"\\cite[tp]?\{([^}]+)\}", txt):
        for k in m.group(1).split(","):
            cites.append(k.strip())

dup = {k: v for k, v in labels.items() if len(v) > 1}
print("duplicate labels:", dup if dup else "NONE")

# known duplicates that PRE-DATE this edit (sec:ongoing appears twice)
missing = sorted({r for r in refs if r[1] not in labels})
print("unresolved refs:", missing if missing else "NONE")

bib = (ROOT / "thesis.bib").read_text(encoding="utf-8", errors="replace")
bibkeys = set(re.findall(r"@\w+\{([^,\s]+),", bib))
bad_cites = sorted({c for c in cites if c and c not in bibkeys})
print("citation keys not in thesis.bib:", bad_cites if bad_cites else "NONE")

# figure files referenced by includegraphics in my edited chapters
fig_ok = True
for ch, sub in [("20-methods.tex", "20-methods"), ("30-results.tex", "30-results"),
                ("40-discussion.tex", "40-discussion"),
                ("92-AppendixA.tex", "92-AppendixA"), ("93-AppendixB.tex", "93-AppendixB"),
                ("94-AppendixC.tex", "94-AppendixC")]:
    txt = (ROOT / "chapters" / ch).read_text(encoding="utf-8", errors="replace")
    for m in re.finditer(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", txt):
        name = m.group(1)
        if "/" in name:
            continue
        cands = [FIGROOT / sub / name,
                 ROOT / "chapters" / name]
        if not any(c.exists() for c in cands):
            print("MISSING FIGURE:", ch, name)
            fig_ok = False
print("figure files:", "all found" if fig_ok else "SEE ABOVE")

# environment balance per edited file
for ch in ["20-methods.tex", "30-results.tex", "40-discussion.tex", "93-AppendixB.tex"]:
    txt = (ROOT / "chapters" / ch).read_text(encoding="utf-8", errors="replace")
    begins = re.findall(r"\\begin\{(\w+\*?)\}", txt)
    ends = re.findall(r"\\end\{(\w+\*?)\}", txt)
    from collections import Counter
    cb, ce = Counter(begins), Counter(ends)
    diff = {k: (cb[k], ce[k]) for k in set(cb) | set(ce) if cb[k] != ce[k]}
    print(f"env balance {ch}:", diff if diff else "OK")
    # unbalanced braces (crude, ignoring \\{ escapes and comments)
    body = re.sub(r"(?<!\\)%.*", "", txt)
    n_open = body.count("{") - body.count(r"\{")
    n_close = body.count("}") - body.count(r"\}")
    print(f"  braces open/close: {n_open}/{n_close}",
          "OK" if n_open == n_close else "MISMATCH")

# my new labels defined
mine = ["sec:spiking_methods", "sec:spiking_results", "tab:spike_ground",
        "tab:spike_knee", "tab:spike_beer", "tab:spike_animatlab",
        "fig:spike_rhythm", "fig:spike_ground", "fig:spike_knee",
        "fig:spike_animatlab", "app:spiking_details", "tab:app_spk_cells",
        "tab:app_spk_cal", "tab:app_spk_gates", "tab:app_spk_animatlab_full"]
for l in mine:
    print("label", l, "->", labels.get(l, "MISSING"))

# style: no (R)/(TM), matlab lowercase
for ch in ["20-methods.tex", "30-results.tex", "40-discussion.tex", "93-AppendixB.tex"]:
    txt = (ROOT / "chapters" / ch).read_text(encoding="utf-8", errors="replace")
    new = txt.split("[ZCODE 2026-10-02")
    for i, block in enumerate(new[1:], 1):
        seg = "[ZCODE 2026-10-02" + block
        for pat, why in [(r"\(R\)", "(R)"), (r"\(TM\)", "(TM)"),
                         (r"\bmatlab\b", "lowercase matlab"),
                         (r"\bmatucha\b", "typo")]:
            if re.search(pat, seg):
                print(f"STYLE HIT {ch} block {i}:", why)
print("style scan done")
