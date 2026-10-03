"""Verify ChatGPT's claim: local ProofFinal (stale) chapter image references all resolve under
local Figures tree. Counts references and reports misses."""
import re
from pathlib import Path

PF = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal")
FIG = PF.parent / "Figures"
EXTS = [".pdf", ".png", ".jpg", ".jpeg", ".eps", ".PNG", ".JPG", ".PDF"]
KNOWN = {e.lower() for e in EXTS}

count, miss = 0, []
for tex in sorted((PF / "chapters").glob("*.tex")):
    text = "\n".join(ln for ln in tex.read_text(encoding="utf-8", errors="replace").splitlines()
                     if not ln.lstrip().startswith("%"))
    dirs = []
    for g in re.findall(r"\\graphicspath\{((?:\{[^}]*\})+)\}", text):
        dirs += re.findall(r"\{([^}]*)\}", g)
    dirs = [d.replace("../Figures/", "") for d in dirs]
    for m in re.finditer(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", text):
        name = m.group(1)
        count += 1
        cands = [name] if Path(name).suffix.lower() in KNOWN else [name + e for e in EXTS]
        ok = any((FIG / d / c).exists() for d in dirs for c in cands) or any(
            (FIG / c).exists() for c in cands)
        if not ok:
            miss.append((tex.name, name, dirs))
print(f"stale ProofFinal includegraphics refs: {count}; missing: {len(miss)}")
for x in miss:
    print("  MISS", x)
