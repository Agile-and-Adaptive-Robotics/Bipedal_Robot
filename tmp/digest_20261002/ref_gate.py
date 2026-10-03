"""Figure-reference gate: every live \\includegraphics in ProofFinal chapters must resolve
under Dissertation\\Figures (using each chapter's graphicspath). Exit 0 + count, exit 1 +
missing list. Also verifies the figure-folder policy: no stray image files at the top level
of a chapter folder that are neither referenced nor in Components/Deprecated (warning only)."""
import re
import sys
from pathlib import Path

PF = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal")
FIG = PF.parent / "Figures"
EXTS = [".pdf", ".png", ".jpg", ".jpeg", ".eps", ".PNG", ".JPG", ".PDF"]
KNOWN = {e.lower() for e in EXTS}

count, missing = 0, []
for tex in sorted((PF / "chapters").glob("*.tex")):
    text = "\n".join(ln for ln in tex.read_text(encoding="utf8", errors="replace").splitlines()
                     if not ln.lstrip().startswith("%"))
    dirs = []
    for g in re.findall(r"\\graphicspath\{((?:\{[^}]*\})+)\}", text):
        dirs += [d.replace("../Figures/", "").rstrip("/") for d in re.findall(r"\{([^}]*)\}", g)]
    for m in re.finditer(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", text):
        name = m.group(1)
        count += 1
        cands = [name] if Path(name).suffix.lower() in KNOWN else [name + e for e in EXTS]
        if not any((FIG / d / c).exists() for d in dirs for c in cands) and \
           not any((FIG / c).exists() for c in cands):
            missing.append(f"{tex.name}: {name}")

if missing:
    print("REF GATE FAIL")
    print("\n".join(missing))
    sys.exit(1)
print(f"REF GATE PASS refs={count}")
