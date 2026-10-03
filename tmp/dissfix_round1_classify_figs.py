"""Classify top-level image files in Figures/<chapter>: referenced / live-epstopdf-byproduct / unreferenced.

A file is REFERENCED if its basename (any extension form) appears in any chapter's
\includegraphics. A `<name>-eps-converted-to.pdf` whose source `<name>` is referenced
is a LIVE BYPRODUCT (epstopdf regenerates it) -> keep. Everything else at top level
that no chapter references -> candidate for Components/ or Deprecated/.
"""
import os
import re
from pathlib import Path

BASE = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation")
PF = BASE / "ProofFinal" / "chapters"
FIG = BASE / "Figures"
EXTS = {".pdf", ".png", ".jpg", ".jpeg", ".eps"}

names = set()          # raw includegraphics arguments
for texf in sorted(PF.glob("*.tex")):
    body = "\n".join(ln for ln in texf.read_text(encoding="utf8", errors="replace").splitlines()
                     if not ln.lstrip().startswith("%"))
    names |= set(re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", body))

basenames = set()
for n in names:
    basenames.add(n)
    basenames.add(os.path.basename(n))
    basenames.add(os.path.splitext(os.path.basename(n))[0])

print(f"total distinct includegraphics names: {len(names)}")
print()
for d in sorted(p for p in FIG.iterdir() if p.is_dir()):
    for f in sorted(d.iterdir()):
        if not f.is_file() or f.suffix.lower() not in EXTS:
            continue
        stem = f.stem  # filename without extension
        if stem in basenames or f.name in basenames:
            status = "REFERENCED"
        elif stem.endswith("-eps-converted-to") and stem[: -len("-eps-converted-to")] in basenames:
            status = "LIVE-BYPRODUCT (source referenced; epstopdf regenerates)"
        else:
            src = None
            if stem.endswith("-eps-converted-to"):
                base_src = stem[: -len("-eps-converted-to")]
                hits = [p.name for p in d.glob(base_src + ".*")]
                src = f"source files present: {hits}" if hits else "source .eps ABSENT"
            status = f"UNREFERENCED {src or ''}"
        if status != "REFERENCED":
            print(f"{d.name}\\{f.name}  ->  {status}")
