"""Sync ProofFinal .tex from the canonical Overleaf zip (Bolen_Dissertation.zip of 2026-10-02),
adapting chapter graphicspath declarations from figs/<X>/ to ../Figures/<X>/.
Local-only files in ProofFinal are left untouched. Writes a .bak? No - git is the backup (clean tree at 89b9cbef)."""
import re
from pathlib import Path

ZIP = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\overleaf_20261002")
PF = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal")

files = [ZIP / "main.tex"] + sorted((ZIP / "chapters").glob("*.tex"))
for src in files:
    rel = src.relative_to(ZIP)
    dst = PF / rel
    text = src.read_text(encoding="utf-8")
    lines = []
    for ln in text.splitlines():
        if "\\graphicspath" in ln:
            ln = ln.replace("figs/", "../Figures/")
        lines.append(ln)
    dst.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print("synced", rel)
