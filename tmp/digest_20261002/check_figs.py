"""Cross-check every \\includegraphics target in the Overleaf zip against (a) the zip's own
figs tree and (b) the local Dissertation\\Figures tree (with figs/ prefix stripped)."""
import re
from pathlib import Path

ZIP = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\overleaf_20261002")
FIG = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures")
EXTS = [".pdf", ".png", ".jpg", ".jpeg", ".eps", ".PNG", ".JPG", ".PDF"]
KNOWN = {e.lower() for e in EXTS}

chapter_gp = {}
all_targets = []
for tex in sorted((ZIP / "chapters").glob("*.tex")) + [ZIP / "main.tex"]:
    text = tex.read_text(encoding="utf-8", errors="replace")
    # strip comment-only lines so commented-out figures are excluded
    live = "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("%"))
    dirs = []
    for g in re.findall(r"\\graphicspath\{((?:\{[^}]*\})+)\}", live):
        dirs += re.findall(r"\{([^}]*)\}", g)
    chapter_gp[tex.name] = dirs
    for m in re.finditer(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", live):
        all_targets.append((tex.name, m.group(1)))

def find(root: Path, name: str, dirs):
    cands = [name] if Path(name).suffix.lower() in KNOWN else [name + e for e in EXTS]
    for d in dirs:
        for c in cands:
            if (root / d / c).exists():
                return f"{d}/{c}"
    hits = sorted({str(p.relative_to(root)).replace("\\", "/") for c in cands for p in root.rglob(c)})
    if len(hits) == 1:
        return hits[0] + "  (project-wide)"
    return f"AMBIG {hits}" if hits else None

print(f"{'src':22s} {'target':40s} {'zip':46s} local")
problems = []
for src, tgt in all_targets:
    zd = chapter_gp.get(src, [])
    z = find(ZIP, tgt, zd + ["."])
    l = find(FIG, tgt, [d.replace("figs/", "") for d in zd] + ["."])
    if not (z and l):
        problems.append((src, tgt, z, l))
    print(f"{src:22s} {tgt:40s} {(z or 'MISSING-IN-ZIP'):46s} {l or 'MISSING-LOCALLY'}")

print(f"\nPROBLEMS: {len(problems)}")
for p in problems:
    print("  ", p)
