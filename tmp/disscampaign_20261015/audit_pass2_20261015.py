# -*- coding: utf-8 -*-
"""Pass 2: comment-stripped include/graphics check, orphan recount, PDF page count,
log summary, external-draft old-path check. Read-only."""
import os, re
from pathlib import Path

BASE = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation")
PF = BASE / "ProofFinal"
FIGS = BASE / "Figures"

def norm(s):
    return s.replace("\r\n", "\n").replace("\r", "\n")

def strip_comments(text):
    out = []
    for line in text.split("\n"):
        l = line.replace("\\%", "\x00")
        i = l.find("%")
        if i >= 0:
            l = l[:i]
        out.append(l.replace("\x00", "\\%"))
    return "\n".join(out)

main_tex = norm((PF / "main.tex").read_text(encoding="utf-8", errors="replace"))
main_code = strip_comments(main_tex)

# ---- include order, comment-stripped ----
includes = re.findall(r"\\include\{([^}]*)\}", main_code)
print("INCLUDE ORDER (comment-stripped):")
for i, inc in enumerate(includes, 1):
    print("  %2d. %s" % (i, inc))

def gp_dirs_of(text):
    dirs = []
    for m in re.finditer(r"\\graphicspath\{((?:\s*\{[^{}]*\}\s*)+)\}", text):
        for d in re.findall(r"\{([^{}]*)\}", m.group(1)):
            dirs.append(d)
    return dirs

EXTS = [".pdf", ".png", ".jpg", ".jpeg", ".eps"]
inc_re = re.compile(r"\\includegraphics\s*(?:\[[^\]]*\])?\s*\{([^}]*)\}")
current = gp_dirs_of(main_code)
resolved, unresolved_g, commented_only = [], [], []
for inc in includes:
    fp = PF / (inc + ".tex")
    t = norm(fp.read_text(encoding="utf-8", errors="replace"))
    tc = strip_comments(t)
    own = gp_dirs_of(tc)
    if own:
        current = own
    for m in inc_re.finditer(tc):
        name = m.group(1).strip()
        root, ext = os.path.splitext(name)
        bases = [name + e for e in EXTS] if not ext else (
            [name, root + "-eps-converted-to.pdf"] if ext.lower() == ".eps" else [name])
        existing = []
        for b in bases:
            if (PF / b).exists():
                existing.append(("cwd", str(PF / b)))
            for d in current:
                if (PF / d / b).exists():
                    existing.append(("gp:" + d, str(PF / d / b)))
        if existing:
            resolved.append((inc, name, existing))
        else:
            unresolved_g.append((inc, name))
    # commented-out includegraphics that do NOT resolve (informational)
    for m in inc_re.finditer(t):
        if m.group(1).strip() in [u[1] for u in unresolved_g]:
            pass
print("\nACTIVE (comment-stripped) \\includegraphics: %d resolved, %d UNRESOLVED" %
      (len(resolved), len(unresolved_g)))
for inc, name in unresolved_g:
    print("   UNRESOLVED: %s {%s}" % (inc, name))
print("active graphicspath at each chapter with figures:")
for inc in includes:
    t = strip_comments(norm((PF / (inc + ".tex")).read_text(encoding="utf-8", errors="replace")))
    n = len(inc_re.findall(t))
    if n:
        print("   %-28s %2d figs  gp=%s" % (inc, n, gp_dirs_of(t)))

# ---- orphan recount on comment-stripped references ----
referenced = set()
for inc, name, existing in resolved:
    for _, abspath in existing:
        referenced.add(Path(abspath).resolve())
extra = set()
for r in referenced:
    if r.suffix.lower() == ".eps":
        c = r.with_name(r.stem + "-eps-converted-to.pdf")
        if c.exists():
            extra.add(c.resolve())
referenced |= extra
allfigs = sorted(p for p in FIGS.rglob("*") if p.is_file())
orphans = [p for p in allfigs if p.resolve() not in referenced]
act = [p for p in orphans if "Deprecated" not in p.parts]
dep = [p for p in orphans if "Deprecated" in p.parts]
print("\nORPHANS (active-tree, comment-stripped refs): %d; Deprecated: %d; total files: %d"
      % (len(act), len(dep), len(allfigs)))
byfolder = {}
for p in act:
    byfolder.setdefault(str(p.parent.relative_to(FIGS)), []).append(p.name)
for k in sorted(byfolder):
    print("  %-55s %2d files" % (k, len(byfolder[k])))

# ---- PDF page count (naive) vs log ----
pdf = PF / "build_wf" / "main.pdf"
data = pdf.read_bytes()
n_page_objs = len(re.findall(rb"/Type\s*/Page[^s]", data))
counts = re.findall(rb"/Count\s+(\d+)", data)
print("\nmain.pdf: %d bytes; naive /Type/Page objects=%d; max /Count=%s" %
      (len(data), n_page_objs, max(int(c) for c in counts) if counts else None))

# ---- main.log scan ----
log = (PF / "build_wf" / "main.log").read_text(encoding="utf-8", errors="replace")
lines = log.split("\n")
print("main.log: %d lines" % len(lines))
err = [l for l in lines if l.startswith("!")]
print("lines starting with '!': %d %s" % (len(err), err[:10]))
und = [l for l in lines if "ndefined" in l]
print("'undefined/Undefined' lines: %d" % len(und))
for l in und[:20]:
    print("   ", l.strip()[:140])
out = [l for l in lines if "Output written" in l]
print("'Output written':", out)
rerun = [l for l in lines if "Rerun" in l]
print("'Rerun' lines: %d" % len(rerun))
for l in rerun[:10]:
    print("   ", l.strip()[:140])
missing = [l for l in lines if "not found" in l.lower() or "File `" in l and "not found" in l]
print("'not found' lines: %d" % len(missing))
for l in missing[:15]:
    print("   ", l.strip()[:140])

# ---- external drafts referencing OLD figure paths ----
print("\nOLD-PATH references (Figures/20-methods|30-results|40-discussion|50-futurework|60-conclusion) in Dissertation .tex outside ProofFinal:")
for t in sorted(BASE.rglob("*.tex")):
    if PF in t.parents:
        continue
    txt = norm(t.read_text(encoding="utf-8", errors="replace"))
    for m in re.finditer(r"Figures/(?:20-methods|30-results|40-discussion|50-futurework|60-conclusion)/[^\s{}'\"]*", txt):
        print("   %s: %s" % (t.relative_to(BASE).as_posix(), m.group(0)))
