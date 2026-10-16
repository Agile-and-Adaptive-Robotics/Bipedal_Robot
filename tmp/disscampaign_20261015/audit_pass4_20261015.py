# -*- coding: utf-8 -*-
"""Pass 4: page count from log; pre-vs-post includegraphics set comparison;
CPG_spinal_section_draft figure reference style."""
import re
from pathlib import Path

BASE = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation")
PF = BASE / "ProofFinal"
OL = Path(r"D:\Github\Bipedal_Robot\tmp\disscampaign_20261015\overleaf_pull")

log = (PF / "build_wf" / "main.log").read_text(encoding="utf-8", errors="replace")
m = re.search(r"Output written on .*?\((\d+) pages?,\s*(\d+) bytes?\)", log)
print("Output written: pages=%s bytes=%s (regex on full log text)" %
      (m.group(1) if m else "NOT-PARSED", m.group(2) if m else "NOT-PARSED"))

def strip_comments(text):
    out = []
    for line in text.replace("\r\n", "\n").split("\n"):
        l = line.replace("\\%", "\x00")
        i = l.find("%")
        if i >= 0:
            l = l[:i]
        out.append(l.replace("\x00", "\\%"))
    return "\n".join(out)

inc_re = re.compile(r"\\includegraphics\s*(?:\[[^\]]*\])?\s*\{([^}]*)\}")

def active_names(files):
    s = []
    for t in files:
        tc = strip_comments(t.read_text(encoding="utf-8", errors="replace"))
        for mm in inc_re.finditer(tc):
            s.append((t.name, mm.group(1).strip()))
    return s

old = active_names(sorted(OL.rglob("*.tex")))
new = active_names(sorted(PF.rglob("*.tex")))
print("\npre-restructure (Overleaf pull) active includegraphics: %d" % len(old))
print("current (ProofFinal) active includegraphics: %d" % len(new))

def basekey(arg):
    """Normalize: strip figs/ or ../Figures/ prefixes and subfolders to basename."""
    return Path(arg).name

old_by_key = {}
for f, a in old:
    old_by_key.setdefault(basekey(a), []).append((f, a))
new_keys = set(basekey(a) for _, a in new)
print("\nOLD-referenced basenames missing from the NEW active set (candidates for 'no longer referenced'):")
missing = 0
for k in sorted(old_by_key):
    if k not in new_keys:
        missing += 1
        print("   %s  (was: %s)" % (k, old_by_key[k]))
print("count:", missing)
print("\nNEW-referenced basenames not in OLD set (newly referenced):")
old_keys = set(old_by_key)
for k in sorted(new_keys):
    if k not in old_keys:
        print("   %s" % k)

# CPG draft figure style
draft = BASE / "Notes" / "CPG_spinal_section_draft.tex"
if draft.exists():
    txt = draft.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
    print("\nCPG_spinal_section_draft.tex: %d includegraphics; first 8 args:" %
          len(inc_re.findall(strip_comments(txt))))
    for mm in list(inc_re.finditer(txt))[:8]:
        print("   {%s}" % mm.group(1))
    gp = re.findall(r"\\graphicspath\{((?:\s*\{[^{}]*\}\s*)+)\}", txt)
    print("   graphicspath:", gp)
else:
    print("\nCPG_spinal_section_draft.tex NOT FOUND at", draft)
