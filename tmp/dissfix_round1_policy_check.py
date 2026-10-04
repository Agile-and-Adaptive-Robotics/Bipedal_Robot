"""Policy checks per the standing dissertation rules:

1. Every figure environment must have an ALT TEXT comment above its \\caption.
2. Figures/<chapter> top level: only image files referenced by the dissertation;
   pieces belong in Components/, superseded in Deprecated/.
3. Report em-dash characters in chapters (candidate new-prose violations).
"""
import os
import re

BASE = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
PF = os.path.join(BASE, "ProofFinal")
FIG = os.path.join(BASE, "Figures")
EXTS = {".pdf", ".png", ".jpg", ".jpeg", ".eps"}

# ---- collect includegraphics per chapter + graphicspath dirs (ref-gate logic)
used = {}  # chapter num -> set of basenames referenced
for texf in sorted(os.listdir(os.path.join(PF, "chapters"))):
    if not texf.endswith(".tex"):
        continue
    text = open(os.path.join(PF, "chapters", texf), encoding="utf8", errors="replace").read()
    lines = [ln for ln in text.splitlines() if not ln.lstrip().startswith("%")]
    body = "\n".join(lines)
    ch = texf.split("-")[0]
    s = used.setdefault(ch, set())
    for m in re.finditer(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", body):
        s.add(os.path.basename(m.group(1)))

# ---- 1. ALT TEXT check
print("=== ALT TEXT check (figures missing an ALT TEXT comment above caption) ===")
alt_missing = 0
for texf in sorted(os.listdir(os.path.join(PF, "chapters"))):
    if not texf.endswith(".tex"):
        continue
    path = os.path.join(PF, "chapters", texf)
    lines = open(path, encoding="utf8", errors="replace").read().splitlines()
    in_fig = 0
    for i, ln in enumerate(lines):
        s = ln.lstrip()
        if s.startswith("%"):
            continue
        if re.match(r"\\begin\{figure\*?\}", s):
            in_fig += 1
            continue
        if re.match(r"\\end\{figure\*?\}", s):
            in_fig = max(0, in_fig - 1)
            continue
        if in_fig and re.match(r"\\caption", s):
            # look upward inside the figure for an ALT TEXT comment line
            window = lines[max(0, i - 30):i]
            # stop at figure begin
            for w in reversed(window):
                ws = w.strip()
                if re.match(r"\\(begin|end)\{figure", ws):
                    break
                if ws.startswith("%") and "ALT TEXT" in ws:
                    break
            else:
                print(f"  {texf}:{i+1}: caption without ALT TEXT comment")
                alt_missing += 1
print(f"  total missing: {alt_missing}")

# ---- 2. top-level figure policy
print("=== Top-level figure-folder policy (image files at top level, referenced?) ===")
for d in sorted(os.listdir(FIG)):
    dp = os.path.join(FIG, d)
    if not os.path.isdir(dp):
        continue
    ch = d.split("-")[0]
    ref = used.get(ch, set())
    for f in sorted(os.listdir(dp)):
        fp = os.path.join(dp, f)
        if os.path.isfile(fp) and os.path.splitext(f)[1].lower() in EXTS:
            base = os.path.splitext(f)[0]
            hit = any(r.startswith(base) or base == os.path.splitext(r)[0] for r in ref)
            if not hit:
                print(f"  Figures/{d}/{f}: NOT referenced by chapter {d}")

# ---- 3. em-dash unicode scan
print("=== Unicode em-dash (\u2014) occurrences in chapters ===")
for texf in sorted(os.listdir(os.path.join(PF, "chapters"))):
    if not texf.endswith(".tex"):
        continue
    path = os.path.join(PF, "chapters", texf)
    for i, ln in enumerate(open(path, encoding="utf8", errors="replace").read().splitlines(), 1):
        if "\u2014" in ln:
            print(f"  {texf}:{i}: {ln.strip()[:100]}")
