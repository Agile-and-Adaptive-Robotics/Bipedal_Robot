"""Static diagnosis of ProofFinal chapters: env balance + includegraphics existence.

Scratch diagnostic for dissertation repair round 1. Writes findings to stdout.
"""
import os
import re
import sys

BASE = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
CH = os.path.join(BASE, "ProofFinal", "chapters")
FIGROOT = os.path.join(BASE, "Figures")

chapters = sorted(f for f in os.listdir(CH) if f.endswith(".tex"))

env_re = re.compile(r"\\(begin|end)\{([^}]*)\}")
inc_re = re.compile(r"\\includegraphics\*?(?:\[[^\]]*\])?\{([^}]*)\}")

fail = False

for name in chapters:
    path = os.path.join(CH, name)
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    # strip comments (unescaped %) to avoid false positives
    stripped = re.sub(r"(?<!\\)%.*", "", text)
    begins, ends = {}, {}
    for kind, env in env_re.findall(stripped):
        if kind == "begin":
            begins[env] = begins.get(env, 0) + 1
        else:
            ends[env] = ends.get(env, 0) + 1
    probs = []
    for env in sorted(set(begins) | set(ends)):
        b, e = begins.get(env, 0), ends.get(env, 0)
        if b != e:
            probs.append(f"env '{env}' begin={b} end={e}")
    if probs:
        fail = True
        print(f"[ENV] {name}: " + "; ".join(probs))
    # includegraphics existence check (graphicspath = ../Figures/, chapter subdirs
    # declared per-chapter via \graphicspath additions -- collect those too)
    for grp in inc_re.findall(stripped):
        rel = grp.replace("\\", "/")
        if not os.path.splitext(rel)[1]:
            # no extension: try common ones
            cands = [rel + ext for ext in (".pdf", ".png", ".jpg", ".eps")]
        else:
            cands = [rel]
        hits = []
        for root, dirs, files in os.walk(FIGROOT):
            for c in cands:
                cand = os.path.join(root, *c.split("/"))
                if os.path.isfile(cand):
                    hits.append(os.path.relpath(cand, FIGROOT))
        if not hits:
            fail = True
            print(f"[MISSING-FIG] {name}: {grp} -> no file under Figures/")
    # per-chapter graphicspath additions
    for gp in re.findall(r"\\graphicspath\{([^}]*)\}", stripped):
        inner = re.findall(r"\{([^}]*)\}", gp)
        for d in inner:
            if not os.path.isdir(os.path.join(BASE, "ProofFinal", d.replace("../", "ProofFinal_XX"))) and not os.path.isdir(os.path.normpath(os.path.join(BASE, "ProofFinal", d))):
                # graphicspath entries are relative to the compile cwd (ProofFinal)
                p = os.path.normpath(os.path.join(BASE, "ProofFinal", d))
                if not os.path.isdir(p):
                    fail = True
                    print(f"[GRAPHICSPATH] {name}: entry '{d}' -> not a dir: {p}")

print("STATIC CHECK DONE fail=", fail)
