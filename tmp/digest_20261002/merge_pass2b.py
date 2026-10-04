"""Merge pass 2 continuation (methods already applied): results, discussion, appendixB, appendixC."""
from pathlib import Path

PF = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal")
REP = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\o808_export")
EE2 = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\ee2_tex")

def read(p):
    return p.read_text(encoding="utf8")

def write(p, t):
    p.write_text(t, encoding="utf8", newline="\n")

def lines_of(p, a, b):
    return "\n".join(read(p).splitlines()[a - 1:b])

def spiking_block(rel):
    ls = read(EE2 / rel).splitlines()
    i0 = next(i for i, l in enumerate(ls) if "Spiking-mirror campaign BEGIN" in l)
    i1 = next(i for i, l in enumerate(ls) if "Spiking-mirror campaign END" in l)
    return "\n".join(ls[i0:i1 + 1])

# ---------- 30-results ----------
r = read(PF / "chapters/30-results.tex")
new_vals = lines_of(REP / "chapters/30-results.tex", 362, 363)
new_redesign = lines_of(REP / "chapters/30-results.tex", 365, 366)
old_vals_start = "The correction terms were re-identified over the full archived test sets"
old_redesign_start = "The route redesign built on these values produced usable flexor and extensor model definitions."
ls = r.splitlines()
oi = next(i for i, l in enumerate(ls) if l.startswith(old_vals_start))
oj = next(i for i, l in enumerate(ls) if l.startswith(old_redesign_start))
assert oi < oj, (oi, oj)
assert "8.9" in ls[oi], "old values paragraph unexpected"
r = "\n".join(ls[:oi] + [new_vals, "", new_redesign] + ls[oj + 1:]) + "\n"
assert "+3.94" in r and "0.158" in r and "flxr77" in r, "new values missing"
assert r.count("\\section{Follow-Up Identification and Route Redesign}") == 1, "section duplicated"
assert r.count("\\label{sec:ongoing}") == 1, "label duplicated"
r = r.rstrip() + "\n\n" + spiking_block("30-results.tex") + "\n"
write(PF / "chapters/30-results.tex", r)
print("results OK")

# ---------- 40-discussion ----------
d = read(PF / "chapters/40-discussion.tex")
rl = read(REP / "chapters/40-discussion.tex").splitlines()
ki = next(i for i, l in enumerate(rl) if "The neural layer is no longer only a proposal" in l)
repo_synth = rl[ki]
dl = d.splitlines()
zi = next(i for i, l in enumerate(dl) if l.startswith("This chain was developed deliberately so that each stage can feed the next."))
assert "What remains missing is the neural layer itself" in dl[zi]
dl[zi] = repo_synth
d = "\n".join(dl) + "\n"
d = d.rstrip() + "\n\n" + spiking_block("40-discussion.tex") + "\n"
write(PF / "chapters/40-discussion.tex", d)
print("discussion OK")

# ---------- 93-AppendixB ----------
b = read(PF / "chapters/93-AppendixB.tex")
if "Spiking-mirror campaign BEGIN" not in b:
    b = b.rstrip() + "\n\n" + spiking_block("93-AppendixB.tex") + "\n"
    write(PF / "chapters/93-AppendixB.tex", b)
print("appendixB OK")

# ---------- 94-AppendixC ----------
c = read(REP / "chapters/94-AppendixC.tex")
cl = []
for ln in c.splitlines():
    if "\\graphicspath" in ln:
        ln = ln.replace("figs/", "../Figures/")
    cl.append(ln)
write(PF / "chapters/94-AppendixC.tex", "\n".join(cl) + "\n")
print("appendixC OK")
