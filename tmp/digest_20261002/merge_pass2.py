"""Merge pass 2026-10-02 (evening, take 2):
- Fix Xi picks in 30-results (flexor front row 77 / extensor flxr77 row 32) using the repo's
  audited paragraphs (o808 export), replacing the superseded pick-1 paragraphs.
- Port repo-only Methods content: GoF decomposition paragraphs + 5%-margin sentence + audit comments.
- Port repo Discussion synthesis paragraph.
- Replace 94-AppendixC with the repo's audited version (paths adapted).
- Append the easteregg2 spiking-mirror fenced blocks to Methods/Results/Discussion/AppendixB.
Every replacement asserts its anchor; nothing is silently skipped."""
from pathlib import Path

PF = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal")
REP = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\o808_export")
EE2 = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\ee2_tex")

def read(p):
    return p.read_text(encoding="utf8")

def write(p, t):
    p.write_text(t, encoding="utf8", newline="\n")

def must_replace(text, old, new, what):
    assert old in text, f"ANCHOR MISSING ({what})"
    assert text.count(old) == 1, f"ANCHOR NOT UNIQUE ({what}): {text.count(old)}"
    return text.replace(old, new)

def lines_of(p, a, b):
    """1-based inclusive line range."""
    return "\n".join(read(p).splitlines()[a - 1:b])

def spiking_block(rel):
    t = read(EE2 / rel)
    ls = t.splitlines()
    i0 = next(i for i, l in enumerate(ls) if "Spiking-mirror campaign BEGIN" in l)
    i1 = next(i for i, l in enumerate(ls) if "Spiking-mirror campaign END" in l)
    return "\n".join(ls[i0:i1 + 1])

# ---------- 20-methods ----------
m = read(PF / "chapters/20-methods.tex")
# (a) GoF decomposition paragraphs before the optimization section
gof = lines_of(REP / "chapters/20-methods.tex", 370, 372)
anchor = "\\section{Optimization Algorithms and Cost Functions}"
assert anchor in m
m = m.replace(anchor, gof + "\n\n" + anchor, 1)
# (b) margin sentence: swap zip 1%/raw sentence for repo 5%-both + audit comment
rep_margin_para = lines_of(REP / "chapters/20-methods.tex", 389, 390)
i = rep_margin_para.find("The flexor target also carries")
assert i > 0
rep_sentence = rep_margin_para[i:].strip()
audit = lines_of(REP / "chapters/20-methods.tex", 390, 390)
old_sentence = ("The flexor target also carries a required torque margin of 1\\% above the human curve, "
                "whereas the extensor target is the raw human curve.")
assert old_sentence in m, "zip margin sentence not found"
m = must_replace(m, old_sentence, rep_sentence + "\n" + audit, "methods margin sentence")
# (c) spiking methods block
m = m.rstrip() + "\n\n" + spiking_block("20-methods.tex") + "\n"
write(PF / "chapters/20-methods.tex", m)
print("methods: GoF paras + margin + spiking block  OK")

# ---------- 30-results ----------
r = read(PF / "chapters/30-results.tex")
new_vals = lines_of(REP / "chapters/30-results.tex", 362, 363)
new_redesign = lines_of(REP / "chapters/30-results.tex", 365, 366)
old_vals_start = "The correction terms were re-identified over the full archived test sets"
old_redesign_start = "The route redesign built on these values produced usable flexor and extensor model definitions."
ls = r.splitlines()
oi = next(i for i, l in enumerate(ls) if l.startswith(old_vals_start))
oj = next(i for i, l in enumerate(ls) if l.startswith(old_redesign_start))
assert oi < oj
# old block = values para line oi .. redesign para line oj (inclusive), each is one long line
r = "\n".join(ls[:oi] + [new_vals, "", new_redesign] + ls[oj + 1:]) + "\n"
r = must_replace(r, "$\\chi_{0} = \\qty{+8.9}{\\mm}$", "$\\chi_{0} = \\qty{+3.94}{\\mm}$", "value sanity check")
r = must_replace(r, "with $\\chi_{3} = 0.621$", "with $\\chi_{3} = 0.158$", "ext chi3 sanity")
r = r.rstrip() + "\n\n" + spiking_block("30-results.tex") + "\n"
write(PF / "chapters/30-results.tex", r)
print("results: pick-77/row-32 paragraphs + spiking block  OK")

# ---------- 40-discussion ----------
d = read(PF / "chapters/40-discussion.tex")
rl = read(REP / "chapters/40-discussion.tex").splitlines()
ki = next(i for i, l in enumerate(rl) if "The neural layer is no longer only a proposal" in l)
repo_synth = rl[ki]
old_synth_start = "This chain was developed deliberately so that each stage can feed the next."
dl = d.splitlines()
zi = next(i for i, l in enumerate(dl) if l.startswith(old_synth_start))
# zip paragraph contains the old 'What remains missing' clause; replace that one line wholesale
assert "What remains missing is the neural layer itself" in dl[zi]
dl[zi] = repo_synth
d = "\n".join(dl) + "\n"
# my cross-engine paragraph follows; keep it. Append spiking discussion section.
d = d.rstrip() + "\n\n" + spiking_block("40-discussion.tex") + "\n"
write(PF / "chapters/40-discussion.tex", d)
print("discussion: synthesis swap + spiking block  OK")

# ---------- 93-AppendixB ----------
b = read(PF / "chapters/93-AppendixB.tex")
b = b.rstrip() + "\n\n" + spiking_block("93-AppendixB.tex") + "\n"
write(PF / "chapters/93-AppendixB.tex", b)
print("appendixB: spiking block  OK")

# ---------- 94-AppendixC: repo version wholesale, path-adapted ----------
c = read(REP / "chapters/94-AppendixC.tex")
cl = []
for ln in c.splitlines():
    if "\\graphicspath" in ln:
        ln = ln.replace("figs/", "../Figures/")
    cl.append(ln)
write(PF / "chapters/94-AppendixC.tex", "\n".join(cl) + "\n")
print("appendixC: repo audited version installed  OK")
