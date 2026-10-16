# -*- coding: utf-8 -*-
"""Post-restructure verification (read-only): conservation + label checks done right."""
import os, re
from collections import Counter

ROOT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
CH = os.path.join(ROOT, "ProofFinal", "chapters")
LOGDIR = os.path.join(ROOT, "Notes", "revision_campaign_20261015")

# donors are deleted; their content is recoverable for checking from the new files
# themselves + the recorded slices. Instead re-check: (a) prior_walker lines present in
# 15-background; (b) all labels across the full chapter set are unique and complete vs
# the pre-restructure label inventory recorded in the phase-A log; (c) include order.

def rd(p):
    with open(p, "r", encoding="utf-8", newline="") as f:
        return f.read()

lab = re.compile(r"\\label\{([^}]+)\}")
all_labels = Counter()
per_file = {}
for name in sorted(os.listdir(CH)):
    if not name.endswith(".tex"):
        continue
    txt = rd(os.path.join(CH, name))
    ls = lab.findall(txt)
    per_file[name] = ls
    all_labels.update(ls)

dups = [l for l, c in all_labels.items() if c > 1]
print("total labels:", sum(all_labels.values()), "unique:", len(all_labels))
print("DUPLICATE labels:", dups if dups else "none")

for want in ["ch:methods", "ch:results", "ch:discussion", "ch:futurework", "ch:conclusion",
             "ch:muscles", "ch:torque", "ch:ankle", "ch:database", "ch:sns", "ch:overall",
             "ch:introduction", "ch:background",
             "sec:prior_walker", "sec:methods_overview", "sec:foundation", "sec:synthesis",
             "sec:force_model", "sec:results_force", "sec:results_torque", "sec:ongoing",
             "sec:balance_testbeds", "sec:results_balance", "sec:preliminary_sim_methods",
             "sec:preliminary_sim_results", "sec:spiking_methods", "sec:spiking_results",
             "sec:disc_error", "sec:disc_pulley", "sec:disc_spiking", "sec:optimization"]:
    where = [n for n, ls in per_file.items() if want in ls]
    print("%-32s -> %s" % (want, where if where else "*** MISSING ***"))

# every \ref target in every chapter resolves to a defined label
ref = re.compile(r"\\(?:ref|eqref|pageref)\{([^}]+)\}")
unresolved = set()
for name in sorted(os.listdir(CH)):
    if not name.endswith(".tex"):
        continue
    for t in ref.findall(rd(os.path.join(CH, name))):
        if t not in all_labels:
            unresolved.add((name, t))
main_txt = rd(os.path.join(ROOT, "ProofFinal", "main.tex"))
print("\nUNRESOLVED \\ref targets in chapters+main:", sorted(unresolved) if unresolved else "none")

# prior_walker lines present in 15-background verbatim
bg = rd(os.path.join(CH, "15-background.tex")).splitlines()
pw = ["\\section{Prior Neuromechanical Simulation: The AnimatLab Bipedal Walker}\\label{sec:prior_walker}"]
hits = [i for i, ln in enumerate(bg, 1) if pw[0] in ln]
print("\nprior_walker heading in 15-background at line(s):", hits)
for i in hits:
    for j in range(i, i + 6):
        print("  bg:%d %s" % (j, bg[j-1][:100]))

# include order
inc = re.findall(r"\\include\{(chapters/[^}]+)\}", main_txt)
print("\nmain.tex include order:")
for i in inc:
    print("  ", i)
