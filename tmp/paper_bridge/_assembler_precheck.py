# Assembler pre-check: cite keys in fragments vs bib; graphics targets vs disk.
# Read-only. Run with D:/Anaconda/python.exe.
import re, os, sys

ROOT = r"D:/GitHub/Bipedal_Robot"
SECDIR = os.path.join(ROOT, "tmp/paper_bridge/sections")
FRAGS = ["10_abstract_intro.tex", "20_methods.tex", "30_results.tex",
         "40_discussion_conclusion.tex", "tables_main.tex"]
FIGDIR = os.path.join(ROOT, "Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures")
BIB = os.path.join(ROOT, "Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/BolenFrontiers22.bib")

def strip_comments(t):
    return "\n".join(l for l in t.splitlines() if not l.lstrip().startswith("%"))

bib = open(BIB, encoding="utf-8", errors="replace").read()
bibkeys = set(re.findall(r"@\w+\{([^,\s]+)\s*,", bib))

allkeys = {}
for f in FRAGS:
    t = strip_comments(open(os.path.join(SECDIR, f), encoding="utf-8").read())
    for m in re.finditer(r"\\(?:cite|citet|citep|citealp|citeauthor)\*?(?:\[[^\]]*\])*\{([^}]*)\}", t):
        for k in m.group(1).split(","):
            k = k.strip()
            if k:
                allkeys.setdefault(k, []).append(f)

missing = {k: v for k, v in allkeys.items() if k not in bibkeys}
print(f"UNIQUE CITE KEYS IN FRAGMENTS: {len(allkeys)}")
print("MISSING FROM BIB:", missing if missing else "NONE")

# graphics targets
targets = []
for f in FRAGS:
    t = strip_comments(open(os.path.join(SECDIR, f), encoding="utf-8").read())
    for m in re.finditer(r"\\includegraphics\[[^\]]*\]\{([^}]*)\}", t):
        targets.append((f, m.group(1)))
print("\nGRAPHICS TARGETS:")
ok = True
for f, g in targets:
    p = os.path.join(os.path.dirname(FIGDIR), g)  # g is Figures/x relative to ProofFinal
    exists = os.path.exists(p)
    ok = ok and exists
    print(f"  {f}: {g} -> {'EXISTS' if exists else 'MISSING'}")
print("ALL GRAPHICS EXIST:", ok)

# labels defined in fragments + refs used (to catch collisions/missing within fragments+tables)
labels, refs = [], []
for f in FRAGS:
    t = strip_comments(open(os.path.join(SECDIR, f), encoding="utf-8").read())
    labels += [(f, m.group(1)) for m in re.finditer(r"\\label\{([^}]*)\}", t)]
    refs += [(f, m.group(1)) for m in re.finditer(r"\\(?:ref|eqref)\{([^}]*)\}", t)]
dupes = {}
for f, l in labels:
    dupes.setdefault(l, []).append(f)
d = {l: v for l, v in dupes.items() if len(v) > 1}
print("\nDUPLICATE LABELS ACROSS FRAGMENTS:", d if d else "NONE")
labelset = set(dupes)
unresolved = [(f, r) for f, r in refs if r not in labelset]
print("REFS NOT DEFINED IN FRAGMENTS (need Main/supp):", unresolved if unresolved else "NONE")
print("\nLABELS:", sorted(labelset))
