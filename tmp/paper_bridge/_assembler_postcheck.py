# Post-assembly verification of the merged Main.tex.
# 1. Environment balance. 2. Cite keys vs bib. 3. Graphics targets vs disk.
# 4. Refs vs labels (Main + supplement). 5. Fragment bodies appear verbatim
# (normalized whitespace, comments stripped) as subsequences of Main.
# 6. Banned-token scan. 7. Display-item count.
import re, os, sys

ROOT = r"D:/GitHub/Bipedal_Robot"
SEC = os.path.join(ROOT, "tmp/paper_bridge/sections")
MAIN = os.path.join(ROOT, "Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex")
SUPP = os.path.join(ROOT, "Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/frontiers_SupplementaryMaterial.tex")
BIB = os.path.join(ROOT, "Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/BolenFrontiers22.bib")
FIGDIR = os.path.join(ROOT, "Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures")

raw = open(MAIN, encoding="utf-8").read()
lines = raw.splitlines()
active = [l for l in lines if not l.lstrip().startswith("%")]
t = "\n".join(active)

fail = []

# 1. environment balance
envs = {}
for m in re.finditer(r"\\(begin|end)\{([^}]*)\}", t):
    envs.setdefault(m.group(2), []).append(m.group(1))
for e, v in sorted(envs.items()):
    b, en = v.count("begin"), v.count("end")
    if b != en:
        fail.append(f"ENV IMBALANCE {e}: begin={b} end={en}")
print("ENVIRONMENTS:", {e: v.count('begin') for e, v in sorted(envs.items())})

# 2. cite keys vs bib
bib = open(BIB, encoding="utf-8", errors="replace").read()
bibkeys = set(re.findall(r"@\w+\{([^,\s]+)\s*,", bib))
keys = set()
for m in re.finditer(r"\\(?:cite|citet|citep|citealp|citeauthor)\*?(?:\[[^\]]*\])*\{([^}]*)\}", t):
    for k in m.group(1).split(","):
        if k.strip():
            keys.add(k.strip())
missing = sorted(keys - bibkeys)
print(f"\nCITE KEYS: {len(keys)} unique; MISSING FROM BIB: {missing if missing else 'NONE'}")
if missing:
    fail.append(f"MISSING CITE KEYS {missing}")

# 3. graphics targets
targets = re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]*)\}", t)
print(f"\nINCLUDEGRAPHICS ({len(targets)}):")
for g in targets:
    p = os.path.join(os.path.dirname(FIGDIR), g)
    ex = os.path.exists(p)
    print(f"  {g} -> {'EXISTS' if ex else 'MISSING'}")
    if not ex:
        fail.append(f"MISSING GRAPHIC {g}")

# 4. refs vs labels (main + supp via xr)
labels = set(re.findall(r"\\label\{([^}]*)\}", t))
supp_raw = open(SUPP, encoding="utf-8", errors="replace").read()
supp_active = "\n".join(l for l in supp_raw.splitlines() if not l.lstrip().startswith("%"))
supp_labels = set(re.findall(r"\\label\{([^}]*)\}", supp_active))
refs = set()
for m in re.finditer(r"\\(?:ref|eqref)\{([^}]*)\}", t):
    refs.add(m.group(1))
unres = sorted(refs - labels - supp_labels)
unused = sorted(labels - refs)
print(f"\nREFS: {len(refs)}; UNRESOLVED (main+supp): {unres if unres else 'NONE'}")
print(f"LABELS NEVER REFERENCED: {unused if unused else 'NONE'}")
if unres:
    fail.append(f"UNRESOLVED REFS {unres}")
dup = [l for l in labels if len(re.findall(r"\\label\{" + re.escape(l) + r"\}", t)) > 1]
print("DUPLICATE LABELS IN MAIN:", dup if dup else "NONE")
if dup:
    fail.append(f"DUPLICATE LABELS {dup}")

# 5. fragment bodies as subsequences (whitespace-normalized, comments stripped)
def norm_lines(txt):
    out = []
    for l in txt.splitlines():
        s = re.sub(r"\s+", " ", l.lstrip()).rstrip()
        if s.startswith("%") or not s:
            continue
        out.append(s)
    return out

main_norm = norm_lines(raw)
FRAGS = ["10_abstract_intro.tex", "20_methods.tex", "30_results.tex",
         "40_discussion_conclusion.tex", "tables_main.tex"]
print("\nFRAGMENT-INCLUSION CHECK (normalized-line subsequence):")
for f in FRAGS:
    frag = norm_lines(open(os.path.join(SEC, f), encoding="utf-8").read())
    # skip known intentionally-different lines: the 40_discussion ". each enlarged" fix
    it = iter(enumerate(main_norm))
    pos = -1
    missing_lines = []
    for fl in frag:
        found = False
        for j in range(pos + 1, len(main_norm)):
            if main_norm[j] == fl:
                pos = j
                found = True
                break
        if not found:
            missing_lines.append(fl[:90])
    status = "VERBATIM-INCLUDED" if not missing_lines else f"{len(missing_lines)} LINE(S) NOT FOUND"
    print(f"  {f}: {status}")
    for ml in missing_lines:
        print(f"    MISSING: {ml}")

# 6. banned tokens / stale references
for tok in ["tab:OptCoeff", "tab:TorqueFit", "fig:Ext_passthru", "$\\phi$\\qty", "$\\phi$\\;",
            "robot robot", "the the ", "offet", "gamultobj", "acheivable", "Function Fit",
            "delve", "crucial role", "Torque Values About Robotic Knee", "0.0.0958",
            "52.0", "41.5 cm"]:
    hits = [i + 1 for i, l in enumerate(lines) if tok in l and not l.lstrip().startswith("%")]
    if hits:
        print(f"TOKEN '{tok}' at active lines {hits}")
        if tok not in ("41.5 cm",):  # 41.5 cm should be gone (now \qty{41.5}{\cm})
            fail.append(f"STALE TOKEN {tok} at {hits}")

# 7. display items
nfig = t.count("\\begin{figure}")
ntab = t.count("\\begin{table}")
print(f"\nDISPLAY ITEMS: {nfig} figures + {ntab} tables = {nfig + ntab} (cap 15)")
if nfig + ntab > 15:
    fail.append("DISPLAY CAP EXCEEDED")

# 8. line endings
data = open(MAIN, "rb").read()
crlf = data.count(b"\r\n")
print(f"\nLINE ENDINGS: CRLF={crlf}, total LF={data.count(chr(10).encode())}")

print("\n" + ("ALL CHECKS PASSED" if not fail else "FAILURES: " + "; ".join(fail)))
