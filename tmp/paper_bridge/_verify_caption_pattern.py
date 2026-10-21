import re, sys
path = r"D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex"
text = open(path, encoding="utf-8").read()
lines = text.splitlines()
# E1 sweep: a value (number with unit or qty{...}) immediately followed by a panel letter (A)-(H),
# i.e. the banned value-then-letter pattern like "48.5 cm (A)" or "\qty{48.5}{\cm} \textbf{(A)}"
banned = re.compile(
    r"(\d\s*(?:cm|mm|m|kPa|MPa|N|N\\cdot m|\\newton\w*)\s*\}?\s*|\d\s+)"  # a value w/ unit (incl. qty tail)
    r"(\\textbf\{\()?[（(]\s*[A-H]\s*[)）]"
)
hits = []
for i, ln in enumerate(lines, 1):
    if "\\caption" in ln or "\\ multicolumn" in ln:
        pass
    for m in banned.finditer(ln):
        hits.append((i, m.group(0), ln[:160]))
print("VALUE-THEN-LETTER HITS:", len(hits))
for h in hits:
    print(h)
# Also list every caption line number for the letter-first audit
for i, ln in enumerate(lines, 1):
    if "\\caption{" in ln:
        # find panel-letter usages with values
        panels = re.findall(r"\\textbf\{\(([A-H])\)\}\s*([^,.;]{0,40})", ln)
        print(f"CAPTION line {i}: panel-prefix usages: {panels[:8]}")
