# Supervisor Gate-1 check for 70-overall.tex revision batch (campaign 2026-10-15).
# Extracts text from Ben's two papers + the 2020 course report, then checks every
# 10-word window of the batch's CHANGED lines (git diff, + side) for verbatim
# overlap, and pulls 'Hz' context from the course report for the 4 Hz claim.
import subprocess, re, sys
from pathlib import Path
from pypdf import PdfReader

REPO = Path(r"D:\Github\Bipedal_Robot")
GIT = r"C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.7\resources\app\git\cmd\git.exe"
TARGET = "Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/70-overall.tex"

PDFS = {
    "bolen_2026_Actuators": REPO / "Documentation/Reports and Papers/Isometric_Force_Test/IsometricForceCharacterizationofBPAs.pdf",
    "bolen_isometric_2025_ms": REPO / "Documentation/ Reports and Papers/Knee_Torque_Test/Torque_Values_About_Robotic_Knee_Using_Braided_Pneumatic_Actuators.pdf",
    "bolen_design_2020_report": REPO / "Documentation/Reports and Papers/CPG_Walker/Design_of_a_Bipedal_Walker_Controlled_by_a_CPG_Network_in_Animatlab.pdf",
}
# fix the accidental space in the second path
PDFS["bolen_isometric_2025_ms"] = REPO / "Documentation/Reports and Papers/Knee_Torque_Test/Torque_Values_About_Robotic_Knee_Using_Braided_Pneumatic_Actuators.pdf"

def words(text):
    return re.findall(r"[A-Za-z0-9]+", text.lower())

corpora = {}
for name, p in PDFS.items():
    if not p.exists():
        print(f"MISSING PDF: {name}: {p}")
        continue
    r = PdfReader(str(p))
    txt = "\n".join((pg.extract_text() or "") for pg in r.pages)
    (REPO / "tmp").mkdir(exist_ok=True)
    (REPO / "tmp" / f"sup70_{name}.txt").write_text(txt, encoding="utf-8")
    corpora[name] = txt
    print(f"extracted {name}: {len(r.pages)} pages, {len(words(txt))} words")

# changed lines from the diff (+ side, minus LaTeX-noise headers)
diff = subprocess.run([GIT, "-C", str(REPO), "diff", "--", TARGET],
                      capture_output=True, text=True, encoding="utf-8", errors="replace").stdout
changed = [l[1:] for l in diff.splitlines() if l.startswith("+") and not l.startswith("+++")]
changed_text = " ".join(changed)
cw = words(changed_text)
print(f"\nchanged lines: {len(changed)}, words: {len(cw)}")

N = 10
paper_streams = {n: words(t) for n, t in corpora.items()}
hits = []
for i in range(len(cw) - N + 1):
    win = cw[i:i+N]
    key = " ".join(win)
    for n, stream in paper_streams.items():
        if key in " ".join(stream):
            hits.append((n, key))
print(f"\n10-word verbatim overlaps vs Ben's papers: {len(hits)}")
for n, h in hits:
    print(f"  HIT [{n}]: {h}")

# 7-word external check vs course report (the only non-Ben-owned-adjacent doc here
# is the course report, co-authored with Morrow -> still Ben's own; kept at 10)
# Hz context in the course report (4 Hz claim source)
rep = corpora.get("bolen_design_2020_report", "")
print("\n--- Hz mentions in bolen_design_2020 course report ---")
for m in re.finditer(r"[^.]*\bHz\b[^.]*\.", rep):
    s = " ".join(m.group(0).split())
    if len(s) < 400:
        print(" *", s)
print("\n--- Hz mentions in 70-overall changed lines ---")
for l in changed:
    if "Hz" in l or "Hz}" in l:
        for m in re.finditer(r"[^.]*Hz[^.]*\.", l):
            print(" *", " ".join(m.group(0).split())[:300])
