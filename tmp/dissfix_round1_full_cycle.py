"""Full gate-sequence verification in the scratch copy (clean state):
pdflatex -> bibtex -> pdflatex -> pdflatex, mirroring compile_gate.py's cycle.
Checks per final log: no '!' lines, no undefined, no multiply, Output written.
"""
import os
import subprocess
import time

SCRATCH = r"D:\Github\Bipedal_Robot\tmp\dissfix_scratch\ProofFinal"
PDFLATEX = r"D:\MiKTeX\miktex\bin\x64\pdflatex.exe"
BIBTEX = r"D:\MiKTeX\miktex\bin\x64\bibtex.exe"

for f in os.listdir(SCRATCH):
    if f.startswith("main.") and os.path.splitext(f)[1] not in (".tex", ".bib"):
        os.remove(os.path.join(SCRATCH, f))
chaux = os.path.join(SCRATCH, "chapters")
for f in os.listdir(chaux):
    if f.endswith(".aux"):
        os.remove(os.path.join(chaux, f))


def pdf_pass(label):
    t0 = time.time()
    r = subprocess.run([PDFLATEX, "-interaction=nonstopmode", "main.tex"],
                       cwd=SCRATCH, capture_output=True, text=True, timeout=420)
    print(f"{label}: exit={r.returncode} wall={time.time()-t0:.1f}s")
    return r.returncode


pdf_pass("pass1")
rb = subprocess.run([BIBTEX, "main"], cwd=SCRATCH, capture_output=True, text=True, timeout=180)
bbl = os.path.isfile(os.path.join(SCRATCH, "main.bbl"))
print(f"bibtex: exit={rb.returncode} bbl_written={bbl}")
if not bbl:  # retry once -- the transient flake observed earlier this session
    rb = subprocess.run([BIBTEX, "main"], cwd=SCRATCH, capture_output=True, text=True, timeout=180)
    bbl = os.path.isfile(os.path.join(SCRATCH, "main.bbl"))
    print(f"bibtex-retry: exit={rb.returncode} bbl_written={bbl}")
pdf_pass("pass2")
pdf_pass("pass3")

log = open(os.path.join(SCRATCH, "main.log"), encoding="latin-1").read()
lines = log.splitlines()
print("'!' lines:", sum(1 for l in lines if l.startswith("!")))
print("undefined lines:", sum(1 for l in lines if "undefined" in l.lower()))
print("multiply lines:", sum(1 for l in lines if "multiply" in l.lower()))
ow = [l for l in lines if l.startswith("Output written on")]
print("Output written:", ow)
