"""Time a cold single-pass pdflatex run of the scratch dissertation copy."""
import glob
import os
import subprocess
import time

SCRATCH = r"D:\Github\Bipedal_Robot\tmp\dissfix_scratch\ProofFinal"
PDFLATEX = r"D:\MiKTeX\miktex\bin\x64\pdflatex.exe"

for f in glob.glob(os.path.join(SCRATCH, "main.*")):
    if f.endswith((".tex", ".bib")) or ("." not in os.path.basename(f)[5:]):
        pass
# remove aux state for a cold run, keep .tex/.bib/.cls/.bst
for f in os.listdir(SCRATCH):
    if f.startswith("main.") and os.path.splitext(f)[1] not in (".tex", ".bib"):
        os.remove(os.path.join(SCRATCH, f))

t0 = time.time()
p = subprocess.run(
    [PDFLATEX, "-interaction=nonstopmode", "-file-line-error", "main.tex"],
    cwd=SCRATCH, capture_output=True, text=True, timeout=560,
)
dt = time.time() - t0
print(f"cold pass exit={p.returncode} wall={dt:.1f}s")
log = open(os.path.join(SCRATCH, "main.log"), encoding="latin-1").read()
print("Output written present:", "Output written on" in log)
