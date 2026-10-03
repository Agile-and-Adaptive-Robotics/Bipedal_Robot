"""Compile gate for the dissertation ProofFinal. Builds in a temp dir (ProofFinal stays
clean), full pdflatex/bibtex cycle, then scans the log. Exit 0 + 'GATE PASS pages=N' on a
clean build; exit 1 + the offending lines otherwise."""
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

PF = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal")
PDFLATEX = r"D:\MiKTeX\miktex\bin\x64\pdflatex.exe"
BIBTEX = r"D:\MiKTeX\miktex\bin\x64\bibtex.exe"

build = Path(tempfile.mkdtemp(prefix="pf_gate_"))
for bib in ("thesis.bib", "BolenFrontiers22.bib", "Frontiers-Harvard.bst"):
    src = PF / bib
    if src.exists():
        shutil.copy(src, build / bib)


def run(cmd, args, timeout=420):
    return subprocess.run([cmd, *args], cwd=str(PF), capture_output=True, text=True,
                          encoding="utf8", errors="replace", timeout=timeout)


fail = None
try:
    r1 = run(PDFLATEX, ["-interaction=nonstopmode", f"-output-directory={build}", "main.tex"])
    rb = subprocess.run([BIBTEX, "main"], cwd=str(build), capture_output=True, text=True,
                        encoding="utf8", errors="replace", timeout=180)
    r2 = run(PDFLATEX, ["-interaction=nonstopmode", f"-output-directory={build}", "main.tex"])
    r3 = run(PDFLATEX, ["-interaction=nonstopmode", f"-output-directory={build}", "main.tex"])
    log = (build / "main.log").read_text(encoding="utf8", errors="replace")
    errors = [l for l in log.splitlines() if l.startswith("!")]
    undef = [l for l in log.splitlines() if "undefined" in l.lower()]
    multi = [l for l in log.splitlines() if "multiply" in l.lower()]
    pages = re.search(r"Output written on .*?\((\d+) page", log, re.S)
    if errors:
        fail = "ERRORS:\n" + "\n".join(errors[:20])
    elif undef:
        fail = "UNDEFINED:\n" + "\n".join(undef[:20])
    elif multi:
        fail = "MULTIPLY-DEFINED:\n" + "\n".join(multi[:20])
    elif not pages:
        fail = "no 'Output written' line found; last exit " + str(r3.returncode)
except subprocess.TimeoutExpired:
    fail = "timeout during compile"

if fail:
    print("GATE FAIL")
    print(fail)
    sys.exit(1)
print(f"GATE PASS pages={pages.group(1)}")
shutil.rmtree(build, ignore_errors=True)
