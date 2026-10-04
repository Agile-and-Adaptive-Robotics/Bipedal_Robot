"""Replicate the compile-gate invocation exactly (cwd=PF, output-directory=empty temp).

Diagnostic only: does NOT evaluate pass/fail the way the gate does; prints the
tail of the resulting log and whether 'Output written' appears, so the failure
mechanism (if any) is visible.
"""
import shutil
import subprocess
import tempfile
from pathlib import Path

PF = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal")
PDFLATEX = r"D:\MiKTeX\miktex\bin\x64\pdflatex.exe"

build = Path(tempfile.mkdtemp(prefix="pf_repro_"))
print("build dir:", build)

r1 = subprocess.run(
    [PDFLATEX, "-interaction=nonstopmode", f"-output-directory={build}", "main.tex"],
    cwd=str(PF), capture_output=True, text=True, encoding="utf8", errors="replace", timeout=420,
)
print("r1 exit:", r1.returncode)
print("--- r1 stdout first 5 lines ---")
print("\n".join(r1.stdout.splitlines()[:5]))
log = build / "main.log"
print("log exists:", log.exists())
if log.exists():
    text = log.read_text(encoding="latin-1", errors="replace")
    lines = text.splitlines()
    bang = [l for l in lines if l.startswith("!")]
    print("'!' lines:", bang[:10])
    print("Output written:", any("Output written on" in l for l in lines))
    print("--- log tail 15 ---")
    print("\n".join(lines[-15:]))
print("--- build dir contents ---")
for p in sorted(build.rglob("*")):
    print(" ", p.relative_to(build))
shutil.rmtree(build, ignore_errors=True)
