"""Build Final/main.tex from easteregg2 in a new directory passed on the command line."""
import os
import argparse
import subprocess
from pathlib import Path

dissertation = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument("build",type=Path)
parser.add_argument("--source",type=Path,default=dissertation / "Final")
arguments=parser.parse_args()
source = arguments.source.resolve()
build = arguments.build.resolve()
build.mkdir(parents=True, exist_ok=False)
bin_dir = Path(r"D:\MiKTeX\miktex\bin\x64")
env = dict(os.environ)
for variable in ["BIBINPUTS", "BSTINPUTS"]:
    env[variable] = str(source) + os.pathsep + env.get(variable, "")
def run(name, args, directory):
    logfile = build / (name + ".stdout")
    with logfile.open("wb") as stream:
        result = subprocess.run(args, cwd=directory, env=env, stdout=stream, stderr=subprocess.STDOUT)
    print(name, result.returncode, logfile, flush=True)
    if result.returncode:
        print(logfile.read_bytes()[-2500:].decode(errors="replace"))
        raise SystemExit(result.returncode)

latex = [bin_dir / "pdflatex.exe", "-interaction=nonstopmode", "-halt-on-error",
         "-file-line-error", "-recorder", "-output-directory="+str(build), "main.tex"]
run("latex1", latex, source)
run("bibliography", [bin_dir / "bibtex.exe", "main"], build)
for number in range(2, 6):
    run("latex"+str(number), latex, source)
    log = (build/"main.log").read_text(encoding="utf-8", errors="replace")
    if number >= 3 and not any(s in log for s in
            ("Rerun to get", "Rerun LaTeX", "Label(s) may have changed", "undefined references")):
        break
