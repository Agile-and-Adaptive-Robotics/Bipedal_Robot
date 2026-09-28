# tex_structure_check.py - environment/brace/dollar balance in the edited tex files
import re
import io
import os

DIR = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters"
FILES = ["94-AppendixC.tex", "30-results.tex", "20-methods.tex"]

pairs = [(r"\\begin\{table\}", r"\\end\{table\}"),
         (r"\\begin\{sidewaystable\}", r"\\end\{sidewaystable\}"),
         (r"\\begin\{tabular\}", r"\\end\{tabular\}"),
         (r"\\begin\{equation\}", r"\\end\{equation\}")]

for f in FILES:
    t = io.open(os.path.join(DIR, f), encoding="utf-8").read()
    # strip LaTeX comments (incl. my provenance notes) before counting
    stripped = []
    for line in t.splitlines():
        m = re.search(r"(?<!\\)%", line)
        stripped.append(line[:m.start()] if m else line)
    body = "\n".join(stripped)
    ok = True
    msgs = []
    for b, e in pairs:
        nb, ne = len(re.findall(b, body)), len(re.findall(e, body))
        if nb != ne:
            ok = False
            msgs.append(f"{b}={nb} vs {e}={ne}")
    dollars = body.count("$")
    braces = body.count("{") - body.count("}")
    print(f"{f}: env={'OK' if ok else 'MISMATCH ' + '; '.join(msgs)} | "
          f"dollar parity {'OK' if dollars % 2 == 0 else 'ODD'} ({dollars}) | brace delta {braces}")
