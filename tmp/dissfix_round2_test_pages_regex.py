"""Round-2 unit test: the patched pages regex in tmp/digest_20261002/compile_gate.py
against the REAL main.log files left behind by the two failed gate runs (the gate only
deletes its build dir on PASS, so both are still in %TEMP%). Regex-only -- does not run
pdflatex or the gate. Also runs the gate's errors/undefined/multiply line checks on the
same logs, and a negative control (log with the Output-written line stripped must NOT match).
"""
import os
import re
from pathlib import Path

PATTERN = re.compile(r"Output written on .*?\((\d+) page", re.S)

logs = [Path(os.environ["TEMP"]) / d / "main.log" for d in ("pf_gate_j_zoe3dv", "pf_gate_lvwgels7")]
for lg in logs:
    log = lg.read_text(encoding="utf8", errors="replace")
    errs = [l for l in log.splitlines() if l.startswith("!")]
    undef = [l for l in log.splitlines() if "undefined" in l.lower()]
    multi = [l for l in log.splitlines() if "multiply" in l.lower()]
    m = PATTERN.search(log)
    print(f"{lg.parent.name}: errors={len(errs)} undef={len(undef)} multi={len(multi)} "
          f"pages_match={bool(m)} pages={m.group(1) if m else None}")

# negative control: strip the Output written sentence -> must NOT match
log2 = logs[0].read_text(encoding="utf8", errors="replace")
log2 = re.sub(r"Output written on.*?\(\d+ pages?[^\n]*\n?", "", log2, flags=re.S)
print("negative control (line stripped):", bool(PATTERN.search(log2)))
