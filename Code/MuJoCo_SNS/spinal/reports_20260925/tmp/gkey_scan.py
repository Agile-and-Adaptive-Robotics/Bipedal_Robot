"""Scan the two variant builders for the params.G / params.TAU keys they
read (goal4 wiring: the curriculum must search functional knobs only).
Also list the desc_* knobs the runner applies to the descending drive.
"""
import re
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

for fn in ("build_network_w2lvar.py", "build_network_syn6.py"):
    src = open(fn, encoding="utf-8").read()
    gkeys = sorted(set(re.findall(r'G\["([a-z0-9_]+)"\]', src)))
    tau = sorted(set(re.findall(r'TAU\["([a-z0-9_]+)"\]', src)))
    print(f"--- {fn}")
    print("G keys read:", gkeys)
    print("TAU keys read:", tau)

rsrc = open("runner.py", encoding="utf-8").read()
gkeys_r = sorted(set(re.findall(r'G\["([a-z0-9_]+)"\]', rsrc)))
print("--- runner.py G keys read:")
print(gkeys_r)
for pat in ("desc_e", "desc_f", "rg_nap_h"):
    hits = [i + 1 for i, ln in enumerate(rsrc.splitlines()) if pat in ln]
    print(f"runner.py lines containing {pat!r}:", hits)
