"""Drive-scale ladder probe on a fixed feedback config (the wean-ladder
structural test): winner feedback knobs HELD, drive_scale swept; records
where bilateral gait dies.

Usage: python wean_probe_20261016.py <winner.json> <s1> <s2> ...
Each scale = fresh subprocess eval (`wean_tune one`, appends to
wean_human/probe.jsonl), 3 concurrent. Progress in probe.log.
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reports_20261016", "wean_human")
PY = r"D:\Anaconda\envs\myo\python.exe"

win = json.loads(open(sys.argv[1], encoding="utf-8").read())
scales = [float(s) for s in sys.argv[2:]]
# the winner's TUNED feedback must ride along (spec["fb"]) — without it
# set_stage loads s3k values and the ladder would probe the wrong gains
import wean_tune_20261016 as _WT
win_fb = {k: float(win["params"].get(k, 0.0)) for k in _WT.FB_KEYS}
jsonl = os.path.join(OUT, "probe.jsonl")
log = open(os.path.join(OUT, "probe.log"), "a", encoding="utf-8",
           buffering=1)

pending = list(scales)
running = []
while pending or running:
    while pending and len(running) < 3:
        s = pending.pop(0)
        spec = dict(tag=f"probe_s{s:g}", knobs=win["knobs"], scale=s,
                    jsonl=jsonl, fb=dict(win_fb))
        sp = os.path.join(OUT, f"_probe_spec_{s:g}.json")
        json.dump(spec, open(sp, "w", encoding="utf-8"))
        wl = open(os.path.join(OUT, f"_probe_s{s:g}.log"), "w",
                  encoding="utf-8", buffering=1)
        pr = subprocess.Popen([PY, "-u", os.path.join(
            HERE, "wean_tune_20261016.py"), "one", sp], stdout=wl,
            stderr=wl, cwd=HERE,
            env={**os.environ, "PYTHONIOENCODING": "utf-8"})
        running.append((s, pr, wl))
        log.write(f"launched scale {s:g} pid {pr.pid}\n")
        time.sleep(5)
    time.sleep(10)
    for item in list(running):
        s, pr, wl = item
        if pr.poll() is not None:
            wl.close()
            running.remove(item)
            log.write(f"scale {s:g} EXIT {pr.returncode}\n")
log.write("PROBE_DONE\n")
log.close()
print("probe rows in " + jsonl)
