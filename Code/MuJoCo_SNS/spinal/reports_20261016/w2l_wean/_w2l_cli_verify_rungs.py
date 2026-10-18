r"""CLI-verify EVERY ladder rung winner in a FRESH process per rung
(2026-10-16): the r020 'survivor' died in a fresh CLI process (campaign
in-process context artifact) -- so the robust wean frontier must be
established by reproducing each rung the way any user/supervisor would
run it. One subprocess per rung = one identical CLI context per rung.

Writes cli_verify_rungs.json + appends rows to wean_trials.csv."""
from __future__ import annotations

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PY = r"D:\Anaconda\envs\myo\python.exe"


def main() -> int:
    lad = json.load(open(os.path.join(HERE, "ladder_summary.json"),
                         encoding="utf-8"))
    out = []
    for r in sorted(lad["rounds"], key=lambda x: -x["dscale"]):
        ds = r["dscale"]
        p = r["params"]
        cmd = [PY, os.path.join(HERE, "_w2l_wean_eval.py"),
               f"--id=cli_verify_r{ds:.3f}", "--dur=30",
               f"--dscale={ds}",
               f"--g-heel={p['g_heel']}", f"--g-ib={p['g_ib']}",
               f"--ibnA={p['ibnA']}", f"--ctn-amp={p['ctn_amp']}",
               f"--cref={p['cref']}", f"--comm={p['comm']}",
               f"--toe-m={p['toe_m']}"]
        res = subprocess.run(cmd, capture_output=True, text=True,
                             encoding="utf-8", errors="replace",
                             cwd=HERE, timeout=600)
        line = next((ln for ln in res.stdout.splitlines()
                     if ln.startswith("RESULT ")), None)
        m = json.loads(line[7:]) if line else dict(error=res.stderr[-300:])
        out.append(dict(dscale=ds, cli_survive=m.get("survive"),
                        cli_metrics=m))
        print(f"rung {ds:.3f}: CLI survive={m.get('survive')} "
              f"cyc {m.get('cycL')}/{m.get('cycR')} "
              f"slope {m.get('slope')} phase {m.get('phase_med')} "
              f"obj {m.get('obj')}", flush=True)
        json.dump(out, open(os.path.join(HERE, "cli_verify_rungs.json"),
                            "w", encoding="utf-8"), indent=1)
    robust = [o["dscale"] for o in out if o["cli_survive"]]
    print("ROBUST RUNGS:", robust, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
