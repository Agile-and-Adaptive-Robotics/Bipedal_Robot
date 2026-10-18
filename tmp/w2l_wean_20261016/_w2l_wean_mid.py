r"""W2L wean midpoint 0.175 (2026-10-16): death boundary is bracketed
0.20 (pass) <-> 0.15 (fail); run ONE sharpening round at 0.175 and update
ladder_summary.json with corrected death semantics: death_scale = the
HIGHEST failing dscale (first rung that dies going down), finalist = the
LOWEST surviving dscale."""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _w2l_wean_campaign as C          # noqa: E402


def main() -> int:
    seed = json.load(open(os.path.join(HERE, "round_best_r020.json"),
                          encoding="utf-8"))["params"]
    r = C.run_round(0.175, seed)
    lad = json.load(open(os.path.join(HERE, "ladder_summary.json"),
                         encoding="utf-8"))
    if not any(abs(x["dscale"] - 0.175) < 1e-9 for x in lad["rounds"]):
        lad["rounds"].append(r)
    surv = [x for x in lad["rounds"] if x["survive30"]]
    fails = [x["dscale"] for x in lad["rounds"] if not x["survive30"]]
    lad["finalist_dscale"] = (min(surv, key=lambda x: x["dscale"])["dscale"]
                              if surv else None)
    lad["death_scale"] = max(fails) if fails else None   # first rung down that dies
    json.dump(lad, open(os.path.join(HERE, "ladder_summary.json"), "w",
                        encoding="utf-8"), indent=1)
    print(f"MID DONE finalist={lad['finalist_dscale']} "
          f"death={lad['death_scale']}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
