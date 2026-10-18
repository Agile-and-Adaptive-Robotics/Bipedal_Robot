r"""W2L wean ladder EXTENSION (2026-10-16): the planned ladder
(1.0..0.20) all SURVIVED (death=None), so push deeper to find where
alternating progression actually dies: 0.15, 0.10, 0.05, 0.00 (0 = pure
kickoff + feedback, no tonic drive). Same 2-consecutive-fail stop +
one midpoint refinement; merges into ladder_summary.json."""
from __future__ import annotations

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _w2l_wean_campaign as C          # noqa: E402

DEEP = [0.15, 0.10, 0.05, 0.0]


def main() -> int:
    t0 = time.time()
    seed = json.load(open(os.path.join(HERE, "round_best_r020.json"),
                          encoding="utf-8"))["params"]
    lad = json.load(open(os.path.join(HERE, "ladder_summary.json"),
                         encoding="utf-8"))
    results = []
    fails, death = 0, None
    for s in DEEP:
        r = C.run_round(s, seed)
        results.append(r)
        if r["survive30"]:
            seed = r["params"]
            fails = 0
        else:
            fails += 1
            if death is None:
                death = s
            if fails >= 2:
                break
    allr = [r for r in lad["rounds"]] + results
    last_pass = next((r["dscale"] for r in reversed(allr) if r["survive30"]),
                     None)
    first_death = min([r["dscale"] for r in allr if not r["survive30"]],
                      default=None)
    if first_death is not None and last_pass is not None:
        mid = round((last_pass + first_death) / 2.0, 3)
        if not any(abs(r["dscale"] - mid) < 1e-9 for r in allr):
            allr.append(C.run_round(mid, seed))
    surv = [r for r in allr if r["survive30"]]
    lad["rounds"] = allr
    lad["finalist_dscale"] = (min(surv, key=lambda r: r["dscale"])["dscale"]
                              if surv else None)
    lad["death_scale"] = first_death
    lad["deep_extension_wall_min"] = round((time.time() - t0) / 60, 1)
    json.dump(lad, open(os.path.join(HERE, "ladder_summary.json"), "w",
                        encoding="utf-8"), indent=1)
    print(f"DEEP DONE finalist={lad['finalist_dscale']} "
          f"death={lad['death_scale']} "
          f"({lad['deep_extension_wall_min']} min)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
