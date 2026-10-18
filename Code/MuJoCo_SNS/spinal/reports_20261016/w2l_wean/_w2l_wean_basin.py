r"""W2L robot WEAN finalist: BASIN GATE + GAIN-SENSITIVITY PROBE
(overnight campaign 2026-10-16, egg; stage "feedback optimizer (robot)").

Reads a round_best_rXXX.json from _w2l_wean_campaign.py, then:

(1) BASIN GATE (the human-side basin_gate.py recipe, robot adaptation):
    12 seeded perturbations, ALL 7 tuned scalars jointly x (1 +- 1%),
    full walk re-run at the finalist drive scale, 30 s. PASS = every
    perturbed run still SURVIVES (finite, no fall, >= 6 primary cycles
    both sides, inter-leg phase 0.35-0.65, fwd > 0.05 m). Margin logged.

(2) GAIN-SENSITIVITY PROBE (Ben's anti-lottery requirement): each tuned
    pathway removed ONE AT A TIME (gain -> 0; comm=0 is the documented
    structural ablation in build_w2l_split_net), plus the all-Ben-defaults
    variant, each at the finalist drive scale, 30 s. REMOVING A REAL GAIN
    MUST DEGRADE THE RESULT -- a knob whose removal does NOT degrade is
    reported as inert (a finding, not a pass).

Usage:
    D:\Anaconda\envs\myo\python.exe _w2l_wean_basin.py round_best_r060.json
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _w2l_wean_eval as EV          # noqa: E402
from _w2l_wean_campaign import BOOT_PARAMS, knobs_of  # noqa: E402

TUNED = ("g_heel", "g_ib", "ibnA", "ctn_amp", "cref", "comm", "toe_m")
NEUTRAL = dict(g_heel=0.0, g_ib=0.0, ibnA=0.0, ctn_amp=0.0, comm=0.0,
               toe_m=0.0)


def mrow(m: dict) -> dict:
    return {k: m[k] for k in EV.METRIC_KEYS}


def main() -> int:
    fin = json.load(open(os.path.join(HERE, sys.argv[1]),
                         encoding="utf-8"))
    ds = float(fin["dscale"])
    p0 = dict(fin["params"])
    path = os.path.join(HERE, "basin_results.json")
    out = dict(finalist=sys.argv[1], dscale=ds, params=p0)

    # anchor: exact finalist params re-run (should match round verify)
    m0 = EV.evaluate(knobs_of(ds, p0, 30.0), tag="finalist_anchor")
    out["anchor"] = mrow(m0)

    # ---- (1) basin gate: +-1% on all 7 scalars jointly, 12 draws
    import numpy as np
    rng = np.random.default_rng(20261016)
    basin = []
    for j in range(12):
        pj = dict(p0)
        for k in TUNED:
            sgn = float(rng.choice((-1.0, 1.0)))
            pj[k] = float(p0[k]) * (1.0 + 0.01 * sgn)
        mj = EV.evaluate(knobs_of(ds, pj, 30.0), tag=f"basin{j:02d}")
        basin.append(mrow(mj))
    out["basin_runs"] = basin
    out["basin_pass"] = bool(all(b["survive"] for b in basin))
    out["basin_margin_cyc"] = min(min(b["cycL"], b["cycR"]) for b in basin)
    out["basin_margin_fwd"] = min(b["fwd"] for b in basin)
    json.dump(out, open(path, "w", encoding="utf-8"), indent=1)

    # ---- (2) sensitivity: one pathway removed at a time + Ben defaults
    # NOTE comm: the aff-net census FORBIDS comm=0 (conditional topology:
    # comm=0 unbuilt commissurals change the coupled-baseline census), so
    # the comm "removal" probe uses 0.001 = 0.24% of the tuned 0.418 --
    # effectively silenced without the structural rebuild.
    sens = {}
    for k in TUNED:
        if k == "cref":
            continue        # not a gain: cref is the encoder saturation
        pk = dict(p0)
        pk[k] = 0.001 if k == "comm" else NEUTRAL[k]
        try:
            sens[f"zero_{k}"] = mrow(EV.evaluate(knobs_of(ds, pk, 30.0),
                                                 tag=f"sens_zero_{k}"))
        except AssertionError as e:
            sens[f"zero_{k}"] = dict(error=str(e)[:120])
        json.dump(out | {"sensitivity": sens},
                  open(path, "w", encoding="utf-8"), indent=1)
    sens["ben_defaults"] = mrow(EV.evaluate(knobs_of(ds, BOOT_PARAMS, 30.0),
                                            tag="sens_ben_defaults"))
    out["sensitivity"] = sens
    out["anchor_obj"] = m0["obj"]

    path = os.path.join(HERE, "basin_results.json")
    json.dump(out, open(path, "w", encoding="utf-8"), indent=1)
    print(f"BASIN {'PASS' if out['basin_pass'] else 'FAIL'} "
          f"(min cyc {out['basin_margin_cyc']}, min fwd "
          f"{out['basin_margin_fwd']:+.3f}); sensitivity objs: "
          + ", ".join(f"{k}={v['obj']:.0f}" for k, v in sens.items())
          + f" | anchor {m0['obj']:.0f} -> {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
