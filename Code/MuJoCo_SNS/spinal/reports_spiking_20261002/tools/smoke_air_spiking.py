"""GATE 4 (SPIKING_MIRROR_PLAN.md): 20 s air smoke, spiking vs
non-spiking (_smoke_nap.py pattern, regime match not bit-exact).

Runs the FULL runner (stand->ramp->air-stepping, deafferented,
interleg off - the stage-1 air preparation) once per network flavor:
  ns      : AARL_NET unset  (the default non-spiking build)
  spiking : AARL_NET=spiking (the hybrid mirror through the selector)
Both use the SAME stage-1 params (best_walk_params_v10 multipliers,
drive 2.5, rg_nap_h 0.35) and write scratch npz files; this script then
prints the regime metrics (RG period, E-duty, knee/hip/ankle ranges,
rises) side by side and, when BOTH runs' jsons exist, the comparison
verdict.

Regime-match rule (documented; the plan only asks for a comparison, not
identity): PASS = spiking run finite AND >= 3 rises AND RG period within
35% of the non-spiking run AND knee and hip swing ranges >= 50% of the
non-spiking run.

Usage: D:\\Anaconda\\envs\\myo\\python.exe smoke_air_spiking.py ns|spiking|compare
"""
import json
import os
import sys
from pathlib import Path

SPINAL = Path(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
os.chdir(SPINAL)
sys.path.insert(0, str(SPINAL))

mode = sys.argv[1] if len(sys.argv) > 1 else "compare"
OUT = SPINAL / "reports_spiking_20261002"

if mode in ("ns", "spiking"):
    if mode == "spiking":
        os.environ["AARL_NET"] = "spiking"
    else:
        os.environ.pop("AARL_NET", None)
    os.environ["AARL_NPZ"] = f"scratch_spk_air_{mode}.npz"
    import _curriculum as C
    import runner as R
    mul = json.loads((SPINAL / "best_walk_params_v10.json")
                     .read_text(encoding="utf-8"))["multipliers"]
    p = {**mul, "drive": 2.5, "rg_nap_h": 0.35, "desc_e": 1.7,
         "desc_f": 1.4, "rg_to_pf": 2.4}
    C.set_stage(1, p)
    R.main(["--no-ground", "--no-afferents", "--no-interleg",
            "--time", "20", "--drive", "2.5"])

    import numpy as np
    z = np.load(SPINAL / os.environ["AARL_NPZ"], allow_pickle=True)
    t, q, neuro = z["t"], z["q"], z["neuro"]
    m = (t >= 5.0) & (t <= 23.0)
    res = dict(net=mode, npz=os.environ["AARL_NPZ"])
    res["finite"] = bool(np.all(np.isfinite(q[m]))
                         and np.all(np.isfinite(neuro[m])))
    knee = q[m, 4]
    rge = neuro[m, 2]
    rgf = neuro[m, 3]
    on = rge > 0.5 * max(rge.max(), 1e-9)
    rises = np.flatnonzero(np.diff(on.astype(int)) == 1)
    res["rises"] = int(len(rises))
    res["knee_min"] = float(knee.min())
    res["knee_range"] = float(knee.ptp())
    res["hip_range"] = float(q[m, 3].ptp())
    res["ankle_range"] = float(q[m, 5].ptp())
    res["e_duty"] = float(on.mean())
    if len(rises) >= 2:
        per = float(np.mean(np.diff(t[m][rises])))
        res["period_s"] = per
    d = (rge - rgf)[rises[0]:rises[-1] + 1] if len(rises) >= 2 else None
    (OUT / f"air_smoke_{mode}.json").write_text(
        json.dumps(res, indent=1), "utf-8")
    print("AIR SMOKE " + json.dumps(res))
    sys.exit(0)

# ---- compare ----
a = json.loads((OUT / "air_smoke_ns.json").read_text(encoding="utf-8"))
b = json.loads((OUT / "air_smoke_spiking.json")
               .read_text(encoding="utf-8"))
print(f"{'metric':14s} {'non-spiking':>12s} {'spiking':>12s}")
for k in ("finite", "rises", "period_s", "e_duty", "knee_min",
          "knee_range", "hip_range", "ankle_range"):
    print(f"{k:14s} {str(a.get(k)):>12s} {str(b.get(k)):>12s}")
T_ns, T_sp = a.get("period_s"), b.get("period_s")
per_ok = T_ns and T_sp and abs(T_sp - T_ns) / T_ns <= 0.35
rng_ok = (b["knee_range"] >= 0.5 * a["knee_range"]
          and b["hip_range"] >= 0.5 * a["hip_range"])
ok = bool(b["finite"] and b["rises"] >= 3 and per_ok and rng_ok)
print(f"finite={b['finite']} rises>=3={b['rises'] >= 3} "
      f"period-35%={bool(per_ok)} ranges-50%={bool(rng_ok)}")
print("GATE4 AIR SMOKE:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
