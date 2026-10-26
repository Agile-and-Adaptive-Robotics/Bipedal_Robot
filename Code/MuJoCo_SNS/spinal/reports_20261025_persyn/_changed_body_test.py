"""CHANGED-BODY TEST (Ben's ruling, 2026-10-26): freeze the tuned neural
controller, perturb the BIOMECHANICS (leg-segment masses scaled +/-%,
mass+inertia together), and check the same conductances still walk. This
replaces the parameter nudge test as the acceptance bar (the literature's
test: Deng torque perturbations, Shinohara belt levels/holes, Di Russo
pushes — all with frozen controllers).

PRE-REGISTERED PREDICTIONS (written BEFORE this ran, per the 2026-10-26
post-hoc correction discipline):
  P1: the w25 winner survives +/-10% leg-mass scaling in at least 5 of 6
      conditions (mass shifts swing timing smoothly; the nudge failures
      were direct conductance-timing disruptions).
  P2: at +/-20%, distal-only (tibia+foot+toe) survives more often than
      femur-only (distal mass changes swing period less than proximal).
  P3: failures, if any, are the same phase-drift mode, not falls.

Conditions: femur +/-10/20%; distal (tibia+foot+toe) +/-10/20%; all-leg
+/-10%. Winner params FROZEN from round_w25.json at wsup=0.25.
"""
from __future__ import annotations

import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
WEAN = os.path.normpath(os.path.join(
    HERE, "..", "reports_20261016", "w2l_wean"))
for p in (WEAN, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
import _w2l_wean_eval as EV  # noqa: E402

# --- monkeypatch: scale leg masses at model load, controller untouched ---
SCALE = {"factor": 1.0, "groups": ()}      # set per condition
GROUPS = {
    "femur": ("femur_L", "femur_R"),
    "distal": ("tibia_L", "tibia_R", "foot_L", "foot_R", "toe_L", "toe_R"),
}
import test_w2l_ground as T  # noqa: E402
_orig_load = T.mujoco.MjModel.from_xml_path


def patched_load(path, *a, **k):
    m = _orig_load(path, *a, **k)
    if getattr(m, "nbody", 0) and str(path).endswith("w2l_ground.xml") \
            and SCALE["factor"] != 1.0:
        names = {b for g in SCALE["groups"] for b in GROUPS[g]}
        for b in range(m.nbody):
            nm = T.mujoco.mj_id2name(m, T.mujoco.mjtObj.mjOBJ_BODY, b)
            if nm in names:
                m.body_mass[b] *= SCALE["factor"]
                m.body_inertia[b] *= SCALE["factor"]
    return m


T.mujoco.MjModel.from_xml_path = patched_load
EV.T.KV["net"] = "bilateralrg2"

d = json.load(open(os.path.join(HERE, "round_w25.json"), encoding="utf-8"))
P, WSUP = d["params"], d["wsup"]

CONDITIONS = [
    ("femur+10", ("femur",), 1.10), ("femur-10", ("femur",), 0.90),
    ("femur+20", ("femur",), 1.20), ("femur-20", ("femur",), 0.80),
    ("distal+10", ("distal",), 1.10), ("distal-10", ("distal",), 0.90),
    ("distal+20", ("distal",), 1.20), ("distal-20", ("distal",), 0.80),
    ("allleg+10", ("femur", "distal"), 1.10),
    ("allleg-10", ("femur", "distal"), 0.90),
]

rows = []
for name, groups, f in CONDITIONS:
    SCALE["factor"], SCALE["groups"] = f, groups
    kn = dict(EV.BASE, dscale=1.0, te=3.0, tf=4.0, wsup=WSUP,
              sagittal=50.0, dur_walk=30.0, **P)
    m = EV.evaluate(kn, tag=f"bodytest_{name}",
                    csv_path=os.path.join(HERE, "bodytest_trials.csv"))
    rows.append(dict(cond=name, factor=f, groups=list(groups),
                     survive=bool(m["survive"]), obj=float(m["obj"]),
                     cycL=int(m["cycL"]), cycR=int(m["cycR"]),
                     phase=m["phase_med"], tilt=float(m["tilt_max"]),
                     pelz=float(m["pelz_min"])))
    print(f"  {name:10s} survive={m['survive']} obj={m['obj']:.1f} "
          f"cyc {m['cycL']}/{m['cycR']} phase {m['phase_med']} "
          f"pelz_min {m['pelz_min']:.2f}", flush=True)
SCALE["factor"] = 1.0
n_ok = sum(r["survive"] for r in rows)
json.dump(dict(pre_registration=open(__file__, encoding="utf-8").read()
               .split('"""')[1], rows=rows, survived=n_ok,
               of=len(rows), winner="round_w25", wsup=WSUP),
          open(os.path.join(HERE, "changed_body_test_w25.json"), "w"),
          indent=1)
print(f"\nCHANGED-BODY TEST w25: {n_ok}/{len(rows)} conditions survive")
