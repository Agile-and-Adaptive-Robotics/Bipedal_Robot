"""Bias + stretch-reflex posture test (2026-10-21, Ben's REFLEXIVE-balance
directive) - the filling-glass experiment.

Method (directive): alpha bias = the per-muscle PRELOAD-FREE static-solve
activations (gravity-only NNLS at the standing pose - the bsolve/static-opt
backsolve pattern, runner --reflex-stand); Ia/II/Ib loops ON with the
repo's literature-based encoder gains (params.AFF, SADb framing; II
setpoint moved to the standing-pose length = the SCONE T3a MuscleReflex
L0/KL analog); rig OFF (--rig-scale 0); gamma coactivation available
coupled to the alpha bias (default 0). Supraspinal BAL_* surrogates OFF -
only spinal reflexes may hold.

Acceptance (filling glass): trunk load stepped 5/10/20 %BW; the stretch
reflex should compensate (per-muscle activation rises with load, COM/ankle
drift stays small). Records COM drift, sway, ankle angle, per-muscle
activation response per load step.

Matrix (12 s runs, load steps from t=3 s, 3 s dwell):
  A bias-only   : --reflex-stand --no-afferents (alpha bias alone, rig 0)
  B bias+reflex : --reflex-stand, default AFF gains
  C/D           : B with ii_gain 3 / 6 (the postural II loop-gain question)
  E             : D + gamma_coact 0.5
  S             : B at rig_scale 1.0 (supported reference: reflex response
                  exists without the free-pelvis problem)
Usage (myo env, cwd=spinal):  python _reflex_stand_20261021.py [A|B|C|D|E|S|ALL]
Results: _reflex_stand_results_20261021.json + per-run npz/log files.
"""
import io
import json
import os
import sys
import time

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

MODE = sys.argv[1].upper() if len(sys.argv) > 1 else "ALL"

os.environ["AARL_NPZ"] = "_reflexstand_boot_20261021.npz"
for _e in ("AARL_NET", "AARL_STANCE_W", "AARL_KNEE_BIAS", "AARL_PELVIS_TY",
           "AARL_POSE", "AARL_DEAFF"):
    os.environ.pop(_e, None)

import numpy as np  # noqa: E402
import mujoco  # noqa: E402
import runner as R  # noqa: E402
import params as P  # noqa: E402
from muscle_map import classify  # noqa: E402

_m = mujoco.MjModel.from_xml_path(str(R.MODEL))
BW = float(_m.body_mass.sum()) * 9.81
ACT_NAMES = [mujoco.mj_id2name(_m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
             for i in range(_m.nu)]
del _m
GROUPS = {}
for j, nm in enumerate(ACT_NAMES):
    mi = classify(nm)
    if mi is not None:
        GROUPS.setdefault(mi.groups[0], []).append(j)

CFG = {
    "A": dict(ii_gain=1.2, gamma=0.0, aff=True,  rig=0.0),
    "B": dict(ii_gain=1.2, gamma=0.0, aff=True,  rig=0.0),
    "C": dict(ii_gain=3.0, gamma=0.0, aff=True,  rig=0.0),
    "D": dict(ii_gain=6.0, gamma=0.0, aff=True,  rig=0.0),
    "E": dict(ii_gain=6.0, gamma=0.5, aff=True,  rig=0.0),
    "S": dict(ii_gain=1.2, gamma=0.0, aff=True,  rig=1.0),
    # rig-weaned ladder (2026-10-21 evening): rig 0 collapsed for every
    # config (A broke at t=1.23, B fell by t=2); rig 1 survives but the
    # springs carry the glass load (feet unload under it). rig 0.3 = the
    # legs carry most weight while the release transient stays survivable;
    # ii_gain swept to the SCONE-KL-mapped postural scale (~30).
    "F": dict(ii_gain=1.2, gamma=0.0, aff=True,  rig=0.3),
    "G": dict(ii_gain=30.0, gamma=0.0, aff=True, rig=0.3),
    "H": dict(ii_gain=30.0, gamma=0.5, aff=True, rig=0.3),
    # J (the headroom experiment): F/G/H showed the antigravity MNs run
    # at 0.78-1.0 (glut_max saturated at 1.0) under the RG-E->PF-E
    # extensor tonic + 6x bias - no headroom for reflex expression, and
    # d_soleus stays +0.002 even at ii_gain 30. J removes the PF->MN
    # tonic (G.pf_to_mn=0) so the per-muscle alpha bias + Ia/II/Ib are
    # the ONLY drive: does bias+reflex then hold, and does the
    # activation respond to the load steps?
    "J": dict(ii_gain=30.0, gamma=0.5, aff=True, rig=0.3),
}
PF_ZERO = {"J"}   # extra: P.G["pf_to_mn"] = 0 for these tags
T0, DWELL, DUR = 3.0, 3.0, 12.0
LOADS = [5.0, 10.0, 20.0]
WATCH = ("soleus_r", "med_gas_r", "vas_lat_r", "glut_max2_r", "ercspn_r",
         "tib_ant_r", "rect_fem_r", "iliacus_r")


def run_cfg(tag):
    c = CFG[tag]
    # per-run full knob reset (leak-proof across in-process runs)
    P.G["posture_to_mn"] = 1.0
    P.G["ii_set_pose"] = 0.0
    P.G["gamma_coact"] = 0.0
    P.G["pf_to_mn"] = 2.0
    P.AFF["ii_gain"] = 1.2
    npz = f"_reflexstand_{tag}_20261021.npz"
    os.environ["AARL_NPZ"] = npz
    args = ["--stand-eval", repr(DUR), "--rig-scale", repr(c["rig"]),
            "--trunk-load", ",".join(str(int(x)) for x in LOADS),
            "--reflex-stand"]
    if c["gamma"] > 0:
        args += ["--gamma-coact", repr(c["gamma"])]
    if tag == "A":
        args += ["--no-afferents"]
    # set AFTER flag defaults so they are the values of record
    P.AFF["ii_gain"] = c["ii_gain"]
    if tag in PF_ZERO:
        P.G["pf_to_mn"] = 0.0
        print(f"[{tag}] pf_to_mn=0: PF tonic removed (headroom test)",
              flush=True)
    t0 = time.time()
    m = R.main(args)
    return m, npz, time.time() - t0


def analyze(tag, npz, m, wall):
    d = np.load(npz)
    t, com, q, ctrl, contact = (d["t"], d["com"], d["q"], d["ctrl"],
                                d["contact"])
    n = min(int(np.sum(t > 0)) + 1, len(t))
    t, com, q, ctrl, contact = (t[:n], com[:n], q[:n], ctrl[:n],
                                contact[:n])
    iw = {jn: i for i, jn in enumerate(R.KEY_JOINTS)}

    def win(lo, hi):
        a = min(int(np.searchsorted(t, lo)), n - 1)
        b = min(int(np.searchsorted(t, hi)), n - 1)
        return slice(a, max(b, a + 1))

    base = win(1.8, T0 - 0.2)
    com_ref = com[base].mean(axis=0)
    ank_r0 = float(q[base, iw["ankle_angle_r"]].mean())
    ank_l0 = float(q[base, iw["ankle_angle_l"]].mean())
    steps = []
    for i, pct in enumerate(LOADS):
        s = win(T0 + i * DWELL + 1.0, T0 + (i + 1) * DWELL)
        cc = com[s]
        steps.append(dict(
            pctBW=pct,
            com_dx=round(float(cc[:, 0].mean() - com_ref[0]), 4),
            com_dy=round(float(cc[:, 1].mean() - com_ref[1]), 4),
            com_dz=round(float(cc[:, 2].mean() - com_ref[2]), 4),
            sway_max=round(float(np.max(np.hypot(
                cc[:, 0] - com_ref[0], cc[:, 1] - com_ref[1]))), 4),
            ankle_r_deg=round(float(q[s, iw["ankle_angle_r"]].mean()
                                    - ank_r0), 2),
            ankle_l_deg=round(float(q[s, iw["ankle_angle_l"]].mean()
                                    - ank_l0), 2),
            feet_pctBW=round(100.0 * 0.5 * float(np.mean(
                contact[s, 0] + contact[s, 1])) / BW, 1),
            act={nm: round(float(ctrl[s, ACT_NAMES.index(nm)].mean()), 3)
                 for nm in WATCH},
            gmean={g: round(float(ctrl[s, idx].mean()), 3)
                   for g, idx in GROUPS.items()},
            act_sum=round(float(np.mean(np.sum(ctrl[s], axis=1))), 1),
        ))
    base_act = {nm: float(ctrl[base, ACT_NAMES.index(nm)].mean())
                for nm in WATCH}
    return dict(
        tag=tag, cfg=CFG[tag], wall_s=round(wall, 1),
        fell=bool(m.get("bal_fell", False)),
        t_end=float(m.get("t_end", t[-1])),
        com_z_min=round(float(np.min(com[:, 2])), 3),
        base_act={k: round(v, 3) for k, v in base_act.items()},
        base_com_z=round(float(com_ref[2]), 3),
        load_steps=steps,
    )


def main():
    tags = list("ABCDE S".replace(" ", "")) if MODE == "ALL" else [MODE]
    results = {}
    out_path = "_reflex_stand_results_20261021.json"
    if os.path.exists(out_path):
        results = json.loads(open(out_path, encoding="utf-8").read())
    for tag in tags:
        print(f"\n===== reflex-stand config {tag}: {CFG[tag]} =====",
              flush=True)
        m, npz, wall = run_cfg(tag)
        a = analyze(tag, npz, m, wall)
        results[tag] = a
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=1, default=str)
        s0 = a["load_steps"][0]
        print(f"[{tag}] fell={int(a['fell'])} comz_min={a['com_z_min']} "
              f"base_comz={a['base_com_z']}", flush=True)
        for st in a["load_steps"]:
            d_soleus = st["act"]["soleus_r"] - a["base_act"]["soleus_r"]
            d_vas = st["act"]["vas_lat_r"] - a["base_act"]["vas_lat_r"]
            print(f"  {st['pctBW']:4.0f}%BW: com_dz={st['com_dz']:+.3f} "
                  f"sway={st['sway_max']:.3f} ank_r={st['ankle_r_deg']:+.1f} "
                  f"feet={st['feet_pctBW']}% d_soleus={d_soleus:+.3f} "
                  f"d_vaslat={d_vas:+.3f} actsum={st['act_sum']}", flush=True)
    print(f"\nwrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
