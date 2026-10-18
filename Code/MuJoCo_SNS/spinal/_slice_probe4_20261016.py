"""Slice probe v4: SLICE-EQUIVALENT numpy references (AFF.i0_ii = 0).

Probes 2/3 ran the full winner network with all afferent ports zero - but
the afferent ENCODER neurons still carry their resting baseline tone
(AFF.i0_ii = 1.0 nA -> II cells sit at ~1 mV -> Sat = 0.2 excitatory
leakage into MNs via ii_to_mn and into PF_E/RG_E/InE via ii_e_central).
The Simulink slice has NO afferent encoders, so its exact reference is the
full net with the II baseline zeroed (II cells at 0 mV -> Sat = 0 ->
conductance exactly 0 = absent).

Re-runs the full verification set under i0_ii = 0 and prints the targets
the MATLAB checks should use:
  base / gamma_r x4 / alpha_r x0.5 / alpha_r x2 / contact cycle 2.1

Run on egg, cwd = Code\\MuJoCo_SNS\\spinal.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WINNER = os.path.join(HERE, "reports_20261016", "wean_human",
                      "round8_best3.json")
MODEL = os.path.join(HERE, "..", "..", "..", "Solid_Models", "OpenSim",
                     "Gait2392_Robotbody", "mjc", "gait2392_simbody",
                     "gait2392_simbody_cvt3.xml")

import wean_tune_20261016 as WT  # noqa: E402

w = json.loads(open(WINNER, encoding="utf-8").read())
knobs = {k: float(w["knobs"][k]) for k in WT.DRIVE_KEYS}
fb = {k: float(w["params"].get(k, 0.0)) for k in WT.FB_KEYS}
WT.apply_config(knobs, scale=float(w["scale"]), fb=fb)

import params as P  # noqa: E402
import runner as R  # noqa: E402
import build_network as bn  # noqa: E402

P.AFF["i0_ii"] = 0.0        # <- the slice-equivalence switch
print("AFF.i0_ii ->", P.AFF["i0_ii"], flush=True)

import mujoco  # noqa: E402

m = mujoco.MjModel.from_xml_path(MODEL)
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        for i in range(m.nu)]
net = bn.build(acts, dt=P.DT, interleg=True)
iport = {nm: i for i, nm in enumerate(net.inputs)}
STATE0 = dict(V=np.array(net.compiled.V, copy=True),
              b_gate=np.array(net.compiled.b_gate, copy=True),
              c_gate=np.array(net.compiled.c_gate, copy=True))

DT = P.DT
NSTEPS = int(round(22.0 / DT))
CYCLE = 1.7284
T0 = 5.0
G = P.G
WATCH = ["RG_E_r", "RG_F_r", "MN_vas_lat_r", "MN_bifemsh_r"]
wcol = {nm: net.idx[nm] for nm in WATCH}


def foot_forces(phi):
    fh = np.zeros_like(phi)
    ft = np.zeros_like(phi)
    fh = np.where(phi < 0.03, 0.95 * phi / 0.03, fh)
    fh = np.where((phi >= 0.03) & (phi < 0.30),
                  0.95 - 0.60 * (phi - 0.03) / 0.27, fh)
    fh = np.where((phi >= 0.30) & (phi < 0.55),
                  0.35 - 0.35 * (phi - 0.30) / 0.25, fh)
    ft = np.where((phi >= 0.30) & (phi < 0.55),
                  0.55 * (phi - 0.30) / 0.25, ft)
    ft = np.where((phi >= 0.55) & (phi < 0.62), 0.55, ft)
    ft = np.where((phi >= 0.62) & (phi < 0.70),
                  0.55 + 0.40 * (phi - 0.62) / 0.08, ft)
    ft = np.where((phi >= 0.70) & (phi < 0.72),
                  0.95 - 0.95 * (phi - 0.70) / 0.02, ft)
    return np.clip(fh, 0, 1.5), np.clip(ft, 0, 1.5)


def run(gamma_r_mult=1.0, alpha_r_mult=1.0, cycle=CYCLE):
    net.compiled.V = np.array(STATE0["V"], copy=True)
    net.compiled.b_gate = np.array(STATE0["b_gate"], copy=True)
    net.compiled.c_gate = np.array(STATE0["c_gate"], copy=True)
    t_arr = np.zeros(NSTEPS + 1)
    tr = {nm: np.zeros(NSTEPS + 1) for nm in WATCH}
    kick_state = {"r": 0.0, "l": 0.0}
    prev_loaded = {"r": False, "l": False}
    decay = float(np.exp(-DT / 0.25))
    u = np.zeros(len(net.inputs))
    for k in range(NSTEPS + 1):
        t = k * DT
        drive, posture = R.drive_posture(t, WT.BASE_DRIVE)
        u[:] = 0.0
        u[iport["DRIVE"]] = drive
        u[iport["POSTURE"]] = posture
        for s2 in ("r", "l"):
            u[iport["GAMMA_c_" + s2]] = drive * G["gamma_" + s2] * (
                gamma_r_mult if s2 == "r" else 1.0)
            u[iport["ALPHA_c_" + s2]] = drive * G["alpha_" + s2] * (
                alpha_r_mult if s2 == "r" else 1.0)
        if t >= T0:
            for s2, lag in (("r", 0.0), ("l", 0.5)):
                phi = ((t - T0) / cycle - lag) % 1.0
                fh, ft = foot_forces(np.array(phi))
                heel_sig = min(float(fh) / 0.35, 1.5)
                toe_sig = min(float(ft) / 0.50, 1.5)
                load_sig = min((float(fh) + float(ft)) / 0.60, 1.5)
                loaded = load_sig > 0.05
                if loaded and not prev_loaded[s2]:
                    kick_state[s2] = 1.0
                elif not loaded and prev_loaded[s2]:
                    kick_state[s2] = -1.0
                else:
                    kick_state[s2] *= decay
                prev_loaded[s2] = loaded
                kick = G["contact_onset"] * kick_state[s2]
                u[iport["HEEL_c_" + s2]] = G["heel_rge"] * heel_sig + kick
                u[iport["TOE_c_" + s2]] = G["toe_rge"] * toe_sig + kick
                u[iport["LOAD_c_" + s2]] = G["ib_rge"] * load_sig
        t_arr[k] = t
        v = net.compiled.V
        for nm in WATCH:
            tr[nm][k] = v[wcol[nm]]
        if k < NSTEPS:
            net.compiled.forward(list(u))
    return t_arr, tr


def bstats(t, v, t0=6.0, t1=15.0):
    msk = (t >= t0) & (t <= t1)
    vv, tt = v[msk], t[msk]
    thr = 0.5 * (vv.min() + vv.max())
    on = vv > thr
    rises = np.where(np.diff(on.astype(int)) == 1)[0]
    if len(rises) < 2:
        return dict(cycles=len(rises), period=None, duty=None)
    return dict(cycles=int(len(rises)),
                period=round(float(np.mean(np.diff(tt[rises]))), 4),
                duty=round(float(on.mean()), 3),
                vmin=round(float(vv.min()), 3),
                vmax=round(float(vv.max()), 3))


for tag, kw in (("base", {}), ("gamma_x4", dict(gamma_r_mult=4.0)),
                ("alpha_x05", dict(alpha_r_mult=0.5)),
                ("alpha_x2", dict(alpha_r_mult=2.0)),
                ("cycle_21", dict(cycle=2.1))):
    t_arr, tr = run(**kw)
    se = bstats(t_arr, tr["RG_E_r"])
    sf = bstats(t_arr, tr["RG_F_r"])
    msk = (t_arr >= 6.0) & (t_arr <= 15.0)
    mne = float(np.mean(tr["MN_vas_lat_r"][msk]))
    mnf = float(np.mean(tr["MN_bifemsh_r"][msk]))
    print("P4 %-9s RG_E cyc=%s per=%s duty=%s v=[%s,%s] | RG_F duty=%s "
          "v=[%s,%s] | MNvas=%.3f MNbif=%.3f"
          % (tag, se["cycles"], se["period"], se["duty"], se["vmin"],
             se["vmax"], sf["duty"], sf["vmin"], sf["vmax"], mne, mnf),
          flush=True)
print("PROBE4_DONE", flush=True)
