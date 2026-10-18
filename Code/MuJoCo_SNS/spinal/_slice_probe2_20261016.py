"""Slice probe v2 (2026-10-16 night): SCRIPTED-CONTACT numpy references.

Probe v1 established that with ZERO afferents the winner RG latches on the
NaP plateau (V ~ 5.1 mV, no bursts) - the s3k-lineage rhythm is
afferent-sustained. This probe drives the full winner network with a
STYLIZED bilateral antiphase foot-contact script at the production cadence
(round8_full3.npz: RG_E_r period 1.7284 s, duty ~0.71), shaped exactly like
the runner feeds it:

    heel_sig = min(Fh/0.35BW, 1.5); toe_sig = min(Ft/0.50BW, 1.5);
    load_sig = min((Fh+Ft)/0.60BW, 1.5)
    u[HEEL_c_s] = heel_rge*heel_sig + contact_onset*kick_s   (same TOE_c)
    u[LOAD_c_s] = ib_rge*load_sig
    kick: +1 at loading edge, -1 at unloading edge, x exp(-dt/0.25)

Outputs: reports_20261016/sns_simulink/slice_ref2.npz + prints.
Run on egg, cwd = Code\\MuJoCo_SNS\\spinal.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reports_20261016", "sns_simulink")
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
T_END = 22.0
NSTEPS = int(round(T_END / DT))
CYCLE = 1.7284          # production RG_E_r period (round8_full3.npz)
T0 = 5.0                # contacts start when drive reaches full walk value
G = P.G

WATCH = ["RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l",
         "PF_E1_r", "PF_E2_r", "PF_F1_r", "PF_F2_r",
         "MN_vas_lat_r", "MN_bifemsh_r",
         "KINH_r", "IBEXC_knee_ext_r", "LBIN_r", "HEEL_r", "TOE_r"]
wcol = {nm: net.idx[nm] for nm in WATCH}


def foot_forces(phi):
    """Stylized stance foot forces, fractions of BW. phi in [0,1), phi=0
    = heel strike; any-contact duty ~0.72 (production ~0.71)."""
    fh = np.zeros_like(phi)
    ft = np.zeros_like(phi)
    # heel: impact spike then decay through foot-flat, off by phi 0.55
    fh = np.where(phi < 0.03, 0.95 * phi / 0.03, fh)
    fh = np.where((phi >= 0.03) & (phi < 0.30),
                  0.95 - 0.60 * (phi - 0.03) / 0.27, fh)
    fh = np.where((phi >= 0.30) & (phi < 0.55),
                  0.35 - 0.35 * (phi - 0.30) / 0.25, fh)
    # toe: loads mid-stance, push-off peak, off by phi 0.72 (toe-off)
    ft = np.where((phi >= 0.30) & (phi < 0.55),
                  0.55 * (phi - 0.30) / 0.25, ft)
    ft = np.where((phi >= 0.55) & (phi < 0.62), 0.55, ft)
    ft = np.where((phi >= 0.62) & (phi < 0.70),
                  0.55 + 0.40 * (phi - 0.62) / 0.08, ft)
    ft = np.where((phi >= 0.70) & (phi < 0.72),
                  0.95 - 0.95 * (phi - 0.70) / 0.02, ft)
    ft = np.clip(ft, 0.0, 1.5)
    fh = np.clip(fh, 0.0, 1.5)
    return fh, ft


def run_scripted(gamma_r_mult=1.0, cycle=CYCLE, heel_scale=1.0):
    net.compiled.V = np.array(STATE0["V"], copy=True)
    net.compiled.b_gate = np.array(STATE0["b_gate"], copy=True)
    net.compiled.c_gate = np.array(STATE0["c_gate"], copy=True)
    t_arr = np.zeros(NSTEPS + 1)
    drv = np.zeros(NSTEPS + 1)
    tr = {nm: np.zeros(NSTEPS + 1) for nm in WATCH}
    heel_arr = np.zeros(NSTEPS + 1)
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
            u[iport["ALPHA_c_" + s2]] = drive * G["alpha_" + s2]
        # scripted contacts (both feet, antiphase), only during the drive
        if t >= T0:
            for s2, lag in (("r", 0.0), ("l", 0.5)):
                phi = ((t - T0) / cycle - lag) % 1.0
                fh, ft = foot_forces(np.array(phi))
                heel_sig = min(float(fh) / 0.35, 1.5) * heel_scale
                toe_sig = min(float(ft) / 0.50, 1.5) * heel_scale
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
                if s2 == "r":
                    heel_arr[k] = u[iport["HEEL_c_r"]]
        t_arr[k], drv[k] = t, drive
        v = net.compiled.V
        for nm in WATCH:
            tr[nm][k] = v[wcol[nm]]
        if k < NSTEPS:
            net.compiled.forward(list(u))
    return t_arr, drv, tr, heel_arr


def burst_stats(t, v, t0=6.0, t1=15.0):
    msk = (t >= t0) & (t <= t1)
    vv, tt = v[msk], t[msk]
    thr = 0.5 * (vv.min() + vv.max())
    on = vv > thr
    rises = np.where(np.diff(on.astype(int)) == 1)[0]
    if len(rises) < 2:
        return dict(cycles=len(rises), period_s=None, duty=None,
                    vmax=float(vv.max()), vmin=float(vv.min()))
    per = float(np.mean(np.diff(tt[rises])))
    return dict(cycles=int(len(rises)), period_s=per, duty=float(on.mean()),
                vmax=float(vv.max()), vmin=float(vv.min()))


results = {}
for tag, kw in (("scripted", {}),
                ("scripted_g2", dict(gamma_r_mult=2.0)),
                ("scripted_gh", dict(gamma_r_mult=0.5))):
    t_arr, drv, tr, heel_arr = run_scripted(**kw)
    results[tag] = (t_arr, drv, tr, heel_arr)
    for nm in ("RG_E_r", "RG_F_r"):
        s = burst_stats(t_arr, tr[nm])
        print("SCRIPTED %-12s %-7s cycles=%s period=%s duty=%s "
              "v=[%.3f,%.3f]"
              % (tag, nm, s["cycles"],
                 ("%.4f" % s["period_s"]) if s["period_s"] else "None",
                 ("%.3f" % s["duty"]) if s["duty"] is not None else "None",
                 s["vmin"], s["vmax"]), flush=True)

t_arr, drv, tr, heel_arr = results["scripted"]
np.savez_compressed(
    os.path.join(OUT, "slice_ref2.npz"),
    t=t_arr, drive=drv, heel_r=heel_arr,
    **{("sc_" + nm): tr[nm] for nm in WATCH},
    **{("g2_" + nm): results["scripted_g2"][2][nm] for nm in WATCH},
    **{("gh_" + nm): results["scripted_gh"][2][nm] for nm in WATCH},
    watch=np.array(WATCH),
    cycle=CYCLE, t0=T0,
)
print("WROTE slice_ref2.npz", flush=True)
print("PROBE2_DONE", flush=True)
