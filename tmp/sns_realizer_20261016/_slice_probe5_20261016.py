"""Slice probe v5: simulate THE SLICE ITSELF (34 neurons from
slice_export.json) with the production inputs, exact toolbox semantics.

Why: probes 2-4 simulate the FULL winner network (46 MN pools/side). The
Renshaw network effect differs between the full net (every RC inhibited by
~44 partner pools -> RC sits low -> MN_vas 3.089 mV) and the 2-pool slice
(RC only sees its pair partner -> more recurrent inhibition -> Simulink
MN_vas 2.247 mV). The slice is a DIFFERENT (smaller) circuit, so its
reference must be computed on the slice topology. This probe implements
the exact SNS_NumpyFixedTau math on the exported 34-neuron subgraph:

  Cm dV/dt = Gm(Vrest-V) + Iapp + sum_k g_k*Sat(Vpre_k)*(E_k - V)
             [+ g_ion*a_inf(V)*h*(e_ion-V)  for the RG NaP somata]
  Sat(Vpre) = clip((Vpre-0)/5, 0, 1);  dh/dt = (b_inf(V)-h)/tau_h (FIXED)
  Euler dt = 2 ms, V(0) = 0, h(0) = h0 from the compiled winner backend

The RG/PF part is cross-anchored to the REAL backend (probe2/3 full-net
stats, which the Simulink model reproduced near-exactly), so the stepper
is validated where it can be; the MN numbers are then the slice-exact
CHECK-C targets.

Run on egg, cwd = Code\\MuJoCo_SNS\\spinal.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reports_20261016", "sns_simulink")

J = json.loads(open(os.path.join(OUT, "slice_export.json"),
                    encoding="utf-8").read())
sched = J["meta"]["schedule"]
kn = sched["knob_G_values"]
h0 = J["meta"]["h0"]
DT = 0.002
NSTEPS = int(round(22.0 / DT))
CYCLE = 1.7284
T0 = 5.0
HEEL_RGE = 1.4827569817170443
TOE_RGE = 0.36000460790520544
IB_RGE = 0.31221919098930606
ONSET = 0.8172727211325486

names = [n["name"] for n in J["neurons"]]
idx = {nm: i for i, nm in enumerate(names)}
N = len(names)
isnap = np.zeros(N, bool)
prm = dict(Cm=np.zeros(N), Gm=np.zeros(N), Vrest=np.zeros(N))
napg = np.zeros(N); nape = np.zeros(N)
mk = np.zeros(N); sm = np.zeros(N); em = np.zeros(N)
kh = np.zeros(N); sh = np.zeros(N); eh = np.zeros(N)
tauh = np.zeros(N)
for n in J["neurons"]:
    i = idx[n["name"]]
    prm["Cm"][i] = n["Cm_nF"]
    prm["Gm"][i] = n["Gm_uS"]
    prm["Vrest"][i] = n["Vrest_mV"]
    if n["kind"] == "nap":
        isnap[i] = True
        q = n["nap"]
        napg[i] = q["g_ion"]; nape[i] = q["e_ion"]
        mk[i] = q["k_m"]; sm[i] = q["slope_m"]; em[i] = q["e_m"]
        kh[i] = q["k_h"]; sh[i] = q["slope_h"]; eh[i] = q["e_h"]
        tauh[i] = q["tau_h_s"] * 1000.0   # ms
src = np.array([idx[s["src"]] for s in J["synapses"]])
dst = np.array([idx[s["dst"]] for s in J["synapses"]])
gsyn = np.array([s["g_uS"] for s in J["synapses"]])
esyn = np.array([s["Esyn_mV"] for s in J["synapses"]])
thr = np.array([s["ThrPre_mV"] for s in J["synapses"]])
slope = np.array([s["SlopePre_mV"] for s in J["synapses"]])
NS = len(src)


def drive_posture(t, wd):
    if t < 3.0:
        return 0.0, 1.0
    if t < 5.0:
        f = (t - 3.0) / 2.0
        return wd * f, 1.0 - 0.4 * f
    if t < 15.0:
        return wd, 0.6
    if t < 17.0:
        f = (t - 15.0) / 2.0
        return wd * (1.0 - f), 0.6 + 0.4 * f
    return 0.0, 1.0


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
    V = np.zeros(N)
    h = np.full(N, h0)
    t_arr = np.zeros(NSTEPS + 1)
    watch = ["RG_E_r", "RG_F_r", "MN_vas_lat_r", "MN_bifemsh_r"]
    tr = {nm: np.zeros(NSTEPS + 1) for nm in watch}
    kick_state = {"r": 0.0, "l": 0.0}
    prev_loaded = {"r": False, "l": False}
    decay = float(np.exp(-DT / 0.25))
    wd = sched["walk_drive_nA"]
    iapp = np.zeros(N)
    for k in range(NSTEPS + 1):
        t = k * DT
        drive, posture = drive_posture(t, wd)
        iapp[:] = 0.0
        iapp[idx["DRIVE"]] = drive
        iapp[idx["POSTURE"]] = posture
        for s2 in ("r", "l"):
            iapp[idx["GAMMA_" + s2]] = drive * kn["gamma_" + s2] * (
                gamma_r_mult if s2 == "r" else 1.0)
            iapp[idx["ALPHA_" + s2]] = drive * kn["alpha_" + s2] * (
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
                kickc = ONSET * kick_state[s2]
                iapp[idx["HEEL_" + s2]] = HEEL_RGE * heel_sig + kickc
                iapp[idx["TOE_" + s2]] = TOE_RGE * toe_sig + kickc
                iapp[idx["LBIN_" + s2]] = IB_RGE * load_sig
        t_arr[k] = t
        for nm in watch:
            tr[nm][k] = V[idx[nm]]
        if k == NSTEPS:
            break
        # synapse conductances from CURRENT V (Euler: V_last)
        sat = np.clip((V[src] - thr) / slope, 0.0, 1.0)
        g = gsyn * sat
        isyn = np.zeros(N)
        np.add.at(isyn, dst, g * esyn)
        np.add.at(isyn, dst, -V[dst] * g)
        dV = (-prm["Gm"] * (V - prm["Vrest"]) + iapp + isyn) * 1000.0 \
            / prm["Cm"]
        if isnap.any():
            ii = np.where(isnap)[0]
            a_inf = 1.0 / (1.0 + mk[ii] * np.exp(sm[ii] * (em[ii] - V[ii])))
            b_inf = 1.0 / (1.0 + kh[ii] * np.exp(sh[ii] * (eh[ii] - V[ii])))
            # backend order (SNS_NumpyFixedTau): update h FIRST (from
            # V_last), then i_ion uses the NEW h
            h[ii] = h[ii] + DT * 1000.0 * (b_inf - h[ii]) / tauh[ii]
            i_nap = np.zeros(N)
            i_nap[ii] = napg[ii] * a_inf * h[ii] * (nape[ii] - V[ii])
            dV = dV + i_nap * 1000.0 / prm["Cm"]
        # dV is in mV/s here (1000/Cm converts nA/nF -> mV/s); DT is in s
        V = V + DT * dV
    return t_arr, tr


def bstats(t, v, t0=6.0, t1=15.0):
    msk = (t >= t0) & (t <= t1)
    vv, tt = v[msk], t[msk]
    thr_ = 0.5 * (vv.min() + vv.max())
    on = vv > thr_
    rises = np.where(np.diff(on.astype(int)) == 1)[0]
    if len(rises) < 2:
        return dict(cycles=len(rises), period=None, duty=None,
                    vmin=round(float(vv.min()), 3),
                    vmax=round(float(vv.max()), 3))
    return dict(cycles=int(len(rises)),
                period=round(float(np.mean(np.diff(tt[rises]))), 4),
                duty=round(float(on.mean()), 3),
                vmin=round(float(vv.min()), 3),
                vmax=round(float(vv.max()), 3))


res = {}
for tag, kw in (("base", {}), ("gamma_x4", dict(gamma_r_mult=4.0)),
                ("alpha_x05", dict(alpha_r_mult=0.5)),
                ("alpha_x2", dict(alpha_r_mult=2.0)),
                ("cycle_21", dict(cycle=2.1))):
    t_arr, tr = run(**kw)
    se = bstats(t_arr, tr["RG_E_r"])
    sf = bstats(t_arr, tr["RG_F_r"])
    msk = (t_arr >= 6.0) & (t_arr <= 15.0)
    mne = round(float(np.mean(tr["MN_vas_lat_r"][msk])), 3)
    mnf = round(float(np.mean(tr["MN_bifemsh_r"][msk])), 3)
    res[tag] = dict(RG_E=se, RG_F=sf, MN_vas=mne, MN_bif=mnf)
    print("P5 %-9s RG_E cyc=%s per=%s duty=%s v=[%s,%s] | RG_F duty=%s "
          "v=[%s,%s] | MNvas=%.3f MNbif=%.3f"
          % (tag, se["cycles"], se["period"], se["duty"], se["vmin"],
             se["vmax"], sf["duty"], sf["vmin"], sf["vmax"], mne, mnf),
          flush=True)

with open(os.path.join(OUT, "slice_ref_targets.json"), "w",
          encoding="utf-8") as f:
    json.dump(res, f, indent=1)
print("WROTE slice_ref_targets.json", flush=True)
print("PROBE5_DONE", flush=True)
