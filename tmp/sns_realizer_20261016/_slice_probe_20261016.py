"""Slice probe for the SNS_Simscape winner realization (2026-10-16 night).

Loads the round8_best3.json winner EXACTLY the way wean_tune_20261016.py
cmd_one does (apply_config: set_stage(4, base+fb) + knobs*scale), builds the
full human spinal network, then

  1. extracts the SLICE subgraph (right-leg RG/PF/heel-toe/KINH/IBEXC/LBIN
     + one knee MN pair with Renshaw + the left RG pair for the L/R
     gamma/alpha drive demonstration + C1 commissurals) with all EFFECTIVE
     parameters (units: mV, nA, uS, nF - SNS toolbox / SNS_Library
     conventions, same mapping export_network_json.py uses);
  2. runs the FULL compiled network (the real SNS_NumpyFixedTau backend,
     dt = 2 ms) with ONLY the slice's input currents nonzero (all afferent
     / balance / load ports zero) -> a numpy REFERENCE trace for the slice
     (RG dynamics are topologically closed under these inputs: no other
     cell projects into RG/PF);
  3. repeats the reference at gamma_r x2 and x0.5 (drive-knob modulation
     reference) and reads the production round8_full3.npz cadence for
     context;
  4. writes reports_20261016/sns_simulink/{slice_export.json,
     slice_ref.npz} + prints a summary.

Run ON egg, cwd = Code\\MuJoCo_SNS\\spinal:
  D:\\Anaconda\\envs\\myo\\python.exe _slice_probe_20261016.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reports_20261016", "sns_simulink")
os.makedirs(OUT, exist_ok=True)
WINNER = os.path.join(HERE, "reports_20261016", "wean_human",
                      "round8_best3.json")
MODEL = os.path.join(HERE, "..", "..", "..", "Solid_Models", "OpenSim",
                     "Gait2392_Robotbody", "mjc", "gait2392_simbody",
                     "gait2392_simbody_cvt3.xml")

# ---------------------------------------------------------------- load winner
import wean_tune_20261016 as WT  # noqa: E402  (cwd = spinal)

w = json.loads(open(WINNER, encoding="utf-8").read())
knobs = {k: float(w["knobs"][k]) for k in WT.DRIVE_KEYS}
fb = {k: float(w["params"].get(k, 0.0)) for k in WT.FB_KEYS}
eff = WT.apply_config(knobs, scale=float(w["scale"]), fb=fb)
print("APPLIED winner knobs ->", json.dumps(eff), flush=True)

import params as P  # noqa: E402
import runner as R  # noqa: E402
import build_network as bn  # noqa: E402
import muscle_map  # noqa: E402

print("EFFECTIVE: descend_to_rg_e=%.6f f=%.6f rg_to_pf=%.6f "
      "rg_mutual_inh=%.6f pf_recip_inh=%.6f pf_to_mn=%.6f "
      "posture_to_mn=%.6f posture_to_rg_e=%.6f rg_nap_h=%.6f "
      "c1_gain=%.6f heel_rge=%.6f toe_rge=%.6f ib_rge=%.6f "
      "ib_e_central=%.6f renshaw=%.6f f1_kneext_inh=%.6f "
      "ib_exc_to_mn=%.6f ib_group_exc=%.6f"
      % (P.G["descend_to_rg_e"], P.G["descend_to_rg_f"], P.G["rg_to_pf"],
         P.G["rg_mutual_inh"], P.G["pf_recip_inh"], P.G["pf_to_mn"],
         P.G["posture_to_mn"], P.G["posture_to_rg_e"], P.TAU["rg_nap_h"],
         P.G["c1_gain"], P.G["heel_rge"], P.G["toe_rge"], P.G["ib_rge"],
         P.G["ib_e_central"], P.G["renshaw"], P.G["f1_kneext_inh"],
         P.G["ib_exc_to_mn"], P.G["ib_group_exc"]), flush=True)

# ---------------------------------------------------------------- build net
import mujoco  # noqa: E402

m = mujoco.MjModel.from_xml_path(MODEL)
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        for i in range(m.nu)]
net = bn.build(acts, dt=P.DT, interleg=True)
n = net.net
print("BUILT: %d neurons, %d synapses, %d inputs" %
      (len(n.populations), len(n.connections), len(n.inputs)), flush=True)

# ---------------------------------------------------------------- slice set
def pick_muscle(group, prefer):
    for a in prefer:
        if a in net.mn_names:
            mi = muscle_map.classify(a)
            if mi and mi.side == "r" and mi.groups[0] == group:
                return a
    for a in net.mn_names:
        mi = muscle_map.classify(a)
        if mi and mi.side == "r" and mi.groups[0] == group:
            return a
    raise SystemExit("no %s right muscle found" % group)

MN_EXT = pick_muscle("knee_ext", ["vas_lat_r", "vas_med_r", "vas_int_r"])
MN_FLX = pick_muscle("knee_flex", ["bifemsh_r", "semiten_r", "semimem_r"])
print("MN pair: ext=%s flex=%s" % (MN_EXT, MN_FLX), flush=True)

SLICE = set(
    ["DRIVE", "POSTURE", "GAMMA_r", "GAMMA_l", "ALPHA_r", "ALPHA_l"]
    + ["RG_E_r", "RG_F_r", "InE_r", "InF_r",
       "RG_E_l", "RG_F_l", "InE_l", "InF_l",
       "CIN_F_r", "CIN_F_l"]
    + ["PF_E1_r", "PF_E2_r", "PF_F1_r", "PF_F2_r",
       "PF_IN_E_r", "PF_IN_F_r"]
    + ["HEEL_r", "TOE_r", "LBIN_r", "KINH_r",
       "HEEL_l", "TOE_l", "LBIN_l",
       "IBEXC_knee_ext_r",
       "MN_%s" % MN_EXT, "MN_%s" % MN_FLX,
       "RC_%s" % MN_EXT, "RC_%s" % MN_FLX]
)
by_name = {p["name"]: p for p in n.populations}

# ---------------------------------------------------------------- dump
from sns_toolbox.neurons import NonSpikingNeuron  # noqa: E402
from sns_toolbox.neurons import \
    NonSpikingNeuronWithPersistentSodiumChannel  # noqa: E402

neurons = []
for nm in sorted(SLICE):
    p = by_name[nm]
    t = p["type"]
    pr = t.params
    rec = {"name": nm, "kind": "rc", "tau_s": None, "Vrest_mV": 0.0,
           "Gm_uS": 1.0, "Cm_nF": None, "bias_nA": 0.0}
    rec["Vrest_mV"] = float(pr["resting_potential"])
    rec["Gm_uS"] = float(pr["membrane_conductance"])
    rec["Cm_nF"] = 1000.0 * float(pr["membrane_capacitance"])
    rec["tau_s"] = float(pr["membrane_capacitance"]) \
        / float(pr["membrane_conductance"])
    rec["bias_nA"] = float(pr.get("bias", 0.0))
    if isinstance(t, NonSpikingNeuronWithPersistentSodiumChannel):
        rec["kind"] = "nap"
        pa, pb = pr["paramsA"], pr["paramsB"]
        rec["nap"] = {
            "g_ion": float(np.atleast_1d(pr["Gion"])[0]),
            "e_ion": float(np.atleast_1d(pr["Eion"])[0]),
            "k_m": float(np.atleast_1d(pa["k"])[0]),
            "slope_m": float(np.atleast_1d(pa["slope"])[0]),
            "e_m": float(np.atleast_1d(pa["reversal"])[0]),
            "k_h": float(np.atleast_1d(pb["k"])[0]),
            "slope_h": float(np.atleast_1d(pb["slope"])[0]),
            "e_h": float(np.atleast_1d(pb["reversal"])[0]),
            "tau_h_s": float(np.atleast_1d(pb["TauMax"])[0]),
        }
    elif not isinstance(t, NonSpikingNeuron):
        raise SystemExit("unexpected neuron type for %s: %r" % (nm, t))
    neurons.append(rec)

synapses = []
for c in n.connections:
    src = n.populations[c["source"]]["name"]
    dst = n.populations[c["destination"]]["name"]
    if src in SLICE and dst in SLICE:
        cp = dict(c["type"].params)
        cp.update(c["params"])
        synapses.append({
            "src": src, "dst": dst,
            "g_uS": float(cp["max_conductance"]),
            "Esyn_mV": float(cp["reversal_potential"]),
            "ThrPre_mV": float(cp.get("e_lo", 0.0)),
            "SlopePre_mV": float(cp.get("e_hi", 5.0))
                           - float(cp.get("e_lo", 0.0)),
        })
print("SLICE: %d neurons, %d synapses" % (len(neurons), len(synapses)),
      flush=True)

inputs = []
for i in range(len(n.inputs)):
    dst = n.populations[n.inputs[i]["destination"]]["name"]
    if dst in SLICE:
        inputs.append({"port": net.inputs[i], "dst": dst})
print("SLICE input ports:", [r["port"] for r in inputs], flush=True)

# initial h of the RG gates: read the compiled backend BEFORE any step
bg = np.asarray(net.compiled.b_gate)
rg_col = net.idx["RG_E_r"]
h0 = float(np.atleast_1d(bg.ravel()[rg_col])
           if bg.size > rg_col else bg.ravel()[0])
print("compiled b_gate shape=%s RG_E_r h0=%.6f (all-cols first 8: %s)"
      % (bg.shape, h0, np.array2string(bg.ravel()[:8])), flush=True)
# snapshot the compiled initial state for exact re-runs
STATE0 = dict(V=np.array(net.compiled.V, copy=True),
              b_gate=np.array(bg, copy=True),
              c_gate=np.array(net.compiled.c_gate, copy=True))

export = {
    "meta": {
        "source": "round8_best3.json (curr_s8b_wean04 trial 35, "
                  "drive_scale %.16g)" % w["scale"],
        "loader": "wean_tune_20261016.apply_config(knobs, scale, fb) "
                  "== cmd_one winner path",
        "dt_s": P.DT,
        "E_HI_mV": 5.0,
        "mn_pair": {"ext": MN_EXT, "flx": MN_FLX},
        "units": "mV, nA, uS, nF (SNS toolbox / SNS_Library conventions)",
        "h0": h0,
        "schedule": {
            "walk_drive_nA": WT.BASE_DRIVE,
            "segments": {k: list(v) for k, v in P.SCHEDULE.items()},
            "posture": "1.0 stand / 0.6 walk (runner.drive_posture)",
            "knob_G_values": {
                "gamma_l": float(P.G["gamma_l"]),
                "gamma_r": float(P.G["gamma_r"]),
                "alpha_l": float(P.G["alpha_l"]),
                "alpha_r": float(P.G["alpha_r"]),
                "alpha_to_mn": float(P.G["alpha_to_mn"]),
            },
            "port_semantics": "u[GAMMA_c_s] = drive(t)*G[gamma_s], "
                              "u[ALPHA_c_s] = drive(t)*G[alpha_s] "
                              "(runner.py:1372-1378)",
        },
        "notes": [
            "Afferent encoders (Ia/II/Ib), BAL_*/VEST cells, pm_* runner "
            "phase machine, contra_* crossed edges (no_cross=1) and the "
            "muJoCo body are OUT of the slice; their input ports are "
            "simply not driven.",
            "POST_<muscle> MN bias ports default 0: production values "
            "come from the runner's MuJoCo standing solve (6.0*act_stand"
            "*stand_frac, ~0.05 floor in walk) - runner machinery, not "
            "circuit; they cannot affect the RG/PF rhythm (MNs project "
            "only to their Renshaw cells).",
        ],
    },
    "neurons": neurons,
    "synapses": synapses,
    "inputs": inputs,
}

# ------------------------------------------------------- numpy references
DT = P.DT
T_END = 22.0
NSTEPS = int(round(T_END / DT))
iport = {nm: i for i, nm in enumerate(net.inputs)}

WATCH = ["RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l",
         "PF_E1_r", "PF_E2_r", "PF_F1_r", "PF_F2_r",
         "MN_%s" % MN_EXT, "MN_%s" % MN_FLX,
         "KINH_r", "IBEXC_knee_ext_r", "LBIN_r", "HEEL_r", "TOE_r"]
wcol = {nm: net.idx[nm] for nm in WATCH}


def run_reference(gamma_r_mult=1.0):
    """Full-net run, only slice-modeled inputs nonzero; returns traces."""
    # exact re-init: restore the compiled backend's initial state
    net.compiled.V = np.array(STATE0["V"], copy=True)
    net.compiled.b_gate = np.array(STATE0["b_gate"], copy=True)
    net.compiled.c_gate = np.array(STATE0["c_gate"], copy=True)
    t_arr = np.zeros(NSTEPS + 1)
    drv = np.zeros(NSTEPS + 1)
    pos = np.zeros(NSTEPS + 1)
    tr = {nm: np.zeros(NSTEPS + 1) for nm in WATCH}
    u = np.zeros(len(net.inputs))
    for k in range(NSTEPS + 1):
        t = k * DT
        drive, posture = R.drive_posture(t, WT.BASE_DRIVE)
        u[:] = 0.0
        u[iport["DRIVE"]] = drive
        u[iport["POSTURE"]] = posture
        for s2 in ("r", "l"):
            gp, ap = "GAMMA_c_" + s2, "ALPHA_c_" + s2
            if gp in iport:
                u[iport[gp]] = drive * P.G["gamma_" + s2] * (
                    gamma_r_mult if s2 == "r" else 1.0)
            if ap in iport:
                u[iport[ap]] = drive * P.G["alpha_" + s2]
        t_arr[k], drv[k], pos[k] = t, drive, posture
        v = net.compiled.V.copy()
        for nm in WATCH:
            tr[nm][k] = v[wcol[nm]]
        if k < NSTEPS:
            net.compiled.forward(list(u))
    return t_arr, drv, pos, tr


def burst_stats(t, v, t0=5.0, t1=15.0):
    msk = (t >= t0) & (t <= t1)
    vv = v[msk]
    tt = t[msk]
    thr = 0.5 * (vv.min() + vv.max())
    on = vv > thr
    rises = np.where(np.diff(on.astype(int)) == 1)[0]
    if len(rises) < 2:
        return dict(cycles=len(rises), period_s=None, duty=None,
                    vmax=float(vv.max()), vmin=float(vv.min()))
    per = float(np.mean(np.diff(tt[rises])))
    return dict(cycles=int(len(rises)), period_s=per,
                duty=float(on.mean()), vmax=float(vv.max()),
                vmin=float(vv.min()))


refs = {}
stats = {}
for tag, mult in (("baseline", 1.0), ("gamma_r_x2", 2.0),
                  ("gamma_r_half", 0.5)):
    t_arr, drv, pos, tr = run_reference(mult)
    refs[tag] = (t_arr, drv, pos, tr)
    stats[tag] = {"RG_E_r": burst_stats(t_arr, tr["RG_E_r"]),
                  "RG_F_r": burst_stats(t_arr, tr["RG_F_r"]),
                  "gamma_r_mult": mult}
    s = stats[tag]["RG_E_r"]
    print("REF %-12s RG_E_r cycles=%s period=%s duty=%.3f "
          "v=[%.3f,%.3f] mV"
          % (tag, s["cycles"],
             ("%.4f" % s["period_s"]) if s["period_s"] else "None",
             float(s["duty"] or 0.0), float(s["vmin"]),
             float(s["vmax"])), flush=True)

# production npz for context (real gait cadence WITH afferents)
prod = {}
pnpz = os.path.join(HERE, "reports_20261016", "wean_human",
                    "round8_full3.npz")
if os.path.exists(pnpz):
    z = np.load(pnpz, allow_pickle=True)
    names = list(z["neuro_names"])
    if "RG_E_r" in names:
        ci = names.index("RG_E_r")
        prod = burst_stats(z["t"], z["neuro"][:, ci])
        print("PRODUCTION round8_full3 RG_E_r cycles=%s period=%s "
              "duty=%.3f" % (prod.get("cycles"),
                             ("%.4f" % prod["period_s"])
                             if prod.get("period_s") else "None",
                             prod.get("duty") or 0.0), flush=True)
else:
    print("PRODUCTION npz not found:", pnpz, flush=True)

export["meta"]["reference_stats"] = {"slice_input_runs": stats,
                                     "production_npz": prod}

with open(os.path.join(OUT, "slice_export.json"), "w",
          encoding="utf-8") as f:
    json.dump(export, f, indent=1)
print("WROTE slice_export.json", flush=True)

np.savez_compressed(
    os.path.join(OUT, "slice_ref.npz"),
    t=refs["baseline"][0], drive=refs["baseline"][1],
    posture=refs["baseline"][2],
    **{("base_" + nm): refs["baseline"][3][nm] for nm in WATCH},
    **{("g2_" + nm): refs["gamma_r_x2"][3][nm] for nm in WATCH},
    **{("gh_" + nm): refs["gamma_r_half"][3][nm] for nm in WATCH},
    watch=np.array(WATCH),
    stats=json.dumps(stats),
)
print("WROTE slice_ref.npz", flush=True)
print("PROBE_DONE", flush=True)
