"""OpenSim 4.6 analysis: add_mag3_r moment arms + torque envelopes about
hip_flexion_r and hip_adduction_r for the two P2 variants, vs the human
reference (stock gait2392 route, identical Thelen parameters).

Human reference: Add_Mag_Mesh_Opt.m lines 93-99 (Adductor Magnus 3, MIF 488 N,
OFL 0.131, TSL 0.249, pennation 0.08726646 rad) on the stock gait2392 route
P1 pelvis (-0.11108,-0.11413,0.04882) -> P2 femur_r (0.007,-0.3837,-0.0266) =
the add_mag3_r path carried by ConnorBipedal.osim (grep-verified: identical
Thelen params in ConnorBipedal / robotbody / stock gait2392_simbody).

Two torque conventions, both reported:
  A) MIF x arm (project convention, MonoMuscleData/Add_Mag_Mesh_Opt):
     tau = |arm| * MIF(488 N). Pure path geometry, no passive-stretch
     inflation -- PRIMARY for Ben's performance ruling.
  B) Thelen2003 isometric equilibrium at activation=1 (tau = arm * F_tendon,
     incl. passive force) -- cross-check. add_mag3_r ONLY is equilibrated
     (all other muscles disabled in memory; full-model equilibrium fails on
     Ben's rerouted semimem_r, irrelevant here).

Two grids:
  PRIMARY = project anatomical hip RoM (Add_Mag_Mesh_Opt.m:41-52):
            flexion -25..+85 deg, adduction -45..+20 deg.
  FULL    = model clamped ranges (all three models +/-120 deg, i.e. includes
            anatomically impossible poses; passive forces saturate up to
            ~14 kN) -- kept as supplementary CSV only.

Run: D:/Anaconda/envs/opensim/python.exe add_mag3_torque_compare.py
"""

import csv
import json
import math
from pathlib import Path

import opensim as osim

HERE = Path(__file__).parent
MODELS = {
    "variant_repo": HERE / "variant_repo.osim",
    "variant_thumb": HERE / "variant_thumb.osim",
    "human_ref": HERE.parent / "ConnorBipedal.osim",
}
MUSCLE = "add_mag3_r"
DOFS = ("hip_flexion_r", "hip_adduction_r")
MIF = 488.0  # N, human Adductor Magnus 3 max isometric force, Add_Mag_Mesh_Opt.m:95

# PRIMARY grid = the project's hip RoM for this muscle (Add_Mag_Mesh_Opt.m:41-52)
PRIMARY_DEG = {"hip_flexion_r": (-25.0, 85.0, 23), "hip_adduction_r": (-45.0, 20.0, 14)}


def coord_range(model, name):
    c = model.getCoordinateSet().get(name)
    try:
        lo, hi = c.getRangeMin(), c.getRangeMax()
        if math.isnan(lo) or math.isnan(hi):
            raise ValueError("no range")
    except Exception:
        raise RuntimeError(f"coordinate {name} has no clamped range")
    return lo, hi


def pose_row(mu, mu_t, state, coords, fv, av):
    """One grid pose: moment arms + both torque conventions for add_mag3_r."""
    coords["hip_flexion_r"].setValue(state, fv)
    coords["hip_adduction_r"].setValue(state, av)
    mu.setActivation(state, 1.0)
    try:
        model_eq_ok = True
        mu_t.computeEquilibrium(state)
    except Exception:
        model_eq_ok = False
    arm_f = mu.computeMomentArm(state, coords["hip_flexion_r"])
    arm_a = mu.computeMomentArm(state, coords["hip_adduction_r"])
    row = {
        "hip_flexion_deg": round(math.degrees(fv), 4),
        "hip_adduction_deg": round(math.degrees(av), 4),
        "l_mt_m": mu.getLength(state),
        "arm_hip_flexion_m": arm_f,
        "arm_hip_adduction_m": arm_a,
        "tauA_hip_flexion_Nm": abs(arm_f) * MIF,
        "tauA_hip_adduction_Nm": abs(arm_a) * MIF,
    }
    tf = None
    if model_eq_ok:
        try:
            tf = mu.getTendonForce(state)
        except Exception:
            tf = None
    row["tendon_force_N"] = tf if tf is not None else ""
    row["tauB_hip_flexion_Nm"] = arm_f * tf if tf is not None else ""
    row["tauB_hip_adduction_Nm"] = arm_a * tf if tf is not None else ""
    row["equilibrium_ok"] = tf is not None
    return row


def analyze(model_path):
    model = osim.Model(str(model_path))
    model.finalizeFromProperties()
    for m in model.getMuscles():          # isolate our muscle for equilibrium
        if m.getName() != MUSCLE:
            m.set_appliesForce(False)
    state = model.initSystem()
    mu = model.getMuscles().get(MUSCLE)
    mu_t = osim.Thelen2003Muscle.safeDownCast(mu)
    coords = {d: model.getCoordinateSet().get(d) for d in DOFS}
    clamped = {d: coord_range(model, d) for d in DOFS}

    # PRIMARY grid: project anatomical RoM on both DOFs
    fl, fh, nf = PRIMARY_DEG["hip_flexion_r"]
    al, ah, na = PRIMARY_DEG["hip_adduction_r"]
    f_lo = math.radians(fl)
    f_hi = math.radians(fh)
    a_lo = math.radians(al)
    a_hi = math.radians(ah)
    flex_vals = [f_lo + (f_hi - f_lo) * i / (nf - 1) for i in range(nf)]
    add_vals = [a_lo + (a_hi - a_lo) * i / (na - 1) for i in range(na)]
    primary = []
    for fv in flex_vals:
        for av in add_vals:
            primary.append(pose_row(mu, mu_t, state, coords, fv, av))

    # FULL grid: model clamped ranges, one DOF at a time (supplementary)
    full = []
    for d in DOFS:
        lo, hi = clamped[d]
        n = 25
        vals = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
        for v in vals:
            fv = v if d == "hip_flexion_r" else 0.0
            av = v if d == "hip_adduction_r" else 0.0
            full.append(pose_row(mu, mu_t, state, coords, fv, av))

    def peaks(rows):
        return {
            "peak_tauA_hip_flexion_Nm": max(r["tauA_hip_flexion_Nm"] for r in rows),
            "peak_tauA_hip_adduction_Nm": max(r["tauA_hip_adduction_Nm"] for r in rows),
            "peak_arm_hip_flexion_m": max(abs(r["arm_hip_flexion_m"]) for r in rows),
            "peak_arm_hip_adduction_m": max(abs(r["arm_hip_adduction_m"]) for r in rows),
            "peak_tauB_hip_flexion_Nm": max(abs(r["tauB_hip_flexion_Nm"]) for r in rows if r["equilibrium_ok"]),
            "peak_tauB_hip_adduction_Nm": max(abs(r["tauB_hip_adduction_Nm"]) for r in rows if r["equilibrium_ok"]),
            "peak_tendon_force_N": max(r["tendon_force_N"] for r in rows if r["equilibrium_ok"]),
            "mean_tauA_hip_flexion_Nm": sum(r["tauA_hip_flexion_Nm"] for r in rows) / len(rows),
            "mean_tauA_hip_adduction_Nm": sum(r["tauA_hip_adduction_Nm"] for r in rows) / len(rows),
        }

    summary = {
        "model": str(model_path),
        "clamped_range_deg": {d: [round(math.degrees(lo), 2), round(math.degrees(hi), 2)]
                              for d, (lo, hi) in clamped.items()},
        "primary_grid": {"hip_flexion_deg": [fl, fh, nf], "hip_adduction_deg": [al, ah, na]},
        "primary": peaks(primary),
        "n_primary_poses": len(primary),
        "n_primary_equilibrium_failed": sum(1 for r in primary if not r["equilibrium_ok"]),
    }
    return summary, primary, full


def main():
    out = {}
    store = {}
    for name, path in MODELS.items():
        summary, primary, full = analyze(path)
        store[name] = (primary, full)
        out[name] = summary
        p = summary["primary"]
        print(f"[{name}] clamped {summary['clamped_range_deg']}; primary grid "
              f"{summary['primary_grid']}; eq-fail {summary['n_primary_equilibrium_failed']}")
        print(f"[{name}] PRIMARY peak |tauA(add-flx envelope, MIF x arm)| "
              f"flex {p['peak_tauA_hip_flexion_Nm']:.2f} N*m, "
              f"add {p['peak_tauA_hip_adduction_Nm']:.2f} N*m | "
              f"peak |arm| flex {p['peak_arm_hip_flexion_m']*100:.2f} cm, "
              f"add {p['peak_arm_hip_adduction_m']*100:.2f} cm | "
              f"mean tauA flex {p['mean_tauA_hip_flexion_Nm']:.2f}, "
              f"add {p['mean_tauA_hip_adduction_Nm']:.2f} | "
              f"peak |tauB(eq, act=1)| flex {p['peak_tauB_hip_flexion_Nm']:.2f}, "
              f"add {p['peak_tauB_hip_adduction_Nm']:.2f}")

    def load_store(name, which):
        return store[name][0 if which == "primary" else 1]

    hu = load_store("human_ref", "primary")
    for name in ("variant_repo", "variant_thumb"):
        vr = load_store(name, "primary")
        for key, tag in (("tauA_hip_flexion_Nm", "flex"), ("tauA_hip_adduction_Nm", "add"),
                         ("tauB_hip_flexion_Nm", "flexB"), ("tauB_hip_adduction_Nm", "addB")):
            ratios = []
            for rv, rh in zip(vr, hu):
                th = abs(rh[key])
                if th > 1.0:
                    ratios.append((abs(rv[key]) / th, (rv["hip_flexion_deg"], rv["hip_adduction_deg"])))
            ratios.sort()
            if ratios:
                out[name][f"min_ratio_{tag}"] = round(ratios[0][0], 4)
                out[name][f"min_ratio_pose_{tag}"] = list(ratios[0][1])
                out[name][f"n_poses_below_human_{tag}"] = sum(1 for r, _ in ratios if r < 1.0)
                out[name][f"n_poses_compared_{tag}"] = len(ratios)
                print(f"[{name}] {tag}: min ratio {ratios[0][0]:.4f} at {ratios[0][1]}, "
                      f"{out[name][f'n_poses_below_human_{tag}']}/{len(ratios)} poses below human")

    with open(HERE / "add_mag3_torque_summary.json", "w") as f:
        json.dump(out, f, indent=2)
    for name in MODELS:
        for which in ("primary", "full"):
            rows = load_store(name, which)
            with open(HERE / f"arm_tau_{name}_{which}.csv", "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
    print("wrote add_mag3_torque_summary.json + arm_tau_*_{primary,full}.csv")

if __name__ == "__main__":
    main()
