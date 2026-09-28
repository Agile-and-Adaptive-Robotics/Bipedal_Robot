"""Validation + profile extraction for the add_mag3_r comparison.

1. Numerical cross-check of computeMomentArm vs central-difference dl/dq.
2. Profile tables on the primary grid at the conjugate DOF = 0:
   - variant_repo hip_flexion arm vs human (the flexion-null question)
   - variant_thumb hip_adduction arm vs human (the adduction-deficit question)
"""

import math
import csv
from pathlib import Path

import opensim as osim

HERE = Path(__file__).parent
EPS = 1e-5  # rad, central difference step
MIF = 488.0


def check_arms(model_path, label, flex_deg, add_deg):
    model = osim.Model(str(model_path))
    model.finalizeFromProperties()
    state = model.initSystem()
    mu = model.getMuscles().get("add_mag3_r")
    cf = model.getCoordinateSet().get("hip_flexion_r")
    ca = model.getCoordinateSet().get("hip_adduction_r")
    f, a = math.radians(flex_deg), math.radians(add_deg)

    def l_mt(fv, av):
        cf.setValue(state, fv)
        ca.setValue(state, av)
        return mu.getLength(state)

    for name, c, v, other in (("flexion", cf, f, a), ("adduction", ca, a, f)):
        if name == "flexion":
            num = (l_mt(f + EPS, a) - l_mt(f - EPS, a)) / (2 * EPS)
        else:
            num = (l_mt(f, a + EPS) - l_mt(f, a - EPS)) / (2 * EPS)
        c.setValue(state, v)
        (ca if name == "flexion" else cf).setValue(state, other)
        osa = mu.computeMomentArm(state, c)
        print(f"  [{label}] arm_{name} @ (flex {flex_deg}, add {add_deg}): "
              f"OpenSim {osa:+.6f} m vs dL/dq {num:+.6f} m "
              f"(|diff| {abs(osa-num):.2e})")


def profile(csv_name, flex_deg, dof_col, tau_col):
    rows = []
    with open(HERE / csv_name) as f:
        for r in csv.DictReader(f):
            if abs(float(r["hip_adduction_deg"]) - (0.0 if dof_col.endswith("flexion_m") else 0.0)) < 1e-6:
                pass
            rows.append(r)
    return rows


def axis_profile(csv_name, axis, fixed):
    """Rows where the conjugate DOF == fixed."""
    out = []
    with open(HERE / csv_name) as f:
        for r in csv.DictReader(f):
            conj = float(r["hip_adduction_deg"]) if axis == "flex" else float(r["hip_flexion_deg"])
            if abs(conj - fixed) < 1e-6:
                out.append(r)
    return out


print("=== moment-arm numerical cross-check (dL/dq, eps=1e-5 rad) ===")
check_arms(HERE / "variant_repo.osim", "repo", 45.0, 5.0)
check_arms(HERE / "variant_repo.osim", "repo", 0.0, 0.0)
check_arms(HERE / "variant_thumb.osim", "thumb", 45.0, 5.0)
check_arms(HERE.parent / "ConnorBipedal.osim", "human", 0.0, 0.0)

print("\n=== flexion-arm profile at adduction=0 deg (repo vs human) ===")
hu = axis_profile("arm_tau_human_ref_primary.csv", "flex", 0.0)
rp = axis_profile("arm_tau_variant_repo_primary.csv", "flex", 0.0)
th = axis_profile("arm_tau_variant_thumb_primary.csv", "flex", 0.0)
print(f"{'flex_deg':>8} | {'arm_repo':>9} {'arm_hum':>9} {'tauA_repo':>9} {'tauA_hum':>9} | {'arm_thumb':>9} {'tauA_thumb':>10}")
for r, h, t in zip(rp, hu, th):
    print(f"{float(r['hip_flexion_deg']):8.1f} | "
          f"{float(r['arm_hip_flexion_m'])*100:9.3f} {float(h['arm_hip_flexion_m'])*100:9.3f} "
          f"{float(r['tauA_hip_flexion_Nm']):9.2f} {float(h['tauA_hip_flexion_Nm']):9.2f} | "
          f"{float(t['arm_hip_flexion_m'])*100:9.3f} {float(t['tauA_hip_flexion_Nm']):10.2f}")

print("\n=== adduction-arm profile at flexion=0 deg (thumb vs human vs repo) ===")
hu = axis_profile("arm_tau_human_ref_primary.csv", "add", 0.0)
rp = axis_profile("arm_tau_variant_repo_primary.csv", "add", 0.0)
th = axis_profile("arm_tau_variant_thumb_primary.csv", "add", 0.0)
print(f"{'add_deg':>8} | {'arm_repo':>9} {'arm_thumb':>9} {'arm_hum':>9} | {'tauA_repo':>9} {'tauA_thumb':>10} {'tauA_hum':>9}")
for r, t, h in zip(rp, th, hu):
    print(f"{float(r['hip_adduction_deg']):8.1f} | "
          f"{float(r['arm_hip_adduction_m'])*100:9.3f} {float(t['arm_hip_adduction_m'])*100:9.3f} "
          f"{float(h['arm_hip_adduction_m'])*100:9.3f} | "
          f"{float(r['tauA_hip_adduction_Nm']):9.2f} {float(t['tauA_hip_adduction_Nm']):10.2f} "
          f"{float(h['tauA_hip_adduction_Nm']):9.2f}")
