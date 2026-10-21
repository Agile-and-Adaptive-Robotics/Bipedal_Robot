"""gen_required_torque_20261020.py -- Section-5 missing input (EDIT_SPEC
20261020 para 5): human required-torque curves tau_req(theta) = OpenSim
moment arm x max isometric force for the four batch sites that lack a
recorded human target (soleus_r, tib_ant_r, med_gas_r, rect_fem_r), plus
bifemsh_r as a cross-check against the recorded mono-mat ctx.humanTorque.

Model: stock human Gait2392 (gait2392_simbody.osim), the same file the
2026-10-03 provenance probes used as 'stock'. Style follows
Testing_Data/2022_02_Festo/OpenSim_*_Results.txt (moment arm x fmax, joint
coordinate swept, other coordinates at default 0). One .txt per
(muscle, joint) into tmp/paper_bridge/reqtorque_20261020/ (scratch) plus a
combined CSV. Writes NOTHING else; records untouched.

Run: D:/Anaconda/envs/myo/python.exe gen_required_torque_20261020.py
"""
import os
import opensim as osim

ROOT = r"D:\GitHub\Bipedal_Robot"
MODEL = os.path.join(ROOT, "Solid_Models", "OpenSim", "Gait2392_Robotbody",
                     "gait2392_simbody.osim")
OUTDIR = os.path.join(ROOT, "tmp", "paper_bridge", "reqtorque_20261020")

# (muscle, coordinate, deg_min, deg_max) -- ranges = the recorded batch
# thetaRangesDeg per site (BiPulley mats), so the requirement is evaluated
# on exactly the robot's design range.
SPECS = [
    ("soleus_r", "ankle_angle_r", -50.0, 20.0),
    ("tib_ant_r", "ankle_angle_r", -50.0, 20.0),
    ("med_gas_r", "knee_angle_r", -120.0, 10.0),
    ("med_gas_r", "ankle_angle_r", -50.0, 20.0),
    ("rect_fem_r", "hip_flexion_r", -15.0, 110.0),
    ("bifemsh_r", "knee_angle_r", -120.0, 10.0),
]
N = 100


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    model = osim.Model(MODEL)
    state = model.initSystem()
    ms = model.getMuscles()
    combined = []
    for muscle, qname, lo, hi in SPECS:
        coord = model.getCoordinateSet().get(qname)
        mu = ms.get(muscle)
        fmax = mu.getMaxIsometricForce()
        rows = []
        for i in range(N):
            deg = lo + (hi - lo) * i / (N - 1)
            coord.setValue(state, deg * 3.141592653589793 / 180.0)
            model.realizePosition(state)
            tau = mu.computeMomentArm(state, coord) * fmax
            rows.append((i, deg, tau))
        tag = "%s_%s" % (muscle, qname.replace("_angle_r", ""))
        txt = os.path.join(OUTDIR, "OpenSim_ReqTorque_%s_20261020.txt" % tag)
        with open(txt, "w") as f:
            f.write("UNKNOWN\nversion=1\nnRows=%d\nnColumns=3\n"
                    "inDegrees=no\nendheader\n" % N)
            f.write("time\t%s\ttau_req_Nm\n" % qname)
            for t, deg, tau in rows:
                f.write("%.8f\t%.8f\t%.8f\n" % (t, deg, tau))
        taus = [r[2] for r in rows]
        print("%-28s fmax=%7.1f N  tau_req [%.3f, %.3f] N*m  -> %s"
              % (tag, fmax, min(taus), max(taus), os.path.basename(txt)))
        combined.append((tag, fmax, rows))
    # combined CSV
    csvp = os.path.join(OUTDIR, "reqtorque_all_20261020.csv")
    with open(csvp, "w") as f:
        f.write("row,angle_deg," + ",".join(c[0] for c in combined) + "\n")
        for i in range(N):
            f.write("%d,%.6f," % (i, combined[0][2][i][1]))
            f.write(",".join("%.8f" % c[2][i][2] for c in combined) + "\n")
    print("combined ->", csvp)


if __name__ == "__main__":
    main()
