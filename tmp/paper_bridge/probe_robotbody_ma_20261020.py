"""probe_robotbody_ma_20261020.py -- diagnostic: moment arms of the FIVE
batch muscles in the robotbody model (gait2392_robotbody.osim -- the file
the batch specs were built from) over each site's recorded joint grid.
Purpose: decide whether the reconstructed BiPam_pulley torques carry a sign
convention flip or whether the robot routes genuinely pull the same
rotational way. Diagnostic only; writes nothing.
"""
import opensim as osim

MODEL = r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_robotbody.osim"
SPECS = [
    ("soleus_r", "ankle_angle_r", -50.0, 20.0),
    ("tib_ant_r", "ankle_angle_r", -50.0, 20.0),
    ("med_gas_r", "knee_angle_r", -120.0, 10.0),
    ("med_gas_r", "ankle_angle_r", -50.0, 20.0),
    ("rect_fem_r", "hip_flexion_r", -15.0, 110.0),
    ("bifemsh_r", "knee_angle_r", -120.0, 10.0),
]
N = 21


def main():
    model = osim.Model(MODEL)
    state = model.initSystem()
    ms = model.getMuscles()
    for muscle, qname, lo, hi in SPECS:
        coord = model.getCoordinateSet().get(qname)
        mu = ms.get(muscle)
        fmax = mu.getMaxIsometricForce()
        vals = []
        for i in range(N):
            deg = lo + (hi - lo) * i / (N - 1)
            coord.setValue(state, deg * 3.141592653589793 / 180.0)
            model.realizePosition(state)
            vals.append(mu.computeMomentArm(state, coord) * fmax)
        print("%-14s %-14s fmax=%7.1f  robotbody tau_req [%.2f, %.2f] N*m  "
              "endpoints %.2f / %.2f"
              % (muscle, qname, fmax, min(vals), max(vals), vals[0], vals[-1]))


if __name__ == "__main__":
    main()
