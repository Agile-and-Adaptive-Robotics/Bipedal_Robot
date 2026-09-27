import math
import opensim as osim

mp = r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\add_mag3_compare_20260926\variant_repo.osim"
model = osim.Model(mp)
model.finalizeFromProperties()
s = model.initSystem()
mu = model.getMuscles().get("add_mag3_r")
cf = model.getCoordinateSet().get("hip_flexion_r")

# default pose: full-model equilibrium vs per-muscle equilibrium
mu.setActivation(s, 1.0)
try:
    model.equilibrateMuscles(s)
    f_full = mu.getTendonForce(s)
    l_full = mu.getFiberLength(s)
except Exception as e:
    f_full = l_full = None
    print("full equilibrium failed at default pose:", e)

mu.setActivation(s, 1.0)
mu.computeInitialFiberEquilibrium(s)
f_per = mu.getTendonForce(s)
l_per = mu.getFiberLength(s)
print(f"full: F={f_full} lce={l_full}")
print(f"per-muscle: F={f_per} lce={l_per}")
assert f_full is None or abs(f_full - f_per) < 1e-9 * max(1.0, f_full)

# sweep hip flexion across a wide range; does per-muscle equilibrium always succeed?
for deg in range(-25, 86, 10):
    s2 = model.getWorkingState() if hasattr(model, "getWorkingState") else s
    cf.setValue(s2, math.radians(deg))
    mu.setActivation(s2, 1.0)
    try:
        mu.computeInitialFiberEquilibrium(s2)
        print(f"flex {deg:+4d} deg: lce {mu.getFiberLength(s2):.4f} m, "
              f"F {mu.getTendonForce(s2):8.1f} N, "
              f"arm_flex {mu.computeMomentArm(s2, cf):+.4f} m")
    except Exception as e:
        print(f"flex {deg:+4d} deg: PER-MUSCLE EQUILIBRIUM FAILED: {e}")
