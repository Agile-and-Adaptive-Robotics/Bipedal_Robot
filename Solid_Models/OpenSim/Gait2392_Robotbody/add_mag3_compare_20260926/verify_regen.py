"""Post-regeneration verification: gait2392_robot.osim (winner applied) must
reproduce the thumb-variant moment arms exactly."""
import math
import opensim as osim

RB = r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_robot.osim"
TH = r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\add_mag3_compare_20260926\variant_thumb.osim"


def arms(path, flex_deg, add_deg):
    m = osim.Model(path)
    m.finalizeFromProperties()
    s = m.initSystem()
    mu = m.getMuscles().get("add_mag3_r")
    cf = m.getCoordinateSet().get("hip_flexion_r")
    ca = m.getCoordinateSet().get("hip_adduction_r")
    cf.setValue(s, math.radians(flex_deg))
    ca.setValue(s, math.radians(add_deg))
    return (mu.getLength(s),
            mu.computeMomentArm(s, cf),
            mu.computeMomentArm(s, ca))


worst = 0.0
for flex, add in ((0, 0), (45, 5), (-25, -45), (85, 20)):
    a = arms(RB, flex, add)
    b = arms(TH, flex, add)
    d = max(abs(x - y) for x, y in zip(a, b))
    worst = max(worst, d)
    print(f"robot.osim @ (flex {flex:+3d}, add {add:+3d}): l_mt {a[0]:.6f} m, "
          f"arm_flex {a[1]:+.6f} m, arm_add {a[2]:+.6f} m | max dev vs thumb variant {d:.2e}")
assert worst < 1e-9, f"regenerated model deviates from thumb variant by {worst}"
print("REGENERATION VERIFIED: robot.osim add_mag3_r == thumb variant on all probes")
