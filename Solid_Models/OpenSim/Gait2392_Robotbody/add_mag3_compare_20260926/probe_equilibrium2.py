import math
import opensim as osim

mp = r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\add_mag3_compare_20260926\variant_repo.osim"
model = osim.Model(mp)
model.finalizeFromProperties()
s = model.initSystem()
mu = model.getMuscles().get("add_mag3_r")
t = osim.Thelen2003Muscle.safeDownCast(mu)
print("Thelen2003Muscle methods:",
      [m for m in dir(t) if "qui" in m.lower() or "Forc" in m and "comput" in m.lower()])

# try: all other muscles appliesForce=false, then equilibrateMuscles
n = 0
for m in model.getMuscles():
    if m.getName() != "add_mag3_r":
        m.set_appliesForce(False)
        n += 1
model.setAllPropertiesEnabled() if False else None
s = model.initSystem()
mu = model.getMuscles().get("add_mag3_r")
mu.setActivation(s, 1.0)
try:
    model.equilibrateMuscles(s)
    print(f"with {n} others disabled: F={mu.getTendonForce(s):.2f} N "
          f"lce={mu.getFiberLength(s):.5f} m")
except Exception as e:
    print("still fails:", str(e)[:200])
