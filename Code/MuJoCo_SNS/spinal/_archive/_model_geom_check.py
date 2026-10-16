import mujoco
import runner

m = mujoco.MjModel.from_xml_path(str(runner.MODEL))
print("MODEL:", runner.MODEL)
for nm in ("pelvis", "femur_r", "tibia_r", "talus_r", "toe_r", "calcn_r"):
    b = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_BODY, nm)
    if b >= 0:
        print(nm, "pos", [round(float(x), 4) for x in m.body_pos[b]])
    else:
        print(nm, "N/A")
print("total mass", round(float(m.body_mass.sum()), 2))
