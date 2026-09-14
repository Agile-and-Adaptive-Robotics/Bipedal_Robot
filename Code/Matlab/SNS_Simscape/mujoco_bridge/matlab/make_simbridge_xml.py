"""Create the sensor-patched bridge model copy (cvt3.xml is never modified).

Writes gait2392_simbody_cvt3_simbridge.xml next to cvt3.xml (so the relative
Geometry\ paths resolve), adding a <sensor> section before </mujoco>:
  - knee_angle_r jointpos / jointvel
  - vas_med_r actuatorpos / actuatorvel / actuatorfrc
Element names verified against mujoco 2.3.7 (see sensor_name_test.py).
The file is regenerated idempotently; a previous <sensor> block injected by
this script is stripped first.
"""
from pathlib import Path

import mujoco

HERE = Path(__file__).resolve().parents[5]          # .../Bipedal_Robot
MODEL_DIR = HERE / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody" / "mjc" / "gait2392_simbody"
SRC = MODEL_DIR / "gait2392_simbody_cvt3.xml"
DST = MODEL_DIR / "gait2392_simbody_cvt3_simbridge.xml"

SENSOR = """
  <sensor>
    <jointpos name="knee_r_pos" joint="knee_angle_r"/>
    <jointvel name="knee_r_vel" joint="knee_angle_r"/>
    <actuatorpos name="vas_med_r_len" actuator="vas_med_r"/>
    <actuatorvel name="vas_med_r_vel" actuator="vas_med_r"/>
    <actuatorfrc name="vas_med_r_frc" actuator="vas_med_r"/>
  </sensor>
"""

text = SRC.read_text(encoding="utf-8")
start = text.find("  <sensor>")
if start != -1 and text.find("simbridge") == -1:
    end = text.find("</sensor>", start)
    text = text[:start] + text[end + len("</sensor>") + 1:]
text = text.replace("</mujoco>", SENSOR + "</mujoco>")
DST.write_text(text, encoding="utf-8")

m = mujoco.MjModel.from_xml_path(str(DST))
names = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_SENSOR, i) for i in range(m.nsensor)]
print(f"wrote {DST.name}: nsensor={m.nsensor} names={names} "
      f"nsensordata={m.nsensordata} timestep={m.opt.timestep}")

# 2 ms production variant: the spinal runner rewrites the option line to
# timestep 0.002 + implicitfast before loading (runner.py line ~160); this
# copy matches those semantics so the bridge plant steps at the SAME 2 ms
# as the SNS network (single-rate co-sim, exact parity).
DST2 = MODEL_DIR / "gait2392_simbody_cvt3_simbridge2.xml"
text2 = text.replace(
    '<option timestep="0.005" collision="predefined"/>',
    '<option timestep="0.002" collision="predefined" integrator="implicitfast"/>')
DST2.write_text(text2, encoding="utf-8")
m2 = mujoco.MjModel.from_xml_path(str(DST2))
print(f"wrote {DST2.name}: nsensor={m2.nsensor} timestep={m2.opt.timestep} "
      f"integrator={m2.opt.integrator}")
