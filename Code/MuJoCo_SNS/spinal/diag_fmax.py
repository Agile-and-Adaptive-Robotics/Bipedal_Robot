"""Fmax audit: which actuators ship with near-zero force capacity, and is
that the converter or the prunes? Print gainprm[2] for every actuator,
flag < 50 N, and show the raw <general> tag for the trunk set."""
import re

import mujoco

import runner

m = mujoco.MjModel.from_xml_path(str(runner.MODEL))
print("actuators with Fmax (gainprm[2]) < 50 N:")
for i in range(m.nu):
    fmax = m.actuator_gainprm[i, 2]
    if fmax < 50.0:
        nm = mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        print(f"  {nm:16s} Fmax {fmax:8.2f} N")

xml = open(str(runner.MODEL), encoding="utf-8").read()
for nm in ("ercspn_r", "ext_hal_r", "soleus_r"):
    mm = re.search(rf'<general[^>]*name="{nm}"[^>]*/>', xml)
    if mm:
        tag = mm.group(0)
        keep = re.search(r'gainprm="[^"]*"', tag)
        print(f"\n{nm}: {keep.group(0) if keep else tag[:160]}")
