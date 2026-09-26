# -*- coding: utf-8 -*-
import mujoco, os
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
schema = mujoco.mj_readSchema(None) if hasattr(mujoco, "mj_readSchema") else None
# fallback: parse the schema text from the docs embedded in _mujoco
import mujoco as mj
txt = None
try:
    from mujoco.schema import _SchemaDef
except Exception:
    pass
# Simplest: try a tiny model with candidate attributes and see which are accepted
cands = {
    "base": '<mujoco><worldbody/><actuator><muscle name="m" tendon="t" force="10" range="0 1" lengthrange="0 1" timeconst="0.01"/></actuator><tendon><spatial name="t"><site site="s"/></spatial></tendon><worldbody><site name="s" pos="0 0 0"/></worldbody></mujoco>',
}
attrs = ["damp", "scale", "shift", "delay", "gear", "lmin", "lmax", "forcerange", "ctrlrange", "actrange", "ctrllimited", "forcelimited", "gaintype", "gainprm", "biasprm", "dynprm"]
for a in attrs:
    xml = ('<mujoco><worldbody><site name="s" pos="0 0 0"/></worldbody>'
           '<tendon><spatial name="t"><site site="s"/></spatial></tendon>'
           '<actuator><muscle name="m" tendon="t" %s="0.05"/></actuator></mujoco>') % a
    try:
        m = mj.MjModel.from_xml_string(xml)
        print("%-12s ACCEPTED (value read: %s)" % (a, m.actuator_biasprm[0] if a == "damp" else ""))
    except Exception as e:
        print("%-12s rejected: %s" % (a, str(e)[:70]))
