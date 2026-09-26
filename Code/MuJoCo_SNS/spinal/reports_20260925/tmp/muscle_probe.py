# -*- coding: utf-8 -*-
"""Probe muscle actuator internals in mujoco 2.3.7."""
import mujoco, os, numpy as np
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")

W = ('<mujoco><worldbody><body name="b" pos="0 0 1"><geom type="sphere" size="0.01" mass="0.1"/>'
     '<site name="s1" pos="0 0 0"/><site name="s2" pos="0 0 -0.05"/></body>'
     '<site name="s0" pos="0 0 0.2"/></worldbody>'
     '<tendon><spatial name="t"><site site="s0"/><site site="s1"/><site site="s2"/></spatial></tendon>'
     '<actuator>%s</actuator></mujoco>')

mus = '<muscle name="m" tendon="t" force="1500" range="0.0811 0.2589" lengthrange="0.0811 0.2589" timeconst="0.01 0.01"/>'
m1 = mujoco.MjModel.from_xml_string(W % mus)
print("muscle shortcut: gainprm =", m1.actuator_gainprm[0][:3], " biasprm =", m1.actuator_biasprm[0][:3])
print("  dyntype:", m1.actuator_dyntype[0], "dyntprm:", m1.actuator_dynprm[0][:3],
      "biastype:", m1.actuator_biastype[0], "gaintype:", m1.actuator_gaintype[0])
print("  gear:", m1.actuator_gear[0][:3])

# general with muscle types + custom biasprm[2] = -damp
gen = ('<general name="m" tendon="t" dyntype="muscle" gaintype="muscle" biastype="muscle" '
       'dyntprm="0.01 0.01" gainprm="1500 0 0" biasprm="0 0 -142.24"/>')
try:
    m2 = mujoco.MjModel.from_xml_string(W % gen)
    print("general muscle: gainprm =", m2.actuator_gainprm[0][:3], " biasprm =", m2.actuator_biasprm[0][:3])
    print("  dyntype:", m2.actuator_dyntype[0], "dynprm:", m2.actuator_dynprm[0][:3],
          "biastype:", m2.actuator_biastype[0], "gaintype:", m2.actuator_gaintype[0])
except Exception as e:
    print("general FAILED:", str(e)[:200])
