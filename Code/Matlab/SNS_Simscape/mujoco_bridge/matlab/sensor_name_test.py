"""Verify MJCF sensor element names against mujoco 2.3.7 (myo env).

Loads a minimal muscle model with jointpos/jointvel/actuatorpos/
actuatorvel/actuatorfrc sensors. If it compiles, those element names
are valid for 2.3.7 and we can use them in the sensor-patched bridge model.
"""
import mujoco

xml = """
<mujoco>
  <worldbody>
    <body name="b"><joint name="j" type="hinge" limited="true" range="-1 1"/>
      <geom type="sphere" size=".1"/></body>
  </worldbody>
  <actuator><muscle name="m1" joint="j"/></actuator>
  <sensor>
    <jointpos name="jp" joint="j"/>
    <jointvel name="jv" joint="j"/>
    <actuatorpos name="al" actuator="m1"/>
    <actuatorvel name="av" actuator="m1"/>
    <actuatorfrc name="af" actuator="m1"/>
  </sensor>
</mujoco>
"""
m = mujoco.MjModel.from_xml_string(xml)
names = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_SENSOR, i) for i in range(m.nsensor)]
print("mujoco", mujoco.__version__, "-> nsensor", m.nsensor, "names", names)
