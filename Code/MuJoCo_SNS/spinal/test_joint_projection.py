"""Prescribed poses must satisfy MuJoCo's own residuals for scalar equalities."""
import os
from pathlib import Path
import sys

os.environ.setdefault("CONDA_PREFIX", str(Path(sys.executable).parent))
import mujoco
import numpy as np

from joint_projection import project_joint_followers


def test_quartic_chains_references_fixed_and_inactive():
    model = mujoco.MjModel.from_xml_string("""
    <mujoco>
      <compiler angle="radian"/>
      <worldbody>
        <body><joint name="driver" type="slide" ref="0.25"/><geom size=".1"/></body>
        <body><joint name="first" type="slide" ref=".2"/><geom size=".1"/></body>
        <body><joint name="child" type="slide" ref="-.1"/><geom size=".1"/></body>
        <body><joint name="fixed" type="slide" ref=".3"/><geom size=".1"/></body>
        <body><joint name="free" type="slide" ref="-.35"/><geom size=".1"/></body>
      </worldbody>
      <equality>
        <joint joint1="child" joint2="first" polycoef=".1 2 0 0 0"/>
        <joint joint1="first" joint2="driver" polycoef="0 .2 -.4 .3 1.25"/>
        <joint joint1="fixed" polycoef="-.125 9 8 7 6"/>
        <joint joint1="free" polycoef=".9 0 0 0 0" active="false"/>
      </equality>
    </mujoco>""")
    data = mujoco.MjData(model)
    for displacement in (0., -.7, .8):
        data.qpos[:] = model.qpos0
        data.qpos[0] += displacement
        data.qpos[-1] = -.47
        project_joint_followers(model, data)
        poly = np.polynomial.polynomial.polyval(displacement, [0, .2, -.4, .3, 1.25])
        np.testing.assert_allclose(data.qpos, [.25 + displacement, .2 + poly, 2*poly, .175, -.47],
                                   atol=1e-14, rtol=0)
        mujoco.mj_forward(model, data)
        rows = data.efc_type == mujoco.mjtConstraint.mjCNSTR_EQUALITY
        assert rows.sum() == 3
        np.testing.assert_allclose(data.efc_pos[rows], 0., atol=1e-14, rtol=0)
    # Cached topology must not cache coefficient values or ignore toggled equalities.
    model.eq_data[1, 4] = .5
    if hasattr(data, "eq_active"):
        data.eq_active[-1] = 1
    else:
        model.eq_active[-1] = 1
    project_joint_followers(model, data)
    np.testing.assert_allclose(data.qpos[-1], -.35 + .9, atol=1e-14, rtol=0)
    mujoco.mj_forward(model, data)
    rows = data.efc_type == mujoco.mjtConstraint.mjCNSTR_EQUALITY
    assert rows.sum() == 4
    np.testing.assert_allclose(data.efc_pos[rows], 0., atol=1e-14, rtol=0)


if __name__ == "__main__":
    test_quartic_chains_references_fixed_and_inactive()
