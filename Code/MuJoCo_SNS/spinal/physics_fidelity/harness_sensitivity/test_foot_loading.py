"""Regression checks for calibrated foot load and the archived doubled sensor convention."""
from pathlib import Path
import sys
import tempfile
import unittest

import mujoco
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from foot_loading import FootGroundLoads
from reference_angle_errors import compare_angles


class ContactLoadTest(unittest.TestCase):
    def test_each_ground_contact_is_counted_once_on_the_correct_side(self):
        names = ('calcn_r', 'toes_r', 'calcn_l', 'toes_l')
        bodies = ''.join(
            f'<body name="{name}"><geom name="g_{name}" type="sphere" size=".03"/></body>' for name in names
        )
        model = mujoco.MjModel.from_xml_string(
            '<mujoco><worldbody><geom name="floor" type="plane" size="1 1 .1"/>'
            + bodies
            + '</worldbody></mujoco>'
        )
        recorder = FootGroundLoads(model, 1)
        ground = model.geom('floor').id
        right_heel, right_toe, left_heel, left_toe = [model.geom('g_' + name).id for name in names]
        upward_normal = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

        def force(value):
            return np.array([value, 0.0, 0.0, 0.0, 0.0, 0.0])

        self.assertEqual(recorder.add(0, ground, right_heel, upward_normal, force(32.0)), ('heel', 'r'))
        np.testing.assert_array_equal(recorder.normal[0], [32.0, 0.0])
        self.assertEqual(recorder.add(0, ground, right_toe, upward_normal, force(25.0)), ('toe', 'r'))
        self.assertEqual(recorder.add(0, left_heel, ground, -upward_normal, force(51.0)), ('heel', 'l'))
        self.assertIsNone(recorder.add(0, right_toe, left_toe, upward_normal, force(300.0)))
        np.testing.assert_array_equal(recorder.normal[0], [57.0, 51.0])
        np.testing.assert_array_equal(recorder.upward[0], [57.0, 51.0])
        np.testing.assert_array_equal(recorder.count[0], [2, 1])

    def test_40_newton_pulses_cross_only_the_historical_threshold(self):
        t = np.arange(0.0, 5.0, 0.01)
        pulse = (np.sin(2 * np.pi * t) > 0.3).astype(float)
        physical = np.column_stack([40.0 * pulse, 40.0 * pulse])
        with tempfile.TemporaryDirectory() as directory:
            trace = Path(directory) / 'trajectory.npz'
            loads = Path(directory) / 'foot_ground.npz'
            np.savez(trace, t=t, q=np.zeros((len(t), 13)), neuro=np.zeros((len(t), 14)), contact=2 * physical)
            np.savez(loads, t=t, upward_N=physical)
            historical = compare_angles(trace, start=0.0, recorded_force_scale=1.0)
            calibrated = compare_angles(trace, start=0.0, physical_loads=loads)
            corrected_archive = compare_angles(trace, start=0.0, recorded_force_scale=0.5)
        self.assertGreaterEqual(historical['sides']['r']['cycles'], 3)
        self.assertEqual(calibrated['sides']['r']['cycles'], 0)
        self.assertEqual(calibrated['sides'], corrected_archive['sides'])


if __name__ == '__main__':
    unittest.main()
