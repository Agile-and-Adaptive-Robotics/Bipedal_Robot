"""Check mass balance and gravity torque against an independent analytic case."""
import json
from pathlib import Path
from contextlib import nullcontext

import mujoco
import numpy as np

from add_torso import add_torso


def main():
    results = {}
    for angle in ('degree', 'radian'):
        fixture = Path(__file__).parent / 'verification_fixtures' / angle
        fixture.mkdir(parents=True, exist_ok=True)
        with nullcontext(fixture) as directory:
            source = Path(directory) / 'lower.xml'
            output = Path(directory) / 'with_torso.xml'
            source.write_text('<mujoco><compiler angle="%s"/>'
                              '<option gravity="0 0 -9.81"/>'
                              '<worldbody><body name="pelvis">'
                              '<geom type="box" size=".1 .1 .1" mass="5"/>'
                              '</body></worldbody></mujoco>' % angle)
            audit = add_torso(source, output, 'pelvis', [0, 0, .1])
            assert np.isclose(audit['spine_mass_kg'], .12)
            assert np.isclose(audit['payload_mass_kg'], 69.88)
            expected = 69.88 * 9.81 * .1
            assert np.isclose(abs(audit['lumbar_gravity_compensation_Nm']), expected)
            model = mujoco.MjModel.from_xml_path(str(output))
            assert np.allclose(model.jnt_range[0], np.deg2rad([-30, 30]))
            results[angle] = audit
    path = Path(__file__).with_name('verification.json')
    path.write_text(json.dumps(results, indent=2))
    print('PASS: degree/radian sources, 75 kg mass balance, and analytic lumbar torque %.6f Nm' % expected)


if __name__ == '__main__':
    main()
