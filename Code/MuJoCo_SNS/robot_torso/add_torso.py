"""Add the approved torso to an already validated MJCF lower-body model.

The pelvis frame must have x anterior, y superior, and z right. The caller
supplies the measured lumbar hinge position in that frame. No muscles or
unmeasured attachment locations are silently synthesized.
"""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np


def values(sequence):
    return ' '.join(format(float(x), '.17g') for x in sequence)


def add_torso(source, destination, pelvis_name, hinge_position, spec=None):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if source == destination:
        raise ValueError('Write a separate model; preserve the source.')
    if source.parent != destination.parent:
        raise ValueError('Keep output beside source to preserve relative mesh/include paths.')
    if spec is None:
        spec = json.loads(Path(__file__).with_name('torso_spec.json').read_text())
    tree = ET.parse(source)
    root = tree.getroot()
    pelvis = root.find(".//body[@name='%s']" % pelvis_name)
    if pelvis is None or root.find(".//body[@name='lumbar_torso']") is not None:
        raise ValueError('Missing pelvis, or torso already present.')
    lower = mujoco.MjModel.from_xml_path(str(source))
    lower_mass = float(lower.body_mass.sum())
    length = float(spec['spine_length_m'])
    width = float(spec['spine_square_width_m'])
    bar_mass = length * width**2 * float(spec['spine_density_kg_m3'])
    payload_mass = float(spec['target_total_mass_kg']) - lower_mass - bar_mass
    if payload_mass <= 0 or min(length, width, spec['payload_sphere_radius_m']) <= 0:
        raise ValueError('Nonpositive payload mass or geometry size.')
    torso = ET.SubElement(pelvis, 'body', name='lumbar_torso', pos=values(hinge_position))
    # Explicit attributes prevent inherited damping/stiffness from masking control needs.
    joint_type=spec.get('lumbar_joint_type','ball')
    ET.SubElement(torso, 'joint', name='lumbar_ball', type=joint_type,
                  limited='false',
                  damping=str(spec['lumbar_damping_Nms_per_rad']), stiffness='0',
                  armature='0', frictionloss='0')
    ET.SubElement(torso, 'geom', name='spine_bar', type='box',
                  pos=values([0, length/2, 0]), size=values([width/2, length/2, width/2]),
                  mass=str(bar_mass), contype='0', conaffinity='0', rgba='0.25 0.25 0.25 1')
    fraction=float(spec.get('torso_payload_fraction',.5))
    if not 0 < fraction <= 1:
        raise ValueError('Torso payload fraction must lie in (0,1].')
    ET.SubElement(torso, 'geom', name='torso_payload', type='sphere',
                  pos=values(spec['payload_center_from_hinge_m']),
                  size=str(spec['payload_sphere_radius_m']), mass=str(payload_mass*fraction),
                  contype='0', conaffinity='0', rgba='0.65 0.65 0.7 0.65')
    if fraction < 1:
        # A rigid pelvis child supplies its own inertia without altering the
        # existing pelvis inertial element or adding a joint.
        ilium=ET.SubElement(pelvis,'body',name='ilium_payload',pos=values(hinge_position))
        ET.SubElement(ilium,'geom',name='ilium_payload_mass',type='sphere',
                      pos=values(spec['pelvis_payload_center_from_joint_m']),
                      size=str(spec['payload_sphere_radius_m']),mass=str(payload_mass*(1-fraction)),
                      contype='0',conaffinity='0',rgba='.5 .65 .7 .4')
    ET.SubElement(torso, 'site', name='spine_top', pos=values([0, length, 0]), size='0.005')
    ET.indent(tree, space='  ')
    tree.write(destination, encoding='utf-8', xml_declaration=True)
    model = mujoco.MjModel.from_xml_path(str(destination))
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    joint = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, 'lumbar_ball')
    dof = model.jnt_dofadr[joint]
    total = float(model.body_mass.sum())
    if not np.isclose(total, spec['target_total_mass_kg'], rtol=0, atol=1e-9):
        raise AssertionError('Compiled total mass does not match target.')
    audit = dict(lower_body_mass_kg=lower_mass, spine_mass_kg=bar_mass,
                 payload_mass_kg=payload_mass, total_mass_kg=total,
                 torso_payload_mass_kg=payload_mass*fraction,ilium_payload_mass_kg=payload_mass*(1-fraction),
                 frame='x anterior, y superior, z right',lumbar_joint_type=joint_type,
                 lumbar_gravity_compensation_Nm=data.qfrc_bias[dof:dof+3].tolist(),
                 neutral_upright_payload_moment_magnitude_Nm=
                 payload_mass*fraction * np.linalg.norm(model.opt.gravity) * abs(spec['payload_center_from_hinge_m'][0]),
                 source=str(source), output=str(destination),
                 status='Compiled mass and static torque audit only; no locomotion validation')
    destination.with_suffix('.torso_audit.json').write_text(json.dumps(audit, indent=2))
    return audit


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--pelvis', required=True)
    parser.add_argument('--hinge', nargs=3, type=float, required=True)
    args = parser.parse_args()
    print(json.dumps(add_torso(args.source, args.destination, args.pelvis, args.hinge), indent=2))
