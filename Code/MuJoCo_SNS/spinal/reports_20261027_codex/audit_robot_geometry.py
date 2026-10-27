"""Audit original import joint frames against the current CAD instance pose."""
import json
from pathlib import Path
import re
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]


def rotation(axis, angle):
    axis = np.array(axis, dtype=float)
    if abs(angle) < 1e-14:
        return np.eye(3)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    cross = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    return np.eye(3) + np.sin(angle)*cross + (1-np.cos(angle))*(cross@cross)


def component_name(text):
    return text.replace('09_BA_001:', '09_BA_001-1/').replace('09_BA_001-1:', '09_BA_001-1/')


def main(directory=None):
    global HERE
    if directory is not None:
        HERE=Path(directory).resolve()
    cad = json.loads((HERE/'robot_mass_geometry.json').read_text())
    components = {x['name']: x for x in cad['components']}
    source = (REPO/'Code/Matlab/SNS_Simscape/mdl_humanoid_lower_ah001_DataFile.m').read_text()
    transforms = {}
    for index, field, value in re.findall(r'smiData\.RigidTransform\((\d+)\)\.(translation|angle|axis|ID) = (.*?);', source):
        fields = transforms.setdefault(int(index), {})
        fields[field] = value.strip('"') if field == 'ID' else np.fromstring(value.strip('[]'), sep=' ')
    by_id = {x['ID']: x for x in transforms.values() if x.get('ID')}
    joints = []
    for kind, index, ident in re.findall(r'smiData\.(RevoluteJoint|CylindricalJoint|SphericalJoint)\((\d+)\)\.ID = "(.*?)";', source):
        if not ident:
            continue
        a, b = ident[1:-1].split(':-:')
        # Restrict moving joints to the actual knee linkages and hip spheres.
        if not (kind == 'SphericalJoint' or all('04_' in x for x in (a, b))):
            continue
        anchor, axes = [], []
        for prefix, name in zip('BF', (a, b)):
            key = prefix + '[' + a.replace('09_BA_001-1:', '09_BA_001:') + ':-:' + b.replace('09_BA_001-1:', '09_BA_001:') + ']'
            frame = by_id[key]
            comp = components[component_name(name)]
            transform = np.array(comp['transform'])
            rot = transform[:9].reshape(3, 3).T
            anchor.append(rot@(frame['translation']/1000) + transform[9:12])
            axes.append(rot@rotation(frame['axis'], float(frame['angle'][0]))[:, 2])
        joints.append(dict(kind=kind, a=component_name(a), b=component_name(b),
                           anchor_a=anchor[0].tolist(), anchor_b=anchor[1].tolist(),
                           axis_a=axes[0].tolist(), axis_b=axes[1].tolist(),
                           separation_m=float(np.linalg.norm(anchor[0]-anchor[1])),
                           absolute_axis_dot=float(abs(np.dot(*axes)))))
    masses = []
    for comp in cad['components']:
        name = comp['name']
        if re.search(r'/0[2-5]_', name):
            density, material = 1200., 'solid Onyx (user specified)'
        elif 'POLYPROPYLENE' in name:
            density, material = 900., 'polypropylene (provisional density)'
        elif 'BRASS' in name:
            density, material = 8500., 'brass (provisional density)'
        else:
            density, material = 7850., 'steel (provisional density; includes fastener simplification)'
        volume = sum(x[3] for x in comp['mass_properties_density_one'])
        masses.append(dict(name=name, volume_m3=volume, density_kg_m3=density,
                           material=material, mass_kg=volume*density))
    total = sum(x['mass_kg'] for x in masses)
    # Generic import constraints are not welds. Define the rigid mechanical
    # members explicitly; keep all four links in each knee mechanism separate.
    membership = {}
    for name in components:
        short = name.split('/')[-1]
        for prefix, group in [('02_01_PE','pelvis'), ('03_01_FB','thigh_r'),
                              ('03_02_FH','thigh_r'), ('04_01_KT','thigh_r'),
                              ('03_03_FB','thigh_l'), ('03_04_FH','thigh_l'),
                              ('04_03_KT','thigh_l'), ('04_02_KB','shank_r'),
                              ('05_01_TI','shank_r'), ('04_04_KB','shank_l'),
                              ('05_02_TI','shank_l')]:
            if short.startswith(prefix):
                membership[name] = group
        if short.startswith(('04_05_BL', '04_06_FL')):
            membership[name] = ('back_' if short.startswith('04_05') else 'front_') + ('r' if short.endswith('-1') else 'l')
        if '6630T310' in name or '92997A108' in name or '90695A038' in name:
            membership[name] = 'thigh_r' if name.endswith('-1') else 'thigh_l'
    adjacency = {name:set() for name in components}
    for ident in by_id:
        if not ident.startswith('B[') or ':-:' not in ident:
            continue
        names = tuple(component_name(x) for x in ident[2:-1].split(':-:'))
        if all(x in components for x in names):
            adjacency[names[0]].add(names[1])
            adjacency[names[1]].add(names[0])
    attachments = []
    for name in components:
        if name in membership:
            continue
        visited, frontier, options = {name}, {name}, set()
        while frontier and not options:
            frontier = set.union(*(adjacency[x] for x in frontier)) - visited
            visited |= frontier
            options = {membership[x] for x in frontier if x in membership}
        if not options:
            raise RuntimeError('Unattached component '+name)
        # Joint pins touch multiple links; choosing which rigid member owns their
        # spin inertia is explicitly recorded as an approximation.
        membership[name] = sorted(options)[0]
        attachments.append(dict(name=name, assigned=membership[name], candidates=sorted(options)))
    groups = {}
    for name in components:
        groups.setdefault(membership[name], []).append(name)
    output = dict(components=masses, joints=joints, lower_body_mass_kg=total,
                  rigid_groups=groups, hardware_attachments=attachments,
                  spine_mass_kg=.12, payload_mass_kg=75-total-.12,
                  warning='CAD mass excludes pneumatic actuators, tubing, and unmodeled hardware; hardware densities provisional')
    (HERE/'robot_geometry_audit.json').write_text(json.dumps(output, indent=2))
    (HERE/'export_groups.json').write_text(json.dumps([dict(name=name, components=members) for name,members in sorted(groups.items(), key=lambda x:x[0]!='pelvis')], indent=2))
    print(json.dumps({k:v for k,v in output.items() if k not in ('components','joints','rigid_groups','hardware_attachments')}, indent=2))
    for joint in joints:
        print(joint['a'].split('/')[-1], joint['b'].split('/')[-1],
              'gap mm', round(1000*joint['separation_m'], 5), 'axis dot', joint['absolute_axis_dot'])


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--directory',type=Path)
    main(parser.parse_args().directory)
