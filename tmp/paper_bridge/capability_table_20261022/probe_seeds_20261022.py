"""Quick probe (2026-10-22): for every robot_seed_preference muscle of
actuator_map_27.json, list the path-point body frames (in order) and the
PathWrap count in gait2392_robotbody.osim -- the two facts that decide
biPulleySpecsFromOpenSim buildability (wraps refused; single-frame refused;
unknown frames refused). Read-only.
"""
import json
import os
import xml.etree.ElementTree as ET

REPO = r'D:\GitHub\Bipedal_Robot'
OSIM = os.path.join(REPO, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody',
                    'gait2392_robotbody.osim')
MAP = os.path.join(REPO, 'Code', 'Matlab', 'Mesh_Optimization',
                   'Routing27_v2', 'actuator_map_27.json')

tree = ET.parse(OSIM)
root = tree.getroot()

muscles = {}
for m in root.iter('Thelen2003Muscle'):
    name = m.get('name')
    gp = m.find('GeometryPath')
    frames = []
    for p in gp.iter('PathPoint'):
        sock = p.find('socket_parent_frame')
        fr = sock.text.split('/')[-1] if sock is not None else '?'
        frames.append(fr)
    nwrap = len(list(gp.iter('PathWrap')))
    muscles[name] = (frames, nwrap)

with open(MAP) as f:
    cmap = json.load(f)

known = {'pelvis', 'femur_r', 'tibia_r', 'talus_r', 'calcn_r', 'toes_r',
         'torso'}
for a in cmap['actuators']:
    print(f"[{a['id']:2d}] {a['master_name']} ({a['primary_dof']})")
    for seed in a['robot_seed_preference']:
        if seed not in muscles:
            print(f"      {seed:12s} MISSING in robotbody")
            continue
        frames, nw = muscles[seed]
        # collapse consecutive duplicates for readability
        seq = []
        for fr in frames:
            if not seq or seq[-1] != fr:
                seq.append(fr)
        bad = [f for f in seq if f not in known]
        flag = ''
        if nw:
            flag += ' WRAPS=%d' % nw
        if bad:
            flag += ' UNKNOWN_FRAMES=%s' % bad
        if len(set(seq)) < 2:
            flag += ' SINGLE_FRAME'
        print(f"      {seed:12s} {seq}{flag}")
print('PROBE DONE')
