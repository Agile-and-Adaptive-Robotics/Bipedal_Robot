"""Probe robot-side routes in gait2392_robotbody.osim for the 27-actuator set.

Pure XML (no OpenSim needed): for each candidate robot-seed muscle, count
path-point types (PathPoint / MovingPathPoint / ConditionalPathPoint),
PathWrap objects, list the frames each point lives on, and the muscle's
max_isometric_force. This sets the campaign spec-builder policy (which
muscles parse clean as plain via-point routes, which need the documented
moving/conditional/wrap simplification).
"""
import xml.etree.ElementTree as ET
import json
import sys

MODEL = (r'D:\Github\Bipedal_Robot\Solid_Models\OpenSim'
         r'\Gait2392_Robotbody\gait2392_robotbody.osim')
OUT = (r'D:\Github\Bipedal_Robot\tmp\routing_opt'
       r'\robot_route_probe.json')

# Robot-side seed preference per campaign actuator (first that parses clean
# wins; the map JSON documents every entry).
CANDIDATES = [
    'ercspn_r', 'intobl_r', 'extobl_r',
    'bifemlh_r', 'semimem_r', 'sar_r', 'tfl_r', 'rect_fem_r', 'grac_r',
    'bifemsh_r', 'vas_med_r', 'vas_int_r', 'vas_lat_r',
    'soleus_r', 'med_gas_r', 'lat_gas_r',
    'tib_post_r', 'tib_ant_r', 'per_brev_r', 'per_long_r', 'per_tert_r',
    'flex_dig_r', 'flex_hal_r', 'ext_dig_r', 'ext_hal_r',
    'glut_max3_r', 'glut_max2_r', 'glut_max1_r',
    'glut_med1_r', 'glut_min1_r',
    'add_mag3_r', 'add_mag2_r', 'add_mag1_r',
    'iliacus_r', 'psoas_r',
]

tree = ET.parse(MODEL)
root = tree.getroot()

report = {}
for mu in root.iter('Thelen2003Muscle'):
    name = mu.get('name')
    if name not in CANDIDATES:
        continue
    gp = mu.find('GeometryPath')
    if gp is None:
        continue
    entry = {'mif': float(mu.find('max_isometric_force').text)}
    counts = {}
    frames = []
    details = []
    for tag in ('PathPoint', 'MovingPathPoint', 'ConditionalPathPoint'):
        pts = gp.findall('.//' + tag)
        counts[tag] = len(pts)
        for p in pts:
            sock = p.find('socket_parent_frame')
            fr = (sock.text if sock is not None else '?').split('/')[-1]
            frames.append(fr)
            if tag != 'PathPoint':
                details.append({'tag': tag, 'frame': fr})
    counts['PathWrap'] = len(gp.findall('.//PathWrap'))
    entry['counts'] = counts
    entry['frames'] = frames
    entry['dynamic_details'] = details
    entry['appliesForce'] = (mu.find('appliesForce').text
                             if mu.find('appliesForce') is not None else '?')
    report[name] = entry

for name in CANDIDATES:
    if name not in report:
        report[name] = {'missing': True}

with open(OUT, 'w') as f:
    json.dump(report, f, indent=1)

# Console summary: clean = plain PathPoints only, >= 2, no wrap/dynamic
print(f'{"muscle":14s} {"PP":>3s} {"Mv":>3s} {"Cd":>3s} {"Wr":>3s} '
      f'{"F":>1s}  frames')
for name in CANDIDATES:
    e = report[name]
    if 'missing' in e:
        print(f'{name:14s} MISSING from model')
        continue
    c = e['counts']
    flag = 'Y' if e['appliesForce'] in ('?', 'true') else 'N'
    clean = (c['MovingPathPoint'] == 0 and c['ConditionalPathPoint'] == 0
             and c['PathWrap'] == 0 and c['PathPoint'] >= 2)
    tag = '' if clean else '  <-- needs simplification'
    print(f'{name:14s} {c["PathPoint"]:3d} {c["MovingPathPoint"]:3d} '
          f'{c["ConditionalPathPoint"]:3d} {c["PathWrap"]:3d}  {flag}  '
          f'{",".join(e["frames"])}{tag}')
print(f'\nwrote {OUT}')
sys.exit(0)
