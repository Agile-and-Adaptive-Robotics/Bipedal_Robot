"""Probe (2026-10-22): ground the manifest claims before generating the
human torque curves.

Checks, against STOCK gait2392_simbody.osim with pinned opensim 4.4.1:
  1. opensim.__version__ == 4.4.1 (pin).
  2. Every human_group muscle of actuator_map_27.json (Routing27_v2) exists
     and its max_isometric_force (echoed for traceability).
  3. Every primary/secondary DOF exists; each coordinate's default value
     (rad->deg) and XML range (deg).
  4. Coordinate.setValue does NOT clamp: set a value outside the model range
     (knee_angle_r at +10 deg, above a 0-deg max if that is the model range)
     and read it back.
Read-only: no file written.
"""
import json
import os
import sys
import xml.etree.ElementTree as ET

import opensim as osim

REPO = r'D:\GitHub\Bipedal_Robot'
MODEL = os.path.join(REPO, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody',
                     'gait2392_simbody.osim')
MAP = os.path.join(REPO, 'Code', 'Matlab', 'Mesh_Optimization',
                   'Routing27_v2', 'actuator_map_27.json')
DEG2RAD = 3.141592653589793 / 180.0

ver = str(getattr(osim, '__version__', 'unknown'))
print(f'opensim.__version__ = {ver}')
if ver != '4.4.1':
    sys.exit('PIN VIOLATION')

with open(MAP) as f:
    cmap = json.load(f)

dofs, muscles = [], []
for a in cmap['actuators']:
    if a['primary_dof'] not in dofs:
        dofs.append(a['primary_dof'])
    for d in a['secondary_dofs']:
        if d not in dofs:
            dofs.append(d)
    for m in a['human_group']:
        if m not in muscles:
            muscles.append(m)
print(f'{len(muscles)} unique muscles, {len(dofs)} unique DOFs in map')

# XML coordinate ranges (generator pattern, gen_gait2392_torque_targets.py:66)
tree = ET.parse(MODEL)
xml_rng = {}
for c in tree.getroot().iter('Coordinate'):
    name = c.get('name')
    lo, hi = c.find('range_min'), c.find('range_max')
    dv = c.find('default_value')
    if name:
        rng = (float(lo.text) / DEG2RAD, float(hi.text) / DEG2RAD) \
            if lo is not None and hi is not None else None
        dflt = float(dv.text) / DEG2RAD if dv is not None else None
        xml_rng[name] = (rng, dflt)

model = osim.Model(MODEL)
state = model.initSystem()
coords = model.getCoordinateSet()
mus = model.getMuscles()

print('\n-- muscles (Fmax N) --')
missing = []
for m in muscles:
    try:
        mu = mus.get(m)
        print(f'{m:14s} Fmax={mu.getMaxIsometricForce():8.1f}')
    except Exception:
        missing.append(m)
        print(f'{m:14s} MISSING')

print('\n-- DOFs (model XML range deg, default deg; live default deg) --')
for d in dofs:
    c = coords.get(d)
    rng, dflt = xml_rng.get(d, (None, None))
    live = c.getValue(state) / DEG2RAD
    print(f'{d:20s} xml_range={rng} xml_default={dflt} live_default={live:.6f}')

print('\n-- clamp probe: knee_angle_r set to +10 deg, read back --')
kc = coords.get('knee_angle_r')
kc.setValue(state, 10.0 * DEG2RAD)
model.realizePosition(state)
back = kc.getValue(state) / DEG2RAD
print(f'set 10.0 deg -> read {back:.6f} deg (clamp would show model max)')

if missing:
    sys.exit(f'MISSING MUSCLES: {missing}')
print('\nPROBE DONE')
