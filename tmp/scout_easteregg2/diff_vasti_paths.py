import xml.etree.ElementTree as ET
import sys

def muscle_paths(fn, names):
    tree = ET.parse(fn)
    root = tree.getroot()
    out = {}
    for mus in root.iter('Thelen2003Muscle'):
        nm = mus.get('name')
        if nm in names:
            pts = []
            gp = mus.find('GeometryPath')
            if gp is None:
                continue
            pps = gp.find('PathPointSet')
            for pp in pps.findall('PathPoint'):
                loc = pp.find('location')
                bod = pp.get('body')
                pts.append((pp.get('name'), bod, loc.text.strip() if loc is not None else None))
            pw = gp.find('PathWrapSet')
            wraps = [w.get('name') for w in pw.findall('PathWrap')] if pw is not None else []
            mif = mus.find('max_isometric_force')
            out[nm] = {'points': pts, 'wraps': wraps,
                       'mif': mif.text.strip() if mif is not None else None}
    return out

names = {'vas_med_r', 'vas_lat_r', 'vas_int_r', 'bifemsh_r', 'add_mag3_r', 'soleus_r', 'tib_ant_r'}
adj = muscle_paths('ConnorBipedal_Vastus_Adjusted.osim', names)
stk = muscle_paths('gait2392_simbody.osim', names)
for nm in sorted(names):
    a = adj.get(nm)
    s = stk.get(nm)
    print('===', nm)
    print('  stock    :', 'None' if s is None else f"{len(s['points'])} pts, MIF {s['mif']}, wraps {s['wraps']}")
    if s:
        for p in s['points']:
            print('     ', p)
    print('  adjusted :', 'None' if a is None else f"{len(a['points'])} pts, MIF {a['mif']}, wraps {a['wraps']}")
    if a:
        for p in a['points']:
            print('     ', p)
