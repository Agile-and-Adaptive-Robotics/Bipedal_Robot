import xml.etree.ElementTree as ET

def muscle_paths(fn, names):
    root = ET.parse(fn).getroot()
    out = {}
    for mus in root.iter('Thelen2003Muscle'):
        nm = mus.get('name')
        if nm not in names:
            continue
        pts = []
        gp = mus.find('GeometryPath')
        pps = gp.find('PathPointSet')
        for pp in pps:
            if pp.tag != 'PathPoint':
                continue
            loc = pp.find('location')
            pts.append((pp.get('name'), pp.get('body'),
                        loc.text.strip() if loc is not None else None))
        pw = gp.find('PathWrapSet')
        wraps = [w.get('name') for w in pw] if pw is not None else []
        mif = mus.find('max_isometric_force')
        out[nm] = (pts, wraps, mif.text.strip())
    return out

names = {'vas_med_r', 'vas_lat_r', 'vas_int_r', 'bifemsh_r', 'add_mag3_r', 'soleus_r', 'tib_ant_r'}
adj = muscle_paths('ConnorBipedal_Vastus_Adjusted.osim', names)
stk = muscle_paths('gait2392_simbody.osim', names)
for nm in sorted(names):
    print('===', nm)
    for label, d in (('stock', stk), ('adjusted', adj)):
        if nm not in d:
            print(f'  {label:9s}: ABSENT')
            continue
        pts, wraps, mif = d[nm]
        print(f'  {label:9s}: {len(pts)} pts, MIF {mif}, wraps {wraps}')
        for p in pts:
            print('      ', p)
