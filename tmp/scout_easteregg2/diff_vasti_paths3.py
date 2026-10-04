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
        for obj in pps.find('objects'):
            loc = obj.find('location')
            frame = obj.find('socket_parent_frame')
            pts.append((obj.tag, obj.get('name'),
                        frame.text if frame is not None else obj.get('body'),
                        loc.text.strip() if loc is not None else 'MOVING/spline'))
        pw = gp.find('PathWrapSet')
        wraps = []
        if pw is not None:
            for w in (pw.find('objects') if pw.find('objects') is not None else pw):
                wraps.append(w.get('name'))
        mif = mus.find('max_isometric_force')
        ofl = mus.find('optimal_fiber_length')
        tsl = mus.find('tendon_slack_length')
        out[nm] = (pts, wraps, mif.text.strip(), ofl.text.strip(), tsl.text.strip())
    return out

names = {'vas_med_r', 'vas_lat_r', 'vas_int_r', 'bifemsh_r', 'add_mag3_r'}
adj = muscle_paths('ConnorBipedal_Vastus_Adjusted.osim', names)
stk = muscle_paths('gait2392_simbody.osim', names)
for nm in sorted(names):
    print('===', nm)
    for label, d in (('stock', stk), ('adjusted', adj)):
        if nm not in d:
            print(f'  {label:9s}: ABSENT')
            continue
        pts, wraps, mif, ofl, tsl = d[nm]
        print(f'  {label:9s}: {len(pts)} pts, MIF {mif}, OFL {ofl}, TSL {tsl}, wraps {wraps}')
        for p in pts:
            print('      ', p)
