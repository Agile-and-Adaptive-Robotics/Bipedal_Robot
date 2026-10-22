"""Per-site OpenSim muscle-path artifact audit (plan Item 7, 2026-10-22).

Directive basis (PULLEY_SWEEP_DIRECTIVE_20261021.md section 6): before any
requirement envelope is trusted, verify per site that the Gait2392 muscle
path shows no via-point jumps (ConditionalPathPoint activation steps) and
no wrap-side flips across the FULL graded range used for the requirement.
Gait2392 represents wrapping with conditional via points; GeometryPath
behavior changed after OpenSim 4.1, so this audit pins and records the
opensim version (4.4.1 asserted; the repo's myo env reports exactly that,
probe _probe_osim_20261022.py this session).

What it does (READ-ONLY; model and map file mtimes are asserted unchanged):
  1. XML inventory per muscle: static PathPoint / ConditionalPathPoint /
     MovingPathPoint / PathWrap counts; each conditional point's coordinate
     and activation range (rad -> deg).
  2. Live API sweep per site actuator and per DOF (primary + secondaries,
     same no-leak pattern as gen_gait2392_torque_targets.py:209-218) over
     the graded RoM (map rom_table_deg) at --step-deg (default 0.5):
     moment arm, muscle length, and the per-point isActive(state) bit
     pattern (API exposure verified by probe this session).
  3. Flags per (site, DOF, muscle):
       cond_in_range  - a conditional point on this muscle keys on THIS
                        DOF with an activation range intersecting the
                        graded RoM (structural via-point-jump risk);
       activation_angles - live adjacent-sample changes of the active
                        point set (the actual via-point jump);
       step_angles    - |d(arm)| between adjacent samples > 2% of the
                        sweep's peak |arm| (moment-arm step);
       wrap_flips     - PathWrap-based flips; 0 by construction on the
                        Gait2392 twins (all PathWrapSets empty, verified:
                        184 sets, 0 <PathWrap> in gait2392_simbody.osim);
                        the check stays for future wrapped models.
     artifact_flag = cond_in_range OR activation_angles OR step_angles.
  4. Bone-penetration proxy = the point-set change signals above; a true
     STL proximity check against the MyoConverter STLs is TODO and never
     silent (not implemented here; directive asks for visual verify too).
  5. --curves additionally scans the RECORDED requirement curves of
     record (OpenSim_Vasti_Results.txt, OpenSim_Bifem_Results.txt) with
     the same step detector.

Output: <out>/path_artifact_audit_<model>_<stamp>.csv with the pinned
opensim_version in '# ' header comments (mirrors the generator's manifest
record at gen_gait2392_torque_targets.py:194).

Gate expectations (plan Item 7 gate): MUST flag vas_med_r and vas_lat_r
(conditional + moving points) and semiten_r (conditional points); MUST
report zero flags for a clean two-point muscle (glut_med1_r).

Run (gate command of record, easteregg2):
  "D:/Anaconda/envs/myo/python.exe" tmp/paper_bridge/opt_improve_20261022/audit_osim_path_artifacts_20261022.py --model gait2392_simbody --out tmp/paper_bridge/opt_improve_20261022/ [--curves]
"""
import argparse
import csv
import datetime
import hashlib
import json
import os
import sys
import xml.etree.ElementTree as ET

import opensim as osim

REPO = r'D:\GitHub\Bipedal_Robot'
MODELS = {
    'gait2392_simbody': os.path.join(
        REPO, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody',
        'gait2392_simbody.osim'),
    'gait2392_robotbody': os.path.join(
        REPO, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody',
        'gait2392_robotbody.osim'),
}
DEFAULT_MAP = os.path.join(
    REPO, 'Code', 'Matlab', 'Mesh_Optimization', 'Routing27_v2',
    'actuator_map_27.json')
DEFAULT_OUT = os.path.join(REPO, 'tmp', 'paper_bridge',
                           'opt_improve_20261022')
CURVES = [
    ('knee extensor requirement of record',
     os.path.join(REPO, 'Testing_Data', '2022_02_Festo',
                 'OpenSim_Vasti_Results.txt')),
    ('knee flexor requirement of record',
     os.path.join(REPO, 'Testing_Data', '2022_02_Festo',
                 'OpenSim_Bifem_Results.txt')),
]
PINNED_VERSION = '4.4.1'
DEG2RAD = 3.141592653589793 / 180.0


def sha1_of(path):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 16), b''):
            h.update(chunk)
    return h.hexdigest()


def stat_sig(path):
    st = os.stat(path)
    return (st.st_mtime_ns, st.st_size)


def short_name(socket_text):
    """'/jointset/knee_r/knee_angle_r' -> 'knee_angle_r'."""
    return socket_text.strip().split('/')[-1]


def muscle_path_inventory(model_path):
    """muscle -> dict(nStatic, nCond, nMoving, nPathWrap, cond=[{coord,
    lo_deg, hi_deg}], movingCoords=[...]) parsed from the osim XML."""
    tree = ET.parse(model_path)
    root = tree.getroot()
    out = {}
    for el in root.iter():
        tag = el.tag.split('}')[-1]
        if not tag.endswith('Muscle'):
            continue
        name = el.get('name')
        if not name:
            continue
        inv = {'nStatic': 0, 'nCond': 0, 'nMoving': 0, 'nPathWrap': 0,
               'cond': [], 'movingCoords': []}
        for sub in el.iter():
            stag = sub.tag.split('}')[-1]
            if stag == 'PathPoint' and sub.get('name'):
                inv['nStatic'] += 1
            elif stag == 'ConditionalPathPoint':
                inv['nCond'] += 1
                coord, rng = None, (None, None)
                for ch in sub:
                    ctag = ch.tag.split('}')[-1]
                    if ctag == 'socket_coordinate':
                        coord = short_name(ch.text or '')
                    elif ctag == 'range':
                        parts = (ch.text or '').split()
                        if len(parts) == 2:
                            rng = (float(parts[0]) / DEG2RAD,
                                   float(parts[1]) / DEG2RAD)
                inv['cond'].append({'coord': coord, 'lo': rng[0],
                                    'hi': rng[1]})
            elif stag == 'MovingPathPoint':
                inv['nMoving'] += 1
                for ch in sub:
                    if ch.tag.split('}')[-1].startswith('socket_') and \
                            ch.text and 'coordinate' in ch.tag:
                        inv['movingCoords'].append(short_name(ch.text))
            elif stag == 'PathWrap' and sub.get('name'):
                inv['nPathWrap'] += 1
        out[name] = inv
    return out


def rom_for(dof, rom_table, xml_ranges):
    tbl = rom_table.get(dof)
    if tbl:
        return float(tbl[0]), float(tbl[1]), 'campaign table'
    rng = xml_ranges.get(dof)
    if rng is None:
        sys.exit(f'no range available for {dof}; set it in rom_table_deg')
    return rng[0], rng[1], 'stock model range'


def angles_for(lo, hi, step):
    n = int(round((hi - lo) / step)) + 1
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def detect_steps(vals, degs, step_pct):
    """Indices (in degs[:-1] terms) where |diff| > step_pct% of peak."""
    peak = max(abs(v) for v in vals) if vals else 0.0
    if peak <= 0:
        return [], 0.0
    flags, worst = [], 0.0
    for i in range(len(vals) - 1):
        d = abs(vals[i + 1] - vals[i])
        worst = max(worst, d)
        if d / peak * 100.0 > step_pct:
            flags.append(round((degs[i] + degs[i + 1]) / 2.0, 2))
    return flags, worst / peak * 100.0


def detect_activation_changes(patterns, degs):
    """Angles where the active path-point set changes between adjacent
    samples (the live via-point jump)."""
    out = []
    for i in range(len(patterns) - 1):
        if patterns[i] != patterns[i + 1]:
            out.append(round((degs[i] + degs[i + 1]) / 2.0, 2))
    return out


def scan_storage_curve(path, step_pct):
    """Read an OpenSim .sto requirement curve; return (angle header,
    {muscle: (n, peak, max_step_pct, step_angles)})."""
    with open(path, errors='replace') as f:
        lines = f.read().splitlines()
    hdr = next(i for i, ln in enumerate(lines) if ln.strip() == 'endheader')
    cols = [c for c in lines[hdr + 1].replace('\t', ' ').split(' ') if c]
    rows = []
    for ln in lines[hdr + 2:]:
        if not ln.strip():
            continue
        vals = [float(v) for v in ln.split()]
        rows.append(vals)
    out = {}
    angle = [r[1] for r in rows]           # col 0 = time idx, col 1 = angle
    for j in range(2, len(cols)):
        name = cols[j]
        series = [r[j] for r in rows]
        steps, worst = detect_steps(series, angle, step_pct)
        out[name] = (len(series), max(abs(v) for v in series), worst, steps)
    return cols[1], out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='gait2392_simbody',
                    choices=sorted(MODELS.keys()))
    ap.add_argument('--model-path', default=None,
                    help='explicit osim path (overrides --model name)')
    ap.add_argument('--map', default=DEFAULT_MAP)
    ap.add_argument('--out', default=DEFAULT_OUT)
    ap.add_argument('--step-deg', type=float, default=0.5)
    ap.add_argument('--step-pct', type=float, default=2.0)
    ap.add_argument('--xml-only', action='store_true',
                    help='inventory only, no API sweep (twin scan mode)')
    ap.add_argument('--curves', action='store_true',
                    help='also step-scan the recorded requirement txts')
    args = ap.parse_args()

    model_path = args.model_path or MODELS[args.model]
    model_tag = args.model if not args.model_path else \
        os.path.splitext(os.path.basename(args.model_path))[0]
    os.makedirs(args.out, exist_ok=True)

    ver = str(getattr(osim, '__version__', 'unknown'))
    print(f'model   : {model_path}')
    print(f'map     : {args.map}')
    print(f'opensim : {ver} (pinned {PINNED_VERSION})')
    if ver != PINNED_VERSION:
        sys.exit(f'OPENSIM VERSION PIN VIOLATION: {ver} != {PINNED_VERSION}; '
                 'requirements must be recomputed on the pinned build')

    sig_model, sig_map = stat_sig(model_path), stat_sig(args.map)
    inv = muscle_path_inventory(model_path)
    print(f'muscles inventoried: {len(inv)}  '
          f'(cond points: {sum(v["nCond"] for v in inv.values())}, '
          f'moving: {sum(v["nMoving"] for v in inv.values())}, '
          f'pathwraps: {sum(v["nPathWrap"] for v in inv.values())})')

    with open(args.map) as f:
        cmap = json.load(f)
    rom_table = cmap.get('rom_table_deg', {})

    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    csv_path = os.path.join(args.out,
                            f'path_artifact_audit_{model_tag}_{stamp}.csv')
    rows = []

    if args.xml_only:
        for m, v in sorted(inv.items()):
            rows.append(dict(actuator_id='', master_name='(inventory)',
                             dof='', kind='xml', muscle=m, rom_lo='', rom_hi='',
                             n_static_pp=v['nStatic'], n_cond_pp=v['nCond'],
                             n_moving_pp=v['nMoving'],
                             n_pathwrap=v['nPathWrap'], cond_in_range='',
                             n_activation_angles=0, activation_angles_deg='',
                             n_step_angles=0, step_angles_deg='',
                             max_step_pct='', wrap_flips=0,
                             artifact_flag=(v['nCond'] + v['nMoving']) > 0))
    else:
        model = osim.Model(model_path)
        state = model.initSystem()
        coords = model.getCoordinateSet()
        coord_names = [coords.get(i).getName()
                       for i in range(coords.getSize())]
        mus_names = set(model.getMuscles().get(i).getName()
                        for i in range(model.getMuscles().getSize()))
        missing = []
        for a in cmap['actuators']:
            for m in a['human_group']:
                if m not in mus_names:
                    missing.append(f'{a["master_name"]}: muscle {m}')
            for d in [a['primary_dof']] + a['secondary_dofs']:
                if d not in coord_names:
                    missing.append(f'{a["master_name"]}: dof {d}')
        if missing:
            for m in missing:
                print(f'MISSING: {m}')
            sys.exit('map references missing muscles/DOFs; fix the map first')

        for a in cmap['actuators']:
            for kind, dofs in (('primary', [a['primary_dof']]),
                               ('secondary', a['secondary_dofs'])):
                for dof in dofs:
                    lo, hi, src = rom_for(dof, rom_table, {})
                    degs = angles_for(lo, hi, args.step_deg)
                    coord = coords.get(dof)
                    state = model.initSystem()      # no-leak fresh state
                    arms, lens, pats = {m: [] for m in a['human_group']}, \
                        {m: [] for m in a['human_group']}, \
                        {m: [] for m in a['human_group']}
                    for deg in degs:
                        coord.setValue(state, deg * DEG2RAD)
                        model.realizePosition(state)
                        for m in a['human_group']:
                            mu = model.getMuscles().get(m)
                            arms[m].append(mu.computeMomentArm(state, coord))
                            lens[m].append(mu.getLength(state))
                            pts = mu.getGeometryPath().getPathPointSet()
                            pats[m].append(tuple(
                                pts.get(i).isActive(state)
                                for i in range(pts.getSize())))
                    site_flag = False
                    for m in a['human_group']:
                        v = inv[m]
                        cond_hit = any(
                            c['coord'] == dof and c['lo'] is not None
                            and not (c['hi'] < lo or c['lo'] > hi)
                            for c in v['cond'])
                        act = detect_activation_changes(pats[m], degs)
                        steps, worst = detect_steps(arms[m], degs,
                                                    args.step_pct)
                        wrap_flips = 0 if v['nPathWrap'] == 0 else len(steps)
                        flag = bool(cond_hit or act or steps)
                        site_flag = site_flag or flag
                        rows.append(dict(
                            actuator_id=a['id'], master_name=a['master_name'],
                            dof=dof, kind=kind, muscle=m,
                            rom_lo=lo, rom_hi=hi,
                            n_static_pp=v['nStatic'], n_cond_pp=v['nCond'],
                            n_moving_pp=v['nMoving'],
                            n_pathwrap=v['nPathWrap'],
                            cond_in_range=int(cond_hit),
                            n_activation_angles=len(act),
                            activation_angles_deg=';'.join(str(x) for x in act),
                            n_step_angles=len(steps),
                            step_angles_deg=';'.join(str(x) for x in steps[:8]),
                            max_step_pct=round(worst, 4),
                            wrap_flips=wrap_flips, artifact_flag=int(flag)))
                    tag = 'FLAG' if site_flag else 'ok'
                    print(f'[{a["id"]:2d}] {a["master_name"][:30]:30s} '
                          f'{dof:18s} {kind:9s} {tag}')

    # Write CSV (deterministic column order; '# ' provenance header).
    cols = ['actuator_id', 'master_name', 'dof', 'kind', 'muscle',
            'rom_lo', 'rom_hi', 'n_static_pp', 'n_cond_pp', 'n_moving_pp',
            'n_pathwrap', 'cond_in_range', 'n_activation_angles',
            'activation_angles_deg', 'n_step_angles', 'step_angles_deg',
            'max_step_pct', 'wrap_flips', 'artifact_flag']
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([f'# opensim_version {ver} (pinned {PINNED_VERSION})'])
        w.writerow([f'# model {model_path} sha1 {sha1_of(model_path)}'])
        w.writerow([f'# map {args.map} sha1 {sha1_of(args.map)}'])
        w.writerow([f'# generated {datetime.datetime.now().isoformat()} '
                    f'step_deg {args.step_deg} step_pct {args.step_pct}'])
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])
    print(f'wrote {csv_path} ({len(rows)} muscle-DOF rows)')

    flagged = [r for r in rows if r['artifact_flag']]
    print(f'\nsite-DOF muscle rows flagged: {len(flagged)} of {len(rows)}')
    for gate_mus in ('vas_med_r', 'vas_lat_r', 'semiten_r', 'glut_med1_r'):
        hits = [r for r in rows if r['muscle'] == gate_mus]
        if not hits:
            print(f'gate {gate_mus:12s}: NOT IN AUDIT SCOPE')
            continue
        fl = any(r['artifact_flag'] for r in hits)
        detail = '; '.join(f'{r["dof"]}:flag={r["artifact_flag"]}'
                           f'(cond={r["cond_in_range"]},act='
                           f'{r["n_activation_angles"]},step='
                           f'{r["n_step_angles"]})' for r in hits)
        print(f'gate {gate_mus:12s}: flagged={fl}  [{detail}]')

    if args.curves:
        print('\nrecorded requirement-curve step scan (same 2%-of-peak rule):')
        for label, path in CURVES:
            if not os.path.exists(path):
                print(f'  {label}: MISSING {path}')
                continue
            ang_col, res = scan_storage_curve(path, args.step_pct)
            print(f'  {label} ({os.path.basename(path)}, angle col '
                  f'{ang_col.split("/")[-2] if "/" in ang_col else ang_col}):')
            for m, (n, peak, worst, steps) in res.items():
                tag = 'FLAG' if steps else 'ok'
                print(f'    {m:12s} n={n} peak={peak:9.4f} '
                      f'max_step={worst:6.3f}% steps@{steps if steps else "-"}'
                      f' -> {tag}')

    if stat_sig(model_path) != sig_model or stat_sig(args.map) != sig_map:
        sys.exit('READ-ONLY VIOLATION: model or map file changed during audit')
    print('read-only check: model and map mtimes/sizes unchanged')
    print('AUDIT DONE')


if __name__ == '__main__':
    main()
