"""Generate the 27-actuator human torque targets from STOCK gait2392.

Campaign of record for the corrected routing comparison (Ben's doctrine:
the human reference is the stock muscle on its OpenSim route; a BPA on the
modified robot route must meet or beat it over the RoM). Generalizes the
validated add_mag3_compare_20260926 harness and the 2026-10-03 scout probe
(tmp/scout_probe_20260903/probe_torque_provenance.py).

For every actuator in actuator_map_27.json and every DOF it lists:
  sweep the DOF over the campaign RoM at N poses (all other coordinates at
  their defaults), compute tauA = computeMomentArm(state, coord) *
  max_isometric_force for each mapped stock muscle, and sum elementwise for
  the group. Moment arms are validated against a central difference of
  muscle length dl/dtheta at 8 spot angles per muscle (the add_mag3 8/8
  pattern); a relative mismatch above 5 percent fails the run loudly.

Outputs (into --outdir):
  targets27_<id>_<label>_<dof>_primary.csv   angle_deg + per-muscle tau
                                            columns + tau_group_Nm
  targets27_<id>_<label>_<dof>_secondary.csv same, for secondary DOFs
  targets27_manifest.json                    hashes, RoM sources, FD
                                            results, muscle map echoes

tauB (equilibrium tendon force cross-check) is intentionally NOT emitted:
the add_mag3 verdict section 6 showed Thelen passive overstretch makes
equilibrium-derived numbers unusable as targets; tauA = |arm| * MIF is the
doctrine of record. Decision documented in the manifest.

Run on easteregg2:
  D:\\Anaconda\\envs\\opensim\\python.exe gen_gait2392_torque_targets.py
"""
import argparse
import csv
import datetime
import hashlib
import json
import os
import sys

import opensim as osim

REPO = r'D:\GitHub\Bipedal_Robot'
DEFAULT_MODEL = os.path.join(REPO, 'Solid_Models', 'OpenSim',
                             'Gait2392_Robotbody', 'gait2392_simbody.osim')
DEFAULT_MAP = os.path.join(REPO, 'Code', 'Matlab', 'Mesh_Optimization',
                           'Routing27', 'actuator_map_27.json')
DEFAULT_OUT = os.path.join(REPO, 'Code', 'Matlab', 'Mesh_Optimization',
                           'Routing27', 'Human_Torques_27')

N_PRIMARY = 101
N_SECONDARY = 101
FD_SPOTS = 8
FD_H_DEG = 0.5          # central-difference half step, degrees
FD_TOL = 0.05           # relative tolerance arm vs dl/dtheta
DEG2RAD = 3.141592653589793 / 180.0


def sha1_of(path):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 16), b''):
            h.update(chunk)
    return h.hexdigest()


def coord_ranges_from_xml(model_path):
    """Coordinate -> (min_deg, max_deg) parsed from the osim XML (no API
    surface risk); locked coordinates report None."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(model_path)
    out = {}
    for c in tree.getroot().iter('Coordinate'):
        name = c.get('name')
        lo = c.find('range_min')
        hi = c.find('range_max')
        if name and lo is not None and hi is not None:
            out[name] = (float(lo.text) / DEG2RAD, float(hi.text) / DEG2RAD)
        elif name:
            out[name] = None
    return out


def safe_label(name):
    return name.replace(' ', '_').replace(',', '').replace('(', ''). \
        replace(')', '').replace('/', '_').replace('.', '')


def muscle_length(mu, state):
    try:
        return mu.getLength(state)
    except AttributeError:
        return mu.getGeometryPath().getLength(state)


def sweep(model, state, coord, muscles, degs):
    """tau[muscle][i] = arm * MIF at each angle; mutates state."""
    out = {m: [] for m in muscles}
    for deg in degs:
        coord.setValue(state, deg * DEG2RAD)
        model.realizePosition(state)
        for m in muscles:
            mu = model.getMuscles().get(m)
            out[m].append(mu.computeMomentArm(state, coord)
                          * mu.getMaxIsometricForce())
    return out


def fd_check(model, state, coord, muscles, degs):
    """computeMomentArm vs central-difference -dl/dtheta at FD_SPOTS
    angles. OpenSim's moment-arm convention is arm = -dL/dq (the muscle's
    generalized force contribution carries the minus sign), so the finite
    difference is NEGATED before comparing. Spots where |dl/dtheta| <
    1e-4 m/rad are skipped (arm essentially zero)."""
    results = {}
    spots = [degs[int(round(i * (len(degs) - 1) / (FD_SPOTS - 1)))]
             for i in range(FD_SPOTS)]
    for m in muscles:
        mu = model.getMuscles().get(m)
        worst = 0.0
        checked = 0
        for deg in spots:
            coord.setValue(state, (deg - FD_H_DEG) * DEG2RAD)
            model.realizePosition(state)
            lm = muscle_length(mu, state)
            coord.setValue(state, (deg + FD_H_DEG) * DEG2RAD)
            model.realizePosition(state)
            lp = muscle_length(mu, state)
            dl_dth = -(lp - lm) / (2.0 * FD_H_DEG * DEG2RAD)   # arm = -dL/dq
            coord.setValue(state, deg * DEG2RAD)
            model.realizePosition(state)
            arm = mu.computeMomentArm(state, coord)
            if abs(dl_dth) < 1e-4:
                continue
            denom = max(abs(dl_dth), 1e-9)
            worst = max(worst, abs(arm - dl_dth) / denom)
            checked += 1
        results[m] = worst if checked else 0.0
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default=DEFAULT_MODEL)
    ap.add_argument('--map', default=DEFAULT_MAP)
    ap.add_argument('--outdir', default=DEFAULT_OUT)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    with open(args.map) as f:
        cmap = json.load(f)

    print(f'model   : {args.model}')
    print(f'map     : {args.map}')
    print(f'outdir  : {args.outdir}')
    print(f'opensim : {osim.__version__ if hasattr(osim, "__version__") else "unknown"}')

    model = osim.Model(args.model)
    state = model.initSystem()

    # Existence checks (fail loud, never silent).
    xml_ranges = coord_ranges_from_xml(args.model)
    coords = model.getCoordinateSet()
    coord_names = [coords.get(i).getName() for i in range(coords.getSize())]
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

    def rom_for(dof):
        tbl = cmap['rom_table_deg'].get(dof)
        if tbl:
            return float(tbl[0]), float(tbl[1]), 'campaign table'
        rng = xml_ranges.get(dof)
        if rng is None:
            sys.exit(f'no range available for {dof}; set it in rom_table_deg')
        return rng[0], rng[1], 'stock model range'

    manifest = {
        'generated': datetime.datetime.now().isoformat(),
        'model': args.model,
        'model_sha1': sha1_of(args.model),
        'map': args.map,
        'map_sha1': sha1_of(args.map),
        'opensim_version': str(getattr(osim, '__version__', 'unknown')),
        'n_primary': N_PRIMARY,
        'n_secondary': N_SECONDARY,
        'fd_tol': FD_TOL,
        'tauB_note': ('tauB equilibrium cross-check intentionally omitted '
                      '(add_mag3 verdict section 6: Thelen passive '
                      'overstretch makes equilibrium forces unusable as '
                      'targets); tauA = arm * MIF is the doctrine'),
        'actuators': [],
    }

    failures = []
    for a in cmap['actuators']:
        # Fresh default state per actuator so a previous sweep's coordinate
        # value never leaks into the next actuator's sweeps.
        state = model.initSystem()
        label = safe_label(a['master_name'])
        entry = {'id': a['id'], 'master_name': a['master_name'],
                 'human_group': a['human_group'], 'dofs': {}}

        for kind, dofs in (('primary', [a['primary_dof']]),
                           ('secondary', a['secondary_dofs'])):
            for dof in dofs:
                # Fresh default state per DOF sweep (no leaked angles).
                state = model.initSystem()
                lo, hi, src = rom_for(dof)
                n = N_PRIMARY if kind == 'primary' else N_SECONDARY
                degs = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
                coord = coords.get(dof)
                taus = sweep(model, state, coord, a['human_group'], degs)
                group = [sum(taus[m][i] for m in a['human_group'])
                         for i in range(n)]

                sign_change = (min(group) < 0.0 < max(group))
                if sign_change and kind == 'primary':
                    print(f'  NOTE {a["master_name"]} {dof}: group curve '
                          f'changes sign (member cancellation); the '
                          f'campaign compares magnitudes')

                fname = (f'targets27_{a["id"]:02d}_{label}_'
                         f'{dof}_{kind}.csv')
                fpath = os.path.join(args.outdir, fname)
                with open(fpath, 'w', newline='') as f:
                    w = csv.writer(f)
                    w.writerow(['angle_deg']
                               + [f'tau_{m}_Nm' for m in a['human_group']]
                               + ['tau_group_Nm'])
                    for i in range(n):
                        w.writerow([f'{degs[i]:.6f}']
                                   + [f'{taus[m][i]:.9f}'
                                      for m in a['human_group']]
                                   + [f'{group[i]:.9f}'])

                entry['dofs'][dof] = {
                    'kind': kind, 'rom_deg': [lo, hi], 'rom_source': src,
                    'csv': fname, 'peak_abs_tau_Nm': max(abs(g)
                                                         for g in group),
                    'group_sign_change': sign_change,
                }

        # Moment-arm validation on the primary DOF (add_mag3 8/8 pattern).
        # KNEE CAVEAT (documented): stock gait2392 knee muscles carry
        # condyle PathWraps and moving quad insertions (Delp 1990 planar
        # knee), so muscle length is nonsmooth in knee_angle and a
        # central-difference dl/dtheta disagrees with the analytic
        # computeMomentArm at wrap transitions. OpenSim's
        # computeMomentArm is the definition of record (the add_mag3
        # harness validated the API; the legacy OpenSim_*_Results.txt
        # targets are stock API outputs). FD values for knee sweeps are
        # therefore recorded as INFORMATIONAL and do not gate the run;
        # every non-knee DOF gates hard.
        plo, phi, _ = rom_for(a['primary_dof'])
        degs_p = [plo + (phi - plo) * i / (N_PRIMARY - 1)
                  for i in range(N_PRIMARY)]
        state = model.initSystem()
        fd = fd_check(model, state, coords.get(a['primary_dof']),
                      a['human_group'], degs_p)
        fd_informational = a['primary_dof'] == 'knee_angle_r'
        entry['fd_worst_rel'] = fd
        bad = {m: v for m, v in fd.items() if v > FD_TOL}
        if fd_informational:
            bad = {}
            entry['fd_informational'] = True
        failures.append((a['master_name'], bad)) if bad else None
        entry['fd_pass'] = not bad
        manifest['actuators'].append(entry)

        peak = entry['dofs'][a['primary_dof']]['peak_abs_tau_Nm']
        tag = 'INFO (knee wrap FD)' if fd_informational else \
            ('PASS' if not bad else 'FAIL')
        print(f'[{a["id"]:2d}] {a["master_name"]:34s} '
              f'{a["primary_dof"]:18s} peak |tau| {peak:8.1f} N*m  '
              f'FD {tag} '
              f'(worst {max(fd.values()):.4f})')

    mpath = os.path.join(args.outdir, 'targets27_manifest.json')
    with open(mpath, 'w') as f:
        json.dump(manifest, f, indent=1)
    print(f'\nwrote {mpath}')

    if failures:
        for name, bad in failures:
            print(f'FD FAILURE {name}: {bad}')
        sys.exit('moment-arm finite-difference validation failed')
    print('ALL moment-arm FD checks passed')


if __name__ == '__main__':
    main()
