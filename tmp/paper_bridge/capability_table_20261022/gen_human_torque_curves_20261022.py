"""Generate the human torque curves for the 27-actuator capability table
(2026-10-22, easteregg2).

Single grading source for BOTH capability-table rows. Doctrine (ask of
2026-10-22 + actuator_map_27.json): the human target is the STOCK
gait2392_simbody muscle on its STOCK route; the group curve is the signed
elementwise sum of (moment arm x Fmax) over human_group, on a 25-point
inclusive grid across rom_table_deg[primary_dof]; secondary DOFs held at
the model default (verified 0 deg for all coordinates this run,
probe_model_20261022.py). Convention reused verbatim from the campaign
generator of record, Code/Matlab/Mesh_Optimization/Routing27_v2/
gen_gait2392_torque_targets.py:95-105 (tau = computeMomentArm *
max_isometric_force; fresh initSystem per actuator so no angle leaks).

Signs: OpenSim native. computeMomentArm returns arm = -dL/dq (gen_...py:108
-128; re-validated live below by central differences), so tau_h > 0 acts in
the coordinate's POSITIVE direction. Knee flexion is NEGATIVE in gait2392
(rom knee_angle_r [-120, 10]), so knee-extensor torques (vasti) are
positive and knee-flexor torques (hamstrings) negative in this CSV.

Artifact caveat: per-actuator string built from
path_artifact_audit_gait2392_simbody_20261022_034030.csv rows whose
(actuator_id, dof, kind) match this actuator's primary DOF and whose
artifact_flag=1. CAVEAT ONLY - curves are never corrected (directive).

Moment-arm validation: FD spot check (arm vs -dl/dtheta, 8 spots, 5% tol)
per the generator's add_mag3 8/8 pattern; knee_angle_r FD is informational
only (condyle wrap makes length nonsmooth; documented at
gen_gait2392_torque_targets.py:254-264); all other DOFs gate hard.

Outputs (this folder):
  human_torque_curves_20261022.csv   one row per actuator: id, name, dof,
                                     theta_deg_01..25, tau_h_Nm_01..25,
                                     artifact_flag_caveat; '# ' provenance
                                     header
  human_torque_curves_manifest_20261022.json   full conventions + hashes
                                     + per-muscle Fmax + FD + caveats

Run: "D:/Anaconda/envs/myo/python.exe" gen_human_torque_curves_20261022.py
"""
import csv
import datetime
import hashlib
import json
import os
import sys

import opensim as osim

REPO = r'D:\GitHub\Bipedal_Robot'
OUTDIR = os.path.join(REPO, 'tmp', 'paper_bridge',
                      'capability_table_20261022')
MODEL = os.path.join(REPO, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody',
                     'gait2392_simbody.osim')
MAP = os.path.join(REPO, 'Code', 'Matlab', 'Mesh_Optimization',
                   'Routing27_v2', 'actuator_map_27.json')
AUDIT = os.path.join(REPO, 'tmp', 'paper_bridge', 'opt_improve_20261022',
                     'path_artifact_audit_gait2392_simbody_'
                     '20261022_034030.csv')

PINNED_VERSION = '4.4.1'
N_GRID = 25                    # inclusive points
FD_SPOTS = 8
FD_H_DEG = 0.5
FD_TOL = 0.05
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


def muscle_length(mu, state):
    try:
        return mu.getLength(state)
    except AttributeError:
        return mu.getGeometryPath().getLength(state)


def load_audit(path):
    """(rows, header_comments) from the dated audit CSV. Keys used:
    actuator_id, dof, kind, muscle, artifact_flag, activation_angles_deg,
    step_angles_deg, max_step_pct, cond_in_range."""
    rows, comments = [], []
    with open(path, newline='') as f:
        for ln in f:
            if ln.startswith('#'):
                comments.append(ln.strip())
                continue
            break
    with open(path, newline='') as f:
        rdr = csv.DictReader(
            (ln for ln in f if not ln.startswith('#')))
        for r in rdr:
            rows.append(r)
    return rows, comments


def caveat_for(rows, act_id, dof):
    """Primary-DOF artifact caveat string for one actuator (never alters
    the curves)."""
    flagged = [r for r in rows
               if r['actuator_id'] == str(act_id) and r['dof'] == dof
               and r['kind'] == 'primary' and r['artifact_flag'] == '1']
    if not flagged:
        return 'none'
    parts = []
    for r in flagged:
        segs = [f"{r['muscle']}"]
        if r['cond_in_range'] == '1':
            segs.append('cond_in_range=1')
        if r['n_activation_angles'] not in ('', '0'):
            segs.append(f"act@{r['activation_angles_deg']}")
        if r['n_step_angles'] not in ('', '0'):
            segs.append(f"steps@{r['step_angles_deg']}")
        if r['max_step_pct']:
            segs.append(f"max_step={r['max_step_pct']}%")
        parts.append('[' + ' '.join(segs) + ']')
    return ' | '.join(parts) + ' (caveat only; curve NOT corrected)'


def fd_check(model, state, coord, muscles, degs):
    """arm vs central-difference -dl/dtheta at FD_SPOTS grid angles
    (gen_gait2392_torque_targets.py:108-138 pattern, re-run live)."""
    results = {}
    spots = [degs[int(round(i * (len(degs) - 1) / (FD_SPOTS - 1)))]
             for i in range(FD_SPOTS)]
    for m in muscles:
        mu = model.getMuscles().get(m)
        worst, checked = 0.0, 0
        for deg in spots:
            coord.setValue(state, (deg - FD_H_DEG) * DEG2RAD)
            model.realizePosition(state)
            lm = muscle_length(mu, state)
            coord.setValue(state, (deg + FD_H_DEG) * DEG2RAD)
            model.realizePosition(state)
            lp = muscle_length(mu, state)
            dl_dth = -(lp - lm) / (2.0 * FD_H_DEG * DEG2RAD)  # arm=-dL/dq
            coord.setValue(state, deg * DEG2RAD)
            model.realizePosition(state)
            arm = mu.computeMomentArm(state, coord)
            if abs(dl_dth) < 1e-4:
                continue
            worst = max(worst, abs(arm - dl_dth) / max(abs(dl_dth), 1e-9))
            checked += 1
        results[m] = worst if checked else 0.0
    return results


def main():
    ver = str(getattr(osim, '__version__', 'unknown'))
    if ver != PINNED_VERSION:
        sys.exit(f'OPENSIM VERSION PIN VIOLATION: {ver} != {PINNED_VERSION}')

    with open(MAP) as f:
        cmap = json.load(f)
    audit_rows, audit_comments = load_audit(AUDIT)

    # The audit must describe THIS model and THIS map (sha1s in its header).
    hdr = ' '.join(audit_comments)
    need = [sha1_of(MODEL), sha1_of(MAP)]
    for h in need:
        if h not in hdr:
            sys.exit(f'audit CSV does not pin current file sha1 {h}; '
                     'regenerate the audit before grading')
    print(f'audit pins model+map sha1: OK ({need[0][:8]}..., {need[1][:8]}...)')

    sig_model, sig_map, sig_audit = (stat_sig(MODEL), stat_sig(MAP),
                                     stat_sig(AUDIT))

    model = osim.Model(MODEL)
    state = model.initSystem()
    coords = model.getCoordinateSet()
    coord_names = set(coords.get(i).getName()
                      for i in range(coords.getSize()))
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
        sys.exit('map references missing muscles/DOFs')

    # Defaults of every DOF the map names (the 'held at default (0)' claim).
    defaults_deg = {}
    for d in sorted({a['primary_dof'] for a in cmap['actuators']}
                    | {d for a in cmap['actuators']
                       for d in a['secondary_dofs']}):
        defaults_deg[d] = coords.get(d).getValue(state) / DEG2RAD

    theta_cols = [f'theta_deg_{i + 1:02d}' for i in range(N_GRID)]
    tau_cols = [f'tau_h_Nm_{i + 1:02d}' for i in range(N_GRID)]
    out_csv = os.path.join(OUTDIR, 'human_torque_curves_20261022.csv')

    manifest = {
        'generated': datetime.datetime.now().isoformat(),
        'purpose': ('single grading source (both capability-table rows): '
                    'human group torque curves for the 27 actuators of '
                    'actuator_map_27.json (Routing27_v2)'),
        'opensim_version': ver,
        'python': sys.executable,
        'model': MODEL,
        'model_sha1': sha1_of(MODEL),
        'map': MAP,
        'map_sha1': sha1_of(MAP),
        'audit_csv': AUDIT,
        'audit_csv_sha1': sha1_of(AUDIT),
        'conventions': {
            'grid': (f'{N_GRID} points, INCLUSIVE endpoints, uniform: '
                     'theta_i = rom_lo + (rom_hi - rom_lo) * i / '
                     f'{N_GRID - 1}, i = 0..{N_GRID - 1}'),
            'rom_source': ('campaign rom_table_deg of actuator_map_27.json '
                           '(all 10 DOFs covered by the table)'),
            'torque_definition': ('tau_h(theta) = sum over human_group of '
                                  'computeMomentArm(theta) [m, OpenSim '
                                  'native] x max_isometric_force [N] '
                                  '(Fmax from the stock model); signed '
                                  'elementwise group sum, NEVER '
                                  'magnitude-only and never corrected'),
            'moment_arm_sign': ('computeMomentArm = -dL/dq (validated live '
                                'by central differences this run, see '
                                'fd_check per actuator), so tau_h > 0 acts '
                                'in the POSITIVE direction of the OpenSim '
                                'coordinate'),
            'coordinate_signs': ('OpenSim gait2392 native: knee flexion '
                                 'NEGATIVE (knee_angle_r rom [-120, 10], so '
                                 'extension is the positive direction: '
                                 'vasti group tau is positive, hamstrings '
                                 'negative); hip flexion POSITIVE '
                                 '(iliopsoas tau positive); ankle '
                                 'dorsiflexion POSITIVE'),
            'secondary_dofs': ('held at the model default pose; verified '
                              'this run: default 0 deg for every DOF in '
                              'the map (defaults_deg below)'),
            'artifact_caveat': ('column carries the primary-DOF artifact '
                                'flags of path_artifact_audit_...csv '
                                '(via-point activation jumps / moment-arm '
                                'steps); CAVEAT ONLY, curves never '
                                'corrected'),
        },
        'coordinate_defaults_deg': defaults_deg,
        'secondary_only_flagged_muscles_note': (
            'med_gas_r (actuator 13) and lat_gas_r (actuator 14) carry '
            'artifact flags on their SECONDARY knee_angle_r sweep in the '
            'audit; their shipped PRIMARY ankle_angle_r curves are '
            'unflagged, so their caveat column reads none'),
        'fd_settings': {'spots': FD_SPOTS, 'h_deg': FD_H_DEG, 'tol': FD_TOL,
                        'knee_informational': (
                            'knee_angle_r FD is informational only: stock '
                            'gait2392 knee paths carry condyle wraps/'
                            'moving points making length nonsmooth '
                            '(documented gen_gait2392_torque_targets.py:'
                            '254-264); every other DOF gates hard')},
        'actuators': [],
    }

    failures = []
    csv_rows = []
    for a in cmap['actuators']:
        state = model.initSystem()          # fresh: no angle leak
        dof = a['primary_dof']
        lo, hi = (float(x) for x in cmap['rom_table_deg'][dof])
        degs = [lo + (hi - lo) * i / (N_GRID - 1) for i in range(N_GRID)]
        coord = coords.get(dof)

        per_mu = {}
        group = [0.0] * N_GRID
        for i, deg in enumerate(degs):
            coord.setValue(state, deg * DEG2RAD)
            model.realizePosition(state)
            got = coord.getValue(state) / DEG2RAD
            assert abs(got - deg) < 1e-6, f'coordinate clamped: {got} != {deg}'
            for m in a['human_group']:
                mu = model.getMuscles().get(m)
                t = (mu.computeMomentArm(state, coord)
                     * mu.getMaxIsometricForce())
                per_mu.setdefault(m, {'fmax_n': mu.getMaxIsometricForce(),
                                      'peak_abs_tau_Nm': 0.0})
                per_mu[m]['peak_abs_tau_Nm'] = max(
                    per_mu[m]['peak_abs_tau_Nm'], abs(t))
                group[i] += t

        fd = fd_check(model, state, coord, a['human_group'], degs)
        knee = dof == 'knee_angle_r'
        bad = {} if knee else {m: v for m, v in fd.items() if v > FD_TOL}
        if bad:
            failures.append((a['master_name'], bad))

        caveat = caveat_for(audit_rows, a['id'], dof)
        fsum = sum(per_mu[m]['fmax_n'] for m in a['human_group'])
        entry = {
            'id': a['id'], 'master_name': a['master_name'],
            'primary_dof': dof, 'human_group': a['human_group'],
            'rom_deg': [lo, hi], 'n_grid': N_GRID,
            'per_muscle_fmax_n': {m: per_mu[m]['fmax_n']
                                  for m in a['human_group']},
            'group_fmax_sum_n': fsum,
            'master_mif_n_from_map': a['master_mif_n'],
            'peak_abs_tau_h_Nm': max(abs(g) for g in group),
            'tau_h_sign_change': min(group) < 0.0 < max(group),
            'artifact_flag_caveat': caveat,
            'fd_worst_rel': fd,
            'fd_status': 'informational (knee wrap)' if knee else
                         ('PASS' if not bad else 'FAIL'),
        }
        manifest['actuators'].append(entry)
        csv_rows.append(
            [a['id'], a['master_name'], dof]
            + [f'{g:.6f}' for g in degs]
            + [f'{g:.9f}' for g in group]
            + [caveat])
        print(f'[{a["id"]:2d}] {a["master_name"]:34s} {dof:18s} '
              f'peak|tau| {entry["peak_abs_tau_h_Nm"]:8.1f} N*m  '
              f'FD {entry["fd_status"]:24s} '
              f'worst {max(fd.values()):.4f}  '
              f'caveat {"FLAG" if caveat != "none" else "none"}')

    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([f'# human_torque_curves_20261022 - grading source for '
                    f'the 27-actuator capability table (both rows)'])
        w.writerow([f'# opensim_version {ver} (pinned {PINNED_VERSION}); '
                    f'python {sys.executable}'])
        w.writerow([f'# model {MODEL} sha1 {sha1_of(MODEL)}'])
        w.writerow([f'# map {MAP} sha1 {sha1_of(MAP)}'])
        w.writerow([f'# audit {AUDIT} sha1 {sha1_of(AUDIT)} '
                    f'(caveat column source)'])
        w.writerow([f'# grid: {N_GRID} inclusive uniform points over '
                    f'rom_table_deg[primary_dof] (campaign table)'])
        # NOTE: comment lines are kept comma-free so csv.writer never
        # quotes them; every '#' line then literally starts with '#'
        # under any reader (naive line-skip or parsed-row filter).
        w.writerow(['# tau_h(theta) = sum over human_group of '
                    'computeMomentArm x Fmax (OpenSim native signs; '
                    'arm = -dL/dq so tau>0 = positive coordinate '
                    'direction; knee flexion NEGATIVE: vasti positive / '
                    'hamstrings negative)'])
        w.writerow(['# secondary DOFs held at model default 0 deg '
                    '(verified); artifact_flag_caveat = via-point flags; '
                    'CAVEAT ONLY; curves never corrected'])
        w.writerow([f'# generated {datetime.datetime.now().isoformat()}'])
        w.writerow(['actuator_id', 'master_name', 'primary_dof']
                   + theta_cols + tau_cols + ['artifact_flag_caveat'])
        for r in csv_rows:
            w.writerow(r)
    print(f'\nwrote {out_csv} ({len(csv_rows)} actuator rows)')

    mpath = os.path.join(OUTDIR, 'human_torque_curves_manifest_20261022.json')
    with open(mpath, 'w') as f:
        json.dump(manifest, f, indent=1)
    print(f'wrote {mpath}')

    if (stat_sig(MODEL), stat_sig(MAP), stat_sig(AUDIT)) != \
            (sig_model, sig_map, sig_audit):
        sys.exit('READ-ONLY VIOLATION: an input changed during the run')

    if failures:
        for name, bad in failures:
            print(f'FD FAILURE {name}: {bad}')
        sys.exit('moment-arm finite-difference validation failed')
    print('ALL non-knee moment-arm FD checks passed '
          '(knee informational per documented wrap caveat)')
    print('DONE')


if __name__ == '__main__':
    main()
