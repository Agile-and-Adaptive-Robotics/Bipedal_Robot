"""Self-check for the SINGLE-MUSCLE human_torque_curves_20261022.csv
(fresh process; 2026-10-22 evening, easteregg2).

Checks, in this fresh process:
  1. CSV structure: '#' provenance headers, 54-column header row
     (id, name, dof + 25 theta + 25 tau + caveat), exactly 27 data rows,
     actuator ids 1..27.
  2. Grid: theta values equal the uniform inclusive 25-point grid over
     rom_table_deg[primary_dof] from actuator_map_27.json (endpoints exact).
  3. Archive divergence: human_torque_curves_group_targets_20261022.csv
     exists; its curves differ from the new single-muscle curves exactly
     on the actuators whose target set changed (5, 11, 24, 25) and are
     identical (to the printed precision) on unchanged single-muscle
     actuators (spot: 1, 4, 12, 16, 26).
  4. Independent recompute (fresh Model load): actuator 5 (semimem_r ALONE
     - ruling drop), 24 (glut_max1+2+3 compartment sum), 16 (tib_ant_r
     duplicate-name fix) recomputed from scratch and compared to the CSV
     values at every one of the 25 points (tol 1e-6 N*m, printed 9 dp).
  5. Caveat column consistency vs the audit CSV: 'none' exactly when no
     chosen muscle has artifact_flag=1 on (actuator, primary dof).

Run: "D:/Anaconda/envs/myo/python.exe" selfcheck_human_torque_curves_singlemuscle_20261022.py
"""
import csv
import json
import os
import sys

import opensim as osim

REPO = r'D:\GitHub\Bipedal_Robot'
OUTDIR = os.path.join(REPO, 'tmp', 'paper_bridge',
                      'capability_table_20261022')
NEW_CSV = os.path.join(OUTDIR, 'human_torque_curves_20261022.csv')
OLD_CSV = os.path.join(OUTDIR,
                       'human_torque_curves_group_targets_20261022.csv')
MANIFEST = os.path.join(OUTDIR,
                        'human_torque_curves_manifest_20261022.json')
MODEL = os.path.join(REPO, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody',
                     'gait2392_simbody.osim')
MAP = os.path.join(REPO, 'Code', 'Matlab', 'Mesh_Optimization',
                   'Routing27_v2', 'actuator_map_27.json')
AUDIT = os.path.join(REPO, 'tmp', 'paper_bridge', 'opt_improve_20261022',
                     'path_artifact_audit_gait2392_simbody_'
                     '20261022_034030.csv')
DEG2RAD = 3.141592653589793 / 180.0
N_GRID = 25
PINNED_VERSION = '4.4.1'

CHOSEN = {5: ['semimem_r'],
          16: ['tib_ant_r'],
          24: ['glut_max1_r', 'glut_max2_r', 'glut_max3_r']}


def read_csv_rows(path):
    header_comments, rows = [], []
    with open(path, newline='') as f:
        rdr = csv.reader((ln for ln in f if not ln.startswith('#')))
        hdr = next(rdr)
        for r in rdr:
            rows.append(r)
    return hdr, rows, header_comments


def main():
    ok = True

    ver = str(getattr(osim, '__version__', 'unknown'))
    print(f'opensim : {ver} (pinned {PINNED_VERSION})')
    if ver != PINNED_VERSION:
        sys.exit('pin violation')

    with open(MAP) as f:
        cmap = json.load(f)
    rom = cmap['rom_table_deg']

    comments = [ln.strip() for ln in path_comments(NEW_CSV)]
    hdr, rows, _ = read_csv_rows(NEW_CSV)

    # 1. structure
    ncol = len(hdr)
    if ncol != 3 + 2 * N_GRID + 1:
        print(f'FAIL header has {ncol} columns; expected {3 + 2 * N_GRID + 1}')
        ok = False
    if hdr[:3] != ['actuator_id', 'master_name', 'primary_dof'] or \
            hdr[3:3 + N_GRID] != [f'theta_deg_{i+1:02d}' for i in range(N_GRID)] \
            or hdr[3 + N_GRID:3 + 2 * N_GRID] != [f'tau_h_Nm_{i+1:02d}' for i in range(N_GRID)] \
            or hdr[-1] != 'artifact_flag_caveat':
        print('FAIL header column names/order wrong')
        ok = False
    if len(rows) != 27:
        print(f'FAIL {len(rows)} data rows; expected 27')
        ok = False
    if [int(r[0]) for r in rows] != list(range(1, 28)):
        print('FAIL actuator ids not 1..27 in order')
        ok = False
    need_hdr_tokens = ['opensim_version 4.4.1', 'TARGET RULE Ben 2026-10-22',
                       'knee flexion NEGATIVE', 'CAVEAT ONLY']
    joined = ' | '.join(comments)
    for tok in need_hdr_tokens:
        if tok not in joined:
            print(f'FAIL provenance header missing token: {tok}')
            ok = False
    if not comments:
        print('FAIL no # provenance header lines found')
        ok = False
    print(f'check 1 structure: {"PASS" if ok else "FAIL"} '
          f'({len(rows)} rows, {ncol} cols, {len(comments)} header lines)')

    # 2. grid vs map rom_table
    grid_ok = True
    for r in rows:
        dof = r[2]
        lo, hi = (float(x) for x in rom[dof])
        thetas = [float(x) for x in r[3:3 + N_GRID]]
        want = [lo + (hi - lo) * i / (N_GRID - 1) for i in range(N_GRID)]
        if any(abs(a - b) > 5e-5 for a, b in zip(thetas, want)):
            print(f'FAIL grid actuator {r[0]} {dof}')
            grid_ok = False
    print(f'check 2 grids vs rom_table_deg: {"PASS" if grid_ok else "FAIL"}')
    ok = ok and grid_ok

    # 3. archive divergence
    _, old_rows, _ = read_csv_rows(OLD_CSV)
    old_by_id = {int(r[0]): r for r in old_rows}
    changed_ids, same_ids = {5, 11, 24, 25}, {1, 4, 12, 16, 26}
    div_ok = True
    for r in rows:
        aid = int(r[0])
        old = old_by_id[aid]
        new_tau = r[3 + N_GRID:3 + 2 * N_GRID]
        old_tau = old[3 + N_GRID:3 + 2 * N_GRID]
        differs = new_tau != old_tau
        if aid in changed_ids and not differs:
            print(f'FAIL actuator {aid}: expected single-muscle curve to '
                  f'differ from the archived group curve')
            div_ok = False
        if aid in same_ids and differs:
            print(f'FAIL actuator {aid}: target set unchanged, curve must '
                  f'match the archived group version at printed precision')
            div_ok = False
    print(f'check 3 archive divergence (changed {sorted(changed_ids)} / '
          f'unchanged spot {sorted(same_ids)}): '
          f'{"PASS" if div_ok else "FAIL"}')
    ok = ok and div_ok

    # 4. independent recompute of spot actuators in this fresh process
    model = osim.Model(MODEL)
    state = model.initSystem()
    coords = model.getCoordinateSet()
    rec_ok = True
    by_id = {int(r[0]): r for r in rows}
    for aid, chosen in CHOSEN.items():
        r = by_id[aid]
        dof = r[2]
        lo, hi = (float(x) for x in rom[dof])
        coord = coords.get(dof)
        worst = 0.0
        for i in range(N_GRID):
            deg = lo + (hi - lo) * i / (N_GRID - 1)
            coord.setValue(state, deg * DEG2RAD)
            model.realizePosition(state)
            tau = 0.0
            for m in chosen:
                mu = model.getMuscles().get(m)
                tau += (mu.computeMomentArm(state, coord)
                        * mu.getMaxIsometricForce())
            csv_val = float(r[3 + N_GRID + i])
            worst = max(worst, abs(tau - csv_val))
        if worst > 1e-6:
            print(f'FAIL recompute actuator {aid} {"+".join(chosen)}: '
                  f'max |diff| {worst:.3e} N*m')
            rec_ok = False
        else:
            print(f'  recompute actuator {aid:2d} {"+".join(chosen):35s} '
                  f'max |diff| {worst:.2e} N*m over all 25 points')
    print(f'check 4 independent recompute: {"PASS" if rec_ok else "FAIL"}')
    ok = ok and rec_ok

    # 5. caveat column vs audit flags for CHOSEN muscles
    audit_rows = []
    with open(AUDIT, newline='') as f:
        rdr = csv.DictReader((ln for ln in f if not ln.startswith('#')))
        audit_rows = list(rdr)
    with open(MANIFEST) as f:
        manifest = json.load(f)
    cav_ok = True
    for entry in manifest['actuators']:
        aid, dof, chosen = entry['id'], entry['primary_dof'], \
            entry['chosen_muscles']
        flagged = any(r['actuator_id'] == str(aid) and r['dof'] == dof
                      and r['kind'] == 'primary' and r['muscle'] in chosen
                      and r['artifact_flag'] == '1' for r in audit_rows)
        col = by_id[aid][-1]
        if flagged != (col != 'none'):
            print(f'FAIL caveat consistency actuator {aid}: audit flagged='
                  f'{flagged} column={col[:40]}')
            cav_ok = False
    print(f'check 5 caveat/audit consistency: {"PASS" if cav_ok else "FAIL"}')
    ok = ok and cav_ok

    flagged_muscles = sorted({r['muscle'] for r in audit_rows
                              if r['kind'] == 'primary'
                              and r['artifact_flag'] == '1'
                              and any(r['muscle'] == m and int(r['actuator_id']) == e['id']
                                      for e in manifest['actuators']
                                      for m in e['chosen_muscles'])})
    print(f'chosen muscles with primary-DOF artifact flags: '
          f'{flagged_muscles}')

    print('SELFCHECK ' + ('PASS' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)


def path_comments(path):
    class _C:
        def __init__(self, p):
            self.p = p
        def __iter__(self):
            with open(self.p, errors='replace') as f:
                for ln in f:
                    if ln.startswith('#'):
                        yield ln.rstrip('\n')
                    elif ln.strip():
                        return
    return _C(path)


if __name__ == '__main__':
    main()
