"""Generate the vas_med_r SINGLE-MUSCLE knee torque surface for study 1
(minimum 20 mm BPAs vs human vastus medialis), 2026-10-27, easteregg2.

Mirrors gen_human_surfaces_2d_20261023.py conventions exactly (pinned
opensim 4.4.1, sha-pinned stock model, computeMomentArm on STOCK routing,
Fmax x arm, OpenSim coordinate direction, knee flexion NEGATIVE) but for
the single muscle vas_med_r on the 1-D graded knee grid:

  knee_angle_r nodes = 5 inclusive uniform nodes over rom_table_deg
  [-120, 10] -> [-120, -87.5, -55, -22.5, 10] deg
  (the same node convention the 5x5 survey grid uses on the knee axis,
   human_surfaces_5x5_grading_20261023.csv act 11 p-nodes verified).

tau_h(node) = arm(node) * Fmax(vas_med_r). All other DOFs at their model
defaults (asserted 0 for the campaign DOFs). A 1-deg dense sweep over
[-120, 10] is also exported and cross-checked against the machinery's own
OpenSim_Vasti_Results.txt column 5 (the extensor campaign's vas_med
target of record; interpolated comparison, reported).

Flag layer (same detectors as the 2-D generator, applied on the 1-D knee
axis): FLAG_VIA (conditional-path-point active-set change inside the
node's bin), FLAG_STEP (adjacent 0.5-deg |d arm| > 2% of surface-wide A*),
FLAG_SLACK / FLAG_SHORT on the 1-deg box. Flagged nodes still count in
the headline (conservative); nothing is corrected.

Outputs (this folder, all NEW, dated _20261027):
  vas_med_surface_5node_20261027.csv   - node rows + flags
  vas_med_surface_dense_20261027.csv   - 1-deg sweep (cross-check/viewing)
  vas_med_surface_manifest_20261027.json

Run:
  "D:/Anaconda/envs/myo/python.exe" tmp/paper_bridge/capability_table_20261022/gen_vas_med_surface_20261027.py
"""
import csv
import hashlib
import json
import math
import os
import sys
import xml.etree.ElementTree as ET

import opensim as osim

REPO = r'D:\GitHub\Bipedal_Robot'
OUTDIR = os.path.join(REPO, 'tmp', 'paper_bridge', 'capability_table_20261022')
MODEL = os.path.join(REPO, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody',
                     'gait2392_simbody.osim')
VASTI = os.path.join(REPO, 'Testing_Data', '2022_02_Festo',
                     'OpenSim_Vasti_Results.txt')
MUSCLE = 'vas_med_r'
KNEE = 'knee_angle_r'
KNEE_ROM = (-120.0, 10.0)     # rom_table_deg (actuator_map_27.json)
N_NODE = 5
BOX_DEG = 1.0
FINE_DEG = 0.5
STEP_PCT = 2.0
SHORT_LCEN = 0.5

PINNED_VERSION = '4.4.1'
PIN_MODEL_SHA = 'b89062c763e6d8fa3d7651c0a4fd2820055d807c'  # same pin as the 2-D generator

DEG2RAD = math.pi / 180.0


def sha1(path):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(65536), b''):
            h.update(blk)
    return h.hexdigest()


def nodes_for(lo, hi, n):
    return [lo + (hi - lo) * k / (n - 1) for k in range(n)]


def angles_for(lo, hi, step):
    n = int(round((hi - lo) / step))
    return [lo + step * k for k in range(n + 1)]


def bin_edges(nodes, lo, hi):
    edges = [lo]
    for k in range(1, len(nodes)):
        edges.append(0.5 * (nodes[k - 1] + nodes[k]))
    edges.append(hi)
    return edges


def cond_moving_from_xml(path, muscle):
    """Conditional/Moving path points of `muscle` from the raw XML
    (same detector family as the 2-D generator)."""
    out = {'cond': [], 'moving': []}
    tree = ET.parse(path)
    for fm in tree.getroot().iter('ForceSet'):
        pass
    root = tree.getroot()
    for mu in root.iter('Millard2012EquilibriumMuscle'):
        if mu.get('name') != muscle:
            continue
        gp = mu.find('GeometryPath')
        if gp is None:
            continue
        for cpp in gp.iter('ConditionalPathPoint'):
            out['cond'].append(cpp.get('name'))
        for mpp in gp.iter('MovingPathPoint'):
            out['moving'].append(mpp.get('name'))
    return out


def main():
    if osim.__version__ != PINNED_VERSION:
        sys.exit(f'opensim {osim.__version__} != pinned {PINNED_VERSION}')
    got = sha1(MODEL)
    if got != PIN_MODEL_SHA:
        sys.exit(f'model sha {got} != pinned {PIN_MODEL_SHA}')
    print(f'opensim {osim.__version__} | model sha {got} OK')

    model = osim.Model(MODEL)
    state = model.initSystem()
    coords = {}
    model_cs = model.getCoordinateSet()
    for i in range(model_cs.getSize()):
        c = model_cs.get(i)
        coords[c.getName()] = c
    if KNEE not in coords:
        sys.exit(f'coordinate {KNEE} missing')

    # all campaign DOFs at default 0 (same guard as the 2-D generator)
    campaign_dofs = ['hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r',
                     'knee_angle_r', 'ankle_angle_r', 'subtalar_angle_r',
                     'mtp_angle_r', 'lumbar_extension', 'lumbar_bending',
                     'lumbar_rotation']
    for d in campaign_dofs:
        if d in coords:
            dv = coords[d].getValue(state) / DEG2RAD
            if abs(dv) > 1e-9:
                sys.exit(f'coordinate {d} default {dv} deg, expected 0')

    mu = model.getMuscles().get(MUSCLE)
    if mu is None:
        sys.exit(f'muscle {MUSCLE} missing')
    try:
        nwrap = mu.getGeometryPath().getWrapSet().getSize()
    except Exception as e:
        sys.exit(f'wrap-set check failed: {e}')
    if nwrap != 0:
        sys.exit(f'WRAP HALT GATE: {MUSCLE} has {nwrap} live wrap objects')
    print('wrap halt gate: empty wrap set OK')

    tsl = mu.get_tendon_slack_length()
    lo = mu.get_optimal_fiber_length()
    phi0 = mu.get_pennation_angle_at_optimal()
    fmax = mu.get_max_isometric_force()
    print(f'{MUSCLE}: Fmax = {fmax:.1f} N, tsl = {tsl:.6f} m, '
          f'lo = {lo:.6f} m, phi0 = {phi0:.4f} rad')

    kcoord = coords[KNEE]

    def arm_at(deg):
        kcoord.setValue(state, deg * DEG2RAD)
        model.realizePosition(state)
        got_deg = kcoord.getValue(state) / DEG2RAD
        assert abs(got_deg - deg) < 1e-6, f'clamped {got_deg} != {deg}'
        a = mu.computeMomentArm(state, kcoord)
        lmt = mu.getLength(state)
        h = lo * math.sin(phi0)
        lcen = math.sqrt(max((lmt - tsl) ** 2 + h * h, 0.0)) / lo
        return a, lmt, lcen

    # ---- 1-deg dense box --------------------------------------------------
    box = angles_for(KNEE_ROM[0], KNEE_ROM[1], BOX_DEG)
    arm_box, lmt_box, lcen_box = [], [], []
    for d in box:
        a, lmt, lc = arm_at(d)
        arm_box.append(a)
        lmt_box.append(lmt)
        lcen_box.append(lc)
    A_star = max(abs(a) for a in arm_box)
    print(f'A* (surface-wide max |arm|) = {A_star:.6f} m')

    # ---- 0.5-deg fine sweep for VIA / STEP --------------------------------
    cmv = cond_moving_from_xml(MODEL, MUSCLE)
    print(f'conditional points: {cmv["cond"]} | moving points: {cmv["moving"]}')
    fine = angles_for(KNEE_ROM[0], KNEE_ROM[1], FINE_DEG)
    arm_fine, active_fine = [], []
    for d in fine:
        kcoord.setValue(state, d * DEG2RAD)
        model.realizePosition(state)
        arm_fine.append(mu.computeMomentArm(state, kcoord))
        gp = mu.getGeometryPath()
        act = []
        for i in range(gp.getPathPointSet().getSize()):
            pp = gp.getPathPointSet().get(i)
            try:
                act.append(bool(pp.isActive(state)))
            except Exception:
                act.append(True)
        active_fine.append(tuple(act))

    nodes = nodes_for(KNEE_ROM[0], KNEE_ROM[1], N_NODE)
    edges = bin_edges(nodes, KNEE_ROM[0], KNEE_ROM[1])

    def bin_flags(i):
        e0, e1 = edges[i], edges[i + 1]
        via = False
        for k in range(1, len(fine)):
            if e0 - 1e-9 <= fine[k] <= e1 + 1e-9:
                if active_fine[k] != active_fine[k - 1]:
                    via = True
        step_flag = False
        for k in range(1, len(arm_fine)):
            if e0 - 1e-9 <= fine[k] <= e1 + 1e-9:
                if abs(arm_fine[k] - arm_fine[k - 1]) > STEP_PCT / 100.0 * A_star:
                    step_flag = True
        slack = False
        short = False
        for k, d in enumerate(box):
            if e0 - 1e-9 <= d <= e1 + 1e-9:
                if lmt_box[k] <= tsl:
                    slack = True
                if lcen_box[k] < SHORT_LCEN:
                    short = True
        return via, step_flag, slack, short

    # ---- node rows --------------------------------------------------------
    rows = []
    for i, nd in enumerate(nodes):
        a, lmt, lc = arm_at(nd)
        via, stp, slk, shrt = bin_flags(i)
        rows.append({
            'node_idx': i + 1, 'knee_angle_deg': nd,
            'arm_m': a, 'fmax_N': fmax, 'tau_h_Nm': a * fmax,
            'l_mt_m': lmt, 'lceN_est': lc,
            'flag_via': int(via), 'flag_step': int(stp),
            'flag_slack': int(slk), 'flag_short': int(shrt),
        })
        print(f'node {nd:8.2f} deg: arm {a:+.6f} m -> tau_h {a * fmax:+8.3f} N*m'
              f'  flags via={int(via)} step={int(stp)} slack={int(slk)} short={int(shrt)}')

    # ---- cross-check vs OpenSim_Vasti_Results.txt col 5 --------------------
    vasti = []
    with open(VASTI) as f:
        for ln in f:
            parts = ln.split()
            if len(parts) == 5:
                try:
                    vasti.append((float(parts[1]), float(parts[4])))
                except ValueError:
                    pass
    xs = [v[0] for v in vasti]
    ys = [v[1] for v in vasti]
    my_dense = [a * fmax for a in arm_box]
    dmax = 0.0
    for d, mine in zip(box, my_dense):
        if d < xs[0] or d > xs[-1]:
            continue
        # linear interpolation of the vasti table
        j = max(0, min(len(xs) - 2, int((d - xs[0]) / (xs[1] - xs[0]))))
        t = (d - xs[j]) / (xs[j + 1] - xs[j])
        ref = ys[j] + t * (ys[j + 1] - ys[j])
        dmax = max(dmax, abs(ref - mine))
    print(f'cross-check vs OpenSim_Vasti_Results.txt col5 (vas_med): '
          f'max |interp diff| = {dmax:.4f} N*m over {len(box)} samples')

    # ---- outputs ----------------------------------------------------------
    out_csv = os.path.join(OUTDIR, 'vas_med_surface_5node_20261027.csv')
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    dens = os.path.join(OUTDIR, 'vas_med_surface_dense_20261027.csv')
    with open(dens, 'w', newline='') as f:
        f.write('knee_angle_deg,arm_m,tau_h_Nm,l_mt_m,lceN_est\n')
        for d, a, lmt, lc in zip(box, arm_box, lmt_box, lcen_box):
            f.write(f'{d:.4f},{a:.8f},{a * fmax:.6f},{lmt:.8f},{lc:.6f}\n')
    man = {
        'created': '2026-10-27', 'muscle': MUSCLE, 'knee_dof': KNEE,
        'knee_rom_deg': list(KNEE_ROM), 'nodes_deg': nodes,
        'fmax_N': fmax, 'tsl_m': tsl, 'lo_m': lo, 'phi0_rad': phi0,
        'A_star_m': A_star, 'opensim': osim.__version__,
        'model_sha1': got, 'model': MODEL,
        'sign_convention': 'OpenSim coordinate direction; knee flexion '
                           'NEGATIVE; computeMomentArm x Fmax (stock routing)',
        'flags': {c: sum(r[c] for r in rows) for c in
                  ('flag_via', 'flag_step', 'flag_slack', 'flag_short')},
        'crosscheck_vasti_col5_max_abs_diff_Nm': dmax,
        'source_script': __file__,
    }
    with open(os.path.join(OUTDIR, 'vas_med_surface_manifest_20261027.json'),
              'w') as f:
        json.dump(man, f, indent=2)
    print(f'WROTE {out_csv}')
    print(f'WROTE {dens}')
    print('VAS_MED SURFACE DONE')


if __name__ == '__main__':
    main()
