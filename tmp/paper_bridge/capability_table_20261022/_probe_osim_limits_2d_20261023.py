"""Live 2-D limits probe for the Gait2392 human surfaces (2026-10-23).

Purpose: ground the opensim_limits_checklist_20261023.md flag rules in LIVE
measurements on the pinned stack (opensim 4.4.1, gait2392_simbody.osim sha1
b89062c7..., actuator_map_27.json sha1 a97afffe...), complementing the 1-D
path-artifact audit of record
(path_artifact_audit_gait2392_simbody_20261022_034030.csv).

READ-ONLY: model and map mtimes/sizes asserted unchanged at exit. Outputs
land only in this folder (probe_osim_limits_2d_20261023.csv).

Per actuator (TARGET_RULE muscle sets copied verbatim from
gen_human_torque_curves_singlemuscle_20261022.py, Ben's 2026-10-22 ruling)
and per (primary x secondary) surface, all other of the actuator's DOFs
pinned at 0 deg (the generator's convention):

  1. live wrap-set size per muscle (expect 0: stock model has no PathWrap);
  2. tendon_slack_length / optimal_fiber_length / pennation at optimum /
     max_isometric_force;
  3. full-box sweep at --step-deg (default 1.0): moment arm of the primary
     DOF, muscle length, per-point isActive pattern;
  4. per-primary-column detectors (same rule as the audit): activation-set
     changes (via-point jumps) and |d(arm)| > step-pct% of that column's
     peak |arm|; boundary-angle span across columns answers "does the
     activation boundary move with the OTHER angle?";
  5. min muscle length over the box, min l_mt/(tsl+lo), and min estimated
     normalized fiber length from the constant-width pennation model
     lce = sqrt((l_mt-tsl)^2 + (lo*sin(phi0))^2), lceN = lce/lo;
  6. fine 0.5 deg 1-D pass at secondary=0 for reconciliation with the
     audit's numbers.

Run (easteregg2):
  "D:/Anaconda/envs/myo/python.exe" tmp/paper_bridge/capability_table_20261022/_probe_osim_limits_2d_20261023.py
"""
import csv
import datetime
import hashlib
import json
import math
import os
import sys
import xml.etree.ElementTree as ET

import opensim as osim

REPO = r'D:\GitHub\Bipedal_Robot'
FOLDER = os.path.join(REPO, 'tmp', 'paper_bridge', 'capability_table_20261022')
MODEL = os.path.join(REPO, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody',
                     'gait2392_simbody.osim')
MAP = os.path.join(REPO, 'Code', 'Matlab', 'Mesh_Optimization', 'Routing27_v2',
                   'actuator_map_27.json')
OUT_CSV = os.path.join(FOLDER, 'probe_osim_limits_2d_20261023.csv')
PINNED_VERSION = '4.4.1'
PIN_MODEL_SHA = 'b89062c763e6d8fa3d7651c0a4fd2820055d807c'
PIN_MAP_SHA = 'a97afffeeae85687810986273db6ba0151eba1a2'
DEG2RAD = math.pi / 180.0
STEP_DEG = 1.0
STEP_PCT = 2.0
FINE_DEG = 0.5

# TARGET_RULE copied verbatim from gen_human_torque_curves_singlemuscle_
# 20261022.py (Ben 2026-10-22 ruling): id -> chosen muscles.
TARGET_RULE = {
    1:  ['ercspn_r'], 2:  ['intobl_r'], 3:  ['extobl_r'],
    4:  ['bifemlh_r'], 5:  ['semimem_r'], 6:  ['sar_r'],
    7:  ['tfl_r'], 8:  ['rect_fem_r'], 9:  ['grac_r'],
    10: ['bifemsh_r'], 11: ['vas_int_r'], 12: ['soleus_r'],
    13: ['med_gas_r'], 14: ['lat_gas_r'], 15: ['tib_post_r'],
    16: ['tib_ant_r'], 17: ['per_brev_r'], 18: ['per_long_r'],
    19: ['per_tert_r'], 20: ['flex_dig_r'], 21: ['flex_hal_r'],
    22: ['ext_dig_r'], 23: ['ext_hal_r'],
    24: ['glut_max1_r', 'glut_max2_r', 'glut_max3_r'],
    25: ['add_mag1_r', 'add_mag2_r', 'add_mag3_r'],
    26: ['iliacus_r'], 27: ['psoas_r'],
}


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
    return (socket_text or '').strip().split('/')[-1]


def cond_points_from_xml(model_path, muscle_name):
    """[(coord, lo_deg, hi_deg)] for one muscle, from the osim XML."""
    tree = ET.parse(model_path)
    root = tree.getroot()
    out = []
    for el in root.iter():
        tag = el.tag.split('}')[-1]
        if not tag.endswith('Muscle') or el.get('name') != muscle_name:
            continue
        for sub in el.iter():
            stag = sub.tag.split('}')[-1]
            if stag != 'ConditionalPathPoint':
                continue
            coord, rng = None, (None, None)
            for ch in sub:
                ctag = ch.tag.split('}')[-1]
                if ctag == 'socket_coordinate':
                    coord = short_name(ch.text)
                elif ctag == 'range':
                    parts = (ch.text or '').split()
                    if len(parts) == 2:
                        rng = (float(parts[0]) / DEG2RAD,
                               float(parts[1]) / DEG2RAD)
            out.append((coord, rng[0], rng[1]))
    return out


def moving_count_from_xml(model_path, muscle_name):
    tree = ET.parse(model_path)
    n = 0
    for el in tree.getroot().iter():
        tag = el.tag.split('}')[-1]
        if tag.endswith('Muscle') and el.get('name') == muscle_name:
            for sub in el.iter():
                if sub.tag.split('}')[-1] == 'MovingPathPoint':
                    n += 1
    return n


def angles_for(lo, hi, step):
    n = int(round((hi - lo) / step)) + 1
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def detect_steps(vals, degs, step_pct):
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
    out = []
    for i in range(len(patterns) - 1):
        if patterns[i] != patterns[i + 1]:
            out.append(round((degs[i] + degs[i + 1]) / 2.0, 2))
    return out


def main():
    ver = str(getattr(osim, '__version__', 'unknown'))
    print(f'opensim : {ver} (pinned {PINNED_VERSION})')
    if ver != PINNED_VERSION:
        sys.exit(f'OPENSIM VERSION PIN VIOLATION: {ver} != {PINNED_VERSION}')
    msha, psha = sha1_of(MODEL), sha1_of(MAP)
    print(f'model   : {MODEL}')
    print(f'  sha1  : {msha}')
    print(f'map     : {MAP}')
    print(f'  sha1  : {psha}')
    if msha != PIN_MODEL_SHA:
        sys.exit(f'MODEL SHA PIN VIOLATION: {msha}')
    if psha != PIN_MAP_SHA:
        sys.exit(f'MAP SHA PIN VIOLATION: {psha}')
    sig_model, sig_map = stat_sig(MODEL), stat_sig(MAP)

    with open(MAP) as f:
        cmap = json.load(f)
    rom = cmap['rom_table_deg']
    by_id = {a['id']: a for a in cmap['actuators']}
    if sorted(TARGET_RULE) != sorted(by_id):
        sys.exit('TARGET_RULE ids do not match the map actuators')

    model = osim.Model(MODEL)
    state = model.initSystem()
    coords = model.getCoordinateSet()
    coord_names = [coords.get(i).getName() for i in range(coords.getSize())]
    for a in cmap['actuators']:
        for d in [a['primary_dof']] + a['secondary_dofs']:
            if d not in coord_names:
                sys.exit(f'map dof {d} missing in model')
        for m in TARGET_RULE[a['id']]:
            try:
                model.getMuscles().get(m)
            except Exception:
                sys.exit(f'target muscle {m} missing in model')

    rows = []
    for aid in sorted(TARGET_RULE):
        a = by_id[aid]
        p_dof = a['primary_dof']
        s_dofs = list(a['secondary_dofs'])
        p_lo, p_hi = rom[p_dof]
        p_degs = angles_for(p_lo, p_hi, STEP_DEG)
        fine_degs = angles_for(p_lo, p_hi, FINE_DEG)
        surfaces = [(sd, rom[sd]) for sd in s_dofs] or [(None, None)]
        for muscle in TARGET_RULE[aid]:
            mu = model.getMuscles().get(muscle)
            # --- static properties + live wrap check ---
            try:
                live_wrap = mu.getGeometryPath().getWrapSet().getSize()
            except Exception as e:
                live_wrap = f'ERR:{e}'
            tsl = mu.get_tendon_slack_length()
            lo = mu.get_optimal_fiber_length()
            fmax = mu.get_max_isometric_force()
            try:
                phi0 = mu.get_pennation_angle_at_optimal()
            except Exception:
                phi0 = float('nan')
            cond_xml = cond_points_from_xml(MODEL, muscle)
            n_moving = moving_count_from_xml(MODEL, muscle)
            cond_str = ';'.join(
                f'{c}:{lo_:.2f}..{hi_:.2f}' for c, lo_, hi_ in cond_xml)

            for s_dof, s_rng in surfaces:
                if s_dof is None:
                    s_degs = [0.0]
                    fine_s = [0.0]
                    s_tag = '(none)'
                else:
                    s_degs = angles_for(s_rng[0], s_rng[1], STEP_DEG)
                    fine_s = [0.0]
                    s_tag = s_dof

                # ---------- full-box sweep (coarse 2-D) ----------
                state = model.initSystem()
                pcoord = coords.get(p_dof)
                act_dofs = [p_dof] + s_dofs
                for d in act_dofs:
                    coords.get(d).setValue(state, 0.0)
                n_steps_2d, max_step_2d = 0, 0.0
                act_boundaries = set()
                n_act_changes_2d = 0
                min_lmt, min_at = None, None
                min_lcen, min_lcen_at = None, None
                max_abs_arm = 0.0
                n_pts = 0
                for s_deg in s_degs:
                    if s_dof is not None:
                        coords.get(s_dof).setValue(state, s_deg * DEG2RAD)
                    arms, pats = [], []
                    for p_deg in p_degs:
                        pcoord.setValue(state, p_deg * DEG2RAD)
                        model.realizePosition(state)
                        arms.append(mu.computeMomentArm(state, pcoord))
                        lmt = mu.getLength(state)
                        n_pts += 1
                        if min_lmt is None or lmt < min_lmt:
                            min_lmt, min_at = lmt, (p_deg, s_deg)
                        lcen = (lmt - tsl) / lo if lo > 0 else float('nan')
                        if not math.isnan(phi0):
                            h = lo * math.sin(phi0)
                            lcen = math.sqrt(
                                max((lmt - tsl) ** 2 + h * h, 0.0)) / lo
                        if min_lcen is None or lcen < min_lcen:
                            min_lcen, min_lcen_at = lcen, (p_deg, s_deg)
                        pts = mu.getGeometryPath().getPathPointSet()
                        pats.append(tuple(
                            pts.get(i).isActive(state)
                            for i in range(pts.getSize())))
                    max_abs_arm = max(
                        [max_abs_arm] + [abs(v) for v in arms])
                    steps, worst = detect_steps(arms, p_degs, STEP_PCT)
                    n_steps_2d += len(steps)
                    max_step_2d = max(max_step_2d, worst)
                    acts = detect_activation_changes(pats, p_degs)
                    n_act_changes_2d += len(acts)
                    for x in acts:
                        act_boundaries.add(x)
                    # reset primary for next column
                    pcoord.setValue(state, p_degs[0] * DEG2RAD)
                ab = sorted(act_boundaries)
                ab_span = (ab[-1] - ab[0]) if len(ab) > 1 else 0.0

                # ---------- fine 1-D pass at secondary = 0 ----------
                state = model.initSystem()
                for d in act_dofs:
                    coords.get(d).setValue(state, 0.0)
                arms_f, pats_f, lmts_f = [], [], []
                for p_deg in fine_degs:
                    pcoord.setValue(state, p_deg * DEG2RAD)
                    model.realizePosition(state)
                    arms_f.append(mu.computeMomentArm(state, pcoord))
                    lmts_f.append(mu.getLength(state))
                    pts = mu.getGeometryPath().getPathPointSet()
                    pats_f.append(tuple(
                        pts.get(i).isActive(state)
                        for i in range(pts.getSize())))
                fine_steps, fine_worst = detect_steps(
                    arms_f, fine_degs, STEP_PCT)
                fine_acts = detect_activation_changes(pats_f, fine_degs)
                fine_min_norm = min(
                    l / (tsl + lo) for l in lmts_f) if (tsl + lo) > 0 else \
                    float('nan')

                rows.append(dict(
                    actuator_id=aid, master_name=a['master_name'],
                    muscle=muscle, primary_dof=p_dof, secondary_dof=s_tag,
                    n_grid_points=n_pts, live_wrap_count=live_wrap,
                    n_cond_pp=len(cond_xml), n_moving_pp=n_moving,
                    cond_points_xml=cond_str,
                    phi0_deg=round(math.degrees(phi0), 2)
                    if not math.isnan(phi0) else 'n/a',
                    tsl_m=round(tsl, 5), lo_m=round(lo, 5),
                    fmax_n=round(fmax, 1),
                    min_lmt_m=round(min_lmt, 5),
                    min_lmt_at_deg=f'{min_at[0]:g};{min_at[1]:g}',
                    min_lmt_norm=round(min_lmt / (tsl + lo), 4),
                    min_lcen_est=round(min_lcen, 4),
                    min_lcen_est_at_deg=(
                        f'{min_lcen_at[0]:g};{min_lcen_at[1]:g}'),
                    primary_only_min_lmt_norm=round(fine_min_norm, 4),
                    n_act_changes_2d=n_act_changes_2d,
                    act_boundaries_deg=';'.join(str(x) for x in ab[:12]),
                    act_boundary_span_deg=round(ab_span, 2),
                    n_steps_2d=n_steps_2d,
                    max_step_pct_2d=round(max_step_2d, 4),
                    fine_steps_sec0_deg=';'.join(
                        str(x) for x in fine_steps[:8]),
                    fine_max_step_pct_sec0=round(fine_worst, 4),
                    fine_act_changes_sec0_deg=';'.join(
                        str(x) for x in fine_acts[:8]),
                    max_abs_arm_m=round(max_abs_arm, 5)))
                print(f'[{aid:2d}] {a["master_name"][:26]:26s} '
                      f'{muscle:12s} {p_dof:16s} x {s_tag:16s} '
                      f'n={n_pts:6d} steps2d={n_steps_2d:3d} '
                      f'maxstep={max_step_2d:6.3f}% '
                      f'actchg={n_act_changes_2d:3d} '
                      f'abspan={ab_span:5.2f} '
                      f'minL={min_lmt:.4f} minN='
                      f'{min_lmt / (tsl + lo):.3f} minLceN={min_lcen:.3f}')

    cols = ['actuator_id', 'master_name', 'muscle', 'primary_dof',
            'secondary_dof', 'n_grid_points', 'live_wrap_count',
            'n_cond_pp', 'n_moving_pp', 'cond_points_xml', 'phi0_deg',
            'tsl_m', 'lo_m', 'fmax_n', 'min_lmt_m', 'min_lmt_at_deg',
            'min_lmt_norm', 'min_lcen_est', 'min_lcen_est_at_deg',
            'primary_only_min_lmt_norm', 'n_act_changes_2d',
            'act_boundaries_deg', 'act_boundary_span_deg', 'n_steps_2d',
            'max_step_pct_2d', 'fine_steps_sec0_deg',
            'fine_max_step_pct_sec0', 'fine_act_changes_sec0_deg',
            'max_abs_arm_m']
    with open(OUT_CSV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([f'# opensim_version {ver} (pinned {PINNED_VERSION})'])
        w.writerow([f'# model {MODEL} sha1 {msha}'])
        w.writerow([f'# map {MAP} sha1 {psha}'])
        w.writerow([f'# generated {datetime.datetime.now().isoformat()} '
                    f'step_deg {STEP_DEG} step_pct {STEP_PCT} '
                    f'fine_deg {FINE_DEG}'])
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])
    print(f'\nwrote {OUT_CSV} ({len(rows)} muscle-surface rows)')

    if stat_sig(MODEL) != sig_model or stat_sig(MAP) != sig_map:
        sys.exit('READ-ONLY VIOLATION: model or map changed during probe')
    print('read-only check: model and map mtimes/sizes unchanged')
    print('PROBE DONE')


if __name__ == '__main__':
    main()
