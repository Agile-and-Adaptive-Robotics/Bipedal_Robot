"""Independent verification of human_surfaces_5x5_grading_20261023.csv
(2026-10-23). Read-only. Three checks:

A. box stats vs probe_osim_limits_2d_20261023.csv: for every (muscle,
   primary, secondary) surface the probe sampled with the SAME axis pair,
   compare min_lcen_est / min_lmt / max_abs_arm and locations (probe
   rounded 4-5 dp). For the five uniarticular-class surfaces the probe
   pinned the secondary at 0 (tag '(none)'); my box contains that line, so
   check mine <= probe min (lengths) and mine >= probe max (arm).
B. re-evaluate 3 CSV cells from scratch (fresh process state, direct API,
   no generator code) and compare tau_h to the CSV values exactly printed.
C. structure: 675 rows, 25 per actuator, envelope exclusions only on
   24-27, v1 outputs untouched (mtimes predate this run).
"""
import csv
import json
import math
import os

import opensim as osim

REPO = r'D:\GitHub\Bipedal_Robot'
FOLDER = os.path.join(REPO, 'tmp', 'paper_bridge', 'capability_table_20261022')
CSV5 = os.path.join(FOLDER, 'human_surfaces_5x5_grading_20261023.csv')
MAN = os.path.join(FOLDER, 'human_surfaces_manifest_20261023.json')
PROBE = os.path.join(FOLDER, 'probe_osim_limits_2d_20261023.csv')
MODEL = os.path.join(REPO, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody',
                     'gait2392_simbody.osim')
DEG2RAD = math.pi / 180.0

fails = []


def check(cond, msg):
    print(('PASS ' if cond else 'FAIL ') + msg)
    if not cond:
        fails.append(msg)


# ---------- load my outputs ----------
man = json.load(open(MAN))
rows = list(csv.DictReader(
    (ln for ln in open(CSV5) if not ln.startswith('#'))))
print(f'CSV rows: {len(rows)}; manifest actuators: {len(man["actuators"])}')

# ---------- C. structure ----------
check(len(rows) == 675, f'row count 675 (got {len(rows)})')
per_act = {}
for r in rows:
    per_act.setdefault(int(r['actuator_id']), []).append(r)
check(sorted(per_act) == list(range(1, 28)), 'actuator ids 1..27 present')
check(all(len(v) == 25 for v in per_act.values()),
      '25 cells per actuator')
env_acts = sorted({int(r['actuator_id']) for r in rows
                   if r['envelope_limited'] == 'yes'})
check(env_acts == [24, 25, 26, 27], f'envelope actuators == 24-27 '
      f'(got {env_acts})')
for aid in env_acts:
    nex = sum(1 for r in per_act[aid] if r['in_graded_domain'] == '0')
    check(nex == 10, f'act {aid}: 10 envelope-excluded cells (got {nex})')
non_env_excl = [r for aid, v in per_act.items() if aid not in env_acts
                for r in v if r['in_graded_domain'] == '0']
check(not non_env_excl, 'no exclusions outside envelope actuators')
flagged_cells = [r for r in rows if r['in_graded_domain'] == '1' and (
    r['flag_via'] == '1' or r['flag_step'] == '1'
    or r['flag_slack'] == '1' or r['flag_short'] == '1')]
check(len(flagged_cells) == man['reporting']['totals'][
    'n_flagged_graded_cells'], 'flagged-cell total matches manifest '
    f'({len(flagged_cells)} vs {man["reporting"]["totals"]["n_flagged_graded_cells"]})')
check(man['opensim_version'] == '4.4.1', 'manifest pins opensim 4.4.1')

# ---------- A. probe cross-check ----------
probe = list(csv.DictReader(
    (ln for ln in open(PROBE) if not ln.startswith('#'))))
pidx = {}
for p in probe:
    pidx[(p['muscle'], p['primary_dof'], p['secondary_dof'])] = p

n_cmp = 0
for e in man['actuators']:
    p_dof = e['axes']['primary']
    s_dof = e['axes']['secondary']
    for m, st in e['box_stats_1deg'].items():
        key = (m, p_dof, s_dof)
        if key in pidx:
            p = pidx[key]
            d_c = abs(st['min_lcen_est'] - float(p['min_lcen_est']))
            loc_me = st['min_lcen_est_at_deg'].split(';')
            loc_p = p['min_lcen_est_at_deg'].split(';')
            d_loc = max(abs(float(loc_me[0]) - float(loc_p[0])),
                        abs(float(loc_me[1]) - float(loc_p[1])))
            d_l = abs(st['min_lmt_m'] - float(p['min_lmt_m']))
            d_a = abs(st['max_abs_arm_m'] - float(p['max_abs_arm_m']))
            ok = (d_c <= 1e-4 and d_l <= 1e-5 and d_a <= 1e-5
                  and d_loc <= 1.0)
            n_cmp += 1
            if not ok:
                check(False, f'act {e["id"]} {m}: lcen '
                      f'{st["min_lcen_est"]}@{st["min_lcen_est_at_deg"]} vs '
                      f'{p["min_lcen_est"]}@{p["min_lcen_est_at_deg"]}; '
                      f'lmt {st["min_lmt_m"]} vs {p["min_lmt_m"]}; '
                      f'arm {st["max_abs_arm_m"]} vs {p["max_abs_arm_m"]}')
        else:
            key0 = (m, p_dof, '(none)')
            if key0 in pidx:
                p = pidx[key0]
                # my 2-D box contains the probe's secondary=0 line
                ok = (st['min_lcen_est'] <= float(p['min_lcen_est']) + 1e-9
                      and st['min_lmt_m'] <= float(p['min_lmt_m']) + 1e-9
                      and st['max_abs_arm_m']
                      >= float(p['max_abs_arm_m']) - 1e-9)
                n_cmp += 1
                if not ok:
                    check(False, f'act {e["id"]} {m}: containment vs '
                          f'(none) probe row violated '
                          f'(lcen {st["min_lcen_est"]} vs '
                          f'{p["min_lcen_est"]}; arm '
                          f'{st["max_abs_arm_m"]} vs {p["max_abs_arm_m"]})')
print(f'A: compared {n_cmp} muscle-surfaces against the probe CSV '
      f'(31 = every muscle-surface of the manifest)')
check(n_cmp == 31, f'probe comparisons cover every muscle-surface '
      f'({n_cmp}==31)')

# ---------- B. independent cell re-evaluation ----------
model = osim.Model(MODEL)
state = model.initSystem()
coords = model.getCoordinateSet()


def tau_direct(muscles, p_dof, s_dof, pd, sd):
    state = model.initSystem()
    coords.get(p_dof).setValue(state, pd * DEG2RAD)
    coords.get(s_dof).setValue(state, sd * DEG2RAD)
    model.realizePosition(state)
    tot = 0.0
    for m in muscles:
        mu = model.getMuscles().get(m)
        tot += mu.computeMomentArm(
            state, coords.get(p_dof)) * mu.getMaxIsometricForce()
    return tot


spots = [(13, -50.0, -55.0), (24, -15.0, 3.75), (1, 6.25, 7.5)]
for aid, pd, sd in spots:
    r = next(r for r in per_act[aid]
             if abs(float(r['p_node_deg']) - pd) < 1e-9
             and abs(float(r['s_node_deg']) - sd) < 1e-9)
    muscles = r['chosen_muscles'].split('+')
    got = tau_direct(muscles, r['primary_dof'], r['secondary_dof'], pd, sd)
    want = float(r['tau_h_Nm'])
    check(abs(got - want) <= 5e-7,
          f'act {aid} tau at ({pd},{sd}): direct {got:.9f} vs CSV '
          f'{want:.9f}')

# ---------- v1 untouched ----------
import datetime
run_start = datetime.datetime(2026, 10, 23, 12, 0).timestamp()
for name in ['human_torque_curves_20261022.csv',
             'human_torque_curves_group_targets_20261022.csv',
             'human_torque_curves_manifest_20261022.json']:
    p = os.path.join(FOLDER, name)
    mt = os.path.getmtime(p)
    check(mt < run_start, f'v1 file untouched today: {name} '
          f'(mtime {datetime.datetime.fromtimestamp(mt)})')

print('\n' + ('ALL CHECKS PASSED' if not fails
              else f'{len(fails)} FAILURES'))
raise SystemExit(1 if fails else 0)
