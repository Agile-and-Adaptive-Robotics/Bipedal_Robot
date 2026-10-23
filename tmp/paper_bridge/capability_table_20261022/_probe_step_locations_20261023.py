"""Supplementary location probe (2026-10-23): WHERE the 2-D step events of
probe_osim_limits_2d_20261023.csv live, distinguishing (a) activation-boundary
jumps, (b) smooth-but-steep gradients (event magnitude scales with sampling
step), and (c) per-column small-peak / sign-crossing normalization artifacts;
plus the tendon-slack contours (l_mt < tendon_slack_length) for the muscles
whose graded boxes dip below slack.

READ-ONLY; pinned 4.4.1; outputs printed to stdout only.
Run: "D:/Anaconda/envs/myo/python.exe" <this file>
"""
import json
import math
import os
import sys

import opensim as osim

REPO = r'D:\GitHub\Bipedal_Robot'
MODEL = os.path.join(REPO, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody',
                     'gait2392_simbody.osim')
MAP = os.path.join(REPO, 'Code', 'Matlab', 'Mesh_Optimization', 'Routing27_v2',
                   'actuator_map_27.json')
PINNED = '4.4.1'
DEG2RAD = math.pi / 180.0

CASES = [
    # (muscle, primary, secondary or None, label)
    ('sar_r', 'knee_angle_r', 'hip_flexion_r', 'sartorius 1764-step case'),
    ('grac_r', 'knee_angle_r', 'hip_flexion_r', 'gracilis 7.06% case'),
    ('per_long_r', 'subtalar_angle_r', 'ankle_angle_r', 'per_long 10.23%'),
    ('per_brev_r', 'subtalar_angle_r', 'ankle_angle_r', 'per_brev 6.72%'),
    ('tib_post_r', 'subtalar_angle_r', 'ankle_angle_r', 'tib_post 817'),
    ('iliacus_r', 'hip_flexion_r', None, 'iliacus jitter regime'),
    ('psoas_r', 'hip_flexion_r', None, 'psoas jitter regime'),
]
SLACK_CASES = [
    ('semimem_r', 'knee_angle_r', 'hip_flexion_r'),
    ('tfl_r', 'knee_angle_r', 'hip_flexion_r'),
    ('med_gas_r', 'ankle_angle_r', 'knee_angle_r'),
    ('lat_gas_r', 'ankle_angle_r', 'knee_angle_r'),
]


def angles_for(lo, hi, step):
    n = int(round((hi - lo) / step)) + 1
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def sweep(model, muscle, p_dof, s_dof, rom, step=1.0):
    state = model.initSystem()
    coords = model.getCoordinateSet()
    mu = model.getMuscles().get(muscle)
    pcoord = coords.get(p_dof)
    p_lo, p_hi = rom[p_dof]
    p_degs = angles_for(p_lo, p_hi, step)
    s_degs = [0.0] if s_dof is None else angles_for(*rom[s_dof], step=step)
    cols = {}
    lmt_min_overall = None
    for s_deg in s_degs:
        if s_dof is not None:
            coords.get(s_dof).setValue(state, s_deg * DEG2RAD)
        arms = []
        for p_deg in p_degs:
            pcoord.setValue(state, p_deg * DEG2RAD)
            model.realizePosition(state)
            arms.append(mu.computeMomentArm(state, pcoord))
        cols[s_deg] = arms
    return p_degs, s_degs, cols


def main():
    ver = str(getattr(osim, '__version__', 'unknown'))
    if ver != PINNED:
        sys.exit(f'PIN VIOLATION {ver}')
    with open(MAP) as f:
        rom = json.load(f)['rom_table_deg']
    model = osim.Model(MODEL)
    model.initSystem()

    for muscle, p_dof, s_dof, label in CASES:
        p_degs, s_degs, cols = sweep(model, muscle, p_dof, s_dof, rom)
        surf_peak = max(max(abs(v) for v in arms) for arms in cols.values())
        print(f'\n=== {muscle} {p_dof} x {s_dof} ({label}) ===')
        print(f'surface peak |arm| = {surf_peak:.5f} m')
        # per-column: count events vs column peak and vs surface peak
        n_col_events, n_surf_events = 0, 0
        worst_col = None
        surf_event_locs = {}
        for s_deg, arms in cols.items():
            cpeak = max(abs(v) for v in arms)
            for i in range(len(arms) - 1):
                d = abs(arms[i + 1] - arms[i])
                loc = round((p_degs[i] + p_degs[i + 1]) / 2.0, 2)
                if cpeak > 0 and d / cpeak * 100.0 > 2.0:
                    n_col_events += 1
                    if worst_col is None or d / cpeak > worst_col[0]:
                        worst_col = (d / cpeak * 100.0, s_deg, loc,
                                     cpeak, d)
                if surf_peak > 0 and d / surf_peak * 100.0 > 2.0:
                    n_surf_events += 1
                    surf_event_locs.setdefault(s_deg, []).append(loc)
        print(f'events vs COLUMN peak (>2%): {n_col_events}')
        if worst_col:
            pct, s_deg, loc, cpeak, d = worst_col
            print(f'  worst: {pct:.2f}% of column peak at secondary='
                  f'{s_deg:g}, primary={loc:g} '
                  f'(col peak {cpeak:.5f}, |d| {d:.5f})')
        print(f'events vs SURFACE peak (>2%): {n_surf_events}')
        for s_deg in list(surf_event_locs)[:6]:
            locs = surf_event_locs[s_deg]
            print(f'  secondary={s_deg:g}: {len(locs)} at primary '
                  + ','.join(str(x) for x in locs[:14]))
        if len(surf_event_locs) > 6:
            print(f'  ... {len(surf_event_locs)} secondary values total')

    print('\n=== tendon-slack contours (l_mt <= tsl) ===')
    for muscle, p_dof, s_dof in SLACK_CASES:
        state = model.initSystem()
        coords = model.getCoordinateSet()
        mu = model.getMuscles().get(muscle)
        tsl = mu.get_tendon_slack_length()
        pcoord = coords.get(p_dof)
        p_degs = angles_for(*rom[p_dof], step=1.0)
        s_degs = angles_for(*rom[s_dof], step=2.0)
        n_slack, first, last = 0, None, None
        s_lo_hit, s_hi_hit = None, None
        for s_deg in s_degs:
            coords.get(s_dof).setValue(state, s_deg * DEG2RAD)
            row_hit = []
            for p_deg in p_degs:
                pcoord.setValue(state, p_deg * DEG2RAD)
                model.realizePosition(state)
                if mu.getLength(state) <= tsl:
                    n_slack += 1
                    row_hit.append(p_deg)
                    if first is None:
                        first = (p_deg, s_deg)
                    last = (p_deg, s_deg)
            if row_hit:
                if s_lo_hit is None:
                    s_lo_hit = s_deg
                s_hi_hit = s_deg
        print(f'{muscle:10s} {p_dof:16s} x {s_dof:16s}: tsl={tsl:.3f} '
              f'slack cells={n_slack} of {len(p_degs)*len(s_degs)} '
              f'(2-deg secondary grid); first@{first} last@{last}; '
              f'secondary span with slack: {s_lo_hit}..{s_hi_hit}')

    print('\nDONE')


if __name__ == '__main__':
    main()
