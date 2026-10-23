"""Generate the 2-D SINGLE-MUSCLE human capability SURFACES for the
27-actuator capability table v2 (2026-10-23, easteregg2).

Implements FIX_SPEC_20261023.md sections 1, 1b and 2 with the flag layer of
opensim_limits_checklist_20261023.md. Supersedes NOTHING on disk: the v1
1-D CSVs (human_torque_curves_20261022*.csv) stay untouched; all outputs are
NEW files dated _20261023.

Surface definition (FIX_SPEC 1): per actuator, tau_h(theta_p, theta_s) =
sum over the CHOSEN single-muscle set (TARGET_RULE, Ben 2026-10-22, copied
verbatim from gen_human_torque_curves_singlemuscle_20261022.py) of
computeMomentArm(about the PRIMARY dof) x Fmax, on a 5x5 grading grid of
INCLUSIVE UNIFORM nodes over rom_table_deg for the TWO DOFs the route
crosses (the BiPam class's two hinges). Grading values are DIRECT OpenSim
evaluations at the nodes - never interpolated. All other DOFs of the
actuator's map entry are pinned at 0 deg and recorded. STANDING TRAP: human
muscles on STOCK gait2392_simbody routing only; no robot/BPA path enters
this computation anywhere.

Axis pairs:
  - map primary_dof + FIRST map secondary where the map lists one
    (1-9, 13-23, 25; remaining map secondaries pinned at 0);
  - the CLASS hinge pair for the uniarticular-mounted actuators whose map
    entry lists no secondary (FIX_SPEC 1: the class "even serves
    uniarticular muscles over multiple DOFs"):
      10 bifemsh, 11 vas_int -> knee_angle_r x hip_flexion_r
        (Bi_Hip class hinges; biPulleySpecsFromMasters.m:121 thetaHipKnee
         [-15,110] x [-120,10]),
      12 soleus -> ankle_angle_r x knee_angle_r
        (Calves class hinges; biPulleySpecsFromMasters.m:122
         thetaKneeAnkle),
      26 iliacus, 27 psoas -> hip_flexion_r x hip_adduction_r
        (Uni_Hip class: the two-hinge pair the class itself grades for
         glut_max (map first secondary) and add_mag (map only secondary);
         FIX_SPEC 1b's "any hip surface that grades the adduction axis"
         anticipates these).

Hip adduction coupling (FIX_SPEC 1b, Ben 2026-10-23; supersedes the static
rectangle AND the 10-22 placeholder): on every surface that grades
hip_adduction_r (24, 25, 26, 27) the graded domain is the COUPLED envelope
- adduction max = 20 - 10*(f+15)/125 deg, abduction max (negative) =
-10 - 50*(f+15)/125 deg, LINEAR in hip flexion f in [-15, 110] (bone
interference, femur-pelvis contact). A cell is in the graded domain iff its
center satisfies abd_min(f_c) <= a_c <= add_max(f_c) (closed). Cells
outside are EXCLUDED from the graded domain for BOTH rows (removed from the
area denominator, never failures); counts are reported. hip_rotation stays
the static rom_table_deg box and is PINNED AT 0 on actuator 24.

Flag layer (checklist 20261023; flagged cells STILL COUNT in the headline;
with- and without-flag numbers both reported; curves NEVER corrected):
  FLAG_VIA  - a ConditionalPathPoint isActive-set change inside the cell's
              angular bin (0.5 deg cross sweeps; boundary is a constant-
              angle line, checklist-verified).
  FLAG_STEP - adjacent-sample |d(arm)| > 2% of A* (A* = per-muscle max
              |arm| over the whole surface - SURFACE-WIDE, never
              per-column) inside the cell's bin (0.5 deg cross sweeps).
  FLAG_SLACK- any 1-deg full-box sample in the cell's bin with
              l_mt <= tendon_slack_length.
  FLAG_SHORT- any 1-deg full-box sample in the cell's bin with
              lceN_est < 0.5 (constant-width pennation estimate,
              lce = sqrt((l_mt-tsl)^2 + (lo sin phi0)^2)).
Cell bins: 5 inclusive nodes per axis; bin of node i = the midpoint-Voronoi
interval clipped to [lo, hi] (outer bins half-width).

Halt gates (checklist section 8): opensim != 4.4.1; model/map sha1 mismatch
(pins below); any live muscle wrap set non-empty; any gate-table flag of
checklist section 7 missing from the graded domain. Knee finite-difference
validation stays INFORMATIONAL only; every other primary DOF gates hard at
5% (v1 behavior).

Outputs (this folder, all NEW):
  human_surfaces_5x5_grading_20261023.csv   - 675 cell rows, single grading
        source of truth for BOTH capability-table rows
  human_surfaces_9x9_20261023/              - per-actuator 9x9 viewing
        surfaces, bilinearly interpolated from the 1-deg OpenSim box
  human_surfaces_manifest_20261023.json     - axes, ranges, pinning,
        envelope, detectors, gates, flag counts, provenance

Run (easteregg2):
  "D:/Anaconda/envs/myo/python.exe" tmp/paper_bridge/capability_table_20261022/gen_human_surfaces_2d_20261023.py
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
OUTDIR = os.path.join(REPO, 'tmp', 'paper_bridge',
                      'capability_table_20261022')
MODEL = os.path.join(REPO, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody',
                     'gait2392_simbody.osim')
MAP = os.path.join(REPO, 'Code', 'Matlab', 'Mesh_Optimization',
                   'Routing27_v2', 'actuator_map_27.json')
AUDIT = os.path.join(REPO, 'tmp', 'paper_bridge', 'opt_improve_20261022',
                     'path_artifact_audit_gait2392_simbody_'
                     '20261022_034030.csv')
OUT_CSV = os.path.join(OUTDIR, 'human_surfaces_5x5_grading_20261023.csv')
OUT_DIR9 = os.path.join(OUTDIR, 'human_surfaces_9x9_20261023')
OUT_MAN = os.path.join(OUTDIR, 'human_surfaces_manifest_20261023.json')

PINNED_VERSION = '4.4.1'
PIN_MODEL_SHA = 'b89062c763e6d8fa3d7651c0a4fd2820055d807c'
PIN_MAP_SHA = 'a97afffeeae85687810986273db6ba0151eba1a2'
DEG2RAD = math.pi / 180.0
N_NODE = 5          # inclusive nodes per axis (5x5 = 25 cells)
N_VIEW = 9          # viewing export nodes per axis
BOX_DEG = 1.0       # full-box sampling (SLACK/SHORT bins, A*, 9x9 source)
FINE_DEG = 0.5      # flag cross sweeps (VIA/STEP)
STEP_PCT = 2.0      # FLAG_STEP threshold, % of surface-wide A*
SHORT_LCEN = 0.5    # FLAG_SHORT threshold
FD_SPOTS = 8
FD_H_DEG = 0.5
FD_TOL = 0.05

# TARGET_RULE copied verbatim from gen_human_torque_curves_singlemuscle_
# 20261022.py (Ben 2026-10-22 ruling).
SAME = 'same-named single stock muscle (map single-member human_group)'
TARGET_RULE = {
    1:  (['ercspn_r'], SAME),
    2:  (['intobl_r'], SAME),
    3:  (['extobl_r'], SAME),
    4:  (['bifemlh_r'], SAME),
    5:  (['semimem_r'],
         'semimem_r ALONE - drop semiten (Ben 2026-10-22 explicit ruling)'),
    6:  (['sar_r'], SAME),
    7:  (['tfl_r'], SAME),
    8:  (['rect_fem_r'], SAME),
    9:  (['grac_r'], SAME),
    10: (['bifemsh_r'], SAME),
    11: (['vas_int_r'],
         'vas_int_r ALONE - drop vas_lat/vas_med (Ben 2026-10-22 ruling)'),
    12: (['soleus_r'], SAME),
    13: (['med_gas_r'], SAME),
    14: (['lat_gas_r'], SAME),
    15: (['tib_post_r'],
         'tib_post_r (the real tibialis-posterior actuator; map ruling)'),
    16: (['tib_ant_r'],
         'tib_ant_r (map duplicate-name resolution: this actuator is '
         'tibialis anterior)'),
    17: (['per_brev_r'], SAME),
    18: (['per_long_r'], SAME),
    19: (['per_tert_r'], SAME),
    20: (['flex_dig_r'], SAME),
    21: (['flex_hal_r'], SAME),
    22: (['ext_dig_r'], SAME),
    23: (['ext_hal_r'], SAME),
    24: (['glut_max1_r', 'glut_max2_r', 'glut_max3_r'],
         'glut_max1_r + glut_max2_r + glut_max3_r ONLY - one anatomical '
         'muscle the model splits into three compartments; drop every '
         'med/min compartment (Ben 2026-10-22 ruling)'),
    25: (['add_mag1_r', 'add_mag2_r', 'add_mag3_r'],
         'add_mag1_r + add_mag2_r + add_mag3_r ONLY - compartment split of '
         'one anatomical muscle; drop add_long/add_brev (Ben 2026-10-22 '
         'ruling)'),
    26: (['iliacus_r'], SAME),
    27: (['psoas_r'], SAME),
}

# Surface axes: id -> (primary, secondary, pinned-at-0 list, provenance).
MAP_1ST = 'map primary_dof + FIRST map secondary (map ordering); '
AXES = {
    1:  ('lumbar_extension', 'lumbar_bending', ['lumbar_rotation'],
         MAP_1ST + 'lumbar_rotation (2nd map secondary) pinned at 0'),
    2:  ('lumbar_extension', 'lumbar_bending', [], MAP_1ST + 'exact'),
    3:  ('lumbar_extension', 'lumbar_bending', [], MAP_1ST + 'exact'),
    4:  ('knee_angle_r', 'hip_flexion_r', [], MAP_1ST + 'exact'),
    5:  ('knee_angle_r', 'hip_flexion_r', [], MAP_1ST + 'exact'),
    6:  ('knee_angle_r', 'hip_flexion_r', [], MAP_1ST + 'exact'),
    7:  ('knee_angle_r', 'hip_flexion_r', ['hip_adduction_r'],
         MAP_1ST + 'hip_adduction_r (2nd map secondary) pinned at 0 '
         '(feasible at every graded flexion under the 1b envelope: 0 is '
         'inside [abd_min, add_max] for all f)'),
    8:  ('knee_angle_r', 'hip_flexion_r', [], MAP_1ST + 'exact'),
    9:  ('knee_angle_r', 'hip_flexion_r', [], MAP_1ST + 'exact'),
    10: ('knee_angle_r', 'hip_flexion_r', [],
         'Bi_Hip class hinge pair (biPulleySpecsFromMasters.m thetaHipKnee '
         '[-15,110]x[-120,10]); map lists no secondary - the class serves '
         'the uniarticular muscle over its two hinges (FIX_SPEC 1)'),
    11: ('knee_angle_r', 'hip_flexion_r', [],
         'Bi_Hip class hinge pair (biPulleySpecsFromMasters.m '
         'thetaHipKnee); map lists no secondary - class two-hinge serving '
         '(FIX_SPEC 1)'),
    12: ('ankle_angle_r', 'knee_angle_r', [],
         'Calves class hinge pair (biPulleySpecsFromMasters.m '
         'thetaKneeAnkle); map lists no secondary - class two-hinge '
         'serving (FIX_SPEC 1); human soleus crosses only the ankle so '
         'the surface is honest-flat in knee'),
    13: ('ankle_angle_r', 'knee_angle_r', [], MAP_1ST + 'exact'),
    14: ('ankle_angle_r', 'knee_angle_r', [], MAP_1ST + 'exact'),
    15: ('subtalar_angle_r', 'ankle_angle_r', [], MAP_1ST + 'exact'),
    16: ('ankle_angle_r', 'subtalar_angle_r', [], MAP_1ST + 'exact'),
    17: ('subtalar_angle_r', 'ankle_angle_r', [], MAP_1ST + 'exact'),
    18: ('subtalar_angle_r', 'ankle_angle_r', [], MAP_1ST + 'exact'),
    19: ('subtalar_angle_r', 'ankle_angle_r', [], MAP_1ST + 'exact'),
    20: ('mtp_angle_r', 'ankle_angle_r', ['subtalar_angle_r'],
         MAP_1ST + 'subtalar_angle_r (2nd map secondary) pinned at 0'),
    21: ('mtp_angle_r', 'ankle_angle_r', ['subtalar_angle_r'],
         MAP_1ST + 'subtalar_angle_r pinned at 0'),
    22: ('mtp_angle_r', 'ankle_angle_r', ['subtalar_angle_r'],
         MAP_1ST + 'subtalar_angle_r pinned at 0'),
    23: ('mtp_angle_r', 'ankle_angle_r', ['subtalar_angle_r'],
         MAP_1ST + 'subtalar_angle_r pinned at 0'),
    24: ('hip_flexion_r', 'hip_adduction_r', ['hip_rotation_r'],
         MAP_1ST + 'hip_rotation_r (2nd map secondary) pinned at 0 '
         '(rotation stays the static rom box per FIX_SPEC 1b); adduction '
         'axis is envelope-coupled'),
    25: ('hip_adduction_r', 'hip_flexion_r', [],
         MAP_1ST + 'exact; adduction (primary) axis is envelope-coupled'),
    26: ('hip_flexion_r', 'hip_adduction_r', [],
         'Uni_Hip class two-hinge pair hip_flexion x hip_adduction (the '
         'pair the class grades for glut_max (map first secondary) and '
         'add_mag (map only secondary)); map lists no secondary - class '
         'two-hinge serving (FIX_SPEC 1); adduction axis envelope-coupled '
         '(1b "any hip surface that grades the adduction axis")'),
    27: ('hip_flexion_r', 'hip_adduction_r', [],
         'Uni_Hip class two-hinge pair (see actuator 26); adduction axis '
         'envelope-coupled'),
}

# Checklist section 7 gate table - hard gates only ("yes"/named boundaries;
# borderline and marginal entries deliberately NOT gated).
GATE_FLAGS = {
    4:  ['SHORT'],
    5:  ['VIA', 'SLACK', 'SHORT'],
    7:  ['SLACK', 'SHORT'],
    8:  ['VIA', 'SHORT'],
    9:  ['VIA', 'STEP'],
    11: ['VIA'],
    13: ['VIA', 'SLACK', 'SHORT'],
    14: ['VIA', 'SLACK', 'SHORT'],
    15: ['STEP'],
    17: ['STEP'],
    18: ['STEP'],
    24: ['SHORT'],
    25: ['SHORT'],
    26: ['VIA', 'STEP', 'SHORT'],
    27: ['VIA', 'STEP', 'SHORT'],
}

# FIX_SPEC 1b coupled envelope (adduction positive; LINEAR in flexion).
ENV_FLEX_LO, ENV_FLEX_HI = -15.0, 110.0


def env_add_max(f):
    return 20.0 - 10.0 * (f - ENV_FLEX_LO) / (ENV_FLEX_HI - ENV_FLEX_LO)


def env_abd_min(f):
    return -10.0 - 50.0 * (f - ENV_FLEX_LO) / (ENV_FLEX_HI - ENV_FLEX_LO)


def sha1_of(path):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 16), b''):
            h.update(chunk)
    return h.hexdigest()


def stat_sig(path):
    st = os.stat(path)
    return (st.st_mtime_ns, st.st_size)


def angles_for(lo, hi, step):
    n = int(round((hi - lo) / step)) + 1
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def nodes_for(lo, hi, n):
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def bin_edges(nodes, lo, hi):
    """Midpoint-Voronoi bin edges of the inclusive node set."""
    edges = [lo]
    for i in range(len(nodes) - 1):
        edges.append((nodes[i] + nodes[i + 1]) / 2.0)
    edges.append(hi)
    return edges


def bin_index(x, edges):
    """Bin index of x under closed bins [e_i, e_{i+1}]; boundary hits flag
    BOTH adjacent bins (conservative)."""
    hits = []
    for i in range(len(edges) - 1):
        if edges[i] <= x <= edges[i + 1]:
            hits.append(i)
    return hits


def short_name(socket_text):
    return (socket_text or '').strip().split('/')[-1]


def cond_moving_from_xml(model_path, muscle_names):
    """{muscle: {'cond': [(coord, lo, hi)...], 'n_moving': int}}."""
    tree = ET.parse(model_path)
    out = {m: {'cond': [], 'n_moving': 0} for m in muscle_names}
    for el in tree.getroot().iter():
        tag = el.tag.split('}')[-1]
        if not tag.endswith('Muscle') or el.get('name') not in out:
            continue
        name = el.get('name')
        for sub in el.iter():
            stag = sub.tag.split('}')[-1]
            if stag == 'MovingPathPoint':
                out[name]['n_moving'] += 1
            elif stag == 'ConditionalPathPoint':
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
                out[name]['cond'].append((coord, rng[0], rng[1]))
    return out


def load_audit(path):
    rows, comments = [], []
    with open(path, newline='') as f:
        for ln in f:
            if ln.startswith('#'):
                comments.append(ln.strip())
                continue
            break
    with open(path, newline='') as f:
        rdr = csv.DictReader((ln for ln in f if not ln.startswith('#')))
        rows.extend(rdr)
    return rows, comments


def caveat_for(rows, act_id, dof, chosen):
    """Primary-DOF 1-D artifact caveat string (v1 machinery, verbatim)."""
    flagged = [r for r in rows
               if r['actuator_id'] == str(act_id) and r['dof'] == dof
               and r['kind'] == 'primary' and r['muscle'] in chosen
               and r['artifact_flag'] == '1']
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
    """arm vs central-difference -dl/dtheta (v1 machinery, re-run live)."""
    results = {}
    spots = [degs[int(round(i * (len(degs) - 1) / (FD_SPOTS - 1)))]
             for i in range(FD_SPOTS)]
    for m in muscles:
        mu = model.getMuscles().get(m)
        worst, checked = 0.0, 0
        for deg in spots:
            coord.setValue(state, (deg - FD_H_DEG) * DEG2RAD)
            model.realizePosition(state)
            lm = mu.getLength(state)
            coord.setValue(state, (deg + FD_H_DEG) * DEG2RAD)
            model.realizePosition(state)
            lp = mu.getLength(state)
            dl_dth = -(lp - lm) / (2.0 * FD_H_DEG * DEG2RAD)
            coord.setValue(state, deg * DEG2RAD)
            model.realizePosition(state)
            arm = mu.computeMomentArm(state, coord)
            if abs(dl_dth) < 1e-4:
                continue
            worst = max(worst, abs(arm - dl_dth) / max(abs(dl_dth), 1e-9))
            checked += 1
        results[m] = worst if checked else 0.0
    return results


def bilinear(box, p_degs, s_degs, tp, ts):
    """box[i][j] = value at (p_degs[i], s_degs[j]); interpolate at (tp, ts)."""

    def frac(vals, t):
        if t <= vals[0]:
            return 0, 0.0
        if t >= vals[-1]:
            return len(vals) - 2, 1.0
        for i in range(len(vals) - 1):
            if vals[i] <= t <= vals[i + 1]:
                w = (t - vals[i]) / (vals[i + 1] - vals[i]) \
                    if vals[i + 1] > vals[i] else 0.0
                return i, w
        return len(vals) - 2, 1.0

    ip, wp = frac(p_degs, tp)
    js, ws = frac(s_degs, ts)
    v00 = box[ip][js]
    v01 = box[ip][js + 1]
    v10 = box[ip + 1][js]
    v11 = box[ip + 1][js + 1]
    return ((1 - wp) * (1 - ws) * v00 + (1 - wp) * ws * v01
            + wp * (1 - ws) * v10 + wp * ws * v11)


def main():
    ver = str(getattr(osim, '__version__', 'unknown'))
    print(f'opensim : {ver} (pinned {PINNED_VERSION})')
    if ver != PINNED_VERSION:
        sys.exit(f'OPENSIM VERSION PIN VIOLATION: {ver} != {PINNED_VERSION}')
    msha, psha = sha1_of(MODEL), sha1_of(MAP)
    print(f'model   : {MODEL}\n  sha1  : {msha}')
    print(f'map     : {MAP}\n  sha1  : {psha}')
    if msha != PIN_MODEL_SHA:
        sys.exit(f'MODEL SHA PIN VIOLATION: {msha}')
    if psha != PIN_MAP_SHA:
        sys.exit(f'MAP SHA PIN VIOLATION: {psha}')
    audit_rows, audit_comments = load_audit(AUDIT)
    hdr = ' '.join(audit_comments)
    for h in (msha, psha):
        if h not in hdr:
            sys.exit(f'audit CSV does not pin current file sha1 {h}')
    print('audit pins model+map sha1: OK')
    sigs = {p: stat_sig(p) for p in (MODEL, MAP, AUDIT)}

    with open(MAP) as f:
        cmap = json.load(f)
    rom = cmap['rom_table_deg']
    by_id = {a['id']: a for a in cmap['actuators']}
    if sorted(TARGET_RULE) != sorted(by_id) or sorted(AXES) != sorted(by_id):
        sys.exit('TARGET_RULE/AXES ids do not match the map actuators')

    model = osim.Model(MODEL)
    state = model.initSystem()
    coords = model.getCoordinateSet()
    coord_names = [coords.get(i).getName() for i in range(coords.getSize())]
    all_muscles = [model.getMuscles().get(i).getName()
                   for i in range(model.getMuscles().getSize())]

    # Existence checks + default-pose claim (v1: every map DOF defaults 0).
    for a in cmap['actuators']:
        aid = a['id']
        for m in TARGET_RULE[aid][0]:
            if m not in all_muscles:
                sys.exit(f'target muscle {m} missing in model')
        p, s, pinned, _ = AXES[aid]
        for d in [p, s] + pinned + a['secondary_dofs']:
            if d not in coord_names:
                sys.exit(f'dof {d} missing in model')
    for d in sorted(set(sum(([AXES[i][0], AXES[i][1]] + AXES[i][2]
                             for i in AXES), []))):
        dv = coords.get(d).getValue(state) / DEG2RAD
        if abs(dv) > 1e-9:
            sys.exit(f'coordinate {d} default is {dv} deg, expected 0')

    all_chosen = sorted({m for v in TARGET_RULE.values() for m in v[0]})
    xmlinfo = cond_moving_from_xml(MODEL, all_chosen)

    # Halt gate: live wrap sets must be EMPTY on every chosen muscle
    # (checklist L3; stock model has zero PathWrap entries).
    for m in all_chosen:
        mu = model.getMuscles().get(m)
        try:
            nwrap = mu.getGeometryPath().getWrapSet().getSize()
        except Exception as e:
            sys.exit(f'wrap-set check failed on {m}: {e}')
        if nwrap != 0:
            sys.exit(f'WRAP HALT GATE: {m} has {nwrap} live wrap objects')
    print('wrap halt gate: all chosen muscles have empty wrap sets')

    for aid in sorted(AXES):
        for axis in (AXES[aid][0], AXES[aid][1]):
            lo, hi = rom[axis]
            if abs((hi - lo) - round(hi - lo)) > 1e-9:
                sys.exit(f'axis {axis} span {hi - lo} not integer degrees')

    csv_rows = []
    man_acts = []
    flagged_acts = []
    gate_failures = []
    fd_failures = []
    tau_boxes = {}      # aid -> tau box on the 1-deg lattice (for 9x9)
    for aid in sorted(AXES):
        a = by_id[aid]
        chosen, ruling = TARGET_RULE[aid]
        p_dof, s_dof, pinned, provenance = AXES[aid]
        p_lo, p_hi = (float(x) for x in rom[p_dof])
        s_lo, s_hi = (float(x) for x in rom[s_dof])
        envelope = 'hip_adduction_r' in (p_dof, s_dof)

        nodes_p = nodes_for(p_lo, p_hi, N_NODE)
        nodes_s = nodes_for(s_lo, s_hi, N_NODE)
        edges_p = bin_edges(nodes_p, p_lo, p_hi)
        edges_s = bin_edges(nodes_s, s_lo, s_hi)
        box_p = angles_for(p_lo, p_hi, BOX_DEG)
        box_s = angles_for(s_lo, s_hi, BOX_DEG)

        state = model.initSystem()
        pcoord = coords.get(p_dof)
        scoord = coords.get(s_dof)
        for d in pinned:
            coords.get(d).setValue(state, 0.0)

        props = {}
        for m in chosen:
            mu = model.getMuscles().get(m)
            props[m] = {
                'tsl': mu.get_tendon_slack_length(),
                'lo': mu.get_optimal_fiber_length(),
                'fmax': mu.get_max_isometric_force(),
                'phi0': mu.get_pennation_angle_at_optimal(),
                'mu': mu,
            }

        def set_pose(pd, sd):
            pcoord.setValue(state, pd * DEG2RAD)
            scoord.setValue(state, sd * DEG2RAD)
            model.realizePosition(state)

        # ---------- 1-deg full box: arms, lengths, A*, slack/short bins --
        arm_box = {m: [[0.0] * len(box_s) for _ in box_p] for m in chosen}
        lmt_box = {m: [[0.0] * len(box_s) for _ in box_p] for m in chosen}
        lcen_box = {m: [[0.0] * len(box_s) for _ in box_p] for m in chosen}
        for ip, pd in enumerate(box_p):
            for js, sd in enumerate(box_s):
                set_pose(pd, sd)
                gp, sc = (pcoord.getValue(state) / DEG2RAD,
                          scoord.getValue(state) / DEG2RAD)
                assert abs(gp - pd) < 1e-6 and abs(sc - sd) < 1e-6, \
                    f'coordinate clamped: ({gp},{sc}) != ({pd},{sd})'
                for m in chosen:
                    pr = props[m]
                    arm = pr['mu'].computeMomentArm(state, pcoord)
                    lmt = pr['mu'].getLength(state)
                    h = pr['lo'] * math.sin(pr['phi0'])
                    lcen = math.sqrt(
                        max((lmt - pr['tsl']) ** 2 + h * h, 0.0)) / pr['lo']
                    arm_box[m][ip][js] = arm
                    lmt_box[m][ip][js] = lmt
                    lcen_box[m][ip][js] = lcen

        # ---------- 0.5-deg cross sweeps: VIA + STEP --------------------
        fine_p = angles_for(p_lo, p_hi, FINE_DEG)
        fine_s = angles_for(s_lo, s_hi, FINE_DEG)

        def sweep(seq_dofs, seq, fixed_dof, fixed_val):
            """Sweep seq_dofs over seq at fixed_dof=fixed_val; per chosen
            muscle return (arms, isActive-patterns)."""
            coords.get(fixed_dof).setValue(state, fixed_val * DEG2RAD)
            out = {m: ([], []) for m in chosen}
            cobj = coords.get(seq_dofs)
            for deg in seq:
                cobj.setValue(state, deg * DEG2RAD)
                model.realizePosition(state)
                for m in chosen:
                    pr = props[m]
                    out[m][0].append(
                        pr['mu'].computeMomentArm(state, pcoord))
                    pts = pr['mu'].getGeometryPath().getPathPointSet()
                    out[m][1].append(tuple(
                        pts.get(i).isActive(state)
                        for i in range(pts.getSize())))
            return out

        via_cells = set()      # (i, j) node indices
        step_cells = set()
        flag_detail = {}       # (i, j) -> set of tokens

        def add_flag(i, j, token, why):
            key = (i, j)
            if token == 'FLAG_VIA':
                via_cells.add(key)
            elif token == 'FLAG_STEP':
                step_cells.add(key)
            flag_detail.setdefault(key, set()).add(f'{token}:{why}')

        # primary-axis sweeps at each secondary node
        for j, snode in enumerate(nodes_s):
            res = sweep(p_dof, fine_p, s_dof, snode)
            for m in chosen:
                arms, pats = res[m]
                astar = max(abs(v) for r in arm_box[m] for v in r)
                for k in range(len(arms) - 1):
                    d = abs(arms[k + 1] - arms[k])
                    mid = (fine_p[k] + fine_p[k + 1]) / 2.0
                    if pats[k] != pats[k + 1]:
                        for i in bin_index(mid, edges_p):
                            add_flag(i, j, 'FLAG_VIA', m)
                    if astar > 0 and d / astar * 100.0 > STEP_PCT:
                        for i in bin_index(mid, edges_p):
                            add_flag(i, j, 'FLAG_STEP', m)
        # secondary-axis sweeps at each primary node
        for i, pnode in enumerate(nodes_p):
            res = sweep(s_dof, fine_s, p_dof, pnode)
            for m in chosen:
                arms, pats = res[m]
                astar = max(abs(v) for r in arm_box[m] for v in r)
                for k in range(len(arms) - 1):
                    d = abs(arms[k + 1] - arms[k])
                    mid = (fine_s[k] + fine_s[k + 1]) / 2.0
                    if pats[k] != pats[k + 1]:
                        for j in bin_index(mid, edges_s):
                            add_flag(i, j, 'FLAG_VIA', m)
                    if astar > 0 and d / astar * 100.0 > STEP_PCT:
                        for j in bin_index(mid, edges_s):
                            add_flag(i, j, 'FLAG_STEP', m)

        # ---------- slack / short bins from the 1-deg box ----------------
        slack_cells, short_cells = set(), set()
        for m in chosen:
            pr = props[m]
            for ip, pd in enumerate(box_p):
                pb = bin_index(pd, edges_p)
                for js, sd in enumerate(box_s):
                    if lmt_box[m][ip][js] <= pr['tsl']:
                        for i in pb:
                            for j in bin_index(sd, edges_s):
                                slack_cells.add((i, j))
                                flag_detail.setdefault((i, j), set()).add(
                                    f'FLAG_SLACK:{m}')
                    if lcen_box[m][ip][js] < SHORT_LCEN:
                        for i in pb:
                            for j in bin_index(sd, edges_s):
                                short_cells.add((i, j))
                                flag_detail.setdefault((i, j), set()).add(
                                    f'FLAG_SHORT:{m}')

        # ---------- grading nodes: direct OpenSim evaluation ------------
        tau_nodes = [[0.0] * N_NODE for _ in range(N_NODE)]
        in_domain = [[False] * N_NODE for _ in range(N_NODE)]
        for i, pd in enumerate(nodes_p):
            for j, sd in enumerate(nodes_s):
                if envelope:
                    f = pd if p_dof == 'hip_flexion_r' else sd
                    add_ = sd if p_dof == 'hip_flexion_r' else pd
                    in_domain[i][j] = (env_abd_min(f) <= add_
                                       <= env_add_max(f))
                else:
                    in_domain[i][j] = True
                set_pose(pd, sd)
                for m in chosen:
                    pr = props[m]
                    tau_nodes[i][j] += (
                        pr['mu'].computeMomentArm(state, pcoord)
                        * pr['fmax'])

        # tau box for the 9x9 viewing interpolation (from cached arms)
        tau_box = [[sum(arm_box[m][ip][js] * props[m]['fmax']
                        for m in chosen) for js in range(len(box_s))]
                   for ip in range(len(box_p))]

        # ---------- FD validation along primary at secondary = 0 --------
        scoord.setValue(state, 0.0)
        model.realizePosition(state)
        fd = fd_check(model, state, pcoord, chosen, fine_p)
        knee = p_dof == 'knee_angle_r'
        bad = {} if knee else {m: v for m, v in fd.items() if v > FD_TOL}
        if bad:
            fd_failures.append((a['master_name'], bad))

        # ---------- envelope cell accounting -----------------------------
        n_graded = sum(1 for i in range(N_NODE) for j in range(N_NODE)
                       if in_domain[i][j])
        n_excluded = N_NODE * N_NODE - n_graded

        # ---------- gate check on the GRADED domain ----------------------
        present = set()
        for (i, j) in via_cells | step_cells | slack_cells | short_cells:
            if in_domain[i][j]:
                if (i, j) in via_cells:
                    present.add('VIA')
                if (i, j) in step_cells:
                    present.add('STEP')
                if (i, j) in slack_cells:
                    present.add('SLACK')
                if (i, j) in short_cells:
                    present.add('SHORT')
        missing = [g for g in GATE_FLAGS.get(aid, []) if g not in present]
        if missing:
            gate_failures.append((aid, a['master_name'], missing))

        # ---------- audit caveat + bookkeeping ---------------------------
        caveat = caveat_for(audit_rows, aid, p_dof, chosen)
        fsum = sum(props[m]['fmax'] for m in chosen)
        n_flagged_graded = 0
        counts = {'VIA': 0, 'STEP': 0, 'SLACK': 0, 'SHORT': 0}
        for i in range(N_NODE):
            for j in range(N_NODE):
                cell_flags = []
                if in_domain[i][j]:
                    if (i, j) in via_cells:
                        cell_flags.append('FLAG_VIA')
                        counts['VIA'] += 1
                    if (i, j) in step_cells:
                        cell_flags.append('FLAG_STEP')
                        counts['STEP'] += 1
                    if (i, j) in slack_cells:
                        cell_flags.append('FLAG_SLACK')
                        counts['SLACK'] += 1
                    if (i, j) in short_cells:
                        cell_flags.append('FLAG_SHORT')
                        counts['SHORT'] += 1
                    if cell_flags:
                        n_flagged_graded += 1
                # every cell is listed (excluded cells carry
                # in_graded_domain=0 and never enter counts/denominators)
                csv_rows.append([
                    aid, a['master_name'], '+'.join(chosen),
                    p_dof, f'{p_lo:g}', f'{p_hi:g}',
                    s_dof, f'{s_lo:g}', f'{s_hi:g}',
                    ';'.join(pinned) if pinned else '(none)',
                    'yes' if envelope else 'no',
                    f'{nodes_p[i]:.6g}', f'{nodes_s[j]:.6g}',
                    '1' if in_domain[i][j] else '0',
                    f'{tau_nodes[i][j]:.9f}', f'{fsum:g}',
                    '1' if (i, j) in via_cells else '0',
                    '1' if (i, j) in step_cells else '0',
                    '1' if (i, j) in slack_cells else '0',
                    '1' if (i, j) in short_cells else '0',
                    '; '.join(sorted(flag_detail.get((i, j), [])))
                    or 'none',
                    caveat])
        if n_flagged_graded:
            flagged_acts.append(a['master_name'])
        tau_boxes[aid] = (tau_box, box_p, box_s)

        box_stats = {}
        for m in chosen:
            minl, at = None, None
            minc, atc = None, None
            maxa = 0.0
            for ip, pd in enumerate(box_p):
                for js, sd in enumerate(box_s):
                    if minl is None or lmt_box[m][ip][js] < minl:
                        minl, at = lmt_box[m][ip][js], (pd, sd)
                    if minc is None or lcen_box[m][ip][js] < minc:
                        minc, atc = lcen_box[m][ip][js], (pd, sd)
                    maxa = max(maxa, abs(arm_box[m][ip][js]))
            box_stats[m] = {
                'min_lmt_m': round(minl, 5),
                'min_lmt_at_deg': f'{at[0]:g};{at[1]:g}',
                'min_lcen_est': round(minc, 4),
                'min_lcen_est_at_deg': f'{atc[0]:g};{atc[1]:g}',
                'max_abs_arm_m': round(maxa, 5),
                'fmax_n': round(props[m]['fmax'], 1),
                'tsl_m': round(props[m]['tsl'], 5),
                'lo_m': round(props[m]['lo'], 5),
                'cond_points_xml': ';'.join(
                    f'{c}:{lo_:.2f}..{hi_:.2f}'
                    for c, lo_, hi_ in xmlinfo[m]['cond']),
                'n_moving_pp': xmlinfo[m]['n_moving'],
            }

        man_acts.append({
            'id': aid, 'master_name': a['master_name'],
            'chosen_muscles': chosen, 'target_rule_ruling': ruling,
            'axes': {
                'primary': p_dof, 'primary_rom_deg': [p_lo, p_hi],
                'secondary': s_dof, 'secondary_rom_deg': [s_lo, s_hi],
                'pinned_at_0_deg': pinned,
                'provenance': provenance,
                'map_secondary_dofs': a['secondary_dofs'],
            },
            'envelope': {
                'applies': envelope,
                'rule': ('coupled adduction envelope (FIX_SPEC 1b): '
                         'add_max(f) = 20 - 10*(f+15)/125 deg, abd_min(f) '
                         '= -10 - 50*(f+15)/125 deg, LINEAR in flexion f '
                         '[-15,110]; cell graded iff center inside (closed)'
                         ) if envelope else
                        'rom_table_deg rectangle (no hip_adduction axis)',
                'n_graded_cells': n_graded,
                'n_excluded_cells': n_excluded,
            },
            'grid': {
                'nodes_primary_deg': [round(x, 6) for x in nodes_p],
                'nodes_secondary_deg': [round(x, 6) for x in nodes_s],
                'bins_primary_deg': [round(x, 6) for x in edges_p],
                'bins_secondary_deg': [round(x, 6) for x in edges_s],
            },
            'chosen_fmax_sum_n': round(fsum, 1),
            'tau_h_definition': ('sum over chosen muscles of '
                                 'computeMomentArm(about primary) x Fmax, '
                                 'direct OpenSim evaluation at each node'),
            'peak_abs_tau_h_Nm': round(
                max(abs(v) for r in tau_nodes for v in r), 6),
            'flag_counts_graded_cells': counts,
            'n_flagged_graded_cells': n_flagged_graded,
            'n_graded_cells': n_graded,
            'gate_table_expected': GATE_FLAGS.get(aid, []),
            'gate_table_missing': missing,
            'fd_worst_rel': {m: round(v, 5) for m, v in fd.items()},
            'fd_status': 'informational (knee)' if knee else
                         ('PASS' if not bad else 'FAIL'),
            'audit_caveat_1d_primary': caveat,
            'box_stats_1deg': box_stats,
        })
        print(f'[{aid:2d}] {a["master_name"][:30]:30s} '
              f'{p_dof:17s} x {s_dof:17s} '
              f'graded {n_graded:2d}/25 env-excl {n_excluded:2d}  '
              f'flags V{counts["VIA"]:2d} S{counts["STEP"]:2d} '
              f'K{counts["SLACK"]:2d} H{counts["SHORT"]:2d}  '
              f'gate {"MISSING " + ",".join(missing) if missing else "ok"}  '
              f'FD {man_acts[-1]["fd_status"]}')

    # ---------------- halt gates ----------------------------------------
    if gate_failures:
        for aid, name, miss in gate_failures:
            print(f'GATE FAILURE act {aid} ({name}): missing {miss}')
        sys.exit('checklist section 7 gate table not satisfied - halting, '
                 'no outputs written')
    if fd_failures:
        for name, bad in fd_failures:
            print(f'FD FAILURE {name}: {bad}')
        sys.exit('moment-arm finite-difference validation failed - '
                 'no outputs written')
    for p, s0 in sigs.items():
        if stat_sig(p) != s0:
            sys.exit(f'READ-ONLY VIOLATION: {p} changed during the run')

    # ---------------- write outputs -------------------------------------
    os.makedirs(OUT_DIR9, exist_ok=True)
    with open(OUT_CSV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['# human_surfaces_5x5_grading_20261023 - 2-D '
                    'SINGLE-MUSCLE human capability surfaces; SINGLE '
                    'grading source of truth for BOTH capability-table '
                    'rows (hand-tuned and optimized); FIX_SPEC_20261023 '
                    'sections 1/1b/2; flags per '
                    'opensim_limits_checklist_20261023'])
        w.writerow([f'# opensim_version {ver} (pinned {PINNED_VERSION}); '
                    f'python {sys.executable}'])
        w.writerow([f'# model {MODEL} sha1 {msha} (STOCK routing; human '
                    f'muscles never routed on robot/BPA paths)'])
        w.writerow([f'# map {MAP} sha1 {psha}'])
        w.writerow([f'# audit {AUDIT} sha1 {sha1_of(AUDIT)} (1-D caveat '
                    f'column source)'])
        w.writerow(['# grid: 5x5 INCLUSIVE uniform nodes per axis over '
                    'rom_table_deg (5 nodes per axis; 25 cells per '
                    'actuator; cell value = DIRECT OpenSim evaluation at '
                    'the node; cell bin = midpoint-Voronoi within the '
                    'axis range)'])
        w.writerow(['# tau_h_Nm = sum over chosen muscles of '
                    'computeMomentArm(about PRIMARY dof) x Fmax (OpenSim '
                    'native signs; arm = -dL/dq so tau>0 = positive '
                    'coordinate direction; knee flexion NEGATIVE)'])
        w.writerow(['# in_graded_domain=0 ONLY on hip-envelope surfaces '
                    '(hip_adduction coupled to flexion; FIX_SPEC 1b): '
                    'excluded cells are OUTSIDE the area denominator for '
                    'BOTH rows - never counted as failures'])
        w.writerow(['# flags: FLAG_VIA / FLAG_STEP (0.5-deg cross sweeps; '
                    'STEP threshold 2% of per-muscle surface-wide peak) / '
                    'FLAG_SLACK / FLAG_SHORT (1-deg box bins; lceN_est '
                    '< 0.5); FLAGGED CELLS STILL COUNT in the headline '
                    'percentage; with- and without-flag numbers are both '
                    'reported; curves are NEVER corrected'])
        w.writerow(['# pinned DOFs are held at 0 deg; non-surface joints '
                    'of 3+-joint routes are pinned at 0 and named in '
                    'pinned_dofs_at_0'])
        w.writerow([f'# generated {datetime.datetime.now().isoformat()}'])
        w.writerow(['actuator_id', 'master_name', 'chosen_muscles',
                    'primary_dof', 'primary_rom_lo_deg', 'primary_rom_hi_deg',
                    'secondary_dof', 'secondary_rom_lo_deg',
                    'secondary_rom_hi_deg', 'pinned_dofs_at_0',
                    'envelope_limited', 'p_node_deg', 's_node_deg',
                    'in_graded_domain', 'tau_h_Nm', 'fmax_sum_n',
                    'flag_via', 'flag_step', 'flag_slack', 'flag_short',
                    'flags_detail', 'audit_caveat_1d_primary'])
        for r in csv_rows:
            w.writerow(r)
    print(f'\nwrote {OUT_CSV} ({len(csv_rows)} cell rows)')

    for entry in man_acts:
        aid = entry['id']
        p_dof = entry['axes']['primary']
        s_dof = entry['axes']['secondary']
        p_lo, p_hi = entry['axes']['primary_rom_deg']
        s_lo, s_hi = entry['axes']['secondary_rom_deg']
        chosen = entry['chosen_muscles']
        tau_box, box_p, box_s = tau_boxes[aid]
        vp = nodes_for(p_lo, p_hi, N_VIEW)
        vs = nodes_for(s_lo, s_hi, N_VIEW)
        envelope = entry['envelope']['applies']
        fname = (f'act{aid:02d}_'
                 + chosen[0].replace('_r', '').replace('-', '')
                 + '_9x9_20261023.csv')
        path9 = os.path.join(OUT_DIR9, fname)
        with open(path9, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow([f'# 9x9 VIEWING surface (interpolated; NOT the '
                        f'grading source) - actuator {aid} '
                        f'{entry["master_name"]}; tau about {p_dof} over '
                        f'{p_dof} x {s_dof}'])
            w.writerow([f'# values bilinearly interpolated from the '
                        f'{BOX_DEG:g}-deg OpenSim box samples; grading '
                        f'values live in human_surfaces_5x5_grading_'
                        f'20261023.csv (direct evaluations)'])
            w.writerow([f'# interpolation smooths via-point discontinuity '
                        f'corners for viewing only; curves never corrected'])
            w.writerow(['p_node_deg', 's_node_deg', 'tau_h_Nm_interp',
                        'in_graded_domain'])
            for tp in vp:
                for ts in vs:
                    ind = 1
                    if envelope:
                        fl = tp if p_dof == 'hip_flexion_r' else ts
                        ad = ts if p_dof == 'hip_flexion_r' else tp
                        ind = 1 if (env_abd_min(fl) <= ad
                                    <= env_add_max(fl)) else 0
                    w.writerow([f'{tp:.6g}', f'{ts:.6g}',
                                f'{bilinear(tau_box, box_p, box_s, tp, ts):.9f}',
                                ind])
    print(f'wrote {len(man_acts)} files in {OUT_DIR9}')

    total_graded = sum(e['n_graded_cells'] for e in man_acts)
    total_flagged = sum(e['n_flagged_graded_cells'] for e in man_acts)
    manifest = {
        'generated': datetime.datetime.now().isoformat(),
        'purpose': ('2-D SINGLE-MUSCLE human capability surfaces (v2 '
                    'metric: coupled-DOF surfaces) for the 27-actuator '
                    'capability table; single grading source of truth for '
                    'BOTH rows (hand-tuned and optimized)'),
        'spec': ['FIX_SPEC_20261023.md sections 1, 1b, 2',
                 'opensim_limits_checklist_20261023.md (flag layer + '
                 'gates)'],
        'supersedes': ('nothing on disk - v1 1-D CSVs '
                       '(human_torque_curves_20261022*.csv) stay '
                       'untouched; this is the v2 surface regeneration'),
        'opensim_version': ver,
        'python': sys.executable,
        'model': MODEL, 'model_sha1': msha,
        'map': MAP, 'map_sha1': psha,
        'audit_csv': AUDIT, 'audit_csv_sha1': sha1_of(AUDIT),
        'probe_csv_crosscheck': os.path.join(
            OUTDIR, 'probe_osim_limits_2d_20261023.csv'),
        'target_rule': {
            'ruling': ('same-named single stock muscle per actuator; no '
                       'cross-muscle group sums; one anatomical muscle '
                       'split into compartments keeps its own compartments '
                       '(glut_max1-3, add_mag1-3); explicit rulings 5/11/'
                       '15/16/24/25 (Ben 2026-10-22)'),
            'per_actuator': {str(k): TARGET_RULE[k][1] for k in TARGET_RULE},
        },
        'surface_definition': {
            'tau': ('tau_h(theta_p, theta_s) = sum over CHOSEN muscles of '
                    'computeMomentArm(about the PRIMARY dof) x Fmax; '
                    'OpenSim native signs (arm = -dL/dq; tau>0 = positive '
                    'coordinate direction; knee flexion NEGATIVE: vasti '
                    'positive / hamstrings negative)'),
            'grid': ('5x5 INCLUSIVE uniform nodes per axis over '
                     'rom_table_deg (theta_i = lo + (hi-lo)*i/4); 25 '
                     'cells per actuator; grading values are DIRECT '
                     'OpenSim evaluations at the nodes (never '
                     'interpolated); cell bin = midpoint-Voronoi interval '
                     'within [lo, hi] (outer bins half-width)'),
            'area_share': ('headline metric is AREA SHARE over the graded '
                           'cells (see FIX_SPEC 1 phrasing style); '
                           'envelope-excluded cells are removed from the '
                           'denominator for BOTH rows and never counted '
                           'as failures'),
            'standing_trap': ('human surfaces computed on STOCK '
                              'gait2392_simbody routing only; robot/BPA '
                              'paths never substituted'),
        },
        'hip_adduction_envelope': {
            'ruling': 'FIX_SPEC 1b (Ben 2026-10-23)',
            'applies_to_actuators': [24, 25, 26, 27],
            'rule': ('adduction max = 20 deg at flexion -15 tapering '
                     'LINEARLY to +10 deg at flexion +110; abduction max '
                     'from -10 deg at full extension to -60 deg at full '
                     'flexion; LINEAR interpolation between extremes '
                     '(supersedes the 2026-10-22 quarter-ellipse '
                     'placeholder)'),
            'physical_origin': ('BONE INTERFERENCE (femur-pelvis contact, '
                                'per Ben) - a real anatomical/hardware '
                                'constraint, not a modeling convenience'),
            'cell_rule': ('cell graded iff its CENTER (p_node, s_node) is '
                          'inside the closed envelope'),
            'per_actuator_excluded_cells': {
                str(e['id']): e['envelope']['n_excluded_cells']
                for e in man_acts if e['envelope']['applies']},
            'hip_rotation': ('stays the static rom_table_deg box until '
                             'Ben rules otherwise; pinned at 0 on '
                             'actuator 24 (not graded anywhere)'),
            'tfl_note': ('actuator 7 grades knee x hip_flexion; its '
                         'hip_adduction (2nd map secondary) is PINNED AT '
                         '0, which is inside the envelope at every graded '
                         'flexion, so no envelope exclusion applies to '
                         'its surface'),
        },
        'pinning_3plus_joint_routes': {
            '1': 'lumbar_rotation pinned at 0 (grades lumbar_extension x '
                 'lumbar_bending)',
            '7': 'hip_adduction_r pinned at 0',
            '20': 'subtalar_angle_r pinned at 0', '21': 'subtalar_angle_r '
                 'pinned at 0', '22': 'subtalar_angle_r pinned at 0',
            '23': 'subtalar_angle_r pinned at 0',
            '24': 'hip_rotation_r pinned at 0',
        },
        'uniarticular_class_serving': {
            'note': ('FIX_SPEC 1: the BiPam class serves uniarticular '
                     'muscles over its two hinges; the five actuators '
                     'whose map entry lists no secondary grade on their '
                     'CLASS hinge pair'),
            '10': 'knee_angle_r x hip_flexion_r (Bi_Hip thetaHipKnee)',
            '11': 'knee_angle_r x hip_flexion_r (Bi_Hip thetaHipKnee)',
            '12': 'ankle_angle_r x knee_angle_r (Calves thetaKneeAnkle; '
                  'human soleus is knee-independent so its surface is '
                  'flat in knee)',
            '26': 'hip_flexion_r x hip_adduction_r (Uni_Hip two-hinge '
                  'pair; envelope-coupled)',
            '27': 'hip_flexion_r x hip_adduction_r (Uni_Hip; envelope-'
                  'coupled)',
        },
        'flag_layer': {
            'policy': ('flagged cells STILL COUNT in the headline meet '
                       'percentage (conservative); with-flag and '
                       'without-flag numbers are BOTH reported; curves '
                       'are NEVER corrected'),
            'FLAG_VIA': ('ConditionalPathPoint isActive-set change inside '
                         'the cell bin; 0.5-deg cross sweeps (primary '
                         'swept at each secondary node and vice versa); '
                         'boundary is a constant-angle line '
                         '(checklist-verified)'),
            'FLAG_STEP': ('adjacent-sample |d(arm)| > 2.0% of the '
                          'per-muscle SURFACE-WIDE peak A* (never '
                          'per-column); 0.5-deg cross sweeps'),
            'FLAG_SLACK': ('any 1-deg full-box sample in the cell bin '
                           'with l_mt <= tendon_slack_length'),
            'FLAG_SHORT': ('any 1-deg full-box sample in the cell bin '
                           'with lceN_est < 0.5 (constant-width pennation '
                           'estimate, lce = sqrt((l_mt-tsl)^2 + (lo sin '
                           'phi0)^2); upper bound on true fiber length '
                           'when tendon is loaded)'),
            'FLAG_BONE': ('NOT APPLIED - bone-penetration check TODO and '
                          'not run (checklist L5); disclosed here, no '
                          'flags invented'),
            'audit_caveats': ('the 1-D audit of record '
                              '(path_artifact_audit_...034030.csv) '
                              'primary-DOF caveat strings are carried '
                              'per actuator in the grading CSV column '
                              'audit_caveat_1d_primary'),
        },
        'halt_gates': {
            'version_pin': f'{ver} == {PINNED_VERSION}: PASS',
            'sha1_pins': 'model+map+audit pins: PASS',
            'live_wrap_sets': 'empty on every chosen muscle: PASS',
            'gate_table': 'checklist section 7 hard gates: PASS',
            'fd': ('non-knee primary moment-arm FD (5% tol): PASS; knee '
                   'informational per documented wrap/moving-point '
                   'nonsmoothness'),
        },
        'viewing_export': {
            'dir': OUT_DIR9,
            'grid': '9x9 inclusive uniform nodes per axis',
            'method': ('bilinear interpolation from the 1-deg OpenSim box '
                       'samples; grading values in the 5x5 CSV are direct '
                       'evaluations; interpolation smooths via-point '
                       'discontinuities for VIEWING only'),
        },
        'reporting': {
            'headline': ('robot-row meet percentages computed against '
                         'tau_h_Nm in the 5x5 CSV count FLAGGED cells '
                         '(conservative)'),
            'without_flag_variant': ('meet percentage over graded cells '
                                     'with no flags (denominator = '
                                     'n_graded_cells - '
                                     'n_flagged_graded_cells per '
                                     'actuator)'),
            'totals': {
                'n_actuators': len(man_acts),
                'n_cells_total': len(csv_rows),
                'n_graded_cells': total_graded,
                'n_excluded_envelope_cells': len(csv_rows) - total_graded,
                'n_flagged_graded_cells': total_flagged,
            },
        },
        'actuators': man_acts,
    }
    with open(OUT_MAN, 'w') as f:
        json.dump(manifest, f, indent=1)
    print(f'wrote {OUT_MAN}')

    print('\nflagged actuators (>=1 flagged graded cell):')
    for n in flagged_acts:
        print(f'  - {n}')
    print(f'\nDONE ({len(man_acts)} actuators; {total_graded} graded '
          f'cells; {total_flagged} flagged graded cells)')


if __name__ == '__main__':
    main()
