"""Assemble the 27-actuator capability-table deliverables (2026-10-22,
second assembly pass, easteregg2). Re-derives every number from the three
CSVs of record in this folder and fails hard on any mismatch:

  inputs (READ-ONLY):
    hand_tuned_metrics_20261022.csv                      (hand-tuned row)
    BiPulley_survey_perMuscle_regraded_20261022_223448.csv (optimized row)
    human_torque_curves_20261022.csv                     (targets + flags)

  outputs (overwritten in place; names are the dated deliverable names):
    per_muscle_capability_20261022.csv
    capability_table_20261022.md
    capability_table_20261022.tex
    capability_table_20261022.xlsx
    CAVEATS_REPORT_20261022.md   (section 11 rewritten by hand after the run)
    assemble_summary_20261022.json

AVERAGING RULES (Ben's ask): hand-tuned row averages the covered muscles
(27 of 27); optimized row averages the muscles that were built AND could be
graded (16; 23 of 27 built, 4 not buildable, 7 built but grade-unable).
Row labels state the count actually averaged.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HAND = os.path.join(HERE, 'hand_tuned_metrics_20261022.csv')
REGRADE = os.path.join(
    HERE, 'BiPulley_survey_perMuscle_regraded_20261022_223448.csv')
CURVES = os.path.join(HERE, 'human_torque_curves_20261022.csv')


def read_commented_csv(path):
    with open(path, newline='') as f:
        rows = list(csv.reader(ln for ln in f if not ln.startswith('#')))
    hdr = rows[0]
    return [dict(zip(hdr, r)) for r in rows[1:]]


def as_float(s):
    s = (s or '').strip()
    return float(s) if s not in ('', '0') and s != '0' else (
        float(s) if s else None)


hand = read_commented_csv(HAND)
reg = read_commented_csv(REGRADE)
cur = read_commented_csv(CURVES)

# ---------------- hard cross-checks ----------------
assert len(hand) == 27, f'hand rows {len(hand)} != 27'
assert len(reg) == 27, f'regrade rows {len(reg)} != 27'
assert len(cur) == 27, f'curve rows {len(cur)} != 27'
for i, (h, r, c) in enumerate(zip(hand, reg, cur)):
    assert h['muscle'] == r['muscle'] == c['master_name'], \
        f'row {i + 1}: {h["muscle"]} | {r["muscle"]} | {c["master_name"]}'
    assert int(c['actuator_id']) == i + 1
    assert h['source'] == 'connor-corrected', h['source']

built = [r for r in reg if r['built'] == 'true']
graded = [r for r in built if int(r['nBPA']) > 0
          and r['converged'] == '1']
unbuilt = [r for r in reg if r['built'] != 'true']
built_ungradeable = [r for r in built if r not in graded]
assert len(built) == 23, f'built {len(built)} != 23'
assert len(graded) == 16, f'graded {len(graded)} != 16'
assert len(unbuilt) == 4, f'unbuilt {len(unbuilt)} != 4'
assert len(built_ungradeable) == 7, \
    f'built-ungradeable {len(built_ungradeable)} != 7'
assert len([r for r in graded if r['converged'] == '1']) == 16

suspects = [h for h in hand if h['suspect'] == 'yes']
assert len(suspects) == 9, f'suspects {len(suspects)} != 9'
flagged = [c for c in cur if c['artifact_flag_caveat'] != 'none']
assert len(flagged) == 6, f'flagged curves {len(flagged)} != 6'
assert [c['master_name'] for c in flagged] == [
    'Semimembranosus', 'Rectus Femoris', 'Gracilis', 'Vastus Intermedius',
    'Iliacus', 'Psoas X']

hand_meet = sum(float(h['meetPct']) for h in hand) / 27.0
hand_act = sum(float(h['activePct']) for h in hand) / 27.0
opt_meet = sum(float(r['meetPct']) for r in graded) / 16.0
opt_act = sum(float(r['activePct']) for r in graded) / 16.0
assert abs(hand_meet - 72.8888888889) < 1e-6, hand_meet
assert abs(hand_act - 98.8148148148) < 1e-6, hand_act
assert abs(opt_meet - 52.0) < 1e-9, opt_meet
assert abs(opt_act - 82.25) < 1e-9, opt_act

evals = sorted(int(r['evalsUsed']) for r in graded)
med = (evals[7] + evals[8]) / 2
capped = [r for r in graded if 'capability-capped' in r['note']]
print(f'checks: built={len(built)} graded={len(graded)} unbuilt={len(unbuilt)}'
      f' built_ungradeable={len(built_ungradeable)} suspects={len(suspects)}'
      f' flagged={len(flagged)}')
print(f'avgs: hand {hand_meet:.2f}/{hand_act:.2f}  opt {opt_meet:.2f}/{opt_act:.2f}')
print(f'graded evals: min {evals[0]} median {med:g} max {evals[-1]}; '
      f'capability-capped rows {len(capped)}')

# ---------------- per-muscle CSV ----------------
caveat_by_id = {int(c['actuator_id']): c['artifact_flag_caveat'] for c in cur}
dof_by_id = {int(c['actuator_id']): c['primary_dof'] for c in cur}
rows_out = []
for i, (h, r) in enumerate(zip(hand, reg)):
    aid = i + 1
    is_graded = r in graded
    opt_note = r['note']
    parts = ['HAND: ' + h['note'], 'OPT: ' + opt_note]
    if caveat_by_id[aid] != 'none':
        parts.append('CURVE-ARTIFACT: ' + caveat_by_id[aid])
    rows_out.append({
        'muscle': h['muscle'],
        'actuator_id': aid,
        'primary_dof': dof_by_id[aid],
        'coverage_source': h['source'],
        'hand_target_version': 'single-muscle curves',
        'opt_target_version': ('group-targets chosen; regraded vs '
                               'single-muscle curves' if is_graded
                               else 'n/a (not graded)'),
        'hand_meet_pct': f"{float(h['meetPct']):.4f}",
        'hand_active_pct': f"{float(h['activePct']):.4f}",
        'opt_meet_pct': f"{float(r['meetPct']):.4f}" if is_graded else '',
        'opt_active_pct': f"{float(r['activePct']):.4f}" if is_graded else '',
        'opt_nBPA': r['nBPA'],
        'opt_G': r['tackleG'],
        'opt_evals_used': r['evalsUsed'],
        'opt_converged': 'yes' if r['converged'] == '1' else 'no',
        'path_artifact_flag': 'yes' if caveat_by_id[aid] != 'none' else 'no',
        'transcription_suspect_flag': h['suspect'],
        'notes': ' | '.join(parts),
    })

per_csv = os.path.join(HERE, 'per_muscle_capability_20261022.csv')
with open(per_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
    w.writeheader()
    w.writerows(rows_out)
print(f'wrote {per_csv}')

# ---------------- markdown table ----------------
md = f"""# 27-actuator capability table (2026-10-22)

Each cell is the row average: the share of the 25-point range-of-motion grid
at the actuator's primary coordinate where the condition holds (every other
joint the muscle crosses pinned at 0 deg). Human targets: single same-named
muscle torque curves, Gait2392 via OpenSim 4.4.1
(`human_torque_curves_20261022.csv`).

| Row | % of RoM torque meets or exceeds human | % of RoM with active force |
|---|---|---|
| Hand-tuned (27 of 27 graded) | {hand_meet:.2f}% | {hand_act:.2f}% |
| Optimized (16 of 27 graded; 23 built) | {opt_meet:.2f}% | {opt_act:.2f}% |

- **Hand-tuned row:** averages all 27 covered actuators (0 excluded).
- **Optimized row:** Ben asked for all 27; **23 of 27 were built**, and the
  average runs over the **16 built muscles that could be graded**. The row
  label states the count actually averaged; the missing are stated, never
  faked:
  - Could not build (4): Erector Spinae, Internal Oblique, External Oblique,
    Vastus Intermedius - the master muscle path crosses no joint (all points
    on one body), so no joint actuator can be formed.
  - Built but grade-unable (7): Rectus Femoris, Tibialis Posterior, Peroneus
    Brevis, Peroneus Longus, Peroneus Tertius, Adductor Magnus, Psoas X -
    built from the seed route, but the regrade spec's hinge chain does not
    model their primary coordinate (knee_angle_r for Rectus Femoris;
    subtalar_angle_r for the three peroneals and Tibialis Posterior;
    hip_adduction_r for Adductor Magnus; hip_flexion_r for Psoas X), so no
    regraded percentage exists. Their opt grade cells are blank in the
    per-muscle CSV.
  - Full reasons, the nine transcription suspects, the six artifact-flagged
    human curves, and the budget disclosure: `CAVEATS_REPORT_20261022.md`.
"""
md_path = os.path.join(HERE, 'capability_table_20261022.md')
with open(md_path, 'w') as f:
    f.write(md)
print(f'wrote {md_path}')

# ---------------- LaTeX table ----------------
tex = r"""% capability_table_20261022.tex -- 27-actuator capability table,
% assembled 2026-10-22 from per_muscle_capability_20261022.csv in this
% folder. Float placement NOT decided: Ben wraps/places it.
% Suggested caption: Isometric capability of the hand-tuned and optimized
% actuator routings against the same-named human muscle torque curves
% (Gait2392, OpenSim 4.4.1). Each cell is the row average over the 25-point
% range-of-motion grid at the actuator's primary coordinate, with every
% non-graded joint the muscle crosses pinned at 0 deg. The hand-tuned row
% averages all 27 covered actuators; the optimized row averages the 16 of
% 23 built actuators that could be graded on the survey hinge chain
% (4 could not be built, 7 built but grade-unable; see
% CAVEATS_REPORT_20261022.md).
\begin{tabular}{lcc}
\toprule
 & Torque meets or exceeds human (\% of RoM) & Active force (\% of RoM) \\
\midrule
Hand-tuned (27 of 27 graded) & {HM}\% & {HA}\% \\
Optimized (16 of 27 graded; 23 built) & {OM}\% & {OA}\% \\
\bottomrule
\end{tabular}
""".replace('{HM}', f'{hand_meet:.2f}').replace('{HA}', f'{hand_act:.2f}') \
   .replace('{OM}', f'{opt_meet:.2f}').replace('{OA}', f'{opt_act:.2f}')
tex_path = os.path.join(HERE, 'capability_table_20261022.tex')
with open(tex_path, 'w') as f:
    f.write(tex)
print(f'wrote {tex_path}')

# ---------------- xlsx ----------------
import openpyxl  # noqa: E402  (verified importable in this env before run)

wb = openpyxl.Workbook()
ws = wb.active
ws.title = 'Capability table'
ws.append(['27-actuator capability table, 2026-10-22'])
ws.append(['Cells are row averages over the 25-point RoM grid at each '
           'actuator\'s primary coordinate; every non-graded joint pinned '
           'at 0 deg.'])
ws.append([])
ws.append(['Row', '% of RoM torque meets or exceeds human',
           '% of RoM with active force', 'Muscles averaged'])
ws.append(['Hand-tuned (27 of 27 graded)', round(hand_meet, 2),
           round(hand_act, 2), 27])
ws.append(['Optimized (16 of 27 graded; 23 built)', round(opt_meet, 2),
           round(opt_act, 2), 16])
ws.append([])
ws.append(['Traceability'])
ws.append(['Hand-tuned meet sum / count',
           sum(float(h['meetPct']) for h in hand), 27])
ws.append(['Hand-tuned active sum / count',
           sum(float(h['activePct']) for h in hand), 27])
ws.append(['Optimized meet sum / count',
           sum(float(r['meetPct']) for r in graded), 16])
ws.append(['Optimized active sum / count',
           sum(float(r['activePct']) for r in graded), 16])
ws.append(['Built of 27', 23,
           '4 could not build; 7 built but grade-unable (see per-muscle '
           'sheet + caveats report)'])
ws.append(['Hand-tuned source', 'hand_tuned_metrics_20261022.csv'])
ws.append(['Optimized source',
           'BiPulley_survey_perMuscle_regraded_20261022_223448.csv'])
ws.append(['Human targets', 'human_torque_curves_20261022.csv '
           '(single-muscle, OpenSim 4.4.1)'])

ws2 = wb.create_sheet('Per-muscle rows')
ws2.append(list(rows_out[0].keys()))
for ro in rows_out:
    ws2.append([ro[k] for k in ro])
xlsx_path = os.path.join(HERE, 'capability_table_20261022.xlsx')
wb.save(xlsx_path)
print(f'wrote {xlsx_path}')

# ---------------- summary json ----------------
summary = {
    'assembled': '2026-10-22 (second pass, easteregg2; re-derived from the '
                 'CSVs of record by assemble_deliverables_20261022.py)',
    'hand_meet_avg': round(hand_meet, 2),
    'hand_act_avg': round(hand_act, 2),
    'opt_meet_avg': round(opt_meet, 2),
    'opt_act_avg': round(opt_act, 2),
    'hand_covered': 27,
    'opt_covered_graded': 16,
    'opt_built': 23,
    'could_not_build': 4,
    'built_grade_unable': 7,
    'transcription_suspects': 9,
    'artifact_flagged_curves': 6,
    'unconverged_surveyed_sites': 0,
    'files': [per_csv, md_path, tex_path, xlsx_path],
}
with open(os.path.join(HERE, 'assemble_summary_20261022.json'), 'w') as f:
    json.dump(summary, f, indent=2)
print('wrote assemble_summary_20261022.json')
print('DONE')
