"""Diff the 2026-10-22 audit CSV of record against the 2026-10-23 re-run."""
import csv

OLD = (r'D:\GitHub\Bipedal_Robot\tmp\paper_bridge\opt_improve_20261022'
       r'\path_artifact_audit_gait2392_simbody_20261022_034030.csv')
NEW = (r'D:\GitHub\Bipedal_Robot\tmp\paper_bridge\capability_table_20261022'
       r'\audit_rerun_20261023\path_artifact_audit_gait2392_simbody_'
       r'20261023_021610.csv')


def load(p):
    with open(p, newline='') as f:
        return list(csv.DictReader(
            (ln for ln in f if not ln.startswith('#'))))


old, new = load(OLD), load(NEW)
print('rows old/new:', len(old), len(new))
diffs = 0
for o, n in zip(old, new):
    key = (o['actuator_id'], o['dof'], o['kind'], o['muscle'])
    assert key == (n['actuator_id'], n['dof'], n['kind'], n['muscle']), key
    for c in ('cond_in_range', 'n_activation_angles', 'activation_angles_deg',
              'n_step_angles', 'step_angles_deg', 'max_step_pct',
              'artifact_flag'):
        if o[c] != n[c]:
            diffs += 1
            print('DIFF', key, c, 'old=', repr(o[c]), 'new=', repr(n[c]))
print('total field diffs:', diffs)
print('old flag=1 rows:', sum(1 for r in old if r['artifact_flag'] == '1'))
print('new flag=1 rows:', sum(1 for r in new if r['artifact_flag'] == '1'))
