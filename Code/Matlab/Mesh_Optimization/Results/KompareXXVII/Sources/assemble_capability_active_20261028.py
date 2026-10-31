"""Assemble Ben's two-metric capability table from canonical frozen-design replay.

Run regrade_capability_active_20261028.m first. No mechanics are implemented here.
Historical mats and grading CSVs are read-only; corrected outputs have new names.
"""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'tmp/torque_comparison_20261028'
CAP = ROOT / 'tmp/paper_bridge/capability_table_20261022'
records = list(csv.DictReader((OUT / 'canonical_replay.csv').open(newline='')))
roster = list(csv.DictReader((CAP / 'BiPulley_survey2_perMuscle_20261023.csv').open(newline='')))
assert len(records) == 54 and len(roster) == 27
lookup = {(r['site'], r['side']): r for r in records}
assert len(lookup) == 54
rows = []
for item in roster:
    o, u = [lookup[item['site'], s] for s in ('original', 'updated')]
    assert o['domain'] == u['domain']
    assert abs(float(u['meet']) - float(item['meetPct'])) < 0.0001
    for r in (o,u):
        assert 0 <= float(r['meet']) <= float(r['active']) <= 100
        assert float(r['max_torque_delta']) < 1e-9
    rows.append(dict(actuator_id=int(item['actuator_id']), muscle=item['master_name'],
        site=item['site'], original_meet_pct=float(o['meet']),
        original_active_pct=float(o['active']), updated_meet_pct=float(u['meet']),
        updated_active_pct=float(u['active']), domain_cells=int(o['domain']),
        original_mat=o['mat'], updated_mat=u['mat'],
        original_nBPA=int(o['nBPA']), original_G=int(o['G']), original_diameter_mm=float(o['diameter']),
        updated_nBPA=int(u['nBPA']), updated_G=int(u['G']), updated_diameter_mm=float(u['diameter'])))
with (OUT / 'goal_table_torque_active_corrected_20261028.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
means={k:sum(r[k] for r in rows)/27 for k in ('original_meet_pct','original_active_pct','updated_meet_pct','updated_active_pct')}
assert abs(means['original_meet_pct']-37.28395061728395)<1e-10
assert abs(means['updated_meet_pct']-65.18518518518519)<1e-10
assert sum(r['updated_meet_pct']>r['original_meet_pct']+1e-9 for r in rows)==23
summary=r'''\begin{table}[htbp]
\centering
\caption{Mean torque and active-force coverage across 27 actuator sites at \qty{620}{\kPa}.}
\label{tab:CapabilitySurvey}
\begin{tabular}{@{}lcc@{}}
\toprule
\textbf{Configuration} & \textbf{Human torque met (\%)} & \textbf{Active force (\%)} \\
\midrule
Baseline & 37.28 & 92.44 \\
Transmission configurations & 65.19 & 98.81 \\
\bottomrule
\end{tabular}
\end{table}
'''
(OUT / 'capability_summary_table.tex').write_text(summary,encoding='utf-8')
prose=r'''Active-force coverage is the fraction of included grid cells with positive insertion tension, excluding slack or infeasible configurations. It does not require the torque to match the human reference. Both coverage measures are calculated for each actuator and then averaged equally across the 27 sites. Table~\ref{tab:CapabilitySurvey} reports the resulting means.

The mean fraction meeting the human torque requirement increases from \qty{37.28}{\percent} for the baseline to \qty{65.19}{\percent} for the saved transmission configurations, improving at 23 sites and tying at four. Mean active-force coverage increases from \qty{92.44}{\percent} to \qty{98.81}{\percent}. The higher active-force coverage shows that the transmission configurations maintain tension over more of the graded domain, while the lower torque coverage shows that positive tension alone is insufficient to meet the human torque requirement.
'''
(OUT / 'capability_results_prose.tex').write_text(prose,encoding='utf-8')
detail=[r'''\begin{table}[p]
\centering\small
\caption{Torque and active-force coverage at each actuator site. Percentages use the included grid cells for that site.}
\label{tab:capability_by_site}
\begin{tabular}{@{}lrrrr@{}}
\toprule
& \multicolumn{2}{c}{Baseline (\%)} & \multicolumn{2}{c}{Transmission (\%)} \\
\cmidrule(lr){2-3}\cmidrule(l){4-5}
\textbf{Site} & \textbf{Torque met} & \textbf{Active force} & \textbf{Torque met} & \textbf{Active force} \\
\midrule
''']
for r in rows:
    vals=[r[k] for k in ('original_meet_pct','original_active_pct','updated_meet_pct','updated_active_pct')]
    detail.append(r['site'].replace('_',r'\_')+' & '+' & '.join(f'{v:.2f}' for v in vals)+r' \\'+'\n')
detail.append(r'''\midrule
Mean & 37.28 & 92.44 & 65.19 & 98.81 \\
\bottomrule
\end{tabular}
\end{table}
''')
(OUT/'capability_by_site_table.tex').write_text(''.join(detail),encoding='utf-8')
files={Path(r['mat']) for r in records}
files.update([ROOT/'Code/Matlab/Mesh_Optimization/BiPam_pulley2_20261022.m',
    CAP/'BiPulley_survey2_perMuscle_20261023.csv',CAP/'human_surfaces_5x5_grading_20261023.csv',
    CAP/'human_surfaces_manifest_20261023.json',Path(__file__),
    Path(__file__).with_name('regrade_capability_active_20261028.m'),
    Path(__file__).with_name('biPulleyActiveForceMask_20261028.m')])
manifest={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
result=dict(means=means,configurations=54,sites=27,max_torque_delta_Nm=max(float(r['max_torque_delta']) for r in records),
    changed_active_rows={s:[dict(site=r['site'],before=float(r['old_active']),after=float(r['active'])) for r in records if r['side']==s and abs(float(r['old_active'])-float(r['active']))>1e-9] for s in ('original','updated')},input_sha256=manifest)
(OUT/'verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='input_sha256'},indent=2))

