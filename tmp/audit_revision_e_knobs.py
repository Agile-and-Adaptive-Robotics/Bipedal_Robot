from pathlib import Path
import ast,json
R=Path(r'D:\Github\Bipedal_Robot'); S=R/'Code/MuJoCo_SNS/spinal'
t=ast.parse((S/'w2l_mujoco/test_w2l_ground.py').read_text(encoding='utf-8'))
node=next(n for n in t.body if isinstance(n,ast.Assign) and any(isinstance(a,ast.Name) and a.id=='KV' for a in n.targets))
known={k.arg for k in node.value.keywords}
p=json.loads((S/'reports_20261025_persyn/round_w25.json').read_text())['params']
r={'checked_commit':'a988508a','method':'AST-only inspection; no simulation import or execution','winner_key_count':len(p),'accepted_by_evaluator':sorted(set(p)&known),'discarded_by_evaluator':sorted(set(p)-known),'reason':'evaluate() retains only keys already in T.KV; runner gains dict also does not forward the individual keys','scope':'Committed local source. A newer source on easteregg2 could differ; recover its hash before assigning these outputs to this code.'}
o=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/review/parameter_plumbing_audit.json';o.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r,indent=2))
