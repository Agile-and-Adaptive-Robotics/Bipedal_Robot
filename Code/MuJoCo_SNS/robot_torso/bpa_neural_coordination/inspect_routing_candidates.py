"""Inventory historical Routing27 dimensions; do not assign them to CAD routes."""
import json,re
from pathlib import Path
from scipy.io import loadmat
import numpy as np
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[2]/'Matlab/Mesh_Optimization/Routing27_v2/Results'
files={}
for p in SOURCE.glob('Routing27_a*.mat'):
    key=int(re.search(r'_a(\d+)_',p.name).group(1))
    if key not in files or p.name>files[key].name:files[key]=p
records=[]
for key,p in sorted(files.items()):
    data=loadmat(p,variable_names=['best','ctx'],simplify_cells=True)
    best=data.get('best',{});ctx=data.get('ctx',{});res=best.get('res',{})
    if not isinstance(res,dict) or 'x' not in res:
        records.append(dict(actuator_id=key,source=str(p),status='No readable selected x',best_keys=list(best)))
        continue
    x=np.asarray(res['x']).ravel()
    records.append(dict(actuator_id=key,source=str(p),configuration=best.get('cfgName'),
        rest_m=float(x[int(ctx['iRest'])-1]),tendon_m=float(x[int(ctx['iTend'])-1]),
        diameter_mm=float(ctx['diameter']),pressure_kpa=float(ctx['pressure']),
        fitting_parameter_m=float(ctx['fitting']),worst_constraint=float(res.get('worstC',np.nan)),
        margin_relative=float(res.get('marginRel',np.nan)),
        warning='Historical robotbody primary-DOF optimization, not verified against retargeted AH001 CAD paths; pulley, fitting, compliance, and ROM conventions require reconciliation.'))
out=HERE/'runs/cad_parameter_audit_20261031';out.mkdir(parents=True,exist_ok=True)
(out/'HISTORICAL_ROUTING_CANDIDATES.json').write_text(json.dumps(dict(source=str(SOURCE),candidates=records,
    absent_master_ids=[i for i in range(1,28) if i not in files],candidate_count=len(records)),indent=2))
print(json.dumps(dict(candidate_count=len(records),records=records),indent=2))
