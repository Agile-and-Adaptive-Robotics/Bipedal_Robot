from pathlib import Path
from scipy.io import loadmat
import json
W=Path(__file__).resolve().parent
R=next(p for p in W.parents if p.name=='Bipedal_Robot')
out={}
for fn in ['minimizeFlxPin10_results.mat','minimizeExt10mmX3_results.mat']:
    m=loadmat(R/'Testing_Data/2022_02_Festo'/fn,simplify_cells=True)
    info={}
    for k,v in m.items():
        if k.startswith('__'):continue
        desc={'type':type(v).__name__,'shape':list(v.shape) if hasattr(v,'shape') else None}
        if isinstance(v,dict):desc['fields']=list(v)
        if k in ['POP','MAXGEN','MODE','ALLBPA','NUMHOLD','SOLVER','PICK','LOCKSRC','TRANSMODE','USE_BRACKET2']:
            desc['value']=v.tolist() if hasattr(v,'tolist') else str(v)
        if k=='filtered_results':desc['adopted_row']=(v[76] if fn.startswith('minimizeFlx') else v[31]).tolist()
        info[k]=desc
    out[fn]=info
(W/'build/identification_inventory.json').write_text(json.dumps(out,indent=2,default=str),encoding='utf-8')
print(json.dumps(out,indent=2,default=str))
