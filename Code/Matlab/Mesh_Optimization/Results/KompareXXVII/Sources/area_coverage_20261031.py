"""Postprocess frozen grid indicators with trapezoidal angular-domain weights.
No mechanical model, fit, simulation, or optimizer is run. Normalize each site
before averaging sites. Sparse-grid quadrature estimates are not exact boundaries.
"""
from pathlib import Path
import csv,hashlib,json,warnings
import numpy as np
from scipy.io import loadmat
ROOT=next(p for p in Path(__file__).resolve().parents if p.name=='Bipedal_Robot')
HERE=Path(__file__).resolve().parent
rows=list(csv.DictReader((HERE/'canonical_replay.csv').open(newline='')))
out=[]
for record in rows:
    path=ROOT/'Code/Matlab/Mesh_Optimization/Results/BiPulley_surveys'/Path(record['mat'].replace('\\','/')).name
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        saved=loadmat(path,simplify_cells=True,variable_names=['dom','g','spec'])
    dom=np.asarray(saved['dom'],bool);g=saved['g'];meet=np.asarray(g['meet'],bool)
    active=np.asarray(g['taut'],bool)&~np.asarray(g['slack'],bool)&~np.asarray(g['infeas'],bool)&np.isfinite(g['tauB'])
    assert np.sum(meet&dom)==int(record['meet_cells'])
    assert np.sum(active&dom)==int(record['active_cells'])
    assert not np.any(meet&dom&~active)
    ndim=int(np.sum(np.ptp(np.asarray(saved['spec']['kin']['thetaRangesDeg']),axis=1)>0))
    if ndim==1:
        for mask in [meet,active]:assert np.all(mask==mask[0,:]) or np.all(mask==mask[:,0,None])
    weights=[]
    for size in dom.shape:
        axis=np.ones(size)
        if size>1:axis[[0,-1]]=.5
        weights.append(axis)
    weight=np.outer(*weights);den=np.sum(weight*dom)
    m=100*np.sum(weight*dom*meet)/den;a=100*np.sum(weight*dom*active)/den
    assert 0<=m<=a<=100
    out.append(dict(site=record['site'],side=record['side'],domain_dimensions=ndim,node_meet_pct=float(record['meet']),area_meet_pct=float(m),area_active_pct=float(a),input_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
means={side:{key:float(np.mean([x[key] for x in out if x['side']==side])) for key in ['area_meet_pct','area_active_pct']} for side in ['original','updated']}
result=dict(means=means,method='Tensor trapezoidal quadrature of sampled domain and satisfaction indicators on original uniform grids; normalize within each site, then equal site mean. One-DoF repeated axis verified and cancels. Coarse-grid estimates, no physical model interpolation or reevaluation.',rows=out)
(HERE/'area_coverage_20261031.json').write_text(json.dumps(result,indent=2),encoding='utf8')
with (HERE/'area_coverage_20261031.csv').open('w',newline='',encoding='utf8') as f:
    writer=csv.DictWriter(f,fieldnames=out[0].keys());writer.writeheader();writer.writerows(out)
print(json.dumps(means,indent=2))
