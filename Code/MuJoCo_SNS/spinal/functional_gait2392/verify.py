"""Topology, anatomical recruitment, live-cache, and backend checks."""
import os,sys,json,copy
from pathlib import Path
HERE=Path(__file__).resolve().parent;S=HERE.parent;sys.path[:0]=[str(S),str(HERE),str(S/'w2l_cpg')]
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import mujoco
import bsolve_ik as B
import network as F
from preallocated_exact_eval import PreallocatedExactForward
from mechanics import JOINTS
OUT=S/'reports_20261028_functional_pf'
m,d=B.build_air_model(False)
acts=[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(m.nu)]
net=F.build(acts,options=dict(axes=True,biarticular=True))
reference=copy.deepcopy(net.compiled);net.compiled=PreallocatedExactForward(net.compiled);net.bind_mujoco(m,d)
# Compare the exact same dynamically configured matrix under the reference backend.
reference.g_max_non[:]=net.compiled.g_max_non
rng=np.random.default_rng(20261028);maximum=0.
for _ in range(200):
    u=net.make_inputs();u[:]=rng.uniform(0,.2,len(u));u[net.inputs.index('DRIVE')]=2.
    reference.forward(u);net.compiled.forward(u)
    maximum=max(maximum,float(np.max(abs(reference.V-net.compiled.V))))
assert maximum==0.
base=np.load(S/'bsolve_out.npz',allow_pickle=True)['qpos'][0]
ai=acts.index('tfl_r');out=[]
for flex in (-20.,100.):
    d.qpos[:]=base
    for side in ('r','l'):
        for coord,val in [('hip_flexion',flex),('hip_adduction',0.),('hip_rotation',0.),('knee_angle',-20.)]:
            d.qpos[net.arm_model.addresses[JOINTS.index(coord+'_'+side)]]=np.deg2rad(val)
    net.arm_model.project(d);mujoco.mj_forward(m,d);net._update_weights(force=True)
    row=JOINTS.index('hip_rotation_r');value=float(net.arm_model.arms[row,ai]);weights=[]
    for cell in net.extra_pairs['HIP_Y_r']:
        post,pre=net.idx['MN_tfl_r'],net.idx[cell]
        edge=np.flatnonzero((net.compiled.post==post)&(net.compiled.pre==pre))[0]
        weights.append(float(net.compiled.gmax[edge]))
        assert net.compiled.g_max_non[post,pre]==net.compiled.gmax[edge]
    assert weights[0]==0. if value<0 else weights[1]==0.
    assert max(weights)>0
    out.append(dict(hip_flexion_deg=flex,rotation_arm_m=value,pf_positive_negative_weights=weights))
assert out[0]['rotation_arm_m']*out[1]['rotation_arm_m']<0
names=[p['name'] for p in net.net.populations];edges=[]
for c in net.net.connections:
    edges.append(dict(source=names[c['source']],target=names[c['destination']],g=float(c['params']['max_conductance']),reversal=float(c['params']['reversal_potential'])))
bi_edges=[e for e in edges if e['source'].startswith('PF_FUNC_BI_') and e['target'].startswith('MN_')]
bi_covered=sorted({e['target'][3:] for e in bi_edges})
expected=sorted(a for a,v in net.routing['actuator_mapping'].items() if v['confirmed_biarticular'])
assert bi_covered==expected,(bi_covered,expected)
for key,cells in net.extra_pairs.items():
    label,side=key.rsplit('_',1)
    for i,cell in enumerate(cells):
        relay=f'PF_IN_FUNC_{label}_{"P" if i==0 else "N"}_{side}'
        assert any(e['source']==cell and e['target']==relay and e['reversal']>0 for e in edges)
        assert any(e['source']==relay and e['target']==cells[1-i] and e['reversal']<0 for e in edges)
result=dict(machine='easteregg2',neurons=int(net.compiled.V.size),inputs=len(net.inputs),pairs=net.extra_pairs,
    biarticular_actuators=bi_covered,reference_fast_max_voltage_difference=maximum,
    verified_pose_dependent_switch=out,interneuron_lamination_verified=True,
    scope='Controller construction, geometry, cache refresh, and numerical neural parity; no gait acceptance claim')
(OUT/'controller_verification.json').write_text(json.dumps(result,indent=2))
(OUT/'functional_edges.json').write_text(json.dumps(edges,indent=2))
print(json.dumps(dict(neurons=result['neurons'],inputs=result['inputs'],biarticular_count=len(bi_covered),switch=out,parity=maximum),indent=2))
