"""Make the legacy MuJoCo/OpenSim axis adapter explicit at neutral posture."""
import os,sys,json
from pathlib import Path
HERE=Path(__file__).resolve().parent;S=HERE.parent;sys.path[:0]=[str(S),str(HERE)]
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import mujoco
import bsolve_ik as B
from mechanics import CoupledArms,JOINTS
m,d=B.build_air_model(False);arms=CoupledArms(m);d.qpos[arms.addresses]=0.;arms.project(d);mujoco.mj_forward(m,d)
records={}
for name,jid in zip(JOINTS,arms.jids):
    v=d.xaxis[jid];isb=np.array([v[0],v[2],-v[1]])
    action='positive native coordinate'
    if name.startswith('hip_adduction'):
        action='adduction' if isb[0]*(1 if name.endswith('_r') else -1)>0 else 'abduction'
    elif name.startswith('hip_rotation'):
        action='internal rotation' if isb[1]*(1 if name.endswith('_r') else -1)>0 else 'external rotation'
    elif name.startswith('subtalar'):
        action='inversion' if isb[0]*(1 if name.endswith('_r') else -1)>0 else 'eversion'
    elif name=='lumbar_rotation':action='left axial rotation' if isb[1]>0 else 'right axial rotation'
    elif name=='lumbar_bending':action='right lateral bending' if isb[0]>0 else 'left lateral bending'
    elif name=='lumbar_extension':action='extension' if isb[2]>0 else 'flexion'
    records[name]=dict(positive_joint_axis_in_isb=isb.tolist(),positive_native_action=action)
out=S/'reports_20261028_functional_pf'
(out/'axis_conventions.json').write_text(json.dumps(dict(frame='x anterior, y superior, z right',
    adapter='legacy MuJoCo world (x,y,z) -> ISB (x,z,-y)',neutral_axes=records,
    limitation='Joint axes rotate with posture. Subtalar inversion/eversion is about its actual oblique axis, not an artificially imposed pure x hinge.'),indent=2))
print(json.dumps({k:v for k,v in records.items() if k.startswith(('hip_rotation','hip_adduction','subtalar','lumbar'))},indent=2))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
z=np.load(out/'routing_samples.npz');names=z['actuator_names'].tolist();rho=z['sweep_moment_arms'];angles=z['hip_flexion_deg']
factor=1 if records['hip_rotation_r']['positive_native_action']=='internal rotation' else -1
plt.rcParams.update({'font.family':'Arial','font.size':11})
fig,ax=plt.subplots(figsize=(8,5),constrained_layout=True)
for name,label,style,color in [('glut_med1_r','Gluteus medius 1','o-','#0077BB'),('tfl_r','Tensor fasciae latae','s--','#CC6677'),('sar_r','Sartorius','D-.','#228833')]:
    ax.plot(angles,1000*factor*rho[:,JOINTS.index('hip_rotation_r'),names.index(name)],style,color=color,label=label,markersize=4)
ax.axhline(0,color='.3',lw=1);ax.set(xlabel='Hip flexion (degrees)',ylabel='Internal-rotation moment arm (mm)',title='Gait2392: hip rotation recruitment changes with posture')
ax.legend();fig.savefig(out/'hip_rotation_routing.png',dpi=220);plt.close(fig)
(out/'hip_rotation_routing.alt.txt').write_text('Three measured Gait2392 hip-rotation moment-arm curves across minus 20 to 100 degrees of hip flexion. Gluteus medius compartment 1, tensor fasciae latae, and sartorius cross zero, demonstrating that fixed agonist/antagonist labels do not cover the sampled range. Internal-rotation sign is derived from the neutral anatomical axis audit.')
