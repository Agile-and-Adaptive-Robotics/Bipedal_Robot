"""Check the new fit's fixed-conductance sparse SNS engine and encoder parity."""
import os,sys,json,time
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
H=Path(__file__).resolve().parent;campaign=os.environ.get('HUNT_CAMPAIGN','fourcell_hunt_20261029_v1')
if Path(campaign).name!=campaign:raise ValueError('Campaign must be a directory name under hunt_training/runs')
root=H/'runs'/campaign
import numpy as np,mujoco
import motor_fit as M,controller as C
real=C.build
def build(*args,**kwargs):
    return real(*args,options=dict(weights_directory=str(root/'weights'),biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True),**kwargs)
C.build=build
z,m,d,a,t,q,v,target,pf=M.setup(1);z2,m2,d2,b,*_=M.setup(1)
a.apply_weights(1.);b.apply_weights(1.)
rng=np.random.default_rng(291026);inputs=np.zeros((500,len(a.inputs)));inputs[:,a.ports['DRIVE']]=2.5
for side in ('r','l'):
    for f in range(6):inputs[:,a.ports[f'HUNT_AF_{f}_{side}']]=rng.uniform(0.,1.,500)
inputs[:50,a.ports['HUNT_INIT_E_r']]=2.;inputs[:50,a.ports['HUNT_INIT_F_l']]=2.
fast=np.empty((len(inputs),len(a.compiled.V)));native=np.empty_like(fast)
start=time.perf_counter()
for k,u in enumerate(inputs):a.compiled.forward(u);fast[k]=a.compiled.V
fast_s=time.perf_counter()-start;start=time.perf_counter()
for k,u in enumerate(inputs):b.compiled.c.forward(u);native[k]=b.compiled.V
native_s=time.perf_counter()-start
error=float(np.max(abs(fast-native)));assert error<1e-10,error
encoder_error=0.
for k in np.linspace(0,len(t)-1,12).astype(int):
    d.qpos[:]=q[k];a.arm_model.project(d);d.qvel[:]=v[k];d.act[:]=target[k];mujoco.mj_forward(m,d)
    u=a.make_inputs();M.afferents(a,u,d)
    ranges=m.actuator_lengthrange;length=np.clip((d.actuator_length-ranges.mean(axis=1))/np.maximum(np.ptp(ranges,axis=1)/2,.001),-1,1)
    velocity=np.clip(d.actuator_velocity/M.P.AFF['ia_vel_ref'],-1,1);force=np.clip(abs(d.actuator_force)/a.Fmax,0,1.5)
    for i,act in enumerate(z['act_names']):
        info=a.muscles[act];stance=float(np.clip(a.compiled.V[a.idx['RG_E_'+info.side]]/5.,0.,1.));ext=info.groups[0] in M.EXTENSOR_STANCE_GROUPS
        drive=float(np.clip(a.compiled.V[a.idx['DRIVE']]/5.,0.,1.));ia=M.P.AFF['ia_gain']*(1+M.P.MOD['ia']*drive)*(.5+.5*stance if ext else 1.)
        ii=M.P.AFF['ii_gain']*(1+M.P.MOD['ii']*drive)*(.3+.7*stance if ext else .5);ib=M.P.AFF['ib_gain']*(1-.7*stance if ext else 1.)
        for name,value in [('Ia',max(ia*velocity[i],0.)),('II',max(M.P.AFF['i0_ii']*(.3+.7*stance)+ii*length[i],0.)),('Ib',.5*ib*force[i])]:encoder_error=max(encoder_error,abs(u[a.ports[name+'_'+act]]-value))
assert encoder_error<1e-12,encoder_error
result=dict(neurons=len(a.compiled.V),steps=len(inputs),max_voltage_error=error,bit_exact=bool(np.array_equal(fast,native)),fast_s=fast_s,native_s=native_s,ratio=native_s/fast_s,encoder_max_error=encoder_error,
    scope='Fixed-conductance neural updates with identical changing inputs; excludes physics, per-step conductance updates, compilation, and logging. This is a single timing measurement, not a controlled hardware benchmark. Encoder check covers midpoint-length reference and gamma_coact=0.')
(root/'acceleration_check.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
