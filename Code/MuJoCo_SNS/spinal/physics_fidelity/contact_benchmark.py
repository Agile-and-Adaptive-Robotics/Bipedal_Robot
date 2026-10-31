"""Explicit ball/plane contact sensitivity; no walker parameter adoption."""
import os,sys,json,time
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np,mujoco
OUT=Path(__file__).resolve().parent/'runs/audit_20261030'
OUT.mkdir(parents=True,exist_ok=True)
variants={
 'baseline':([.02,1.],[.9,.95,.001,.5,2.]),
 'lessviscous_existing':([.01,.7],[.9,.95,.001,.5,2.]),
 'legacy_nonlinear_actual':([.02,1.],[.9,.9,.001,.003,2.]),
 'depth_width_corrected':([.02,1.],[.9,.95,.003,.5,2.]),
 'constant_impedance_control':([.02,1.],[.9,.9,.003,.5,2.])}
rows={};traces={}
for name,(ref,imp) in variants.items():
 xml=f'<mujoco><option timestep=".002" integrator="implicitfast"/><worldbody><geom name="floor" type="plane" size="1 1 .1"/><body pos="0 0 .15"><freejoint/><geom name="ball" type="sphere" size=".05" mass="1"/></body></worldbody><contact><pair geom1="floor" geom2="ball" solref="{" ".join(map(str,ref))}" solimp="{" ".join(map(str,imp))}"/></contact></mujoco>'
 m=mujoco.MjModel.from_xml_string(xml);d=mujoco.MjData(m);force=np.zeros(6);trace=np.zeros((1000,4));start=time.perf_counter()
 for k in range(1000):
  mujoco.mj_step(m,d);f=0.;penetration=0.
  for j in range(d.ncon):
   mujoco.mj_contactForce(m,d,j,force);f+=force[0];penetration=max(penetration,-d.contact[j].dist)
  trace[k]=[d.time,d.qpos[2],f,penetration]
 rows[name]=dict(solref=ref,solimp=imp,wall_s=time.perf_counter()-start,finite=bool(np.all(np.isfinite(trace))),peak_contact_force_N=float(trace[:,2].max()),peak_penetration_m=float(trace[:,3].max()),final_height_m=float(trace[-1,1]))
 traces[name]=trace
 assert rows[name]['finite']
delta=float(np.max(abs(traces['legacy_nonlinear_actual']-traces['constant_impedance_control'])))
assert delta<1e-10,delta
rows['constant_impedance_trace_max_difference']=delta
(OUT/'CONTACT_RESULTS.json').write_text(json.dumps(rows,indent=2));np.savez_compressed(OUT/'contact_traces.npz',**traces)
print(json.dumps(rows,indent=2))
