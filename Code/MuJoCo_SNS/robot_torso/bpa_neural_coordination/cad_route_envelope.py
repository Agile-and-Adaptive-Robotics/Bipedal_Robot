"""CAD route feasibility from existing motion; no dynamic replay or hardware guesses."""
import os,sys,json,hashlib,argparse,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROBOT=HERE.parent
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent);sys.path.insert(0,str(ROBOT))
import numpy as np,mujoco
from dynamic_routes import DynamicRoutes

def main(out):
 source=ROBOT/'runs/verified65_exact_fast.npz';modelpath=ROBOT/'ah001_neural.xml';routespath=ROBOT/'robotbody_routes.json'
 z=np.load(source);q=z['qpos'];m=mujoco.MjModel.from_xml_path(str(modelpath));d=mujoco.MjData(m)
 if q.shape[1]!=m.nq:raise ValueError('Historical motion and current CAD nq mismatch: refuse reinterpretation')
 routes=json.loads(routespath.read_text());dynamic=DynamicRoutes(m,routes);indices=np.linspace(0,len(q)-1,min(401,len(q)),dtype=int);lengths=np.empty((len(indices),m.nu));names=[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(m.nu)];tid=m.actuator_trnid[:,0]
 for k,i in enumerate(indices):d.qpos[:]=q[i];mujoco.mj_forward(m,d);dynamic.update(d);lengths[k]=d.ten_length[tid]
 minimum=lengths.min(0);maximum=lengths.max(0);e=.25;budget=(minimum-(1-e)*maximum)/e
 rows=[dict(name=n,route_min_m=float(a),route_max_m=float(b),route_excursion_m=float(b-a),direct_1to1_max_contraction_assumed=e,
  maximum_fixed_nonmuscle_length_m=float(c),feasible_with_zero_cable=bool(c>=0),feasible_if_two_25mm_fittings_assumed=bool(c>=.05),
  assumed_series_extension_m=0.,scope='Zero-series-stretch length-envelope feasibility only. No sleeve/cable assignment, force-dependent compliance, force capacity, mounting, or pressure dynamics identified.') for n,a,b,c in zip(names,minimum,maximum,budget)]
 out.mkdir(parents=True,exist_ok=True);np.savez_compressed(out/'CAD_ROUTE_ENVELOPES.npz',lengths=lengths,indices=indices,names=np.array(names))
 result=dict(machine='easteregg2',model=str(modelpath),source_motion=str(source),source_motion_scope='Historical hybrid allocator/reference-tension motion used solely as a geometric pose library; not a neural-only BPA trial',poses=len(indices),actuators=m.nu,
  source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (modelpath,routespath,source,Path(__file__))},
  zero_cable_length_feasible=int(np.sum(budget>=0)),two_25mm_fitting_length_feasible=int(np.sum(budget>=.05)),routes=rows,
  formula='At zero series stretch, Lroute=Lfixed+Lrest*(1-e), e in[0,.25], gives Lfixed <= (Lmin−.75Lmax)/.25. Finite force-dependent cable stretch changes this bound. Does not prove full compliance, force, or packaging feasibility.',
  next='Use per-route resting length/cable/parallel-count/compliance assignments with this envelope. Do not silently transfer historical Routing27 designs or human optimal fiber lengths.')
 (out/'CAD_ROUTE_ENVELOPE_REPORT.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k not in ('routes','source_sha256')},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args();main(Path(a.out))
