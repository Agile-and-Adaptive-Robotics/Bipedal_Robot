"""Bounded full-Gait neural-family parity prerequisite, not a walking trial."""
import os,sys,time,json,copy,hashlib,argparse
from pathlib import Path
HERE=Path(__file__).resolve().parent;SPINAL=HERE.parents[1]/'spinal';REPO=HERE.parents[3]
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
sys.path[:0]=[str(SPINAL),str(SPINAL/'w2l_cpg')]
import numpy as np,mujoco
import build_network as B
from exact_fast_eval import ExactFastForward
from sparse_pairwise_eval import SparsePairwiseForward

def run(family,out,steps=256,sparse_existing=False):
 os.environ['AARL_NET']=family
 modelpath=REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
 m=mujoco.MjModel.from_xml_path(str(modelpath));acts=[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(m.nu)]
 tick=time.perf_counter();net=B.build(acts,dt=.002,interleg=True);buildwall=time.perf_counter()-tick
 if sparse_existing:
  saved=np.load(out/(family+'_parity.npz'));u=saved['inputs'];reference=saved['reference'];sparse=SparsePairwiseForward(copy.deepcopy(net.compiled));warm=SparsePairwiseForward(copy.deepcopy(net.compiled));warm.forward(u[0]);tick=time.perf_counter();observed=np.empty_like(reference)
  for i,x in enumerate(u):sparse.forward(x);observed[i]=sparse.V
  elapsed=time.perf_counter()-tick;error=float(np.max(abs(observed-reference)));old=json.loads((out/(family+'_parity.json')).read_text());old.update(sparse_pairwise_wall_s=elapsed,sparse_pairwise_max_voltage_error_mV=error,sparse_pairwise_bit_exact=error==0.,sparse_pairwise_speed_ratio=old['reference_neural_wall_s']/elapsed)
  assert np.isfinite(observed).all() and error<1e-10
  (out/(family+'_parity.json')).write_text(json.dumps(old,indent=2));print(json.dumps({k:v for k,v in old.items() if 'sparse_pairwise' in k},indent=2));return old
 original=copy.deepcopy(net.compiled);fast=ExactFastForward(copy.deepcopy(original))
 rng=np.random.default_rng(100031);u=rng.uniform(0.,.2,size=(steps,len(net.inputs)))
 for key,value in [('DRIVE',2.5),('POSTURE',.15)]:
  if key in net.inputs:u[:,net.inputs.index(key)]=value
 # Warm only an independent clone; reference and fast comparison start identically.
 warm=ExactFastForward(copy.deepcopy(original));tick=time.perf_counter();warm.forward(u[0]);coldwall=time.perf_counter()-tick
 errors=[];reference_states=np.empty((steps,len(original.V)));fast_states=np.empty_like(reference_states)
 tick=time.perf_counter()
 for i in range(steps):original.forward(u[i]);reference_states[i]=original.V
 refwall=time.perf_counter()-tick;tick=time.perf_counter()
 for i in range(steps):fast.forward(u[i]);fast_states[i]=fast.V
 fastwall=time.perf_counter()-tick
 error=float(np.max(np.abs(reference_states-fast_states)));finite=bool(np.isfinite(reference_states).all() and np.isfinite(fast_states).all())
 # Inspect existing corresponding neural weights; this does not fit either side.
 names=[None]*len(net.compiled.V)
 for name,index in net.idx.items():
  if names[index] is None:names[index]=name
 matrix=net.compiled.g_max_non;asym=[];lookup={n:i for i,n in enumerate(names)};seen=set()
 for post,pre in zip(*np.nonzero(matrix)):
  a,b=names[pre],names[post]
  def mirror(name):
   for x,y in (('_r','_l'),('_l','_r'),('_R','_L'),('_L','_R')):
    if name.endswith(x):return name[:-2]+y
   return name
  ma,mb=mirror(a),mirror(b)
  if ma in lookup and mb in lookup:
   pair=tuple(sorted(((int(post),int(pre)),(lookup[mb],lookup[ma]))))
   if pair in seen:continue
   seen.add(pair)
   delta=float(matrix[post,pre]-matrix[lookup[mb],lookup[ma]])
   if abs(delta)>1e-12:asym.append(dict(source=a,target=b,mirrored_source=ma,mirrored_target=mb,conductance_difference=delta))
 result=dict(machine='easteregg2',family=family,actuators=len(acts),neurons=len(original.V),inputs=len(net.inputs),nonzero_connections=int(np.count_nonzero(matrix)),steps=steps,dt_s=.002,
  finite=finite,max_voltage_error_mV=error,bit_exact=error==0.,parity_pass=finite and error<1e-10,reference_neural_wall_s=refwall,fast_neural_wall_s=fastwall,warm_neural_speed_ratio=refwall/fastwall,
  build_wall_s=buildwall,cold_fast_first_forward_wall_s=coldwall,seed=100031,model=str(modelpath),model_sha256=hashlib.sha256(modelpath.read_bytes()).hexdigest(),
  mirrored_nonzero_weight_differences=len(asym),first_weight_differences=asym[:12],weight_audit_definition='Unequal mirrored pairs counted once, including one-sided zero pathways; both suffix directions swapped.',
  scope='Reference/accelerated neuronal-state comparison on full92-actuator Gait names with seeded synthetic currents. No physical stepping, standing, balance, force-law, or full-trajectory acceptance is implied.')
 out.mkdir(parents=True,exist_ok=True);(out/(family+'_parity.json')).write_text(json.dumps(result,indent=2));np.savez_compressed(out/(family+'_parity.npz'),inputs=u,reference=reference_states,fast=fast_states,neuron_names=np.array(names),actuator_names=np.array(acts))
 print(json.dumps(result,indent=2),flush=True)
 if not result['parity_pass']:raise RuntimeError('Do not tune this family until acceleration parity is repaired')
 return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--family',choices=['s3k','w2lvar','syn6'],required=True);p.add_argument('--out',required=True);p.add_argument('--steps',type=int,default=256);p.add_argument('--sparse-existing',action='store_true');a=p.parse_args();run(a.family,Path(a.out),a.steps,a.sparse_existing)
