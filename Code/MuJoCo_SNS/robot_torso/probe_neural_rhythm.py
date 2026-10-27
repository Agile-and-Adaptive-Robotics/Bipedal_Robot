"""Staged neural-only rhythm check before further robot stepping trials."""
import os,sys,json
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import mujoco
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'spinal'));sys.path.insert(0,str(HERE.parent/'spinal/w2l_cpg'))
import build_network as bn
import params
from exact_fast_eval import ExactFastForward
def main():
    model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_neural.xml'))
    names=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(model.nu)]
    results=[]
    for drive in (.05,.1,.2,.4,.8,1.5,3.):
        params.G.update(c1_gain=.15,v3_gain=0.,heel_rge=.15,toe_rge=.05)
        net=bn.build(names,dt=.002);net.compiled=ExactFastForward(net.compiled)
        net.compiled.V[net.idx['RG_F_r']]=5.
        net.compiled.V[net.idx['RG_E_l']]=5.
        u=net.make_inputs();u[net.input_index('DRIVE')]=drive;u[net.input_index('POSTURE')]=.15
        trace=[]
        for step in range(6000):
            net.step(u)
            if step>=1000 and step%5==0:
                v=net.compiled.V
                trace.append([v[net.idx['RG_F_'+s]]-v[net.idx['RG_E_'+s]] for s in ('r','l')])
        trace=np.array(trace)
        crossings=[int(np.sum((x[1:]>0)&(x[:-1]<=0))) for x in trace.T]
        result=dict(drive=drive,flexor_crossings=crossings,difference_excursion_mV=np.ptp(trace,axis=0).tolist(),
                    correlation=float(np.corrcoef(trace.T)[0,1]) if np.min(np.std(trace,axis=0))>1e-10 else None)
        results.append(result);print(result,flush=True)
    (HERE/'neural_rhythm_probe.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':main()
