"""Small sensory/command ablations, with actual neural traces and no body claims."""
import argparse,json,os,sys,time
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
from replay_adapter import Replay
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    cfg=json.loads(a.manifest.read_text());r=Replay(cfg)
    u=r.u
    stats={name:dict(minimum=float(u[:,i].min()),maximum=float(u[:,i].max()),rms=float(np.sqrt(np.mean(u[:,i]**2))))for name,i in r.net.command_runtime.ports.items() if np.any(u[:,i])}
    (a.output/'INPUT_STATS.json').write_text(json.dumps(stats,indent=2))
    original=dict(command_F=cfg['command_F'],command_E=cfg['command_E'])
    cases=[('no_sensory_no_command',0.,0.,0.,0.,0.),('no_sensory_command',.15,.15,0.,0.,0.),('reference',.15,.15,1.,1.,1.)]
    for aff in range(3):
        for gain in (8.,32.,96.):
            values=[0.,0.,0.];values[aff]=gain
            cases.append((('Ia','II','Ib')[aff]+'_'+str(int(gain)),.15,.15,*values))
    cases += [('all_32',.15,.15,32.,32.,32.),('flex_drive_only',.15,0.,1.,1.,1.)]
    rows=[]
    for tag,f,e,ia,ii,ib in cases:
        r.manifest.update(command_F=f,command_E=e)
        parameters={**cfg.get('fixed_parameters',{}),'Ia_RG':ia,'II_RG':ii,'Ib_RG':ib}
        result=r.evaluate(parameters,time.time()+180.)
        row=dict(case=tag,parameters=parameters,command_F=f,command_E=e,**result)
        rows.append(row);(a.output/'RESULTS.json').write_text(json.dumps(rows,indent=2))
        np.savez_compressed(a.output/(tag+'.npz'),t=r.t,rg=r.rg_trace,pf=r.pf_trace,mn=r.output,desired=r.desired,rg_names=np.array([r.names[i]for i in r.rg]),pf_names=np.array([r.names[i]for i in r.pf]),acts=np.array(r.acts))
        print(json.dumps(dict(case=tag,**result)),flush=True)
    (a.output/'MANIFEST.json').write_text(json.dumps(cfg,indent=2))
if __name__=='__main__':main()
