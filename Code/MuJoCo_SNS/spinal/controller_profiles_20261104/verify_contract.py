"""Graph-only strict contact fault injection for selected Gait2392 profiles."""
import argparse,json,os,sys
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import numpy as np
from controller_profiles_20261104 import profile as P
from controller_profiles_20261104.prepare import audit

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    root=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104')
    folders=dict(hunt=root/'recalibrated_v2/hunt',s3k=root/'intrinsic_seed_v3/s3k',w2lvar=root/'commissural_fit_v4/w2lvar',syn6=root/'commissural_fit_v4/syn6')
    rows=[]
    for family,folder in folders.items():
        manifest=json.loads((folder/'training/B/MANIFEST.json').read_text());acts=np.load(manifest['targets'])['act_names'].tolist();n=P.build(family,acts,options=manifest.get('hunt_options'),config=manifest['controller_profile'])
        baseline=audit(n)
        connection=next(c for c,source,target in P.connections(n)if source=='HEEL_r'and target=='InE_r');rev=connection['params']['reversal_potential'];connection['params']['reversal_potential']=-8.
        try:single=any(f['issue']=='inhibitory_contact_entry'and f.get('target')=='InE_r'for f in P.contact_flags(n)['flags'])
        finally:connection['params']['reversal_potential']=rev
        import build_network as B
        n.net.add_connection(B._syn(.1,True),'TOE_r','RG_F_r');injected=n.net.connections[-1]
        try:toe=any(f['issue']=='toe_contact_at_RG'for f in P.contact_flags(n)['flags'])
        finally:n.net.connections.remove(injected)
        heel_to_wrong=any(source=='HEEL_r'and target=='InF_r'for _,source,target in P.connections(n))
        result=dict(family=family,static_pass=baseline['static_pass'],single_inhibitory_contact_detected=single,toe_RG_detected=toe,ipsilateral_heel_to_InF_present=heel_to_wrong,final_contact_flags=P.contact_flags(n)['flags'],neural_steps=0,physics_steps=0)
        result['passed']=bool(baseline['static_pass']and single and toe and not heel_to_wrong and not result['final_contact_flags']);rows.append(result);(a.output/'RESULTS.json').write_text(json.dumps(rows,indent=2));print(json.dumps(result),flush=True)
    if not all(row['passed']for row in rows):raise SystemExit(2)
if __name__=='__main__':main()
