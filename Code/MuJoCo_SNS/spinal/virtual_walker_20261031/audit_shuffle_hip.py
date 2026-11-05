"""Read saved neural voltages to attribute hip-MN synaptic current."""
import argparse,json
from pathlib import Path
import numpy as np
from native_kernel import NativeController
from recruitment import apply_recruitment
def main():
 p=argparse.ArgumentParser();p.add_argument('folder',type=Path);a=p.parse_args();meta=json.loads((a.folder.parent/'IMPORT_AUDIT.json').read_text());result=json.loads((a.folder/'RESULT.json').read_text());z=np.load(a.folder/'trajectory.npz');mask=z['t']>=2.;v=z['neuro'][mask]
 c=NativeController(meta['source_snapshot'],meta['muscles'],ankle_pf_crossed=result['ankle_pf_crossed'],symmetric_rg_ext_tonic=result['symmetric_rg_ext_tonic_nA'],neural_time_scale=result['neural_time_scale'],contact_pf_scale=result['contact_pf_scale']);apply_recruitment(c,result['pf_mn_recruitment_scale'],result.get('toe_mn_gmax'),result['hip_pf_mn_scale'],result.get('hip_flexor_pf_mn_scale',1.))
 rows=[]
 for pre,post,g,threshold,saturation,reversal in c.edges:
  pre,post=int(pre),int(post)
  if c.names[post]not in('L Hip MN flx','L Hip MN ext','R Hip MN flx','R Hip MN ext'):continue
  conductance=g*np.clip((v[:,pre]-threshold)/(saturation-threshold),0.,1.);current=conductance*(reversal-v[:,post]);rows.append(dict(source=c.names[pre],target=c.names[post],gmax_uS=float(g),mean_current_nA=float(current.mean()),mean_abs_current_nA=float(abs(current).mean()),maximum_abs_current_nA=float(abs(current).max())))
 rows.sort(key=lambda x:x['mean_abs_current_nA'],reverse=True);(a.folder/'HIP_MN_CURRENT_AUDIT.json').write_text(json.dumps(dict(rows=rows,scope='Synaptic currents evaluated at saved interval-end voltages; direct adapter/bias/ionic currents are not included; not exact previous substep currents'),indent=2));print(json.dumps(rows[:16],indent=2))
if __name__=='__main__':main()
