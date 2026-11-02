"""Source-graph fan-in and late current audit at frozen standing finalists."""
import os,json
from pathlib import Path
import numpy as np
import gait_standing_family as F
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def main():
 F.P.G['ii_set_pose']=1.;F.P.G['posture_to_mn']=0.;rows=[]
 m=F.mujoco.MjModel.from_xml_path(str(F.REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'));acts=[m.actuator(i).name for i in range(m.nu)]
 for family in ('s3k','w2lvar'):
  os.environ['AARL_NET']=family;best=json.loads((ROOT/('rustuna_'+family+'_v1/BEST.json')).read_text());params=best['trial']['params'];n=F.B.build(acts,dt=.002,interleg=True);F.tie_weights(n);c=n.compiled;names=F.canonical_neuron_names(n)
  for post,pre in zip(*np.nonzero(c.g_max_non)):
   if names[post].startswith('MN_') and names[pre].startswith('PF_'):c.g_max_non[post,pre]*=params['pf_mn_gain']
   if names[post].startswith('MN_') and names[pre].startswith('IBEXC'):c.g_max_non[post,pre]*=params['ib_group_gain']
  path=Path(best['trial']['diagnostics']['result']).with_suffix('').with_suffix('.apparatus.npz');z=np.load(path);v=z['neuro'][-250:];assert v.shape[1]==len(names);totals={};fan=[]
  for post,name in enumerate(names):
   if not name.startswith('MN_'):continue
   sources=np.flatnonzero(c.g_max_non[post]);source_rows=[]
   for pre in sources:
    cls=names[pre].split('_')[0];fraction=np.clip((v[:,pre]-c.e_lo[post,pre])/(c.e_hi[post,pre]-c.e_lo[post,pre]),0,1);current=float(np.mean(c.g_max_non[post,pre]*fraction*(c.del_e[post,pre]-v[:,post])));totals[cls]=totals.get(cls,0)+current/len(acts)
    if cls=='IIIN':source_rows.append(dict(source=names[pre],max_conductance_uS=float(c.g_max_non[post,pre]),mean_current_nA=current))
   if source_rows:fan.append(dict(target=name,IIIN_incoming_edges=len(source_rows),IIIN_total_gmax_uS=sum(r['max_conductance_uS'] for r in source_rows),IIIN_mean_current_nA=sum(r['mean_current_nA'] for r in source_rows),mean_motor_voltage_mV=float(v[:,post].mean()),sources=source_rows))
  rows.append(dict(family=family,parameters=params,mean_per_MN_signed_current_nA=totals,IIIN_fan_in=sorted(fan,key=lambda r:r['IIIN_mean_current_nA']),nonzero_edges=int(np.count_nonzero(c.g_max_non)),scope='Analytical synaptic currents at saved late-standing states, rebuilt matching reflex flags and mirrored pooling. External POST input and membrane currents are separate. Source-derived current attribution, not a literature claim or intervention.'))
 report=dict(machine='easteregg2',families=rows,interpretation='Turning PF and aggregate-Ib branches down is not sufficient for w2lvar: its multi-muscle antagonist-II inhibitory fan-in remains. Any pooling/gain intervention must be opt-in, bilateral, and separately labeled; this audit does not silently alter connectivity.')
 (ROOT/'OTHER_FAMILY_CURRENT_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps([{k:v for k,v in r.items() if k!='IIIN_fan_in'} for r in rows],indent=2))

if __name__=='__main__':main()
