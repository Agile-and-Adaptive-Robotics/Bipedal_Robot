"""Original full-gain standing PF/current attribution from preserved neural states."""
import os,json
import numpy as np
import gait_standing_family as F
from muscle_map import classify
from pathlib import Path
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');OUT=ROOT/'group_diagnosis'

def main():
 F.P.G['ii_set_pose']=1.;F.P.G['posture_to_mn']=0.;m=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));acts=[m.actuator(i).name for i in range(m.nu)];rows=[]
 for family in ('s3k','syn6','w2lvar'):
  p=ROOT/'standing'/(family+'_baseline_v2');r=json.loads(p.with_suffix('.result.json').read_text());matched=json.loads((OUT/(family+'_'+family+'_v1.result.json')).read_text());targets=r.get('standing_activation_targets',matched['standing_activation_targets']);os.environ['AARL_NET']=family;n=F.B.build(acts,dt=.002,interleg=True);F.tie_weights(n);c=n.compiled;names=F.canonical_neuron_names(n)
  with np.load(p.with_suffix('.apparatus.npz')) as z:before=z['neuro'][-251:-1].copy();after=z['neuro'][-250:].copy()
  motor=[]
  for post,name in enumerate(names):
   if not name.startswith('MN_'):continue
   totals={};series=np.zeros(250)
   for pre in np.flatnonzero(c.g_max_non[post]):
    I=c.g_max_non[post,pre]*np.clip((before[:,pre]-c.e_lo[post,pre])/(c.e_hi[post,pre]-c.e_lo[post,pre]),0,1)*(c.del_e[post,pre]-before[:,post]);key=names[pre].split('_')[0];totals[key]=totals.get(key,0)+float(I.mean());series+=I
   inferred=(after[:,post]-before[:,post])/c.time_factor_membrane[post]+c.g_m[post]*(before[:,post]-c.V_rest[post])-c.i_b[post]-series;expected=6*targets[acts.index(name[3:])]*r.get('standing_bias_scale',1.);error=float(np.max(abs(inferred-expected)))
   if error>1e-9:raise ValueError('Original inferred tonic does not match runner standing law: '+str((family,name,error)))
   motor.append(dict(actuator=name[3:],primary_group=classify(name[3:]).groups[0],currents_nA=totals,tonic_input_inferred_nA=float(inferred.mean()),expected_tonic_nA=expected,tonic_law_error_nA=error,mean_voltage_mV=float(after[:,post].mean())))
  rows.append(dict(family=family,motor_neurons=motor))
 report=dict(machine='easteregg2',families=rows,scope='Preserved original all-global-gain1 standing runs. Signed synaptic currents evaluated on actual previous voltages. External tonic current is inferred from the membrane balance then independently checked against the runner6×standing-target bias law, not presented as directly logged input. No body force or neural voltage is modified.')
 (OUT/'ORIGINAL_FULL_GAIN_MOTOR_CURRENTS.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print('Original full-gain motor-current and tonic-law audit passed')

if __name__=='__main__':main()
