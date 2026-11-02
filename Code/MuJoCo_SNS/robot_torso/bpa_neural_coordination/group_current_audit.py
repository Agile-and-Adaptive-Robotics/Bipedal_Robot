"""Actual input-aligned pre-update currents, per homologous motor neuron."""
import os,json
from pathlib import Path
import numpy as np
import gait_standing_family as F
from muscle_map import classify
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');OUT=ROOT/'group_diagnosis'

def main():
 F.P.G['ii_set_pose']=1.;F.P.G['posture_to_mn']=0.;m=F.mujoco.MjModel.from_xml_path(str(F.REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'));acts=[m.actuator(i).name for i in range(m.nu)];rows=[];targets=[]
 for family,label in (('s3k','s3k'),('syn6','syn6'),('w2lvar','w2lvar'),('w2lvar','fanin_mean')):
  p=OUT/(family+'_'+label+'_v1');r=json.loads(p.with_suffix('.result.json').read_text());params=json.loads((ROOT/('rustuna_'+family+'_v1/BEST.json')).read_text())['trial']['params'];os.environ['AARL_NET']=family;n=F.B.build(acts,dt=.002,interleg=True);F.tie_weights(n);c=n.compiled;names=F.canonical_neuron_names(n)
  for post,pre in zip(*np.nonzero(c.g_max_non)):
   if names[post].startswith('MN_') and names[pre].startswith('PF_'):c.g_max_non[post,pre]*=params['pf_mn_gain']
   if names[post].startswith('MN_') and names[pre].startswith('IBEXC'):c.g_max_non[post,pre]*=params['ib_group_gain']
   if names[post].startswith('MN_'):
    key=names[pre].split('_')[0]+':'+classify(names[post][3:]).groups[0]
    if key in (r.get('pathway_group_scales') or {}):c.g_max_non[post,pre]*=r['pathway_group_scales'][key]
  with np.load(p.with_suffix('.apparatus.npz')) as z:previous=z['neuro'][-251:-1].copy();following=z['neuro'][-250:].copy()
  with np.load(p.with_suffix('.inputs.npz')) as z:
   if z['names'].tolist()!=n.inputs:raise ValueError('Input identity mismatch')
   inputs=z['inputs'][-250:].copy()
  motor=[];error=0.
  for post,name in enumerate(names):
   if not name.startswith('MN_'):continue
   act=name[3:];sources=np.flatnonzero(c.g_max_non[post]);currents={};counts={}
   for pre in sources:
    cls=names[pre].split('_')[0];fraction=np.clip((previous[:,pre]-c.e_lo[post,pre])/(c.e_hi[post,pre]-c.e_lo[post,pre]),0,1);I=c.g_max_non[post,pre]*fraction*(c.del_e[post,pre]-previous[:,post]);currents[cls]=currents.get(cls,0.)+float(np.mean(I));counts[cls]=counts.get(cls,0)+1
   applied=inputs@c.input_connectivity[post];bias=float(c.i_b[post]);leak=-c.g_m[post]*(previous[:,post]-c.V_rest[post]);syn=sum(currents.values());predict=previous[:,post]+c.time_factor_membrane[post]*(leak+bias+applied+syn);error=max(error,float(np.max(abs(predict-following[:,post]))))
   port=n.inputs.index('POST_'+act);tonic=inputs[:,port]*c.input_connectivity[post,port]
   motor.append(dict(actuator=act,primary_group=classify(act).groups[0],all_source_functional_groups=classify(act).groups,mean_voltage_mV=float(following[:,post].mean()),mean_normalized_motor_command=float(np.clip(following[:,post]/F.R.E_HI,0,1).mean()),synaptic_current_nA=currents,incoming_edges_by_class=counts,external_tonic_POST_nA=float(tonic.mean()),external_total_nA=float(applied.mean()),intrinsic_bias_nA=bias,leak_current_nA=float(leak.mean())))
  # The exact instantaneous membrane current also depends on V. Averaged
  # current above is descriptive; verify the full time-resolved balance below.
  exact_error=0.
  for row in motor:
   post=n.idx['MN_'+row['actuator']];Isyn=np.zeros(250)
   for pre in np.flatnonzero(c.g_max_non[post]):Isyn+=c.g_max_non[post,pre]*np.clip((previous[:,pre]-c.e_lo[post,pre])/(c.e_hi[post,pre]-c.e_lo[post,pre]),0,1)*(c.del_e[post,pre]-previous[:,post])
   pred=previous[:,post]+c.time_factor_membrane[post]*(-c.g_m[post]*(previous[:,post]-c.V_rest[post])+c.i_b[post]+inputs@c.input_connectivity[post]+Isyn);exact_error=max(exact_error,float(np.max(abs(pred-following[:,post]))))
  if exact_error>1e-9:raise ValueError('Input/voltage temporal alignment or missing motor current: '+str(exact_error))
  groups={}
  for group in sorted({a['primary_group'] for a in motor}):
   members=[a for a in motor if a['primary_group']==group];classes=sorted({k for a in members for k in a['synaptic_current_nA']});groups[group]=dict(muscles=len(members),mean_voltage_mV=float(np.mean([a['mean_voltage_mV'] for a in members])),mean_tonic_nA=float(np.mean([a['external_tonic_POST_nA'] for a in members])),mean_synaptic_nA={k:float(np.mean([a['synaptic_current_nA'].get(k,0) for a in members])) for k in classes})
  targets.append(np.array(r['standing_activation_targets']));rows.append(dict(family=family if label==family else 'w2lvar_meanIIIN',motor_neurons=motor,groups=groups,membrane_balance_max_error_mV=exact_error,standing_solver_diagnostics=r['standing_solver_diagnostics'],parameters=params))
 report=dict(machine='easteregg2',families=rows,shared_static_targets_bit_exact=all(np.array_equal(targets[0],a) for a in targets[1:]),scope='Same corrected-pose static gravity/contact fit; actual250late input/previous-voltage pairs. Per-class currents retain signed conductance/reversal interactions. Matched homolog names across families, with time-resolved membrane-balance gate. Joint torque comparison requires an equality-aware common mechanical projection; this report does not treat current as torque.')
 (OUT/'PER_MOTOR_CURRENT_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(shared_static_targets_bit_exact=report['shared_static_targets_bit_exact'],families=[dict(family=r['family'],balance_error=r['membrane_balance_max_error_mV'],groups=r['groups']) for r in rows]),indent=2))

if __name__=='__main__':main()
