"""Causal ankle-block follow-up and explicitly labeled fan-in normalization."""
import argparse,os,json
from pathlib import Path
import numpy as np
import gait_standing_family as F
from muscle_map import classify
from gait_standing_objective import score
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');OUT=ROOT/'group_diagnosis'

def main(case):
 params=json.loads((ROOT/'rustuna_w2lvar_v1/BEST.json').read_text())['trial']['params'];scales={'IIIN:knee_ext':0.,'IIIN:hip_ext':0.};derivation=None
 if case=='ankle_block':scales.update({'IIIN:ankle_pf':0.,'IIIN:ankle_df':0.})
 elif case=='fanin_mean':
  os.environ['AARL_NET']='w2lvar';F.P.G['ii_set_pose']=1.;F.P.G['posture_to_mn']=0.;m=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));names=[m.actuator(i).name for i in range(m.nu)];n=F.B.build(names,dt=.002,interleg=True);F.tie_weights(n);ordered=F.canonical_neuron_names(n);bygroup={}
  for post,target in enumerate(ordered):
   if not target.startswith('MN_'):continue
   sources=[pre for pre in np.flatnonzero(n.compiled.g_max_non[post]) if ordered[pre].startswith('IIIN_')]
   if sources:bygroup.setdefault(classify(target[3:]).groups[0],[]).append(len(sources))
  if any(len(set(v))!=1 for v in bygroup.values()):raise ValueError('A group has nonuniform IIIN fan-in; do not silently average')
  derivation={group:values[0] for group,values in bygroup.items()};scales={'IIIN:'+g:1./count for g,count in derivation.items()}
 else:raise ValueError(case)
 r=F.run('w2lvar',OUT,tag=case+'_v1',gain=params['pf_mn_gain'],ib_group_gain=params['ib_group_gain'],standing_bias_scale=params['standing_bias_scale'],pathway_group_scales=scales,capture_neural_inputs=True);prefix=OUT/('w2lvar_'+case+'_v1');report=dict(case=case,assessment=score(r,prefix)[1],result=str(prefix.with_suffix('.result.json')),scales=scales,source_fanin=derivation,scope='Distinct opt-in standing mapping variant. Fan-in mean preserves every nonzero IIIN pathway but divides its target-group edge gains by actual number of incoming antagonist muscles. Engineering response to expanded muscle pools, not asserted author intent or unchanged literature replication. Ankle block is a causal ablation atop the two-group best; global gains/mechanics unchanged.')
 prefix.with_suffix('.assessment.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case',required=True);main(p.parse_args().case)
