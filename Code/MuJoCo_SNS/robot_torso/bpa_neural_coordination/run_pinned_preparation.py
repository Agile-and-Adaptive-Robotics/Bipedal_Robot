"""Opt-in prepared actuator-state probes; historical main warmup is unchanged."""
import argparse,json,copy
from pathlib import Path
import gait_standing_family as F
from gait_standing_objective import score

ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def main(case):
 out=ROOT/'pinned_preparation';out.mkdir(exist_ok=True)
 nominal=ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json'
 reference=out/'PREPARATION_REFERENCE.json';cfg=copy.deepcopy(json.loads(nominal.read_text()));cfg['pressure_history_capacity_s']=.016
 if reference.exists() and json.loads(reference.read_text())!=cfg:raise ValueError('Preserve prepared-state reference')
 reference.write_text(json.dumps(cfg,indent=2),encoding='utf-8')
 params=json.loads((ROOT/'rustuna_syn6_v1/BEST.json').read_text())['trial']['params']
 allowed={'nominal','nominal_fresh','nominal_replay','harsh_live','harsh_replay','harsh_late_live','harsh_late_replay','harsh_zero_live','harsh_zero_replay','harsh_cold_live','harsh_cold_replay','bio_nominal','bio_nominal_fresh','rot0p8','rot0p5','lat0p8','vertical0p8','selected_fresh','selected20'}
 if case not in allowed:raise ValueError('Unknown bounded case')
 tag='prepared_'+case+'_v1';path=out/('syn6_'+tag+'.result.json')
 if path.exists():raise FileExistsError('Completed case already preserved: '+str(path))
 actual=ROOT/'virtual_bpa_design/HELDOUT_SLOW_WEAK_MEMORY.json' if case.startswith('harsh') else nominal
 replay=out/'syn6_prepared_nominal_v1.inputs.npz' if case.endswith('replay') else None
 kwargs=dict(pinned_actuator_settling=True,settling_reference_configuration=reference,response_onset_physics_s=2. if 'late' in case else 0.)
 if case.startswith('bio_nominal'):actual=None;kwargs['settling_reference_configuration']=None
 if case.startswith(('harsh_zero','harsh_cold')):kwargs['release_pressure_reset']=True
 if case.startswith('harsh_cold'):kwargs['release_history_reset']=True
 if case.startswith('rot'):kwargs['rig_rotational_scale']=.8 if '0p8'in case else .5
 if case=='lat0p8':kwargs.update(rig_rotational_scale=.5,rig_lateral_scale=.8)
 if case=='vertical0p8':kwargs.update(rig_rotational_scale=.5,rig_lateral_scale=.8,rig_vertical_scale=.8)
 if case.startswith('selected'):
  scales=json.loads((out/'SELECTED_RUNG.json').read_text())['rig_scales'];kwargs.update(rig_rotational_scale=scales['rotational'],rig_lateral_scale=scales['lateral'],rig_vertical_scale=scales['vertical'])
 r=F.run('syn6',out,duration=20. if case=='selected20' else 8.,tag=tag,gain=params['pf_mn_gain'],ib_group_gain=params['ib_group_gain'],standing_bias_scale=params['standing_bias_scale'],bpa_configuration=actual,sensory_replay=replay,**kwargs)
 assessment=score(r,path.with_suffix('').with_suffix(''))[1]
 ledger=out/'SUMMARY.json';previous=json.loads(ledger.read_text()) if ledger.exists() else dict(machine='easteregg2',controller_parameters=params,cases={},scope='Common nominal pinned actuator-state preparation. Harsh response changes begin at release or2s of free physics. Live/replay preparation is identical; all initialized states and modeled spring energies must match before comparison. Historical main warmup unchanged; rig and knee protection retained.')
 previous['cases'][case]=dict(result=str(path),assessment=assessment,release_state=r.get('release_state_file'),failure=r['failure'],upward_rig_bw=r.get('mean_upward_rig_support_bw'),downward_rig_bw=r.get('mean_downward_rig_force_bw'),mesh_range_m=r.get('actual_visual_foot_mesh_clearance_range_m'),knee_max_deg=r.get('max_knee_angle_deg'))
 ledger.write_text(json.dumps(previous,indent=2),encoding='utf-8');print(json.dumps(previous['cases'][case],indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case',required=True);main(p.parse_args().case)
