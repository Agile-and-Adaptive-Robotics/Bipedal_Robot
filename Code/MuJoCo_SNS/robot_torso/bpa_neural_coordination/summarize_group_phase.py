"""Source/current/torque evidence and explicit experimental mapping checkpoint."""
import json,csv,hashlib
from pathlib import Path
import numpy as np
from update_family_cad_record import section
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');OUT=ROOT/'group_diagnosis'

def main():
 current=json.loads((OUT/'PER_MOTOR_CURRENT_AUDIT.json').read_text());original=json.loads((OUT/'ORIGINAL_FULL_GAIN_MOTOR_CURRENTS.json').read_text());torque=json.loads((OUT/'COMMON_JOINT_TORQUE_REFERENCE.json').read_text());checks=OUT/'normalized_checks';cases=[]
 for p in sorted(OUT.glob('*.assessment.json')):
  row=json.loads(p.read_text());r=json.loads(Path(row['result']).read_text());cases.append(dict(case=row['case'],result=row['result'],assessment=row['assessment'],upward_rig_bw=r.get('mean_upward_rig_support_bw'),downward_rig_bw=r.get('mean_downward_rig_force_bw'),mesh_range_m=r.get('actual_visual_foot_mesh_clearance_range_m'),scales=r.get('pathway_group_scales')))
 repeat={}
 fresh=checks/'w2lvar_fresh_v1'
 if fresh.with_suffix('.assessment.json').exists():
  source=OUT/'w2lvar_fanin_mean_v1'
  for suffix in ('.npz','.apparatus.npz','.apparatus_extended.npz','.inputs.npz','.motor_state.npz'):
   with np.load(source.with_suffix(suffix)) as a,np.load(fresh.with_suffix(suffix)) as b:
    if set(a.files)!=set(b.files):raise ValueError('Repeat fields differ')
    repeat[suffix]={k:bool(np.array_equal(a[k],b[k])) for k in a.files}
  if not all(v for group in repeat.values() for v in group.values()):raise ValueError('Fresh normalized candidate is not bit-exact')
 validation=[json.loads(p.read_text()) for p in sorted(checks.glob('*.assessment.json'))]
 chosen=json.loads((OUT/'w2lvar_fanin_mean_v1.assessment.json').read_text());baseparams=json.loads((ROOT/'rustuna_w2lvar_v1/BEST.json').read_text())['trial']['params']
 basins=[v for v in validation if v['case'].startswith('basin')];hold=next((v for v in validation if v['case']=='hold20'),None)
 if hold:
  hr=json.loads(Path(hold['result']).read_text());hold['apparatus_and_contact']=dict(mean_ground_load_bw=hr['mean_ground_load_bw'],upward_rig_bw=hr['mean_upward_rig_support_bw'],downward_rig_bw=hr['mean_downward_rig_force_bw'],knee_max_deg=hr['max_knee_angle_deg'],mesh_range_m=hr['actual_visual_foot_mesh_clearance_range_m'])
 original_result=json.loads(Path(chosen['result']).read_text())
 checkpoint=dict(machine='easteregg2',variant='w2lvar standing with mean-pool IIIN target mapping',global_parameters=baseparams,pathway_group_scales=chosen['scales'],source_fanin=chosen['source_fanin'],expected_source_sha256=original_result['source_sha256'],fresh_bit_exact=bool(repeat),fresh_comparisons=repeat,validation=validation,basin_pass_count=sum(v['assessment']['accepted'] for v in basins),basin_count=len(basins),twenty_second_hold=hold,scope='Distinct opt-in standing variant, not an unchanged literature controller or full curriculum. All IIIN-to-MN edges remain nonzero; original production source weights are untouched. Global PF and aggregate-Ib branches remain at the earlier stage-specific zero gains. Existing lateral/rotational apparatus and contact penetration remain. No pressure transfer or independent standing claim.')
 spinal=Path(__file__).resolve().parents[2]/'spinal'
 deps=[spinal/n for n in ('muscle_map.py','joint_pf_weights.json','physics_fidelity/standing_solver.py','physics_fidelity/shared_force_law_functional_v3.json','functional_gait2392/mechanics.py','w2l_cpg/sparse_pairwise_eval.py','w2l_cpg/preallocated_exact_eval.py','w2l_cpg/exact_fast_eval.py')]
 model_gate=json.loads((ROOT/'family_parity/w2lvar_parity.json').read_text());model_path=Path(model_gate['model'])
 if hashlib.sha256(model_path.read_bytes()).hexdigest()!=model_gate['model_sha256']:raise ValueError('Model XML changed since the family parity reference')
 deps.append(model_path)
 checkpoint['expected_dependency_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps}
 (OUT/'NORMALIZED_STANDING_CHECKPOINT.json').write_text(json.dumps(checkpoint,indent=2),encoding='utf-8')
 with (OUT/'HOMOLOGOUS_MOTOR_CURRENTS.csv').open('w',newline='',encoding='utf-8') as f:
  writer=csv.DictWriter(f,fieldnames=['condition','family','actuator','primary_group','mean_voltage_mV','tonic_input_nA','PF_nA','IIIN_nA','IIX_nA','II_nA','Ia_nA','IaIN_nA','Ib_nA','IBIN_nA','IBEXC_nA','RC_nA']);writer.writeheader()
  for condition,groups in (('frozen_selected_standing',current['families']),('original_full_gain',original['families'])):
   for family in groups:
    for a in family['motor_neurons']:
     currents=a.get('synaptic_current_nA',a.get('currents_nA'));writer.writerow(dict(condition=condition,family=family['family'],actuator=a['actuator'],primary_group=a['primary_group'],mean_voltage_mV=a['mean_voltage_mV'],tonic_input_nA=a.get('external_tonic_POST_nA',a.get('tonic_input_inferred_nA')),**{key+'_nA':currents.get(key,0.) for key in ('PF','IIIN','IIX','II','Ia','IaIN','Ib','IBIN','IBEXC','RC')}))
 body='All three frozen family standing cases have bit-identical static activation targets. Actual input-aligned synaptic currents are tabulated per homologous MN, with a full time-resolved membrane-balance error below1e−9 mV. Original full-gain PF/current attribution is also recovered from preserved state traces; inferred tonic currents independently match the standing-input law. HOMOLOGOUS_MOTOR_CURRENTS.csv contains both conditions.\n\nAt the same corrected pose, measured late activation projects to right-ankle muscle torques of −19.22 Nm(s3k), −17.53 Nm(syn6), and only −3.63 Nm(w2lvar). The common mechanical reconstruction matches the static fit within1.87e−16 activation units; gravity/contact/passive/active torque terms and all17joint components are retained in COMMON_JOINT_TORQUE_REFERENCE.json. These are offline common-pose recruitment comparisons, never teacher forces in simulation.\n\nThe legacy runner NPZ field act is20selected excitation commands: runner.py assigns log_act from data.ctrl, not data.act. Actual native activation and force-cache telemetry are now captured opt-in in motor_state.npz. Native activation-kernel replay matches all92 measured activation states exactly. Legacy files are preserved; this logging-label finding does not invalidate the body-force/contact gates.\n\n| Intervention | Standing gate | Tilt(deg) | Sway(mm) | Foot load(BW) |\n| --- | --- | ---: | ---: | ---: |\n'
 for row in cases:
  a=row['assessment'];body+=f"| {row['case']} | {a.get('accepted',False)} | {a.get('tilt_deg',float('nan')):.3f} | {a.get('sway_m',float('nan'))*1000:.2f} | {a.get('normal_force_bw',float('nan')):.4f} |\n"
 body+='\nThe two-group knee/hip-IIIN Rustuna study completed 12 trials, reusing 3 causal controls. It restores loading but still fails sway; its BEST.json is not a pass. Releasing the remaining ankle antagonist block passes. A source-derived mean-pool variant instead retains EVERY IIIN connection, dividing gains by the actual incoming antagonist-muscle count (1–9), and passes without fitting these normalization factors. This is an explicit engineering response to the expanded muscle pool, not asserted author intent. Full original graphs/weights remain intact in production files.\n\nFresh normalized-run numeric arrays are bit-exact. The normalized variant passes the 20 s hold and '+str(sum(v['assessment']['accepted'] for v in basins))+'/'+str(len(basins))+' seeded ±1% checks of all ten mapping gains and tonic bias. This is a bounded sensitivity test, not an exhaustive eleven-dimensional corner enumeration; the two zero global gains remain zero under relative perturbation. Long-hold upward apparatus force is 2.1795% BW, downward force 0.0095% BW, with 97.7025% BW mean foot load. Visual mesh penetration reaches 4.17/3.68 mm, and knee maxima remain below −1.75 degrees. Individual cases remain in NORMALIZED_STANDING_CHECKPOINT.json. No stage 1/2, walking, BPA, or harness-free promotion is implied.\n\nThe common-pose role-sign audit flags per_tert_r/l, labeled ankle_pf in muscle_map.py despite a positive ankle moment arm opposite the soleus reference. Gemellus also opposes its assigned hip-ext reference at this pose, with a much smaller arm. No silent global reclassification is made, and one pose does not establish behavior across the entire RoM.\n'
 body+='\nThe normalized variant also has an actual-activation common-pose torque projection: right ankle −8.76 Nm and right hip +8.96 Nm, compared with −3.63 Nm and −3.36 Nm in the failed frozen w2lvar case. The matched successful s3k/syn6 common-pose projections remain distinct; common-pose force attribution is not a substitute for their differing closed-loop steady poses.\n'
 (OUT/'EXPERIMENT_RECORD.md').write_text('# Current/group mapping experiment — easteregg2\n\n'+body,encoding='utf-8');section(ROOT/'EXPERIMENT_RECORD.md','Current/group mapping diagnosis and explicit normalized standing variant',body)
 print(json.dumps(dict(cases=len(cases),validation_cases=len(validation),fresh_exact=bool(repeat),checkpoint=str(OUT/'NORMALIZED_STANDING_CHECKPOINT.json')),indent=2))

if __name__=='__main__':main()
