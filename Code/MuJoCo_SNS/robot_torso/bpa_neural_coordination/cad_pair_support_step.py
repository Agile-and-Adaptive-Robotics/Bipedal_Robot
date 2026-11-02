"""One preregistered 0.90→0.85BW feedforward-support step; no height change."""
import argparse,json,hashlib,copy,shutil
from pathlib import Path
import numpy as np
import cad_pair_nominal8 as N

BASE=N.OUT;ROOT=BASE.parent/'support085';BW=735.75;SUPPORT=.85
def read(p):return json.loads(p.read_text())
def arrays(p):
 with np.load(p) as z:return {k:z[k] for k in z.files}
def law(a,support):
 q=a['qpos_full'];v=a['qvel_full'];spring=800*(q[0,1]-q[:,1]);damping=-120*v[:,1];prediction=np.clip(support*BW+spring+damping,0,1.2*BW)
 return dict(feedforward_N=support*BW,spring=spring,damping=damping,predicted=prediction,max_error_N=float(np.max(np.abs(prediction-a['rig'][:,N.C.mujoco.MjModel.from_xml_path(str(N.C.ROBOT/'ah001_neural.xml')).body('pelvis').id,1]))))

def baseline_audit():
 ROOT.mkdir(exist_ok=True);data={}
 for mode in ('reference','mixed'):
  a=arrays(BASE/mode/mode/'ACTUATOR_BRIDGE.npz');proof=law(a,.9);data[mode]=dict(max_law_error_N=proof['max_error_N'],initial_root_y_m=float(a['qpos_full'][0,1]),initial_qvel_max=float(np.max(np.abs(a['qvel_full'][0]))),mean_spring_N=float(proof['spring'].mean()),mean_damping_N=float(proof['damping'].mean()),initial_qpos_sha256=hashlib.sha256(a['qpos_full'][0].tobytes()).hexdigest())
 result=dict(machine='easteregg2',mass_kg=75.,gravity_m_s2=9.81,source=str(N.C.R.__file__),source_law='clip(support*BW + 800*(initial_root_y-current_root_y) - 120*root_vy, 0, 1.2*BW)',passed=all(row['max_law_error_N']<1e-9 for row in data.values()),cases=data)
 (ROOT/'SUPPORT_LAW_AUDIT.json').write_text(json.dumps(result,indent=2));return result

def protocol():
 root=N.OUT;root.mkdir(exist_ok=True);base=read(BASE/'PREREGISTRATION.json')
 assert read(BASE/'FRESH_REPEAT_AUDIT.json')['passed']
 for path,digest in base['source_sha256'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest,'Frozen dependency changed:'+path
 proof=read(ROOT/'SUPPORT_LAW_AUDIT.json');assert proof['passed']
 p=copy.deepcopy(base);p.update(experiment='Single vertical feedforward support reduction, not harness-height change',parent_accepted_nominal8=str(BASE),nominal_support_before_bw=.9,nominal_support_after_bw=SUPPORT,nominal_feedforward_reduction_N=.05*BW,initial_root_y_unchanged=True,vertical_stiffness_N_m=800.,vertical_damping_Ns_m=120.,vertical_force_cap_bw=1.2,other52_algorithm_unchanged_but_measured_load_response_allowed=True,apparatus='Feedforward85%BW + unchangedheight spring800N/m and damping120Ns/m, cap120%BW; lateral250/damping60 and angular300/damping60 frozen',actual_support_claim='Measured comparison; spring response can offset some of the nominal5percentagepoint reduction',support_law_audit=str(ROOT/'SUPPORT_LAW_AUDIT.json'),same_gates_as_parent=True,no_posture_or_controller_adjustment_after_failure=True,only_one_support_step=True,additional_audits=['actual consumed CLI support=.85','same initial qpos and qvel as parent','reconstructed measuredverticalforce law'],driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 path=root/'PREREGISTRATION.json'
 if path.exists():assert read(path)==p,'Registered support study changed'
 else:path.write_text(json.dumps(p,indent=2))
 return p

def run(mode,fresh=False):
 if not (ROOT/'SUPPORT_LAW_AUDIT.json').exists():assert baseline_audit()['passed']
 N.OUT=ROOT/('fresh' if fresh else '');N.protocol=protocol
 old=N.L.compile_opt_in
 def compile_support():
  fn,parser,code,edits=old();parse=parser.parse_args
  def reduced(args=None,*a,**kw):
   args=list(args);i=args.index('--support');assert float(args[i+1])==.9;args[i+1]=str(SUPPORT);return parse(args,*a,**kw)
  parser.parse_args=reduced;return fn,parser,code,edits
 N.L.compile_opt_in=compile_support
 score=N.run(mode);case=N.OUT/mode;folder=case/mode
 # Correct only inherited experiment descriptions, preserving exact legacytext.
 corrections=[]
 for path,updates in [(case/'RECRUITMENT_PREREGISTRATION.json',dict(nominal_support_bw=SUPPORT)),(folder/'LIVE_PREREGISTRATION.json',dict(same_controller='Two literal freeSNS MNs;other52explicitexternalhybridcontext')),(folder/'cad.json',dict(controller='Layered gradedSNS; twofreeMNs, other52externalhybrid allocator/currentclamp',controller_scope='Literal SNS recruitment of vas_med bilateralpair;zero POSTclamp/no pairfloor. Other52unchangedhybridcontext.')),(folder/'BRIDGE_SUMMARY.json',dict(scope='Two freeSNS pairroutes with'+('virtualpressureBPAs' if mode=='mixed' else 'idealreferenceforces')+';other52hybrid; single85%BWfeedforwardstep;notautonomouswholebodycontrol'))]:
  d=read(path);prior={k:d.get(k) for k in updates};d['legacy_support_step_descriptions']=prior;d.update(updates);path.write_text(json.dumps(d,indent=2));corrections.append(dict(path=str(path),old=prior,new=updates,numerical_results_changed=False))
 (case/'METADATA_CORRECTIONS.json').write_text(json.dumps(corrections,indent=2))
 # Audit after physics without changing the predeclared safety/settling decisions.
 traces=arrays(case/'REPORT_CURVES.npz');base=arrays(BASE/mode/'REPORT_CURVES.npz');params=read(folder/'LIVE_PREREGISTRATION.json')['arguments'];audit=dict(consumed_support=float(params['support']),unchanged_initial_root_y=bool(traces['root_position'][0,1]==base['root_position'][0,1]),same_initial_root_position=bool(np.array_equal(traces['root_position'][0],base['root_position'][0])),measured_mean_upward_change_percentage_points=100*(score['metrics']['mean_upward_bw']-read(BASE/mode/'ASSESSMENT.json')['metrics']['mean_upward_bw']),COM_final_displacement_mm=(1000*(traces['COM'][-1]-traces['COM'][0])).tolist(),COM_min_anterior_displacement_mm=float(1000*np.min(traces['COM'][:,0]-traces['COM'][0,0])),mean_other52_output_sum=float(traces['other52_activation_sum'].mean()),tail_other52_output_sum=float(traces['other52_activation_sum'][traces['t']>=6].mean()) if np.any(traces['t']>=6) else None)
 if (folder/'ACTUATOR_BRIDGE.npz').exists():
  data=arrays(folder/'ACTUATOR_BRIDGE.npz');prior=arrays(BASE/mode/mode/'ACTUATOR_BRIDGE.npz');reconstruction=law(data,SUPPORT);audit.update(vertical_force_law_max_error_N=reconstruction['max_error_N'],initial_full_qpos_exact=bool(np.array_equal(data['qpos_full'][0],prior['qpos_full'][0])),initial_full_qvel_exact=bool(np.array_equal(data['qvel_full'][0],prior['qvel_full'][0])),mean_vertical_spring_N=float(reconstruction['spring'].mean()),mean_vertical_damping_N=float(reconstruction['damping'].mean()))
 audit['consumed_support_correct']=audit['consumed_support']==SUPPORT
 (case/'SUPPORT_COMPARISON.json').write_text(json.dumps(audit,indent=2));return score

def report():
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
 outcomes={};fig,axes=plt.subplots(4,2,figsize=(12,10),sharex=True,constrained_layout=True)
 for col,mode in enumerate(('reference','mixed')):
  score=read(ROOT/mode/'ASSESSMENT.json');audit=read(ROOT/mode/'SUPPORT_COMPARISON.json');outcomes[mode]=dict(assessment=score,audit=audit)
  d=arrays(ROOT/mode/'REPORT_CURVES.npz');b=arrays(BASE/mode/'REPORT_CURVES.npz');t=d['t'];axes[0,col].set_title(mode+': '+('settling pass' if score['accepted'] else 'failed gates'))
  for leg in range(2):axes[0,col].plot(t,d['knee_deg'][:,leg],('-','--')[leg],color=('#0077BB','#EE7733')[leg],label=('Right','Left')[leg])
  axes[1,col].plot(t,d['rig_force_N'][:,1]/BW,label='85% feedforward');axes[1,col].plot(b['t'],b['rig_force_N'][:,1]/BW,'--',label='90% parent')
  for j in range(3):axes[2,col].plot(t,1000*(d['COM'][:,j]-d['COM'][0,j]),('-','--',':')[j],label=('Anterior','Superior','Right')[j])
  axes[3,col].plot(t,d['pair_Qnorm_Nm'],label='Pair');axes[3,col].plot(t,d['other52_Qnorm_Nm'],'--',label='Other52')
  for row,label in enumerate(('Knee angle (°)','Actual vertical harness force / BW','Whole-model COM displacement (mm)','Internal generalized-force norm (N m)')):
   axes[row,col].set_ylabel(label);axes[row,col].grid(alpha=.2);axes[row,col].legend(frameon=False,fontsize=8);axes[row,col].axvspan(6,8,color='.85',alpha=.4);axes[row,col].set_xlim(0,8)
  axes[3,col].set_xlabel('Time (s)')
 fig.suptitle('One frozen CAD support step: 90% → 85% BW feedforward\nHeight target, vertical spring/damping, lateral/rotation guides, six SNS weights, and other52 hybrid algorithm unchanged',fontsize=12)
 for ext in ('png','pdf','svg'):fig.savefig(ROOT/('SUPPORT_STEP.'+ext),dpi=200)
 plt.close(fig)
 fresh={}
 for mode in ('reference','mixed'):
  path=ROOT/'fresh'/mode/'ASSESSMENT.json'
  if not path.exists():continue
  a=read(ROOT/mode/'ASSESSMENT.json');b=read(path);rows=[]
  for item in read(ROOT/mode/'TELEMETRY_MANIFEST.json')['archives']:
   if not item['full_trace_retained']:continue
   source=Path(item['original_path']);target=ROOT/'fresh'/mode/source.relative_to(ROOT/mode);x=arrays(source);y=arrays(target)
   for key in x:
    if x[key].dtype.kind in 'biufc':rows.append(dict(archive=str(source.relative_to(ROOT/mode)),array=key,exact=x[key].shape==y[key].shape and x[key].dtype==y[key].dtype and x[key].tobytes()==y[key].tobytes()))
  fresh[mode]=dict(assessment_exact=a==b,accepted=b['accepted'],all_numeric_arrays_byte_exact=all(r['exact'] for r in rows),arrays=rows)
 (ROOT/'RESULTS.json').write_text(json.dumps(dict(outcomes=outcomes,fresh=fresh),indent=2))
 if fresh:(ROOT/'FRESH_REPEAT_AUDIT.json').write_text(json.dumps(fresh,indent=2))
 lines=['# Single reduced vertical-assistance test','','Machine: easteregg2. This changes only the consumed `--support` argument from 0.90 to 0.85. The 75 kg model therefore loses 36.7875 N of prescribed vertical feedforward force. Actual force also includes the unchanged 800 N/m height spring and 120 Ns/m damper, clipped to [0,1.2BW]. Initial height and the height target are unchanged; actual measured support is reported rather than assumed to drop exactly five percentage points.','','Six mirrored SNS recruitment weights, assumed BPA geometry/pressure dynamics, other52hybrid algorithm, lateral/rotation guides, sourcepose and initialization remain frozen. The other52algorithm uses the measured unsupported load, so its resulting recruitment can change. All parent8s gates are unchanged. No postfailure tuning.','','| Case | Completed s | Gates pass | Actual upward support %BW | Change, percentagepoints | Footload %BW | Peak pressure kPa |','| --- | ---: | --- | ---: | ---: | ---: | ---: |']
 for mode,row in outcomes.items():
  a=row['assessment'];m=a['metrics'];v=row['audit'];lines.append(f"| {mode} | {m['duration_s']:.3f} | {a['accepted']} | {100*m['mean_upward_bw']:.3f} | {v['measured_mean_upward_change_percentage_points']:.3f} | {100*m['mean_ground_load_bw']:.3f} | {m['pressure_peak_kpa']:.3f} |")
  lines.extend(['',f"{mode}: failed gates {[k for k,v in a['gates'].items() if not v]}; final COM displacement xyz {v['COM_final_displacement_mm']} mm, peak posterior displacement {v['COM_min_anterior_displacement_mm']} mm. Tailother52outputsum {v['tail_other52_output_sum']}; pair/other52internalforcenorms {m.get('tail_pair_Qnorm_mean_Nm')}/{m.get('tail_other52_Qnorm_mean_Nm')} Nm. Actualforce-law and initialization audit: {v}.",''])
 lines+=['','Fresh checks: '+str({k:{name:value for name,value in row.items() if name!='arrays'} for k,row in fresh.items()}),'','Scope: two literal graded SNS-controlled routes, only those two pressureactuated in mixedmode, surrounded by 52hybrid idealreference routes and strongguidance. No pair-autonomousbalance, hardwarecalibration, locomotion, robustness, or endurance claim. Prescribedforce is not adjusted by secretly moving the harness height.','','Protocol: `'+str(ROOT/'PREREGISTRATION.json')+'`. Baseline law audit: `'+str(ROOT/'SUPPORT_LAW_AUDIT.json')+'`. Per-case comparison and rawcurves are saved under reference/ and mixed/. All newly inherited controller-only descriptions are corrected with originaltext retained in per-case METADATA_CORRECTIONS.json; numericresults are not changed.']
 (ROOT/'EXPERIMENT_RECORD.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');snapshot=ROOT/'source_snapshot';snapshot.mkdir(exist_ok=True);shutil.copy2(__file__,snapshot/Path(__file__).name)
 print(json.dumps({k:dict(pass_=v['assessment']['accepted'],support=v['assessment']['metrics']['mean_upward_bw'],fresh=fresh.get(k,{}).get('all_numeric_arrays_byte_exact')) for k,v in outcomes.items()},indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case',choices=['reference','mixed']);p.add_argument('--fresh',action='store_true');p.add_argument('--audit',action='store_true');p.add_argument('--report',action='store_true');a=p.parse_args()
 if a.audit:print(json.dumps(baseline_audit(),indent=2));N.OUT=ROOT;print(json.dumps(protocol(),indent=2))
 if a.case:run(a.case,a.fresh)
 if a.report:report()
