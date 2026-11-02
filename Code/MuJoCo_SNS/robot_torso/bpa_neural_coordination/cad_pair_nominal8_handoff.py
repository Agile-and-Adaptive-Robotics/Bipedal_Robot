"""Report/fresh-check the frozen nominal8 CAD pair without changing its study."""
import argparse,json,hashlib,shutil
from pathlib import Path
import numpy as np
import cad_pair_nominal8 as N

ROOT=N.OUT
def load(p):return json.loads(p.read_text())

def compare():
 result={}
 for mode in ('reference','mixed'):
  source=ROOT/mode;repeat=ROOT/'fresh'/mode
  a=load(source/'ASSESSMENT.json');b=load(repeat/'ASSESSMENT.json');details=[]
  manifest=load(source/'TELEMETRY_MANIFEST.json')
  for entry in manifest['archives']:
   if not entry['full_trace_retained']:continue
   path=Path(entry['original_path']);other=repeat/path.relative_to(source)
   with np.load(path) as left,np.load(other) as right:
    assert set(left.files)==set(right.files)
    for key in left.files:
     x,y=left[key],right[key]
     if x.dtype.kind not in 'biufc':continue
     same=x.dtype==y.dtype and x.shape==y.shape and x.tobytes()==y.tobytes()
     details.append(dict(archive=str(path.relative_to(source)),array=key,exact=bool(same)))
  with np.load(source/'REPORT_CURVES.npz') as left,np.load(repeat/'REPORT_CURVES.npz') as right:
   for key in left.files:details.append(dict(archive='REPORT_CURVES.npz',array=key,exact=left[key].dtype==right[key].dtype and left[key].shape==right[key].shape and left[key].tobytes()==right[key].tobytes()))
  result[mode]=dict(original_pass=a['accepted'],fresh_pass=b['accepted'],assessment_exact=a==b,arrays=details,all_numeric_arrays_byte_exact=all(row['exact'] for row in details))
 result['passed']=all(v['original_pass'] and v['fresh_pass'] and v['assessment_exact'] and v['all_numeric_arrays_byte_exact'] for k,v in result.items() if k!='passed')
 (ROOT/'FRESH_REPEAT_AUDIT.json').write_text(json.dumps(result,indent=2));return result

def report():
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 colors=['#0077BB','#EE7733'];fig,axes=plt.subplots(5,2,figsize=(12,12),sharex=True,constrained_layout=True)
 results={}
 for col,mode in enumerate(('reference','mixed')):
  folder=ROOT/mode;score=load(folder/'ASSESSMENT.json');results[mode]=score
  with np.load(folder/'REPORT_CURVES.npz') as z:d={k:z[k] for k in z.files}
  t=d['t'];axes[0,col].set_title(('Ideal reference tensions' if mode=='reference' else 'Two virtual BPAs + 52 reference routes')+'\n'+('steady-gate pass' if score['accepted'] else 'steady-gate failure'))
  for leg in range(2):
   label=('Right','Left')[leg];style=('-','--')[leg]
   axes[0,col].plot(t,d['pair_command'][:,leg],style,color=colors[leg],label=label)
   axes[1,col].plot(t,d['pair_force_N'][:,leg],style,color=colors[leg])
   axes[2,col].plot(t,d['knee_deg'][:,leg],style,color=colors[leg])
  axes[0,col].legend(frameon=False,ncol=2);axes[0,col].set_ylabel('Direct pair MN output (0–1)');axes[1,col].set_ylabel('Pair tension (N)');axes[2,col].set_ylabel('Knee angle (°)')
  for j in range(3):axes[3,col].plot(t,d['rig_force_N'][:,j]/735.75,('-', '--', ':')[j],color=('#009988','#CC3311','#AA3377')[j],label=('Anterior','Superior','Right')[j])
  axes[3,col].set_ylabel('Signed harness force / BW');axes[3,col].legend(frameon=False,ncol=3,fontsize=8)
  axes[4,col].plot(t,d['pair_Qnorm_Nm'],color=colors[0],label='Pair');axes[4,col].plot(t,d['other52_Qnorm_Nm'],'--',color=colors[1],label='Other 52');axes[4,col].set_ylabel('Internal generalized-force norm (N m)');axes[4,col].set_xlabel('Time (s)');axes[4,col].legend(frameon=False,ncol=2)
  for ax in axes[:,col]:ax.axvspan(6,8,color='0.85',alpha=.45);ax.grid(alpha=.2);ax.set_xlim(0,8)
 fig.suptitle('CAD model: frozen free-SNS vastus-medialis pair, strongly supported nominal standing\nAll cells graded; other 52 routes retain the external hybrid allocator',fontsize=12)
 for suffix in ('png','pdf','svg'):fig.savefig(ROOT/('CAD_PAIR_NOMINAL8.'+suffix),dpi=200)
 plt.close(fig)
 # Pressure, COM drift, and compensation shown explicitly, without interpreting norms as net joint torques.
 fig,axes=plt.subplots(3,2,figsize=(12,8),sharex=True,constrained_layout=True)
 for col,mode in enumerate(('reference','mixed')):
  with np.load(ROOT/mode/'REPORT_CURVES.npz') as z:d={k:z[k] for k in z.files}
  t=d['t'];axes[0,col].set_title(mode)
  for leg in range(2):
   axes[0,col].plot(t,d['pressure_kpa'][:,leg],('-','--')[leg],color=colors[leg],label=('Right','Left')[leg])
  for j in range(3):axes[1,col].plot(t,1000*(d['COM'][:,j]-d['COM'][0,j]),('-', '--', ':')[j],label=('Anterior','Superior','Right')[j])
  axes[2,col].plot(t,d['other52_activation_sum'],color=colors[0],label='Sum of other 52 outputs');axes[2,col].plot(t,d['other52_floor_fraction']*52,'--',color=colors[1],label='Equivalent count at 5% floor')
  for row,label in enumerate(('Actual pair pressure (kPa)','Whole-model COM displacement (mm)','Other 52 output summary')):axes[row,col].set_ylabel(label);axes[row,col].legend(frameon=False,fontsize=8);axes[row,col].grid(alpha=.2);axes[row,col].axvspan(6,8,color='0.85',alpha=.45);axes[row,col].set_xlim(0,8)
  axes[2,col].set_xlabel('Time (s)')
 fig.suptitle('Nominal8 pressure, body drift, and external-controller context',fontsize=12)
 for suffix in ('png','pdf'):fig.savefig(ROOT/('CAD_PAIR_NOMINAL8_CONTEXT.'+suffix),dpi=200)
 plt.close(fig)
 fresh=load(ROOT/'FRESH_REPEAT_AUDIT.json') if (ROOT/'FRESH_REPEAT_AUDIT.json').exists() else {'passed':False,'status':'Not completed'}
 (ROOT/'RESULTS.json').write_text(json.dumps(dict(cases=results,fresh_repeat=fresh,joint_nominal8_accepted=all(v['accepted'] for v in results.values()) and fresh.get('passed',False)),indent=2))
 lines=['# Frozen CAD pair: nominal8 extension','','Machine: easteregg2. The six mirrored PF/IBEXC→MN weights, virtual geometry, pressure lag, source physics, and strong rig are unchanged from the accepted 2 s pair. Gray plot shading marks the preregistered 6–8 s settling window.','','Only the vas_med right/left pair uses literal, unclamped SNS recruitment, with zero external POST current and no voltage/activation floor. The remaining 52 routes retain the original hybrid torque allocator and 5% output floor. The mixed case contains two assumed pressure-actuated BPAs, not 54 calibrated hardware BPAs. All neurons are graded.','','| Case | Completed (s) | All gates | Upward harness (% BW) | Mean foot load (% BW) | Max knee (°), R/L | Peak pressure (kPa) |','| --- | ---: | --- | ---: | ---: | --- | ---: |']
 for mode,row in results.items():
  m=row['metrics'];lines.append(f"| {mode} | {m['duration_s']:.3f} | {row['accepted']} | {100*m['mean_upward_bw']:.3f} | {100*m['mean_ground_load_bw']:.3f} | {np.rad2deg(m['knee_max_rad']).round(3).tolist()} | {m['pressure_peak_kpa']:.3f} |")
 lines+=['',f"Fresh-process comparison: {fresh.get('passed',False)}. Byte equality includes every numeric full archive and every selected reporting curve; preparation/scheduling differences are not hidden.",'','IBEXC is the existing source-neuron name. Its dominant input in this condition is additive RG excitation; it is not a multiplicative load gate. The six calibrated outgoing PF/IBEXC connections preserve their original ratios and bilateral equality.','','## Stationarity and attribution','']
 for mode,row in results.items():
  m=row['metrics'];lines += [f"- {mode}: failed gates {[k for k,v in row['gates'].items() if not v]}; tail COM slopes (m/s) {m.get('tail_COM_slope_xyz_m_s')}; knee slopes (°/s) {m.get('tail_knee_slope_deg_s')}; knee ranges (°) {m.get('tail_knee_ptp_deg')}; internal velocity RMS maximum {m.get('tail_max_internal_velocity_rms_rad_s')} rad/s; pair mean tension {m.get('tail_pair_mean_force_N')} N; pair/other52 internal generalized-force norms {m.get('tail_pair_Qnorm_mean_Nm')}/{m.get('tail_other52_Qnorm_mean_Nm')} N m. These norms do not establish torque cancellation or causal independence."]
 lines+=['','The 8 s duration and tail gates establish a bounded nominal settling test only. They do not establish endurance, robustness, harness-free balance, locomotion, or whole-body neural autonomy. No harness reduction or postfailure tuning was performed. Signed harness forces/torques and actual moving-route forces were recorded. Minimum visual mesh height is a penetration diagnostic, not evidence that the meshes define the collision geometry.','','## IO recovery record','','Two physically completed reference runs lost RAM histories due to readout errors (wrong core time key, then NumPy Boolean JSON serialization). Their original summaries/source/protocol versions remain in reference_readout_error and reference_serialization_error. They are not numerical failures or independent scientific repetitions. The reader was tested on explicitly synthetic tiled-data fixtures in repo tmp, which are not experimental results. Automatic error-time RAM archive persistence now prevents further readout errors from requiring physics repeats. The final reference rerun recovers the missing stationarity data without changing physical parameters or gate thresholds.','','## Files and reproduction','',f'- Protocol: `{ROOT / "PREREGISTRATION.json"}`',f'- Cases: `{ROOT / "reference"}` and `{ROOT / "mixed"}`',f'- Full-time selected curves: each case `REPORT_CURVES.npz`; full archives retained only for gate passers.',f'- Program: `{Path(N.__file__)}`',f'- Plot/report and fresh-case program: `{Path(__file__)}`','- Run a fresh named comparison with this report program `--fresh-case reference`, then in another process `--fresh-case mixed`. Existing completed assessments are reused, not overwritten.','- Build the comparison and plots with `--report`. Source model/controller originals and accepted 2 s figures remain unchanged.']
 (ROOT/'EXPERIMENT_RECORD.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 (ROOT/'CAD_PAIR_NOMINAL8.alt.txt').write_text('Two-column eight-second comparison of a free graded SNS vastus-medialis pair using ideal reference tension versus virtual BPA pressures. Rows show pair neural output, actual pair tension, both knee angles, signed harness forces, and internal generalized-force norms for the pair and externally controlled other52 routes. The final two seconds are shaded for preregistered settling checks. The strong harness and other52 external hybrid controller remain active.',encoding='utf-8')
 sources=ROOT/'source_snapshot';sources.mkdir(exist_ok=True);ledger={}
 for file in (Path(N.__file__),Path(__file__)):
  digest=hashlib.sha256(file.read_bytes()).hexdigest();target=sources/file.name
  if target.exists():assert hashlib.sha256(target.read_bytes()).hexdigest()==digest,'Existing source snapshot differs'
  else:shutil.copy2(file,target)
  ledger[str(file)]=dict(sha256=digest,snapshot=str(target))
 (sources/'MANIFEST.json').write_text(json.dumps(dict(new_programs=ledger,frozen_scientific_dependencies=str(ROOT/'PREREGISTRATION.json'),prior_controller_snapshot=str(N.BASE/'source_snapshot/MANIFEST.json')),indent=2))
 return results

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--fresh-case',choices=['reference','mixed']);p.add_argument('--report',action='store_true');a=p.parse_args()
 if a.fresh_case:N.OUT=ROOT/'fresh';N.OUT.mkdir(exist_ok=True);N.run(a.fresh_case)
 if a.report:
  if all((ROOT/'fresh'/m/'ASSESSMENT.json').exists() for m in ('reference','mixed')):compare()
  report()
