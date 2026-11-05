"""Source108 closed-loop, four-parameter supported-ground calibration adapter."""
import hashlib,json,re
from pathlib import Path
import numpy as np
from import_virtual_walker import HERE
from run_virtual_walker import run
from assess_trial import assess
SOURCE=HERE.parents[3]/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer_Standalone.asim'
EXPECTED='01ece6edcbcf17c274410b23f1ce7603eab288aadb946d6b548c5102ed5c1154'
CONFIG=HERE/'CANONICAL108_NEARFLAT_ANKLE_PHASE.json'

ROOT=None

def handle_request(request):
 global ROOT
 op=request['op']
 if op=='preflight':
  ROOT=Path(request['manifest']['output_dir'])/'body_trials'
  baseline=HERE/'runs/W2L_VirtualWalker_Newer/canonical108_ground_nearflat_ankle_crossed_fresh_v1'
  a=json.loads((baseline/'ASSESSMENT.json').read_text());parity=json.loads((baseline/'FRESH_PARITY.json').read_text());source_ok=hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED
  ready=source_ok and a['alternating_supported_stepping_gate'] and parity['all_numeric_arrays_bit_exact']
  return dict(ok=True,ready_for_optimization=ready,source_sha256=EXPECTED,neurons=108,chemical_edges=216,baseline_gate=a['alternating_supported_stepping_gate'],fresh_byte_parity=parity['all_numeric_arrays_bit_exact'],reason='Actual108 near-flat supported stepping and source hash verified; not independent walking.',objective='Strict physical stepping gate first; then upward support + gentle p95/peak impact and relative-belt progress costs. Net downward forces cannot cancel lifting support in the objective.')
 if op!='evaluate':raise ValueError(op)
 params=request['params'];assert set(params)=={'pf_mn_recruitment_scale','toe_mn_gmax','hip_pf_mn_scale','neural_time_scale'}
 assert .8<=params['pf_mn_recruitment_scale']<=1.5 and 0<=params['toe_mn_gmax']<=.05
 assert .5<=params['hip_pf_mn_scale']<=1.3 and .85<=params['neural_time_scale']<=1.15
 assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED,'Canonical source changed mid-study'
 cfg=json.loads(CONFIG.read_text());tag='rustuna108_'+re.sub('[^A-Za-z0-9_.-]','_',str(request['evaluation_id']))
 r=run(newer=True,duration=10.,tag=tag,output_root=ROOT,**cfg,**params);folder=ROOT/'W2L_VirtualWalker_Newer'/tag;a=assess(folder);z=np.load(folder/'trajectory.npz');bw=r['robot_mass_kg']*9.81;mask=z['t']>=1.
 if not np.any(mask):mask=np.ones(len(z['t']),dtype=bool)
 support={}
 for label,key in [('rails','rail_force_world'),('guides','guide_force_world'),('harness','harness_force')]:
  f=z[key][mask,1];support[label]=dict(net_BW=float(np.mean(f)/bw),upward_BW=float(np.mean(np.maximum(f,0))/bw),downward_BW=float(np.mean(np.maximum(-f,0))/bw))
 upward=sum(x['upward_BW'] for x in support.values());foot=np.sum(z['foot_vertical'][mask],axis=1)/bw;p95=float(np.percentile(foot,95));peak=float(np.max(foot));speed=r['root_forward_relative_to_belt_m']/r['duration_completed_s'];knees=r['anatomical_knee_min_deg'];passed=a['alternating_supported_stepping_gate'] and r['failure'] is None and min(knees)>0 and r['root_y_min_m']>.90
 if passed:value=upward+.02*p95+.002*peak+.05*max(0.,(.03-speed)/.03)
 else:value=100.+10.*sum(max(0,3-a['legs'][s]['unloaded_clear_swings']) for s in ('L','R'))+20.*max(0,.8-a['alternating_transition_fraction'])+(10000. if r['failure'] is not None or min(knees)<=0 else 0.)
 assert np.isfinite(value)
 diagnostics=dict(folder=str(folder),physical_gate_passed=passed,source_sha256=EXPECTED,swings={s:a['legs'][s]['unloaded_clear_swings'] for s in ('L','R')},alternation=a['alternating_transition_fraction'],upward_assistance_BW=upward,net_assistance_BW=r['mean_net_vertical_assistance_BW'],support_components=support,foot_force_p95_BW=p95,foot_force_peak_BW=peak,relative_belt_progress_m=r['root_forward_relative_to_belt_m'],relative_belt_speed_m_s=speed,knee_min_deg=knees,source_condition='Actual108; near-flat physical home; explicit4edgePF→ankleMN phase swap; source30mm passiveguides/rails−20mm; runtimeclock1; unchangedIa/II/Ib/c1/V3/heelpaths; no allocator.')
 out=dict(ok=True,value=float(value),diagnostics=diagnostics);(folder/'OPTIMIZATION_EVALUATION.json').write_text(json.dumps(dict(request=request,result=out),indent=2));return out
