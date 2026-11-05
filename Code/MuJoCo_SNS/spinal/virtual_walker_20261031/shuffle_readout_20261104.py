"""Second shuffle search: explicit neural-output adapter calibration, same task gates."""
import json,hashlib,re,contextlib
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
from import_virtual_walker import HERE
from run_virtual_walker import run
from w2l_ground_optimization import SOURCE,EXPECTED
ROOT=None

def measure(folder):
 z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());t=z['t'];mask=t>=2.;dt=r['physics_dt'];bw=r['robot_mass_kg']*9.81;load=z['foot_vertical']/bw;loaded=load>.05;both=mask&loaded.all(axis=1)
 edges=np.diff(np.r_[False,both,False].astype(int));starts=np.flatnonzero(edges==1);ends=np.flatnonzero(edges==-1);dur=(ends-starts)*dt;sustained=float(np.sum(dur[dur>=.04])/(mask.sum()*dt));flight=float(np.mean(~loaded[mask].any(axis=1)));legs={};events=[]
 for si,s in enumerate(('L','R')):
  clear=np.min(z['foot_clearance'][:,2*si:2*si+2],axis=1);swing=mask&(clear>.001)&(load[:,si]<.025);e=np.diff(np.r_[False,swing,False].astype(int));aa=np.flatnonzero(e==1);bb=np.flatnonzero(e==-1);qualified=[]
  for begin,end in zip(aa,bb):
   if (end-begin)*dt>=.04 and np.max(clear[begin:end])>=.003 and np.mean(load[begin:end,1-si]>.1)>=.5:
    qualified.append(float(t[begin]));events.append((float(t[begin]),s))
  position=z['foot_body_position'][:,si,0]-z['root'][:,0];excursion=float(np.percentile(position[mask],95)-np.percentile(position[mask],5))
  legs[s]=dict(clear_unloaded_swings=len(qualified),event_times_s=qualified,period_s=float(np.median(np.diff(qualified)))if len(qualified)>1 else None,foreaft_excursion_m=excursion,peak_whole_foot_clearance_m=float(np.max(clear[mask])))
 events.sort();alternation=sum(a[1]!=b[1]for a,b in zip(events[:-1],events[1:]))/max(1,len(events)-1)
 upward=sum(float(np.mean(np.maximum(z[k][mask,1],0.))/bw)for k in('rail_force_world','guide_force_world','harness_force'));downward=sum(float(np.mean(np.maximum(-z[k][mask,1],0.))/bw)for k in('rail_force_world','guide_force_world','harness_force'))
 knees=r['anatomical_knee_min_deg'];complete=r['duration_completed_s']>=15.99;periods=[v['period_s']for v in legs.values()];ranges=[v['foreaft_excursion_m']for v in legs.values()]
 gates=dict(complete=complete,finite=r['failure']is None,knees=min(knees)>0,swings=min(v['clear_unloaded_swings']for v in legs.values())>=5,alternation=alternation>=.8,sustained_double_support=sustained>=.08,no_flight=flight<=.02,short_steps=max(ranges)<=.20,slow_cycles=all(p is not None and 1.2<=p<=2.8 for p in periods),assistance=upward<=.30,no_downward_bracing=downward<=.02)
 passed=all(gates.values());penalty=sum(not v for v in gates.values());value=(0.5*upward+sum(ranges)+abs(sustained-.16)+flight)if passed else 10*penalty+10*max(0,.08-sustained)+10*flight+sum(ranges)+upward
 return dict(value=float(value),passed=passed,gates=gates,legs=legs,alternation=float(alternation),double_support_fraction=float(np.mean(both[mask])),sustained_double_support_fraction=sustained,double_support_bouts_s=dur.tolist(),airborne_fraction=flight,upward_assistance_BW=upward,downward_assistance_BW=downward,knee_min_deg=knees,root_x_change_m=r['root_x_change_m'],measured_belt_x_velocity_m_s=r['belt_world_x_change_m']/max(r['duration_completed_s']-.5,.001),simulated_s_per_wall_s=r['simulated_s_per_wall_s'])

def handle_request(request):
 global ROOT
 if request['op']=='preflight':
  ROOT=Path(request['manifest']['output_dir'])/'body_trials';good=hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED
  return dict(ok=True,ready_for_optimization=good,source_sha256=EXPECTED,scope='Previously verified W2L108 port; baseline does not pass the new shuffle task. Task gates predeclared in adapter, not relaxed after search.')
 if request['op']!='evaluate':raise ValueError(request['op'])
 params=request['params'];assert set(params)=={'belt_speed_m_s','hip_excitation_scale','knee_ext_excitation_scale','neural_time_scale'}
 assert .02<=params['belt_speed_m_s']<=.06 and .5<=params['hip_excitation_scale']<=1. and .5<=params['knee_ext_excitation_scale']<=1. and 1.8<=params['neural_time_scale']<=2.4
 cfg=json.loads((HERE/'CANONICAL108_NEARFLAT_ANKLE_PHASE.json').read_text());cfg.update(pf_mn_recruitment_scale=1.1152227449218066,toe_mn_gmax=.04884415411663129,belt_motor_end_s=16.,hip_pf_mn_scale=1.0375630740749227,symmetric_rg_ext_tonic=.09891086315354553,**params);tag=re.sub('[^A-Za-z0-9_.-]','_',str(request['evaluation_id']))
 run(newer=True,duration=16.,tag=tag,output_root=ROOT,**cfg);folder=ROOT/'W2L_VirtualWalker_Newer'/tag;result=measure(folder);result['folder']=str(folder);result['configuration']=cfg;(folder/'SHUFFLE_ASSESSMENT.json').write_text(json.dumps(result,indent=2))
 return dict(ok=True,value=result['value'],diagnostics=result)
