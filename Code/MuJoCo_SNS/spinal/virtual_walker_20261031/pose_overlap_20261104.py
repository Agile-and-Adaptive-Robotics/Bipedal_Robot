"""Bounded paired pelvis-pose/stance-overlap diagnostics, with actual foot loads."""
import argparse,json,contextlib
from pathlib import Path
import numpy as np
from import_virtual_walker import HERE,mujoco
from run_virtual_walker import run
from assess_trial import assess
from render_trial import render

def metrics(folder):
 z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());a=assess(folder);t=z['t'];mask=t>=1.;bw=r['robot_mass_kg']*9.81;loaded=z['foot_vertical']>.05*bw;both=loaded.all(axis=1)&mask
 edges=np.diff(np.r_[False,both,False].astype(int));starts=np.flatnonzero(edges==1);ends=np.flatnonzero(edges==-1);dur=(ends-starts)*r['physics_dt']
 # Native free-root quaternion rotates local superior Y into the world frame.
 q=z['q'][:,3:7];w,x,y,zz=q.T;up_x=2*(x*y-w*zz);up_y=1-2*(x*x+zz*zz);tilt=np.degrees(np.arctan2(up_x,up_y))
 swing_angles=[]
 for si,side in enumerate(('L','R')):
  for event in a['legs'][side]['swing_intervals']:
   keep=(t>=event['start_s'])&(t<event['end_s']);swing_angles.extend(z['angles'][keep,1+3*si].tolist())
 upward=sum(np.mean(np.maximum(z[key][mask,1],0.))/bw for key in ('rail_force_world','guide_force_world','harness_force'))
 return dict(passed=bool(a['alternating_supported_stepping_gate']and min(r['anatomical_knee_min_deg'])>0),double_support_fraction=float(np.mean(both[mask])),double_support_bouts_s=dur.tolist(),mean_double_support_bout_s=float(np.mean(dur))if len(dur)else 0.,airborne_fraction=float(np.mean(~loaded[mask].any(axis=1))),mean_anterior_pelvis_tilt_deg=float(np.mean(tilt[mask])),pelvis_tilt_range_deg=[float(tilt[mask].min()),float(tilt[mask].max())],mean_swing_knee_flexion_deg=float(np.mean(swing_angles))if swing_angles else 0.,straight_knee_swing_fraction=float(np.mean(np.array(swing_angles)<10.))if swing_angles else 1.,upward_assistance_BW=float(upward),foot_force_peak_BW=float(np.max(z['foot_vertical'].sum(axis=1))/bw),swings={s:a['legs'][s]['unloaded_clear_swings']for s in ('L','R')},alternation=a['alternating_transition_fraction'],minimum_knee_flexion_deg=r['anatomical_knee_min_deg'],failure=r['failure'])

def main():
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--feedback-ladder',action='store_true');p.add_argument('--posture-ladder',action='store_true');p.add_argument('--posture-refine',action='store_true');p.add_argument('--posture-damping',action='store_true');p.add_argument('--rail-ladder',action='store_true');p.add_argument('--repeat-case');a=p.parse_args()
 if a.repeat_case:
  saved=json.loads((a.root/'COMPLETE.json').read_text());selected=next(r for r in saved['cases']if r['case']==a.repeat_case)
  with(a.root/'fresh_repeat.log').open('w',encoding='utf-8')as log,contextlib.redirect_stdout(log):
   run(newer=True,duration=10.,tag='repeat',output_root=a.root/'fresh_repeat',**selected['configuration']);folder=a.root/'fresh_repeat/W2L_VirtualWalker_Newer/repeat';measure=metrics(folder);old=np.load(Path(selected['folder'])/'trajectory.npz');fresh=np.load(folder/'trajectory.npz');parity={key:old[key].shape==fresh[key].shape and old[key].tobytes()==fresh[key].tobytes()for key in old.files if old[key].dtype.kind in 'fiu'}
   (a.root/'FRESH_REPEAT.json').write_text(json.dumps(dict(case=selected['case'],all_numeric_arrays_byte_exact=all(parity.values()),arrays=parity,metrics=measure,folder=str(folder)),indent=2))
   if not(all(parity.values())and measure['passed']):raise RuntimeError('Fresh-repeat gate failed')
   render(folder,5.,8.2,title='W2L108 | stance-gated hip PF posture',detail_caption=f"Double support {100*measure['double_support_fraction']:.1f}% | Upward assistance {100*measure['upward_assistance_BW']:.1f}% BW (after 1 s)")
  print(json.dumps(dict(repeated=all(parity.values()),metrics=measure)));return
 a.root.mkdir(parents=True,exist_ok=False)
 parent=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\LiDeng_W2L_20261104');best=json.loads((parent/'SUMMARY.json').read_text())['best'];cfg=json.loads((HERE/'CANONICAL108_NEARFLAT_ANKLE_PHASE.json').read_text());cfg.update(best['params']);m=mujoco.MjModel.from_xml_path(str(parent/'fresh_repeat/W2L_VirtualWalker_Newer/repeat/model.xml'))
 cases=[('baseline',0.,None,1.),('tilt2',2.,None,1.),('tilt4',4.,None,1.),('ext005',0.,.05,1.),('ext01',0.,.1,1.),('ext02',0.,.2,1.),('tilt2_ext01',2.,.1,1.),('tilt4_ext01',4.,.1,1.),('hipflex08',0.,None,.8),('tilt2_hipflex08',2.,None,.8)]
 cases=[(*case,{})for case in cases]
 if a.feedback_ladder:
  cases=[('baseline',0.,None,1.,{})]+[(label,0.,None,1.,extra)for label,extra in [('pf125',dict(contact_pf_scale=1.25)),('pf15',dict(contact_pf_scale=1.5)),('pf2',dict(contact_pf_scale=2.)),('v3half',dict(v3_output_scale=.5)),('v3up',dict(v3_output_scale=1.25)),('knee12',dict(pf_knee_flex_scale=1.2)),('knee14',dict(pf_knee_flex_scale=1.4)),('slower',dict(neural_time_scale=1.3))]]
 if a.posture_ladder:
  cases=[('baseline',0.,None,1.,{})]+[(f'posture_{target}_{kp}',0.,None,1.,dict(pelvis_pf_posture=dict(target_deg=target,kp=kp,kd=.15,limit_nA=.6)))for target in(0.,4.)for kp in(1.,2.,4.)]
 if a.posture_refine:
  cases=[('baseline',0.,None,1.,{})]+[(f'posture_refine_{kp}',0.,None,1.,dict(pelvis_pf_posture=dict(target_deg=0.,kp=kp,kd=.15,limit_nA=.6)))for kp in(1.2,1.4,1.6,1.8)]
 if a.posture_damping:
  cases=[('baseline',0.,None,1.,{})]+[(f'pd_{kp}_{kd}',0.,None,1.,dict(pelvis_pf_posture=dict(target_deg=2.,kp=kp,kd=kd,limit_nA=.3)))for kp in(1.2,1.4,1.8)for kd in(.5,.9)]
 if a.rail_ladder:
  cases=[('baseline',0.,None,1.,{})]+[(f'rail_{offset}',0.,None,1.,dict(rail_height_offset=offset))for offset in(-.005,-.01,-.02,-.03)]+[(f'rail_posture_{offset}',0.,None,1.,dict(rail_height_offset=offset,pelvis_pf_posture=dict(target_deg=0.,kp=1.,kd=.15,limit_nA=.6)))for offset in(-.01,-.02)]
 manifest=dict(source_study=str(parent),base_parameters=best['params'],cases=cases,scope='Initial anterior pelvis tilt with compensating hip coordinates, not maintained tilt control; extensor tonic bias is a tested stance-duration hypothesis; fixed rail/guide/muscle physics',double_support_definition='Both feet superior GRF >0.05 BW, after1s; actual timestep samples; report bouts to expose brief impact contacts',selection='Existing strict stepping/knee gate, DS at least0.03 above baseline, assistance no more than baseline+0.02BW, mean swing knee flexion no less than baseline-5deg; then favor DS near0.20 and low assistance/impact. Engineering comparison targets, not human norms.')
 (a.root/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));rows=[]
 for label,tilt,tonic,hip,extra in cases:
  c=dict(cfg);c.update(extra);q=np.array(c['initial_qpos']);theta=np.radians(tilt);q[3:7]=[np.cos(theta/2),0.,0.,-np.sin(theta/2)]
  for side in ('L','R'):
   j=m.joint('hip_'+side);q[int(j.qposadr[0])]+=theta*float(m.jnt_axis[j.id,2])
  c['initial_qpos']=q.tolist();c['hip_flexor_pf_mn_scale']=hip
  if tonic is not None:c['symmetric_rg_ext_tonic']=tonic
  with(a.root/(label+'.log')).open('w',encoding='utf-8')as log,contextlib.redirect_stdout(log):
   run(newer=True,duration=10.,tag=label,output_root=a.root/'trials',**c);folder=a.root/'trials/W2L_VirtualWalker_Newer'/label;measure=metrics(folder)
  row=dict(case=label,initial_anterior_tilt_deg=tilt,extensor_tonic_nA=tonic,hip_flexor_scale=hip,folder=str(folder),configuration=c,metrics=measure);rows.append(row);(a.root/'STATUS.json').write_text(json.dumps(rows,indent=2))
 baseline=rows[0]['metrics'];eligible=[r for r in rows[1:]if r['metrics']['passed']and r['metrics']['double_support_fraction']>=baseline['double_support_fraction']+.03 and r['metrics']['upward_assistance_BW']<=baseline['upward_assistance_BW']+.02 and r['metrics']['mean_swing_knee_flexion_deg']>=baseline['mean_swing_knee_flexion_deg']-5.]
 score=lambda r:abs(r['metrics']['double_support_fraction']-.2)+r['metrics']['upward_assistance_BW']+.002*r['metrics']['foot_force_peak_BW']
 selected=min(eligible,key=score)if eligible else None
 (a.root/'COMPLETE.json').write_text(json.dumps(dict(cases=rows,eligible=len(eligible),selected=selected),indent=2))
 if selected:
  with(a.root/'render.log').open('w')as log,contextlib.redirect_stdout(log):render(Path(selected['folder']),5.,8.2,title='W2L108 | measured stance-overlap test')
 print(json.dumps(dict(cases=len(rows),eligible=len(eligible),selected=None if selected is None else selected['case'])))
if __name__=='__main__':main()
