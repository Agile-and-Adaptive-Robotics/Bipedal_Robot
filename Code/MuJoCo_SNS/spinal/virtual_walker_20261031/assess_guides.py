"""Declared stance dwell/strike metrics and guide perturbation diagnostics."""
import json,math
from pathlib import Path
import numpy as np
from import_virtual_walker import mujoco
from assess_trial import assess
def stance_states(force,dt,weight):
 # Offline analysis only; the neural controller receives no filtered signal.
 width=max(1,round(.01/dt));f=np.convolve(force,np.ones(width)/width,mode='full')[:len(force)];on=.025*weight;off=.0125*weight;non=max(1,round(.04/dt));noff=max(1,round(.02/dt));state=False;count=0;states=np.zeros(len(force),bool);strikes=[]
 for i,v in enumerate(f):
  if not state:
   count=count+1 if v>=on else 0
   if count>=non:state=True;start=i-non+1;strikes.append(start);states[start:i+1]=True;count=0
  else:
   states[i]=True;count=count+1 if v<=off else 0
   if count>=noff:state=False;states[i-noff+1:i+1]=False;count=0
 return states,strikes
def analyze(folder):
 folder=Path(folder);z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());basic=assess(folder);t=z['t'];dt=r['physics_dt'];weight=r.get('robot_mass_kg',41.919002)*9.81;mask=t>=1.;states=[];feet={}
 # The actual belt translation follows its source negative world-X slider.
 import xml.etree.ElementTree as E
 model=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));j=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'Treadmill_Slide');assert j>=0;belt=-z['q'][:,model.jnt_qposadr[j]]
 for si,side in enumerate(('L','R')):
  state,strikes=stance_states(z['foot_force_world'][:,si,1],dt,weight);states.append(state);strikes=[k for k in strikes if t[k]>=1.];locations=[]
  for k in strikes:
   stop=min(len(t),k+max(1,round(.04/dt)));cop=z['foot_cop_world'][k:stop,si];finite=np.isfinite(cop[:,0]);copmean=np.mean(cop[finite],axis=0) if np.any(finite) else np.full(3,np.nan);locations.append(dict(time_s=float(t[k]),contact_position_world_m=copmean.tolist(),foot_center_world_m=z['foot_body_position'][k,si].tolist(),belt_displacement_m=float(belt[k])))
  qualified=[];used=set()
  for swing in basic['legs'][side]['swing_intervals']:
   candidates=[(i,p) for i,p in enumerate(locations) if i not in used and swing['end_s']-.02<=p['time_s']<=swing['end_s']+.35]
   if candidates:i,p=candidates[0];used.add(i);qualified.append(p)
  strides=[qualified[i+1]['foot_center_world_m'][0]-qualified[i]['foot_center_world_m'][0]-(qualified[i+1]['belt_displacement_m']-qualified[i]['belt_displacement_m']) for i in range(len(qualified)-1)];feet[side]=dict(contact_onsets_including_recontacts=locations,qualified_strikes_after_wholefoot_swing=qualified,relative_belt_stride_lengths_m=strides,relative_belt_stride_mean_m=float(np.nanmean(strides)) if strides else None,confirmed_stance_fraction=float(np.mean(state[mask])),raw_load_above_threshold_fraction=float(np.mean(z['foot_force_world'][mask,si,1]>=.025*weight)))
 raw_double=np.all(z['foot_force_world'][:,:,1]>=.025*weight,axis=1);confirmed=states[0]&states[1];pulse=r.get('yaw_perturb_start_s');windows={}
 if pulse is not None:
  for name,a,b in [('before',pulse-.5,pulse),('early_response',pulse,pulse+.5),('late_response',pulse+.5,min(pulse+2.,t[-1]))]:
   win=(t>=a)&(t<b)
   if np.any(win):windows[name]=dict(guide_contact_fraction=float(np.mean(z['guide_contact_count'][win]>0)),yaw_range_deg=float(np.ptp(z['pelvis_euler_XYZ_deg'][win,1])),yaw_omega_abs_peak_rad_s=float(np.max(abs(z['pelvis_omega_world'][win,1]))),rotational_energy_peak_J=float(np.max(z['pelvis_rotational_energy_J'][win])),guide_yaw_torque_abs_peak_Nm=float(np.max(abs(z['guide_torque_about_pelvis'][win,1]))))
 gap=None
 if r.get('guide_geometry'):
  values=[];roof=r['guide_geometry']['cap_inner_height_m']
  for q,p in zip(z['pelvis_quaternion'][mask],z['root'][mask]):
   mat=np.empty(9);mujoco.mju_quat2Mat(mat,q);R=mat.reshape(3,3);tops=[p[1]+(R@np.array([0,-.04,s*.25]))[1]+np.abs(R[1])@np.array([.06,.0125,.15]) for s in (-1,1)];values.append(roof-max(tops))
  gap=dict(median_m=float(np.median(values)),minimum_m=float(np.min(values)),scope='Roof-to-whole-wing bounding-box vertical separation, independent of lateral overlap; not an exact Euclidean collision distance.')
 relative=z['root'][:,0]-belt;steady_progress=float(relative[-1]-relative[np.flatnonzero(mask)[0]])
 out=dict(source=str(folder.resolve()),neural_time_scale=r['neural_time_scale'],guide_clearance_m=r.get('guide_clearance_m'),harness_scale=r['sagittal_harness_scale'],failure=r['failure'],mean_net_vertical_assistance_BW=r.get('mean_net_vertical_assistance_BW'),root_forward_relative_to_belt_m=r['root_forward_relative_to_belt_m'],evaluation_window_relative_belt_progress_m=steady_progress,stance_protocol=dict(load_on_BW=.025,load_off_BW=.0125,causal_average_window_s=.01,on_dwell_s=.04,off_dwell_s=.02,scope='Offline state estimation only; raw forces retained and not changed in controller or primary unloaded-clearance gate.'),raw_simultaneous_contact_fraction=float(np.mean(raw_double[mask])),dwell_qualified_double_support_fraction=float(np.mean(confirmed[mask])),feet=feet,pelvis_XYZ_range_deg=r['pelvis_euler_XYZ_range_deg'],pelvis_rotational_energy_max_J=r['pelvis_rotational_energy_max_J'],guide_contact_fraction=r['guide_contact_fraction'],guide_force_rms_N=r['guide_force_rms_N'],guide_torque_rms_Nm=r['guide_torque_rms_Nm'],cap_vertical_gap=gap,perturbation_windows=windows,stepping_gate=basic['alternating_supported_stepping_gate'],unloaded_clear_swings=[basic['legs'][s]['unloaded_clear_swings'] for s in ('L','R')],alternating_fraction=basic['alternating_transition_fraction'])
 (folder/'GUIDE_ASSESSMENT.json').write_text(json.dumps(out,indent=2));return out
