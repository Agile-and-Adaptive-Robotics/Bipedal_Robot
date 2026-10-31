"""Longer matched duration, explicitly extending belt motor command to20s."""
import json,sys
import numpy as np
from run_virtual_walker import run
from import_virtual_walker import HERE,mujoco
from assess_guides import analyze
rows=[]
for scale in (2.,4.,6.):
 tag=f'guides_long_v1_t{scale:g}_c0.03_h0_p5';folder=HERE/'runs/W2L_VirtualWalker'/tag;params=dict(neural_time_scale=scale,guide_clearance=.03,sagittal_harness=0.,rail_height_offset=-.02,yaw_perturb_nm=5.,belt_motor_end_s=20.)
 result=json.loads((folder/'RESULT.json').read_text()) if '--analyze-only' in sys.argv else run(duration=20.,tag=tag,**params);assessment=analyze(folder);z=np.load(folder/'trajectory.npz');mask=z['t']>=1.;model=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));beltbody=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Treadmill_Belt');beltjoint=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'Treadmill_Slide');assert beltbody>0 and beltjoint>=0;initial_data=mujoco.MjData(model);mujoco.mj_forward(model,initial_data);beltx=initial_data.xpos[beltbody,0]-z['q'][:,model.jnt_qposadr[beltjoint]];footmargin=float(np.min(2.-np.abs(z['foot_body_position'][:,:,0]-beltx[:,None])-.15));wingmargin=[]
 railbody=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Left_support_rail');railx=model.body_pos[railbody,0]
 for q,p in zip(z['pelvis_quaternion'],z['root']):
  mat=np.empty(9);mujoco.mju_quat2Mat(mat,q);R=mat.reshape(3,3)
  for side in (-1,1):
   center=p+R@np.array([0,-.04,side*.25]);extent=np.abs(R[0])@np.array([.06,.0125,.15]);wingmargin.append(2.-abs(center[0]-railx)-extent)
 from scipy.signal import find_peaks
 names=z['neuron_names'].tolist();rg=z['neuro'][mask,names.index('L RG flx')];peaks,_=find_peaks(rg,prominence=1.,distance=round(.15/result['physics_dt']));period=float(np.median(np.diff(z['t'][mask][peaks]))) if len(peaks)>1 else None
 e=z['pelvis_rotational_energy_J'][mask];rotomega=z['pelvis_omega_world'][mask];roll=z['pelvis_euler_XYZ_deg'][mask,0];yaw=z['pelvis_euler_XYZ_deg'][mask,1];arena=dict(conservative_wholefoot_inside_finite_belt_margin_min_m=footmargin,wholewing_inside_original_rail_xspan_margin_min_m=float(np.min(wingmargin)),carriage_limit_margin_min_m=float(np.min(1.9-abs(z['guide_carriage_x']))),scope='Foot centers use a conservative150mm envelope; wing axis-aligned projected bounds use actual quaternion and exact source wing dimensions; no arena lengthening or body following.');assessment.update(arena=arena,RG_peak_count=int(len(peaks)),RG_median_period_s=period,pelvis_rotational_energy_mean_J=float(np.mean(e)),pelvis_rotational_energy_p95_J=float(np.percentile(e,95)),pelvis_rotational_energy_peak_poststartup_J=float(np.max(e)),yaw_omega_rms_rad_s=float(np.sqrt(np.mean(rotomega[:,1]**2))),roll_rms_deg=float(np.sqrt(np.mean(roll**2))),yaw_rms_deg=float(np.sqrt(np.mean(yaw**2))))
 (folder/'GUIDE_ASSESSMENT.json').write_text(json.dumps(assessment,indent=2));(folder/'PARAMETERS.json').write_text(json.dumps(dict(contact_si=False,pf_knee_flex_scale=1.,ankle_pf_crossed=False,symmetric_rg_ext_tonic=None,**params),indent=2));rows.append(dict(name=tag,result=result,assessment=assessment));(HERE/'runs/GUIDE_LONG_COMPARISON_V1.json').write_text(json.dumps(rows,indent=2))
