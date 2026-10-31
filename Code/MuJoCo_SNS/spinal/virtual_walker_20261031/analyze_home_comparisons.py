"""Read-only saved trajectory/FK comparison; no dynamics or fitting."""
import json,math
import numpy as np
from import_virtual_walker import HERE,mujoco
base=HERE/'runs/W2L_VirtualWalker_Newer'
tags=['canonical108_ground_baseline_v1','canonical108_ground_nearflat_source_phase_v1','canonical108_ground_nearflat_ankle_crossed_v1','canonical108_ground_nearflat_ankle_crossed_20s_v1']
tags+=['canonical108_ground_nearflat_crossed_'+x+'_v1' for x in ['contact_RG_off','contact_PF_off','contact_MN_off','all_contact_outputs_off']]
rows=[]
for tag in tags:
 folder=base/tag;z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());a=json.loads((folder/'ASSESSMENT.json').read_text());m=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));d=mujoco.MjData(m);bw=r['robot_mass_kg']*9.81;mask=z['t']>=1.;support={}
 for name,key in [('rails','rail_force_world'),('guides','guide_force_world'),('harness','harness_force')]:
  f=z[key][mask,1];support[name]=dict(net_BW=float(np.mean(f)/bw),upward_BW=float(np.mean(np.maximum(f,0))/bw),downward_BW=float(np.mean(np.maximum(-f,0))/bw))
 foot=z['foot_vertical'][mask]/bw;contact={name:float(np.mean(z['contact_sensor_binary'][mask,i])) for i,name in enumerate(['heel_L','toe_L','heel_R','toe_R'])};pitch=[[],[]];loaded_pitch=[[],[]];cap_frames=0;ids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'foot_'+s) for s in 'LR'];capids={g for g in range(m.ngeom) if 'cap' in (mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g) or '').lower()};arena=[];belt=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_Treadmill_Belt')
 for i in np.flatnonzero(mask)[::5]:
  d.qpos[:]=z['q'][i];mujoco.mj_forward(m,d);cap_frames+=int(any(c.geom1 in capids or c.geom2 in capids for c in d.contact));edge=d.geom_xpos[belt,0]+m.geom_size[belt,0]
  arena.append(float(edge-max(z['foot_body_position'][i,:,0])))
  for s,bid in enumerate(ids):
   R=d.xmat[bid].reshape(3,3);p=math.degrees(math.atan2(R[1,0],R[0,0]));pitch[s].append(p)
   if z['foot_vertical'][i,s]>.05*bw:loaded_pitch[s].append(p)
 row=dict(tag=tag,gate=a['alternating_supported_stepping_gate'],swings={s:a['legs'][s]['unloaded_clear_swings'] for s in 'LR'},alternation=a['alternating_transition_fraction'],knee_min_deg=r['anatomical_knee_min_deg'],root_pitch_range_deg=[float(np.min(z['pelvis_euler_XYZ_deg'][mask,2])),float(np.max(z['pelvis_euler_XYZ_deg'][mask,2]))],root_forward_relative_belt_m=r['root_forward_relative_to_belt_m'],root_world_dx_m=r['root_x_change_m'],support=support,upward_assistance_BW=sum(x['upward_BW'] for x in support.values()),contact_sensor_duty=contact,total_foot_force_p95_BW=float(np.percentile(np.sum(foot,axis=1),95)),total_foot_force_peak_N=float(np.max(np.sum(foot,axis=1))*bw),loaded_foot_pitch_median_deg=[float(np.median(x)) if x else None for x in loaded_pitch],loaded_foot_pitch_range_deg=[[min(x),max(x)] if x else None for x in loaded_pitch],cap_contact_frames_at_10ms_FK_sampling=cap_frames,belt_front_margin_to_mainfoot_centers_min_m=min(arena),scope='Contact-geometry sampling every10ms, not force convergence; foot-front tips are farther anterior than mainfootcenters. Up/down split is each measured net-channel sign, not individualcontact decomposition.')
 (folder/'GROUND_COMPARISON.json').write_text(json.dumps(row,indent=2));rows.append(row)
(HERE/'runs/CANONICAL108_GROUND_COMPARISON.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
