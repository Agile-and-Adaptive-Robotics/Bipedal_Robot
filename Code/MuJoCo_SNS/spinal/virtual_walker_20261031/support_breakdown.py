"""Read saved reactions and replay contact geometry; no dynamics integration."""
import json
import numpy as np
from import_virtual_walker import HERE,mujoco
folder=HERE/'runs/W2L_VirtualWalker_Newer_GUI1nA/newer_v3_diag_v1_ground_guides';z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());weight=r['robot_mass_kg']*9.81;out={}
for label,mask in [('whole_trial',np.ones(len(z['t']),bool)),('after_1s',z['t']>=1.)]:
 row={}
 for name,field in [('rails','rail_force_world'),('guides','guide_force_world'),('harness','harness_force')]:
  fy=z[field][mask,1];row[name]=dict(net_mean_BW=float(np.mean(fy)/weight),upward_net_channel_mean_BW=float(np.mean(np.maximum(fy,0))/weight),downward_net_channel_mean_BW=float(np.mean(np.maximum(-fy,0))/weight),peak_upward_net_channel_N=float(max(0,fy.max())),peak_downward_net_channel_N=float(max(0,-fy.min())))
 row['net_total_mean_BW']=sum(v['net_mean_BW'] for v in row.values());out[label]=row
m=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));d=mujoco.MjData(m);capids={i for i in range(m.ngeom) if (mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,i) or '').startswith('guide_cap_')};frames=0;contacts=0
for q in z['q']:
 d.qpos[:]=q;mujoco.mj_forward(m,d);n=sum(c.geom1 in capids or c.geom2 in capids for c in d.contact);contacts+=n;frames+=int(n>0)
out.update(cap_contact_geometry_frames=frames,cap_contact_geometry_count=contacts,frames_audited=len(z['q']),scope='Saved measured reaction channels; positive/negative split is per-time net channel, not sum of individual contacts. Cap geometry audited at every saved q using forward kinematics/contact detection only, no dynamics steps or controller rerun.');(folder/'SUPPORT_BREAKDOWN.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
