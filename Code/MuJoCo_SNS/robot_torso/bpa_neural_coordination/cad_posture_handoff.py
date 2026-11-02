"""Read-only handoff/figures for existing CAD balance and rejected IK proposal."""
import json,hashlib,shutil
from pathlib import Path
import numpy as np
from scipy.spatial import ConvexHull
import cad_posture_audit as P

def jsonsafe(value):
 if isinstance(value,dict):return {key:jsonsafe(item) for key,item in value.items()}
 if isinstance(value,list):return [jsonsafe(item) for item in value]
 if isinstance(value,float) and not np.isfinite(value):return None
 return value

def main():
 m,dyn=P.model();a=P.load(P.S.ROOT/'mixed/mixed/ACTUATOR_BRIDGE.npz');saved=P.load(P.OUT/'OFFLINE_ANTERIOR_REFERENCE.npz');candidate=json.loads((P.OUT/'OFFLINE_ANTERIOR_REFERENCE.json').read_text());base=json.loads((P.OUT/'POSTURE_GEOMETRY_AUDIT.json').read_text());pelvis=m.body('pelvis').id;rows={};points={}
 for name,q in [('initial',a['qpos_full'][0]),('accepted85tail',a['qpos_full'][-1]),('candidate',saved['candidate_qpos'])]:
  d=P.mujoco.MjData(m);d.qpos[:]=q
  if name=='accepted85tail':d.qvel[:]=a['qvel_full'][-1];d.qfrc_applied[:]=a['qfrc'][-2];d.xfrc_applied[:]=a['rig'][-2]
  dyn.update(d);g=P.geometry(m,d,dyn);closures={};headings={};heights={};sites=[]
  for side in ('r','l'):
   e=m.equality('knee_loop_'+side).id;one,two=m.eq_obj1id[e],m.eq_obj2id[e];p1=d.xpos[one]+d.xmat[one].reshape(3,3)@m.eq_data[e,:3];p2=d.xpos[two]+d.xmat[two].reshape(3,3)@m.eq_data[e,3:6];closures[side]=float(np.linalg.norm(p1-p2));bid=m.body('foot_'+side).id;forward=d.xmat[bid].reshape(3,3)[:,0];local=d.xmat[pelvis].reshape(3,3).T@forward;headings[side]=dict(world_deg=float(np.rad2deg(np.arctan2(forward[2],forward[0]))),pelvis_relative_deg=float(np.rad2deg(np.arctan2(local[2],local[0]))))
   values=[]
   for gid in range(m.ngeom):
    if not m.geom(gid).name.startswith('foot_'+side+'_'):continue
    mid=m.geom_dataid[gid]
    if mid>=0:
     vertices=m.mesh_vert[m.mesh_vertadr[mid]:m.mesh_vertadr[mid]+m.mesh_vertnum[mid]];values.append(float((vertices@d.geom_xmat[gid].reshape(3,3)[1]+d.geom_xpos[gid,1]).min()))
   heights[side]=min(values)
  for c in d.contact:
   names=[m.geom(int(k)).name for k in (c.geom1,c.geom2)]
   if 'ground' in names and any(n.startswith('foot_') for n in names):sites.append(c.pos[[0,2]].copy())
  rows[name]=dict(knee_loop_anchor_separation_m=closures,foot_heading=headings,foot_visual_mesh_min_y_m=heights,COM=g['COM'].tolist(),CoP_xz=g['CoP_xz'].tolist(),trunk_forelean_deg=float(g['trunk_forelean_deg']),pelvis_forelean_deg=float(g['pelvis_forelean_deg']));points[name]=np.asarray(sites)
 tail=base['cases'];toe_fractions={}
 for key,value in tail.items():
  loads=np.array(value['tail_mean_contact_loads']);toe_fractions[key]=float(loads[:,1].sum()/loads.sum())
 summary=dict(scope='Existing-state geometry audit and rejected offlineproposal; no timeintegration',pose_details=rows,toe_load_fraction=toe_fractions,candidate_rejected=True,candidate_required_CoP_margin_m=candidate['required_CoP_polygon_margin_m'],candidate_static_capacity_feasible=candidate['capacity_feasible'],current_hardware_ROM_certified=False)
 (P.OUT/'POSTURE_INTERPRETATION.json').write_text(json.dumps(jsonsafe(summary),indent=2,allow_nan=False))
 geometryfile=P.OUT/'POSTURE_GEOMETRY_AUDIT.json';geometryfile.write_text(json.dumps(jsonsafe(base),indent=2,allow_nan=False))
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False});fig,axes=plt.subplots(1,2,figsize=(11,5),constrained_layout=True)
 for ax,name in zip(axes,('accepted85tail','candidate')):
  pts=points[name];hull=ConvexHull(pts);v=pts[hull.vertices];v=np.vstack((v,v[0]));ax.fill(1000*v[:,0],1000*v[:,1],color='.9',label='Geometric contact polygon');ax.plot(1000*pts[:,0],1000*pts[:,1],'kx',label='Contact points');com=np.array(rows[name]['COM']);ax.plot(1000*com[0],1000*com[2],'o',color='#0077BB',label='COM projection')
  cop=np.array(rows[name]['CoP_xz']) if name=='accepted85tail' else np.array(candidate['required_static_CoP_xz']);ax.plot(1000*cop[0],1000*cop[1],'s',color='#EE7733',label='Actual CoP' if name=='accepted85tail' else 'Required static CoP');ax.set_title('Accepted85%tail: supported' if name=='accepted85tail' else '+10mm anterior/uprighttrunk: REJECTED');ax.set_xlabel('Anterior coordinate (mm)');ax.set_ylabel('Right coordinate (mm)');ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=8);ax.set_aspect('equal',adjustable='datalim')
 fig.suptitle('Ground support geometry: anterior COM movement alone does not establish feasible standing\nThe unchanged orientation harness makes the new proposal require CoP outside the contact polygon',fontsize=11)
 for ext in ('png','pdf','svg'):fig.savefig(P.OUT/('CONTACT_POLYGON_AND_REJECTED_REFERENCE.'+ext),dpi=200)
 plt.close(fig)
 record=P.OUT/'EXPERIMENT_RECORD.md';text=record.read_text(encoding='utf-8')
 if '## Offline anterior reference: rejected' not in text:
  text+='\n\n## Offline anterior reference: rejected\n\nOne new proposal requested10mm anterior COM movement, actual savedfootbodyposes fixed, uprighttrunk, and localrotationcomponentchanges bounded±5°. This was an offlinekinematic/capacitycheck, not an appliedcontrollerreference or simulatedbehavior. The IK meets the geometric targets (footpositionerror1.17e−15m, orientationerror1.28e−13rad); knees remain−8.65/−11.22°. However, unchangedorientationguidance would exert+20.25Nm sagittal and−9.63Nm yaw torque. RequiredstaticCoP is26.22mm outside the availablegeometriccontactpolygon, and the54route equilibrium LP fails. Pairpressureinversion remainsvalid at497.47/432.94kPa. Thus this reference is REJECTED before a closedlooprun; no changedneuralweights or gainfit was used. This does not prove every moreanteriorreference impossible.\n\nThe earlier draft mistakenly fixed knee-projectedfootposes ratherthan the actualsavedfootposes. Exactfourbarclosure projection shiftsfeet5.69/6.69mm; the draft/source/protocol are preserved and excluded. The correctedproposal retainsactualsavedfootposes. Actualsavedkneeloopanchor gaps and footheadings are listed in POSTURE_INTERPRETATION.json; no physicalconstraint or geometry was modified. No finiteCADhardwareROM is certified.\n\nToe vertical-loadfraction improves from'+f"{100*toe_fractions['parent_mixed']:.2f}% to{100*toe_fractions['support085_mixed']:.2f}%"+' under reducedsupport, despite the larger posteriorCOMexcursion. This distinction prevents interpreting one scalar supportscore or a trunkangle alone as improvedbalance. Existingneuralcontactsensors still receiveonlytheir originalverticalheel/toechannels; full forcecoordinates here are offlineanalysis, not addedneuralinputs.\n\nNext mechanicalcondition should be selected explicitly: a footplacement/posturereference and apparatusorientationreference compatiblewith staticCoP and routecapacity. Footheading/initialflightpose/softfourbarclosure deserve review before more neuralgain or verticalsupport searches. No furthercondition, harnessheightchange, or manuscriptedit was performed.\n'
 if '## Mechanical baseline limits' not in text:
  text+='\n\n## Mechanical baseline limits\n\nThe savedinitialpose has11.63/11.70mm knee-loop connect-anchor gaps. At the selected accepted85%tail snapshot these remain5.52/6.50mm. These are actual connect-anchor separations, not merely visualmesh gaps. The old settlinggates did not test loopclosure. Left/right footheadings are approximately−30.38/+29.13° at that snapshot; sourceheadingfeedbackgain is0. Initialleftvisualfoot is8.13mm aboveground and right is−1.64mm. Thus the sourceframe1069 reference is neither an exactly closedmechanism nor a squared, groundedstandingIK target. Previously acceptedsettling results retain their original measuredvalues/gates, with this newly identified physics limitation; no rigidCADhardwarefidelity is established. A new mechanical/posturalbaseline should be established before more weight-support or neuralgain searches. Undefined unloadedCoP values are represented as null in JSON and NaN in numericalarchives.\n'
 record.write_text(text,encoding='utf-8');snap=P.OUT/'source_snapshot';snap.mkdir(exist_ok=True)
 files=[Path(P.__file__),Path(__file__)];ledger={}
 for file in files:target=snap/file.name;shutil.copy2(file,target);ledger[str(file)]=dict(sha256=hashlib.sha256(file.read_bytes()).hexdigest(),snapshot=str(target))
 (snap/'MANIFEST.json').write_text(json.dumps(ledger,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
