"""Read-only stage3 reference, support geometry, and PT unit audit.

Uses saved actual model/state/forces; forward-position evaluations only.
No integration, controller step, fit, or change to historical acceptance.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,hashlib
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np,mujoco
from scipy.spatial import ConvexHull,QhullError


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def stats(x):
    x=np.asarray(x)
    return {k:np.asarray(v).tolist() for k,v in dict(mean=np.mean(x,axis=0),std=np.std(x,axis=0),minimum=np.min(x,axis=0),maximum=np.max(x,axis=0),peak_to_peak=np.ptp(x,axis=0)).items()}


def hull_summary(points,com):
    try:
        h=ConvexHull(np.unique(np.asarray(points)[:,:2],axis=0));vertices=np.asarray(points) # use unique point indexing below
        p=h.points[h.vertices];cross=p[:,0]*np.roll(p[:,1],-1)-np.roll(p[:,0],-1)*p[:,1];den=3*cross.sum()
        center=np.sum((p+np.roll(p,-1,axis=0))*cross[:,None],axis=0)/den
        margin=np.min(-(h.equations[:,:2]@com[:2]+h.equations[:,2])/np.linalg.norm(h.equations[:,:2],axis=1))
        return {'centroid_xy_m':center.tolist(),'x_span_m':[float(p[:,0].min()),float(p[:,0].max())],'y_span_m':[float(p[:,1].min()),float(p[:,1].max())],'COM_inside':bool(margin>=0),'COM_edge_margin_m':float(margin),'area_m2':float(h.volume)}
    except (QhullError,ValueError):return None


def hip_followers(modelpath):
    raw=Path(__file__).resolve().parents[4]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
    m=mujoco.MjModel.from_binary_path(str(modelpath));rows=[]
    for e in ET.parse(raw).getroot().find('equality').findall('joint'):
        driver=e.get('joint2','')
        if not driver.startswith(('hip_flexion_','hip_adduction_')):continue
        follower=e.get('joint1');a=np.fromstring(e.get('polycoef','0 1 0 0 0'),sep=' ')
        match=[i for i in range(m.neq)if m.eq_type[i]==int(mujoco.mjtEq.mjEQ_JOINT)and m.joint(m.eq_obj1id[i]).name==follower and m.joint(m.eq_obj2id[i]).name==driver]
        if len(match)!=1:raise AssertionError('Source hip follower not uniquely matched')
        b=m.eq_data[match[0],:5];rows.append({'driver':driver,'follower':follower,'raw_polycoef':a.tolist(),'runtime_polycoef':b.tolist(),'coefficients_exact':bool(np.array_equal(a,b)),'nonzero_odd_coefficients':bool(np.any(a[[1,3]]!=0.))})
    return {'source_raw_XML':str(raw),'source_sha256':sha(raw),'runtime_model':str(modelpath),'rows':rows,'count':len(rows),'nonzero_odd_polynomials':sum(r['nonzero_odd_coefficients']for r in rows),'scope':'If q_native=-q_OpenSim, polynomial composition f(-q) negates odd coefficients only. The inspected hip-driven pathpoint polynomials are even, so they are invariant under that sign reversal. This structural check does not validate the full muscle route, force law, or every IK adapter.'}


def evaluate(cfgpath):
    cfg=json.loads(cfgpath.read_text());result=json.loads(cfgpath.with_suffix('.result.json').read_text());base=cfgpath.with_suffix('');z=np.load(str(base)+'.npz',allow_pickle=True);caps=Path(cfg['verification_capture']);meta=json.loads((caps/'ALIGNED_TRACE_METADATA.json').read_text());cap=np.load(caps/'aligned_trace.npz');modelpath=caps/meta['model_file'];m=mujoco.MjModel.from_binary_path(str(modelpath));d=mujoco.MjData(m)
    ref=np.array(result['balance_calibration']['reference_native_xy_m']);t=z['t'];com=z['com'];windows={}
    assist=np.load(str(base)+'.assistance.npz',allow_pickle=True);acols=assist['columns'].tolist();asamp=assist['samples'];at=asamp[:,acols.index('t')];BW=result['measured']['body_weight_N']
    for start in (.6,4.,5.,7.):
        keep=(t>=start)&(t<8.);c=com[keep,:2];delta=c-ref;center=c.mean(axis=0);centered=c-center;ak=(at>=start)&(at<8.);aa=asamp[ak]
        foot=(aa[:,acols.index('right_vertical_GRF_N')]+aa[:,acols.index('left_vertical_GRF_N')])/BW;spring=aa[:,acols.index('spring_support_N')]/BW;explicit=aa[:,acols.index('explicit_support_N')]/BW;rig=aa[:,[acols.index('harness_Fx_N'),acols.index('harness_Fy_N')]]
        windows[str(start)]={'COM_native_xy_m':stats(c),'mean_offset_from_fixed_talus_reference_xy_m':delta.mean(axis=0).tolist(),'mean_offset_radius_m':float(np.linalg.norm(delta.mean(axis=0))),'maximum_distance_from_talus_reference_m':float(np.max(np.linalg.norm(delta,axis=1))),'RMS_distance_from_talus_reference_m':float(np.sqrt(np.mean(np.sum(delta**2,axis=1)))),'mean_centered_RMS_radius_m':float(np.sqrt(np.mean(np.sum(centered**2,axis=1)))),'mean_centered_maximum_radius_m':float(np.max(np.linalg.norm(centered,axis=1))),'first_last_halfsecond_mean_COM_change_m':(c[-min(250,len(c)):].mean(axis=0)-c[:min(250,len(c))].mean(axis=0)).tolist(),'once_counted_foot_load_BW':stats(foot),'spring_support_BW':stats(spring),'explicit_support_BW':stats(explicit),'combined_vertical_load_BW_mean':float(np.mean(foot+spring+explicit)),'horizontal_rig_force_native_xy_N':stats(rig),'horizontal_rig_RMS_N':float(np.sqrt(np.mean(np.sum(rig**2,axis=1))))}
    replica=windows['0.6']['maximum_distance_from_talus_reference_m']
    if abs(replica-result['assessment']['sway_m'])>1e-12:raise AssertionError('Ankle-reference metric reconstruction mismatch')
    owners={m.body(b+'_'+s).id for s in ('r','l')for b in ('talus','calcn','toes')};foot_geoms=[];mesh=[]
    for g in range(m.ngeom):
        b=int(m.geom_bodyid[g])
        while b and b not in owners:b=int(m.body_parentid[b])
        if not b:continue
        foot_geoms.append(g)
        if m.geom_type[g]==int(mujoco.mjtGeom.mjGEOM_MESH):
            mid=int(m.geom_dataid[g]);a=int(m.mesh_vertadr[mid]);n=int(m.mesh_vertnum[mid]);mesh.append((g,m.mesh_vert[a:a+n].copy()))
    plane=m.geom('ground-plane').id;geometric_rows=[];contact_rows=[];talus=[];moving_com=[];fk_errors=[];moving_ankle_errors=[]
    for i in range(0,len(cap['t_state']),10):
        d.qpos[:]=cap['qpos'][i];mujoco.mj_fwdPosition(m,d);time=float(cap['t_state'][i]);c=d.subtree_com[0];vertices=np.vstack([v@d.geom_xmat[g].reshape(3,3).T+d.geom_xpos[g]for g,v in mesh]);h=hull_summary(vertices,c)
        if h:h['time_s']=time;geometric_rows.append(h)
        points=[ct.pos.copy()for ct in d.contact[:d.ncon] if ((ct.geom1==plane and ct.geom2 in foot_geoms)or(ct.geom2==plane and ct.geom1 in foot_geoms))and ct.dist<=0.]
        h=hull_summary(points,c)if len(points)>=3 else None
        if h:h['time_s']=time;contact_rows.append(h)
        footcenter=np.mean([d.xpos[m.body('talus_'+s).id,:2]for s in ('r','l')],axis=0)
        talus.append(footcenter);moving_com.append(d.subtree_com[m.body('pelvis').id].copy())
        raw_index=int(round((time-m.opt.timestep)/m.opt.timestep))
        if not np.array_equal(cap['qpos'][i],z['qfull'][raw_index]):raise AssertionError('Saved raw/aligned pose pairing differs')
        fk_errors.append(float(np.max(abs(c-z['com'][raw_index]))))
        if time>=5.:moving_ankle_errors.append((c[:2]-footcenter).copy())
    def summarize_hulls(rows):
        selected=[r for r in rows if r['time_s']>=5.]
        return {'samples':len(selected),'COM_inside_fraction':float(np.mean([r['COM_inside']for r in selected])),'COM_edge_margin_m':stats([r['COM_edge_margin_m']for r in selected]),'centroid_native_xy_m':stats([r['centroid_xy_m']for r in selected]),'x_span_m':stats([r['x_span_m']for r in selected]),'y_span_m':stats([r['y_span_m']for r in selected])}if selected else {'samples':0}
    initial={};q0=z['qfull'][0].copy();d.qpos[:]=q0;mujoco.mj_fwdPosition(m,d)
    orientation={'calcn_to_toes_world_m':{s:(d.xpos[m.body('toes_'+s).id]-d.xpos[m.body('calcn_'+s).id]).tolist()for s in ('r','l')},'pelvis_local_plusX_world':d.xmat[m.body('pelvis').id].reshape(3,3)[:,0].tolist(),'torso_local_plusX_world':d.xmat[m.body('torso').id].reshape(3,3)[:,0].tolist(),'compiled_pelvis_body_quaternion_wxyz':m.body_quat[m.body('pelvis').id].tolist(),'runtime_hip_axes_local':{s:m.jnt_axis[m.joint('hip_flexion_'+s).id].tolist()for s in ('r','l')},'initial_COM_FK_minus_logged_m':(d.subtree_com[0]-z['com'][0]).tolist(),'initial_talus_midpoint_native_xy_m':np.mean([d.xpos[m.body('talus_'+s).id,:2]for s in ('r','l')],axis=0).tolist(),'conclusion':'Both forefoot landmarks and pelvis/torso anatomical+X point toward world+X. No global180deg facing reversal; actual runtime hip axes must be distinguished from raw XML.'}
    for side in ('r','l'):
        hip=m.joint('hip_flexion_'+side);knee=m.joint('knee_angle_'+side)
        def thigh_pitch(q):
            d.qpos[:]=q;mujoco.mj_kinematics(m,d);v=d.xanchor[knee.id]-d.xanchor[hip.id];return np.arctan2(v[0],-v[2])
        q1=q0.copy();q1[hip.qposadr[0]]+=1e-5;slope=(thigh_pitch(q1)-thigh_pitch(q0))/1e-5
        initial[side]={'hip_native_deg':float(np.degrees(q0[hip.qposadr[0]])),'knee_native_deg':float(np.degrees(q0[knee.qposadr[0]])),'ankle_native_deg':float(np.degrees(q0[m.joint('ankle_angle_'+side).qposadr[0]])),'thigh_world_pitch_anterior_from_vertical_deg':float(np.degrees(thigh_pitch(q0))),'d_anterior_thigh_pitch_d_positive_hip_q':float(slope),'positive_hip_q_anatomical_action':'flexion'if slope>0 else'extension'}
    bal=np.load(str(base)+'.balance.npz');b=bal['state_and_command'];bt=bal['t'];bk=(bt>=5.)&(bt<8.)
    regional=np.load(str(base)+'.regional_contacts.npz');rk=(regional['t']>=5.)&(regional['t']<8.)
    mech=np.load(str(base)+'.mechanics.npz',allow_pickle=True);mk=(mech['t']>=5.)&(mech['t']<8.);names=mech['names'].tolist();pt={name:stats(mech['voltage'][mk,names.index(name)])for name in names if name.startswith(('PT_BAL_','PT_IN_'))}
    return {'configuration':str(cfgpath),'source_hashes':{'configuration':sha(cfgpath),'result':sha(cfgpath.with_suffix('.result.json')),'captured_model':sha(modelpath)},'historical_assessment_unchanged':result['assessment'],'frame':'native x anterior,y left,z superior; independently checked below','global_orientation_proof':orientation,'PT_fixed_reference_native_xy_m':ref.tolist(),'windows':windows,'source_initial_pose':initial,'late_joint_angles_native_deg':{m.joint(j).name:stats(np.degrees(z['qfull'][z['t']>=5.,m.jnt_qposadr[j]]))for j in range(m.njnt)if m.joint(j).name in [a+'_'+s for s in ('r','l')for a in ('hip_flexion','knee_angle','ankle_angle')]},'actual_ground_alignment':result['compiled']['ground_alignment'],'geometry':{'scope':'20ms sampled FK of actual saved states and exact model. All-foot mesh horizontal projection is a geometric extent,not a load-bearing plantar polygon. Native overlapping foot-plane contacts define a second geometric candidate hull; no contact forces are recomputed and no CoP is inferred.','foot_mesh_geometries':len(mesh),'whole_foot_projected_mesh_hull':summarize_hulls(geometric_rows),'native_overlapping_contact_point_hull':summarize_hulls(contact_rows),'initial_projected_mesh_hull':geometric_rows[0],'talus_midpoint_native_xy_m':stats(talus),'moving_pelvis_subtree_COM_m':stats(moving_com),'late_COM_relative_to_actual_moving_talus_xy_m':stats(moving_ankle_errors),'maximum_FK_minus_logged_COM_m':max(fk_errors),'raw_aligned_qpos_pairing_exact':True,'COM_log_timing_note':'Runner/capture body kinematics retain the last solver evaluation while saved qpos is post-integration. FK evaluates the final qpos; the measured substep discrepancy is reported,not silently retimed.'},'actual_late_regional_vertical_load_N':{'order':['R_hindfoot','R_forefoot','L_hindfoot','L_forefoot'],'mean':regional['force_N'][rk].mean(axis=0).tolist()},'PT_actual_late_delayed_task_displacement_m':stats(b[bk,:2]),'PT_actual_late_command_native_xy_nA':stats(b[bk,6:8]),'PT_recorded_source_relay_voltages_mV':pt,'weight_normalization':{'legacy_total_compiled_BW_N':BW,'actual_pelvis_subtree_BW_N':meta['body_weight_N'],'fixed_body_mass_included_in_legacy_normalization_kg':float(m.body_mass.sum()-m.body_subtreemass[m.body('pelvis').id])},'physics_steps':0,'neural_steps':0,'new_fit':False}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--conditions',nargs='+',default=['source_baseline_s3','candidate_ii02_s3','candidate_weak_s3'])
    a=p.parse_args();conditions=a.conditions
    record={'schema':'stage3_reference_unit_audit_v1',
        'conditions':{c:evaluate(a.root/'configurations'/(c+'.json'))for c in conditions},
        'hip_moving_point_polynomial_sign_composition':hip_followers(a.root/'figures/source_baseline_s3/model.mjb'),
        'S2_S3_saved_PT_calibration':{c:json.loads((a.root/'configurations'/(c+'.result.json')).read_text())['balance_calibration']for c in ('source_baseline_s2','source_baseline_s3','candidate_ii02_s2','candidate_ii02_s3')},
        'unit_contract':{'kp':'nA/m','kd':'nA*s/m','ka':'nA*s^2/m',
            'basis':'Dimensionless activation vector calibrated against a 1 m/s^2 native fixed-pose acceleration target.',
            'source':'TaskBalance.update injects the PD result as applied current; TaskBalance.bind maps basis b to 2*b microSiemens excitation and 0.25*b_opposite microSiemens inhibition; runner maps MN voltage/5 mV to clipped activation. Consequently kp/kd are NOT directly desired-acceleration gains.'},
        'acceptance_unchanged':True,'new_body_or_training':False,'audit_source_sha256':sha(__file__)}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps({'output':str(a.out),'new_body_steps':0,'conditions':conditions}))
