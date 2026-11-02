"""Frozen sagittal reference-envelope gates; never clip simulation states.

These are the saved retargeted IK trial's bounds, not universal human RoMs.
They supplement clearance/contact gates and do not establish normal gait.
"""
import os,sys,json,hashlib,argparse
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
HERE=Path(__file__).resolve().parent;S=HERE.parent;sys.path[:0]=[str(S),str(S/'functional_gait2392')]
import numpy as np,mujoco
from tendon_replay_campaign import build,SOURCE
from mechanics import CoupledArms,JOINTS

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def make_specs(folder):
    folder.mkdir(parents=True,exist_ok=True);target=S/'hunt_training/runs/targets/targets.npz';meta=S/'hunt_training/runs/targets/metadata.json';z=np.load(target,allow_pickle=True);md=json.loads(meta.read_text());base=np.load(SOURCE.with_suffix('.npz'));cfg=json.loads(SOURCE.with_suffix('.json').read_text());m,d=build(base['qfull'][0],cfg,1.);ar=CoupledArms(m)
    if z['qpos'].shape[1]!=m.nq:raise ValueError('Target/native model coordinate dimensions differ')
    bounds={}
    for side in ('r','l'):
        for joint in ('hip_flexion','knee_angle','ankle_angle'):
            name=joint+'_'+side;adr=int(m.joint(name).qposadr[0]);values=z['qpos'][:,adr];lo=float(values.min());hi=float(values.max())
            if joint=='knee_angle':hi=min(hi,0.)
            bounds[name]=dict(qpos_address=adr,native_bounds_rad=[lo,hi],native_bounds_deg=np.degrees([lo,hi]).tolist(),native_to_anatomical_sign=md['native_signs'][name])
    spec=dict(format='aarl_ik_envelope_v1',scope='Exact per-side saved retargeted IK sagittal envelope, no added angular margin; knee upper bound intersected with strict0deg. Necessary reference-shape check, not sufficient normal-gait validation.',source_target=str(target),source_target_sha256=sha(target),metadata_sha256=sha(meta),window_neural_s=[5.,13.],numeric_epsilon_rad=1e-10,joints=bounds,other_axes_scope='Frontal separation and strict per-physics-state knee protection remain separate; unidentifiable subtalar/MTP/lumbar patterns are not assigned fictitious gait envelopes.')
    path=folder/'IK_ENVELOPE_SPEC_V1.json'
    if path.exists():raise RuntimeError('Do not overwrite frozen acceptance specification')
    path.write_text(json.dumps(spec,indent=2))
    names=[m.actuator(i).name for i in range(m.nu)];samples=[];indices=np.linspace(0,len(z['qpos'])-1,17,dtype=int)
    for k in indices:
        d.qpos[:]=z['qpos'][k];ar.project(d);mujoco.mj_forward(m,d);samples.append(ar.update(d).copy())
    samples=np.array(samples);roles={}
    for group,joint,sign in (('knee_ext','knee_angle',1.),('ankle_pf','ankle_angle',-1.)):
        roles[group]={}
        for side in ('r','l'):
            row=JOINTS.index(joint+'_'+side)
            for i,name in enumerate(names):
                if not name.endswith('_'+side):continue
                values=samples[:,row,i];roles[group][name]=dict(min_moment_arm_m=float(values.min()),max_moment_arm_m=float(values.max()),desired_sign=sign,valid_extension_route=bool(np.min(sign*values)>.001))
    role_spec=dict(format='aarl_ik_functional_roles_v1',source_envelope_sha256=sha(path),source_target_sha256=sha(target),pose_indices=indices.tolist(),moment_arm_scope='Equality-aware torque per positive pulling force in native joint coordinates; knee extension positive, ankle plantarflexion negative. Role accepted only if correct by>1mm at every sampled IK pose.',roles=roles)
    (folder/'FUNCTIONAL_ROLES_V1.json').write_text(json.dumps(role_spec,indent=2));print(json.dumps(dict(envelope=str(path),ankles={n:v for n,v in bounds.items() if n.startswith('ankle')},per_tert={k:{n:v for n,v in values.items() if n.startswith('per_tert')} for k,values in roles.items()}),indent=2))

def assess(folder,spec_path):
    spec=json.loads(spec_path.read_text());z=np.load(folder/'controller.npz',allow_pickle=True);t=z['t'];keep=(t>=spec['window_neural_s'][0])&(t<spec['window_neural_s'][1]);rows={};dt=float(np.median(np.diff(t))) if len(t)>1 else .002;covered=bool(len(t)>1 and t[0]<=spec['window_neural_s'][0] and t[-1]>=spec['window_neural_s'][1]-dt-1e-10);accepted=bool(np.any(keep) and covered);cost=0.
    if not np.any(keep):
        result=dict(accepted_reference_envelope=False,cost=None,specification_sha256=sha(spec_path),scope=spec['scope'],joints={},failure='No recorded states in the declared scoring window')
        (folder/'BIOMECHANICAL_ENVELOPE.json').write_text(json.dumps(result,indent=2));return result
    for name,row in spec['joints'].items():
        q=z['qfull'][keep,row['qpos_address']];lo,hi=row['native_bounds_rad'];low=np.maximum(lo-q,0.);high=np.maximum(q-hi,0.);bad=(low+high)>spec['numeric_epsilon_rad'];accepted=accepted and not np.any(bad);excess=np.degrees(low+high)
        rows[name]=dict(observed_native_deg=np.degrees([q.min(),q.max()]).tolist(),accepted=bool(not np.any(bad)),fraction_outside=float(bad.mean()),maximum_excess_deg=float(excess.max()),rms_excess_deg=float(np.sqrt(np.mean(excess**2))))
        cost+=float(np.mean(excess**2))+.1*float(excess.max())**2
    result=dict(accepted_reference_envelope=bool(accepted),complete_scoring_window=covered,cost=cost,specification_sha256=sha(spec_path),scope=spec['scope'],joints=rows)
    (folder/'BIOMECHANICAL_ENVELOPE.json').write_text(json.dumps(result,indent=2));return result

def phase_diagnostics(folder,spec_path):
    """Classify reference-envelope excursions by measured load and RG phase."""
    spec=json.loads(spec_path.read_text());z=np.load(folder/'controller.npz');a=np.load(folder/'controller.assistance.npz');r=json.loads((folder/'controller.result.json').read_text());cols=a['columns'].tolist();samples=a['samples'];t=z['t'];keep=(t>=5)&(t<13);bw=r['measured']['body_weight_N'];result={}
    if (folder/'controller.neurons.npz').exists():v=np.load(folder/'controller.neurons.npz');cells=v['cells'].tolist();phase_source='Recorded RG_E and RG_F voltages'
    else:v=np.load(folder/'source_ib_signals.npz');phase_source='Recorded clipped RG_E-RG_F permission: equality is included in not-extensor-dominant'
    for side,label in [('r','right'),('l','left')]:
        if 'voltage' in v:e=np.interp(t,v['t'],v['voltage'][:,cells.index('RG_E_'+side)])-np.interp(t,v['t'],v['voltage'][:,cells.index('RG_F_'+side)])
        else:e=np.interp(t,v['samples'][:,0],v['samples'][:,2+(side=='l')])
        load=np.interp(t,samples[:,0],samples[:,cols.index(label+'_vertical_GRF_N')])/bw
        masks=dict(loaded=keep&(load>=.025),unloaded=keep&(load<.025),RG_extensor_dominant=keep&(e>0.),RG_not_extensor_dominant=keep&(e<=0.))
        for joint in ('knee_angle','ankle_angle'):
            name=joint+'_'+side;row=spec['joints'][name];q=z['qfull'][:,row['qpos_address']];lo,hi=row['native_bounds_rad'];bad=(q<lo-1e-10)|(q>hi+1e-10);states={}
            for state,mask in masks.items():
                if np.any(mask):states[state]=dict(samples=int(mask.sum()),native_angle_deg=np.degrees([q[mask].min(),q[mask].max()]).tolist(),outside_fraction=float(bad[mask].mean()))
            result[name]=states
    record=dict(scope='Descriptive phase/load split, not evidence of causal necessity.',phase_source=phase_source,joint_states=result);(folder/'PHASE_ENVELOPE_DIAGNOSTICS.json').write_text(json.dumps(record,indent=2));return record

def joint_phase_diagnostics(folder,spec_path):
    """Describe intraleg coordination against the saved reference, without a new gate."""
    from scipy.signal import find_peaks
    spec=json.loads(spec_path.read_text());ref=np.load(spec['source_target']);z=np.load(folder/'controller.npz');t=z['t'];keep=(t>=5)&(t<13);dt=float(np.median(np.diff(t)));result={}
    for side in ('r','l'):
        h=spec['joints']['hip_flexion_'+side]['qpos_address'];k=spec['joints']['knee_angle_'+side]['qpos_address'];r=ref['qpos'][:,[h,k]];q=z['qfull'][:,[h,k]];scale=np.ptp(r,axis=0)
        distance=np.sqrt(np.min(np.sum(((q[keep,None,:]-r[None,:,:])/scale[None,None,:])**2,axis=2),axis=1));peak,_=find_peaks(-q[:,1],prominence=np.radians(10.),distance=max(1,round(.4/dt)));peak=peak[keep[peak]];reference_peak=int(np.argmin(r[:,1]))
        result[side]=dict(reference_at_maximum_knee_flexion=dict(native_knee_deg=float(np.degrees(r[reference_peak,1])),anatomical_hip_flexion_deg=float(-np.degrees(r[reference_peak,0]))),observed_knee_flexion_peaks=[dict(logical_time_s=float(t[i]),native_knee_deg=float(np.degrees(q[i,1])),anatomical_hip_flexion_deg=float(-np.degrees(q[i,0]))) for i in peak],nearest_reference_hip_knee_curve_normalized_distance=dict(mean=float(distance.mean()),p95=float(np.percentile(distance,95)),maximum=float(distance.max())))
    record=dict(scope='Read-only intraleg coordination diagnostic. Positive anatomical hip flexion uses the saved native-to-OpenSim sign; knee flexion remains negative native coordinates. Nearest-reference distance is normalized by each coordinate range. No acceptance threshold or phase warp is added.',reference=str(spec['source_target']),legs=result);(folder/'INTRALEG_REFERENCE_DIAGNOSTICS.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2));return record

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--make-specs',type=Path);p.add_argument('--folder',type=Path);p.add_argument('--spec',type=Path);p.add_argument('--joint-phase',action='store_true');a=p.parse_args()
    if a.make_specs:make_specs(a.make_specs)
    elif a.folder and a.spec:
        if a.joint_phase:joint_phase_diagnostics(a.folder,a.spec)
        else:
            print(json.dumps(assess(a.folder,a.spec),indent=2))
            if (a.folder/'controller.neurons.npz').exists() or (a.folder/'source_ib_signals.npz').exists():print(json.dumps(phase_diagnostics(a.folder,a.spec),indent=2))
    else:p.error('Provide --make-specs, or --folder and --spec')
