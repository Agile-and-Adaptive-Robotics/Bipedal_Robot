"""Conservative offline whole-foot clearance/load gate from one aligned capture.

Uses exact archived MJB and recorded qpos. Calls mj_kinematics only: no physics
integration, constraint solve, force recomputation, controller step, or pose fit.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import mujoco


def runs(mask):
    boundaries=np.r_[0,np.flatnonzero(mask[1:]!=mask[:-1])+1,len(mask)]
    return[(int(a),int(b))for a,b in zip(boundaries[:-1],boundaries[1:])if mask[a]]


def evaluate(source,output,start=5.,end=13.,clearance=.001,peak=.01,duration=.08,unload_bw=.025):
    from ground_contact_check import check_samples
    if(start,end)!=(5.,13.):raise ValueError('The preserved independent contact protocol scores5–13s; use that same window for aligned clearance')
    if not np.isfinite([clearance,peak,duration,unload_bw]).all()or not(0<clearance<peak and duration>0 and 0<unload_bw<1):raise ValueError('Invalid declared clearance/load thresholds')
    source=Path(source).resolve();source=source/'aligned_trace.npz'if source.is_dir()else source
    meta=json.loads((source.parent/'ALIGNED_TRACE_METADATA.json').read_text())
    modelpath=source.parent/meta['model_file']
    if meta['schema']!='aarl_aligned_verification_v1':raise ValueError('Actual aligned capture required')
    if mujoco.__version__!=meta['mujoco_version']:raise ValueError('Use the recorded MuJoCo version')
    if hashlib.sha256(modelpath.read_bytes()).hexdigest()!=meta['model_sha256']:raise ValueError('Exact model hash mismatch')
    z=np.load(source);q=z['qpos'];t=z['t_state'];begin=z['force_interval_start'];finish=z['force_interval_end'];clock=z['physics_clock_s'];n=len(t)
    if n<2 or not np.isfinite(q).all()or not np.all(np.diff(t)>0):raise ValueError('Invalid/truncated/nonfinite state trace')
    if not np.array_equal(finish,t):raise ValueError('Forces are not labelled by recorded state interval end')
    pair_ok=np.zeros(n,bool);pair_ok[1:]=np.isclose(begin[1:],t[:-1],rtol=0,atol=1e-9)&np.isclose(np.diff(clock),np.diff(t),rtol=0,atol=1e-8)
    fixed_step=bool(np.allclose(np.diff(t),.002,rtol=0,atol=1e-9))
    if not fixed_step:raise ValueError('Existing independent contact gate requires the actual 2 ms protocol')
    m=mujoco.MjModel.from_binary_path(str(modelpath));data=mujoco.MjData(m)
    if q.shape[1]!=m.nq:raise ValueError('Saved joint coordinates do not match the archived model')
    plane=m.geom('ground-plane').id
    if int(m.geom_type[plane])!=int(mujoco.mjtGeom.mjGEOM_PLANE):raise ValueError('Expected actual ground-plane geometry')
    owners={m.body(body+'_'+side).id:si for si,side in enumerate(('r','l'))for body in('talus','calcn','toes')}
    geometry=[[],[]];manifest=[]
    for g in range(m.ngeom):
        if int(m.geom_type[g])!=int(mujoco.mjtGeom.mjGEOM_MESH):continue
        body=int(m.geom_bodyid[g])
        while body and body not in owners:body=int(m.body_parentid[body])
        if not body:continue
        mesh=int(m.geom_dataid[g]);a=int(m.mesh_vertadr[mesh]);count=int(m.mesh_vertnum[mesh]);vertices=m.mesh_vert[a:a+count].copy()
        geometry[owners[body]].append((g,vertices))
        manifest.append(dict(side=('r','l')[owners[body]],geom_index=g,geom_name=m.geom(g).name,vertices=count,
                             mesh_index=mesh,includes_visual_noncolliding_meshes=True))
    if any(not shapes for shapes in geometry):raise ValueError('Actual bilateral whole-foot meshes are required')
    heights=np.empty((n,2));normals=np.empty((n,3));up=-np.asarray(m.opt.gravity)
    if np.linalg.norm(up)==0:up=np.asarray(meta['native_to_anatomical'])[1]
    for i in range(n):
        data.qpos[:]=q[i];mujoco.mj_kinematics(m,data)
        normal=data.geom_xmat[plane].reshape(3,3)[:,2].copy()
        if normal@up<0:normal=-normal
        normals[i]=normal
        for side,shapes in enumerate(geometry):
            heights[i,side]=min(float(np.min(v@(data.geom_xmat[g].reshape(3,3).T@normal))+(data.geom_xpos[g]-data.geom_xpos[plane])@normal)for g,v in shapes)
    bw=float(meta['body_weight_N']);force=z['foot_force_anatomical_N'][:,:,1]/bw
    if not np.isfinite(force).all():raise ValueError('Nonfinite actual foot forces cannot qualify')
    conservative=np.full_like(heights,np.nan);conservative[1:]=np.minimum(heights[1:],heights[:-1])
    force_small=abs(force)<unload_bw
    force_pair=np.zeros_like(force_small);force_pair[1:]=force_small[1:]&force_small[:-1]
    scored=(begin>=start)&(finish<end)&pair_ok
    qualified=(conservative>clearance)&force_pair&scored[:,None]
    feet=[]
    for side in(0,1):
        episodes=[]
        for a,b in runs(qualified[:,side]):
            elapsed=float(finish[b-1]-begin[a]);maximum=float(conservative[a:b,side].max())
            if elapsed+1e-12>=duration and maximum>=peak:
                episodes.append(dict(interval_start_s=float(begin[a]),interval_end_s=float(finish[b-1]),duration_s=elapsed,
                                     conservative_peak_clearance_m=maximum,maximum_abs_load_bw=float(abs(force[max(0,a-1):b,side]).max()),intervals=b-a))
        feet.append(dict(side=('r','l')[side],sustained_clear_unloaded_swings=len(episodes),episodes=episodes,
                         minimum_mesh_height_m=float(heights[scored,side].min())if np.any(scored)else None,
                         maximum_mesh_height_m=float(heights[scored,side].max())if np.any(scored)else None))
    samples=np.zeros((n,9));samples[:,0]=finish;samples[:,1]=z['harness_wrench_anatomical'][:,1];samples[:,2]=np.nan_to_num(z['explicit_support_anatomical_N'][:,1]);samples[:,7:9]=force*bw
    contact=check_samples(samples,bw)
    knee=[]
    for side in('r','l'):
        maximum=float(np.degrees(q[:,int(m.joint('knee_angle_'+side).qposadr[0])]).max())
        knee.append(dict(side=side,maximum_native_coordinate_deg=maximum,hyperextension=maximum>0.))
    stage=meta.get('provenance',{}).get('configuration',{}).get('stage');ground_protocol=stage in(4,5)
    clearance_gate=all(f['sustained_clear_unloaded_swings']>=3 for f in feet)
    timing_gate=bool(np.all(pair_ok[1:]));knee_gate=not any(k['hyperextension']for k in knee)
    accepted=bool(ground_protocol and timing_gate and clearance_gate and contact['accepted']and knee_gate)
    original=None;configuration=meta.get('provenance',{}).get('configuration_file')
    if configuration:
        resultpath=Path(configuration).with_suffix('.result.json')
        if resultpath.exists():original=dict(path=str(resultpath),assessment=json.loads(resultpath.read_text())['assessment'])
    result=dict(machine='easteregg2',source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                model=str(modelpath),model_sha256=meta['model_sha256'],samples=n,stage=stage,
                accepted=accepted,applicable_to_supported_ground_stepping=ground_protocol,
                timing_gate=timing_gate,protected_knee_gate=knee_gate,knees=knee,clearance_gate=clearance_gate,
                feet=feet,independent_once_counted_contact_gate=contact,original_assessment_preserved=original,
                protocol=dict(window_s=[start,end],unload_bw=unload_bw,minimum_clearance_m=clearance,peak_clearance_m=peak,
                              minimum_duration_s=duration,minimum_swings_per_side=3,geometry_sample_interval_s=.002),
                geometry=manifest,frame=meta['display_frame'],body_weight_N=bw,
                timing_scope='Final qpos at each interval end. Accepted-interval forces are the last solver evaluation, not recomputed at the final pose. Both adjacent endpoint clearances AND both adjacent force observations must qualify; no interpolation or one-step shifted pairing.',
                conservatism='Endpoint minimum and unfiltered absolute superior loads prevent a one-sample load/height transition from granting clearance. This is a discrete engineering gate, not a continuous-collision proof.',
                scope='Independent load/whole-foot-mesh/knee check; original movement/frontal/controller criteria remain separately required. Air or standing protocols do not receive an alternating-ground pass.',
                force_recomputed=False,physics_steps=0,neural_steps=0,kinematics_evaluations=n)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output/'ALIGNED_CLEARANCE.npz',t_state=t,force_interval_start=begin,force_interval_end=finish,
                        whole_foot_clearance_m=heights,conservative_interval_clearance_m=conservative,
                        superior_foot_load_bw=force,qualified_unloaded_clear_interval=qualified,temporal_pair_valid=pair_ok,plane_normal_world=normals)
    (output/'ALIGNED_CLEARANCE_CHECK.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    (output/'ALIGNED_CLEARANCE_CHECK.md').write_text(
        '# Independent aligned clearance/load check\n\n'
        f'Source: `{source}`\n\nStage: {stage}; applicable supported-ground protocol: {ground_protocol}.\n\n'
        f'Independent combined gate: **{accepted}**. Clearance gate: {clearance_gate}; once-counted contact gate: {contact["accepted"]}; knee gate: {knee_gate}; timing gate: {timing_gate}.\n\n'
        f'Qualified swings R/L: {feet[0]["sustained_clear_unloaded_swings"]}/{feet[1]["sustained_clear_unloaded_swings"]}.\n\n'
        'All actual talus/calcaneus/toes mesh vertices are evaluated at every saved 2 ms state using the exact compiled model and forward kinematics only. Both interval endpoints must exceed 1 mm and both adjacent unfiltered superior-force observations must be below 0.025 BW in magnitude, sustained 80 ms with a conservative 10 mm peak. Three episodes per foot are required. Existing once-counted load/alternation and protected-knee gates are additional.\n\n'
        'The last-solver-stage force timing is explicit. No force is recomputed, no new dynamics is integrated, and no source assessment is overwritten. This remains a discrete conservative measurement, not a claim of full gait control or continuous collision clearance.\n',encoding='utf-8')
    print(json.dumps(dict(accepted=accepted,ground_protocol=ground_protocol,swings=[f['sustained_clear_unloaded_swings']for f in feet],
                          timing_gate=timing_gate,knee_gate=knee_gate,physics_steps=0,output=str(output))))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True);p.add_argument('--output',required=True)
    p.add_argument('--start',type=float,default=5.);p.add_argument('--end',type=float,default=13.)
    a=p.parse_args();evaluate(a.source,a.output,a.start,a.end)
