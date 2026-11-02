"""Opt-in realization of runner.patch_xml's existing damping specifications.

Materialize inherited baseline attributes before the original surgical patch;
do not duplicate its damping/spring tables or change its default behavior.
"""
from __future__ import annotations
import argparse
from contextlib import contextmanager
import hashlib
import json
import numbers
from pathlib import Path
import re
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
RIG_NAMES=tuple('pelvis_'+s for s in ('tx','ty','tz','tilt','list','rotation'))+('lumbar_extension','lumbar_bending','lumbar_rotation','hip_rotation_r','hip_rotation_l')


def default_policy():
    return {'schema':'aarl_materialize_inherited_damping_v1','materialize_inherited_damping':True}


def materialize_damping(xml_text,model,policy=None):
    if policy is None:return xml_text,{'enabled':False,'added_damping_attributes':0}
    if policy!=default_policy():raise ValueError('Unknown or unused damping-repair policy fields')
    section=re.search(r'<worldbody\b[^>]*>.*?</worldbody>',xml_text,re.S)
    if section is None:raise ValueError('No unique worldbody section')
    names=set();rows=[]
    def replace(match):
        tag=match.group(0);name_match=re.search(r'\bname="([^"]+)"',tag)
        if name_match is None:raise ValueError('Named joints required for explicit damping binding')
        name=name_match.group(1)
        if name in names:raise ValueError('Duplicate named joint')
        names.add(name);joint=model.joint(name);jid=joint.id
        count=6 if int(model.jnt_type[jid])==0 else 3 if int(model.jnt_type[jid])==1 else 1
        values=model.dof_damping[int(joint.dofadr[0]):int(joint.dofadr[0])+count]
        if not np.isfinite(values).all() or not np.all(values==values[0]):raise ValueError('Scalar XML damping cannot represent this compiled baseline')
        existing=bool(re.search(r'\bdamping\s*=',tag))
        if not existing:tag=tag[:-2]+' damping="'+format(float(values[0]),'.17g')+'"/>'
        rows.append({'joint':name,'baseline_damping':float(values[0]),'attribute_added':not existing})
        return tag
    body=re.sub(r'<joint\b[^>]*/>',replace,section.group(0))
    if len(names)!=model.njnt:raise ValueError('Not every compiled joint was matched; unsupported freejoint or XML structure')
    return xml_text[:section.start()]+body+xml_text[section.end():],{'enabled':True,'policy':policy,'added_damping_attributes':sum(r['attribute_added'] for r in rows),'joints':rows}


def patch_xml(xml_text,model,data,*args,policy=None,original_patch=None,**kwargs):
    """Return original patch output unless an explicit repair policy is supplied."""
    if original_patch is None:
        import runner
        original_patch=runner.patch_xml
    prepared,_=materialize_damping(xml_text,model,policy)
    return original_patch(prepared,model,data,*args,**kwargs)


@contextmanager
def enabled(policy=None):
    """Scoped, single-process setup hook for existing runner.apply_harness callers."""
    if policy is None:
        yield {'enabled':False,'calls':[]}
        return
    if policy!=default_policy():raise ValueError('Unknown damping-repair policy')
    import runner
    original=runner.patch_xml
    if getattr(original,'_damping_repair',False):raise RuntimeError('Nested damping-repair contexts are unsupported')
    record={'enabled':True,'calls':[]}
    def wrapped(xml_text,model,data,*args,**kwargs):
        prepared,ledger=materialize_damping(xml_text,model,policy);record['calls'].append(ledger)
        return original(prepared,model,data,*args,**kwargs)
    wrapped._damping_repair=True
    runner.patch_xml=wrapped
    try:yield record
    finally:
        if runner.patch_xml is not wrapped:raise RuntimeError('Another callback replaced the active patch; refusing to overwrite it')
        runner.patch_xml=original


def model_fields(model,exclude_damping=True):
    result={}
    def take(obj,prefix):
        for name in dir(obj):
            if name.startswith('_') or (exclude_damping and prefix=='' and name=='dof_damping'):continue
            value=getattr(obj,name)
            if isinstance(value,np.ndarray):
                result[prefix+name]={'dtype':str(value.dtype),'shape':list(value.shape),'sha256':hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()}
            elif isinstance(value,numbers.Number):result[prefix+name]=float(value)
            elif isinstance(value,(str,bytes)):result[prefix+name]=str(value)
    take(model,'');take(model.opt,'opt.');take(model.stat,'stat.')
    return result


def acceleration_jacobians(model,report):
    """Include dM/dq times nonzero affine acceleration; no equilibrium assumed."""
    import mujoco
    import gait_mechanical_scaling as G
    arms=G.CoupledArms(model);d=mujoco.MjData(model);full=np.zeros((model.nv,model.nv));rows=[]
    def acceleration(pose):
        d.qpos[:]=pose;d.qvel[:]=0.;d.act[:]=0.;d.ctrl[:]=0.;arms.project(d);mujoco.mj_forward(model,d);arms.update(d)
        E=arms.E;mujoco.mj_fullM(model,full,d.qM);M=E.T@full@E;r=E.T@(d.qfrc_bias-d.qfrc_passive-d.qfrc_actuator)
        return np.linalg.solve(M,-r)
    for sample in report['samples']:
        pose=model.qpos0.copy()
        for name,value in sample['fixed_pelvis_positions_native'].items():pose[model.joint(name).qposadr[0]]=value
        pose[arms.addresses]=np.radians([sample['internal_angles_deg'][n] for n in G.JOINTS]);a0=acceleration(pose);J=np.zeros((17,17))
        for col,address in enumerate(arms.addresses):
            minus=pose.copy();plus=pose.copy();minus[address]-=1e-5;plus[address]+=1e-5
            J[:,col]=(acceleration(plus)-acceleration(minus))/2e-5
        rows.append({'position_acceleration_jacobian':J.tolist(),'affine_acceleration_rad_s2':a0.tolist(),'scope':'d[-M(q)^-1*r(q)]/dq with equality reduction; includes inertia derivative at non-equilibrium posture'})
    return rows


def modal_audit(report,local_jacobians=None):
    result=[]
    for index,row in enumerate(report['samples']):
        M=np.array(row['mass_kg_m2']);D=np.array(row['joint_damping_Nm_s_rad']);K=np.array(row['total_zero_activation_restoring_stiffness_Nm_rad']);n=len(M)
        A=np.block([[np.zeros((n,n)),np.eye(n)],[-np.linalg.solve(M,K),-np.linalg.solve(M,D)]])
        if local_jacobians is not None:A[n:,:n]=local_jacobians[index]['position_acceleration_jacobian']
        eig=np.linalg.eigvals(A);eig=sorted(eig,key=lambda x:(float(x.real),float(x.imag)))
        result.append({'phase_percent':row['phase_percent'],'A_per_s':A.tolist(),'eigenvalues_per_s':[[float(v.real),float(v.imag)] for v in eig],
                       'maximum_real_eigenvalue_per_s':max(float(v.real) for v in eig),'unstable_poles':int(sum(v.real>1e-7 for v in eig)),
                       'oscillatory_modes':[{'real_per_s':float(v.real),'imaginary_per_s':float(v.imag),'oscillation_hz':float(v.imag/(2*np.pi)),'spectral_damping_ratio':float(-v.real/abs(v))} for v in eig if v.imag>1e-7]})
    return result


def rig_assistance(model,data,baseline_model):
    """Signed named-coordinate forces/powers; no stepping or state modification."""
    import mujoco
    up=-model.opt.gravity/np.linalg.norm(model.opt.gravity);moving_mass=float(model.body_mass[model.body_weldid!=0].sum());force=np.zeros(3);rows=[]
    for name in RIG_NAMES:
        j=model.joint(name);b=baseline_model.joint(name);q=int(j.qposadr[0]);v=int(j.dofadr[0]);bv=int(b.dofadr[0]);bq=int(b.qposadr[0])
        position=float(data.qpos[q]);velocity=float(data.qvel[v]);k=float(model.jnt_stiffness[j.id]);d=float(model.dof_damping[v]);k0=float(baseline_model.jnt_stiffness[b.id]);d0=float(baseline_model.dof_damping[bv]);reference=float(model.qpos_spring[q]);reference0=float(baseline_model.qpos_spring[bq])
        spring=-k*(position-reference);damper=-d*velocity;source_spring=-k0*(position-reference0);source_damper=-d0*velocity;extra=spring+damper-source_spring-source_damper
        slide=int(model.jnt_type[j.id])==int(mujoco.mjtJoint.mjJNT_SLIDE)
        if slide:force+=extra*data.xaxis[j.id]
        rows.append({'joint':name,'coordinate_unit':'m' if slide else 'rad','generalized_force_unit':'N' if slide else 'Nm','q':position,'qvel':velocity,
                     'spring_force':spring,'damper_force':damper,'spring_power_W':spring*velocity,'damper_power_W':damper*velocity,
                     'source_baseline_force':source_spring+source_damper,'rig_increment_force':extra,'rig_increment_power_W':extra*velocity,
                     'incremental_spring_energy_J':.5*k*(position-reference)**2-.5*k0*(position-reference0)**2})
    return {'rows':rows,'incremental_root_translation_force_world_N':force.tolist(),'incremental_vertical_force_N':float(force@up),
            'incremental_vertical_BW':float(force@up/(moving_mass*np.linalg.norm(model.opt.gravity))),
            'incremental_all_rig_power_W':sum(r['rig_increment_power_W'] for r in rows),
            'incremental_positive_power_W':sum(max(0,r['rig_increment_power_W']) for r in rows),'incremental_negative_power_W':sum(min(0,r['rig_increment_power_W']) for r in rows),
            'normalization_moving_mass_kg':moving_mass,'scope':'All11 declared rig coordinates,including lumbar/hiprotation. Increment is relative to explicit raw baseline coefficients; vertical BW alone omits rotational assistance.'}


def probe(output_dir):
    sys.path[:0]=[str(HERE),str(HERE.parent),str(HERE.parent/'w2l_cpg')]
    import gait_mechanical_scaling as G
    import runner as R
    import mujoco
    inputs=(Path(__file__),HERE/'gait_mechanical_scaling.py',G.MODEL,G.REFERENCE,G.OSIM,G.FORCE,G.SPINAL/'runner.py')
    frozen={str(p):G.sha(p) for p in inputs}
    raw=mujoco.MjModel.from_xml_path(str(G.MODEL));data=mujoco.MjData(raw);R.apply_start_pose(raw,data);text=G.MODEL.read_text(encoding='utf-8');original=R.patch_xml
    disabled=patch_xml(text,raw,data,policy=None)
    if disabled!=original(text,raw,data):raise AssertionError('Disabled policy is not XML-identical')
    with enabled(None):
        if R.patch_xml is not original:raise AssertionError('Disabled context installed callback')
    try:
        with enabled(default_policy()):raise LookupError('Deliberate callback-cleanup fixture')
    except LookupError:pass
    if R.patch_xml is not original:raise AssertionError('Exception did not restore original patch callback')
    try:materialize_damping(text,raw,{**default_policy(),'unused':True});unused_rejected=False
    except ValueError:unused_rejected=True
    if not unused_rejected:raise AssertionError('Unused policy was accepted')
    old,old_construction=G.load_variant('runtime_patch_v3')
    old_report=G.audit('runtime_patch_v3')
    with enabled(default_policy()) as ledger:
        new,new_construction=G.load_variant('runtime_patch_v3');new_report=G.audit('runtime_patch_v3')
    if R.patch_xml is not original:raise AssertionError('Scoped hook was not restored')
    a=model_fields(old);b=model_fields(new);different=[k for k in set(a)|set(b) if a.get(k)!=b.get(k)]
    if different:raise AssertionError('Damping repair changed non-damping model fields: '+str(different))
    for before,after in zip(old_report['samples'],new_report['samples']):
        for key in ('mass_kg_m2','gravity_stiffness_Nm_rad','total_zero_activation_restoring_stiffness_Nm_rad'):
            if not np.array_equal(before[key],after[key]):raise AssertionError('Static non-damping mechanics changed: '+key)
    changed=[]
    for j in range(new.njnt):
        name=mujoco.mj_id2name(new,mujoco.mjtObj.mjOBJ_JOINT,j);v=int(new.jnt_dofadr[j]);before=float(old.dof_damping[v]);after=float(new.dof_damping[v])
        if before!=after:changed.append({'joint':name,'legacy_damping':before,'repaired_damping':after,'role':'RIG' if name in RIG_NAMES else 'declared passive restraint surrogate'})
    # All intended internal branch values are tested,not copied into repair logic.
    for name in ('lumbar_extension','lumbar_bending','lumbar_rotation','hip_rotation_r','hip_rotation_l'):
        if new.dof_damping[new.joint(name).dofadr[0]]!=25.:raise AssertionError('Intended rig damping was not consumed')
    for side in ('r','l'):
        for base in ('ankle_angle','subtalar_angle','mtp_angle','hip_adduction'):
            if new.dof_damping[new.joint(base+'_'+side).dofadr[0]]!=3.:raise AssertionError('Intended restraint damping was not consumed')
    # Static log fixture: chosen small q/v offsets exercise every named rig,
    # but no body integration or physically achieved state is claimed.
    fixture=mujoco.MjData(new);fixture.qpos[:]=new.qpos_spring
    for name in RIG_NAMES:
        j=new.joint(name);fixture.qpos[j.qposadr[0]]+=.001;fixture.qvel[j.dofadr[0]]=.01
    mujoco.mj_forward(new,fixture);rig=rig_assistance(new,fixture,raw)
    if any(r['damper_power_W']>0 for r in rig['rows']):raise AssertionError('Viscous damping injects energy')
    zero_vertical=mujoco.MjData(new);zero_vertical.qpos[:]=new.qpos_spring
    for name in RIG_NAMES[6:]:
        j=new.joint(name);zero_vertical.qpos[j.qposadr[0]]+=.001;zero_vertical.qvel[j.dofadr[0]]=.01
    mujoco.mj_forward(new,zero_vertical);zero_vertical_rig=rig_assistance(new,zero_vertical,raw)
    if zero_vertical_rig['incremental_vertical_BW']!=0 or zero_vertical_rig['incremental_all_rig_power_W']==0:raise AssertionError('Internal rig assistance was lost at zero vertical BW')
    if frozen!={str(p):G.sha(p) for p in inputs}:raise RuntimeError('Audit input changed')
    local_jacobians=acceleration_jacobians(old,old_report)
    result={'schema':'aarl_damping_repair_static_v1','source_sha256':frozen,'policy':default_policy(),'default_off_XML_exact':True,'context_restored':True,
            'non_damping_compiled_fields_compared':len(a),'all_non_damping_compiled_fields_exact':True,'static_M_K_exact_all_four_poses':True,'unused_policy_rejected':unused_rejected,'exception_context_cleanup_verified':True,'non_damping_field_hashes':a,'changed_damping':changed,
            'materialization':ledger['calls'][0],'legacy_construction':old_construction,'repaired_construction':new_construction,
            'legacy_frozen_inertia_modal':modal_audit(old_report),'repaired_frozen_inertia_modal':modal_audit(new_report),
            'non_equilibrium_acceleration_jacobians':local_jacobians,'legacy_modal':modal_audit(old_report,local_jacobians),'repaired_modal':modal_audit(new_report,local_jacobians),
            'rig_logging_fixture':rig,'zero_vertical_internal_rig_fixture':zero_vertical_rig,
            'body_trials':0,'Hunt_training':False,'curriculum':False,'global_release':False,
            'scope':'Versioned optional damping-realization candidate. Four static non-equilibrium old walking postures,fixed pelvis,zero activation,native v3 static surrogate. Both frozen-inertia impedance poles and true reduced acceleration tangents including dM/dq*a0 are saved. These are not achieved-gait or equilibrium-stability proof; contacts,soft equality compliance,and activation/fiber states are excluded.'}
    output_dir.mkdir(parents=True,exist_ok=True)
    (output_dir/'DAMPING_REPAIR_STATIC.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    for label,report in (('legacy',old_report),('repaired',new_report)):
        (output_dir/('DAMPING_REPAIR_MECHANICS_'+label+'.json')).write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({'non_damping_fields_exact':len(a),'changed_damping_joints':len(changed),'legacy_max_real':[s['maximum_real_eigenvalue_per_s'] for s in result['legacy_modal']],'repaired_max_real':[s['maximum_real_eigenvalue_per_s'] for s in result['repaired_modal']],'output':str(output_dir)}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);args=p.parse_args();probe(args.output_dir)
