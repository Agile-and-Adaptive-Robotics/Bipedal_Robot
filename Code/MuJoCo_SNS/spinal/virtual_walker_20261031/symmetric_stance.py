"""Offline symmetric-voltage stance targets; never allocates forces online."""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint, linprog, least_squares, brentq
from import_virtual_walker import mujoco
from native_kernel import NativeController


def solve_target(model,data,qpos,bank,params,dofs,feet,record,base,minimum_voltage=-.1):
    """Tie six functional motor voltages; solve exact clipped static Hill laws.

    Binary variables represent the source law's zero-tension clipping. The
    source maximum-force cap is enforced by an activation upper bound, rather
    than approximating a saturated muscle as linear. Contact forces obey an
    inscribed friction diamond, and each foot carries exactly half body weight.
    """
    source=Path(record.get('source',Path(__file__).resolve().parents[4]/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer_Standalone.asim'))
    if not source.exists():
        source=Path('D:/GitHub/Bipedal_Robot/Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer_Standalone.asim')
    assert hashlib.sha256(source.read_bytes()).hexdigest()==record['source_sha256']
    core=NativeController(source,bank)
    roles={'hip_flexor':'Hip MN flx','hip_extensor':'Hip MN ext',
           'knee_flexor':'Knee MN flx','knee_extensor':'Knee MN ext',
           'ankle_dorsiflexor':'Ank MN ext','ankle_plantarflexor':'Ank MN flx'}
    act=[i for i,b in enumerate(bank) if b['type']=='LinearHillMuscle']
    bindings={}
    for pre,idx,coeff in core.outputs:
        if idx not in act:continue
        assert np.array_equal(coeff,[0.,0.,1.,0.]),'Explicit inverse adapter needed'
        bindings[core.names[pre]]=idx
    role_index={};group=[];voltage_pairs={}
    for j,(role,tail) in enumerate(roles.items()):
        pair=[bindings[s+' '+tail] for s in 'LR'];voltage_pairs[role]=pair
        assert np.array_equal(params[pair[0],[6,8]],params[pair[1],[6,8]]),'Unmatched sigmoid midpoint/slope'
        for i in pair:role_index[i]=j
    assert set(role_index)==set(act)
    data.qpos[:]=qpos;data.qvel[:]=0;mujoco.mj_forward(model,data)
    tids=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_TENDON,'route_'+bank[i]['name']) for i in act]
    moment=-data.ten_J.reshape(model.ntendon,model.nv)[np.ix_(tids,dofs)].T
    jac=np.empty((3,model.nv));angular=jac.copy();contacts=[];points=[];sides=[]
    for k,(geom,corners,body) in enumerate(feet):
        vertices=data.geom_xpos[geom]+corners@data.geom_xmat[geom].reshape(3,3).T
        for point in vertices[vertices[:,1]<=vertices[:,1].min()+.00005]:
            mujoco.mj_jac(model,data,jac,angular,point,body)
            contacts.append(jac[:,dofs].T.copy());points.append(point.tolist());sides.append('L' if k<2 else 'R')
    required=(data.qfrc_bias-data.qfrc_passive)[dofs].copy()
    bw=record['robot_mass_kg']*9.81;nm=len(act);nc=len(contacts)
    # a[6], T[12], contact forces[3*n], source clipping switches[12].
    ti=6;fi=ti+nm;zi=fi+3*nc;n=zi+nm
    lo=np.full(n,-np.inf);hi=np.full(n,np.inf);lo[:6]=0.;hi[:6]=1.
    lo[ti:fi]=0.;lo[zi:]=0.;hi[zi:]=1.;rows=[];lower=[];upper=[];muscles=[]
    for role,pair in voltage_pairs.items():
        j=role_index[pair[0]];A,C=params[pair[0],[6,8]]
        lo[j]=1/(1+np.exp(C*(A-minimum_voltage)));hi[j]=1/(1+np.exp(C*A))
    for k,(i,tid) in enumerate(zip(act,tids)):
        l0,width,kse,kpe,damping,maximum,A,B,C,D,*_=params[i]
        assert D>=0,'Additional active-force clipping requires explicit treatment'
        length=float(data.ten_length[tid]);fl=max(0.,1.-((length-l0)/width)**2)
        offset=(kpe*(length-l0)+D*fl)/(1+kpe/kse);slope=B*fl/(1+kpe/kse);j=role_index[i]
        assert slope>0.,'Muscle outside active force-length support'
        hi[j]=min(hi[j],(maximum-offset)/slope);hi[ti+k]=maximum
        muscles.append(dict(bank_index=i,muscle=bank[i]['name'],role=list(roles)[j],length_m=length,length_factor=fl,offset_N=float(offset),activation_slope_N=float(slope),maximum_N=float(maximum)))
    assert np.all(lo[:6]<=hi[:6])
    for k,x in enumerate(muscles):
        j=role_index[x['bank_index']];b=x['offset_N'];s=x['activation_slope_N'];big=max(x['maximum_N'],abs(b+s*lo[j]),abs(b+s*hi[j]))+1
        row=np.zeros(n);row[j]=s;row[ti+k]=-1;rows.append(row);lower.append(-np.inf);upper.append(-b)
        row=np.zeros(n);row[j]=-s;row[ti+k]=1;row[zi+k]=big;rows.append(row);lower.append(-np.inf);upper.append(b+big)
        row=np.zeros(n);row[ti+k]=1;row[zi+k]=-big;rows.append(row);lower.append(-np.inf);upper.append(0.)
    eq=np.zeros((len(dofs),n));eq[:,ti:fi]=moment;eq[:,fi:zi]=np.column_stack(contacts)
    rows.extend(eq);lower.extend(required);upper.extend(required)
    mu=min(.5,min(float(model.geom_friction[g,0]) for g,_,_ in feet))
    for k in range(nc):
        lo[fi+3*k:fi+3*k+3]=[-bw,0,-bw];hi[fi+3*k:fi+3*k+3]=[bw,bw,bw]
        for sx in (-1,1):
            for sz in (-1,1):
                row=np.zeros(n);row[fi+3*k:fi+3*k+3]=[sx,-mu,sz];rows.append(row);lower.append(-np.inf);upper.append(0.)
    for side in 'LR':
        row=np.zeros(n)
        for k,s in enumerate(sides):
            if s==side:row[fi+3*k+1]=1.
        rows.append(row);lower.append(.5*bw);upper.append(.5*bw)
    cost=np.zeros(n);cost[ti:fi]=[1/max(1,x['activation_slope_N']) for x in muscles];cost[:6]=1e-7
    integ=np.zeros(n);integ[zi:]=1
    result=milp(cost,integrality=integ,bounds=Bounds(lo,hi),constraints=LinearConstraint(np.array(rows),lower,upper),options={'time_limit':20.,'mip_rel_gap':1e-9})
    report=dict(machine='easteregg2',baseline=str(base),source_path=str(source),source_sha256=record['source_sha256'],qpos=np.asarray(qpos).tolist(),feasible=bool(result.success),solver_status=result.message,body_weight_N=bw,friction_inscribed_diamond_mu=mu,balanced_foot_normal_load_BW=.5,contact_points=points,contact_sides=sides,minimum_motor_voltage_V=minimum_voltage,muscles=muscles,scope='Offline shared functional L/R motor-voltage target with exact source static Hill clipping and passive joint forces; no rail/guide reaction, no online force allocation, and no dynamic/neural pass.')
    if not result.success:
        # Diagnose symmetry separately from mechanical and foot-load feasibility.
        independent_A=eq[:,ti:zi];independent_bounds=[]
        for k,item in enumerate(muscles):
            j=role_index[item['bank_index']];b=item['offset_N'];s=item['activation_slope_N']
            independent_bounds.append((max(0,b+s*lo[j]),min(item['maximum_N'],max(0,b+s*hi[j]))))
        independent_bounds.extend(zip(lo[fi:zi],hi[fi:zi]))
        contact_rows=np.array(rows)[3*nm+len(dofs):3*nm+len(dofs)+4*nc,ti:zi]
        balance_rows=np.array(rows)[-2:,ti:zi]
        for label,balanced in [('independent_with_balanced_feet',True),('independent_free_load_split',False)]:
            Ae=np.vstack([independent_A,balance_rows]) if balanced else independent_A
            be=np.r_[required,[.5*bw,.5*bw]] if balanced else required
            trial=linprog(cost[ti:zi],A_eq=Ae,b_eq=be,A_ub=contact_rows,b_ub=np.zeros(len(contact_rows)),bounds=independent_bounds,method='highs')
            report[label]=dict(feasible=bool(trial.success),status=trial.message)
            if trial.success:report[label].update(tensions_N=trial.x[:nm].tolist(),normal_load_N={s:float(trial.x[nm:].reshape(nc,3)[np.array(sides)==s,1].sum()) for s in 'LR'})
        report['functional_torque_map']={item['muscle']:{str(mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_JOINT,int(model.dof_jntid[dof]))):float(moment[row,k]) for row,dof in enumerate(dofs) if abs(moment[row,k])>1e-8} for k,item in enumerate(muscles)}
        # Distinguish exact voltage equality from a shared functional primitive:
        # one nonnegative torque magnitude per role, inverted through each
        # actuator's measured-in-model moment arm and source force law offline.
        na=6+nm+3*nc;Ae=np.zeros((len(dofs)+nm+2,na));be=np.r_[required,np.zeros(nm),[.5*bw,.5*bw]]
        Ae[:len(dofs),6:]=independent_A;torque_bounds=[(0,None)]*6
        for k,item in enumerate(muscles):
            j=role_index[item['bank_index']];arm=float(np.max(abs(moment[:,k])))
            assert np.count_nonzero(abs(moment[:,k])>1e-7)==1,'Expected uniarticular actuator'
            Ae[len(dofs)+k,j]=-1;Ae[len(dofs)+k,6+k]=arm;item['anatomical_moment_arm_m']=arm
        Ae[-2:,6:]=balance_rows;cc=np.r_[np.full(6,1e-7),cost[ti:zi]]
        torque=linprog(cc,A_eq=Ae,b_eq=be,A_ub=np.column_stack([np.zeros((len(contact_rows),6)),contact_rows]),b_ub=np.zeros(len(contact_rows)),bounds=torque_bounds+independent_bounds,method='highs')
        primitive=dict(feasible=bool(torque.success),status=torque.message,scope='Shared six functional torque magnitudes; each source muscle voltage derived offline with its actual moment arm, passive offset, and sigmoid. No independently fitted L/R gain and no online allocator.')
        if torque.success:
            commands=[]
            for k,item in enumerate(muscles):
                i=item['bank_index'];T=float(torque.x[6+k]);activation=max(lo[role_index[i]],(T-item['offset_N'])/item['activation_slope_N']);A,C=params[i,[6,8]];voltage=float(A-np.log(1/activation-1)/C)
                commands.append(dict(muscle=item['muscle'],bank_index=i,role=item['role'],MN_name=next(name for name,b in bindings.items() if b==i),target_tension_N=T,target_motor_voltage_V=voltage,target_activation=activation))
            forces=torque.x[6+nm:].reshape(nc,3)
            primitive.update(shared_role_torque_Nm={role:float(torque.x[j]) for j,role in enumerate(roles)},motor_commands=commands,contact_forces_N=forces.tolist(),max_equilibrium_residual=float(max(abs(Ae@torque.x-be))),normal_load_by_side_N={s:float(forces[np.array(sides)==s,1].sum()) for s in 'LR'})
        report['shared_functional_torque_primitive']=primitive
        return report
    x=result.x;commands={}
    for role,pair in voltage_pairs.items():
        j=role_index[pair[0]];A,C=params[pair[0],[6,8]];activation=float(x[j]);voltage=float(A-np.log(1/activation-1)/C)
        commands[role]=dict(shared_motor_voltage_V=voltage,shared_activation=activation,muscles=[bank[i]['name'] for i in pair],MN_names=[s+' '+roles[role] for s in 'LR'])
    tensions=x[ti:fi];normal=x[fi:zi].reshape(nc,3)
    for k,item in enumerate(muscles):
        item['target_tension_N']=float(tensions[k]);item['raw_equilibrium_tension_N']=float(item['offset_N']+item['activation_slope_N']*x[role_index[item['bank_index']]])
    report.update(shared_commands=commands,contact_forces_N=normal.tolist(),normal_load_by_side_N={s:float(normal[np.array(sides)==s,1].sum()) for s in 'LR'},max_equilibrium_residual=float(np.max(abs(eq@x-required))),exact_source_law_residual_N=float(max(abs(np.clip(v['raw_equilibrium_tension_N'],0,v['maximum_N'])-v['target_tension_N']) for v in muscles)),normalized_tension_objective=float(cost@x),robot_dofs=dofs,required_bias_minus_passive=required.tolist())
    return report


def mechanical_hold(target_path,parameters_path,duration=2.,target_override=None,preserve_supports=False,output_prefix='MECHANICAL_HOLD'):
    """Independent constant-voltage check with actual Hill dynamics and contacts.

    Both cold and equilibrium-preloaded muscle states are explicit conditions.
    Rail/guide contacts are disabled; no body force, pose clamp, or inverse force
    is applied. Only short, decimated traces are retained in one compressed file.
    """
    from native_kernel import hill_step
    from contact_guides import Measurements
    target_path=Path(target_path);target=json.loads(target_path.read_text()) if target_override is None else target_override;primitive=target['shared_functional_torque_primitive'];assert primitive['feasible']
    base=Path(target['baseline']);bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];table=np.load(parameters_path);assert table['muscle_names'].tolist()==[b['name'] for b in bank];params=table['pooled']
    summaries={};traces={}
    for prepared in (False,True):
        model=mujoco.MjModel.from_xml_path(str(base/'model.xml'));data=mujoco.MjData(model);data.qpos[:]=target['qpos'];disabled=[]
        for g in range(model.ngeom):
            name=mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_GEOM,g) or ''
            if not preserve_supports and ('support_rail' in name or 'Guide_' in name):
                model.geom_contype[g]=0;model.geom_conaffinity[g]=0;disabled.append(name)
        assert model.npair==0,'Explicit contact pairs require separate removal'
        mujoco.mj_forward(model,data);root=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Root');measure=Measurements(model,root)
        tids=np.array([mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_TENDON,'route_'+b['name']) for b in bank]);act=[(v['bank_index'],mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_ACTUATOR,v['muscle'])) for v in primitive['motor_commands']]
        vm=np.full(len(bank),-.15)
        for command in primitive['motor_commands']:
            # Zero tension is an upper voltage bound, not a reason to excite an
            # otherwise slack antagonist to its exact recruitment threshold.
            vm[command['bank_index']]=command['target_motor_voltage_V'] if command['target_tension_N']>1e-8 else target['minimum_motor_voltage_V']
        length=data.ten_length[tids].copy();last=length.copy();internal=np.zeros(len(bank));force=internal.copy();ia=internal.copy();ii=internal.copy();ib=internal.copy();velocity=internal.copy();dt=float(model.opt.timestep)
        if prepared:
            for i,p in enumerate(params):
                l0,width,kse,kpe,damping,maximum,A,B,C,D,*_=p;fl=max(0,1-((length[i]-l0)/width)**2)
                internal[i]=(kpe*(length[i]-l0)+max(0,B/(1+np.exp(min(700,C*(A-vm[i]))))+D)*fl)/(1+kpe/kse)
        hill_step(length,velocity,vm,params,internal.copy(),force.copy(),ia.copy(),ii.copy(),ib.copy(),dt)
        rows=[];failure=None;start=time.perf_counter();knees=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'knee_'+s) for s in 'LR'];kmin=np.full(2,np.inf);support=[]
        for k in range(round(duration/dt)):
            mujoco.mj_step1(model,data);length=data.ten_length[tids];velocity[:]=(length-last)/dt;last[:]=length;hill_step(length,velocity,vm,params,internal,force,ia,ii,ib,dt)
            for b,a in act:data.ctrl[a]=force[b]
            mujoco.mj_step2(model,data);mujoco.mj_forward(model,data);measured=measure.sample(data);angles=np.array([10,55])+np.degrees(data.qpos[model.jnt_qposadr[knees]]);kmin=np.minimum(kmin,angles)
            support.append([measured['rail_force_world'][1],measured['guide_force_world'][1]])
            if k%5==0:rows.append((k*dt,data.qpos.copy(),force.copy(),measured['foot_force_world'][:,1].copy(),measured['rail_force_world'].copy(),measured['guide_force_world'].copy()))
            if not np.isfinite(data.qpos).all():failure='nonfinite'
            elif angles.min()<0.:failure='knee_hyperextension'
            elif data.xpos[root,1]<.75:failure='root_height_below_.75m'
            if failure:break
        label='equilibrium_preloaded' if prepared else 'cold_start';foot=np.array([z[3] for z in rows]);rail=np.array([z[4] for z in rows]);guide=np.array([z[5] for z in rows]);bw=target['body_weight_N']
        summaries[label]=dict(completed_s=(k+1)*dt,failure=failure,warmed_wall_s=time.perf_counter()-start,knee_min_deg=kmin.tolist(),root_y_final_m=float(data.xpos[root,1]),root_world_displacement_m=(data.qpos[:3]-np.array(target['qpos'][:3])).tolist(),mean_foot_load_BW=(foot.mean(axis=0)/bw).tolist(),maximum_rail_force_N=float(abs(rail).max()),maximum_guide_force_N=float(abs(guide).max()),disabled_support_geometries=disabled,motor_voltage_V=vm.tolist(),scope='Constant motor voltages and actual dynamic source Hill law; no neural controller, no supporting apparatus contact, stationary belt. Preloaded and cold conditions distinct.')
        summaries[label].update(preserve_source_apparatus=preserve_supports,upward_channels_BW=float(np.maximum(support,0).sum(axis=1).mean()/bw),downward_channels_BW=float(np.maximum(-np.array(support),0).sum(axis=1).mean()/bw),target_provenance=target.get('provenance'))
        summaries[label]['terminal_qvel']=data.qvel.tolist();summaries[label]['steady_state_established']=False
        if preserve_supports:summaries[label]['scope']='Constant source motor voltages with actual dynamic Hill law and unchanged physical guide/rail apparatus. Stationary belt. No neural controller, desired-force replay, or body/trajectory clamp. Cold and equilibrium-preloaded conditions distinct; static LP support is not imposed.'
        traces.update({label+'_t':np.array([z[0] for z in rows]),label+'_q':np.array([z[1] for z in rows]),label+'_force':np.array([z[2] for z in rows]),label+'_foot_N':foot})
    out=target_path.parent;np.savez_compressed(out/(output_prefix+'.npz'),**traces);(out/(output_prefix+'.json')).write_text(json.dumps(summaries,indent=2),encoding='utf-8');return summaries


def fit_neural_recruitment(target_path,transient_protected=False):
    """Fit tied stance chemical weights from actual baseline spike readouts.

    This is a frozen-background current fit, not coupled validation. Original
    transient PF→MN axons stay untouched. A new phase/load-qualified spiking
    stance relay supplies tonic excitatory/inhibitory recruitment through two
    fixed chemical reversals, with each axon bounded by the source PF envelope.
    """
    import inspect
    from spiking_central_runtime import SpikingCentralController as C
    from recruitment import apply_recruitment
    from spiking_sensory import encode_step
    target_path=Path(target_path);target=json.loads(target_path.read_text());primitive=target['shared_functional_torque_primitive'];assert primitive['feasible'];base=Path(target['baseline']);record=json.loads((base/'SPIKING_COMPARISON.json').read_text());bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];z=np.load(base/'trajectory.npz')
    C.configuration=record['configuration'];C.sensory_spiking=True;C.balance_configuration=None;C.clock_scale=1.;C.reference_calculations=False;C.audit_state=False
    allowed=set(inspect.signature(C.__bases__[0]).parameters)-{'source','bank'};cfg=record['source_configuration'];c=C(cfg['source'],bank,**{k:v for k,v in cfg.items() if k in allowed});apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale']);c.apply_output_gain()
    dt=float(record['result']['physics_dt']);bw=target['body_weight_N'];N=len(z['t']);phase=np.empty((N,2))
    for k,side in enumerate('LR'):
        F=c.names.index(side+' RG flx');E=c.names.index(side+' RG ext');f=np.maximum(0,z['neuro'][:,F]-c.params[F,0]);e=np.maximum(0,z['neuro'][:,E]-c.params[E,0]);phase[:,k]=e/(f+e+1e-12)
    load=np.zeros(2);v=np.tile(np.linspace(.01,.99,16),(2,1));ref=np.zeros((2,16));decoded=np.zeros(2);spikes=np.zeros((2,16),dtype=np.int8);counts=np.zeros(2,dtype=np.int64);target_activity=np.zeros(2);readout=np.empty((N,2));tau=np.full(2,.025)
    for k in range(N):
        forces=np.zeros(2) if k==0 else z['foot_vertical'][k-1]/bw;load+=(np.clip((forces-.05)/.20,0,1)-load)*(-np.expm1(-dt/.02));target_activity[:]=.2*load*phase[k];a=0.
        for sub in range(round(dt/c.dt)):
            encode_step(target_activity,v,ref,decoded,spikes,counts,c.dt,200.,.005,tau,.001);a=a+np.clip(decoded*5.,0,1)
        readout[k]=a/round(dt/c.dt)
    settings=dict(pathway_kind='phase_stance_recruitment',input_current_nA=4.,off_load_BW=.05,full_load_BW=.25,load_filter_tau_s=.02,blend=1.,offline_target_sha256=hashlib.sha256(target_path.read_bytes()).hexdigest());fits=[]
    if transient_protected:settings['transient_quiet_PF_threshold']=.2
    for role in primitive['shared_role_torque_Nm']:
        observations=[];caps=[]
        for command in primitive['motor_commands']:
            if command['role']!=role:continue
            name=command['MN_name'];post=c.names.index(name);side=0 if name.startswith('L') else 1;assert c.params[post,3]==0 and c.params[post,2]==1
            gs=np.ones(N);ge=np.full(N,c.params[post,0]+c.tonic[post]);matched=[e for e in c.output_edits if e['target']==name];assert len(matched)==1;caps.append(matched[0]['before_gmax_uS'])
            for e in c.edges[c.edges[:,1]==post]:
                active=e[2]*np.clip((z['neuro'][:,int(e[0])]-e[3])/(e[4]-e[3]),0,1);gs+=active;ge+=active*e[5]
            assert not any(i==post for i,*_ in c.inputs),'Direct physical MN input requires explicit reconstruction'
            activity=readout[:,side].copy()
            if transient_protected:
                original=[e for e in c.edges if int(e[1])==post and ' PF ' in c.names[int(e[0])] and not c.names[int(e[0])].endswith(' IN')];assert len(original)==1;e=original[0];pf_active=np.clip((z['neuro'][:,int(e[0])]-e[3])/(e[4]-e[3]),0,1);activity*=np.clip((.2-pf_active)/.2,0,1)
            mask=(z['t']>=2)&(z['foot_vertical'][:,side]>.1*bw)&(phase[:,side]>.6)&(activity>.05)
            if mask.sum()<=10:raise ValueError('Insufficient loaded PF-quiet samples: '+json.dumps(dict(role=role,side='LR'[side],samples=int(mask.sum()),threshold=settings.get('transient_quiet_PF_threshold'))))
            observations.append((gs[mask],ge[mask],activity[mask],command,side))
        assert len(observations)==2 and np.isclose(caps[0],caps[1]);cap=min(caps);bounds=np.array([cap,cap*1.5])
        def errors(weights):
            errors=[]
            for gs,ge,activity,command,side in observations:
                pred=(ge+activity*(-10*weights[0]-100*weights[1]))/(gs+activity*weights.sum());expected=1000*command['target_motor_voltage_V'];residual=pred-expected
                if command['target_tension_N']<1e-8:residual=np.maximum(0,residual)
                errors.extend(residual/np.sqrt(max(1,len(residual))))
            return np.r_[errors,.025*weights/np.maximum(bounds,1e-9)]
        fit=least_squares(errors,bounds*.1,bounds=(np.zeros(2),bounds),xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=100)
        if all(command['target_tension_N']<1e-8 and np.all(ge/gs<=1000*command['target_motor_voltage_V']) for gs,ge,activity,command,side in observations):fit.x[:]=0.
        settings[role+'_exc_uS']=float(fit.x[0]);settings[role+'_inh_uS']=float(fit.x[1]);rows=[]
        for gs,ge,activity,command,side in observations:
            before=ge/gs;after=(ge+activity*(-10*fit.x[0]-100*fit.x[1]))/(gs+activity*fit.x.sum());wanted=1000*command['target_motor_voltage_V'];rows.append(dict(side='LR'[side],samples=len(gs),target_tension_N=command['target_tension_N'],target_voltage_mV=wanted,zero_tension_is_voltage_ceiling=command['target_tension_N']<1e-8,before_mean_mV=float(before.mean()),after_mean_mV=float(after.mean()),after_RMSE_mV=float(np.sqrt(np.mean((np.maximum(0,after-wanted) if command['target_tension_N']<1e-8 else after-wanted)**2)))))
        fits.append(dict(role=role,exc_gmax_uS=float(fit.x[0]),inh_gmax_uS=float(fit.x[1]),source_pre_attenuation_PF_gmax_uS=cap,individual_axon_current_envelope_nA=cap*90,rows=rows))
    out=target_path.parent.parent;suffix='_PROTECTED' if transient_protected else '';(out/('STANCE_RECRUITMENT'+suffix+'_CONFIGURATION.json')).write_text(json.dumps(settings,indent=2),encoding='utf-8');report=dict(source_sha256=target['source_sha256'],baseline=str(base),target=str(target_path),configuration=settings,fits=fits,scope='Frozen actual baseline chemical conductance snapshots; equilibrium voltage approximation ignores transient derivatives and feedback changes. Shared L/R weights, original PF→MN transient axons untouched. New tonic phase/load spiking relay must pass isolated state/current checks, air, and coupled ground tests. Source-cap envelope is pre pf_output_scale attenuation, not measured baseline current.');(out/('OFFLINE_NEURAL'+suffix+'_FIT.json')).write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def audit_dynamic_support(root):
    """Necessary whole-body wrench feasibility of the actual saved gait.

    Muscle forces are internal and cannot repair a missing external root wrench
    while the complete trajectory is held fixed. An outer friction pyramid is
    deliberately optimistic: infeasibility is not caused by using the earlier
    conservative inscribed static-contact approximation. This audit does not
    claim that changing a controller cannot change the gait trajectory.
    """
    from scipy.signal import savgol_filter
    root=Path(root);fit=json.loads((root/'OFFLINE_NEURAL_FIT.json').read_text());base=Path(fit['baseline']);z=np.load(base/'trajectory.npz');record=json.loads((base/'RESULT.json').read_text());model=mujoco.MjModel.from_xml_path(str(base/'model.xml'));data=mujoco.MjData(model);dt=record['physics_dt'];bw=record['robot_mass_kg']*9.81;q=z['q'];velocity=np.zeros((len(q),model.nv))
    for k in range(1,len(q)-1):mujoco.mj_differentiatePos(model,velocity[k],2*dt,q[k-1],q[k+1])
    velocity[0]=velocity[1];velocity[-1]=velocity[-2];velocity=savgol_filter(velocity,31,3,axis=0);acceleration=savgol_filter(velocity,31,3,deriv=1,delta=dt,axis=0);M=np.empty((model.nv,model.nv));robot=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Root');joint=int(model.body_jntadr[robot]);rd=np.arange(model.jnt_dofadr[joint],model.jnt_dofadr[joint]+6);belt=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Treadmill_Belt');bodies={mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,name):name for name in ('foot_L','toe_L','foot_R','toe_R')};jac=np.empty((3,model.nv));angular=jac.copy();rows=[]
    for seconds in np.linspace(2.2,7.8,24):
        k=int(round(seconds/dt));data.qpos[:]=q[k];data.qvel[:]=velocity[k];mujoco.mj_forward(model,data);mujoco.mj_fullM(model,M,data.qM);required=(M@acceleration[k]+data.qfrc_bias-data.qfrc_passive)[rd];foot=[];rails=[];points=[]
        for con in data.contact:
            b1=int(model.geom_bodyid[con.geom1]);b2=int(model.geom_bodyid[con.geom2]);footbody=b1 if b1 in bodies else b2 if b2 in bodies else None
            if footbody is not None and belt in (b1,b2):
                side=0 if bodies[footbody].endswith('L') else 1
                if z['foot_vertical'][k,side]<.05*bw:continue
                mujoco.mj_jac(model,data,jac,angular,con.pos,footbody);foot.append((jac[:,rd].T.copy(),float(con.friction[0])));points.append(con.pos.tolist())
            else:
                names=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_GEOM,int(g)) or '' for g in (con.geom1,con.geom2)]
                if not any('support_rail' in n for n in names):continue
                railfirst='support_rail' in names[0];body=b2 if railfirst else b1;normal=con.frame[:3]*(1 if railfirst else -1)
                if normal[1]<.9:continue
                mujoco.mj_jac(model,data,jac,angular,con.pos,body);rails.append((jac[:,rd].T.copy(),float(con.friction[0])))
        def solve(contacts,nfoot,inner=False):
            if not contacts:return None
            A=np.column_stack([a for a,mu in contacts]);ineq=[];n=len(contacts);cost=np.zeros(3*n);bounds=[(-4*bw,4*bw),(0,4*bw),(-4*bw,4*bw)]*n
            for i,(J,mu) in enumerate(contacts):
                if i>=nfoot:cost[3*i+1]=1/bw
                if inner:
                    for sx in (-1,1):
                        for sz in (-1,1):
                            row=np.zeros(3*n);row[3*i:3*i+3]=[sx,-mu,sz];ineq.append(row)
                else:
                    for axis in (0,2):
                        for sign in (-1,1):
                            row=np.zeros(3*n);row[3*i+axis]=sign;row[3*i+1]=-mu;ineq.append(row)
            return linprog(cost,A_eq=A,b_eq=required,A_ub=np.array(ineq),b_ub=np.zeros(len(ineq)),bounds=bounds,method='highs')
        outer=solve(foot,len(foot));inner=solve(foot,len(foot),True);assisted=solve(foot+rails,len(foot));com=data.subtree_com[robot].copy();actual_up=sum(max(0,float(z[key][k,1])) for key in ('rail_force_world','guide_force_world','harness_force'))/bw;actual_down=sum(max(0,-float(z[key][k,1])) for key in ('rail_force_world','guide_force_world','harness_force'))/bw
        rows.append(dict(t_s=float(z['t'][k]),foot_points=len(foot),rail_points=len(rails),foot_only_outer_feasible=bool(outer is not None and outer.success),foot_only_inner_feasible=bool(inner is not None and inner.success),with_actual_rail_points_feasible=bool(assisted is not None and assisted.success),minimum_rail_upward_BW=float(assisted.fun) if assisted is not None and assisted.success else None,actual_upward_channels_BW=actual_up,actual_downward_channels_BW=actual_down,actual_foot_load_BW=(z['foot_vertical'][k]/bw).tolist(),com_world_m=com.tolist(),contact_world_m=points,required_root_wrench=required.tolist(),root_pitch_deg=float(np.degrees(np.arctan2(data.xmat[robot].reshape(3,3)[0,1],data.xmat[robot].reshape(3,3)[1,1])))))
    report=dict(baseline=str(base),source_sha256=record['source_sha256'],velocity_estimation='Central mj_differentiatePos followed by31sample/62ms cubic Savitzky-Golay smoothing; acceleration is derivative of smoothed velocity.',samples=rows,foot_only_outer_feasible_count=sum(r['foot_only_outer_feasible'] for r in rows),foot_only_inner_feasible_count=sum(r['foot_only_inner_feasible'] for r in rows),with_rail_feasible_count=sum(r['with_actual_rail_points_feasible'] for r in rows),scope='Necessary whole-body force/moment feasibility at24actual gait poses and smoothed accelerations, using actual geometric loaded-foot/rail contact points. Outer friction pyramid optimistically includes the true cone; inner diamond gives sufficient friction margin for root wrench only. Ignores joint-muscle bounds/history, so foot-only success is not full neuromechanical feasibility. Failure means these same poses/accelerations/contact locations need changed kinematics or external support, not that any controller must fail. Not an online allocator.')
    report['interpretation_status']='Preliminary smoothed-kinematics3D diagnostic: apparatus reconstruction incomplete. Do not use this feasibility count as the validated support conclusion.';report['authoritative_followup']=str(root/'MEASURED_WRENCH_SUPPORT_FEASIBILITY.json')
    (root/'DYNAMIC_SUPPORT_FEASIBILITY.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def audit_measured_support(root,baseline=None,output_name='MEASURED_WRENCH_SUPPORT_FEASIBILITY.json'):
    """Sagittal necessary support test using recorded forces, not differentiated q.

    Includes an own-foot-wrench reconstruction control and a bounded moment
    uncertainty from noncoplanar contact heights. This avoids interpreting an
    acceleration-smoothing or omitted-guide mismatch as mechanical failure.
    """
    root=Path(root);fit=json.loads((root/'OFFLINE_NEURAL_FIT.json').read_text());base=Path(fit['baseline'] if baseline is None else baseline);z=np.load(base/'trajectory.npz');r=json.loads((base/'RESULT.json').read_text());m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);rid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root');belt=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Treadmill_Belt');feet={mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,n):i//2 for i,n in enumerate(('foot_L','toe_L','foot_R','toe_R'))};bw=r['robot_mass_kg']*9.81;rows=[]
    for seconds in np.linspace(2.2,7.8,24):
        k=int(round(seconds/r['physics_dt']));d.qpos[:]=z['q'][k];mujoco.mj_forward(m,d);origin=d.xipos[rid].copy();contacts=[];rails=[];side_heights=[[],[]];side_mu=np.zeros(2)
        for con in d.contact:
            b1,b2=[int(m.geom_bodyid[g]) for g in (con.geom1,con.geom2)];body=b1 if b1 in feet else b2 if b2 in feet else None;point=con.pos.copy();mu=float(con.friction[0]);A=np.array([[1.,0],[0,1.],[-(point[1]-origin[1]),point[0]-origin[0]]])
            if body is not None and belt in (b1,b2):contacts.append((A,mu));side_heights[feet[body]].append(point[1]);side_mu[feet[body]]=max(side_mu[feet[body]],mu)
            else:
                names=[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,int(g)) or '' for g in (con.geom1,con.geom2)]
                if any('support_rail' in n for n in names):
                    first='support_rail' in names[0];normal=con.frame[:3]*(1 if first else -1)
                    if normal[1]>.9:rails.append((A,mu))
        foot_wrench=np.zeros(3);uncertainty=.005
        for side in (0,1):
            F=z['foot_force_world'][k,side];cop=z['foot_cop_world'][k,side]
            if np.isfinite(cop).all():
                foot_wrench+=np.r_[F[:2],np.cross(cop-origin,F)[2]]
                if side_heights[side]:uncertainty+=side_mu[side]*max(0,F[1])*max(abs(np.array(side_heights[side])-cop[1]))
        rail_wrench=np.r_[z['rail_force_world'][k,:2],z['rail_torque_about_pelvis'][k,2]];guide_wrench=np.r_[z['guide_force_world'][k,:2],z['guide_torque_about_pelvis'][k,2]];harness=np.r_[z['harness_force'][k,:2],z['harness_torque'][k,2]]
        def solve(target,items,nfoot):
            if not items:return None
            A=np.column_stack([a for a,mu in items]);n=len(items);ineq=[];bounds=[(-100*bw,100*bw),(0,100*bw)]*n;cost=np.zeros(2*n)
            for j,(a,mu) in enumerate(items):
                if j>=nfoot:cost[2*j+1]=1/bw
                for sign in (-1,1):
                    row=np.zeros(2*n);row[2*j:2*j+2]=[sign,-mu];ineq.append(row)
            ineq.extend([A[2],-A[2]]);return linprog(cost,A_eq=A[:2],b_eq=target[:2],A_ub=ineq,b_ub=np.r_[np.zeros(2*n),target[2]+uncertainty,-target[2]+uncertainty],bounds=bounds,method='highs')
        own=solve(foot_wrench,contacts,len(contacts));all_foot=solve(foot_wrench+rail_wrench+guide_wrench+harness,contacts,len(contacts));keep_guide=solve(foot_wrench+rail_wrench,contacts,len(contacts));with_rails=solve(foot_wrench+rail_wrench,contacts+rails,len(contacts));row=dict(t_s=float(z['t'][k]),foot_contact_points=len(contacts),upward_rail_points=len(rails),moment_uncertainty_Nm=float(uncertainty),own_foot_wrench_control_feasible=bool(own is not None and own.success),foot_only_all_support_removed_feasible=bool(all_foot is not None and all_foot.success),foot_only_rails_removed_guides_preserved_feasible=bool(keep_guide is not None and keep_guide.success),rail_reconstruction_control_feasible=bool(with_rails is not None and with_rails.success),minimum_required_upward_rail_BW=float(with_rails.fun) if with_rails is not None and with_rails.success else None,measured_rail_BW=float(rail_wrench[1]/bw),measured_guide_BW=float(guide_wrench[1]/bw),foot_wrench=foot_wrench.tolist(),rail_wrench=rail_wrench.tolist(),guide_wrench=guide_wrench.tolist())
        required=foot_wrench+rail_wrench+guide_wrench+harness;beltg=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_Treadmill_Belt');floor_y=float(d.geom_xpos[beltg,1]+m.geom_size[beltg,1]);signs=np.array(list(__import__('itertools').product((-1,1),repeat=3)));actual_x=[];full_x=[];heels=[]
        for A,mu in contacts:actual_x.append(float(A[2,1]+origin[0]))
        for side_index,side in enumerate('LR'):
            if z['foot_vertical'][k,side_index]<=.05*bw:continue
            for name in ('foot_'+side,'toe_'+side):
                g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+name);R=d.geom_xmat[g].reshape(3,3);vertices=d.geom_xpos[g]+(signs*m.geom_size[g])@R.T;full_x.extend(vertices[:,0].tolist())
                if name.startswith('foot_'):heels.append(float(vertices[signs[:,0]==-np.sign(R[0,0]),1].min()-floor_y))
        if required[1]>1e-6 and actual_x:
            desired=origin[0]+(required[2]+(floor_y-origin[1])*required[0])/required[1];row.update(required_CoP_x_m=float(desired),actual_contact_x_span_m=[min(actual_x),max(actual_x)],required_CoP_gap_m=float(desired-np.clip(desired,min(actual_x),max(actual_x))),required_CoP_uncertainty_m=float(uncertainty/required[1]))
            if full_x:row.update(loaded_foot_full_projection_m=[min(full_x),max(full_x)],required_CoP_inside_full_foot_projection=bool(min(full_x)<=desired<=max(full_x)),loaded_heel_gaps_m=heels)
        rows.append(row)
    valid=[row for row in rows if row['own_foot_wrench_control_feasible'] and row['rail_reconstruction_control_feasible']];report=dict(baseline=str(base),source_sha256=r['source_sha256'],samples=rows,valid_reconstruction_samples=len(valid),all_support_removed_feasible_on_valid=sum(row['foot_only_all_support_removed_feasible'] for row in valid),rails_removed_guides_preserved_feasible_on_valid=sum(row['foot_only_rails_removed_guides_preserved_feasible'] for row in valid),scope='Sagittal necessary external-wrench feasibility at24actual gait poses. Required force/moment is from measured foot/rail/guide reactions, not differentiated kinematics. Foot moment uses measured normal-force CoP, with conservative tangential-height moment uncertainty. Own-foot and rail reconstruction controls must pass before interpreting each frame. Frontal/yaw balance and joint-muscle limits are relaxed; foot-only infeasibility remains a necessary sagittal obstruction to preserving this exact motion without support. It does not rule out changed motions/controllers. No online force allocator.')
    assert Path(output_name).name==output_name;(root/output_name).write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def landing_ik_targets(root,body_pitch_range_deg=0.,body_roll_range_deg=0.,assisted=False):
    """Bounded offline landing-pose references, never injected into simulation."""
    import itertools
    root=Path(root);audit=json.loads((root/'MEASURED_WRENCH_SUPPORT_FEASIBILITY.json').read_text());base=Path(audit['baseline']);z=np.load(base/'trajectory.npz');m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);signs=np.array(list(itertools.product((-1,1),repeat=3)));belt=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_Treadmill_Belt');rows=[];poses=[]
    for sample in audit['samples']:
        if sample['foot_only_all_support_removed_feasible'] or 'loaded_foot_full_projection_m' not in sample:continue
        k=int(round(sample['t_s']/float(z['t'][1]-z['t'][0])));side_index=int(np.argmax(z['foot_vertical'][k]));side='LR'[side_index];q=z['q'][k].copy();d.qpos[:]=q;mujoco.mj_forward(m,d);floor_y=float(d.geom_xpos[belt,1]+m.geom_size[belt,1]);g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_foot_'+side);tg=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_toe_'+side);joint_names=[s+'_'+side for s in ('hip','knee','ankle','toe')];jids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,n) for n in joint_names];adr=m.jnt_qposadr[jids];initial=np.r_[q[1],q[adr]];lower=[];upper=[]
        for j,n in zip(jids,joint_names):
            limits=m.jnt_range[j].copy() if m.jnt_limited[j] else np.array([-.7,.7])
            if n.startswith('knee_'):limits[0]=max(limits[0],np.radians(3-(10 if side=='L' else 55)))
            lower.append(limits[0]);upper.append(limits[1])
        def passive_toe(joint):
            address=int(m.jnt_qposadr[joint]);dof=int(m.jnt_dofadr[joint])
            def balance(v):
                temp=q.copy();temp[address]=v;d.qpos[:]=temp;d.qvel[:]=0;mujoco.mj_forward(m,d);return d.qfrc_bias[dof]-d.qfrc_passive[dof]
            return brentq(balance,-.15,.15,xtol=1e-14)
        initial[4]=passive_toe(jids[3]);lower[3]=max(lower[3],initial[4]-np.radians(5));upper[3]=min(upper[3],initial[4]+np.radians(5));bounds=(np.r_[q[1],lower],np.r_[q[1]+.04,upper]);desired_heel=sample['required_CoP_x_m']-.020
        if body_pitch_range_deg:
            initial=np.r_[initial,0.];bounds=(np.r_[bounds[0],-np.radians(body_pitch_range_deg)],np.r_[bounds[1],np.radians(body_pitch_range_deg)])
        other_toe=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'toe_'+('R' if side=='L' else 'L'));other_adr=int(m.jnt_qposadr[other_toe]);other_dof=int(m.jnt_dofadr[other_toe])
        wing=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+('Left' if side=='R' else 'Right')+'_hip_wing');rail=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+('Left' if side=='R' else 'Right')+'_support_rail')
        if body_roll_range_deg:
            other_knee=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'knee_'+('R' if side=='L' else 'L'));other_knee_adr=int(m.jnt_qposadr[other_knee]);kr=m.jnt_range[other_knee].copy();kr[0]=max(kr[0],np.radians(3-(55 if side=='L' else 10)))
            assert body_pitch_range_deg;initial=np.r_[initial,0.,passive_toe(other_toe),q[other_knee_adr]];bounds=(np.r_[bounds[0],-np.radians(body_roll_range_deg),-.7,kr[0]],np.r_[bounds[1],np.radians(body_roll_range_deg),.7,kr[1]])
        initial=np.clip(initial,bounds[0]+1e-9,bounds[1]-1e-9)
        def evaluate(x):
            candidate=q.copy();candidate[1]=x[0];candidate[adr]=x[1:5]
            if body_pitch_range_deg:
                rot=np.array([np.cos(x[5]/2),0,0,np.sin(x[5]/2)]);mujoco.mju_mulQuat(candidate[3:7],rot,q[3:7])
            if body_roll_range_deg:
                rot=np.array([np.cos(x[6]/2),np.sin(x[6]/2),0,0]);old=candidate[3:7].copy();mujoco.mju_mulQuat(candidate[3:7],rot,old);candidate[other_adr]=x[7];candidate[other_knee_adr]=x[8]
            d.qpos[:]=candidate;d.qvel[:]=0;mujoco.mj_forward(m,d);R=d.geom_xmat[g].reshape(3,3);bottom_axis=int(np.argmax(abs(R[1])));long_axis=max((i for i in range(3) if i!=bottom_axis),key=lambda i:abs(R[0,i]));vertices=d.geom_xpos[g]+(signs*m.geom_size[g])@R.T;rear=vertices[signs[:,long_axis]==-np.sign(R[0,long_axis])];anterior=vertices[signs[:,long_axis]==np.sign(R[0,long_axis])];heel=rear[np.argmin(rear[:,1])];front=anterior[np.argmin(anterior[:,1])];toe=(d.geom_xpos[tg]+(signs*m.geom_size[tg])@d.geom_xmat[tg].reshape(3,3).T)[:,1].min();return candidate,heel,front,toe
        def residual(x):
            candidate,heel,front,toe=evaluate(x);result=1000*np.array([heel[0]-desired_heel,heel[1]-floor_y,max(0,floor_y-front[1]),toe-floor_y])
            if body_roll_range_deg:
                rid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root');com=d.subtree_com[rid];linez=heel[2]+(com[0]-heel[0])/(front[0]-heel[0])*(front[2]-heel[2]);lateral_error=com[2]-linez
                if assisted:
                    vertices=d.geom_xpos[wing]+(signs*m.geom_size[wing])@d.geom_xmat[wing].reshape(3,3).T;lateral_error=vertices[:,1].min()-(d.geom_xpos[rail,1]+m.geom_size[rail,1])
                other_min=np.inf
                for name in ('foot_'+('R' if side=='L' else 'L'),'toe_'+('R' if side=='L' else 'L')):
                    gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+name);other_min=min(other_min,float((d.geom_xpos[gid]+(signs*m.geom_size[gid])@d.geom_xmat[gid].reshape(3,3).T)[:,1].min()))
                result=np.r_[result,1000*max(0,floor_y+.010-other_min),1000*lateral_error,100*(d.qfrc_bias[other_dof]-d.qfrc_passive[other_dof])]
            return result
        original_heel=evaluate(initial)[1].copy();fit=least_squares(residual,initial,bounds=bounds,xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=300);candidate,heel,front,toe=evaluate(fit.x);error=np.array([heel[0]-desired_heel,heel[1]-floor_y,max(0,floor_y-front[1]),toe-floor_y]);valid=bool(max(abs(error))<.0001);pose_number=len(poses);poses.append(candidate);old_angles=z['angles'][k,[0,1,2] if side=='L' else [3,4,5]];new_angles=old_angles+np.degrees(fit.x[1:4]-q[adr[:3]])
        apparatus_penetration=max([-float(con.dist) for con in d.contact if any('Guide_' in (mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,int(geom)) or '') for geom in (con.geom1,con.geom2))]+[0.])
        row=dict(t_s=sample['t_s'],side=side,feasible=valid,optimizer_success=bool(fit.success),optimizer_status=int(fit.status),maximum_foot_constraint_error_m=float(max(abs(error))),foot_constraint_error_m=error.tolist(),stance_toe_passive_residual_Nm=float(d.qfrc_bias[m.jnt_dofadr[jids[3]]]-d.qfrc_passive[m.jnt_dofadr[jids[3]]]),source_hip_knee_ankle_deg=old_angles.tolist(),target_hip_knee_ankle_deg=new_angles.tolist(),root_height_change_m=float(fit.x[0]-q[1]),world_z_body_rotation_change_deg=float(np.degrees(fit.x[5])) if body_pitch_range_deg else 0.,world_x_body_rotation_change_deg=float(np.degrees(fit.x[6])) if body_roll_range_deg else 0.,guide_penetration_m=float(apparatus_penetration),required_heel_shift_m=float(desired_heel-original_heel[0]),target_root_y_m=float(fit.x[0]),pose_index=pose_number)
        if body_roll_range_deg:row['target_swing_knee_flexion_deg']=float((55 if side=='L' else 10)+np.degrees(fit.x[8]))
        if body_roll_range_deg:
            check=residual(fit.x);row['opposite_wing_rail_gap_m' if assisted else 'lateral_CoM_line_error_m']=float(check[-2]/1000);row['other_toe_passive_residual_Nm']=float(check[-1]/100);row['feasible']=bool(valid and abs(check[-2]/1000)<1e-4 and abs(check[-1]/100)<1e-5 and apparatus_penetration<.0005)
        rows.append(row)
    report=dict(machine='easteregg2',baseline=str(base),source_sha256=audit['source_sha256'],landmark_version=3,foot_landmark_rule='Actual lowest world-Y corner on each posterior/anterior endcap, with anterior axis inferred independently. Source foot boxes use localZthickness/localYwidth. Includes lateral tilt rather than assuming a flat localYface.',joint_bound_rule='Exact nominal source ranges, no extra .5degree IK margin; knee>=3degrees. Source model limits unchanged.',body_pitch_range_deg=body_pitch_range_deg,targets=rows,feasible_count=sum(r['feasible'] for r in rows),scope='Offline ROM-respecting geometric references only. Source root x/z retained; root y may rise0–40mm, optional bounded world-Z body rotation explicitly recorded, and four stance-leg joints may change within source limits. Heel/forefoot/toe at belt, required sagittal CoP20mm inside heel. No geometry/ROM edits, no pose injection, no mechanical/neural walking claim. Geometric success requires<.1mm residual independently of solver termination status; muscle/contact feasibility still required.')
    report['body_roll_range_deg']=body_roll_range_deg
    if body_roll_range_deg:report['scope']+=' Additional bounded world-X body rotation aligns whole-body CoM with the actual foot contact edge; swing-toe coordinate satisfies its passive gravity/spring equilibrium. These are offline references, not commanded body torques or guaranteed controllable lateral motions.'
    report['assisted']=assisted
    report['foot_contact_rule']='Heel and toe solid contact, with forefoot nonpenetration rather than forced flatness. Verified standing home has~1.024mm forefoot lift and is therefore incompatible with artificial exact-flatness constraint.'
    report['toe_branch_rule']='Stance toe initialized at native gravity/spring equilibrium and restricted to±5degrees around it for IK branch selection only; source joint remains unlimited. Actual contact forces must balance its passive torque. Swing toe equilibrium and10mm swing clearance checked separately.'
    if assisted:report['scope']+=' Assisted variant replaces unsupported lateral alignment with actual opposite-wing/unchanged-rail contact. Explicit minimum-support force solve still required.'
    filename='LANDING_IK_ASSISTED_TARGETS' if assisted else 'LANDING_IK_BALANCE_TARGETS' if body_roll_range_deg else 'LANDING_IK_BODY_TARGETS' if body_pitch_range_deg else 'LANDING_IK_TARGETS';extra={}
    if (root/(filename+'.json')).exists():
        prior=json.loads((root/(filename+'.json')).read_text());history=prior.get('prior_nonaccepted_diagnostics',[]);history.append({k:prior[k] for k in ('landmark_version','targets','feasible_count')});report['prior_nonaccepted_diagnostics']=history;old=np.load(root/(filename+'.npz'));extra={k:old[k] for k in old.files if k!='qpos'};extra['prior_qpos_'+str(len(history))]=old['qpos']
    np.savez_compressed(root/(filename+'.npz'),qpos=np.array(poses),**extra);(root/(filename+'.json')).write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def landing_force_feasibility(root,filename='LANDING_IK_BALANCE_TARGETS'):
    """Full14-DoF gravity/contact/source-Hill feasibility for geometric teachers."""
    import itertools
    root=Path(root);targets=json.loads((root/(filename+'.json')).read_text());qs=np.load(root/(filename+'.npz'))['qpos'];base=Path(targets['baseline']);bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];params=np.load(base/'POOLED_FORCE_PARAMETERS.npz')['pooled'];m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);rid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root');robot=set()
    for b in range(1,m.nbody):
        a=b
        while a and a!=rid:a=int(m.body_parentid[a])
        if a==rid:robot.add(b)
    dofs=[i for i in range(m.nv) if int(m.jnt_bodyid[m.dof_jntid[i]]) in robot];act=[i for i,b in enumerate(bank) if b['type']=='LinearHillMuscle'];tids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+bank[i]['name']) for i in act];signs=np.array(list(itertools.product((-1,1),repeat=3)));jac=np.empty((3,m.nv));angular=jac.copy();rows=[];floor=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_Treadmill_Belt');bw=float(m.body_subtreemass[rid])*9.81
    for target,q in zip(targets['targets'],qs):
        d.qpos[:]=q;d.qvel[:]=0;mujoco.mj_forward(m,d);top=d.geom_xpos[floor,1]+m.geom_size[floor,1];contacts=[];contact_mu=[];penetration=0.
        for name in ('foot_L','toe_L','foot_R','toe_R'):
            g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+name);vertices=d.geom_xpos[g]+(signs*m.geom_size[g])@d.geom_xmat[g].reshape(3,3).T;penetration=max(penetration,top-float(vertices[:,1].min()))
            for point in vertices[abs(vertices[:,1]-top)<.0001]:mujoco.mj_jac(m,d,jac,angular,point,int(m.geom_bodyid[g]));contacts.append(jac[:,dofs].T.copy());contact_mu.append(min(float(m.geom_friction[g,0]),float(m.geom_friction[floor,0])))
        foot_contact_count=len(contacts);rail_contact_count=0
        if targets.get('assisted',False):
            for side in ('Left','Right'):
                g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+side+'_hip_wing');rg=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+side+'_support_rail');rail_top=d.geom_xpos[rg,1]+m.geom_size[rg,1];vertices=d.geom_xpos[g]+(signs*m.geom_size[g])@d.geom_xmat[g].reshape(3,3).T
                for point in vertices[abs(vertices[:,1]-rail_top)<.0001]:
                    if any(abs(point[axis]-d.geom_xpos[rg,axis])>m.geom_size[rg,axis] for axis in (0,2)):continue
                    mujoco.mj_jac(m,d,jac,angular,point,int(m.geom_bodyid[g]));contacts.append(jac[:,dofs].T.copy());contact_mu.append(float(m.geom_friction[g,0]));rail_contact_count+=1
        bounds=[];law=[]
        for i,tid in zip(act,tids):
            l0,width,kse,kpe,damping,maximum,A,B,C,D,*_=params[i];L=float(d.ten_length[tid]);fl=max(0,1-((L-l0)/width)**2);off=kpe*(L-l0)/(1+kpe/kse);slope=B*fl/(1+kpe/kse)
            def tension(v):return float(np.clip(off+(B/(1+np.exp(min(700,C*(A-v))))+D)*fl/(1+kpe/kse),0,maximum))
            bounds.append((tension(-.1),tension(0)));law.append((i,L,fl,off,slope))
        row=dict(t_s=target['t_s'],side=target['side'],geometric_feasible=target['feasible'],solid_foot_penetration_m=max(0,penetration),contact_points=foot_contact_count,rail_contact_points=rail_contact_count,full_14dof_feasible=False)
        if contacts:
            Aeq=np.column_stack([-d.ten_J.reshape(m.ntendon,m.nv)[np.ix_(tids,dofs)].T,*contacts]);required=(d.qfrc_bias-d.qfrc_passive)[dofs];bounds.extend([(-10*bw,10*bw),(0,10*bw),(-10*bw,10*bw)]*len(contacts));ineq=[]
            for k,mu in enumerate(contact_mu):
                for sx in (-1,1):
                    for sz in (-1,1):
                        a=np.zeros(len(bounds));a[len(act)+3*k:len(act)+3*k+3]=[sx,-mu,sz];ineq.append(a)
            cost=np.r_[1e-5/np.maximum(1,np.array([v[1] for v in bounds[:len(act)]])),np.zeros(3*len(contacts))]
            for k in range(foot_contact_count,len(contacts)):cost[len(act)+3*k+1]=1/bw
            lp=linprog(cost,A_eq=Aeq,b_eq=required,A_ub=ineq,b_ub=np.zeros(len(ineq)),bounds=bounds,method='highs');row['solver_status']=lp.message;row['full_14dof_feasible']=bool(lp.success and target['feasible'] and penetration<.0005 and target['guide_penetration_m']<.0005)
            if lp.success:
                commands=[]
                for value,(i,L,fl,off,slope) in zip(lp.x,law):
                    A,B,C,D=params[i,6:10];activation=np.clip((value-off)/(slope+1e-300),1e-12,1-1e-12);vm=float(A-np.log(1/activation-1)/C) if value>1e-8 else -.1;commands.append(dict(bank_index=i,muscle=bank[i]['name'],tension_N=float(value),motor_voltage_V=vm,length_m=L,length_factor=fl))
                Fn=lp.x[len(act)+1::3];row.update(max_equilibrium_residual=float(max(abs(Aeq@lp.x-required))),muscle_tensions_N={bank[i]['name']:float(v) for i,v in zip(act,lp.x)},motor_commands=commands,normal_foot_load_BW=float(Fn[:foot_contact_count].sum()/bw),upward_rail_support_BW=float(Fn[foot_contact_count:].sum()/bw),toe_passive_torque_Nm={s:float(d.qfrc_passive[m.jnt_dofadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'toe_'+s)]]) for s in 'LR'})
            else:
                nd=len(required);slack=linprog(np.r_[cost,np.full(2*nd,1e3)],A_eq=np.column_stack([Aeq,np.eye(nd),-np.eye(nd)]),b_eq=required,A_ub=np.column_stack([ineq,np.zeros((len(ineq),2*nd))]),b_ub=np.zeros(len(ineq)),bounds=bounds+[(0,None)]*(2*nd),method='highs')
                if slack.success:row['minimum_equilibrium_residual_by_dof']={str(i)+':'+str(mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,int(m.dof_jntid[j]))):float(v) for i,(j,v) in enumerate(zip(dofs,Aeq@slack.x[:len(bounds)]-required)) if abs(v)>1e-7}
        rows.append(row)
    report=dict(target_file=str(root/(filename+'.json')),source_sha256=targets['source_sha256'],results=rows,full_14dof_feasible_count=sum(r['full_14dof_feasible'] for r in rows),scope='Independent full14robotDoF static gravity/source-Hill/contact feasibility, including BOTH passive toe joints and actual pooled source force-length/sigmoid bounds−100..0mV. Actual transformed solid-foot corners within.1mm of floor; inscribed μcontact diamond. Rail reactions ONLY in explicitly assisted variant at actual geometric wing/rail contact, minimized and separately reported; no guide/cap wrench, body force, or applied torque. Reject geometry/guide penetration>.5mm. Static force witness only; not dynamic standing, stepping, neural fit, or proof unactuated lateral posture can be tracked.')
    (root/(filename+'_FORCE_FEASIBILITY.json')).write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def mirrored_landing_phase(root):
    """Mechanical mirror check, not independent left/right neural fitting."""
    root=Path(root);name='LANDING_IK_ASSISTED_TARGETS';ik=json.loads((root/(name+'.json')).read_text());forces=json.loads((root/(name+'_FORCE_FEASIBILITY.json')).read_text());index=next(i for i,r in enumerate(forces['results']) if r['full_14dof_feasible']);original=np.load(root/(name+'.npz'))['qpos'][index];m=mujoco.MjModel.from_xml_path(str(Path(ik['baseline'])/'model.xml'));d=mujoco.MjData(m);mirror=original.copy();mirror[2]*=-1;rotation=np.empty(9);mujoco.mju_quat2Mat(rotation,original[3:7]);S=np.diag([1.,1.,-1.]);mujoco.mju_mat2Quat(mirror[3:7],(S@rotation.reshape(3,3)@S).ravel())
    for joint,offset in [('hip',30.),('knee',45.),('ankle',0.),('toe',0.)]:
        ids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,joint+'_'+s) for s in 'LR'];adr=m.jnt_qposadr[ids];mirror[adr[0]]=original[adr[1]]+np.radians(offset);mirror[adr[1]]=original[adr[0]]-np.radians(offset)
    stance='L' if ik['targets'][index]['side']=='R' else 'R';other='R' if stance=='L' else 'L';toe=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'toe_'+other);address=int(m.jnt_qposadr[toe]);dof=int(m.jnt_dofadr[toe])
    def equilibrium(value):
        d.qpos[:]=mirror;d.qpos[address]=value;d.qvel[:]=0;mujoco.mj_forward(m,d);return d.qfrc_bias[dof]-d.qfrc_passive[dof]
    mirror[address]=brentq(equilibrium,-.15,.15,xtol=1e-14);d.qpos[:]=mirror;mujoco.mj_forward(m,d);penetration=max([-float(c.dist) for c in d.contact if any('Guide_' in (mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,int(g)) or '') for g in (c.geom1,c.geom2))]+[0.]);row=dict(ik['targets'][index],side=stance,pose_index=1,guide_penetration_m=float(penetration),mirror_of_pose_index=index,scope='Exact sagittal reflection plus source-native anatomical joint-reference offsets. Swing toe reset to its own actual gravity/passive equilibrium, not a fitted gain.')
    report=dict(machine='easteregg2',baseline=ik['baseline'],source_sha256=ik['source_sha256'],assisted=True,source_reference=str(root/(name+'.json')),source_index=index,targets=[dict(ik['targets'][index],pose_index=0),row],scope='One accepted mechanical phase and its anatomical mirror; no neural parameters fitted independently and no model geometry changed. Independent force/solid-contact gate is authoritative.')
    np.savez_compressed(root/'LANDING_PHASE_PAIR_TARGETS.npz',qpos=np.array([original,mirror]));(root/'LANDING_PHASE_PAIR_TARGETS.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return landing_force_feasibility(root,'LANDING_PHASE_PAIR_TARGETS')


def map_phase_pair_to_neurons(root):
    """Actual source action/current map for the mechanically valid phase pair.

    One signed chemical projection per functional pool, shared across sides.
    This is a constrained frozen-current prototype, not a trained walker or a
    complete-cycle teacher. No joint torque/force is commanded online.
    """
    import inspect
    from stance_recruitment import ROLES
    from spiking_central_runtime import SpikingCentralController as C
    from recruitment import apply_recruitment
    root=Path(root);pair_path=root/'LANDING_PHASE_PAIR_TARGETS_FORCE_FEASIBILITY.json';pair=json.loads(pair_path.read_text());assert all(r['full_14dof_feasible'] for r in pair['results']);meta=json.loads((root/'LANDING_PHASE_PAIR_TARGETS.json').read_text());base=Path(meta['baseline']);rec=json.loads((base/'SPIKING_COMPARISON.json').read_text());cfg=rec['source_configuration'];bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];params=np.load(base/'POOLED_FORCE_PARAMETERS.npz')['pooled'];z=np.load(base/'trajectory.npz');C.configuration=rec['configuration'];C.balance_configuration=None;C.sensory_spiking=True;C.clock_scale=1.;C.audit_state=False;C.reference_calculations=False;allowed=set(inspect.signature(C.__bases__[0]).parameters)-{'source','bank'};c=C(cfg['source'],bank,**{k:v for k,v in cfg.items() if k in allowed});apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale']);c.apply_output_gain();binding={i:c.names[pre] for pre,i,coeff in c.outputs};reverse={tail:role for role,tail in ROLES.items()};rows=[]
    for phase in pair['results']:
        time=phase['t_s']+(.456 if phase['side']=='L' else 0.);k=int(round(time/rec['result']['physics_dt']));voltage=z['neuro'][k]
        for command in phase['motor_commands']:
            i=command['bank_index'];name=binding[i];post=c.names.index(name);role=reverse[name[2:]];is_stance=name.startswith(phase['side']);vtarget=1000*command['motor_voltage_V'];zero=command['tension_N']<1e-8;p=params[i];L=command['length_m'];fl=command['length_factor'];offset=p[3]*(L-p[0])/(1+p[3]/p[2]);slope=p[7]*fl/(1+p[3]/p[2]);ceiling=None
            if zero:
                activation=-offset/max(1e-300,slope)
                if 0<activation<1:ceiling=float(1000*(p[6]-np.log(1/activation-1)/p[8]));vtarget=max(-100,min(0,ceiling))
            assert c.params[post,2]==1 and c.params[post,3]==0 and not any(j==post for j,*_ in c.inputs)
            current=c.params[post,0]-vtarget+c.tonic[post];incoming=[]
            for e in c.edges[c.edges[:,1]==post]:
                active=e[2]*np.clip((voltage[int(e[0])]-e[3])/(e[4]-e[3]),0,1);I=float(active*(e[5]-vtarget));current+=I;incoming.append(dict(source=c.names[int(e[0])],g_active_uS=float(active),reversal_mV=float(e[5]),current_at_target_nA=I))
            needed=-float(current)
            if zero:needed=0. if ceiling is None or voltage[post]<=vtarget else min(0.,needed)
            excitation=max(0.,needed)/max(1e-6,-10-vtarget);inhibition=max(0.,-needed)/max(1e-6,vtarget+100);old=next(v for v in c.output_edits if v['target']==name);cap=old['before_gmax_uS']*90
            rows.append(dict(mechanical_stance_side=phase['side'],role=role,phase='stance' if is_stance else 'swing',MN=name,physical_route=command['muscle'],target_tension_N=command['tension_N'],target_motor_voltage_mV=1000*command['motor_voltage_V'],effective_zero_force_ceiling_mV=ceiling,neural_budget_voltage_mV=vtarget,reference_time_s=time,reference_MN_mV=float(voltage[post]),required_additional_current_nA=needed,exc_gmax_uS=float(excitation),inh_gmax_uS=float(inhibition),source_pre_attenuation_current_envelope_nA=float(cap),required_current_envelope_ratio=float(abs(needed)/max(1e-9,cap)),incoming=sorted(incoming,key=lambda v:abs(v['current_at_target_nA']),reverse=True)))
    settings=json.loads((root/'STANCE_RECRUITMENT_CONFIGURATION.json').read_text());settings['offline_target_sha256']=hashlib.sha256(pair_path.read_bytes()).hexdigest()
    for role in ROLES:
        observations=[v for v in rows if v['role']==role and v['phase']=='stance'];assert len(observations)==2
        negligible=max(abs(v['required_additional_current_nA'])/v['source_pre_attenuation_current_envelope_nA'] for v in observations)<.01
        for kind in ('exc','inh'):settings[role+'_'+kind+'_uS']=0. if negligible else float(np.mean([v[kind+'_gmax_uS'] for v in observations]))
    amplitudes=[]
    for phase in pair['results']:
        time=phase['t_s']+(.456 if phase['side']=='L' else 0.);k=int(round(time/rec['result']['physics_dt']));i=c.names.index(phase['side']+' RG ext');amplitudes.append(float(z['neuro'][k,i]-c.params[i,0]))
    settings['rg_ext_reference_mV']=float(np.mean(amplitudes))
    ready=all(0<=settings[k]<=15 for k in settings if k.endswith('_uS')) and all(v['required_current_envelope_ratio']<=1 for v in rows if v['phase']=='stance')
    report=dict(machine='easteregg2',source_sha256=meta['source_sha256'],mechanical_teacher=str(pair_path),mechanical_teacher_complete_cycle=False,baseline=str(base),rows=rows,shared_stance_configuration=settings,ready_for_isolated_current_probe=bool(ready),scope='Two mirror-related mechanically verified STATIC phase anchors. Desired muscle force inverted through actual source sigmoid/passive length law OFFLINE; source chemical currents evaluated at desired voltage using actual representative spike-proxy snapshots. Zero-force targets are voltage ceilings, not forced excitation. One signed chemical drive per pool, then mirrored conductances averaged; corrections below1%source-current envelope omitted. Original spike-decoded RG-E amplitude shapes the new PF drive, normalized at the verified phase, rather than an immediate full plateau. Baseline RG/PF/intrinsic/force parameters unchanged. No neural training, complete-cycle teacher, online allocation, or walking success claimed.')
    (root/'NEURAL_PHASE_CURRENT_MAP.json').write_text(json.dumps(report,indent=2),encoding='utf-8');(root/'LANDING_PRIMITIVE_NEURAL_CONFIGURATION.json').write_text(json.dumps(settings,indent=2),encoding='utf-8');return report


def audit_hold_drift(root):
    """Saved10ms traces: report drift without rerunning a completed hold."""
    root=Path(root);fit=json.loads((root/'OFFLINE_NEURAL_FIT.json').read_text());m=mujoco.MjModel.from_xml_path(str(Path(fit['baseline'])/'model.xml'));joints=[n+'_'+s for s in 'LR' for n in ('hip','knee','ankle','toe')];jids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,n) for n in joints];qa=m.jnt_qposadr[jids];va=m.jnt_dofadr[jids];rows=[]
    for folder,prefix in [(root/'functional_target','MECHANICAL_HOLD'),(root,'LANDING_PHASE_MECHANICAL_HOLD'),(root,'LANDING_MIRROR_MECHANICAL_HOLD')]:
        z=np.load(folder/(prefix+'.npz'))
        for condition in ('cold_start','equilibrium_preloaded'):
            t=z[condition+'_t'];q=z[condition+'_q'];dv=np.empty(m.nv);mujoco.mj_differentiatePos(m,dv,float(t[-1]-t[-2]),q[-2],q[-1]);i=int(np.searchsorted(t,t[-1]-.1));row=dict(source=str(folder/(prefix+'.npz')),condition=condition,last_saved_t_s=float(t[-1]),terminal_velocity_estimation='Finite difference of final two saved10msposes; not exact terminal simulator qvel.',root_velocity_m_s=dv[:3].tolist(),terminal_joint_velocity_deg_s={n:float(v) for n,v in zip(joints,np.degrees(dv[va]))},last_segment_duration_s=float(t[-1]-t[i]),last_segment_root_displacement_m=(q[-1,:3]-q[i,:3]).tolist(),last_segment_joint_change_deg={n:float(v) for n,v in zip(joints,np.degrees(q[-1,qa]-q[i,qa]))},last_segment_mean_foot_BW=(z[condition+'_foot_N'][i:].mean(axis=0)/(m.body_subtreemass[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root')]*9.81)).tolist(),steady_state_established=False,scope='Equilibrium-preloaded means initial muscle internal tension only, not a dynamically converged body state. Short0.5slanding holds are feasibility observations; source2shome hold also drifts.')
            rows.append(row)
    report=dict(machine='easteregg2',results=rows,scope='Saved-data drift audit only; no duplicate simulation. Do not call short holds steady equilibria or neural balance/stepping results.');(root/'MECHANICAL_HOLD_DRIFT_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def audit_precontact_motor_forcing(root,selected_path=None):
    """Replay source muscle ODE at recorded pre-step poses; gate force parity first."""
    import inspect
    from native_kernel import polynomial
    root=Path(root);cases=json.loads((root/'INTERLOCK_GATES.json').read_text())['cases'];rows=[];output=root/'PRECONTACT_MOTOR_FORCING_AUDIT.json'
    if selected_path is not None and output.exists():rows=[r for r in json.loads(output.read_text())['results'] if r['path']!=str(Path(selected_path))]
    for case in cases:
        if case['mode']=='air':continue
        if selected_path is not None and Path(case['path'])!=Path(selected_path):continue
        base=Path(case['path']);rec=json.loads((base/'SPIKING_COMPARISON.json').read_text());cfg=rec['source_configuration'];bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];p=np.load(base/'POOLED_FORCE_PARAMETERS.npz')['pooled'];z=np.load(base/'trajectory.npz');m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);allowed=set(inspect.signature(NativeController).parameters)-{'source','bank'};c=NativeController(cfg['source'],bank,**{k:v for k,v in cfg.items() if k in allowed});N=len(z['t']);nb=len(bank);tids=np.array([mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+b['name']) for b in bank]);vm=np.full((N,nb),-.15);binding={}
        for pre,i,coeff in c.outputs:vm[:,i]=polynomial(coeff,z['neuro'][:,z['neuron_names'].tolist().index(c.names[pre])]*.001);binding[c.names[pre]]=i
        L=np.empty((N,nb));moment=np.empty((N,nb));last=None;parts=np.zeros((3,nb));trace=np.empty((N,3,nb));force=np.empty((N,nb));slope=np.empty((N,nb));dt=rec['result']['physics_dt'];counter=np.zeros((2,nb));counterforces=np.empty((N,2,nb));rest=np.array(cfg['initial_qpos']);decay=1-dt*p[:,2]/p[:,4]*(1+p[:,3]/p[:,2])
        for k in range(N):
            d.qpos[:]=rest if k==0 else z['q'][k-1];mujoco.mj_forward(m,d);L[k]=d.ten_length[tids];velocity=np.zeros(nb) if last is None else (L[k]-last)/dt;last=L[k].copy();fl=np.maximum(0,1-((L[k]-p[:,0])/p[:,1])**2);sig=1/(1+np.exp(np.clip(p[:,8]*(p[:,6]-vm[k]),-700,700)));active=np.maximum(0,p[:,7]*sig+p[:,9])*fl;forcing=np.array([active,p[:,3]*(L[k]-p[:,0]),p[:,4]*velocity]);parts=parts*decay+dt*p[:,2]/p[:,4]*forcing;trace[k]=parts;force[k]=np.minimum(p[:,5],np.maximum(0,parts.sum(axis=0)));slope[k]=p[:,7]*p[:,8]*sig*(1-sig)*fl*.001
            for j,offset in enumerate((-.005,.005)):
                act=np.maximum(0,p[:,7]/(1+np.exp(np.clip(p[:,8]*(p[:,6]-vm[k]-offset),-700,700)))+p[:,9])*fl;counter[j]=counter[j]*decay+dt*p[:,2]/p[:,4]*(forcing[1]+forcing[2]+act);counterforces[k,j]=np.minimum(p[:,5],np.maximum(0,counter[j]))
            J=d.ten_J.reshape(m.ntendon,m.nv)
            for name,i in binding.items():
                joint=('hip' if 'Hip' in name else 'knee' if 'Knee' in name else 'ankle')+'_'+name[0];dof=m.jnt_dofadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,joint)];moment[k,i]=-J[tids[i],dof]
        err=float(np.max(np.abs(force-z['forces'])));assert err<1e-7,('Saved muscle-force parity failed',str(base),err)
        groups=[]
        for side in 'LR':
            events=[v for v in case['landings'] if v['side']==side];mask=np.zeros(N,dtype=bool)
            for event in events:mask|=(z['t']>=event['touchdown_s']-.15)&(z['t']<event['touchdown_s'])
            for tail in ('Hip MN ext','Hip MN flx','Knee MN ext','Knee MN flx','Ank MN ext','Ank MN flx'):
                name=side+' '+tail;i=binding[name];passive=np.minimum(p[i,5],np.maximum(0,trace[:,1:,i].sum(axis=1)));groups.append(dict(MN=name,route=bank[i]['name'],sigmoid_A_V=float(p[i,6]),sigmoid_B_N=float(p[i,7]),sigmoid_C_perV=float(p[i,8]),sample_count=int(mask.sum()),precontact_mean_voltage_mV=float(vm[mask,i].mean()*1000),precontact_mean_actual_tension_N=float(force[mask,i].mean()),precontact_mean_signed_torque_Nm=float((force[:,i]*moment[:,i])[mask].mean()),precontact_mean_passive_only_counterfactual_N=float(passive[mask].mean()),precontact_mean_raw_internal_components_N=trace[mask,:,i].mean(axis=0).tolist(),mean_instantaneous_active_slope_N_per_mV=float(slope[mask,i].mean()),frozen_motion_minus5mV_mean_tension_N=float(counterforces[mask,0,i].mean()),frozen_motion_plus5mV_mean_tension_N=float(counterforces[mask,1,i].mean())))
        rows.append(dict(path=str(base),max_all_muscle_force_replay_error_N=err,precontact=groups))
    report=dict(results=rows,scope='No new coupled trial. All source Hill forces reconstructed from actual recorded pre-step poses, finite-difference length velocity, motor voltages, zero initial internal state, and exact source Euler ODE. Force parity required<1e−7N. Raw active/elastic/velocity internal components may be signed and are not unique additive clipped tensions. Passive-only and±5mV counterfactuals hold recorded motion fixed, so diagnose actuator response, not accepted neural control or mechanical feasibility. Precontact=last150ms before identical force-qualified touchdown.');(root/'PRECONTACT_MOTOR_FORCING_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('target');parser.add_argument('parameters',nargs='?');parser.add_argument('--duration',type=float,default=2.);parser.add_argument('--fit-neural',action='store_true');parser.add_argument('--protect-transient',action='store_true');parser.add_argument('--dynamic-support',action='store_true');args=parser.parse_args()
    result=audit_measured_support(args.target) if args.dynamic_support else fit_neural_recruitment(args.target,args.protect_transient) if args.fit_neural else mechanical_hold(args.target,args.parameters,args.duration);print(json.dumps(result,indent=2))
