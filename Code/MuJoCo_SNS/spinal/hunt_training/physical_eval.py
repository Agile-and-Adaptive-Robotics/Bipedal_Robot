"""Fresh-process corrected curriculum with measured harness assistance."""
import os,sys,json,time,hashlib
from pathlib import Path
H=Path(__file__).resolve().parent;S=H.parent
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent);os.environ['AARL_NET']='w2lvar'
sys.path[:0]=[str(H),str(S),str(S/'functional_gait2392')]
path=Path(sys.argv[1]).resolve();cfg=json.loads(path.read_text());out=path.with_suffix('');os.environ['AARL_NPZ']=str(out)+'.npz'
import numpy as np,mujoco
import controller as C,runner as R,build_network as BN
from mechanics import CoupledArms,JOINTS
from curriculum_metrics import assess
from frontal_metrics import assess_frontal
stage=cfg['stage'];capture={};extra=np.zeros((8500,21));balance_log=np.zeros((8500,8));count=0
source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),H/'controller.py',H/'balance.py')}
oldpose=R.apply_start_pose
def pose(m,d):
    oldpose(m,d);targets=np.load(H/'runs/targets/targets.npz',allow_pickle=True)
    ar=CoupledArms(m);mean=targets['qpos'][:,ar.addresses].mean(axis=0)
    for k in range(7):mean[k]=mean[k+7]=.5*(mean[k]+mean[k+7])
    mean[[1,8]]=.04;mean[[2,9,5,12,6,13]]=0.
    for label,row in (('hip',0),('knee',3),('ankle',4)):
        if 'start_'+label+'_deg' in cfg:mean[[row,row+7]]=np.radians(cfg['start_'+label+'_deg'])
    d.qpos[ar.addresses]=mean;ar.project(d);mujoco.mj_forward(m,d)
    slides=[m.joint('pelvis_'+axis).id for axis in ('tx','ty','tz')]
    vertical=max(slides,key=lambda j:abs(d.xaxis[j,2]));ad=m.jnt_qposadr[vertical]
    d.qpos[ad]+=cfg.get('height_offset',0.)/d.xaxis[vertical,2]
    ar.project(d);mujoco.mj_forward(m,d)
R.apply_start_pose=pose
original_harness=R.apply_harness
def harness(*args,**kwargs):
    m=original_harness(*args,**kwargs);d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    # Restrict only out-of-plane coordinates. Sagittal leg motion stays free.
    restriction=cfg.get('restriction',0.)
    for name in ('pelvis_list','pelvis_rotation','hip_adduction_r','hip_adduction_l','hip_rotation_r','hip_rotation_l','subtalar_angle_r','subtalar_angle_l','lumbar_bending','lumbar_rotation'):
        j=m.joint(name);m.jnt_stiffness[j.id]+=restriction*800.;m.dof_damping[j.dofadr[0]]+=restriction*40.
    for axis in ('tx','ty','tz'):
        j=m.joint('pelvis_'+axis);vertical=abs(d.xaxis[j.id,2])>.8
        if vertical:
            m.jnt_stiffness[j.id]*=cfg.get('vertical_stiffness_scale',1.)
            m.dof_damping[j.dofadr[0]]*=np.sqrt(cfg.get('vertical_stiffness_scale',1.))
            if cfg.get('force_controlled_support',False):
                m.jnt_stiffness[j.id]=0.;m.dof_damping[j.dofadr[0]]=0.
        elif abs(d.xaxis[j.id,1])>.8:
            scale=cfg.get('lateral_stiffness_scale',1.);m.jnt_stiffness[j.id]*=scale;m.dof_damping[j.dofadr[0]]*=np.sqrt(scale)
    return m
R.apply_harness=harness
def build(acts,dt,interleg):
    options=dict(use_timing_afferents=False if stage==1 else cfg.get('timing_afferents',True),axis_feedback=0. if stage==1 else cfg.get('axis_feedback',0.),balance=None if stage==1 else cfg.get('balance'),weights_directory=cfg.get('weights_directory',str(H/'runs/fitted')),timing_scale=cfg.get('timing_scale',1.),reciprocal_gain=cfg.get('reciprocal_gain',0.),commissural_gain=0. if stage==1 else cfg.get('commissural_gain',0.),biarticular_pair_gain=cfg.get('biarticular_pair_gain',0.))
    if 'hip_abduction_goal' in cfg:options['hip_abduction_goal']=cfg['hip_abduction_goal']
    options['cross_contact_gain']=cfg.get('cross_contact_gain',0.)
    n=C.build(acts,dt=dt,options=options)
    C.P.G['contact_onset']=cfg.get('contact_onset',0.)
    if cfg.get('backend','sparse')=='native':n.compiled.forward=n.compiled.c.forward
    capture['neural_step_wall_s']=0.
    n.ff_gain*=cfg.get('motor_scale',1.);n.reflex_gain*=cfg.get('reflex_scale',1.) if stage!=1 else 0.
    for i in n.pf_edges:
        act=n.names[n.compiled.post[i]][3:];info=n.muscles[act]
        n.ff_gain[i]*=cfg.get('motor_scales',{}).get(info.groups[0]+'_'+info.side,1.)
    if stage==5:
        for i,(pre,post) in enumerate(zip(n.compiled.pre,n.compiled.post)):
            if n.names[post].startswith('PF_') and n.names[pre].startswith(('HEEL_','TOEDF_')):
                n.compiled.gmax[i]=n.original_g[i]*cfg.get('contact_pf_scale',0.)
    native_step=n.step
    def advance(u):
        global count
        for side in ('r','l'):
            u[n.ports['HUNT_INIT_E_'+side]]=2. if side=='r' and count<50 and stage!=3 else 0.
            u[n.ports['HUNT_INIT_F_'+side]]=2. if side=='l' and count<50 and stage!=3 else 0.
        # Optional single command-onset initialization; never a stride clock.
        if 2.<=count*dt<2.2 and cfg.get('walk_kick',0.)>0.:
            u[n.ports['HUNT_INIT_E_r']]=cfg['walk_kick'];u[n.ports['HUNT_INIT_F_l']]=cfg['walk_kick']
        m,d=n.live;pb=m.body('pelvis').id
        # A finite, explicitly declared 100 ms pelvis impulse, never a phase clock.
        d.xfrc_applied[pb,0]=cfg.get('push_N',0.) if cfg.get('push_time',7.)<=count*dt<cfg.get('push_time',7.)+.1 else 0.
        d.xfrc_applied[pb,1]=cfg.get('lateral_push_N',0.) if cfg.get('push_time',7.)<=count*dt<cfg.get('push_time',7.)+.1 else 0.
        d.xfrc_applied[pb,2]=cfg.get('support_fraction',0.)*m.body_mass.sum()*9.81
        step_start=time.perf_counter();native_step(u);capture['neural_step_wall_s']+=time.perf_counter()-step_start
        if hasattr(n,'balance') and count<len(balance_log):balance_log[count]=np.r_[n.balance.last_task,n.balance.last_command]
        support=0.;wrench=np.zeros(6)
        for axis in ('tx','ty','tz'):
            j=m.joint('pelvis_'+axis);ad=j.qposadr[0];vd=j.dofadr[0]
            scalar=-m.jnt_stiffness[j.id]*(d.qpos[ad]-m.qpos_spring[ad])-m.dof_damping[vd]*d.qvel[vd]
            wrench[:3]+=scalar*d.xaxis[j.id]
        support=wrench[2]
        for name in ('pelvis_tilt','pelvis_list','pelvis_rotation'):
            j=m.joint(name);ad=j.qposadr[0];vd=j.dofadr[0]
            scalar=-m.jnt_stiffness[j.id]*(d.qpos[ad]-m.qpos_spring[ad])-m.dof_damping[vd]*d.qvel[vd]
            wrench[3:]+=scalar*d.xaxis[j.id]
        loads=np.zeros(2);grf=np.zeros((2,3));force=np.zeros(6)
        for ci in range(d.ncon):
            contact=d.contact[ci];g1,g2=contact.geom1,contact.geom2
            if m.geom_type[g1]!=mujoco.mjtGeom.mjGEOM_PLANE and m.geom_type[g2]!=mujoco.mjtGeom.mjGEOM_PLANE:continue
            body=m.geom_bodyid[g2 if m.geom_type[g1]==mujoco.mjtGeom.mjGEOM_PLANE else g1]
            while body:
                name=m.body(body).name
                if name in ('talus_r','calcn_r','toes_r','talus_l','calcn_l','toes_l'):
                    mujoco.mj_contactForce(m,d,ci,force)
                    side=0 if name.endswith('_r') else 1;world=contact.frame.reshape(3,3).T@force[:3]
                    if world[2]<0:world=-world
                    grf[side]+=world;loads[side]+=world[2];break
                body=m.body_parentid[body]
        if count<len(extra):extra[count]=[count*dt,support,d.xfrc_applied[pb,2],np.mean(d.act**2),np.mean(d.ctrl**2),d.ncon,np.max(abs(n.compiled.V)),*loads,*wrench,*grf.ravel()]
        count+=1;return n.compiled.V
    n.step=advance;capture.update(network=n,neurons=len(n.compiled.V),inputs=len(n.inputs));return n
BN.build=build
args=['--eval','--drive','2.5']
if stage in (1,2):args+=['--no-ground']
if stage==1:args+=['--no-afferents','--no-interleg']
if stage==3:args=['--stand-eval','8']
start=time.perf_counter();metrics=R.main(args)
assessment=assess(str(out)+'.npz',metrics,5 if stage in (4,5) else stage)
assessment.update(assess_frontal(str(out)+'.npz'))
assessment['accepted']=bool(assessment['accepted'] and assessment['frontal_gate'])
n=capture.pop('network');m,d=n.live;BW=float(m.body_mass.sum()*9.81);late=extra[:count,0]>=5
z=np.load(str(out)+'.npz',allow_pickle=True);ct=extra[:count,7:9]/BW;window=(extra[:count,0]>=5)&(extra[:count,0]<13)
late_window_available=bool(np.any(window))
if not np.any(window):window=extra[:count,0]>=0.
if not np.any(late):late=extra[:count,0]>=0.
measured=dict(mean_spring_support_bw=float(extra[:count][late,1].mean()/BW),explicit_support_bw=cfg.get('support_fraction',0.),
    mean_activation_squared=float(extra[:count][late,3].mean()),mean_excitation_squared=float(extra[:count][late,4].mean()),
    mean_foot_normal_bw=np.mean(ct[window],axis=0).tolist(),body_weight_N=BW,
    ground_enabled=stage>=3,contact_samples=int(np.count_nonzero(np.sum(ct[window],axis=1)>.1)))
measured.update(lateral_harness_rms_bw=float(np.sqrt(np.mean(extra[:count][window,10]**2))/BW),
    lateral_harness_peak_bw=float(np.max(abs(extra[:count][window,10]))/BW),
    harness_wrench_frame='Native world x anterior, y left, z superior. Anatomical right force = -world y. Torque is pelvis joint-restraint generalized torque projected on instantaneous axes; not a physical single attachment load cell.')
measured['horizontal_harness_rms_bw']=float(np.sqrt(np.mean(np.sum(extra[:count][window,9:11]**2,axis=1)))/BW)
measured['late_window_available']=late_window_available
# Independently count each physical contact once. The legacy runner contact
# logger assigns the same force to both heel and toe; its totals are doubled.
assessment['ground_load_verified']=bool(late_window_available and stage>=3 and np.mean(ct[window].sum(axis=1))>.2 and np.all(np.max(ct[window],axis=0)>.1))
if stage in (4,5):assessment['accepted']=bool(assessment['accepted'] and assessment['ground_load_verified'])
result=dict(stage=stage,configuration=cfg,assessment=assessment,metrics=metrics,measured=measured,compiled=capture,
    wall_s=time.perf_counter()-start,trained_weights={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(cfg.get('weights_directory',H/'runs/fitted')).glob('*weights.json')},
    source_sha256=source_hashes,
    effort_scope='Dimensionless activation/excitation squared, not metabolic watts; only a secondary objective after movement gates.')
result['simulated_s']=count*.002
result['backend']=cfg.get('backend','sparse')
if hasattr(n,'balance'):
    result['balance_calibration']=n.balance.calibration
    np.savez_compressed(str(out)+'.balance.npz',t=extra[:count,0],state_and_command=balance_log[:count])
np.savez_compressed(str(out)+'.assistance.npz',samples=extra[:count],columns=['t','spring_support_N','explicit_support_N','activation_squared','excitation_squared','contacts','max_voltage','right_vertical_GRF_N','left_vertical_GRF_N','harness_Fx_N','harness_Fy_N','harness_Fz_N','harness_Mx_Nm','harness_My_Nm','harness_Mz_Nm','right_GRFx_N','right_GRFy_N','right_GRFz_N','left_GRFx_N','left_GRFy_N','left_GRFz_N'])
Path(str(out)+'.result.json').write_text(json.dumps(result,indent=2,default=float));print(json.dumps(result,default=float),flush=True)
if cfg.get('export_graph'):
    c=n.compiled;standing=c.gmax.copy();n.apply_weights(1.);walking=c.gmax.copy();c.gmax[:]=standing;n.sync()
    graph=dict(nodes=[dict(id='n'+str(i),label=n.names[i]) for i in range(len(c.V))],
        edges=[dict(source=int(pre),target=int(post),reversal=float(c.reversal[i]),standing=float(standing[i]),walking=float(walking[i]),inherited=float(n.original_g[i])) for i,(pre,post) in enumerate(zip(c.pre,c.post))],
        scope='Recorded end-state conductances and a separate gait-expression-gating projection at the same pose. The projection is not proof of walking.',inputs=n.inputs)
    Path(str(out)+'.graph.json').write_text(json.dumps(graph))
