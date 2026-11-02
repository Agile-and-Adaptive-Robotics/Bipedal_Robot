"""CoM-constrained standing initialization and static Hill-force feasibility.

No neural or physical parameter is changed, and static forces are not applied
to the subsequent SNS simulation.
"""
import argparse,itertools,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares,linprog,lsq_linear
from scipy.spatial import ConvexHull
from import_virtual_walker import mujoco
from native_kernel import bank_parameters

p=argparse.ArgumentParser();p.add_argument('baseline');p.add_argument('output');p.add_argument('--plot-only',action='store_true');p.add_argument('--recheck-static-only',action='store_true');p.add_argument('--symmetric-voltage-target',action='store_true');p.add_argument('--force-parameters');p.add_argument('--minimum-motor-voltage-V',type=float,default=-.1);a=p.parse_args();base=Path(a.baseline);out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
if (out/'STANDING_HOME.json').exists() and not (a.plot_only or a.recheck_static_only):raise FileExistsError('Preserving completed standing-home audit; use a new output folder')
r=json.loads((base/'RESULT.json').read_text());m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);q0=np.array(r['initial_qpos']);d.qpos[:]=q0;mujoco.mj_forward(m,d);rid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root');signs=np.array(list(itertools.product((-1,1),repeat=3)));feet=[]
for side in 'LR':
 for name in ('foot_'+side,'toe_'+side):
  g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+name);R=d.geom_xmat[g].reshape(3,3);axis=int(np.argmax(abs(R[1])));lower=-np.sign(R[1,axis]);corners=signs[signs[:,axis]==lower]*m.geom_size[g];feet.append((g,corners,int(m.geom_bodyid[g])))
floor=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Treadmill_Belt');fg=next(i for i in range(m.ngeom) if m.geom_bodyid[i]==floor);top=float(np.max((d.geom_xpos[fg]+(signs*m.geom_size[fg])@d.geom_xmat[fg].reshape(3,3).T)[:,1]));qa={name:int(m.jnt_qposadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,name)]) for name in ['hip_L','knee_L','ankle_L','toe_L','hip_R','knee_R','ankle_R','toe_R']}
def geometry(q):
 d.qpos[:]=q;d.qvel[:]=0;mujoco.mj_forward(m,d);points=np.vstack([d.geom_xpos[g]+corners@d.geom_xmat[g].reshape(3,3).T for g,corners,_ in feet]);com=d.subtree_com[rid].copy();hull=ConvexHull(points[:,[0,2]]);margins=-(hull.equations[:,:2]@com[[0,2]]+hull.equations[:,2]);return points,com,float(margins.min())
if a.plot_only:
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42});saved=json.loads((out/'STANDING_HOME.json').read_text());fig,axes=plt.subplots(1,2,figsize=(7.2,6.3),layout='constrained')
 for ax,q,title,color,marker in zip(axes,[q0,np.array(saved['qpos'])],['(A) Earlier near-flat home','(B) CoM-constrained home'],['.3','#0072B2'],['o','s']):
  points,com,margin=geometry(q);origin=d.xpos[rid,0];anchor=[]
  for name in ('hip_L','knee_L','ankle_L'):
   j=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,name);anchor.append(d.xanchor[j,[0,1]]-[origin,top])
  anchor=np.array(anchor);ax.plot(anchor[:,0],anchor[:,1],color=color,lw=2,marker=marker,ms=5)
  for g,_,_ in feet[:2]:
   vertices=d.geom_xpos[g]+(signs*m.geom_size[g])@d.geom_xmat[g].reshape(3,3).T;xy=vertices[:,[0,1]]-[origin,top];hull=ConvexHull(xy);poly=xy[np.r_[hull.vertices,hull.vertices[0]]];ax.plot(poly[:,0],poly[:,1],color=color,lw=1)
  xmin,xmax=points[:,0].min()-origin,points[:,0].max()-origin;cx,cy=com[[0,1]]-[origin,top];ax.axhline(0,color='.65',lw=1);ax.plot([xmin,xmax],[-.006,-.006],lw=4,color=color,alpha=.5);ax.plot([cx,cx],[0,cy],color=color,ls='--',lw=1);ax.plot(cx,cy,marker=marker,color=color,ms=8);ax.annotate('Whole-model CoM',(cx,cy),xytext=(cx+.025,cy-.08),arrowprops={'arrowstyle':'-','color':color},fontsize=10)
  ax.annotate('',(xmin,-.028),(cx,-.028),arrowprops={'arrowstyle':'<->','color':color});ax.text((xmin+cx)/2,-.047,f'{1000*margin:+.0f} mm',ha='center',va='top',fontsize=10,color=color);ax.set(title=title,xlabel='Anterior position relative to pelvis (m)',xlim=(-.09,.29),ylim=(-.075,1.01));ax.set_aspect('equal');ax.grid(alpha=.15)
 axes[0].set_ylabel('Height above treadmill (m)');fig.suptitle('Virtual W2L initial-pose audit: one leg shown; whole-model CoM',fontsize=11);fig.savefig(out/'COM_ALIGNMENT.png',dpi=220);fig.savefig(out/'COM_ALIGNMENT.pdf');plt.close(fig);(out/'COM_ALIGNMENT_ALT.txt').write_text('Sagittal pose audit using actual model joint/foot geometry and whole-body center of mass. Earlier near-flat pose: CoM about69mm behind the full foot support projection. Corrected source-ROM-respecting pose: CoM20mm inside. One leg shown; supports omitted for geometric clarity. Static feasibility is not a neural or dynamic balance pass. No inverse forces are applied to the neural tests.');print(out/'COM_ALIGNMENT.png');sys.exit(0)
def pose(x):
 y,pitch,hip,knee,ankle,toeL,toeR=x;q=q0.copy();q[1]=y;q[3:7]=[np.cos(pitch/2),0,0,np.sin(pitch/2)]
 for side,h0,k0 in [('L',-10,10),('R',20,55)]:q[qa['hip_'+side]]=np.radians(hip-h0);q[qa['knee_'+side]]=np.radians(knee-k0);q[qa['ankle_'+side]]=np.radians(ankle)
 q[qa['toe_L']]=toeL;q[qa['toe_R']]=toeR;return q
def residual(x):
 points,com,margin=geometry(pose(x));springs=[i for i in range(m.ntendon) if (mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_TENDON,i) or '').startswith('spring_')];return np.r_[[100000*(block[:,1].min()-top) for block in points.reshape(4,4,3)],1000*max(0.,.020-margin),.02*np.degrees(x[1]),.05*(x[3]-8),.02*(x[4]-1),10000*(d.ten_length[springs]-m.tendon_lengthspring[springs,0])]
if a.symmetric_voltage_target:q=q0.copy()
elif a.recheck_static_only:q=np.array(json.loads((out/'STANDING_HOME.json').read_text())['qpos'])
else:fit=least_squares(residual,[.955,0,1.7,4.2,2.49,0,0],bounds=([.85,np.radians(-15),-17.5,3.5,-11.5,-.7,-.7],[1.02,np.radians(5),19.5,25,2.5,.7,.7]),max_nfev=500,xtol=1e-12,ftol=1e-12,gtol=1e-12);q=pose(fit.x);points,com,margin=geometry(q)
bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];params=bank_parameters(bank);act=[i for i,b in enumerate(bank) if b['type']=='LinearHillMuscle'];tids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+bank[i]['name']) for i in act];robot=set()
if a.force_parameters:
 table=np.load(a.force_parameters);assert table['muscle_names'].tolist()==[b['name'] for b in bank];params=table['pooled'].copy()
for body in range(m.nbody):
 ancestor=body
 while ancestor and ancestor!=rid:ancestor=int(m.body_parentid[ancestor])
 if ancestor==rid:robot.add(body)
dofs=[i for i in range(m.nv) if int(m.jnt_bodyid[m.dof_jntid[i]]) in robot]

if a.symmetric_voltage_target:
 from symmetric_stance import solve_target
 target=out/'SYMMETRIC_VOLTAGE_TARGET.json';assert not target.exists(),'Preserve completed target'
 result=solve_target(m,d,q,bank,params,dofs,feet,r,base,a.minimum_motor_voltage_V)
 target.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ('muscles','contact_forces_N','qpos')},indent=2));sys.exit(0)
def force_feasibility(q):
 points,com,margin=geometry(q);moment=-d.ten_J.reshape(m.ntendon,m.nv)[np.ix_(tids,dofs)].T;contact=[];jac=np.empty((3,m.nv));angular=np.empty((3,m.nv))
 for g,corners,body in feet:
  vertices=d.geom_xpos[g]+corners@d.geom_xmat[g].reshape(3,3).T
  for point in vertices[vertices[:,1]<=vertices[:,1].min()+.00005]:mujoco.mj_jac(m,d,jac,angular,point,body);contact.append(jac[:,dofs].T.copy())
 A=np.column_stack([moment,*contact]);required=(d.qfrc_bias-d.qfrc_passive)[dofs].copy();bounds=[];capacities=[]
 for i,tid in zip(act,tids):
  l0,width,kse,kpe,damping,maximum,va,B,C,D,*_=params[i];length=d.ten_length[tid];fl=max(0.,1.-((length-l0)/width)**2)
  def tension(vm):return float(np.clip((kpe*(length-l0)+max(0.,B/(1+np.exp(min(700.,C*(va-vm))))+D)*fl)/(1+kpe/kse),0,maximum))
  low,high=tension(a.minimum_motor_voltage_V),tension(0.);bounds.append((low,high));capacities.append(dict(muscle=bank[i]['name'],minimum_command_tension_N=low,maximum_full_command_N=high,length_factor=fl))
 bw=r['robot_mass_kg']*9.81;bounds.extend([(-bw,bw),(0,bw),(-bw,bw)]*len(contact));cost=np.r_[1/np.maximum(1,np.array([b[1] for b in bounds[:len(act)]])),np.tile([0,1e-8,0],len(contact))];ineq=[]
 for i in range(len(contact)):
  for axis in (0,2):
   for direction in (-1,1):
    row=np.zeros(len(cost));row[len(act)+3*i+axis]=direction;row[len(act)+3*i+1]=-.5;ineq.append(row)
 lp=linprog(cost,A_eq=A,b_eq=required,A_ub=np.array(ineq),b_ub=np.zeros(len(ineq)),bounds=bounds,method='highs',options={'time_limit':10.});result=dict(feasible=bool(lp.success),status=lp.message,sagittal_posterior_margin_m=float(com[0]-points[:,0].min()),polygon_min_margin_m=margin,friction_pyramid_mu=.5,contact_points=len(contact),capacities=capacities)
 if lp.success:
  result.update(max_equilibrium_residual=float(np.max(abs(A@lp.x-required))),muscle_tensions_N={bank[i]['name']:float(t) for i,t in zip(act,lp.x)},normal_contact_total_N=float(lp.x[len(act)+1::3].sum()))
 else:
  slack=linprog(np.r_[cost,np.full(2*len(required),1e6)],A_eq=np.column_stack([A,np.eye(len(required)),-np.eye(len(required))]),b_eq=required,A_ub=np.column_stack([np.array(ineq),np.zeros((len(ineq),2*len(required)))]),b_ub=np.zeros(len(ineq)),bounds=bounds+[(0,None)]*(2*len(required)),method='highs')
  if slack.success:result['friction_constrained_residual_by_dof']={str(i)+':'+str(mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,int(m.dof_jntid[j]))):float(value) for i,(j,value) in enumerate(zip(dofs,A@slack.x[:len(cost)]-required))}
  lower=np.array([b[0] for b in bounds]);upper=np.array([b[1] for b in bounds]);free=upper>lower;fixed=np.where(free,0,lower);least=lsq_linear(A[:,free],required-A@fixed,bounds=(lower[free],upper[free]),max_iter=1000,tol=1e-12);estimate=fixed.copy();estimate[free]=least.x;error=A@estimate-required
  result.update(bounded_equilibrium_residual_by_dof={str(i)+':'+str(mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,int(m.dof_jntid[j]))):float(value) for i,(j,value) in enumerate(zip(dofs,error))},bounded_max_residual=float(np.max(abs(error))),unbounded_matrix_rank=int(np.linalg.matrix_rank(A)),gravity_required_by_dof=required.tolist(),bounded_force_estimate=estimate.tolist())
 return result
if a.recheck_static_only:
 result=dict(minimum_motor_voltage_V=a.minimum_motor_voltage_V,force_parameters=a.force_parameters,original_home=force_feasibility(q0),CoM_constrained_home=force_feasibility(q),scope='Saved pose unchanged. Static contact-friction and source Hill equilibrium checked over the declared motor-voltage range, not a dynamic or neural pass. No forces applied to SNS.');target=out/('STATIC_RECHECK_POOLED.json' if a.force_parameters else 'STATIC_RECHECK_ORIGINAL.json');assert not target.exists(),'Preserve completed static recheck';target.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));sys.exit(0)
old=force_feasibility(q0);new=force_feasibility(q);points,com,margin=geometry(q);errors=np.array([block[:,1].min()-top for block in points.reshape(4,4,3)]);summary=dict(machine='easteregg2',baseline=str(base),source_sha256=r['source_sha256'],qpos=q.tolist(),optimization_success=bool(fit.success),parameters=dict(root_y_m=float(fit.x[0]),root_pitch_deg=float(np.degrees(fit.x[1])),hip_flexion_deg=float(fit.x[2]),knee_flexion_deg=float(fit.x[3]),ankle_dorsiflexion_deg=float(fit.x[4]),toe_q_rad=fit.x[5:].tolist()),maximum_sole_height_error_m=float(np.max(abs(errors))),COM_world_m=com.tolist(),support_polygon_min_margin_m=margin,original_home=old,CoM_constrained_home=new,scope='Symmetric hip/knee/ankle pose, existing geometry/masses/ROM. At least 0.5 degree hip/ankle ROM margin, knee >=3.5 degrees. Foot and toe boxes may touch on their lowest edges; no sole penetration. Passive toe springs preserved. CoM inside standing support hull. Static feasibility uses source Hill force bounds and contact friction pyramid mu=0.5, with no guide/rail reaction. No static forces applied to SNS trials. Not a dynamic balance pass.')
np.savez_compressed(out/'standing_home.npz',qpos=q);(out/'STANDING_HOME.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k not in ('qpos','original_home','CoM_constrained_home')},indent=2));print('Original static feasible:',old['feasible'],'New static feasible:',new['feasible'])
assert fit.success and summary['maximum_sole_height_error_m']<.0005 and margin>.019 and new['feasible'],'Home or static feasibility checks failed; do not promote'
