"""Audit full 3D motion, model constraints, and native muscle authority."""
import os,sys,json,hashlib,xml.etree.ElementTree as ET
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parent;S=H.parent;sys.path[:0]=[str(S),str(S/'functional_gait2392')]
import numpy as np
from scipy.signal import find_peaks
import mujoco,bsolve_ik as B
from mechanics import CoupledArms,JOINTS
O=H/'runs/ik_audit';O.mkdir(parents=True,exist_ok=True)
SOURCE=S.parents[2]/'Solid_Models/OpenSim/Gait2392_Robotbody'
def read_motion(path):
    lines=path.read_text(encoding='utf-8',errors='replace').splitlines()
    end=next(i for i,l in enumerate(lines) if l.strip().lower()=='endheader')
    names=lines[end+1].split();values=np.array([[float(x) for x in l.split()[:len(names)]] for l in lines[end+2:] if l.strip()])
    unit=next((l.split('=',1)[1].strip().lower() for l in lines[:end] if l.lower().startswith('indegrees=')),None)
    return values[:,0],names[1:],values[:,1:],unit
def locks(path):
    root=ET.parse(path).getroot()
    return {e.attrib['name']:dict(locked=e.findtext('locked'),prescribed=e.findtext('prescribed'),
        default=e.findtext('default_value'),range=e.findtext('range')) for e in root.iter('Coordinate')}
paths=[SOURCE/'subject01_walk1_ik.mot',SOURCE/'Tutorial1/normal.mot']
paths+=list((H/'runs/ik_sources/falisse/Results').glob('*/motion.mot'))
records=[];curves={}
for path in paths:
    t,names,v,unit=read_motion(path);name=path.parent.name+'_'+path.stem
    factor=180/np.pi if unit=='no' else 1.
    coordinates={}
    for joint in JOINTS:
        if joint not in names:
            coordinates[joint]=dict(status='missing');continue
        y=v[:,names.index(joint)]*factor
        span=float(np.ptp(y));std=float(np.std(y))
        coordinates[joint]=dict(status='constant' if span<1e-3 else 'varying',range_deg=span,
            mean_deg=float(y.mean()),min_deg=float(y.min()),max_deg=float(y.max()))
        # Kinematic dependency, not a claim about neural PF identification.
        if span>1e-3:
            side=joint[-1] if joint.endswith(('_r','_l')) else 'r'
            predictors=[j for j in ('hip_flexion_'+side,'knee_angle_'+side,'ankle_angle_'+side) if j in names and j!=joint]
            X=np.column_stack([np.ones(len(t))]+[v[:,names.index(j)] for j in predictors])
            keep=np.arange(len(t))%5!=0
            beta=np.linalg.lstsq(X[keep],y[keep],rcond=None)[0];error=y[~keep]-X[~keep]@beta
            coordinates[joint]['heldout_sagittal_linear_R2']=float(1.-np.sum(error**2)/max(np.sum((y[~keep]-np.mean(y[~keep]))**2),1e-12))
            coordinates[joint]['dependency_scope']='Interleaved sample holdout within one trajectory, not subject-level validation'
        if joint.startswith('lumbar'):
            pelvis={'lumbar_extension':'pelvis_tilt','lumbar_bending':'pelvis_list','lumbar_rotation':'pelvis_rotation'}[joint]
            if pelvis in names:
                py=v[:,names.index(pelvis)]*factor
                coordinates[joint]['std_lumbar_plus_pelvis_deg']=float(np.std(y+py))
    records.append(dict(name=name,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        samples=len(t),duration_s=float(t[-1]-t[0]),rotation_units='degrees' if unit=='yes' else 'radians' if unit=='no' else 'unspecified',coordinates=coordinates))
    curves[name]=(t,names,v,factor)
models={str(p):locks(p) for p in [SOURCE/'subject01_simbody.osim',H/'runs/ik_sources/falisse/OpenSimModel/new_model/Model/new_model_scaled_FK_contacts.osim'] if p.exists()}
m,d=B.build_air_model(False);arms=CoupledArms(m)
source=np.load(S/'bsolve_out.npz',allow_pickle=True)
maxarm=np.zeros((len(JOINTS),m.nu));axis=[]
for q in source['qpos'][::3]:
    d.qpos[:]=q;arms.project(d);mujoco.mj_forward(m,d)
    maxarm=np.maximum(maxarm,abs(arms.update(d)))
authority={joint:dict(muscles_over_1mm=int(np.sum(maxarm[k]>.001)),max_moment_arm_m=float(maxarm[k].max()),
    muscles=[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in np.flatnonzero(maxarm[k]>.001)]) for k,joint in enumerate(JOINTS)}
result=dict(cached_references=len(list((S/'gait_refs').glob('*.npz'))),cached_coordinate_scope='Bilateral sagittal hip/knee/ankle only',
    full_motion_sources=records,model_coordinates=models,muscle_authority=authority,
    interpretation='Motion patterns and dependencies do not identify neural PF populations. Missing/constant columns are not evidence of absent human motion.')
(O/'audit.json').write_text(json.dumps(result,indent=2))
np.savez_compressed(O/'muscle_authority.npz',max_moment_arms=maxarm,joints=JOINTS,actuators=source['act_names'])
for record in records:
    print(record['name'])
    for joint,info in record['coordinates'].items():
        if joint.startswith(('lumbar','mtp','subtalar','hip_adduction','hip_rotation')):
            print(joint,info.get('status'),round(info.get('range_deg',0),4),round(info.get('heldout_sagittal_linear_R2',0),3))
print('AUTHORITY',[(j,v['muscles_over_1mm']) for j,v in authority.items()])
