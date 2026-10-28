"""Separate kinematic PCA primitives from muscle NMF synergies, with block tests."""
import os,sys,json,hashlib,warnings
from pathlib import Path
HERE=Path(__file__).resolve().parent;S=HERE.parent;sys.path.insert(0,str(HERE))
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import mujoco
from sklearn.decomposition import NMF
from scipy.optimize import nnls
from mechanics import JOINTS
OUT=S/'reports_20261028_functional_pf';OUT.mkdir(exist_ok=True)
d=np.load(S/'bsolve_out.npz',allow_pickle=True)
model=mujoco.MjModel.from_xml_path(str(S.parents[2]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'))
addresses=[model.jnt_qposadr[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,n)] for n in JOINTS]
Q=d['qpos'][:,addresses];velocity=np.gradient(Q,d['t'],axis=0)
X=np.column_stack((Q,velocity));A=np.maximum(d['acts'],0.)
def vaf(y,p):return float(1-np.sum((y-p)**2)/max(np.sum(y*y),1e-12))
records=[]
for rank in (2,4,6,8):
    for block in range(3):
        test=np.array_split(np.arange(len(X)),3)[block];train=np.setdiff1d(np.arange(len(X)),test)
        mean=X[train].mean(axis=0);scale=np.maximum(X[train].std(axis=0),.005)
        x=(X-mean)/scale;_,singular,V=np.linalg.svd(x[train],full_matrices=False)
        pc=x@V[:rank].T;design=np.column_stack((np.ones(len(X)),pc))
        factor=NMF(n_components=rank,init='nndsvda',random_state=20261028,max_iter=3000,tol=1e-5)
        H=factor.fit_transform(A[train]);W=factor.components_
        oracle=np.array([nnls(W.T,row)[0] for row in A[test]])@W
        ridge=np.linalg.solve(design[train].T@design[train]+.01*np.eye(rank+1),design[train].T@H)
        prediction=np.maximum(design[test]@ridge,0.)@W
        records.append(dict(rank=rank,block=block,heldout_frames=test.tolist(),
            muscle_reconstruction_vaf=vaf(A[test],oracle),ik_to_muscle_prediction_vaf=vaf(A[test],prediction),
            kinematic_training_variance=float(np.sum(singular[:rank]**2)/np.sum(singular**2))))
# Learn temporal movement waveforms independently from the reference library.
curves=[];origins=[]
for path in sorted((S/'gait_refs').glob('*.npz')):
    if not path.name.startswith(('ong_','subject','falisse','Falisse')):continue
    z=np.load(path,allow_pickle=True)
    if not all(side+'_'+j in z for side in ('r','l') for j in ('hip','knee','ankle')):continue
    for side in ('r','l'):
        # The degenerate left cycle is explicitly excluded by the repository ruling.
        if path.stem=='ong_selfsel_Init200' and side=='l':continue
        for j in ('hip','knee','ankle'):
            curve=np.asarray(z[side+'_'+j],float)
            if not np.isfinite(curve).all() or np.ptp(curve)<1.:continue
            curve=np.interp(np.linspace(0,1,100),np.linspace(0,1,len(curve)),curve)
            curves.append(curve);origins.append(path.stem)
C=np.array(curves);centers=C.mean(axis=1,keepdims=True);amplitudes=np.maximum(np.std(C,axis=1,keepdims=True),1.)
normalized=(C-centers)/amplitudes
_,singular,temporal=np.linalg.svd(normalized,full_matrices=False)
temporal_tests=[];unique=sorted(set(origins))
for fold in range(4):
    held=set(unique[fold::4]);test=np.array([o in held for o in origins]);train=~test
    _,_,basis=np.linalg.svd(normalized[train],full_matrices=False)
    for rank in (2,4,6,8):
        prediction=(normalized[test]@basis[:rank].T)@basis[:rank]
        temporal_tests.append(dict(fold=fold,rank=rank,heldout_archives=sorted(held),
            normalized_rmse=float(np.sqrt(np.mean((prediction-normalized[test])**2)))))
np.savez_compressed(OUT/'primitive_decomposition.npz',time=d['t'],joint_names=JOINTS,
    joint_coordinates_rad=Q,joint_velocity_rad_s=velocity,temporal_movement_basis=temporal[:8],
    temporal_singular_values=singular)
report=dict(machine='easteregg2',paired_frames=len(X),paired_duration_s=float(np.ptp(d['t'])),
    paired_source='Existing quasi-static IK + measured-GRF back-solved activations, not measured EMG',
    block_cross_validation=records,temporal_kinematic_archive_count=len(unique),temporal_curves=len(C),
    temporal_cross_validation=temporal_tests,
    limits=['Two-second paired record is insufficient to identify independent biological modules.',
            'Oracle NMF reconstruction uses held-out activation data to solve coefficients; IK prediction does not.',
            'Kinematic PCA and activation NMF are separate mathematical representations, not proof of separate neural circuits.',
            'Reference-archive holdouts may share subjects and speeds; no unseen-subject generalization claim.'],
    source_sha256=hashlib.sha256((S/'bsolve_out.npz').read_bytes()).hexdigest())
(OUT/'primitive_analysis.json').write_text(json.dumps(report,indent=2))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':11})
fig,axes=plt.subplots(1,3,figsize=(16,4.8),constrained_layout=True)
ranks=[2,4,6,8]
for key,label,style,color in [('muscle_reconstruction_vaf','Activation reconstruction','o-','#0077BB'),('ik_to_muscle_prediction_vaf','IK-only activation prediction','s--','#CC6677')]:
    means=[np.mean([r[key] for r in records if r['rank']==k]) for k in ranks]
    axes[0].plot(ranks,means,style,color=color,label=label)
axes[0].axhline(0,color='.7',lw=.8);axes[0].set(xlabel='Retained components',ylabel='Held-out uncentered VAF',title='A  Muscle-space fit and kinematic prediction');axes[0].legend(fontsize=10)
for k,style in zip(range(4),('-','--','-.',':')):axes[1].plot(np.arange(100),temporal[k],style,label='Movement mode '+str(k+1))
axes[1].set(xlabel='Normalized cycle (%)',ylabel='Temporal basis (arbitrary units)',title='B  Kinematic temporal primitives');axes[1].legend(fontsize=10)
axes[2].plot(ranks,[np.mean([r['normalized_rmse'] for r in temporal_tests if r['rank']==k]) for k in ranks],'D-',color='#228833')
axes[2].set(xlabel='Temporal movement modes',ylabel='Held-out normalized RMSE',title='C  Holdout by reference archive')
fig.savefig(OUT/'motor_primitives_vs_synergies.png',dpi=200);plt.close(fig)
(OUT/'motor_primitives_vs_synergies.alt.txt').write_text('Three panels distinguish held-out muscle reconstruction from IK-only activation prediction, show four independent kinematic temporal basis curves, and measure temporal-basis reconstruction on held-out gait-reference archives. The paired activation data cover only two seconds and are back-solved, not measured EMG.')
print(json.dumps(dict(paired_frames=len(X),archives=len(unique),means={str(k):{key:float(np.mean([r[key] for r in records if r['rank']==k])) for key in ('muscle_reconstruction_vaf','ik_to_muscle_prediction_vaf')} for k in ranks}),indent=2))
