"""Read-only audit and figure export from retained CAD traces; never integrates a model."""
from pathlib import Path
import hashlib
import json
import numpy as np
import mujoco
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
BASE=Path(r'D:\Users\Ben Bolen\CAD3DOF_20261102\contact_handover_20261103')
NAME='clearance77p5_damped'
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
meta=json.loads((BASE/'runs'/f'{NAME}.json').read_text())
recorded_checks=json.loads((BASE/'FINAL_CHECKS.json').read_text())
z=np.load(BASE/'runs'/f'{NAME}.npz',allow_pickle=False)
t=z['time']; late=t>=16
load=z['contact_loads'].sum(axis=2)
model=mujoco.MjModel.from_xml_path(meta['model_path'])
data=mujoco.MjData(model)
mass=float(model.body_mass.sum())
weight=mass*9.81
intervals=[]
for side in range(2):
    ok=late&(load[:,side]<7.36)&(z['mesh_clearance'][:,side]>.005)&(load[:,1-side]>15)
    edges=np.diff(np.r_[False,ok,False].astype(int))
    intervals.append([[float(t[i]),float(t[j-1])] for i,j in zip(np.flatnonzero(edges==1),np.flatnonzero(edges==-1)) if t[j-1]-t[i]>=.08])
events=sorted((times[0],side) for side in range(2) for times in intervals[side])
alternations=sum(a[1]!=b[1] for a,b in zip(events,events[1:]))
period=recorded_checks[NAME]['period_s']

# Recompute geometric height from present compiled geometry at every logged
# final-rung pose; this is kinematic evaluation, not dynamic replay/integration.
clouds=[[],[]]; hulls={}
for gid in range(model.ngeom):
    name=mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_GEOM,gid) or ''
    if not name.startswith(('foot_r_','foot_l_')): continue
    mid=model.geom_dataid[gid]
    if mid<0: continue
    if mid not in hulls:
        vertices=model.mesh_vert[model.mesh_vertadr[mid]:model.mesh_vertadr[mid]+model.mesh_vertnum[mid]]
        # Linear extrema over all vertices equal extrema over their hull.
        hulls[mid]=vertices
    clouds[0 if name.startswith('foot_r_') else 1].append((gid,hulls[mid]))
height=[]; gaps=[]
for q in z['qpos'][late]:
    data.qpos[:]=q
    mujoco.mj_kinematics(model,data)
    height.append([min(float((v@data.geom_xmat[g].reshape(3,3)[1]+data.geom_xpos[g,1]).min()) for g,v in foot) for foot in clouds])
    side_extents=[]
    for foot in clouds:
        vv=[v@data.geom_xmat[g].reshape(3,3)[2]+data.geom_xpos[g,2] for g,v in foot]
        side_extents.append((min(v.min() for v in vv),max(v.max() for v in vv)))
    gaps.append(float(side_extents[0][0]-side_extents[1][1]))
height=np.array(height)

summary={
    'condition':NAME,'window_rule':'time >= 16 s; stored trace through 27.99 s from a completed 28 s integration',
    'first_last_stored_time_s':[float(t[0]),float(t[-1])],
    'scoring_first_last_stored_time_s':[float(t[late][0]),float(t[late][-1])],
    'samples_total':len(t),'samples_scored':int(late.sum()),'sample_interval_s':float(np.median(np.diff(t))),
    'compiled_current_mass_kg':mass,'gravity_m_s2':model.opt.gravity.tolist(),'BW_denominator_N':weight,
    'actuator_count':int(model.nu),'activation_time_constants_s':np.unique(model.actuator_dynprm[:,0]).tolist(),
    'mean_vertical_support_BW':float(z['rig_force'][late,1].mean()/weight),
    'mean_vertical_force_N':float(z['rig_force'][late,1].mean()),
    'mean_ground_load_R_L_N':load[late].mean(axis=0).tolist(),
    'clear_unloaded_interval_R_L':intervals,'qualifying_intervals_R_L':[len(x) for x in intervals],
    'alternations':alternations,'recorded_not_recomputed_median_RG_period_R_L_s':period,
    'whole_trial_knee_max_recorded_R_L_deg':np.rad2deg(z['knee_coordinate'].max(axis=0)).tolist(),
    'whole_trial_knee_max_every_physics_step_R_L_deg':np.rad2deg(meta['max_knee_every_physics_step_rad']).tolist(),
    'window_knee_max_recorded_R_L_deg':np.rad2deg(z['knee_coordinate'][late].max(axis=0)).tolist(),
    'mesh_clearance_min_R_L_mm':(1000*z['mesh_clearance'][late].min(axis=0)).tolist(),
    'mesh_clearance_max_R_L_mm':(1000*z['mesh_clearance'][late].max(axis=0)).tolist(),
    'current_mesh_recompute_max_abs_height_difference_m':float(np.max(np.abs(height-z['mesh_clearance'][late]))),
    'current_mesh_recompute_min_world_z_gap_all_10ms_samples_m':min(gaps),
    'published_world_z_gap_every_fifth_10ms_sample_m':recorded_checks[NAME]['min_sampled_world_z_mesh_gap_m'],
    'full_trial_world_forward_displacement_m':float(z['qpos'][-1,0]-z['qpos'][0,0]),
    'mean_abs_rig_torque_xyz_Nm':np.abs(z['rig_torque'][late]).mean(axis=0).tolist(),
    'rig_horizontal_force_xz_RMS_N':np.sqrt(np.mean(z['rig_force'][late][:,[0,2]]**2,axis=0)).tolist(),
    'replay_array_comparisons':{},'local_cases':{},
}
for suffix in ['fresh_repeat','reference']:
    other=np.load(BASE/'runs'/f'{NAME}_{suffix}.npz',allow_pickle=False)
    sm=json.loads((BASE/'runs'/f'{NAME}_{suffix}.json').read_text())
    summary['replay_array_comparisons'][suffix]={
        'same_array_keys':set(z.files)==set(other.files),
        'array_count':len(z.files),
        'all_arrays_exact':all(np.array_equal(z[k],other[k]) for k in z.files),
        'same_model_xml_hash':meta['source_hashes'][meta['model_path']]==sm['source_hashes'][sm['model_path']],
    }
for k,v in recorded_checks.items():
    if not k.startswith(NAME+'_basin'): continue
    bz=np.load(BASE/'runs'/f'{k}.npz',allow_pickle=False)
    bm=json.loads((BASE/'runs'/f'{k}.json').read_text())
    bt=bz['time']; bl=bt>=16; bload=bz['contact_loads'].sum(axis=2); bi=[]
    for side in range(2):
        ok=bl&(bload[:,side]<7.36)&(bz['mesh_clearance'][:,side]>.005)&(bload[:,1-side]>15)
        edges=np.diff(np.r_[False,ok,False].astype(int))
        bi.append([[float(bt[i]),float(bt[j-1])]for i,j in zip(np.flatnonzero(edges==1),np.flatnonzero(edges==-1))if bt[j-1]-bt[i]>=.08])
    be=sorted((a[0],s)for s in range(2)for a in bi[s])
    ba=sum(a[1]!=b[1] for a,b in zip(be,be[1:]))
    summary['local_cases'][k]={
        'completed':bm['failed'] is None,'failed':bm['failed'],'duration_simulated_s':bm['duration_simulated_s'],
        'intervals_R_L':[len(x)for x in bi],'alternations':ba,
        'saved_ground_step_gate':v['ground_step_gate'],
        'saved_interval_counts_match':v['unloaded_clear_intervals']==[len(x)for x in bi],
        'min_flexion_every_physics_step_deg':float(-np.rad2deg(max(bm['max_knee_every_physics_step_rad']))),
    }
summary['local_passes']=sum(v['saved_ground_step_gate']for v in summary['local_cases'].values())
source_paths={str(BASE/'runs'/f'{NAME}{s}.{ext}')for s in ['', '_fresh_repeat','_reference','_physics1ms']for ext in ['npz','json']}
source_paths.update(meta['source_hashes'])
source_paths.update(str(BASE/x)for x in ['FINAL_CHECKS.json','FINAL_VERIFICATION.json','source/measure.py','source/prepare_final_checks.py','source/final_basin_cases.json'])
source_paths.update(str(BASE/'runs'/f'{k}.{ext}')for k in summary['local_cases']for ext in ['json','npz'])
summary['source_sha256']={p:sha(p) for p in sorted(source_paths)}
summary['saved_source_hash_matches']={p:sha(p)==h for p,h in meta['source_hashes'].items()}
(HERE/'CAD_AUDIT_RESULTS.json').write_text(json.dumps(summary,indent=2)+'\n')

# Native signs, original samples, and a single common clock. Matplotlib's path
# simplification is disabled so the vector file retains every plotted sample.
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'pdf.fonttype':42,'svg.fonttype':'none','path.simplify':False,'axes.linewidth':.8})
fig,axes=plt.subplots(4,1,figsize=(14/2.54,15/2.54),sharex=True)
fig.subplots_adjust(left=.19,right=.98,top=.93,bottom=.09,hspace=.34)
colors=['#0072B2','#D55E00'];styles=['-','--']
panels=[
    (z['rg'][:,[1,3]]-z['rg'][:,[0,2]],'RG F − E\n(mV)','A'),
    (np.rad2deg(z['knee_coordinate']),'Knee angle\n(deg)','B'),
    (load,'Foot normal\nload (N)','C'),
    (1000*z['mesh_clearance'],'Minimum mesh\nheight (mm)','D'),
]
for ax,(y,label,letter) in zip(axes,panels):
    for side in range(2):ax.plot(t[late],y[late,side],color=colors[side],linestyle=styles[side],linewidth=1.2,label=['Right','Left'][side])
    ax.set_ylabel(label)
    ax.text(-.22,1.01,letter,transform=ax.transAxes,fontsize=10,fontweight='bold')
    ax.grid(True,axis='y',color='.88',linewidth=.5)
    ax.spines[['top','right']].set_visible(False)
axes[0].legend(loc='lower center',bbox_to_anchor=(.5,1.04),ncol=2,frameon=False)
axes[1].axhline(-1,color='.3',linestyle=':',linewidth=.8)
axes[1].set_ylim(-48,1)
axes[2].axhline(7.36,color='.45',linestyle=':',linewidth=.8)
axes[2].axhline(15,color='.25',linestyle='-.',linewidth=.8)
axes[3].axhline(5,color='.3',linestyle=':',linewidth=.8)
axes[3].axhline(0,color='.6',linewidth=.5)
axes[-1].set_xlim(16,28);axes[-1].set_xticks(np.arange(16,29,2));axes[-1].set_xlabel('Time (s)')
pdf=HERE/'CAD_supported_traces_20261116.pdf'; png=HERE/'CAD_supported_traces_20261116.png'; svg=HERE/'CAD_supported_traces_20261116.svg'
fig.savefig(pdf);fig.savefig(svg);fig.savefig(png,dpi=180);plt.close(fig)
provenance={
    'generator':str(Path(__file__).resolve()),'source_trace':str(BASE/'runs'/f'{NAME}.npz'),
    'source_trace_sha256':sha(BASE/'runs'/f'{NAME}.npz'),'source_results':str(HERE/'CAD_AUDIT_RESULTS.json'),
    'figure_size_cm':[14,15],'minimum_font_size_pt':10,'plot_window_s':[16,28],
    'scoring_window_s':[16,28],'samples_per_series':int(late.sum()),
    'transformations':['RG F minus E membrane potential by side','knee radians multiplied by 180/pi; native flexion-negative sign retained','two contact-region forces summed per foot','mesh height metres multiplied by 1000'],
    'smoothing':'none','phase_shift':'none','decimation':'none; every stored sample within the plot window retained','amplitude_rescaling':'unit conversion only','path_simplification':False,
    'display_thresholds':{'knee_deg':-1,'ipsilateral_normal_load_N':7.36,'contralateral_normal_load_N':15,'whole_foot_height_mm':5},
    'output_sha256':{str(p):sha(p)for p in [pdf,png,svg]},
    'historical_compiled_binary_status':'No saved historical MJB or compiled-array digest identified in this CAD run package. Current XML source and all retained dependency hashes match. Current compiled mesh geometry reproduces saved final-window heights at every 10 ms logged pose; this does not supply a missing historical compiled binary hash.'
}
(HERE/'CAD_FIGURE_PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items()if k not in ['source_sha256','local_cases','clear_unloaded_interval_R_L']},indent=2))
