"""Saved-array verification and review media. No integration or fitting."""
import json,sys,hashlib,shutil
from pathlib import Path
import numpy as np
import summarize_gait_dynamic as H
X=H.C.X;F=H.C.F;R=H.R
from PIL import Image,ImageDraw,ImageFont

def verify():
    p=R/'physics/syn6_central_v4_8p0';b=R/'physics/syn6_central_v4_8p0_fresh';checks={}
    for ext in ('.npz','.central.npz','.inputs.npz','.release_state.npz','.apparatus.npz','.apparatus_extended.npz','.motor_state.npz'):
        with np.load(p.with_suffix(ext)) as a,np.load(b.with_suffix(ext)) as z:
            assert a.files==z.files
            checks[ext]={key:bool(np.array_equal(a[key],z[key],equal_nan=True) if a[key].dtype.kind in 'fc' else np.array_equal(a[key],z[key])) for key in a.files}
            assert all(checks[ext].values()),checks[ext]
    result=json.loads(b.with_suffix('.result.json').read_text());original=json.loads(b.with_suffix('.comparison.json').read_text());total=sum(len(v) for v in checks.values())
    report=dict(fresh_all_arrays_exact=True,checked_array_count=total,checks=checks,standing_accepted=original['original_standing_assessment']['accepted'],source_model=str(F.R.MODEL),source_model_sha256=X.digest(F.R.MODEL),neural_configuration=str(R/'TIMING_CANDIDATE.json'),support=dict(upward_bw=result['mean_upward_rig_support_bw'],downward_bw=result['mean_downward_rig_force_bw'],net_bw=result['mean_net_rig_force_bw'],feet_bw=result['mean_ground_load_bw']),no_new_simulation=True)
    X.write(R/'PHYSICS_FRESH_VERIFICATION.json',report);print(json.dumps(report),flush=True);return report

def render(record,out,stem,lines,start,end):
    m=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));d=F.mujoco.MjData(m);m.vis.global_.offwidth=900;m.vis.global_.offheight=800
    plane=m.geom('ground-plane').id;m.geom_matid[plane]=-1;m.geom_rgba[plane]=[.90,.91,.92,1.]
    for body in range(m.nbody):
        name=m.body(body).name or ''
        if name.endswith(('_r','_l')):
            for g in range(m.body_geomadr[body],m.body_geomadr[body]+m.body_geomnum[body]):m.geom_rgba[g]=[.25,.46,.66,1] if name.endswith('_r') else [.88,.43,.18,1]
    z=H.load(record);assert z['qfull'].shape[1]==m.nq
    candidates=np.flatnonzero((z['t']>=start)&(z['t']<end));stride=max(1,round(.10/float(np.median(np.diff(z['t'])))));ids=candidates[::stride];frames=[]
    renderer=F.mujoco.Renderer(m,height=570,width=420);opt=F.mujoco.MjvOption();opt.flags[F.mujoco.mjtVisFlag.mjVIS_TENDON]=False;font=ImageFont.truetype(r'C:\Windows\Fonts\arial.ttf',17)
    for i in ids:
        d.qpos[:]=z['qfull'][i];F.mujoco.mj_forward(m,d);frame=Image.new('RGB',(840,660),'white');draw=ImageDraw.Draw(frame)
        for j,az in enumerate((0.,45.)):
            cam=F.mujoco.MjvCamera();cam.lookat[:]=[d.xpos[m.body('pelvis').id,0],d.xpos[m.body('pelvis').id,1],.88];cam.distance=2.8;cam.azimuth=az;cam.elevation=-5
            renderer.update_scene(d,camera=cam,scene_option=opt);frame.paste(Image.fromarray(renderer.render().copy()),(420*j,90))
        for j,line in enumerate(lines):draw.text((10,5+25*j),line,font=font,fill='black')
        draw.text((10,650-22),f'Recorded poses | 1x playback | t={z["t"][i]:.2f} s | no pose correction',font=font,fill='black');frames.append(frame)
    del renderer
    gif=out/(stem+'.gif');still=out/(stem+'.png');frames[0].save(gif,save_all=True,append_images=frames[1:],duration=100,loop=0);frames[len(frames)//2].save(still)
    return dict(record=str(record),record_sha256=X.digest(record),source_model=str(F.R.MODEL),source_model_sha256=X.digest(F.R.MODEL),gif=str(gif),still=str(still),frame_times_s=z['t'][ids].tolist(),captions=lines,rendering_only=True,no_pose_correction=True)

def main():
    fresh=verify();out=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\Review_20261101');out.mkdir(exist_ok=True)
    support=fresh['support']
    previous=json.loads((out/'GAIT2392_REVIEW_MANIFEST.json').read_text()) if '--metadata-only' in sys.argv else None
    standing=previous['standing'] if previous else render(R/'physics/syn6_central_v4_8p0_fresh.npz',out,'Gait2392_central_spiking_standing',['Gait2392 supported standing | 256 central spiking cells; 778 graded',f'Vertical rig: up {100*support["upward_bw"]:.3f}% BW; down {100*support["downward_bw"]:.3f}%; net {100*support["net_bw"]:.3f}%',f'Feet {100*support["feet_bw"]:.2f}% BW | original standing pass; fresh exact | rig retained'],.6,7.9)
    root=X.REPO/'Code/MuJoCo_SNS/spinal/hunt_training/runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/aligned_support_ladder_v1';p=root/'repeat_0p5.npz';checkpoint=json.loads((root/'CHECKPOINT.json').read_text())
    stepping=previous['stepping'] if previous else render(p,out,'Gait2392_Hunt_supported_stepping',['Gait2392 Hunt-PF supported stepping | graded controller; in place','Measured vertical support 63.73% BW | sustained cleared swings: 5 right / 4 left','Protected knees | earlier standing stage fails; lower-support clearance fails'],5.,12.8)
    stepping.update(configuration=str(root/'packaged_configuration.json'),trial_configuration=str(root/'repeat_0p5.json'),checkpoint=str(root/'CHECKPOINT.json'),result=str(root/'repeat_0p5.result.json'),limitations=checkpoint['scope'],sustained_swings=checkpoint['sustained_swings'],maximum_knee_extension_deg=checkpoint['maximum_knee_extension_deg'],knee_excursion_deg=[86.47128271184286,86.35804172995154],fresh_array_max_differences=checkpoint['fresh_repeat_max_abs_diff'])
    result=json.loads((root/'repeat_0p5.result.json').read_text());assistance=H.load(root/'repeat_0p5.assistance.npz');a=assistance['samples'];w=(a[:,0]>=5)&(a[:,0]<13);fz=(a[w,11]+a[w,2])/result['measured']['body_weight_N'];stepping['measured_support_5_13s_bw']=dict(upward=float(np.mean(np.maximum(fz,0))),downward=float(np.mean(np.maximum(-fz,0))),net=float(np.mean(fz)),spring_net=float(np.mean(a[w,11])/result['measured']['body_weight_N']),explicit_net=float(np.mean(a[w,2])/result['measured']['body_weight_N']));assert abs(float(fz.mean())-checkpoint['measured_assistance_bw'])<1e-12
    stepping['support_up_down_note']='Actual spring/damping world-Fz plus explicitly applied pelvis world-Fz, matching the5–13s scoring window; positive and negative parts of their total separately averaged. Not a nominal support knob or a spiking-case substitution.'
    standing.update(configuration=str(R/'TIMING_CANDIDATE.json'),controller_file=str(Path(H.__file__)),controller_class='TimingZeroSFA',run_file=str(Path(__file__).parent/'gait_burst_validation.py'),record=str(R/'EXPERIMENT_RECORD.md'),fresh_proof=str(R/'PHYSICS_FRESH_VERIFICATION.json'),result=str(R/'physics/syn6_central_v4_8p0_fresh.result.json'),support=fresh['support'])
    standing['duration_scope']='8 s protocol:0.6s pinned preparation plus7.4s free physical integration under the retained support apparatus; not8s harness-free integration.'
    extended=H.load(R/'physics/syn6_central_v4_8p0_fresh.apparatus_extended.npz');a=extended['samples'];bw=json.loads(Path(standing['result']).read_text())['body_weight_N'];standing['retained_rig_detail']=dict(horizontal_force_RMS_bw=float(np.sqrt(np.mean(np.sum(a[:,1:3]**2,axis=1)))/bw),pelvis_rotational_restraint_torque_RMS_Nm=float(np.sqrt(np.mean(np.sum(a[:,4:7]**2,axis=1)))),pelvis_rotational_restraint_torque_peak_Nm=float(np.max(np.linalg.norm(a[:,4:7],axis=1))),joint_parameters=json.loads(str(extended['joint_parameters_json'])))
    cp=json.loads((R/'CENTRAL_CHECKPOINT.json').read_text());cp.update(status='Versioned burst-envelope v2 neural gates accepted after documented scorer review. Actual8s biological standing passes unchanged scorer;82saved arrays fresh-exact. Original peak-scorer negative retained.',physical_checkpoint=str(R/'PHYSICAL_CHECKPOINT.json'),physical_fresh_verification=fresh,physical_trials=3,numerical_campaigns_stopped=True)
    X.write(R/'CENTRAL_CHECKPOINT.json',cp)
    record=(R/'EXPERIMENT_RECORD.md').read_text();banner='\n**Current physical checkpoint:** reviewed burst-envelope metric v2 passes independent controls; actual8s central-spiking biological Gait2392 standing passes the unchanged physical scorer. All82saved physical/neural/release arrays repeat exactly. Direct measured vertical rig forces:0.5937%BW upward,0.0440%BW downward,0.5497%BW net; feet99.3228%BW. Rotational/lateral apparatus remains. The legacy scorer rig-share estimate0.7193%BW is not substituted for measured force. No gait/harness-free claim. All campaign jobs stopped.\n'
    if '**Current physical checkpoint:**' not in record:record=record.replace('\n', '\n'+banner,1)
    record=record.replace('actual8s central-spiking biological Gait2392 standing','actual 8 s central-spiking biological Gait2392 standing protocol (0.6 s pinned preparation plus 7.4 s integration under the retained apparatus)').replace('All82saved','All 82 saved').replace('forces:0.5937%BW upward,0.0440%BW downward,0.5497%BW net; feet99.3228%BW','forces: 0.5937% BW upward, 0.0440% BW downward, 0.5497% BW net; feet 99.3228% BW')
    if '## Earlier neural-only checkpoint' not in record:record=record.replace('Machine: **easteregg2**. One timing-calibrated','## Earlier neural-only checkpoint and scoring history\n\nThe following retained sections describe the earlier neural-only stage, before the physical trials listed later.\n\nMachine: **easteregg2**. One timing-calibrated',1)
    (R/'EXPERIMENT_RECORD.md').write_text(record,encoding='utf-8')
    # Preserve the already-curated source archive byte-for-byte. This review
    # helper changes only media/metadata; its current hash is recorded separately.
    cp['source_snapshot']=dict(path=str(R/'SOURCE_SNAPSHOT.zip'),sha256=X.digest(R/'SOURCE_SNAPSHOT.zip'));X.write(R/'CENTRAL_CHECKPOINT.json',cp)
    current=dict(status='Supersedes the frozen numerical candidate historical review-pending/physical0 metadata; its parameters and hash remain unchanged.',approved_neural_metric=str(R/'BURST_METRIC_V2.json'),approved_neural_gate_pass=True,current_checkpoint=str(R/'CENTRAL_CHECKPOINT.json'),physical_checkpoint=str(R/'PHYSICAL_CHECKPOINT.json'),physical_standing_gate_pass=True,fresh_verification=str(R/'PHYSICS_FRESH_VERIFICATION.json'),fresh_exact_arrays=82,physical_duration='8s protocol:0.6s pinned preparation+7.4s integration under retained rig',all_numerical_jobs_stopped=True)
    manifest=dict(machine='easteregg2',current_validation=current,source_frame='Legacy converted Gait z-up with explicit frame adapter; actual human Gait2392, not W2L or CAD.',standing=standing,stepping=stepping,review_helper=dict(path=str(Path(__file__)),sha256=X.digest(__file__),source_snapshot_preserved=True),scope='Two different controllers and experiments; not a matched training comparison. Saved-pose rendering only; all numerical campaign jobs stopped.')
    cp['current_validation']=current;X.write(R/'CENTRAL_CHECKPOINT.json',cp)
    X.write(out/'GAIT2392_REVIEW_MANIFEST.json',manifest)
    text='# Gait2392 review media\n\nMachine: **easteregg2**. These are two different controllers and experiments, both actual human Gait2392. No new simulation was run to make these media.\n\n'
    for label,row in [('Central-spiking standing',standing),('Hunt-PF supported stepping',stepping)]:
        text+='## '+label+'\n\n'
        for key in ('gif','still','source_model','configuration','record','result','checkpoint','fresh_proof'):
            if key in row:text+='- '+key+': `'+str(row[key])+'`\n'
        text+='\n'+'; '.join(row['captions'])+'.\n\n'
    text+='Stepping knee excursions are approximately86.47/86.36 degrees, near the protected90-degree flexion range. Its knees do not hyperextend. The selected trial has5/4 sustained whole-mesh-cleared swings; lower-support aligned cases fail clearance, and standing has not passed for that controller. Earlier43.61% support results used a penetrating initial pose and are excluded. This is supported stepping in place, not completed curriculum or harness-free locomotion.\n\nMeasured stepping vertical harness split (5–13s, BW fractions): '+json.dumps(stepping['measured_support_5_13s_bw'])+'.\n\nStanding uses the exact four neural parameters in TIMING_CANDIDATE.json and the accepted syn6 source three gains. Replay source fixture is bilateral in-phase, not walking. Actual8s standing passes the unchanged native scorer and all82 saved arrays match fresh. Old local-peak negative is retained; approved versioned burst-envelope validation is separate. No force allocator or MN override was added. The frozen numerical JSON retains its historical review-pending annotation; BURST_METRIC_V2.json and CENTRAL_CHECKPOINT.json carry the subsequent acceptance and actual physical proof.\n'
    text+='\nStanding duration: 8 s requested protocol, comprising 0.6 s pinned preparation and 7.4 s subsequent integration under the retained rig. Horizontal rig-force RMS is 2.635% BW; pelvis rotational-restraint torque RMS is 5.02 N m (peak 26.96 N m). Small vertical assistance does not mean harness-free balance.\n'
    text=text.replace('approximately86.47/86.36','approximately 86.47/86.36').replace('protected90-degree','protected 90-degree').replace('has5/4','has 5/4').replace('Earlier43.61%','Earlier 43.61%').replace('Actual8s','The 8 s protocol').replace('all82 saved','all 82 saved')
    (out/'GAIT2392_REVIEW.md').write_text(text,encoding='utf-8')
    X.write(R/'FILE_MANIFEST.json',dict(files=[dict(path=str(q),bytes=q.stat().st_size,sha256=X.digest(q)) for q in sorted(R.iterdir()) if q.is_file() and q.name!='FILE_MANIFEST.json']))
    print(json.dumps(dict(review=str(out/'GAIT2392_REVIEW.md'),manifest=str(out/'GAIT2392_REVIEW_MANIFEST.json'),standing_gif=standing['gif'],stepping_gif=stepping['gif'])),flush=True)

if __name__=='__main__':main()
