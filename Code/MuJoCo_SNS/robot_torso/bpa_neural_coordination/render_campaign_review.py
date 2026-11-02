"""Render saved scientific trajectories only; never integrate physics."""
import os
import sys
from pathlib import Path
os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)
import argparse
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np
import mujoco
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
ROBOT = HERE.parent
sys.path.insert(0, str(ROBOT))
from dynamic_routes import DynamicRoutes
from robot_kinematics import camera_axes
CAMPAIGN = Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031')
CAD_RUN = CAMPAIGN/'cad_bpa/pair_neural_recruitment/normalized_unclamped/support085'
OUT = CAMPAIGN/'review_20261101/cad'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cad(fps=10, metadata_only=False):
    OUT.mkdir(parents=True, exist_ok=True)
    case = CAD_RUN/'mixed/mixed'
    assessment = json.loads((CAD_RUN/'mixed/ASSESSMENT.json').read_text())
    protocol = json.loads((case/'LIVE_PREREGISTRATION.json').read_text())
    source_summary = json.loads((case/'cad.json').read_text())
    config = CAD_RUN/'mixed/ASSUMED_PAIR_CONFIGURATION.json'
    route_config = ROBOT/'robotbody_routes.json'
    pair = ('vas_med_r', 'vas_med_l')
    with np.load(case/'ACTUATOR_BRIDGE.npz') as archive:
        keys = ('t', 'qpos_full', 'qvel_full', 'pressure', 'neural_output', 'rig')
        trace = {key: archive[key] for key in keys}
    source = Path(protocol['arguments']['model'])
    original = mujoco.MjModel.from_xml_path(str(source))
    tree = ET.parse(source)
    root = tree.getroot()
    root.find('compiler').set('meshdir', str(ROBOT/'meshes'))
    visual = root.find('visual/global')
    visual.set('offwidth', '1600'); visual.set('offheight', '1600')
    ET.SubElement(root.find('asset'), 'texture', name='review_sky', type='skybox', builtin='gradient',
                  rgb1='1 1 1', rgb2='.93 .95 .98', width='512', height='3072')
    root.find("worldbody/geom[@name='ground']").set('rgba', '.94 .95 .96 1')
    for name, position in [('Anterior', [2.4, .95, 0]), ('Oblique', [1.8, 1.05, 1.8])]:
        position = np.asarray(position, float)
        axes = camera_axes(position, np.array([-.05, .66, 0.]))
        ET.SubElement(root.find('worldbody'), 'camera', name='review_'+name,
                      pos=' '.join(map(str, position)),
                      xyaxes=' '.join(map(str, np.r_[axes[:,0], axes[:,1]])), fovy='37')
    model = mujoco.MjModel.from_xml_string(ET.tostring(root, encoding='unicode'))
    for attribute in ('body_mass','body_inertia','qpos0','dof_damping','eq_solref','eq_solimp'):
        assert np.array_equal(getattr(model,attribute), getattr(original,attribute)), attribute
    # Match the saved mixed-actuation mode; kinematic replay supplies no fictitious native force.
    model.actuator_gainprm[:]=0; model.actuator_biasprm[:]=0
    model.actuator_gaintype[:]=0; model.actuator_biastype[:]=0
    data = mujoco.MjData(model)
    dynamic = DynamicRoutes(model, json.loads(route_config.read_text()))
    for gid in range(model.ngeom):
        name = model.body(int(model.geom_bodyid[gid])).name or ''
        if name.endswith('_r'): model.geom_rgba[gid] = [.267,.467,.667,1]
        elif name.endswith('_l'): model.geom_rgba[gid] = [.80,.56,.25,1]
    for aid in range(model.nu):
        tid = int(model.actuator_trnid[aid,0])
        highlight = model.actuator(aid).name in pair
        model.tendon_rgba[tid] = [.05,.57,.29,.95] if highlight else [.45,.45,.48,.13]
        model.tendon_width[tid] = .0026 if highlight else .0008
    renderers = {}
    arial = 'C:/Windows/Fonts/arial.ttf'
    mean_support = 100*assessment['metrics']['mean_upward_bw']
    def render(index, scale=1):
        data.qpos[:] = trace['qpos_full'][index]
        data.qvel[:] = trace['qvel_full'][index]
        data.time = float(trace['t'][index])
        dynamic.update(data)  # Exact source moving/conditional route geometry; no mj_step.
        if scale not in renderers:
            renderers[scale] = mujoco.Renderer(model, width=600*scale, height=670*scale)
        renderer = renderers[scale]
        canvas = Image.new('RGB', (1200*scale, 900*scale), 'white')
        draw = ImageDraw.Draw(canvas)
        font = lambda size: ImageFont.truetype(arial, size*scale)
        draw.text((20*scale,12*scale), 'Supported CAD standing: pressure-actuated neural pair', font=font(27), fill='#222222')
        draw.text((20*scale,49*scale), '2 modeled BPAs + 52 hybrid reference-tension routes', font=font(24), fill='#444444')
        for column, name in enumerate(('Anterior','Oblique')):
            renderer.update_scene(data, 'review_'+name)
            canvas.paste(Image.fromarray(renderer.render()), (600*scale*column, 105*scale))
            draw.text(((20+600*column)*scale,84*scale), ('A  ' if column==0 else 'B  ')+name+' view', font=font(23), fill='#333333')
        pressure = trace['pressure'][index]
        draw.text((20*scale,791*scale), f't = {data.time:4.2f} s   |   Right / left BPA pressure: {pressure[0]:.0f} / {pressure[1]:.0f} kPa', font=font(23), fill='#333333')
        draw.text((20*scale,831*scale), f'Mean upward harness support: {mean_support:.2f}% BW  |  Lateral and rotational guides', font=font(22), fill='#333333')
        draw.text((20*scale,867*scale), 'Green: actual vastus-medialis routes   |   Saved 8 s supported-standing record', font=font(20), fill='#3D604A')
        return canvas
    times = trace['t']
    selected = [int(np.argmin(abs(times-t))) for t in np.arange(0.,8.,1/fps)]
    still_index = int(np.argmin(abs(times-7.5)))
    still = OUT/'CAD_BPA_PAIR_SUPPORT085_still.png'
    gif = OUT/'CAD_BPA_PAIR_SUPPORT085.gif'
    if not metadata_only:
        render(still_index,2).save(still,dpi=(300,300))
        frames=[]
        for number,index in enumerate(selected):
            frames.append(render(index).convert('P',palette=Image.Palette.ADAPTIVE,colors=192))
            if number%20==0: print(f'Rendered {number}/{len(selected)} saved frames',flush=True)
        ticks=np.rint(np.arange(len(frames)+1)*100/fps).astype(int)
        frames[0].save(gif,save_all=True,append_images=frames[1:],duration=(10*np.diff(ticks)).tolist(),loop=0,optimize=True)
    else:
        assert gif.exists() and still.exists()
        with Image.open(gif) as existing: assert existing.n_frames==len(selected)
    renderers.clear()  # MuJoCo2.3.7 Renderer has no public close() method.
    manifest = dict(machine='easteregg2', replay_only=True, physics_steps_executed=0,
                    gif=str(gif), still=str(still), base_model=str(source),
                    run_summary=str(case/'cad.json'), runtime_protocol=str(case/'LIVE_PREREGISTRATION.json'),
                    exact_consumed_arguments=protocol['arguments'], assumed_BPA_configuration=str(config),
                    initial_qpos_archive=protocol['initial_pose_source'], initial_pose_index=1069,
                    controller_entry=str(HERE/'cad_pair_support_step.py'), controller_case_arguments=['--case','mixed'],
                    pair_normalization_entry=str(HERE/'cad_pair_unclamped.py'),
                    pair_recruitment_config=str(CAD_RUN.parent/'RECRUITMENT_PREREGISTRATION.json'),
                    network_builder=str(ROBOT.parent/'spinal/build_network.py'),
                    network_neurons=source_summary['neuron_count'], neural_inputs=source_summary['neural_inputs'],
                    controller='Original graded 262-cell SNS with six mirrored PF/IBEXC output weights normalized only for vas_med_r/l; selected pair receives literal SNS currents with POST clamp zero. Other52 retain explicit external hybrid allocator and 5% reference-force floor.',
                    actual_pressure_routes=list(pair), other_routes='52 filtered ideal reference-tension channels; not BPAs',
                    route_configuration=str(route_config), dynamic_route_runtime=str(ROBOT/'dynamic_routes.py'),
                    measured_metrics=assessment['metrics'], original_assessment=str(CAD_RUN/'mixed/ASSESSMENT.json'),
                    fresh_exact_repeat_record=str(CAD_RUN/'FRESH_REPEAT_AUDIT.json'),
                    main_record=str(CAD_RUN/'EXPERIMENT_RECORD.md'),
                    limitations=['Historical accepted 8s settling condition, not endurance or independent balance.',
                                 'Measured upward support87.34%BW plus AP/lateral springs and rotational guidance.',
                                 'Only2 assumed virtual pressureBPAs; dimensions/pressure dynamics are not hardware calibrated.',
                                 'Legacy toe-out about+29/-30deg, imported loop gaps11.63/11.70mm; loaded tail about5.52/6.50mm.',
                                 'Legs visibly cross/overlap in frontal view by7.5s; historical settling acceptance did not require uncrossed stance.',
                                 'Posterior COM settlement and substantial apparatus moment dependence; later mechanical audit does not validate this as closed-fourbar fidelity.',
                                 'Legacy generic scope strings in assumed-config/LIVE legacy-description fields predate free-pair override; numerical assumptions preserved, actual selectedpair freeSNS verified by PAIR_NEURAL_LOG and assessment.'],
                    posture_audit=str(CAD_RUN.parent/'posture_balance_audit/POSTURE_GEOMETRY_AUDIT.json'),
                    source_trace=str(case/'ACTUATOR_BRIDGE.npz'), frame_time_range_s=[float(times[selected[0]]),float(times[selected[-1]])],
                    still_time_s=float(times[still_index]), fps=fps, frame_count=len(selected),
                    metadata_recovery_only=metadata_only,
                    anatomical_frame='x anterior, y superior, z right',
                    visual_only_changes='Cameras, pale background, bilateral colors, pairroute highlighting, annotation. Saved qpos unchanged; no stepping or pressure integration.',
                    hashes={str(path):digest(path) for path in (source,route_config,config,case/'ACTUATOR_BRIDGE.npz',gif,still,Path(__file__))})
    (OUT/'REPRODUCTION_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    flags = json.dumps(protocol['arguments'],indent=2)
    text = f'''# CAD pressure-pair review assets

Machine: **easteregg2**. Rendered from an existing saved 8 s trajectory; no new simulation or parameter search.

- GIF: `{gif}`
- Still (2400 × 1800 px, 300 dpi): `{still}`
- Complete model/configuration/flag/provenance manifest: `{OUT/'REPRODUCTION_MANIFEST.json'}`

**Model and test.** Original AH001 CAD meshes with split payload, 75 kg total model mass, ball-and-socket lumbar spine, moving/conditional muscle paths, and two modeled pressure BPAs (`vas_med_r`, `vas_med_l`). Their six mirrored SNS recruitment weights are normalized with clamp input zero. The other 52 routes remain explicitly hybrid ideal-tension actuators. This is a supported-standing settling demonstration, not full-BPA actuation, locomotion, independent balance, calibrated hardware, or endurance.

Base XML: `{source}`. This file alone does not reproduce the pressure-pair experiment. Required pressure configuration: `{config}`. Dynamic route configuration: `{route_config}`. Initial pose: `{protocol['initial_pose_source']}`, qpos row 1069, initial velocity zero. Recruitment configuration: `{CAD_RUN.parent/'RECRUITMENT_PREREGISTRATION.json'}`. Runtime opt-in entry: `{HERE/'cad_pair_support_step.py'}` with `--case mixed`, which applies the documented pair and pressure overrides; exact historical source snapshots are retained in the experiment folders. The consumed flags below are authoritative alongside those opt-in modules.

**Support.** Prescribed vertical feedforward is85%BW; actual mean upward support is{mean_support:.3f}%BW after the unchanged800N/m height spring and120Ns/m damper. There are also AP/lateral guides and angular guidance (300Nm/rad,60Nms/rad). Mean foot load is{100*assessment['metrics']['mean_ground_load_bw']:.3f}%BW, and peak pair pressure is{assessment['metrics']['pressure_peak_kpa']:.3f}kPa. Full signed measured apparatus forces and torques remain in the source trajectory and experiment record.

**Limitations for editors.** The legs visibly cross/overlap in frontal view by7.5s; the historical settling acceptance did not gate uncrossed stance. Preserve that limitation alongside toe-out, posterior COM settlement, legacy four-bar loop gaps, ideal52 context, and uncalibrated pressure assumptions in captions/discussion. Do not relabel this historical eight-second accepted settling run as the later near-closed rigid-link initialization. Its imported gap was about11.6/11.7mm, with tail gaps about5.5/6.5mm. The later mechanical audit remains separate. Green traces are the actual named virtual muscle routes, not visually fabricated BPA cylinders or confirmed hardware mountings. Some inherited metadata scope strings precede the free-pair override; current controller identification comes from the actual neural-current and assessment audits.

Experiment record: `{CAD_RUN/'EXPERIMENT_RECORD.md'}`. Native assessment: `{CAD_RUN/'mixed/ASSESSMENT.json'}`. Exact fresh repeat: `{CAD_RUN/'FRESH_REPEAT_AUDIT.json'}`. Detailed posture/contact audit: `{CAD_RUN.parent/'posture_balance_audit/EXPERIMENT_RECORD.md'}`.

Render-only reproduction: `D:\\Anaconda\\envs\\myo\\python.exe {Path(__file__)} --cad`. The renderer uses saved qpos with source moving/conditional route updates; it never calls `mj_step`.

## Exact consumed historical runtime flags

```json
{flags}
```
'''
    (OUT/'README.md').write_text(text,encoding='utf-8')
    print(gif,still,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--cad',action='store_true');parser.add_argument('--fps',type=int,default=10);parser.add_argument('--metadata-only',action='store_true')
    args=parser.parse_args()
    if args.cad:cad(args.fps,args.metadata_only)
