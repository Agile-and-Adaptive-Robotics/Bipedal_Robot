"""Isolated geometry sensitivity and w2lvar timing probe, not a gait trial.

The proximal foot translation is a hypothetical mount offset. It does not
shorten tibial meshes, change their inertias, or establish a hardware design.
The inherited MN clamp remains active: w2lvar affects RG timing and sensory
central processing, but its direct PF/MN contributions are canceled.
"""
import os,sys,json,argparse,hashlib
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET
HERE=Path(__file__).resolve().parent
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import numpy as np
import mujoco

def prepare():
    source=HERE/'ah001_abduction_diagnostic_v2.xml'
    tree=ET.parse(source);root=tree.getroot();records={}
    model=mujoco.MjModel.from_xml_path(str(source));data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    for side in ('r','l'):
        node=root.find(".//body[@name='foot_"+side+"']")
        original=np.fromstring(node.get('pos'),sep=' ')
        direction=original/np.linalg.norm(original)
        shift=-.060*direction
        node.set('pos',' '.join('%.15g'%x for x in original+shift))
        foot=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'foot_'+side)
        shank=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'shank_'+side)
        records[side]=dict(original_foot_offset_in_shank_m=original.tolist(),
            proximal_translation_in_shank_m=shift.tolist(),
            shank_origin_to_mount_m=float(np.linalg.norm(original)),
            old_mount_to_sole_y_m=float(data.xpos[foot,1]-min(
                (model.mesh_vert[model.mesh_vertadr[mid]:model.mesh_vertadr[mid]+model.mesh_vertnum[mid]]@data.geom_xmat[g].reshape(3,3)[1]+data.geom_xpos[g,1]).min()
                for g in range(model.ngeom) if model.geom_bodyid[g]==foot
                for mid in [model.geom_dataid[g]] if mid>=0)))
    path=HERE/'ah001_foot_proximal60_hypothesis.xml';tree.write(path,encoding='utf-8',xml_declaration=True)
    modified=mujoco.MjModel.from_xml_path(str(path))
    assert np.array_equal(model.body_mass,modified.body_mass)
    assert np.array_equal(model.body_inertia,modified.body_inertia)
    assert model.nu==modified.nu==54
    result=dict(machine='easteregg2',source=str(source),output=str(path),sides=records,
        mass_kg=float(modified.body_mass.sum()),translation_m=.060,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        output_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        scope='Hypothetical proximal foot-mount translation, not tibia shortening or a validated mating datum. Foot meshes, local inertia, tibial meshes, knee linkage, source route mapping, and total mass unchanged; foot-attached points follow the moved foot.')
    (HERE/'foot_mount_sensitivity_spec.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

def run(network,shifted):
    import run_neural_robot as runner
    settings=json.loads((HERE/'runs/abd_defense_exact.json').read_text())['run_options']
    label=('w2lvar' if network=='w2lvar' else 's3k')+('_foot60' if shifted else '_foot0')
    settings.update(output=HERE/'runs'/label,
        model=HERE/('ah001_foot_proximal60_hypothesis.xml' if shifted else 'ah001_abduction_diagnostic_v2.xml'))
    capture={};original=runner.bn.build
    def build(acts,dt,interleg):
        if network=='w2lvar':
            import build_network_w2lvar as variant
            net=variant.build(acts,dt=dt,interleg=interleg)
            capture['variant_gains']=dict(variant.GAINS)
            capture['variant_source_sha256']=hashlib.sha256(Path(variant.__file__).read_bytes()).hexdigest()
        else:net=original(acts,dt=dt,interleg=interleg)
        capture.update(network=network,true_neuron_count=int(net.compiled.V.size),
            inputs=len(net.inputs),pf_names=[name for name in net.idx if name.startswith('PF_')],
            named_index_count=len(net.idx))
        return net
    runner.bn.build=build
    summary=runner.run(SimpleNamespace(**settings))
    summary.update(experimental_network_probe=capture,
        probe_adapter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='Hybrid timing comparison only: motor-target clamp cancels incoming PF/reflex MN synaptic currents. Not a test of full w2lvar neural muscle-pattern control.',
        controller_network=network,foot_mount_translation_m=.060 if shifted else 0.)
    settings['output'].with_suffix('.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--network',choices=['s3k','w2lvar'],default='s3k');p.add_argument('--shifted',action='store_true');args=p.parse_args()
    if args.prepare:prepare()
    else:run(args.network,args.shifted)
