"""Source-specific mechanical import; unsupported features are reported, not dropped."""
import os,sys,json,importlib.util,shutil,xml.etree.ElementTree as E
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1]
SOURCE=HERE/'Biarticular_net_Standalone.asim'
root=E.parse(SOURCE).getroot()
types=sorted({e.findtext('Type') for e in root.iter('RigidBody')});jtypes=sorted({e.findtext('Type') for e in root.iter('Joint')})
unsupported=[x for x in types if x not in ('Box','Plane','Attachment','LinearHillMuscle','LinearHillStretchReceptor','Spring')]
unsupported_joints=[x for x in jtypes if x not in ('Hinge','Prismatic')]
out=HERE/'build';out.mkdir(exist_ok=True)
preflight=dict(source=str(SOURCE),body_types=types,joint_types=jtypes,unsupported_bodies=unsupported,unsupported_joints=unsupported_joints,
    length_units_note='Serialized physical lengths match .aproj Actual SI values, despite the display distance unit Decimeters; mass values are grams.',
    frozen_source_root=root.find('.//Organism/RigidBody/Freeze').text=='True',
    scope='Mechanical import only. Native linear-Hill and neural signals require separate coupled runtime validation.')
(out/'PREFLIGHT.json').write_text(json.dumps(preflight,indent=2))
if unsupported or unsupported_joints:raise RuntimeError('Unsupported source features: '+str(unsupported+unsupported_joints))
spec=importlib.util.spec_from_file_location('rat_source_importer',REPO/'Code/MuJoCo_SNS/spinal/virtual_walker_20261031/import_virtual_walker.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);mod.HERE=out
path=mod.build(source_override=SOURCE)
# Importer's W2L root is free by design; preserve this rat source's actual freeze.
tree=E.parse(path);mj=tree.getroot();body=mj.find('worldbody/body')
if preflight['frozen_source_root']:
 fj=body.find('freejoint')
 if fj is not None:body.remove(fj)
E.indent(tree);final=HERE/'Deng_Biarticular_Rat.xml';tree.write(final,encoding='utf-8',xml_declaration=True)
import mujoco
model=mujoco.MjModel.from_xml_path(str(final));data=mujoco.MjData(model);mujoco.mj_forward(model,data)
preflight.update(xml=str(final),compiled=True,nq=model.nq,nv=model.nv,nu=model.nu,source_frozen_root_preserved=True,
    deviations=['Physics target 2 ms; source 50 microseconds, so convergence is not yet established.',
                'Contact law and collision handling use explicit MuJoCo assumptions; source neural/force trajectory equivalence not claimed.',
                'Imported from source muscles/routes; no third muscle added or dropped.'])
(HERE/'MUJOCO_IMPORT_STATUS.json').write_text(json.dumps(preflight,indent=2));print(json.dumps(preflight,indent=2))
