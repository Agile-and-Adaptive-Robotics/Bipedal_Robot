"""Create an explicit, separate collision/stop diagnostic; preserve baseline."""
import json,xml.etree.ElementTree as ET
from run_neural_robot import *
source=HERE/'ah001_neural.xml';destination=HERE/'ah001_abduction_diagnostic_v2.xml'
tree=ET.parse(source);root=tree.getroot()
base=mujoco.MjModel.from_xml_path(str(source))
changed=[]
for loop in root.findall('equality/connect'):
    loop.set('solref','.004 1');loop.set('solimp','.999 .9999 .0001')
for body in root.findall('.//body'):
    name=body.get('name','')
    if not name.endswith(('_r','_l')):continue
    side=name[-1];own,other=(4,8) if side=='r' else (8,4)
    for geom in body.findall('geom'):
        geom.set('contype',str(int(geom.get('contype','0'))|own))
        geom.set('conaffinity',str(int(geom.get('conaffinity','0'))|other))
        changed.append(geom.get('name'))
for side in ('r','l'):
    stop=ET.SubElement(root.find('tendon'),'fixed',name='relative_knee_stop_'+side,
        limited='true',range='-1.2217304763960306 0.05235987755982989',solreflimit='.02 1',solimplimit='.9 .95 .001')
    for part in ('back','shank'):
        joint=part+'_'+side+'_joint';jid=mujoco.mj_name2id(base,mujoco.mjtObj.mjOBJ_JOINT,joint)
        ET.SubElement(stop,'joint',joint=joint,coef=str(float(base.jnt_axis[jid,2])))
ET.indent(tree,space='  ');tree.write(destination,encoding='utf-8',xml_declaration=True)
model=mujoco.MjModel.from_xml_path(str(destination))
assert np.array_equal(model.body_mass,base.body_mass)
assert model.nu==base.nu==54 and model.nq==base.nq and model.nv==base.nv
assert np.array_equal(model.actuator_gainprm,base.actuator_gainprm)
audit=dict(source=str(source),output=str(destination),mass_kg=float(model.body_mass.sum()),
    modified_interleg_collision_geoms=len(changed),actuators=model.nu,
    relative_knee_stop_deg=[-70.,3.],
    assumptions='Adjustable diagnostic passive stop surrogates relative to the CAD reference; not measured hardware stops. Opposite-leg CAD mesh contacts enabled; ipsilateral joint overlaps remain excluded by collision masks.',
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),output_sha256=hashlib.sha256(destination.read_bytes()).hexdigest())
(HERE/'abduction_diagnostic_model_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))
