"""Frozen numericalvariant validation on bothfixture sides/signedknownloads."""
import json,copy,xml.etree.ElementTree as ET
import numpy as np
from scipy.linalg import null_space
import cad_fourbar_fixture as F
OUT=F.OUT;mujoco=F.mujoco

class Left(F.Fixture):
 def __init__(self,text,tc=None,imp=None,dt=.002,iterations=100):
  self.m=mujoco.MjModel.from_xml_string(text);self.m.opt.timestep=dt;self.m.opt.iterations=iterations
  if tc is not None:self.m.eq_solref[0]=[tc,1.]
  if imp is not None:self.m.eq_solimp[0,:2]=imp
  self.d=mujoco.MjData(self.m);self.c=np.array([self.m.jnt_axis[self.m.joint('back_l_joint').id,2],self.m.jnt_axis[self.m.joint('shank_l_joint').id,2],0.]);self.W=null_space(self.c[None,:]);self.q0=self.d.qpos.copy();self.base();self.J=self.d.efc_J.reshape(self.d.nefc,self.m.nv).copy();self.rows=np.flatnonzero(self.d.efc_type==int(mujoco.mjtConstraint.mjCNSTR_EQUALITY));self.sourceK=self.localK()

def left_text():
 root=ET.fromstring((OUT/'PINNED_SOURCE_FIXTURE.xml').read_text());original=ET.parse(F.C.ROBOT/'ah001_neural.xml').getroot();left={node.get('name'):node for node in original.iter('joint')}
 for parent in root.iter():
  for child in list(parent):
   if child.tag=='joint':parent.remove(child)
 for part in ('back','shank','front'):body=next(node for node in root.iter('body') if node.get('name')==part+'_l');body.insert(0,copy.deepcopy(left[part+'_l_joint']))
 eq=root.find('equality');eq.clear();eq.append(copy.deepcopy(next(node for node in original.find('equality') if node.get('name')=='knee_loop_l')));root.set('model','Pinned source-inertia leftfourbar fixture');text=ET.tostring(root,encoding='unicode');(OUT/'PINNED_SOURCE_FIXTURE_LEFT.xml').write_text(text);return text

def main():
 sel=json.loads((OUT/'DERIVED_VARIANT.json').read_text());texts={'right':(OUT/'PINNED_SOURCE_FIXTURE.xml').read_text(),'left':left_text()};out={}
 for side,text in texts.items():
  cls=F.Fixture if side=='right' else Left;source=cls(text);selected=cls(text,sel['timeconst_s'],sel['constant_impedance'],sel['timestep_s'],500);cases=[]
  for load in (-9.,-8.,0.,8.,9.):
   s,_,_=source.equilibrium(load);v,_,_=selected.equilibrium(load);cases.append(dict(load_Nm=load,source=s,selected=v))
  out[side]=dict(cases=cases,source_details=source.details(),selected_details=selected.details(),mass_inertia_same=bool(np.array_equal(source.m.body_mass,selected.m.body_mass) and np.array_equal(source.m.body_inertia,selected.m.body_inertia)))
 result=dict(scope='Frozenvariant signedload/bilateralfixture validation, no parameterrefit orrobottrial',variant=sel,cases=out,gate=all(row['selected']['actual_gap_m']<=.0001 and row['selected']['qacc_max']<1e-5 for side in out.values() for row in side['cases']));(OUT/'SIGNED_BILATERAL_CHECK.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(gate=result['gate'],cases={side:[dict(load=r['load_Nm'],sourcegap_mm=r['source']['actual_gap_m']*1000,selectedgap_um=r['selected']['actual_gap_m']*1e6,selectedacc=r['selected']['qacc_max']) for r in value['cases']] for side,value in out.items()}),indent=2))

if __name__=='__main__':main()
