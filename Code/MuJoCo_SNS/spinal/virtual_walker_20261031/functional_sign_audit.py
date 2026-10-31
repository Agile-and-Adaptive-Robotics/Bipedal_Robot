"""Trace PF→MN→actual muscle ID→signed anatomical torque on both sides."""
import json
import numpy as np
from import_virtual_walker import build,mujoco,safe
from native_kernel import NativeController
def audit(newer=False,source=None):
 path=build(newer,source);m=mujoco.MjModel.from_xml_path(str(path));d=mujoco.MjData(m);mujoco.mj_forward(m,d);meta=json.loads((path.parent/'IMPORT_AUDIT.json').read_text());bank=meta['muscles'];net=NativeController(meta['source_snapshot'],bank,contact_si=True);rows=[]
 for pre,mi,coeff in net.outputs:
  b=bank[mi];muscle=b['name'];side='L' if '_L_' in muscle else 'R';joint=muscle.split('_')[0]+'_'+side;qi=m.jnt_qposadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,joint)];ti=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+safe(muscle));d.qpos[qi]=1e-5;mujoco.mj_forward(m,d);plus=d.ten_length[ti];d.qpos[qi]=-1e-5;mujoco.mj_forward(m,d);minus=d.ten_length[ti];d.qpos[qi]=0;arm=float(-(plus-minus)/2e-5);name=net.names[pre];expected=1 if ' flx' in name else -1;slope=3*coeff[0]*(-.06)**2+2*coeff[1]*(-.06)+coeff[2]
  incoming=[dict(source=net.names[int(e[0])],source_guid=next(k for k,v in net.ids.items() if v==int(e[0])),gmax_uS=float(e[2]),reversal_mV=float(e[5])) for e in net.edges if int(e[1])==pre and ' PF ' in net.names[int(e[0])]]
  rows.append(dict(mn=name,mn_guid=next(k for k,v in net.ids.items() if v==pre),physical_muscle_name=muscle,physical_muscle_guid=b['id'],adapter_polynomial_ABCD=coeff.tolist(),adapter_slope_at_minus60mV=float(slope),joint=joint,moment_arm_m=arm,torque_per_increased_MN_voltage_sign=int(np.sign(slope*arm)),label_expected_torque_sign=expected,label_matches_function=bool(np.sign(slope*arm)==expected),pf_inputs=incoming))
 out=dict(source=meta['source'],source_sha256=meta['source_sha256'],convention='Positive source-relative joint coordinate increases anatomical hip/knee flexion and ankle dorsiflexion, checked against native moving-marker geometry; tensile force moment arm = −dL/dq.',rows=rows,mismatches=[r['mn'] for r in rows if not r['label_matches_function']],scope='An actuator-name mismatch is distinguished from actual end-to-end MN wiring. No route, neuron signal, or weight is changed by this audit.')
 (path.parent/'FUNCTIONAL_SIGN_AUDIT.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));return out
if __name__=='__main__':audit()
