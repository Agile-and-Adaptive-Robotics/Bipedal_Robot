"""Fresh-repeat the slow candidate and render measured direction arrows."""
import argparse,json,contextlib
from pathlib import Path
import numpy as np
import mujoco
from run_virtual_walker import run
from shuffle_optimization_20261104 import measure
from render_trial import render
def main():
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();best=json.loads((a.root/'SUMMARY.json').read_text())['best'];cfg=best['diagnostics']['configuration']
 with(a.root/'fresh_repeat.log').open('w',encoding='utf-8')as log,contextlib.redirect_stdout(log):
  run(newer=True,duration=16.,tag='repeat',output_root=a.root/'fresh_repeat',**cfg);folder=a.root/'fresh_repeat/W2L_VirtualWalker_Newer/repeat';result=measure(folder);old=np.load(Path(best['diagnostics']['folder'])/'trajectory.npz');z=np.load(folder/'trajectory.npz');parity={key:old[key].shape==z[key].shape and old[key].tobytes()==z[key].tobytes()for key in old.files if old[key].dtype.kind in 'fiu'}
  assert all(parity.values())and result['gates']['knees']and result['gates']['complete'];m=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));d=mujoco.MjData(m);d.qpos[:]=z['q'][0];mujoco.mj_forward(m,d);joint=m.joint('Treadmill_Slide');qdot=np.gradient(z['q'][:,int(joint.qposadr[0])],z['t']);keep=z['t']>=2.;direction=dict(toe_minus_foot_world={s:(d.xpos[m.body('toe_'+s).id]-d.xpos[m.body('foot_'+s).id]).tolist()for s in('L','R')},mean_belt_velocity_world=(d.xaxis[joint.id]*qdot[keep].mean()).tolist(),scope='Green arrow = model-forward +X. Blue arrow follows measured belt slide direction; arrow length does not encode speed.')
  proof=dict(selected_trial=best['number'],configuration=cfg,all_numeric_arrays_byte_exact=all(parity.values()),arrays=parity,assessment=result,belt_direction=direction,scope='Partial slow-stepping candidate; full short-shuffle gate is not passed');(a.root/'FRESH_VERIFICATION.json').write_text(json.dumps(proof,indent=2));(a.root/'BELT_DIRECTION.json').write_text(json.dumps(direction,indent=2))
  render(folder,6.,12.,title='W2L108 | slow stepping; shuffle target not yet met',guidance_caption='Green: forward (+X); blue: measured belt motion',detail_caption=f"Upward assistance {100*result['upward_assistance_BW']:.1f}% BW | {sum(x['clear_unloaded_swings']for x in result['legs'].values())} clear swings / 16 s",direction_arrows=True)
 print(json.dumps(dict(repeated=all(parity.values()),shuffle_pass=result['passed'],assessment={k:result[k]for k in('upward_assistance_BW','sustained_double_support_fraction','airborne_fraction')},folder=str(folder))))
if __name__=='__main__':main()
