"""Inspect recorded robot motion using the fixed OpenSim/ISB y-up camera."""
import argparse,json,time
from pathlib import Path
import mujoco
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from dynamic_routes import DynamicRoutes

HERE=Path(__file__).resolve().parent
def main(args):
    summary=json.loads(args.run.with_suffix('.json').read_text())
    model=mujoco.MjModel.from_xml_path(summary.get('model_path',str(HERE/'ah001_neural.xml')))
    config=json.loads((HERE/'robotbody_routes.json').read_text())
    dynamic=DynamicRoutes(model,config);data=mujoco.MjData(model)
    trace=np.load(args.run)
    summary=json.loads(args.run.with_suffix('.json').read_text())
    floor=100*summary.get('run_options',{}).get('activation_floor',0.)
    shifted=bool(summary.get('run_options',{}).get('lateral_shift',0.))
    camera=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_CAMERA,'isb_preview')
    if args.view:
        from mujoco import viewer as mj_viewer
        with mj_viewer.launch_passive(model,data) as viewer:
            viewer.cam.type=mujoco.mjtCamera.mjCAMERA_FIXED;viewer.cam.fixedcamid=camera
            index=0
            while viewer.is_running():
                begin=time.perf_counter()
                data.qpos[:]=trace['qpos'][index];data.act[:]=trace['activation'][index]
                dynamic.update(data);viewer.sync()
                index=(index+1)%len(trace['time'])
                time.sleep(max(0.,.01-(time.perf_counter()-begin)))
        return
    renderer=mujoco.Renderer(model,height=660,width=800)
    font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',19)
    small=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',16)
    def frame(index):
        data.qpos[:]=trace['qpos'][index];data.act[:]=trace['activation'][index]
        dynamic.update(data)
        for i,tid in enumerate(model.actuator_trnid[:,0]):
            activation=float(np.clip(trace['activation'][index,i],0,1))
            model.tendon_rgba[tid]=[.7+.25*activation,.35-.2*activation,.3-.2*activation,.25+.65*activation]
            model.tendon_width[tid]=.0015
        renderer.update_scene(data,camera)
        picture=Image.fromarray(renderer.render())
        draw=ImageDraw.Draw(picture)
        draw.rectangle([0,0,800,76],fill=(250,250,250))
        draw.text((14,7),'SNS robot | 75 kg | 27 actuators per side',fill=(20,20,20),font=font)
        support=100*trace['rig_force'][index,1]/(75*9.81)
        draw.text((14,33),'t = %.2f s | actual vertical support = %.1f%% BW'%(trace['time'][index],support),fill=(20,20,20),font=small)
        assistance='vertical, lateral and angular pelvis assistance' if shifted else 'vertical and angular pelvis assistance'
        draw.text((14,54),'%.0f%% activation floor | '%floor+assistance,fill=(20,20,20),font=small)
        draw.rectangle([0,634,800,660],fill=(250,250,250))
        draw.text((14,638),'x anterior, y superior, z right | hybrid SNS; reference tension actuators',fill=(20,20,20),font=small)
        return picture
    if args.snapshot:
        index=int(np.argmin(abs(trace['time']-args.time)))
        frame(index).save(args.snapshot)
        print('Saved snapshot',args.snapshot)
    if args.gif:
        indexes=[int(np.argmin(abs(trace['time']-t))) for t in np.arange(args.start_time,trace['time'][-1],.25)]
        frames=[frame(i).convert('P',palette=Image.Palette.ADAPTIVE,colors=128) for i in indexes]
        frames[0].save(args.gif,save_all=True,append_images=frames[1:],duration=250,loop=0,optimize=True)
        print('Saved recorded-motion GIF',args.gif)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,default=HERE/'runs/stand_20s_exact.npz')
    parser.add_argument('--snapshot',type=Path)
    parser.add_argument('--time',type=float,default=18.)
    parser.add_argument('--start-time',type=float,default=0.,help='Start GIF after this recorded time; raw simulation remains unchanged')
    parser.add_argument('--gif',type=Path)
    parser.add_argument('--view',action='store_true')
    main(parser.parse_args())
