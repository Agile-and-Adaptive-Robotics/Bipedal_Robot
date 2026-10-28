"""Render recorded, verified robot physics for a defense GIF and print still."""
import argparse,json,hashlib,xml.etree.ElementTree as ET
from pathlib import Path
import os,sys
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import mujoco
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from dynamic_routes import DynamicRoutes
from robot_kinematics import camera_axes
HERE=Path(__file__).resolve().parent

def main(args):
    summary=json.loads(args.run.with_suffix('.json').read_text());trace=np.load(args.run)
    source=Path(summary.get('model_path',HERE/'ah001_neural.xml'))
    tree=ET.parse(source);root=tree.getroot()
    root.find('compiler').set('meshdir',str(HERE/'meshes'))
    global_visual=root.find('visual/global');global_visual.set('offwidth','1600');global_visual.set('offheight','1600')
    ET.SubElement(root.find('asset'),'texture',name='defense_sky',type='skybox',builtin='gradient',
        rgb1='1 1 1',rgb2='.92 .94 .97',width='512',height='3072')
    root.find("worldbody/geom[@name='ground']").set('rgba','.94 .95 .96 1')
    cameras=[]
    for name,position in [('Anterior',np.array([2.4,.95,0.])),('Oblique',np.array([1.8,1.05,1.8]))]:
        target=np.array([-.04,.66,0.]);axes=camera_axes(position,target)
        ET.SubElement(root.find('worldbody'),'camera',name='defense_'+name,pos=' '.join(map(str,position)),
            xyaxes=' '.join(map(str,np.r_[axes[:,0],axes[:,1]])),fovy='37')
        cameras.append('defense_'+name)
    model=mujoco.MjModel.from_xml_string(ET.tostring(root,encoding='unicode'))
    physical=mujoco.MjModel.from_xml_path(str(source))
    assert np.array_equal(model.body_mass,physical.body_mass) and np.array_equal(model.qpos0,physical.qpos0)
    data=mujoco.MjData(model);dynamic=DynamicRoutes(model,json.loads((HERE/'robotbody_routes.json').read_text()))
    for gid in range(model.ngeom):
        name=mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_BODY,int(model.geom_bodyid[gid])) or ''
        if name.endswith('_r'):model.geom_rgba[gid]=[.267,.467,.667,1.]
        elif name.endswith('_l'):model.geom_rgba[gid]=[.80,.56,.25,1.]
    font_path='C:/Windows/Fonts/arial.ttf'
    args.output.mkdir(parents=True,exist_ok=True)
    renderers={}
    def render(index,scale=1):
        data.qpos[:]=trace['qpos'][index];data.act[:]=trace['activation'][index];dynamic.update(data)
        for aid,tid in enumerate(model.actuator_trnid[:,0]):
            name=mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,aid)
            emphasized=name.startswith(('tfl_','sar_','glut_max3_'))
            a=float(trace['activation'][index,aid])
            model.tendon_rgba[tid]=([.13,.53,.27,.35+.6*a] if emphasized else [.60,.30,.28,.16+.45*a])
            model.tendon_width[tid]=.0018 if emphasized else .0011
        width,height=600*scale,700*scale
        if scale not in renderers:renderers[scale]=mujoco.Renderer(model,height=height,width=width)
        renderer=renderers[scale]
        views=[]
        for camera in cameras:
            renderer.update_scene(data,mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_CAMERA,camera))
            views.append(Image.fromarray(renderer.render()))
        canvas=Image.new('RGB',(1200*scale,800*scale),'white')
        for i,view in enumerate(views):canvas.paste(view,(600*scale*i,50*scale))
        draw=ImageDraw.Draw(canvas)
        title=ImageFont.truetype(font_path,27*scale);label=ImageFont.truetype(font_path,25*scale)
        draw.text((20*scale,12*scale),args.title,fill='#333333',font=title)
        for i,name in enumerate(('A  Anterior view','B  Oblique view')):
            draw.text(((20+600*i)*scale,60*scale),name,fill='#333333',font=label)
        support=100*trace['rig_force'][index,1]/(float(model.body_mass.sum())*9.81)
        footer='t = %.2f s | vertical support %.0f%% BW | minimum activation %.0f%%'%(trace['time'][index],support,100*summary['run_options']['activation_floor'])
        draw.text((20*scale,748*scale),footer,fill='#333333',font=label)
        draw.text((20*scale,777*scale),args.scope_label,fill='#555555',font=ImageFont.truetype(font_path,21*scale))
        return canvas
    begin=max(args.start,float(trace['time'][0]));end=min(args.end,float(trace['time'][-1]))
    indexes=[int(np.argmin(abs(trace['time']-t))) for t in np.arange(begin,end,1/args.fps)]
    # Keep one OpenGL context alive throughout the animation render.
    frames=[]
    for number,index in enumerate(indexes):
        frames.append(render(index).convert('P',palette=Image.Palette.ADAPTIVE,colors=192))
        if number%32==0:print('Rendered',number,'/',len(indexes),flush=True)
    gif=args.output/'ah001_abductor_control.gif'
    ticks=np.rint(np.arange(len(frames)+1)*100/args.fps).astype(int)
    durations=(10*np.diff(ticks)).tolist()
    frames[0].save(gif,save_all=True,append_images=frames[1:],duration=durations,loop=0,optimize=True)
    still_index=int(np.argmin(abs(trace['time']-args.still_time)))
    still=args.output/'ah001_abductor_control_still.png'
    render(still_index,2).crop((0,100,2400,1500)).save(still,dpi=(300,300))
    renderers.clear()
    metadata=dict(source_run=str(args.run),source_model=str(source),source_summary=summary,
        gif_start_s=begin,gif_end_s=end,gif_frames=len(frames),gif_fps=args.fps,
        still_time_s=float(trace['time'][still_index]),still_pixels=[2400,1400],
        visual_changes_only='Light background, camera positions, right/left colors, and selected route highlighting; recorded qpos and activations unchanged',
        outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (gif,still)})
    (args.output/'asset_provenance.json').write_text(json.dumps(metadata,indent=2));print('Saved',gif,still,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--start',type=float,default=20.);p.add_argument('--end',type=float,default=32.)
    p.add_argument('--still-time',type=float,default=24.);p.add_argument('--fps',type=int,default=8)
    p.add_argument('--title',default='Supported CAD robot with increased abductor control')
    p.add_argument('--scope-label',default='Hybrid control; reference tensions; lateral/angular assistance and knee-stop surrogates')
    main(p.parse_args())
