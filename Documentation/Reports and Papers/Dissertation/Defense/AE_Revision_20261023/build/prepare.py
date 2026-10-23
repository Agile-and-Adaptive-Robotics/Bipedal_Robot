from pathlib import Path
import shutil, json, hashlib, zipfile
import fitz
from PIL import Image,ImageDraw
root=Path(__file__).resolve().parents[1]
defense=root.parent
rev=defense.parent/'Notes/advisor_revision_20261022/Revision_C'
fig=defense/'deck_builder/deck_assets/figs'
paper=defense.parent.parent/'Knee_Torque_Test/Figures'
sources={}
def asset(key,path):
    path=Path(path)
    dest=root/'assets'/f'{key}.png'
    if path.suffix.lower()=='.pdf':
        doc=fitz.open(path); doc[0].get_pixmap(matrix=fitz.Matrix(2.5,2.5),alpha=False).save(dest)
    else:
        im=Image.open(path); im.convert('RGB').save(dest)
    sources[key]={'source':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'asset':str(dest)}
for key,name in {'jig':'testJigs_force.png','max10':'MaxForce10.png','max20':'MaxForce20.png','surface10':'FStar10.png','surface20':'FStar20.png','comparison':'Comparison.png','knee_jig':'testJigs_knee.png','hybrid':'b_flxpin_preopt.png','corrected':'b_flxpin_opt.png','bracket':'b_bktframe_pin.png','oldflex':'b_flxgrp_trim.png','oldext':'b_ext10mm_trim.png','robot':'c_mujoco_converted_robot.png','initialbody':'c_animatlab_body_gui.png','initialtrace':'c_animatlab_phase1_joint_angles.png','bilateral':'c_bilateralrg_sns_diagram.png','traces':'c_sns_w2l_traces.png','library':'c_simulink_sns_library_grid.png','spiking':'c_simulink_knee_spiking_twin.png','supported':'c_bpa_nonspiking_joints.png'}.items(): asset(key,fig/name)
for key,name in {'hardware':'cam_extensor_hardware.jpg','camroute':'extensor_cam_route.pdf','camresult':'ExtensorRoute_redesigned.pdf','extphysical':'extensor_20mm_physical.png','flexphysical':'flexor_20mm_physical.png','pulley':'pulley_margin_resolved.pdf'}.items(): asset(key,rev/'figs/30-torque'/name)
for key,name in {'clusters':'sadb_clusters_app.png','pivot':'sadb_pivot_app.png'}.items(): asset(key,rev/'figs/50-database'/name)
for key,name in {'reflex':'KneeReflex_circuit.pdf','pipeline':'pipeline.pdf'}.items(): asset(key,rev/'figs/60-sns'/name)
for key,name in {'human':'Flx20mm_humanTarget.pdf','rigid':'Flx20mm_rigidVsCorrected.pdf','meets':'Flx20mm_meetsHuman.pdf'}.items(): asset(key,paper/name)
for key,name in {'forceeq':'force.png','fstareq':'fstar.png'}.items(): asset(key,defense/'deck_builder/deck_assets/eqs'/name)
(root/'build/asset_manifest.json').write_text(json.dumps(sources,indent=2),encoding='utf-8')
for start in range(0,len(sources),12):
    canvas=Image.new('RGB',(1600,1020),'#ddd'); d=ImageDraw.Draw(canvas)
    for j,(key,info) in enumerate(list(sources.items())[start:start+12]):
        im=Image.open(info['asset']); im.thumbnail((395,305)); x=j%4*400;y=j//4*340
        canvas.paste(im,(x,y+25));d.text((x+5,y+5),key,fill='black')
    canvas.save(root/f'build/assets_{start//12+1}.jpg')
print('Prepared',len(sources),'source-traceable assets')
