"""Readable 14 cm workflow schematic; redraws the prior vector source layout.

No experimental data or model parameters are involved. All text is >=10.5 pt
at the declared 14 cm manuscript width. Output remains in Notes for review.
"""
from pathlib import Path
import hashlib, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import fitz

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[4]
SOURCE=OUT.parent/'ch6_revision_20261116/figscripts/prepare_chapter_assets.py'
plt.rcParams.update({'font.family':'Arial','font.size':10.5,'pdf.fonttype':42,
                     'svg.fonttype':'none','path.simplify':False})
fig,ax=plt.subplots(figsize=(14/2.54,14/2.54))
fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.99)
ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
ink,blue,green,grey='#202A31','#4477AA','#228833','#667078'

def box(x,y,w,h,text,fill='#EAF1F7',edge=blue):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.006,rounding_size=0.012',
                               facecolor=fill,edgecolor=edge,linewidth=1.1))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',color=ink,fontsize=10.5,linespacing=1.22)

def arrow(points,color=ink):
    for a,b in zip(points[:-2],points[1:-1]):
        ax.plot([a[0],b[0]],[a[1],b[1]],color=color,lw=1.2)
    ax.add_patch(FancyArrowPatch(points[-2],points[-1],arrowstyle='-|>',
                                mutation_scale=10,color=color,lw=1.2))

def heading(y,text):
    ax.text(.015,y,text,fontsize=10.5,fontweight='bold',va='top',color=ink)

heading(.99,'PREPARE AND CHECK THE MECHANICS')
box(.02,.815,.445,.115,'OpenSim\nConversion and repair\nMuJoCo model')
box(.535,.815,.445,.115,'SOLIDWORKS CAD\nImport or separate assembly\nSimscape / MuJoCo model')
arrow([(.245,.809),(.245,.755),(.843,.755),(.843,.651)],blue)
arrow([(.755,.809),(.755,.755)],blue)
heading(.717,'DEFINE AND EXECUTE THE COUPLED LOOP')
box(.02,.49,.275,.15,'SNS controller\nGraphical blocks\nor scripted network')
box(.362,.49,.275,.15,'Actuator interface\nActivation / pressure\nForce calculation')
box(.705,.49,.275,.15,'Mechanical model\nSimscape or\nMuJoCo')
arrow([(.301,.55),(.356,.55)])
arrow([(.643,.55),(.699,.55)])
ax.text(.328,.657,'command',ha='center',va='bottom',fontsize=10.5,color=grey)
ax.text(.671,.657,'tension',ha='center',va='bottom',fontsize=10.5,color=grey)
box(.27,.30,.46,.11,'Sensory encoding\nLength, velocity, force, contact',fill='#EAF5EE',edge=green)
arrow([(.843,.484),(.843,.355),(.736,.355)],green)
arrow([(.264,.355),(.157,.355),(.157,.484)],green)
heading(.248,'CALIBRATE, TEST, AND INSPECT')
box(.02,.035,.96,.155,'Functional targets and parameter calibration\nPhysical curriculum, interventions, and fresh replays\nAligned neural and mechanical recordings',fill='#F4F4F0',edge=grey)

pdf=OUT/'BUILD_workflow_20261116.pdf'
png=pdf.with_suffix('.png')
fig.savefig(pdf);fig.savefig(png,dpi=220);plt.close(fig)
alt=('Two mechanical-preparation routes feed the same experimental loop: SNS controller, '
     'actuator interface, mechanical model, and sensory encoding returning to the controller. '
     'The OpenSim route uses conversion and repair; the CAD route uses import or separate '
     'assembly. Functional calibration, physical tests, interventions, fresh replays, and '
     'aligned recordings assess the loop. No automatic Simscape-to-MuJoCo translation is implied.')
pdf.with_suffix('.alt.txt').write_text(alt,encoding='utf-8')
doc=fitz.open(pdf);page=doc[0]
spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
result={'source_generator':str(SOURCE),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'generator':str(Path(__file__)),'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'output':str(pdf),'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),
        'size_cm':[14,14],'minimum_source_font_pt':min(s['size'] for s in spans),
        'minimum_font_at_14cm_pt':min(s['size'] for s in spans)*(14/2.54*72/page.rect.width),
        'experimental_data':None,'transforms':'Schematic layout condensed from prior vector source; no data transformation.',
        'text_and_topology':'Same two preparation routes, three forward elements, sensory return, and evaluation block.'}
doc.close()
(OUT/'BUILD_workflow_provenance.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
