"""Render review sheets from actual saved project page objects and drawing geometry."""
from pathlib import Path
from lxml import etree as E
import json,numpy as np,os,textwrap
os.environ['OPENBLAS_NUM_THREADS']='1'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch,Ellipse,Rectangle,FancyBboxPatch
W=Path(__file__).resolve().parents[1];R=W.parents[1];B=W/'build/subnetwork_rebuild_20261031'
plt.rcParams.update({'font.family':'Arial','svg.fonttype':'none'})
def color(s,default='#eeeeee'):
 try:return '#'+format(int(s)&0xffffff,'06x')
 except:return default
def main(folder,name):
 d=R/'Neuromechanical_Models'/folder;O=d/'subnetwork_review';O.mkdir(exist_ok=True)
 root=E.parse(str(d/(name+'.aproj'))).getroot();pages=[n for n in root.findall('.//NervousSystem//Node') if n.findtext('ClassName')=='AnimatGUI.DataObjects.Behavior.Nodes.Subsystem']
 summary=[]
 def render(p,ax):
  ns={n.findtext('ID'):n for n in p.findall('Nodes/Node')};ls=p.findall('Links/Link');flow=E.fromstring(p.findtext('DiagramXml').encode()).find('.//AddFlow');positions={z.findtext('Tag'):(float(z.get('Left')),float(z.get('Top')),float(z.get('Width')),float(z.get('Height'))) for z in flow if z.tag=='Node'}
  for l in ls:
   a,b=l.findtext('OriginID'),l.findtext('DestinationID');x,y,w,h=positions[a];u,v,m,n=positions[b];style=l.findtext('ArrowDestination/Style');end='-|>' if style not in ('Circle','CircleFilled') else '-['
   ax.add_patch(FancyArrowPatch((x+w/2,y+h/2),(u+m/2,v+n/2),arrowstyle=end,color=color(l.findtext('DrawColor'),'#777777'),lw=.7,mutation_scale=8,shrinkA=18,shrinkB=18,alpha=.7,zorder=1))
  for i,n in ns.items():
   x,y,w,h=positions[i];klass=n.findtext('ClassName') or '';alias='OffPage' in klass;sub='Subsystem' in klass
   shape=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=1,rounding_size=8',facecolor='#dceaf5' if sub else '#fff0ce' if alias else color(n.findtext('FillColor')),edgecolor='#333',linestyle='--' if alias else '-',lw=.8,zorder=2)
   ax.add_patch(shape);label=n.findtext('Text') or '(unnamed)';ax.text(x+w/2,y+h/2,'\n'.join(textwrap.wrap(label,21)),ha='center',va='center',fontsize=8,zorder=3)
  if positions:
   ax.set_xlim(min(t[0] for t in positions.values())-35,max(t[0]+t[2] for t in positions.values())+35);ax.set_ylim(max(t[1]+t[3] for t in positions.values())+40,min(t[1] for t in positions.values())-40)
  ax.set_aspect('equal');ax.axis('off');ax.set_title(p.findtext('Text')+' — '+str(len(ns))+' objects, '+str(len(ls))+' links',loc='left',fontsize=11)
 for i,p in enumerate(pages):
  fig,ax=plt.subplots(figsize=(12,10),layout='constrained');render(p,ax);fig.suptitle(name+' | saved AnimatLab page',fontsize=14);fig.text(.02,.01,'Amber dashed boxes are OffPage references, not additional neurons. Original numerical circuit retained.',fontsize=9)
  slug='page_'+format(i,'02d')+'_'+''.join(c if c.isalnum() else '_' for c in p.findtext('Text'));fig.savefig(O/(slug+'.svg'));plt.close(fig);summary.append({'page':p.findtext('Text'),'svg':str(O/(slug+'.svg'))})
 chosen=[p for p in pages if p.findtext('Text')=='Neural Subsystem' or p.findtext('Text').startswith('LH')]
 fig,axes=plt.subplots(4,2,figsize=(20,28),layout='constrained')
 for ax,p in zip(axes.flat,chosen):render(p,ax)
 for ax in list(axes.flat)[len(chosen):]:ax.axis('off')
 fig.suptitle(name+' — actual saved LH subnetwork pages',fontsize=22);fig.savefig(O/'LH_contact_sheet.png',dpi=150);plt.close(fig)
 before=np.load(B/folder/'verification_before.npz',allow_pickle=True)['data'];after=np.load(d/'verification.npz',allow_pickle=True)['data'];results=json.loads((d/'TEST_RESULTS.json').read_text());labels=json.loads((d/'diagnostic_columns.json').read_text());delta=np.abs(before-after)
 record={'model':name,'native_exit':results['exit_code'],'samples':results['samples'],'finite':bool(np.isfinite(after).all()),'native_standalone_identical_to_before':bool(np.array_equal(before,after)),'native_max_abs_difference':float(delta.max()),'root_y':results['signals']['root_y'],'net_root_x_m':results['signals']['root_x']['last']-results['signals']['root_x']['first'],'GUI_verification':'blocked: CUA GetCursorPos Access denied; existing dirty user AnimatLab window left untouched','root_y_scope':'world Y telemetry already supplied; preserved. Stable root Y with passive rail support is not free standing/walking or measured support force.','review_pages':summary}
 (d/'SUBNETWORK_VERIFICATION.json').write_text(json.dumps(record,indent=2));print(json.dumps({k:v for k,v in record.items() if k!='review_pages'},indent=2))
for args in [('W2L_VirtualWalker_20261031','W2L_VirtualWalker'),('W2L_VirtualWalker_Newer_20261031','W2L_VirtualWalker_Newer')]:main(*args)
