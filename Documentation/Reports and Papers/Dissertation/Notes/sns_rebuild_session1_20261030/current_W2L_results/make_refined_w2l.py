from pathlib import Path
import json,hashlib,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Polygon,FancyBboxPatch,FancyArrowPatch,Circle
R=Path(r'D:\GitHub\Bipedal_Robot');P=R/'Code/MuJoCo_SNS/spinal/virtual_walker_20261031/runs/GENTLE108_EFFECTIVE_CONNECTOME.json';D=json.loads(P.read_text());E=D['edges'];N={n['name']:n for n in D['neurons']};OUT=Path(__file__).parent/'refined';OUT.mkdir(exist_ok=True);audit=[];glyphs={}
plt.rcParams.update({'font.family':'Arial','font.size':10.5,'svg.fonttype':'none','pdf.fonttype':42})
fig=plt.figure(figsize=(9,5.5));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,9);ax.set_ylim(0,5.5);ax.axis('off');ax.set_aspect('equal')
def txt(x,y,s,**kw):ax.text(x,y,s,fontsize=10.5,**kw)
def panel(x,y,w,h,title):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.025,rounding_size=.07',facecolor='#f5f8fb',edgecolor='#9ba6b2',lw=.7,zorder=0));txt(x+.1,y+h-.17,title,fontweight='bold',ha='left',va='center')
def node(key,name,x,y,label=None,off=False,w=.53,h=.33):
 if off:patch=Polygon([(x-w/2,y+h/2),(x+w/2,y+h/2),(x+w/2,y-h*.20),(x,y-h/2),(x-w/2,y-h*.20)],facecolor='#ffdf66',edgecolor='#333333',lw=.7,zorder=4)
 else:patch=Ellipse((x,y),w,h,facecolor='#9dd3ef',edgecolor='#333333',lw=.7,zorder=4)
 ax.add_patch(patch);txt(x,y,label or name,ha='center',va='center',zorder=5,linespacing=1.0);glyphs[key]={'name':name,'xy':np.array([x,y]),'w':w,'h':h,'off':off}
def endpoint(g,other):
 v=np.asarray(other)-g['xy'];k=min([g['w']/2/abs(v[0])if abs(v[0])>1e-9 else np.inf,g['h']/2/abs(v[1])if abs(v[1])>1e-9 else np.inf]);return g['xy']+v*k
def edge(a,b,rad=0,via=None):
 A=glyphs[a];B=glyphs[b];found=[e for e in E if e['source']==A['name']and e['target']==B['name']];assert len(found)==1,(A['name'],B['name'],len(found));e=found[0];audit.append(dict(e,figure_source=a,figure_target=b))
 # Reversal/type defines the explicit graded synapse sign convention.
 inhibitory=('Inhibit'in e['reference_synapse_type_name']or e['reversal_mV']<=-60)
 if via:
  pts=[endpoint(A,via[0])]+[np.array(v)for v in via]+[endpoint(B,via[-1])];ax.plot(*np.array(pts).T,color='#444444',lw=.8,zorder=2);start,end=pts[-2],pts[-1]
 else:start=endpoint(A,B['xy']);end=endpoint(B,A['xy'])
 if inhibitory:
  ax.add_patch(FancyArrowPatch(start,end,connectionstyle=f'arc3,rad={rad}',arrowstyle='-',lw=.85,color='#444444',zorder=2));ax.add_patch(Circle(end,.025,facecolor='#222222',edgecolor='#222222',zorder=3))
 else:ax.add_patch(FancyArrowPatch(start,end,connectionstyle=f'arc3,rad={rad}',arrowstyle='->',mutation_scale=9,lw=.85,color='#444444',zorder=2))
 return e
panel(.1,2.85,4.3,2.5,'(A) Rhythm generation and commissural input')
node('rge','L RG ext',1.65,3.35,'RG-E');node('rgf','L RG flx',2.85,3.35,'RG-F');node('rgiE','L RG ext IN',1.65,4.45,'I-E');node('rgiF','L RG flx IN',2.85,4.45,'I-F')
for a,b,rad in [('rge','rgiE',0),('rgf','rgiF',0),('rgiE','rgf',0),('rgiF','rge',0),('rge','rgf',.45),('rgf','rge',.45)]:edge(a,b,rad)
node('contraC','R c1',.62,3.35,'R c1',True);node('contraV','R V3',3.77,4.45,'R V3',True);edge('contraC','rgf',via=[(.62,3.0),(2.85,3.0)]);edge('contraV','rgiE',via=[(3.77,4.91),(1.65,4.91)])
txt(.3,3.72,'Left RG shown; right circuit retained.',ha='left',va='center')
panel(4.6,2.85,4.3,2.5,'(B) Hip PF, motor neurons, and Renshaw cells')
for key,name,x,y,label in [('pfiE','L Hip PF ext IN',5.9,4.67,'I-E'),('pfiF','L Hip PF flx IN',7.6,4.67,'I-F'),('pfe','L Hip PF ext',5.9,4.08,'PF-E'),('pff','L Hip PF flx',7.6,4.08,'PF-F'),('mne','L Hip MN ext',5.9,3.52,'MN-E'),('mnf','L Hip MN flx',7.6,3.52,'MN-F'),('ree','L Hip MN ext RE',5.9,3.03,'RE-E'),('ref','L Hip MN flx RE',7.6,3.03,'RE-F')]:node(key,name,x,y,label)
for a,b,rad in [('pfe','pfiE',0),('pff','pfiF',0),('pfiE','pff',0),('pfiF','pfe',0),('pfe','mne',0),('pff','mnf',0),('mne','ree',0),('mnf','ref',0),('ree','mne',.6),('ref','mnf',-.6),('ree','mnf',.1),('ref','mne',-.1)]:
 if any(e['source']==glyphs[a]['name']and e['target']==glyphs[b]['name']for e in E):edge(a,b,rad)
node('rgEref','L RG ext',5.04,4.08,'RG-E',True);node('rgFref','L RG flx',8.46,4.08,'RG-F',True);edge('rgEref','pfe');edge('rgFref','pff')
panel(.1,.15,4.3,2.47,'(C) Homonymous stretch-II feedback: all 12 channels')
txt(1.9,2.2,'Left',ha='center');txt(3.4,2.2,'Right',ha='center')
for j,joint in enumerate(('Hip','Knee','Ank')):
 y0=1.82-j*.62;txt(.24,y0-.14,'Ankle'if joint=='Ank'else joint,ha='left',va='center')
 for side,xx in(('L',1.5),('R',3.0)):
  for token,dy in(('ext',0),('flx',-.29)):
   target=f'{side} {joint} MN {token}';es=[e for e in E if e['target']==target and' II relay'in e['source']];assert len(es)==1,es;source=es[0]['source'];suffix='E'if token=='ext'else'F';a=f'II_{side}_{joint}_{token}';b=a+'_MN';node(a,source,xx,y0+dy,'II-'+suffix,w=.46,h=.23);node(b,target,xx+.65,y0+dy,'MN-'+suffix,True,w=.50,h=.25);edge(a,b)
panel(4.6,.15,4.3,2.47,'(D) Contact feedback: left pathway and opposite RG')
node('heel','L heel contact',5.07,1.78,'Heel L',w=.58,h=.34);node('toe','L toe contact',5.07,.75,'Toe L',w=.58,h=.34);node('opRG','R RG flx',6.36,2.12,'R RG-F',True,w=.73);edge('heel','opRG',via=[(5.07,2.12)])
targets=[('L RG ext','RG-E',6.0,1.44),('L Hip PF ext','Hip\nPF-E',7.04,1.44),('L Knee PF ext','K/A\nPF-E',8.11,1.44),('L Hip MN ext','Hip\nMN-E',6.0,.63),('L Knee MN ext','Knee\nMN-E',7.04,.63),('L Ank MN ext','Ankle\nMN-E',8.11,.63)]
for i,(name,label,x,y)in enumerate(targets):
 k='contactTarget'+str(i);node(k,name,x,y,label,True,w=.68,h=.48);edge('heel',k,via=[(5.53,1.78),(x,1.78)])if i<3 else edge('heel',k,via=[(5.53,1.78),(5.53,1.01),(x,1.01)])
edge('toe','contactTarget5',via=[(5.07,.34),(8.11,.34)])
txt(4.78,.19,'Right-side counterparts retained; mapping lists exact IDs.',ha='left',va='bottom')
for ext in('pdf','svg','png'):fig.savefig(OUT/('W2L108_MAIN_SUBNETWORKS_REFINED.'+ext),dpi=210)
plt.close(fig);record={'source_connectome':str(P),'source_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'canvas_inches':[9,5.5],'minimum_font_pt':10.5,'displayed_edges':audit,'displayed_edge_count':len(audit),'unique_connection_IDs':len({e['source_connection_ID']for e in audit}),'scope':'Reflowed selected subnetworks. All displayed biological edges match executed endpoints and signed graded synapses. Right-side RG, hip PF/MN, and contact counterparts, Ia/Ib reflexes, and4retainedphaseII pathways are outside the expanded views; no complete108-node claim.','model_files_changed':False};(OUT/'EDGE_AUDIT.json').write_text(json.dumps(record,indent=2));print(OUT)
