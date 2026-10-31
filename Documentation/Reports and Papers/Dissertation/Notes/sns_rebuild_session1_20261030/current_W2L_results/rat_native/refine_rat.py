from pathlib import Path
from lxml import etree as E
import json,csv,hashlib,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Circle,Polygon,FancyArrowPatch
R=Path(r'D:\GitHub\Bipedal_Robot');b=R/'Documentation/Reports and Papers/Dissertation/Notes/sns_rebuild_session1_20261030/current_W2L_results';out=b/'rat_native/refined';out.mkdir(exist_ok=True)
P=R/'Neuromechanical_Models/Deng_Biarticular_Rat/Biarticular_net.aproj';T=E.parse(str(P));alln=T.findall('.//NervousSystem//Node');pages={n.findtext('Text').strip():n for n in alln if n.find('DiagramXml')is not None};types={x.findtext('ID'):x for x in T.findall('.//SynapseTypes/Link')};audit=[]
plt.rcParams.update({'font.family':'Arial','font.size':10.5,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(9,5.5));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,9);ax.set_ylim(0,5.5);ax.axis('off')
def frame(x,y,w,h,title):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.02,rounding_size=.10',facecolor='#fafcfd',edgecolor='#697984',lw=.8,zorder=0));ax.text(x+.10,y+h-.10,title,weight='bold',ha='left',va='top',fontsize=11)
def node(x,y,label,port=False):
 if port:
  w,h=.62,.40
  ax.add_patch(Polygon([(x-w/2,y+h/2),(x+w/2,y+h/2),(x+w/2,y-h*.28),(x,y-h/2),(x-w/2,y-h*.28)],closed=True,facecolor='#fff1a8',edgecolor='#233746',lw=.9,zorder=5))
 else:ax.add_patch(Circle((x,y),.205,facecolor='#c9e4f3',edgecolor='#233746',lw=.95,zorder=5))
 ax.text(x,y,label,ha='center',va='center',fontsize=10.5,zorder=6,linespacing=.95)
def endpoint(pos,target,port):
 x,y=pos;v=np.array(target)-np.array(pos);norm=np.linalg.norm(v)
 if not port:return np.array(pos)+v/norm*.208
 scale=min(.32/abs(v[0]) if abs(v[0])>1e-8 else 1e8,.205/abs(v[1]) if abs(v[1])>1e-8 else 1e8)
 return np.array(pos)+v*scale
def line(a,b,inh=False,pts=None,rad=0):
 if pts is not None:
  ps=[a]+pts+[b]
  for v,w in zip(ps[:-2],ps[1:-1]):ax.plot([v[0],w[0]],[v[1],w[1]],color='#233746',lw=.8,zorder=2)
  start=ps[-2]
 else:start=a
 if inh:
  if rad:
   ar=FancyArrowPatch(start,b,connectionstyle=f'arc3,rad={rad}',arrowstyle='-',lw=.8,color='#233746',zorder=2);ax.add_patch(ar)
  else:ax.plot([start[0],b[0]],[start[1],b[1]],color='#233746',lw=.8,zorder=2)
  ax.add_patch(Circle(b,.035,facecolor='#233746',edgecolor='#233746',zorder=4))
 else:ax.add_patch(FancyArrowPatch(start,b,connectionstyle=f'arc3,rad={rad}',arrowstyle='-|>',mutation_scale=9,lw=.8,color='#233746',zorder=2))
def page_data(name):
 pg=pages[name];ns={n.findtext('ID'):n for n in pg.findall('Nodes/Node')};return ns,pg.findall('Links/Link')
def draw_page(name,pos,labels,curves={}):
 ns,links=page_data(name)
 for e in links:
  si,ti=e.findtext('OriginID'),e.findtext('DestinationID')
  if si not in pos or ti not in pos:continue
  st=types[e.findtext('SynapticTypeID')];rev=float(st.find('EquilibriumPotential').get('Actual'))*1000;inh=rev<=-60
  n1,n2=ns[si],ns[ti];p1,p2=pos[si],pos[ti];po1='OffPage'in n1.findtext('ClassName','');po2='OffPage'in n2.findtext('ClassName','')
  source_label=labels[si].replace('\n',' ');target_label=labels[ti].replace('\n',' ')
  a=endpoint(p1,p2,po1);z=endpoint(p2,p1,po2)
  line(a,z,inh,rad=curves.get((source_label,target_label),0))
  def val(k):v=st.find(k);return None if v is None else float(v.get('Actual'))
  g=e.find('SynapticConductance')
  audit.append({'page':name,'connection_ID':e.findtext('ID'),'source_ID':si,'target_ID':ti,'source':n1.findtext('Text'),'target':n2.findtext('Text'),'type_ID':st.findtext('ID'),'type_name':st.findtext('Name'),'equilibrium_mV':rev,'conductance_uS':float(g.get('Actual'))*1e6 if g is not None else val('MaxSynapticConductance')*1e6,'threshold_mV':val('PreSynapticThreshold')*1000,'saturation_mV':val('PreSynapticSaturationLevel')*1000,'display_sign':'inhibitory'if inh else'excitatory'})
 for i,p in pos.items():node(*p,labels[i],'OffPage'in ns[i].findtext('ClassName',''))
def by_label(name):return {n.findtext('Text').strip():i for i,n in page_data(name)[0].items()}
frame(.10,2.87,4.28,2.49,'(A) Rhythm-generation subnetwork')
frame(4.60,2.87,4.28,2.49,'(B) Hip pattern-formation subnetwork')
frame(.10,.18,4.28,2.44,'(C) Biarticular integration subnetwork')
frame(4.60,.18,4.28,2.44,'(D) Imported model and notation')
# Native left RG: four intrinsic neurons and every connected sensory port.
m=by_label('LH_RG');pos={};labels={}
for text,p,lab in [('LH_RG IN EXT',(2.15,4.63),'I-E'),('LH_RG IN FLX',(3.43,4.63),'I-F'),('LH_RG EXT',(2.15,3.65),'RG-E'),('LH_RG FLX',(3.43,3.65),'RG-F'),('LH_IP-    FLX Ia',(.62,4.64),'Ia-F'),('LH_HipZ-FLX II',(.62,3.73),'II-F'),('LH_SO- Ext Ib',(1.45,3.19),'Ib-E'),('LH_HipZ-Ext II',(3.93,3.19),'II-E')]:pos[m[text]]=p;labels[m[text]]=lab
curves={('RG-E','RG-F'):.18,('RG-F','RG-E'):.18,('I-E','RG-F'):-.08,('I-F','RG-E'):-.08,('Ia-F','RG-E'):.10,('II-F','RG-E'):-.08}
draw_page('LH_RG',pos,labels,curves)
# Hip PF: omit the three native off-page objects that have no connections on this page.
m=by_label('LH_Hip');pos={};labels={}
for text,p,lab in [('LH_Hip  PF IN EXT',(6.69,4.63),'I-E'),('LH_Hip  PF IN FLX',(8.01,4.63),'I-F'),('LH_Hip  PF EXT',(6.69,3.58),'PF-E'),('LH_Hip  PF FLX',(8.01,3.58),'PF-F'),('LH_IP-    FLX Ia',(5.15,4.63),'Ia-F'),('LH_HipZ-FLX II',(5.15,3.64),'II-F'),('LH_HipZ-Ext II',(8.12,3.12),'II-E')]:pos[m[text]]=p;labels[m[text]]=lab
draw_page('LH_Hip',pos,labels,{('I-E','PF-F'):-.06,('I-F','PF-E'):-.06})
# Biar source contains two integration interneurons and four PF inputs. Do not invent an antagonistic four-cell circuit.
m=by_label('Biar');pos={};labels={}
for text,p,lab in [('LH_Hip  PF EXT',(.68,1.79),'Hip\nPF-E'),('LH_Knee PF FLX',(1.73,1.79),'Knee\nPF-F'),('LH_Hip  PF FLX',(2.77,1.79),'Hip\nPF-F'),('LH_Knee PF EXT',(3.82,1.79),'Knee\nPF-E'),('LH_BFP IN EXT',(1.20,.77),'BFP\nI-E'),('LF_RF IN EXT',(3.30,.77),'RF\nI-E')]:pos[m[text]]=p;labels[m[text]]=lab
draw_page('Biar',pos,labels)
# Model image is context; explicit legend distinguishes native ports from neuron pools.
img=plt.imread(b/'rat_native/RAT_SOURCE_STATIC_BODY.png');ax.imshow(img,extent=(4.78,6.66,.47,2.28),zorder=1)
ax.text(6.84,2.16,'Bilateral hindlimbs\n120 neurons\n16 muscles\nFixed pelvis',ha='left',va='top',fontsize=10.5,linespacing=1.3)
ax.add_patch(Circle((6.94,.92),.095,facecolor='#c9e4f3',edgecolor='#233746',lw=.9));ax.text(7.13,.92,'Neuron',va='center',fontsize=10.5)
ax.add_patch(Polygon([(6.84,.64),(7.04,.64),(7.04,.48),(6.94,.43),(6.84,.48)],facecolor='#fff1a8',edgecolor='#233746',lw=.9));ax.text(7.13,.54,'Off-page input',va='center',fontsize=10.5)
ax.text(.20,.035,'Arrowhead: excitation     Filled circle: inhibition     E/F: extensor/flexor',fontsize=10.5,va='bottom')
for ext in ['pdf','svg','png']:fig.savefig(out/('DENG_RAT_NATIVE_SUBNETWORKS_REFINED.'+ext),dpi=190)
plt.close(fig)
with(out/'DISPLAYED_EDGES.csv').open('w',encoding='utf-8',newline='')as f:w=csv.DictWriter(f,fieldnames=list(audit[0]));w.writeheader();w.writerows(audit)
(out/'AUDIT.json').write_text(json.dumps({'source':str(P),'source_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'page_names':['LH_RG','LH_Hip','Biar'],'displayed_edges':len(audit),'minimum_print_font_pt':10.5,'canvas_inches':[9,5.5],'omitted_unconnected_page_objects':['LH_Hip: RG EXT','LH_Hip: RG FLX','LH_Hip: SO Ext Ib'],'layout':'Publication reflow; native source endpoints and parameters retained, positions changed for readability.','source_or_results_modified':False},indent=2),encoding='utf-8')
(out/'CAPTION.tex').write_text(r'''Imported public biarticular hindlimb model: \textbf{(A)} left rhythm-generation subnetwork, \textbf{(B)} hip pattern-formation subnetwork, \textbf{(C)} native BFP/RF biarticular integration, and \textbf{(D)} model context and symbol key. Blue circles are neurons, yellow terminals identify inputs from other native pages, arrowheads denote excitation, and filled circles denote inhibition. The publication layout separates the subnetworks while preserving their source connections. BFP and RF denote posterior biceps femoris and rectus femoris; E/F denote extensor/flexor. The public precursor retains a fixed pelvis.''',encoding='utf-8')
print(out);print('Audited displayed edges:',len(audit))
