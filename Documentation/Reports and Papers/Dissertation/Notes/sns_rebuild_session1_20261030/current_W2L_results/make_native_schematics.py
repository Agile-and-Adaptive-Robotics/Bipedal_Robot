"""Publication views of saved native AnimatLab drawings and an executed connectome.

Read-only with respect to the .aproj, .asim, and simulation artifacts.
"""
from pathlib import Path
from collections import Counter,defaultdict
import csv,json,hashlib,re,math,shutil,textwrap,sys
import numpy as np
from lxml import etree as E
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Rectangle,Polygon,FancyArrowPatch
from matplotlib.backends.backend_pdf import PdfPages
R=Path(r'D:\GitHub\Bipedal_Robot');W=R/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031';OUT=Path(__file__).parent/'native_schematic';OUT.mkdir(exist_ok=True)
aproj=W/'W2L_VirtualWalker_Newer.aproj';cp=R/'Code/MuJoCo_SNS/spinal/virtual_walker_20261031/runs/GENTLE108_EFFECTIVE_CONNECTOME.json';executed=json.loads(cp.read_text());assert len(executed['neurons'])==108 and len(executed['edges'])==216
tree=E.parse(str(aproj));nodes={n.findtext('ID'):n for n in tree.findall('.//NervousSystem//Node')};pages={n.findtext('Text').strip():n for n in nodes.values()if n.find('DiagramXml')is not None};edges={e['source_connection_ID']:e for e in executed['edges']};aliases={i:n.findtext('LinkedNodeID')for i,n in nodes.items()if'OffPage'in n.findtext('ClassName','')};types={s.findtext('ID'):s for s in tree.findall('.//SynapseTypes/Link')}
native=E.parse(str(W/'W2L_VirtualWalker_Newer_Standalone.asim'));native_connections={e.findtext('ID'):e for e in native.findall('.//Connexions/Connexion')}
for edge in edges.values():edge['reference_connection_G_uS']=float(native_connections[edge['source_connection_ID']].findtext('G'))
plt.rcParams.update({'font.family':'Arial','font.size':12.5,'axes.titlesize':13.5,'svg.fonttype':'none','pdf.fonttype':42})
audit={'source_project':str(aproj),'source_project_sha256':hashlib.sha256(aproj.read_bytes()).hexdigest(),'executed_connectome':str(cp),'executed_connectome_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),'neurons':108,'chemical_edges':216,'landscape_inches':[11,8.5],'minimum_font_pt':10,'publication_labels':{},'page_views':[],'missing_native_chemical_heads_completed':[],'target_overlays':[],'scope':'Native node positions, sizes, shapes, and explicit route control points are preserved within each panel. Physiological abbreviations replace long text only. Four executed target changes overlay the saved GUI routes; source files are not changed.'}
def resolve(i):
 seen=set()
 while i in aliases and i not in seen:seen.add(i);i=aliases[i]
 return i
def color(s,default='#ffffff'):
 try:return '#'+format(int(s)&0xffffff,'06x')
 except:return default
def label(n,keep_side=False):
 text=n.findtext('Text','').strip();s=re.sub(r'^[LR] ','',text);s=s.replace('extensor','ext').replace('flexor','flx');s=re.sub(r'\bext\b','E',s);s=re.sub(r'\bflx\b','F',s)
 cls=n.findtext('ClassName','')
 if cls.endswith('.RigidBody')and text=='6':return'Toe\nsensor'
 if cls.endswith('.RigidBody')and text=='HeelContact':return'Heel\nsensor'
 if 'PhysicalToNodeAdapter'in cls or'NodeToPhysicalAdapter'in cls:return'map'if'PhysicalToNode'in cls else'drive'
 if s in('L c1','R c1'):return s
 if s.startswith(('LH_','RH_')):return s.replace('LH_','L ').replace('RH_','R ').replace('K&A Pattern Formation','K/A PF').replace('K&A Synergy','K/A').replace('Contact feedback','contact').replace('HipZ','hip')
 if 'homonymous'in s.lower():s=s.replace('homonymous','')
 if ' II relay'in s:
  e=next((e for e in edges.values()if e['source_ID']==n.findtext('ID')and' MN 'in e['target']),None)
  if e:return 'II-'+('E'if e['target'].endswith('ext')else'F')
  m=re.search(r'\b(E|F)\b',s);return 'II-'+m.group(1)if m else'II'
 if 'Hip F II'in s:return'II-F\nphase'
 if 'Knee F II'in s:return'II-F\nphase'
 if 'Knee'in s and'PF'in s:s=s.replace('Knee','K/A')
 for joint in('Hip ','Knee ','Ank '):s=s.replace(joint,'')
 if 'StretchReceptor'in cls:return'M-'+('E'if re.search(r'\bE\b',s)else'F')
 if s in('E Ia','F Ia','E Ib','F Ib','E Ia IN','F Ia IN'):
  ss=s.split();s=ss[1]+'-'+ss[0]+ ('\nIN'if len(ss)>2 else'')
 if s in('RG E','RG F','PF E','PF F','MN E','MN F','Ia E','Ia F','Ib E','Ib F'):s=s.replace(' ','-')
 else:s=s.replace(' ext','-E').replace(' flx','-F').replace(' E','-E').replace(' F','-F')
 s=s.replace('MN-E RE','RE-E').replace('MN-F RE','RE-F').replace('RG-E IN','RG-E\nIN').replace('RG-F IN','RG-F\nIN').replace('PF-E IN','PF-E\nIN').replace('PF-F IN','PF-F\nIN').replace('Ia-E IN','Ia-E\nIN').replace('Ia-F IN','Ia-F\nIN')
 if s=='Commissural coordination':s='Commissural\ncoordination'
 if text in('LH','RH'):s=text
 if s=='heel contact':s='heel\ncontact'
 if s=='toe contact':s='toe\ncontact'
 if s=='K/A PF-E':s='K/A\nPF-E'
 if s=='K/A PF-F':s='K/A\nPF-F'
 if 'OffPage'in cls and' MN 'in text:
  group='Hip'if' Hip 'in text else'Knee'if' Knee 'in text else'Ankle'if' Ank 'in text else''
  if group:s=group+'\n'+s
 if 'OffPage'in cls and' Hip PF 'in text:s='Hip\n'+s
 if keep_side and re.match(r'^[LR] ',text):s=text[0]+' '+s
 return s
def render(name,ax,title,selected=None,execution=True,publish_rg_inputs=False,ii_targets=False):
 page=pages[name];flow=E.fromstring(page.findtext('DiagramXml').encode()).find('.//AddFlow');seq=[x for x in flow if x.tag in('Node','Link')];local={n.findtext('ID'):n for n in page.findall('Nodes/Node')};links={l.findtext('ID'):l for l in page.findall('Links/Link')};shapes={s.findtext('Tag'):s for s in seq if s.tag=='Node'}
 for s in shapes.values():
  x,y,w,h=[float(s.get(k))for k in('Left','Top','Width','Height')]
  cls=local[s.findtext('Tag')].findtext('ClassName','')
  scale=1.0 if 'Adapter' in cls or cls.endswith('.RigidBody') else 1.4
  s.attrib.update({'Left':str(x-(scale-1)*w/2),'Top':str(y-(scale-1)*h/2),'Width':str(w*scale),'Height':str(h*scale)})
 audit['native_node_centers_retained']=True;audit['publication_node_glyph_scale']=1.4;audit['minimum_font_pt']=12.5;audit['print_font_pt_at9inch_width']=12.5*9/11
 if selected is None:selected=set(shapes)
 def center(s):return np.array([float(s.get('Left'))+float(s.get('Width'))/2,float(s.get('Top'))+float(s.get('Height'))/2])
 def bound(s,other):
  c=center(s);d=np.asarray(other)-c;h=np.array([float(s.get('Width'))/2,float(s.get('Height'))/2]);k=min([h[i]/abs(d[i])for i in(0,1)if abs(d[i])>1e-10]or[0]);return c+k*d
 points=[];drawn=[]
 for s in seq:
  if s.tag!='Link':continue
  lid=s.findtext('Tag');l=links[lid];a=seq[int(s.get('Org'))];b=seq[int(s.get('Dst'))];ai=a.findtext('Tag');bi=b.findtext('Tag')
  if ai not in selected or bi not in selected:continue
  edge=edges.get(lid);changed=edge is not None and edge['target_changed']
  if changed and execution:
   dest=next((i for i in selected if resolve(i)==edge['target_ID']),None)
   assert dest,(name,lid,edge['target']);oldb=b;b=shapes[dest];bi=dest
  pts=[center(a)]+[np.array([float(q.get('X')),float(q.get('Y'))])for q in s.findall('Point')]+[center(b)];pts[0]=bound(a,pts[1]);pts[-1]=bound(b,pts[-2]);points.extend(pts)
  co='#7733AA'if changed and execution else color(s.findtext('DrawColor')or l.findtext('DrawColor'),'#444444');lw=1.15 if changed else .75
  ax.plot(*np.asarray(pts).T,color=co,lw=lw,zorder=1)
  arrow=s.find('ArrowDst');head=l.findtext('ArrowDestination/Style','')if arrow is None else arrow.get('Head','');filled= l.findtext('ArrowDestination/Filled')=='True'if arrow is None else arrow.get('Filled')=='True'
  if edge and head in('','None'):
   head='Fork'if edge['reversal_mV']>-40 else'Circle';audit['missing_native_chemical_heads_completed'].append({'edge':lid,'page':name,'reversal_mV':edge['reversal_mV'],'publication_head':head})
  if head=='Circle':ax.add_patch(Ellipse(pts[-1],6,6,facecolor=co if filled else'white',edgecolor=co,lw=.8,zorder=2))
  elif head not in('None',''):ax.add_patch(FancyArrowPatch(pts[-2],pts[-1],arrowstyle='->'if head=='Fork'else'-|>',mutation_scale=9,color=co,lw=.8,zorder=2))
  if changed and execution:audit['target_overlays'].append({'page':name,'edge_ID':lid,'source':edge['source'],'GUI_target_ID':edge['reference_target_ID'],'executed_target':edge['target'],'preserved_route_points':len(s.findall('Point'))})
  drawn.append(lid)
 for i,s in shapes.items():
  if i not in selected:continue
  n=local[i];x,y,w,h=[float(s.get(k))for k in('Left','Top','Width','Height')];fill=color(s.findtext('FillColor')or n.findtext('FillColor'));edgecolor=color(s.findtext('DrawColor')or n.findtext('DrawColor'),'#000000');kind=s.find('Shape').get('Style','')if s.find('Shape')is not None else''
  if kind in('Termination','Ellipse'):patch=Ellipse((x+w/2,y+h/2),w,h,facecolor=fill,edgecolor=edgecolor,lw=.8,zorder=3)
  elif 'OffPage'in n.findtext('ClassName',''):patch=Polygon([(x,y),(x+w,y),(x+w,y+h*.7),(x+w/2,y+h),(x,y+h*.7)],facecolor=fill,edgecolor=edgecolor,lw=.8,zorder=3)
  else:patch=Rectangle((x,y),w,h,facecolor=fill,edgecolor=edgecolor,lw=.8,zorder=3)
  ax.add_patch(patch);text=label(n,keep_side=name=='Commissural coordination')
  if name.endswith('_RG'):text=text.replace('RG-E\nIN','I-E').replace('RG-F\nIN','I-F').replace('RG-E','E').replace('RG-F','F')
  if publish_rg_inputs:text=text.replace('PF-E\nIN','I-E').replace('PF-F\nIN','I-F')
  ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=12.5,zorder=4,linespacing=1.0);points.extend([[x,y],[x+w,y+h]]);audit['publication_labels'][i]={'saved_name':n.findtext('Text'),'publication_label':text,'page':name,'resolved_neuron_ID':resolve(i),'class':n.findtext('ClassName')}
 def stub(x,y,w,h,text):
  ax.add_patch(Polygon([(x,y),(x+w,y),(x+w,y+h*.7),(x+w/2,y+h),(x,y+h*.7)],facecolor='#ffdf00',edgecolor='#000000',lw=.8,zorder=3));ax.text(x+w/2,y+h*.43,text,fontsize=12.5,ha='center',va='center',linespacing=1.0,zorder=4);points.extend([[x,y],[x+w,y+h]])
 if ii_targets:
  for i in selected:
   edge=next(e for e in edges.values()if e['source_ID']==i and' MN 'in e['target']);s=shapes[i];c=center(s);sw=105;sh=55;target_x=c[0]+(-22 if edge['target'].endswith('ext') else 22);xx=target_x-sw/2;yy=float(s.get('Top'))+float(s.get('Height'))+32;stub(xx,yy,sw,sh,'MN-E'if edge['target'].endswith('ext')else'MN-F');a=bound(s,[target_x,yy]);b=np.array([target_x,yy]);ax.add_patch(FancyArrowPatch(a,b,arrowstyle='->',mutation_scale=10,color='#333333',lw=.9,zorder=2));audit.setdefault('abstracted_crosspage_edges',[]).append({'edge_ID':edge['source_connection_ID'],'source_ID':i,'target_ID':edge['target_ID'],'source':edge['source'],'target':edge['target'],'view':'II-to-homonymous-MN off-page stub'})
 if publish_rg_inputs:
  for direction in('ext','flx'):
   e=next(e for e in edges.values()if e['source']=='L RG '+direction and e['target']=='L Hip PF '+direction);s=shapes[e['target_ID']];c=center(s);sw=80;sh=44;xx=c[0]-155 if direction=='ext'else c[0]+75;yy=c[1]-sh/2;stub(xx,yy,sw,sh,'RG-E'if direction=='ext'else'RG-F');a=np.array([xx+sw if direction=='ext'else xx,c[1]]);b=bound(s,a);ax.add_patch(FancyArrowPatch(a,b,arrowstyle='->',mutation_scale=10,color='#333333',lw=.9,zorder=2));audit.setdefault('abstracted_crosspage_edges',[]).append({'edge_ID':e['source_connection_ID'],'source_ID':e['source_ID'],'target_ID':e['target_ID'],'source':e['source'],'target':e['target'],'view':'RG-to-hip-PF off-page input stub'})
 points=np.asarray(points);lo=points.min(axis=0)-25;hi=points.max(axis=0)+25;ax.set_xlim(lo[0],hi[0]);ax.set_ylim(hi[1],lo[1]);ax.set_aspect('equal');ax.axis('off');ax.set_title(title,fontsize=13.5,fontweight='bold',loc='left',pad=10)
 audit['page_views'].append({'page':name,'title':title,'selected_nodes':list(selected),'drawn_connections':drawn,'preserved_geometry':True})
def newfig():return plt.figure(figsize=(11,8.5),layout='constrained')
figures=[]
def save(fig,stem):
 if stem=='W2L108_NATIVE_BILATERAL_OVERVIEW':
  fig.set_size_inches(9,5.5)
  for text in fig.findobj(matplotlib.text.Text):text.set_fontsize(max(10,text.get_fontsize()*.8))
  fig.get_layout_engine().set(h_pad=.02,w_pad=.02,hspace=.025,wspace=.025)
 for ext in('pdf','svg','png'):fig.savefig(OUT/(stem+'.'+ext),dpi=180)
 if stem=='W2L108_NATIVE_BILATERAL_OVERVIEW':
  for ext in('pdf','svg','png'):shutil.copyfile(OUT/(stem+'.'+ext),OUT/('W2L108_MAIN_SUBNETWORKS.'+ext))
  (OUT/'MAIN_LAYOUT.json').write_text(json.dumps({'width_in':9,'height_in':5.5,'minimum_font_pt':10,'native_internal_route_points_unchanged':True,'main_only_reflow':True},indent=2))
  if '--main-only'in sys.argv:print(OUT/'W2L108_MAIN_SUBNETWORKS.pdf');sys.exit(0)
 figures.append((stem,fig))
fig=newfig();gs=fig.add_gridspec(3,3,height_ratios=[3.4,1.5,1.7],width_ratios=[.95,1.15,1.2]);sub=gs[0,0].subgridspec(2,1)
render('LH_RG',fig.add_subplot(sub[0]),'(A) Left RG');render('RH_RG',fig.add_subplot(sub[1]),'Right RG')
hip=pages['LH_HipZ'];core={n.findtext('ID')for n in hip.findall('Nodes/Node')if'.Neurons.'in n.findtext('ClassName','')and any(x in n.findtext('Text','')for x in(' PF ',' MN '))}
render('LH_HipZ',fig.add_subplot(gs[0,1]),'(B) Hip PF/MN core',selected=core,publish_rg_inputs=True)
sub=gs[0,2].subgridspec(3,2)
for row,(joint,pg)in enumerate((('Hip','HipZ'),('Knee','Knee Motoneuron'),('Ankle','Anklez Motoneuron'))):
 for col,side in enumerate(('L','R')):
  name=side+'H_'+pg;sel={n.findtext('ID')for n in pages[name].findall('Nodes/Node')if' II relay'in n.findtext('Text','')and'.Neurons.'in n.findtext('ClassName','')};render(name,fig.add_subplot(sub[row,col]),('(C) 'if row==0 and col==0 else'')+side+' '+joint,selected=sel,ii_targets=True)
render('Commissural coordination',fig.add_subplot(gs[1,:]),'(D) Commissural pathways and contralateral heel input')
render('LH_Contact feedback',fig.add_subplot(gs[2,:]),'(E) Heel/toe contact subnetwork (left shown; right counterpart retained)')
fig.supxlabel('Retained phase-II (L/R): hip-F II → hip PF-F and RG-F; knee-F II → K/A PF-F and RG-F.',fontsize=12.5)
save(fig,'W2L108_NATIVE_BILATERAL_OVERVIEW')
fig=newfig();ax=fig.add_subplot(111);render('LH_HipZ',ax,'Expanded left hip: PF, MN, Ia/Ib, Renshaw, and separate II-E/F relays');save(fig,'W2L108_NATIVE_LEFT_HIP_DETAIL')
fig,axs=plt.subplots(1,2,figsize=(11,8.5),layout='constrained');render('LH_Knee Motoneuron',axs[0],'(A) Left knee: inherited MN/reflex subsystem');render('LH_Anklez Motoneuron',axs[1],'(B) Left ankle: executed PF target overlay');save(fig,'W2L108_NATIVE_LEFT_KNEE_ANKLE')
fig,axs=plt.subplots(2,1,figsize=(11,8.5),layout='constrained');render('LH_Contact feedback',axs[0],'(A) Left heel/toe contact conversion and ipsilateral targets');render('RH_Contact feedback',axs[1],'(B) Right heel/toe contact conversion and ipsilateral targets');save(fig,'W2L108_NATIVE_BILATERAL_CONTACT')
for joint,pg in(('HIP','HipZ'),('KNEE','Knee Motoneuron'),('ANKLE','Anklez Motoneuron')):
 fig,axs=plt.subplots(1,2,figsize=(11,8.5),layout='constrained')
 for side,ax in zip(('L','R'),axs):
  name=side+'H_'+pg;p=pages[name];local={n.findtext('ID'):n for n in p.findall('Nodes/Node')};sel={i for i,n in local.items()if' II relay'in n.findtext('Text','')and'.Neurons.'in n.findtext('ClassName','')};assert len(sel)==2
  for e in executed['edges']:
   if e['source_ID']in sel:sel.add(e['target_ID'])
  for repeat in range(2):
   for l in p.findall('Links/Link'):
    if l.findtext('DestinationID')in sel and'Adapter'in nodes[l.findtext('OriginID')].findtext('ClassName',''):sel.add(l.findtext('OriginID'))
  render(name,ax,'('+side+') '+joint.lower()+' II-E and II-F to homonymous MN',selected=sel)
 save(fig,'W2L108_NATIVE_'+joint+'_BILATERAL_II')
with PdfPages(OUT/'W2L108_NATIVE_SCHEMATIC_ATLAS.pdf')as pdf:
 for stem,fig in figures:pdf.savefig(fig);plt.close(fig)
with(OUT/'EXECUTED_EDGES.csv').open('w',newline='',encoding='utf-8')as f:
 fields=list(next(iter(edges.values())).keys());wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader();wr.writerows(edges.values())
groups=defaultdict(list)
for e in edges.values():groups[(e['reference_synapse_type_ID'],e['reference_synapse_type_name'],e['gmax_uS'],e['reversal_mV'],e['threshold_mV'],e['saturation_mV'])].append(e)
rows=[]
for k,ee in sorted(groups.items(),key=lambda kv:(kv[0][1],kv[0][2])):
 t=types[k[0]];g=t.find('MaxSynapticConductance');rows.append({'TypeID':k[0],'NativeName':k[1],'NativeMaxConductance_uS':None if g is None else float(g.get('Actual'))*1e6,'ExecutedMaxConductance_uS':k[2],'Equilibrium_mV':k[3],'Threshold_mV':k[4],'Saturation_mV':k[5],'EdgeCount':len(ee),'ConnectionIDs':';'.join(e['source_connection_ID']for e in ee)})
with(OUT/'EXECUTED_NAMED_TYPES.csv').open('w',newline='',encoding='utf-8')as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
(OUT/'SCHEMATIC_PROVENANCE.json').write_text(json.dumps(audit,indent=2));(OUT/'EXECUTED_NAMED_TYPES.json').write_text(json.dumps(rows,indent=2));print(OUT)
for src,dst in(('W2L108_NATIVE_BILATERAL_OVERVIEW','W2L108_MAIN_SUBNETWORKS'),('W2L108_NATIVE_LEFT_HIP_DETAIL','W2L108_LEFT_SUBNETWORKS')):
 for ext in('pdf','svg','png'):shutil.copyfile(OUT/(src+'.'+ext),OUT/(dst+'.'+ext))
with(OUT/'NEURON_NAMES.csv').open('w',newline='',encoding='utf-8')as f:
 wr=csv.DictWriter(f,fieldnames=['ID','NativeName','PublicationLabel']);wr.writeheader()
 for n in executed['neurons']:wr.writerow({'ID':n['ID'],'NativeName':n['name'],'PublicationLabel':label(nodes[n['ID']])})
with PdfPages(OUT/'EXECUTED_NAMED_TYPES.pdf')as pdf:
 for start in range(0,len(rows),18):
  chunk=rows[start:start+18];fig,ax=plt.subplots(figsize=(11,8.5));ax.axis('off');data=[]
  for v in chunk:data.append(['\n'.join(textwrap.wrap(v['NativeName'],33)),f"{v['NativeMaxConductance_uS']:.6g}",f"{v['ExecutedMaxConductance_uS']:.6g}",f"{v['Equilibrium_mV']:g}",f"{v['Threshold_mV']:g}",f"{v['Saturation_mV']:g}",str(v['EdgeCount'])])
  table=ax.table(cellText=data,colLabels=['Saved named type','GUI max G\n(uS)','Executed G\n(uS)','E rev\n(mV)','Threshold\n(mV)','Saturation\n(mV)','Edges'],colWidths=[.40,.11,.11,.10,.10,.10,.08],loc='center',cellLoc='left');table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,2.15)
  for (rr,cc),cell in table.get_celld().items():
   cell.set_edgecolor('#bbbbbb');cell.set_linewidth(.4)
   if rr==0:cell.set_facecolor('#eeeeee')
  ax.set_title('Executed named synapse types: GUI baseline and current working gains',fontsize=12,pad=12);pdf.savefig(fig,bbox_inches=None);plt.close(fig)
import runpy
runpy.run_path(str(Path(__file__).with_name('verify_native_schematics.py')))
