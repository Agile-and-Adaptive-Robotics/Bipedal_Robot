from pathlib import Path
import json,hashlib,re
import numpy as np
from lxml import etree as E
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Rectangle,Polygon,FancyArrowPatch
import mujoco
from PIL import Image
R=Path(r'D:\GitHub\Bipedal_Robot');W=R/'Neuromechanical_Models/Deng_Biarticular_Rat';OUT=Path(__file__).parent/'rat_native';OUT.mkdir(exist_ok=True);P=W/'Biarticular_net.aproj';tree=E.parse(str(P));nodes={n.findtext('ID'):n for n in tree.findall('.//NervousSystem//Node')};pages={n.findtext('Text').strip():n for n in nodes.values()if n.find('DiagramXml')is not None};audit={'source_project':str(P),'project_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'source_xml':str(W/'Deng_Biarticular_Rat.xml'),'source_xml_sha256':hashlib.sha256((W/'Deng_Biarticular_Rat.xml').read_bytes()).hexdigest(),'source_scope':'Public2020 bilateral-hindlimb precursor, not verified final2023. Frozen pelvis. Static imported mechanics plus actual native diagram pages; not a new coupled simulation.','labels':{},'pages':[],'native_font_pt':12.5,'print_font_pt_at9inch':12.5*9/11}
plt.rcParams.update({'font.family':'Arial','font.size':12.5,'axes.titlesize':13.5,'svg.fonttype':'none','pdf.fonttype':42})
def color(t,default='#ffffff'):
 try:return '#'+format(int(t)&0xffffff,'06x')
 except:return default
def short(n,page=''):
 s=re.sub(r'\s+',' ',n.findtext('Text','').strip());s=s.replace('LH_','').replace('RH_','').replace('LH ','').replace('RH ','');s=s.replace('EXT','E').replace('FLX','F').replace('Ext','E').replace('Flx','F').replace('HipZ','Hip').replace(' IN ',' I-').replace(' PF ',' PF-');s=s.replace('RG IN E','RG I-E').replace('RG IN F','RG I-F');s=s.replace(' IN',' I');s=s.replace('BFP&RF','BFP/RF').replace('BFA&IP','BFA/IP').replace('VA&GA','VA/GA').replace('SO&TA','SO/TA')
 if n.findtext('ClassName','').endswith('Subsystem'):return s
 if 'Adapter'in n.findtext('ClassName',''):return'map'
 if'OffPage'in n.findtext('ClassName',''):
  if'Hip'in s and'II'in s:return'II-E'if'-E'in s else'II-F'
  if'IP'in s and'Ia'in s:return'Ia-F'
  if'SO'in s and'Ib'in s:return'Ib-E'
  if'RG'in s:return'E'if'E'in s else'F'
  if'PF'in s:return('Hip'if'Hip'in s else'Knee')+'\nPF-'+('E'if s.endswith('E')else'F')
 if page=='LH_RG':return'I-E'if'I-E'in s else'I-F'if'I-F'in s else'E'if s.endswith('E')else'F'
 if page=='LH_Hip':return'I-E'if'PF-I E'in s or'PF-I-E'in s else'I-F'if'PF-I F'in s or'PF-I-F'in s else'PF-E'if s.endswith('E')else'PF-F'
 if page=='Biar':return'BFP I-E'if'BFP'in s else'RF I-E'
 return s
def draw(name,ax,title):
 p=pages[name];flow=E.fromstring(p.findtext('DiagramXml').encode()).find('.//AddFlow');seq=[s for s in flow if s.tag in('Node','Link')];local={n.findtext('ID'):n for n in p.findall('Nodes/Node')};links={l.findtext('ID'):l for l in p.findall('Links/Link')};shapes={s.findtext('Tag'):s for s in seq if s.tag=='Node'};points=[]
 for i,s in shapes.items():
  n=local[i];cls=n.findtext('ClassName','');scale=1.4 if '.Neurons.'in cls or 'OffPage'in cls else 1.0
  if'OffPage'in cls and'RG'in n.findtext('Text',''):scale=1.0
  x,y,w,h=[float(s.get(k))for k in('Left','Top','Width','Height')];s.attrib.update({'Left':str(x-(scale-1)*w/2),'Top':str(y-(scale-1)*h/2),'Width':str(w*scale),'Height':str(h*scale)})
 def center(s):return np.array([float(s.get('Left'))+float(s.get('Width'))/2,float(s.get('Top'))+float(s.get('Height'))/2])
 def bound(s,other):
  c=center(s);v=np.array(other)-c;hs=[float(s.get('Width'))/2,float(s.get('Height'))/2];k=min([hs[j]/abs(v[j])for j in(0,1)if abs(v[j])>1e-9]or[0]);return c+k*v
 for s in seq:
  if s.tag!='Link':continue
  l=links[s.findtext('Tag')];a=seq[int(s.get('Org'))];b=seq[int(s.get('Dst'))];pts=[center(a)]+[np.array([float(q.get('X')),float(q.get('Y'))])for q in s.findall('Point')]+[center(b)];pts[0]=bound(a,pts[1]);pts[-1]=bound(b,pts[-2]);points.extend(pts);co=color(s.findtext('DrawColor')or l.findtext('DrawColor'),'#333333');ax.plot(*np.array(pts).T,color=co,lw=.8,zorder=1);h=s.find('ArrowDst');head=l.findtext('ArrowDestination/Style','')if h is None else h.get('Head','');filled=l.findtext('ArrowDestination/Filled')=='True'if h is None else h.get('Filled')=='True'
  if head=='Circle':ax.add_patch(Ellipse(pts[-1],6,6,facecolor=co if filled else'white',edgecolor=co,lw=.8,zorder=2))
  elif head not in('','None'):ax.add_patch(FancyArrowPatch(pts[-2],pts[-1],arrowstyle='->'if head=='Fork'else'-|>',mutation_scale=9,color=co,lw=.8,zorder=2))
 for i,s in shapes.items():
  n=local[i];x,y,w,h=[float(s.get(k))for k in('Left','Top','Width','Height')];fill=color(s.findtext('FillColor')or n.findtext('FillColor'));co=color(s.findtext('DrawColor')or n.findtext('DrawColor'),'#000000');kind=s.find('Shape').get('Style','')if s.find('Shape')is not None else''
  if kind in('Ellipse','Termination'):patch=Ellipse((x+w/2,y+h/2),w,h,facecolor=fill,edgecolor=co,lw=.8,zorder=3)
  elif 'OffPage'in n.findtext('ClassName',''):patch=Polygon([(x,y),(x+w,y),(x+w,y+h*.7),(x+w/2,y+h),(x,y+h*.7)],facecolor=fill,edgecolor=co,lw=.8,zorder=3)
  else:patch=Rectangle((x,y),w,h,facecolor=fill,edgecolor=co,lw=.8,zorder=3)
  ax.add_patch(patch);ss=short(n,name);ax.text(x+w/2,y+h/2,ss,fontsize=12.5,ha='center',va='center',zorder=4);points.extend([[x,y],[x+w,y+h]]);audit['labels'][i]={'source_name':n.findtext('Text'),'display_label':ss,'page':name}
 pp=np.array(points);lo=pp.min(0)-30;hi=pp.max(0)+30;ax.set_xlim(lo[0],hi[0]);ax.set_ylim(hi[1],lo[1]);ax.set_aspect('equal');ax.axis('off');ax.set_title(title,loc='left',fontsize=13.5,fontweight='bold');audit['pages'].append({'name':name,'node_count':len(local),'link_count':len(links),'native_node_centers_and_control_points_retained':True,'publication_neural_glyph_scale':1.4,'publication_scope':'Source positions/routes preserved; some glyphs enlarged around fixed centers for12.5pt text.'})
# Static view of the imported source mechanics; no mj_step call.
mx=E.parse(str(W/'Deng_Biarticular_Rat.xml'));m=mujoco.MjModel.from_xml_path(str(W/'Deng_Biarticular_Rat.xml'));d=mujoco.MjData(m);mujoco.mj_forward(m,d);ids=[i for i in range(m.ngeom)if(mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,i)or'').startswith(('geom_Pelvis','geom_LH','geom_RH'))];c=d.geom_xpos[ids].mean(0);extent=max(np.linalg.norm(d.geom_xpos[i]-c)+max(m.geom_size[i])for i in ids);eye=c+np.array([2.1,.9,2.2])*extent;f=c-eye;f/=np.linalg.norm(f);right=np.cross(f,[0,1,0]);right/=np.linalg.norm(right);up=np.cross(right,f);wb=mx.getroot().find('worldbody');E.SubElement(wb,'camera',name='publication_world',pos=' '.join(map(str,eye)),xyaxes=' '.join(map(str,np.r_[right,up])),fovy='40')
for g in mx.findall('.//geom'):
 name=g.get('name','')
 if name.startswith('geom_LH'):g.set('rgba','.0 .467 .733 1')
 elif name.startswith('geom_RH'):g.set('rgba','.933 .467 .2 1')
 elif name=='geom_Pelvis':g.set('rgba','.25 .28 .32 1')
for tendon in mx.findall('.//tendon/spatial'):tendon.set('width','.0003')
m=mujoco.MjModel.from_xml_string(E.tostring(mx,encoding='unicode'));d=mujoco.MjData(m);mujoco.mj_forward(m,d);rr=mujoco.Renderer(m,480,480);opt=mujoco.MjvOption();opt.sitegroup[:]=0;rr.update_scene(d,camera='publication_world',scene_option=opt);body=rr.render().copy();Image.fromarray(body).save(OUT/'RAT_SOURCE_STATIC_BODY.png');del rr
fig=plt.figure(figsize=(11,8.5),layout='constrained');gs=fig.add_gridspec(2,3);a=fig.add_subplot(gs[0,0]);a.imshow(body);a.axis('off');a.set_title('(A) Bilateral hindlimb model',fontsize=13.5,fontweight='bold');draw('LH',fig.add_subplot(gs[0,1:]),'(B) Native left-side subnetwork hierarchy');draw('LH_RG',fig.add_subplot(gs[1,0]),'(C) Rhythm generation');draw('LH_Hip',fig.add_subplot(gs[1,1]),'(D) Hip pattern formation');draw('Biar',fig.add_subplot(gs[1,2]),'(E) Native biarticular coupling')
fig.set_size_inches(9,5.5)
for text in fig.findobj(matplotlib.text.Text):text.set_fontsize(max(10,text.get_fontsize()*.8))
fig.get_layout_engine().set(h_pad=.02,w_pad=.02,hspace=.025,wspace=.025)
audit['native_font_pt']=10;audit['print_font_pt_at9inch']=10;audit['canvas_inches']=[9,5.5]
for ext in('pdf','svg','png'):fig.savefig(OUT/('DENG_RAT_NATIVE_SUBNETWORKS.'+ext),dpi=190)
plt.close(fig);audit['source_project_unchanged']=hashlib.sha256(P.read_bytes()).hexdigest()==audit['project_sha256'];audit['source_xml_unchanged']=hashlib.sha256((W/'Deng_Biarticular_Rat.xml').read_bytes()).hexdigest()==audit['source_xml_sha256'];(OUT/'PROVENANCE.json').write_text(json.dumps(audit,indent=2));print(OUT)
