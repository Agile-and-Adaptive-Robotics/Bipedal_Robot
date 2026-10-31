from pathlib import Path
from lxml import etree as E
import math,json,csv,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Rectangle,Polygon,FancyBboxPatch,FancyArrowPatch
R=Path(r'D:\GitHub\Bipedal_Robot');P=R/'Neuromechanical_Models/Walker_2_Layer_CPG_BilateralRG/Walker_2_Layer_CPG_BilateralRG.aproj';OUT=R/'Documentation/Reports and Papers/Dissertation/Notes/sns_rebuild_session1_20261030/user_selected_subnetwork';OUT.mkdir(exist_ok=True)
t=E.parse(str(P));pg=t.find('.//NervousSystem//DiagramXml').getparent();nodes={n.findtext('ID'):n for n in pg.findall('Nodes/Node')};links={n.findtext('ID'):n for n in pg.findall('Links/Link')};flow=E.fromstring(pg.findtext('DiagramXml').encode()).find('.//AddFlow');seq=[n for n in flow if n.tag in ('Node','Link')];shapes={n.findtext('Tag'):n for n in seq if n.tag=='Node'};types={x.findtext('ID'):x for x in t.findall('.//SynapseTypes/Link')}
sel=set()
for i,s in shapes.items():
 n=nodes[i];name=n.findtext('Text','').strip();x,y=float(s.get('Left')),float(s.get('Top'))
 if name.startswith(('L Hip','L RG'))or name in ['L c1','R V3','L heel contact']or (x<800 and y<300):sel.add(i)
# Keep the two motor-output adapters between the displayed MN and hip muscles.
for e in links.values():
 a,b=e.findtext('OriginID'),e.findtext('DestinationID')
 if a in sel and 'NodeToPhysicalAdapter'in nodes.get(b,E.Element('Node')).findtext('ClassName',''):sel.add(b)

def center(s):return np.array([float(s.get('Left'))+float(s.get('Width'))/2,float(s.get('Top'))+float(s.get('Height'))/2])
coords={i:center(shapes[i])for i in sel}
# A printed-label map removes only repeated L/Hip prefixes; full saved labels accompany the drawing.
def label(i):
 n=nodes[i];s=n.findtext('Text','').strip();cls=n.findtext('ClassName','')
 if s=='HeelContact':return 'Heel'
 if s=='L heel contact':return 'Heel\ncontact'
 if s=='9':return '9'
 if 'Adapter'in cls:
  if s in ['2','3']:return s
  if 'NodeToPhysical'in cls:return 'Drive'
  return 'II' if ' II 'in s else 'Ia' if ' Ia 'in s else 'Ib'
 if s in ['R V3','L c1']:return s
 if s.startswith('L RG '):
  direction='E'if ' ext'in s else'F'
  return ('RG\nI-'if s.endswith(' IN')else'RG-')+direction
 if s.startswith('L Hip '):
  s=s[6:].replace('ext','E').replace('flx','F')
  if s.startswith('PF '):return 'PF\nI-'+('E'if 'E'in s else'F')if s.endswith(' IN')else'PF-'+s[-1]
  if s.startswith('MN '):return 'RE-'+s.split()[1]if s.endswith(' RE')else'MN-'+s.split()[1]
  if s in ['E','F']:return 'Hip '+s
  if s.endswith(' IN'):return s.split()[1]+'-'+s.split()[0]+'\nIN'
  if len(s.split())==2:return s.split()[1]+'-'+s.split()[0]
 return s
plt.rcParams.update({'font.family':'Arial','font.size':10.5,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(6.5,6.5));ax=fig.add_axes([.035,.035,.93,.93]);ax.axis('off')
# Extent follows the supplied selection. No edge to an omitted node is drawn.
allc=np.array(list(coords.values()));lo=allc.min(0)-43;hi=allc.max(0)+43
ax.set_xlim(lo[0],hi[0]);ax.set_ylim(hi[1],lo[1]);ax.set_aspect('equal')
size={}
for i in sel:
 cls=nodes[i].findtext('ClassName','');size[i]=(46,43)if 'Adapter'in cls else (65,45)if cls.endswith('.Muscle')or cls.endswith('.StretchReceptor')or cls.endswith('.RigidBody')else(57,49)
for i in sel:
 if nodes[i].findtext('Text','').strip()=='L heel contact':size[i]=(78,52)
def boundary(i,other):
 c=coords[i];v=np.array(other)-c;w,h=size[i];cls=nodes[i].findtext('ClassName','')
 if '.Neurons.'in cls:
  k=1/math.sqrt((v[0]/(w/2))**2+(v[1]/(h/2))**2);return c+k*v
 k=min((w/2)/abs(v[0])if abs(v[0])>1e-8 else 1e8,(h/2)/abs(v[1])if abs(v[1])>1e-8 else 1e8);return c+k*v
records=[];omitted=[]
for q in seq:
 if q.tag!='Link':continue
 e=links[q.findtext('Tag')];a,b=e.findtext('OriginID'),e.findtext('DestinationID')
 if a not in sel or b not in sel:
  if a in sel or b in sel:omitted.append(e.findtext('ID'))
  continue
 pts=[coords[a]]+[np.array([float(z.get('X')),float(z.get('Y'))])for z in q.findall('Point')]+[coords[b]]
 rerouted=any(np.any(z<lo)|np.any(z>hi)for z in pts)
 if rerouted:pts=[coords[a],np.array([coords[b][0],coords[a][1]]),coords[b]]
 # Remove a control point concealed within an enlarged endpoint glyph.
 def inside(i,z):c=coords[i];w,h=size[i];return abs(z[0]-c[0])<w/2 and abs(z[1]-c[1])<h/2
 while len(pts)>2 and inside(a,pts[1]):pts.pop(1)
 while len(pts)>2 and inside(b,pts[-2]):pts.pop(-2)
 pts[0]=boundary(a,pts[1]);pts[-1]=boundary(b,pts[-2])
 if len(pts)>2:ax.plot(*np.array(pts[:-1]).T,color='#222222',lw=.68,zorder=1)
 arrow=q.find('ArrowDst');head=e.findtext('ArrowDestination/Style','')if arrow is None else arrow.get('Head','');fill=e.findtext('ArrowDestination/Filled')=='True'if arrow is None else arrow.get('Filled')=='True'
 if head=='Circle':
  ax.plot([pts[-2][0],pts[-1][0]],[pts[-2][1],pts[-1][1]],color='#222222',lw=.68,zorder=1);ax.add_patch(Ellipse(pts[-1],9,9,facecolor='#222222'if fill else'white',edgecolor='#222222',lw=.7,zorder=4))
 elif head not in ['','None']:ax.add_patch(FancyArrowPatch(pts[-2],pts[-1],arrowstyle='->'if head=='Fork'else'-|>',mutation_scale=8,color='#222222',lw=.68,zorder=2,shrinkA=0,shrinkB=0))
 else:ax.plot([pts[-2][0],pts[-1][0]],[pts[-2][1],pts[-1][1]],color='#222222',lw=.68,zorder=1)
 st=types.get(e.findtext('SynapticTypeID'));record={'connection_ID':e.findtext('ID'),'source_ID':a,'target_ID':b,'source':nodes[a].findtext('Text'),'target':nodes[b].findtext('Text'),'native_route_points':len(q.findall('Point')),'rerouted_inside_selection':bool(rerouted)}
 if st is not None:
  record['type_name']=st.findtext('Name');record['type_ID']=st.findtext('ID')
  for key in ['EquilibriumPotential','MaxSynapticConductance','PreSynapticThreshold','PreSynapticSaturationLevel']:
   z=st.find(key);record[key]=None if z is None else z.get('Actual')
 records.append(record)
for i in sel:
 n=nodes[i];cls=n.findtext('ClassName','');x,y=coords[i];w,h=size[i]
 if '.Neurons.'in cls:patch=Ellipse((x,y),w,h,facecolor='#bce3f8',edgecolor='#222',lw=.8,zorder=5)
 elif 'Adapter'in cls:patch=Polygon([(x-w*.4,y-h/2),(x+w*.35,y-h/2),(x+w/2,y),(x+w*.35,y+h/2),(x-w*.4,y+h/2),(x-w/2,y)],facecolor='#f5a092',edgecolor='#222',lw=.8,zorder=5)
 elif cls.endswith('.RigidBody'):patch=Polygon([(x-w/2,y-h/2),(x+w/2,y-h/2),(x,y+h/2)],facecolor='#c7ef96',edgecolor='#222',lw=.8,zorder=5)
 else:patch=FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0,rounding_size=20',facecolor='#f6eadc',edgecolor='#222',lw=.8,zorder=5)
 ax.add_patch(patch);ax.text(x,y,label(i),fontsize=10.5,ha='center',va='center',linespacing=.93,zorder=6)
for ext in ['pdf','svg','png']:fig.savefig(OUT/('w2l_selected_left_subnetwork.'+ext),dpi=220)
plt.close(fig)
fields=list(dict.fromkeys(k for r in records for k in r))
with(OUT/'INTERNAL_CONNECTIONS.csv').open('w',encoding='utf-8',newline='')as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(records)
with(OUT/'LABEL_MAP.csv').open('w',encoding='utf-8',newline='')as f:w=csv.DictWriter(f,fieldnames=['ID','saved_name','display_name']);w.writeheader();w.writerows({'ID':i,'saved_name':nodes[i].findtext('Text'),'display_name':label(i).replace('\n',' ')}for i in sorted(sel))
(OUT/'PROVENANCE.json').write_text(json.dumps({'machine':'easteregg2','source':str(P),'source_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'nodes_shown':len(sel),'internal_connections':len(records),'connections_to_omitted_nodes_removed':len(omitted),'outside_route_excursions_rerouted':sum(r['rerouted_inside_selection']for r in records),'positions':'Native centers retained; symbols resized to hold readable labels.','scope':'Exact native source layout corresponding to the user supplied left-RG/hip crop. Not a display of every neuron in the later 108-neuron execution configuration.','image_icons':'No raster muscle/spindle art; simple labeled vector capsules.','font_pt':10.5,'size_inches':[6.5,6.5]},indent=2),encoding='utf-8')
print(OUT);print('Nodes',len(sel),'internal links',len(records),'external links removed',len(omitted))

