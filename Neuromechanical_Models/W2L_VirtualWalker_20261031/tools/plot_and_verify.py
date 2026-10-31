from pathlib import Path
from lxml import etree as E
import numpy as np,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
base=Path('Neuromechanical_Models');fig,axs=plt.subplots(3,2,figsize=(12,9),sharex='col');details=[]
for c,(folder,title) in enumerate([('W2L_VirtualWalker_20261031','Original controller'),('W2L_VirtualWalker_Newer_20261031','II + contact + bilateral RG')]):
 w=base/folder;z=np.load(w/'verification.npz');a=z['data'];labels=list(z['labels']);t=a[:,1];meta=json.loads((w/'BUILD.json').read_text());get=lambda name:a[:,labels.index(name)]
 root=np.column_stack([get('root_'+x) for x in 'xyz']);left=np.column_stack([get('Left hip wing_'+x) for x in 'xyz']);right=np.column_stack([get('Right hip wing_'+x) for x in 'xyz'])
 Ry=(root-(left+right)/2)/.04;Rz=(left-right)/.5;Rx=np.cross(Ry,Rz);height=(np.abs(Rx[:,1])*.12+np.abs(Ry[:,1])*.025+np.abs(Rz[:,1])*.30)/2
 wl=left[:,1]-height;wr=right[:,1]-height;threshold=meta['rail_top_y']+.0005;touch=t[np.minimum(wl,wr)<=threshold]
 contacts={n:float(t[np.where(get(n)>0)[0][0]]) if np.any(get(n)>0) else None for n in labels if n.endswith('_contacts')}
 preservation={}
 for source,sha in meta['source'].items():
  src=Path(source);out=w/('W2L_VirtualWalker_Newer' if c else 'W2L_VirtualWalker');out=out.with_suffix('.aproj') if src.suffix=='.aproj' else Path(str(out)+'_Standalone.asim')
  rs=E.parse(str(src));ro=E.parse(str(out));same={}
  for tag in ['NeuralModules','NervousSystem','Adapters','BehavioralSystem']:
   ns=rs.findall('.//'+tag);no=ro.findall('.//'+tag)
   if ns:same[tag]=len(ns)==len(no) and all(E.tostring(x)==E.tostring(y) for x,y in zip(ns,no))
  preservation[str(src)]=dict(source_unchanged=hashlib.sha256(src.read_bytes()).hexdigest()==sha,payloads=same)
  assert all(same.values()) and preservation[str(src)]['source_unchanged']
 summary=dict(model=title,root_y_range_m=[float(get('root_y').min()),float(get('root_y').max())],end_root_y_m=float(get('root_y')[-1]),root_x_displacement_m=float(get('root_x')[-1]-get('root_x')[0]),first_approximate_wing_contact_s=float(touch[0]) if len(touch) else None,first_foot_sensor_contacts_s=contacts,wing_surface_min_y_m=[float(wl.min()),float(wr.min())],rail_top_y_m=meta['rail_top_y'],preservation=preservation)
 (w/'VERIFICATION.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');details.append(summary)
 support_y=meta['root_initial_xyz'][1]-meta['wing_to_rail_gap']
 ax=axs[0,c];ax.plot(t,get('root_y'),color='#4477AA',label='Root height');ax.axhline(support_y,color='#777777',ls='--',label='Level support height');ax.set(title=title,ylabel='Root height (m)',ylim=(support_y-.01,meta['root_initial_xyz'][1]+.005));ax.legend(loc='upper right',fontsize=8)
 ax=axs[1,c]
 for name,col,style in [('hip_L_angle','#4477AA','-'),('hip_R_angle','#4477AA','--'),('knee_L_angle','#EE6677','-'),('knee_R_angle','#EE6677','--')]:ax.plot(t,np.rad2deg(get(name)),color=col,ls=style,label=name.replace('_angle','').replace('_',' '))
 ax.set(ylabel='Joint change from initial pose (deg)');ax.legend(ncol=2,fontsize=8)
 ax=axs[2,c]
 for name,col,style in [('L Hip MN flx','#4477AA','-'),('R Hip MN flx','#4477AA','--'),('L Hip MN ext','#CCBB44','-'),('R Hip MN ext','#CCBB44','--')]:ax.plot(t,get(name)*1000,color=col,ls=style,label=name)
 ax.set(xlabel='Time (s)',ylabel='Hip MN voltage (mV)');ax.legend(ncol=2,fontsize=8)
 for ax in axs[:,c]:ax.grid(alpha=.2);ax.set_xlim(0,10)
fig.suptitle('Passive hip-wing support: Root is free, rails provide mechanical assistance',fontsize=14)
fig.tight_layout(rect=(0,0,1,.97));fig.savefig(base/'W2L_VirtualWalker_20261031/verification_comparison.png',dpi=160)
print(json.dumps([{k:v for k,v in d.items() if k!='preservation'} for d in details],indent=2))
