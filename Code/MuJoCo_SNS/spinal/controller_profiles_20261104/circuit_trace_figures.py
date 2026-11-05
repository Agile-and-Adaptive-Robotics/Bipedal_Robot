"""Reconstruct a fitted D circuit and draw its routes beside saved neural traces."""
import argparse,csv,hashlib,json,os,sys,time
from pathlib import Path
for key in('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
HERE=Path(__file__).resolve().parent;sys.path[:0]=[str(HERE),str(HERE.parent)]
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch,Circle
import motor_adapter as M
from muscle_map import classify,GROUPS
from diagram_routes import render_routes
COLORS=['#4477AA','#EE6677','#228833','#CCBB44','#66CCEE','#AA3377']

def main():
 p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--ledger-output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True);a.ledger_output.mkdir(parents=True,exist_ok=True)
 final=json.loads((a.fit/'FINAL.json').read_text());dm=json.loads((a.fit/'MANIFEST.json').read_text());cf=Path(dm['fit']);cm=json.loads((cf/'MANIFEST.json').read_text());family=cm['configuration']['family']
 ready=M.handle_request(dict(op='preflight',fit=str(cf),fit_sha256=hashlib.sha256((cf/'FINAL.json').read_bytes()).hexdigest()))
 if not ready['ready_for_optimization']:raise ValueError('C reconstruction failed')
 checked=M.handle_request(dict(op='evaluate',params=final['parameters']))
 if abs(checked['value']-final['best']['value'])>1e-10:raise ValueError('D reconstruction differs from saved result')
 r=M.INSTANCE;c=r.c;names=r.names;runtime=r.net.command_runtime;slots={int(slot):j for j,slot in enumerate(runtime.slots)}
 saved_trace=np.load(a.fit/'BEST_TRACE.npz');trace_match={key:bool(np.array_equal(saved_trace[key],value))for key,value in(('rg',r.rg_trace),('pf',r.pf_trace),('mn',r.output),('desired',r.desired))}
 if not all(trace_match.values()):raise ValueError('Circuit reconstruction does not reproduce saved traces exactly')
 edges=[]
 for i,(pre,post)in enumerate(zip(c.pre,c.post)):
  row=dict(source=names[pre],target=names[post],snapshot_gmax_uS=float(c.gmax[i]),inherited_C_gmax_uS=float(r.base_g[i]),reversal_mV=float(c.reversal[i]),threshold_mV=float(c.elo[i]),saturation_mV=float(c.elo[i]+c.erange[i]),sign='exc'if c.reversal[i]>0 else'inh')
  if i in slots:
   j=slots[i];row.update(dynamic_central_feedback=True,afferent_group=runtime.groups[j],reference_budget_gmax_uS=float(runtime.base[j]),shared_gain=float(runtime.group_values[runtime.group_indices[j]]),stance_gated=bool(runtime.ib[j]),gate='alpha-dependent presynaptic gain; Ib additionally multiplied by measured stance load')
  edges.append(row)
 def active(e):return e['inherited_C_gmax_uS']>0 or e.get('reference_budget_gmax_uS',0)>0
 edges=[e for e in edges if active(e)]
 pfs=[n for n in names if n.endswith('_r')and n.startswith('PF_')and not n.startswith(('PF_IN','PF_HUNT'))]
 acts=[act for act in r.acts if act.endswith('_r')];acts.sort(key=lambda act:(GROUPS.index(classify(act).groups[0]),act))
 matrix=np.zeros((len(acts),len(pfs)))
 for e in edges:
  if e['source']in pfs and e['target'].startswith('MN_')and e['target'][3:]in acts:matrix[acts.index(e['target'][3:]),pfs.index(e['source'])]+=e['snapshot_gmax_uS']
 groups=[g for g in GROUPS if any(classify(act).groups[0]==g for act in acts)];aggregated=np.zeros((len(groups),len(pfs)))
 for i,act in enumerate(acts):aggregated[groups.index(classify(act).groups[0])]+=matrix[i]
 incoming={pf:[e for e in edges if e['target']==pf and e['source'].startswith(('RG_','Ia_','II_','Ib_'))]for pf in pfs}
 # Local autogenic Ib route: direct or through an Ib relay, never infer it
 # from a flexor/extensor phase label or traverse recurrent MN feedback.
 feedback=[]
 for act in acts:
  source='Ib_'+act;mn='MN_'+act;local=[]
  for first in edges:
   if first['source']!=source:continue
   if first['target']==mn:local.append([first])
   elif first['target'].startswith(('IBIN_','IBEXC_')):
    for second in edges:
     if second['source']==first['target']and second['target']==mn:local.append([first,second])
  central=[e for e in edges if e['source']==source and e['target'].startswith(('RG_','PF_'))]
  feedback.append(dict(muscle=act,groups=list(classify(act).groups),autogenic_routes=[dict(path=[source]+[e['target']for e in path],sign='exc'if np.prod([1 if e['sign']=='exc'else-1 for e in path])>0 else'inh',snapshot_gmax_uS=[e['snapshot_gmax_uS']for e in path])for path in local],central_targets=[dict(target=e['target'],sign=e['sign'],stance_gated=e.get('stance_gated',False))for e in central]))
 mirror_lookup={(e['source'],e['target']):e['snapshot_gmax_uS']for e in edges};mirror_differences=[];missing_mirrors=[]
 for e in edges:
  if e['source']in pfs and e['target'].startswith('MN_')and e['target'].endswith('_r'):
   other=(e['source'][:-1]+'l',e['target'][:-1]+'l')
   if other in mirror_lookup:mirror_differences.append(abs(e['snapshot_gmax_uS']-mirror_lookup[other]))
   else:missing_mirrors.append(other)
 summary=dict(family=family,fit=str(a.fit),C_fit=str(cf),fresh_D_objective_matches=True,saved_traces_exact=trace_match,D_parameters=final['parameters'],PF_cells=pfs,PF_incoming=incoming,Ib_per_muscle=feedback,PF_MN_mirror_max_difference_uS=max(mirror_differences,default=0.),PF_MN_missing_mirrors=missing_mirrors,scope='Actual fitted offline D circuit. Snapshot conductances are at final replay time; zero instantaneous central Ib gain is gating, not missing topology. Right side displayed; all bilateral edges in external ledger.')
 ledger=a.ledger_output/'FITTED_EDGE_LEDGER.json'
 if ledger.exists()and json.loads(ledger.read_text())!=edges:raise ValueError('Existing ledger differs; preserve it and use a new owned output folder')
 ledger.write_text(json.dumps(edges,indent=2));(a.output/'CIRCUIT_SUMMARY.json').write_text(json.dumps(summary,indent=2))
 for filename,array,labels in(('PF_TO_MUSCLE.csv',matrix,acts),('PF_TO_GROUP.csv',aggregated,groups)):
  with(a.output/filename).open('w',newline='',encoding='utf-8')as f:
   writer=csv.writer(f);writer.writerow(['target']+pfs)
   for label,row in zip(labels,array):writer.writerow([label]+row.tolist())
 plt.rcParams.update({'font.family':'Arial','font.size':11,'svg.fonttype':'none'})
 fig=plt.figure(figsize=(20,13));grid=fig.add_gridspec(2,2,width_ratios=[1.1,1.],height_ratios=[1.1,1.],hspace=.35,wspace=.25)
 ax=fig.add_subplot(grid[0,0]);ax.set_xlim(-.5,10.5);ax.set_ylim(-.5,10.);ax.axis('off')
 def box(x,y,text,color='#EEEEEE',w=1.25,h=.55,size=10):
  bounds[(x,y)]=(w/2+.08,h/2+.08)
  ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.08',facecolor=color,edgecolor='#333333',lw=1.,zorder=4));ax.text(x,y,text,ha='center',va='center',fontsize=size,zorder=5)
 def arrow(start,end,exc=True,color='#333333',rad=0.,lw=1.):
  pending.append((start,end,exc,color,rad,lw))
 bounds={};pending=[]
 rgpos={'RG_E_r':(1.4,8.9),'InE_r':(3.7,8.9),'InF_r':(6.3,8.9),'RG_F_r':(8.6,8.9)}
 for name,pos in rgpos.items():box(*pos,name[:-2],color='#DDEEFF'if name.startswith('RG')else'#FFF0BB')
 for e in edges:
  if e['source']in rgpos and e['target']in rgpos:arrow(rgpos[e['source']],rgpos[e['target']],e['sign']=='exc',rad=.2 if e['sign']=='inh'else 0.)
 pfpos={pf:(float(x),6.7)for pf,x in zip(pfs,np.linspace(.5,9.5,len(pfs)))}
 for i,(pf,pos)in enumerate(pfpos.items()):
  box(*pos,pf[:-2],np.array(matplotlib.colors.to_rgb(COLORS[i%len(COLORS)]))*.3+.7,w=1.25,size=10)
  sources=[e for e in edges if e['target']==pf and e['source']in rgpos]
  for e in sources:arrow((rgpos[e['source']][0],8.5),(pos[0],7.05),e['sign']=='exc',lw=1.4)
 relays={n:(3.4 if '_E_'in n else 6.6,5.45)for n in names if n in('PF_IN_E_r','PF_IN_F_r')}
 for name,pos in relays.items():box(*pos,name[:-2],color='#FFF0BB',w=1.5,size=9)
 for e in edges:
  if e['source']in pfpos and e['target']in relays:arrow((pfpos[e['source']][0],6.35),relays[e['target']],True,color='#888888',lw=.7)
  if e['source']in relays and e['target']in pfpos:arrow(relays[e['source']],(pfpos[e['target']][0],6.35),False,color='#888888',rad=.12,lw=.7)
 box(5.,3.7,'MN pools: distinct projections listed in the matrix below',w=9.4,h=.65,size=11)
 for i,pf in enumerate(pfs):
  x=pfpos[pf][0];n=int(np.count_nonzero(matrix[:,i]))
  if n:arrow((x,6.3),(x,4.1),True,color=COLORS[i%len(COLORS)],rad=.04,lw=1.5)
  ax.text(x,3.05,f'{n} muscles',ha='center',fontsize=9,color=COLORS[i%len(COLORS)])
 ax.text(5,2.6,'Shared MN box is a container, not all-to-all wiring.',ha='center',fontsize=10)
 for x,aff in((1.4,'Ia'),(5.,'II'),(8.6,'Ib')):
  targets={e['target']for e in edges if e['source'].startswith(aff+'_')and e['target']in pfs}
  count=sum(1 for e in edges if e['source'].startswith(aff+'_')and e['target']in pfs)
  box(x,1.15,f'{aff}: {count} PF edges\n{len(targets)} target PF cells',w=2.7,h=.8,size=10)
 ax.text(5,.25,'Feedback: Ia/II → muscle-carrying PF channels; Ib → selected stance E channels.\nExact afferent-source / PF targets and signs are listed in the feedback panel.',ha='center',va='center',fontsize=10)
 ax.set_title('(A) Fitted right-side RG → PF → MN architecture',loc='left',fontsize=13)
 architecture_routes=render_routes(ax,bounds,pending)
 heat=fig.add_subplot(grid[1,0]);im=heat.imshow(aggregated,cmap='Blues',aspect='auto');heat.set_xticks(range(len(pfs)),[n[:-2]for n in pfs],rotation=25,ha='right');heat.set_yticks(range(len(groups)),groups)
 for row in range(len(groups)):
  for col in range(len(pfs)):
   if aggregated[row,col]>0:heat.text(col,row,f'{aggregated[row,col]:.2g}',ha='center',va='center',fontsize=9,color='white'if aggregated[row,col]>.6*aggregated.max()else'black')
 heat.set_title('(B) PF → MN: summed fitted conductance by primary muscle group\nµS; conductance is not muscle force or torque',loc='left',fontsize=12);fig.colorbar(im,ax=heat,shrink=.65,label='Summed gmax (µS)')
 trace=np.load(a.fit/'BEST_TRACE.npz');t=trace['t'];keep=t>=t[-1]-3*cm['configuration']['target_period_s'];x=t[keep]-t[keep][0]
 sub=grid[:,1].subgridspec(4,1,hspace=.24);axes=[fig.add_subplot(sub[k])for k in range(4)]
 for name,style in(('RG_E_r','-'),('RG_F_r','--')):axes[0].plot(x,trace['rg'][keep,trace['rg_names'].tolist().index(name)],style,label=name)
 for i,pf in enumerate(pfs):axes[1].plot(x,trace['pf'][keep,trace['pf_names'].tolist().index(pf)],color=COLORS[i%len(COLORS)],label=pf[:-2])
 for name,style in(('iliacus_r','-'),('vas_int_r','--'),('bifemsh_r','-.'),('tib_ant_r',':')):
  i=trace['acts'].tolist().index(name);axes[2].plot(x,trace['mn'][keep,i],style,label=name);axes[3].plot(x,trace['desired'][keep,i],style,label=name)
 for axis,label in zip(axes,('RG (mV)','PF (mV)','MN / 5 mV','Target activation')):axis.set_ylabel(label);axis.grid(alpha=.2);axis.legend(fontsize=9,ncol=2);axis.set_xlim(x[0],x[-1])
 axes[0].set_title('(C) Saved D replay: three reference cycles',loc='left',fontsize=13);axes[-1].set_xlabel('Time in displayed window (s)')
 fig.suptitle(f'{family}: fitted circuit, muscle allocation, and neural activity\nOffline expected-afferent replay; not a physical walking run. Local reflex / commissural details are expanded separately.',fontsize=15)
 fig.subplots_adjust(top=.91,bottom=.07,left=.09,right=.97)
 for ext in('png','pdf','svg'):fig.savefig(a.output/(family+'_D_NEURAL_SCHEMATIC.'+ext),dpi=180)
 plt.close(fig)
 # All-muscle allocation: no threshold or top-N pruning hides weak routes.
 fig,ax=plt.subplots(figsize=(10,max(9,.25*len(acts))),constrained_layout=True);im=ax.imshow(matrix,cmap='Blues',aspect='auto');ax.set_xticks(range(len(pfs)),[n[:-2]for n in pfs],rotation=25,ha='right');ax.set_yticks(range(len(acts)),[act[:-2]+'  ['+', '.join(classify(act).groups)+']'for act in acts],fontsize=9)
 for row in range(len(acts)):
  for col in range(len(pfs)):
   if matrix[row,col]>0:ax.text(col,row,f'{matrix[row,col]:.2g}',ha='center',va='center',fontsize=7,color='white'if matrix[row,col]>.6*matrix.max()else'black')
 ax.set_title(f'{family}: exact fitted right PF → muscle-MN conductances\nBlank = no positive projection; all {len(acts)} muscle pools retained, including zero-drive rows.');fig.colorbar(im,ax=ax,label='gmax (µS)')
 for ext in('png','pdf','svg'):fig.savefig(a.output/(family+'_PF_MUSCLE_MAP.'+ext),dpi=180)
 plt.close(fig)
 # Exact afferent->PF routing, displayed separately from motor allocation.
 sensory=np.zeros((3,len(acts),len(pfs)))
 for e in edges:
  for k,aff in enumerate(('Ia','II','Ib')):
   if e['source'].startswith(aff+'_')and e['source'][len(aff)+1:]in acts and e['target']in pfs:
    sensory[k,acts.index(e['source'][len(aff)+1:]),pfs.index(e['target'])]=e.get('reference_budget_gmax_uS',e['snapshot_gmax_uS'])*e.get('shared_gain',1.)
 fig,axes=plt.subplots(1,3,figsize=(16,max(10,.25*len(acts))),sharey=True,constrained_layout=True)
 for k,aff in enumerate(('Ia','II','Ib')):
  im=axes[k].imshow(sensory[k]>0,cmap='Blues',vmin=0,vmax=1,aspect='auto');axes[k].set_title(aff+' → PF (positive entries only)');axes[k].set_xticks(range(len(pfs)),[pf[:-2]for pf in pfs],rotation=30,ha='right');axes[k].set_yticks(range(len(acts)),[act[:-2]for act in acts],fontsize=9)
 fig.suptitle(f'{family}: exact sensory-muscle origin → PF channel map\nFilled = configured excitatory pathway; Ib is load-gated. Blank = absent, not an assumed zero current.')
 for ext in('png','pdf','svg'):fig.savefig(a.output/(family+'_PF_FEEDBACK_MAP.'+ext),dpi=180)
 plt.close(fig)
 with(a.output/'AFFERENT_TO_PF.csv').open('w',newline='',encoding='utf-8')as f:
  cols=['source','target','sign','snapshot_gmax_uS','reference_budget_gmax_uS','shared_gain','stance_gated','reversal_mV','threshold_mV','saturation_mV'];w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader()
  for e in edges:
   if e['source'].startswith(('Ia_','II_','Ib_'))and e['target']in pfs:w.writerow(e)
 # A representative hip-flexor subnetwork shows actual local Ib topology.
 act='iliacus_r';mn='MN_'+act;parents=[pf for pf in pfs if any(e['source']==pf and e['target']==mn for e in edges)]
 fig,ax=plt.subplots(figsize=(12,10));ax.set_xlim(-.5,10.5);ax.set_ylim(-.8,10.2);ax.axis('off');bounds={};pending=[]
 positions={pf:(float(x),9.)for pf,x in zip(parents,np.linspace(1.,9.,len(parents)))};positions[mn]=(5.,6.8)
 for aff,x in(('Ia',1.),('II',5.),('Ib',9.)):positions[aff+'_'+act]=(x,2.5)
 raw=set(positions);relays_local=[]
 for e in edges:
  if e['source']in raw and e['target'].startswith(('IIX_','IBIN_','IBEXC_'))and any(x['source']==e['target']and x['target']==mn for x in edges):relays_local.append(e['target'])
 for relay in set(relays_local):positions[relay]=(9. if relay.startswith('IB')else 5.,4.7)
 for name,pos in positions.items():
  color='#DDEEFF'if name.startswith(('PF_','MN_'))else'#FFF0BB'if name in relays_local else'#E2F0DF';box(*pos,name[:-2].replace('_iliacus',' iliacus'),color,w=1.65,h=.6,size=10)
 for e in edges:
  if e['source']in positions and e['target']in positions:
   is_central=e['target']in parents;arrow(positions[e['source']],positions[e['target']],e['sign']=='exc',color='#228833'if is_central else'#333333',rad=(-.17 if e['source'].startswith('Ia_')else .14)if is_central else 0.,lw=1.2)
 box(5.,.65,'Expected iliacus signals: velocity (Ia), length (II), force (Ib)\nReplayed inputs in this D test; no mechanical feedback loop is integrated.',w=10.,h=1.,size=11)
 for aff in('Ia','II','Ib'):arrow((positions[aff+'_'+act][0],1.2),(positions[aff+'_'+act][0],2.15),True,color='#777777')
 has_central=any(e['source']=='Ib_'+act and e['target'].startswith(('RG_','PF_'))for e in edges)
 ax.text(5.,9.9,'PF → MN projections and Ia/II → PF routes are extracted from this fitted graph.',ha='center',fontsize=10)
 ax.text(5.,-.55,('Iliacus Ib → RG/PF present.'if has_central else'No iliacus Ib → RG or PF pathway in this checkpoint.')+'\nTriangle = excitation; filled circle = inhibition.',ha='center',fontsize=10)
 ax.set_title(f'{family}: iliacus flexor feedback subnetwork\nAntagonist and commissural paths omitted here; this panel shows homonymous routes and carrying PF channels.',fontsize=13)
 flexor_routes=render_routes(ax,bounds,pending)
 for ext in('png','pdf','svg'):fig.savefig(a.output/(family+'_FLEXOR_IB_SUBNETWORK.'+ext),dpi=180,bbox_inches='tight')
 plt.close(fig)
 (a.output/'ROUTE_VALIDATION.json').write_text(json.dumps(dict(architecture=architecture_routes,flexor_subnetwork=flexor_routes,total_neuron_intersections=0,policy='All synapses route orthogonally outside every neuron box; crossings between synapses do not denote junctions'),indent=2))
 print(json.dumps(dict(family=family,PF_cells=pfs,flexor_Ib=[row for row in feedback if row['muscle']in('iliacus_r','bifemsh_r','tib_ant_r')],ledger=str(a.ledger_output))),flush=True)
if __name__=='__main__':main()
