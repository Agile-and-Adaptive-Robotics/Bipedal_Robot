"""Editable, source-audited schematics of the preserved live 704/80 W2L run.

No simulation, model/editor mutation, or browser storage access. Population
readouts are explicit computational interfaces, not invented neuron synapses.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_key]='1'
import json,csv,hashlib,inspect,collections,textwrap,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from import_virtual_walker import mujoco
import native_kernel as K
from spiking_central_runtime import SpikingCentralController
from recruitment import apply_recruitment
from physics_variants import pool_inherited_hip_law

REPO=Path(r'D:\GitHub\Bipedal_Robot')
BASE=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\W2L_VirtualWalker_Newer\pace_0p2_supported10_v1')
OUT=REPO/'Neuromechanical_Models/Mujoco_SNS_models/W2L_704_80_baseline'
HERE=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(path,obj):path.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')

def extract():
    record=json.loads((BASE/'SPIKING_COMPARISON.json').read_text());cfg=record['source_configuration'];source=Path(cfg['source']);assert digest(source)==record['source_sha256']
    bank=json.loads((BASE.parent/'IMPORT_AUDIT.json').read_text())['muscles'];SpikingCentralController.configuration=record['configuration'];SpikingCentralController.sensory_spiking=True;SpikingCentralController.balance_configuration=None;SpikingCentralController.clock_scale=1.;SpikingCentralController.audit_state=False
    allowed=set(inspect.signature(K.NativeController).parameters)-{'source','bank'};c=SpikingCentralController(source,bank,**{k:v for k,v in cfg.items() if k in allowed})
    m=mujoco.MjModel.from_xml_path(str(BASE/'model.xml'));d=mujoco.MjData(m);d.qpos[:]=cfg['initial_qpos'];mujoco.mj_forward(m,d);p=K.bank_parameters(bank);law=pool_inherited_hip_law(c,bank,p,m,d);assert np.array_equal(p,np.load(BASE/'POOLED_FORCE_PARAMETERS.npz')['pooled'])
    recruitment=apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale']);c.apply_output_gain()
    root=ET.parse(source).getroot();module=next(x for x in root.iter('NeuralModule') if x.findtext('ModuleName')=='IntegrateFireSim');native=list(module.iter('Neuron'));types={t.findtext('ID'):t for t in module.findall('Synapses/NonSpikingSynapses/SynapseType')};connections=list(module.iter('Connexion'));assert len(c.v)==108 and len(c.edges)==len(connections)==216
    central=set(c.encoded.tolist());sensory=set(c.sensors.indices.tolist());nodes=[];cells=[];display_names=list(c.names);aliases=[]
    for i,name in enumerate(c.names):
        if c.names.count(name)==1:continue
        neighbours=[int(e[1]) if int(e[0])==i else int(e[0]) for e in c.edges if int(e[0])==i or int(e[1])==i];sides=collections.Counter(c.names[j][0] for j in neighbours if c.names[j].startswith(('L ','R ')));side=sides.most_common(1)[0][0]
        if side!=name[0]:display_names[i]=side+name[1:];aliases.append(dict(index=i,source_label=name,display_role=display_names[i],reason='Source duplicate L-labelled Renshaw cell has right-side actual neighbours. Display alias only; source and numerical indices unchanged.'))
    assert len(set(display_names))==108
    for i,name in enumerate(c.names):
        kind='CaLIF' if i in central else 'LIF_sensor' if i in sensory else 'graded';count=16 if i in central else 32 if i in sensory else 1;guid=native[i].findtext('ID');node=dict(index=i,id=guid,name=display_names[i],source_name=name,kind=kind,population=count,source_parameters=c.params[i].tolist(),tonic_current_nA=float(c.tonic[i]));nodes.append(node)
        for j in range(count):cells.append(dict(id=guid+':'+str(j),logical_channel=i,kind=kind))
    edges=[];variants={}
    for k,(e,original) in enumerate(zip(c.edges,connections)):
        pre,post=int(e[0]),int(e[1]);type_id=original.findtext('SynapseTypeID');type_name=types[type_id].findtext('Name');signature=(type_name,*map(float,e[2:]));variants.setdefault(signature,'T'+str(len(variants)+1).zfill(2));eid=original.findtext('ID') or 'source-row-'+str(k)
        row=dict(index=k,id=eid,source_index=pre,target_index=post,source=display_names[pre],target=display_names[post],source_label=c.names[pre],target_label=c.names[post],source_type_id=type_id,source_type_name=type_name,type_key=variants[signature],gmax_uS=float(e[2]),threshold_mV=float(e[3]),saturation_mV=float(e[4]),reversal_mV=float(e[5]),source_target_before_runtime_edit=c.names[c.ids[original.findtext('TargetID')]],sign_category='exc' if float(e[5])>-60 else 'inh');edges.append(row)
    ty=[dict(key=key,source_name=sig[0],gmax_uS=sig[1],threshold_mV=sig[2],saturation_mV=sig[3],reversal_mV=sig[4]) for sig,key in variants.items()]
    inputs=[]
    for target,b,kind,coeff,body in c.inputs:
        inputs.append(dict(target=c.names[target],target_index=int(target),source_kind=kind,physical_route=bank[b]['name'] if b>=0 else body,route_index=int(b),polynomial_ABCD_A_per_raw=coeff.tolist(),runtime_current_conversion_nA_per_A=1e9,raw_units='N' if kind=='Tension' else 'contact count (geometric0/1)' if kind in ('ContactCount','Contacts') else 'source receptor output: gain*series/parallel length or tension; not assumed Hz',output_units='nA'))
    outputs=[]
    for pre,b,coeff in c.outputs:
        if bank[b]['type']!='LinearHillMuscle':continue
        joint=('hip' if 'Hip' in c.names[pre] else 'knee' if 'Knee' in c.names[pre] else 'ankle')+'_'+c.names[pre][0];jid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,joint);tid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+bank[b]['name']);arm=-float(d.ten_J.reshape(m.ntendon,m.nv)[tid,m.jnt_dofadr[jid]]);outputs.append(dict(source=c.names[pre],source_index=int(pre),physical_route=bank[b]['name'],route_id=bank[b]['id'],polynomial_ABCD_V_per_V=coeff.tolist(),neuron_input_units='V =mV*0.001',actuator_input_units='V',home_moment_arm_m=arm,positive_coordinate='hip flexion, knee flexion, ankle dorsiflexion',source_Hill_parameters=p[b].tolist()))
    counts=dict(logical_channels=len(nodes),chemical_pathways=len(edges),RG_populations=sum(' RG ' in n['name'] and n['kind']=='CaLIF' for n in nodes),PF_populations=sum(' PF ' in n['name'] and n['kind']=='CaLIF' for n in nodes),sensory_populations=len(sensory),central_spiking_cells=16*len(central),sensory_spiking_cells=32*len(sensory),graded_cells=len(nodes)-len(central)-len(sensory),actual_membranes=len(cells),population_readouts=len(central)+len(sensory))
    assert counts==dict(logical_channels=108,chemical_pathways=216,RG_populations=4,PF_populations=8,sensory_populations=16,central_spiking_cells=192,sensory_spiking_cells=512,graded_cells=80,actual_membranes=784,population_readouts=28)
    ledger=dict(machine='easteregg2',run=str(BASE),source=str(source),source_sha256=digest(source),configuration=record['configuration'],source_configuration=cfg,census=counts,nodes=nodes,physical_membranes=cells,chemical_edges=edges,effective_type_variants=ty,physical_inputs=inputs,motor_outputs=outputs,compiled_edge_columns=['source_index','target_index','gmax_uS','threshold_mV','saturation_mV','reversal_mV'],compiled_edge_matrix=c.edges.tolist(),compiled_parameter_columns=['rest_mV','tau_ms','source_relative_size','GmaxCa_uS','Ca_activation_midpoint_mV','Ca_activation_slope','Ca_activation_tau_ms','Ca_inactivation_midpoint_mV','Ca_inactivation_slope','Ca_inactivation_tau_ms','Ca_reversal_mV'],recruitment_edits=recruitment,force_law_pooling=law,central_output_edits=c.output_edits,ankle_runtime_target_edits=c.ankle_phase_edits,scope='ACTUAL LIVE704spiking/80graded controller, not motor-wave replay.216logical chemical pathways use population-mean filtered spike readouts; source and target populations broadcast shared activation. Do not invent784x784independent synapses or divide conductance without preserving the shared saturating decoder.12PF/RG calcium-LIF populations and16sensory LIF populations replace their graded channel outputs. Population proxies are not soma voltages. Named type variants retain exact effective G/Erev/threshold/saturation after runtime calibration. Source sign labels are shorthand; actual current=g(Erev−Vpost).')
    ledger['explicit_display_aliases']=aliases
    ledger['time_steps_s']=dict(neural=float(c.dt),physics=float(record['result']['physics_dt']))
    ledger['population_readout_contract']=dict(central=dict(cells=16,reference_rate_Hz=float(c.settings[4]),filter_tau_s=float(c.settings[3]),refractory_s=.001,proxy='Vrest+20*filtered_mean_spike_activity; NOT soma voltage'),sensory=dict(cells=32,maximum_target_rate_Hz=200.,membrane_tau_s=.005,refractory_s=.001,filter_tau_s_by_name={c.names[i]:float(tau) for i,tau in zip(c.sensors.indices,c.sensors.tau)},target='max(0,I_external_nA/(20*sourceSize))',proxy='Vrest+20*filtered_mean_spike_activity'),shared_gate='Per logical chemical path: G*clip((proxy−threshold)/(saturation−threshold),0,1), then broadcast to each target membrane; global saturation prevents arbitrary replacement by independently clipped all-to-all synapses.')
    ledger['source_Hill_parameter_columns']=['l0_m','width_m','Kse_N_per_m','Kpe_N_per_m','damping_Ns_per_m','maximum_tension_N','sigmoid_A_V','sigmoid_B_N','sigmoid_C_per_V','sigmoid_D_N','parallel_rest_percent','minimum_parallel_percent','Ia_discharge_gain','II_discharge_gain','Ib_discharge_gain']
    assert record['spiking']['actual_spiking_cells']==704 and record['spiking']['remaining_graded_cells']==80
    assert [c.names[i] for i in c.encoded]==record['spiking']['spiking_names']
    assert [c.names[i] for i in c.sensors.indices]==record['spiking']['sensory_spiking_names']
    for actual,saved in zip(c.output_edits,record['spiking']['output_gain_edits']):assert actual==saved,(actual,saved)
    return ledger,record,c

def node_type(name):
    if ' RG ' in name and not name.endswith(' IN'):return 'HC-RG-E' if name.endswith(' ext') else 'HC-RG-F'
    if ' PF ' in name and not name.endswith(' IN'):return 'HC-PF-E' if name.endswith(' ext') else 'HC-PF-F'
    if name.endswith(' RE'):return 'RC'
    if ' MN ' in name:return 'MN'
    if name.endswith(' Ia IN'):return 'IN-IaIN'
    if ' PF ' in name:return 'IN-PF'
    if ' RG ' in name:return 'IN-InE' if 'ext' in name else 'IN-InF'
    if 'V3' in name:return 'IN-V3'
    if 'c1' in name.lower():return 'IN-C'
    if 'heel contact' in name:return 'SN-heel'
    if 'toe contact' in name:return 'SN-toe'
    return 'SN-II' if 'II' in name else 'SN-Ib' if 'Ib' in name else 'SN-Ia'

def group(name):
    side=name[0] if name.startswith(('L ','R ')) else ''
    if 'V3' in name or 'c1' in name.lower():return 'Commissural'
    if ' RG ' in name:return side+' RG'
    if ' PF ' in name:return side+(' Hip PF' if 'Hip' in name else ' Knee-ankle PF')
    if 'heel contact' in name or 'toe contact' in name:return side+' contact'
    return side+(' hip motor/reflex' if 'Hip' in name else ' knee motor/reflex' if 'Knee' in name else ' ankle motor/reflex')

def tag(e):return f"{e['type_key']} {e['source_type_name']} | G={e['gmax_uS']:.17g}uS; E={e['reversal_mV']:.17g}mV; threshold={e['threshold_mV']:.17g}mV; saturation={e['saturation_mV']:.17g}mV | source-row {e['index']}"

def layout(names):
    """Historical LH ordering: inhibitory cells above E/F, MN, then Renshaw."""
    positions={}
    for name in names:
        ext=(' ext' in name or 'extensor' in name or 'dorsiflexor' in name);column=0 if ext else 1
        if group(name)=='Commissural':positions[name]=(300+650*name.startswith('R '),180+300*('V3' in name));continue
        if ' RG ' in name or ' PF ' in name:
            positions[name]=(440+330*column,160 if name.endswith(' IN') else 330)
        elif name.endswith(' RE'):positions[name]=(440+330*column,760)
        elif ' MN ' in name:positions[name]=(440+330*column,550)
        elif name.endswith(' Ia IN'):positions[name]=(440+330*column,330)
        elif name.endswith(' II relay'):positions[name]=(180 if ext else 1030,140)
        elif 'Ib' in name:positions[name]=(180 if ext else 1030,550)
        elif ' Ia' in name:positions[name]=(180 if ext else 1030,330)
        elif 'II' in name:positions[name]=(180 if ext else 1030,760)
        else:positions[name]=(200+350*('toe' in name),240)
    return positions

def spec(ledger,selected,include_boundary=True):
    lookup={n['name']:n for n in ledger['nodes']};selected=set(selected);pos=layout(selected);nodes=[];labels={}
    for name in sorted(selected):
        n=lookup[name];label=name+(f" [{n['population']}spk]" if n['population']>1 else '')+(' [source L label]' if n['source_name']!=name else '');labels[name]=label;nd=dict(id=n['id'],type=node_type(name),label=label,x=pos[name][0],y=pos[name][1],grp='actual')
        if n['population']>1:
            P=n['population'];subnodes=[dict(type='PORT-load',label='Common population input (not a neuron)',x=80,y=350),dict(type='PORT-load',label='Mean spike filter (not a neuron)',x=1060,y=350)];subedges=[]
            for j in range(P):
                cell=name+' cell '+str(j+1).zfill(2);subnodes.append(dict(type=node_type(name),label=cell,x=300+170*(j%4),y=120+110*(j//4)));subedges.extend([dict(**{'from':'Common population input (not a neuron)','to':cell},sign='exc',gain=1.,tag='Common input broadcast; NOT an additional chemical synapse',pts=[]),dict(**{'from':cell,'to':'Mean spike filter (not a neuron)'},sign='exc',gain=1/P,tag=f'Population mean weight1/{P}; spike filter/normalization, NOTuS',pts=[])])
            nd['sub']=dict(nodes=subnodes,edges=subedges,_scope=f'Actual{P}dynamic spiking membranes plus2nonneural interface blocks; original population chemical pathways remain in parent graph. No invented cell-to-cell chemical matrix.')
        nodes.append(nd)
    edges=[];boundary={}
    for e in ledger['chemical_edges']:
        a,b=e['source'],e['target']
        if a in selected and b in selected:pass
        elif include_boundary and ((a in selected) != (b in selected)):
            other=b if a in selected else a;direction='OUT' if a in selected else 'IN';key=(direction,group(other))
            if key not in boundary:
                label=direction+' '+group(other)+' [boundary]';boundary[key]=label;nodes.append(dict(type='PORT-load',label=label,x=0 if direction=='IN' else 1300,y=80+140*sum(k[0]==direction for k in boundary),grp='boundary'))
            labels[other]=boundary[key]
        else:continue
        aa,bb=labels[a],labels[b];na=next(n for n in nodes if n['label']==aa);nb=next(n for n in nodes if n['label']==bb);lane=min(na['y'],nb['y'])-45-12*(len(edges)%3);pts=[dict(x=na['x'],y=lane),dict(x=nb['x'],y=lane)] if na['x']!=nb['x'] else []
        edges.append(dict(**{'from':aa,'to':bb},sign=e['sign_category'],gain=e['gmax_uS'],tag=tag(e),pts=pts))
    return dict(nodes=nodes,edges=edges,groups=dict(actual='Actual channels/populations',boundary='Boundary references; exact row IDs in tags'),_scope='Each edge tag preserves its actual individual source row and full four effective synaptic fields. Boundary duplicates are display references, not added neurons or weights.')

def editor_exports(ledger):
    groups=collections.defaultdict(list)
    for n in ledger['nodes']:groups[group(n['name'])].append(n['name'])
    models=[];overview=[]
    for k,(name,names) in enumerate(sorted(groups.items())):
        sub=spec(ledger,names);models.append(dict(name=name,data=sub));overview.append(dict(type='SUB',label=name,x=200+(k%4)*420,y=150+(k//4)*350,sub=sub))
    bundle=collections.Counter((group(e['source']),group(e['target']),e['sign_category']) for e in ledger['chemical_edges'] if group(e['source'])!=group(e['target']));ov_edges=[]
    for (a,b,sign),count in sorted(bundle.items()):ov_edges.append(dict(**{'from':a,'to':b},sign=sign,gain=count,tag=f"DISPLAY BUNDLE: {count} chemical pathways; gain field is COUNT, notuS. Exact216edges in ledger/full tab.",pts=[]))
    hierarchy=dict(nodes=overview,edges=ov_edges,_scope='Each SUB opens its actual induced circuit with exact boundary row references; interSUBedges count pathways, not summed conductance.')
    models.insert(0,dict(name='704spiking +80graded hierarchy',data=hierarchy));dump(OUT/'W2L_704_80_hierarchy_connectome.json',hierarchy);full=spec(ledger,[n['name'] for n in ledger['nodes']],False)
    for n in full['nodes']:
        source_name=n['label'].split(' [')[0];g=sorted(groups).index(group(source_name));n['x']+=1600*(g%4);n['y']+=1100*(g//4)
    bylabel={n['label']:n for n in full['nodes']}
    for e in full['edges']:
        a,b=bylabel[e['from']],bylabel[e['to']];mid=(a['x']+b['x'])/2;e['pts']=[dict(x=mid,y=a['y']),dict(x=mid,y=b['y'])]
    models.append(dict(name='Exact108channel216pathway graph',data=full));dump(OUT/'W2L_704_80_workspace.json',dict(models=models,active=0));dump(OUT/'W2L_704_80_exact_connectome.json',full)
    return groups,models

def setup_plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})
    return plt

def savefig(fig,name):
    for suffix in ('svg','pdf','png'):fig.savefig(OUT/(name+'.'+suffix),dpi=220,facecolor='white',bbox_inches='tight')

def orthogonal_routes(positions,edges):
    """Deterministic obstacle-aware Manhattan routing, no source graph changes."""
    import heapq
    scale=10;P={k:(round(x*scale),round(y*scale)) for k,(x,y) in positions.items()};blocked=set()
    for x,y in P.values():blocked.update((a,b) for a in range(x-4,x+5) for b in range(y-4,y+5))
    xmin=min(x for x,y in P.values())-12;xmax=max(x for x,y in P.values())+12;ymin=min(y for x,y in P.values())-12;ymax=max(y for x,y in P.values())+12;occupied=collections.Counter();result=[]
    for e in edges:
        a,b=P[e['source']],P[e['target']];dx=b[0]-a[0];dy=b[1]-a[1]
        if abs(dx)>abs(dy):direction=(1 if dx>0 else -1,0)
        else:direction=(0,1 if dy>0 else -1)
        if a[0]==b[0] and a[1]==b[1]:continue
        start=(a[0]+5*direction[0],a[1]+5*direction[1]);end=(b[0]-5*direction[0],b[1]-5*direction[1]);queue=[(0,0,start,None)];best={(start,None):0};previous={};terminal=None
        while queue:
            _,cost,p,last=heapq.heappop(queue)
            if cost!=best.get((p,last)):continue
            if p==end:terminal=(p,last);break
            for step in ((1,0),(-1,0),(0,1),(0,-1)):
                q=(p[0]+step[0],p[1]+step[1])
                if q in blocked or not(xmin<=q[0]<=xmax and ymin<=q[1]<=ymax):continue
                nc=cost+1+(2 if last is not None and last!=step else 0)+2*occupied[q];state=(q,step)
                if nc<best.get(state,float('inf')):best[state]=nc;previous[state]=(p,last);heapq.heappush(queue,(nc+abs(q[0]-end[0])+abs(q[1]-end[1]),nc,q,step))
        if terminal is None:raise RuntimeError('Orthogonal route failed: '+str((e['source'],e['target'])))
        points=[];state=terminal
        while state in previous:points.append(state[0]);state=previous[state]
        points.append(start);points.reverse()
        for p in points:occupied[p]+=1
        simplified=[points[0]]
        for j in range(1,len(points)-1):
            if (points[j][0]-points[j-1][0],points[j][1]-points[j-1][1])!=(points[j+1][0]-points[j][0],points[j+1][1]-points[j][1]):simplified.append(points[j])
        simplified.append(points[-1]);route=[(a[0]/scale+.38*direction[0],a[1]/scale+.38*direction[1])]+[(x/scale,y/scale) for x,y in simplified]+[(b[0]/scale-.38*direction[0],b[1]/scale-.38*direction[1])];result.append((e,route))
    return result

def circuit_figure(ledger,name,positions,labels,title,caption):
    plt=setup_plot();from matplotlib.patches import Circle,FancyArrowPatch
    lookup={n['name']:n for n in ledger['nodes']};edges=[e for e in ledger['chemical_edges'] if e['source'] in positions and e['target'] in positions];routes=orthogonal_routes(positions,edges);fig,ax=plt.subplots(figsize=(6.5,7));ax.set_aspect('equal');ax.axis('off');ax.set_title(title,loc='left',fontsize=11,pad=15);placed=[(4.,7.5,8.,8.5),(1.9,3.6,3.4,4.),(7.3,10.2,3.4,4.)] if 'left_RG_PF' in name else []
    for e,route in routes:
        line,=ax.plot(*np.array(route).T,color='#333333',lw=.85,zorder=1);line.set_gid('chemical_source_row_'+str(e['index']))
        if e['sign_category']=='inh':ax.add_patch(Circle(route[-1],.045,fc='#222222',ec='#222222',zorder=5))
        else:ax.add_patch(FancyArrowPatch(route[-2],route[-1],arrowstyle='-|>',mutation_scale=9,fc='white',ec='#222222',lw=.85,zorder=5))
        lengths=np.linalg.norm(np.diff(np.array(route),axis=0),axis=1);chosen=None
        for j in np.argsort(lengths)[::-1]:
            if lengths[j]<.4:continue
            for fraction in (.5,.27,.73,.12,.88):
                point=(1-fraction)*np.array(route[j])+fraction*np.array(route[j+1]);bounds=(point[0]-.34,point[0]+.34,point[1]-.19,point[1]+.19)
                if any(not(bounds[1]<old[0] or old[1]<bounds[0] or bounds[3]<old[2] or old[3]<bounds[2]) for old in placed):continue
                if any(abs(point[0]-x)<.7 and abs(point[1]-y)<.55 for x,y in positions.values()):continue
                chosen=(point,bounds);break
            if chosen is not None:break
        if chosen is not None:
            point,bounds=chosen;placed.append(bounds);ax.text(*point,e['type_key'],fontsize=10,ha='center',va='center',bbox=dict(fc='white',ec='none',pad=.5),zorder=4)
    for node,(x,y) in positions.items():
        n=lookup[node];fill='#fff3bf' if n['kind']=='LIF_sensor' else '#87cefa';patch=Circle((x,y),.38,fc=fill,ec='#222222',lw=1.1,zorder=6);patch.set_gid('channel_'+n['id']);ax.add_patch(patch)
        if n['population']>1:ax.add_patch(Circle((x,y),.43,fc='none',ec='#222222',lw=.65,zorder=6))
        ax.text(x,y,labels.get(node,node),ha='center',va='center',fontsize=10,zorder=7)
    if 'left_RG_PF' in name:
        ax.text(5.75,8.25,'Rhythm generation',ha='center',fontsize=11);ax.text(2.75,3.7,'Hip PF',ha='center',fontsize=11);ax.text(8.75,3.7,'Knee–ankle PF',ha='center',fontsize=11)
    xs,ys=np.array(list(positions.values())).T;ax.set_xlim(min(xs)-.8,max(xs)+.8);ax.set_ylim(min(ys)-.8,max(ys)+1.1);fig.subplots_adjust(bottom=.26,top=.91,left=.03,right=.98);fig.text(.04,.105,'\n'.join(textwrap.fill(line,85) for line in caption.splitlines()),fontsize=10,va='bottom');fig.text(.04,.025,textwrap.fill('Double ring: spiking population (cell count shown). Open triangle: excitatory type; filled circle: inhibitory. Current = g(Erev − Vpost).',85),fontsize=10,va='bottom')
    savefig(fig,name);plt.close(fig)
    # A clean functional-editor tab uses the exact same orthogonal routing.
    subset=spec(ledger,positions,False);labels_by_name={n['name']:n['name']+(f" [{n['population']}spk]" if n['population']>1 else '')+(' [source L label]' if n['source_name']!=n['name'] else '') for n in ledger['nodes']}
    for nd in subset['nodes']:
        original=nd['label'].split(' [')[0];x,y=positions[original];nd['x']=round(120*(x+2));nd['y']=round(120*(11-y))
    for ee,(e,route) in zip(subset['edges'],routes):
        assert ee['tag']==tag(e);ee['pts']=[dict(x=round(120*(x+2)),y=round(120*(11-y))) for x,y in route[1:-1]]
    dump(OUT/(name+'_connectome.json'),subset)
    return dict(name=title,data=subset)

def detail_figures(ledger):
    positions={};labels={}
    for motif,x,y in [('RG',4.5,6.),('Hip PF',1.5,1.5),('Knee PF',7.5,1.5)]:
        for action,xx in [('ext',x),('flx',x+2.5)]:
            for inhibitory,yy in [(False,y),(True,y+1.5)]:
                node='L '+motif+('' if motif=='RG' else '')+' '+action+(' IN' if inhibitory else '');positions[node]=(xx,yy);labels[node]=('IN '+('E' if action=='ext' else 'F')) if inhibitory else (('RG ' if motif=='RG' else 'PF ')+('E' if action=='ext' else 'F')+'\n16')
    core=circuit_figure(ledger,'W2L_704_80_left_RG_PF',positions,labels,'(B) Left rhythm and pattern-formation subnetworks','Actual induced circuit: 6 populations ×16 Ca/LIF cells, plus6 graded inhibitory interneurons.\nKnee PF also supplies ankle MNs; it is not a separate ankle PF. External afferent/commissural routes remain in the full216-pathway graph.')
    positions={};labels={}
    for action,x,side in [('ext',3.,0),('flx',7.,1)]:
        short='E' if action=='ext' else 'F';outer=.5 if side==0 else 9.5
        for tail,yy,label in [('PF '+action,9.,'PF '+short+'\n16'),('MN '+action,3.,'MN '+short),('MN '+action+' RE',1.,'RE '+short),(action+' Ia IN',5.,'Ia IN\n'+short)]:
            node='L Knee '+tail;positions[node]=(x,yy);labels[node]=label
        for tail,yy,label in [(action+' Ia',5.,'Ia '+short),(action+' Ib',3.,'Ib '+short),(('extensor' if side==0 else 'flexor')+' II relay',7.,'II '+short+'\n32')]:
            node='L Knee '+tail;positions[node]=(outer,yy);labels[node]=label
    motor=circuit_figure(ledger,'W2L_704_80_left_knee_pair',positions,labels,'(C) Left knee motor and reflex subnetwork','Actual induced circuit; MNs, Ia/Ib-labelled cells, Ia interneurons, and Renshaw cells remain graded.\nPhysical Ia adapters enter the Ia-IN cells in this source, rather than the separate Ia-labelled cells. Full physical input polynomials/units accompany the export.')
    return [core,motor]

def coarse_class(n):
    if n['kind']=='CaLIF':return 'RG' if ' RG ' in n['name'] else 'PF'
    if n['kind']=='LIF_sensor':return 'Spiking inputs'
    if ' MN ' in n['name'] and not n['name'].endswith(' RE'):return 'MN'
    if n['name'].endswith((' Ia',' Ib',' II')):return 'Graded afferents'
    return 'Graded INs'

def overview_figure(ledger):
    plt=setup_plot();from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
    fig,ax=plt.subplots(figsize=(6.5,5.3));ax.axis('off');ax.set_xlim(0,12);ax.set_ylim(0,7);ax.set_title('(A) W2L controller: 704 spiking and 80 graded membranes',loc='left',fontsize=11)
    classes=['Spiking inputs','Graded afferents','RG','PF','Graded INs','MN'];counts=collections.Counter()
    for n in ledger['nodes']:counts[coarse_class(n)]+=n['population']
    lookup={n['name']:coarse_class(n) for n in ledger['nodes']};matrix=np.zeros((6,6),int)
    for e in ledger['chemical_edges']:matrix[classes.index(lookup[e['source']]),classes.index(lookup[e['target']])]+=1
    assert matrix.sum()==216 and sum(counts.values())==784
    text=[('RG',(4.2,5.65),'RG: 4 ×16 Ca/LIF\n64 spiking membranes'),('PF',(4.2,3.8),'PF: 8 ×16 Ca/LIF\n128 spiking membranes'),('MN',(8.4,3.8),'12 graded motor neurons\n12 monoarticular routes'),('Spiking inputs',(1.,1.75),'12 II +4 heel/toe inputs\n16 ×32 =512 LIF cells'),('Graded afferents',(5.,1.75),'28 graded source\nIa / Ib / legacy II cells'),('Graded INs',(8.4,5.65),'40 graded cells\nlocal/commissural INs\nRenshaw cells')]
    for key,(x,y),label in text:
        height=1. if key=='Graded INs' else .8;patch=FancyBboxPatch((x-.05,y-height/2),3.05,height,boxstyle='round,pad=.06',fc='#fff3bf' if key=='Spiking inputs' else '#d9effb',ec='#333333',lw=1.);ax.add_patch(patch);ax.text(x+1.48,y,label,ha='center',va='center',fontsize=10)
    overview_paths={('RG','PF'):[(5.5,5.2),(5.5,4.25)],('PF','MN'):[(7.3,3.8),(8.3,3.8)],('RG','Graded INs'):[(7.3,5.85),(8.3,5.85)],('Graded INs','RG'):[(8.3,5.45),(7.3,5.45)],('PF','Graded INs'):[(7.3,4.),(7.75,4.),(7.75,5.65),(8.3,5.65)],('Graded INs','PF'):[(9.9,5.2),(9.9,4.55),(5.95,4.55),(5.95,4.25)],('Spiking inputs','MN'):[(2.5,2.2),(2.5,2.65),(10.2,2.65),(10.2,3.35)],('Graded afferents','MN'):[(6.5,2.2),(6.5,2.95),(9.4,2.95),(9.4,3.35)]}
    for (a,b),route in overview_paths.items():
        n=matrix[classes.index(a),classes.index(b)]
        if not n:continue
        ax.plot(*np.array(route).T,lw=.85,color='#555555',zorder=0);ax.add_patch(FancyArrowPatch(route[-2],route[-1],arrowstyle='->',mutation_scale=10,lw=.85,color='#555555',zorder=1))
    ax.text(.2,.65,'Inputs: length-related outputs, tendon tension (N), heel/toe contact (0/1).\nOutput: MN mV → V → source LinearHill muscle tension (N).\nActual GRF/support are recorded; no added GRF-to-balance pathway.',fontsize=10,va='center')
    fig.subplots_adjust(bottom=.14,left=.025,right=.98,top=.92);fig.text(.04,.025,textwrap.fill('Layer overview; selected arrows only. All 216 chemical pathways and 44 physical input adapters remain in the editable graph and ledger.',85),fontsize=10,va='bottom');savefig(fig,'W2L_704_80_overview');plt.close(fig)
    dump(OUT/'LAYER_CONNECTION_COUNTS.json',dict(order=classes,matrix_source_rows_target_columns=matrix.tolist(),chemical_pathways=int(matrix.sum()),physical_membrane_counts=dict(counts),scope='All216paths; overview arrows selectively illustrate hierarchy and must not be treated as the full graph.'))

def type_tables(ledger):
    plt=setup_plot();from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(OUT/'W2L_704_80_synapse_types.pdf') as pdf:
        for page,start in enumerate(range(0,len(ledger['effective_type_variants']),14)):
            rows=ledger['effective_type_variants'][start:start+14];fig,ax=plt.subplots(figsize=(6.5,9));ax.axis('off');ax.set_title('Effective synaptic types: full fields (page '+str(page+1)+')',loc='left',fontsize=11)
            cells=[[r['key'],textwrap.fill(r['source_name'],26),f"{r['gmax_uS']:.7g}",f"{r['reversal_mV']:g}",f"{r['threshold_mV']:g}",f"{r['saturation_mV']:g}"] for r in rows];table=ax.table(cellText=cells,colLabels=['Key','Source type name','Gmax\n(µS)','Erev\n(mV)','Pre-thr\n(mV)','Pre-sat\n(mV)'],colWidths=[.06,.36,.16,.12,.15,.15],loc='upper left',cellLoc='left');table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,2.7)
            fig.text(.04,.035,textwrap.fill('Full precision remains in JSON/CSV and editor edge tags. Runtime calibration can create several effective variants of one source-named type.',85),fontsize=10);pdf.savefig(fig,bbox_inches='tight');fig.savefig(OUT/('W2L_704_80_synapse_types_'+str(page+1)+'.svg'),bbox_inches='tight');plt.close(fig)

def four_cycle_figure(record):
    plt=setup_plot();from scipy.signal import find_peaks
    z=np.load(BASE/'trajectory.npz');central=np.load(BASE/'central_spikes.npz');t=z['t'];nt=central['t'];names=central['names'].tolist();decoded=central['decoded_activation'];trace=z['neuro'];nnames=z['neuron_names'].tolist();signal=decoded[:,names.index('L RG flx')];indices=np.flatnonzero(nt>=3.);peaks,_=find_peaks(signal[indices],distance=round(.5/np.median(np.diff(nt))),prominence=.25*np.ptp(signal[indices]));peaks=indices[peaks];assert len(peaks)>=5;start,end=map(float,nt[peaks[[0,4]]]);mask=(t>=start)&(t<=end);nm=(nt>=start)&(nt<=end);bw=record['result']['robot_mass_kg']*9.81
    fig,axes=plt.subplots(9,1,figsize=(6.5,10.5),sharex=True,gridspec_kw={'hspace':.16});blue='#0072B2';red='#D55E00';styles=[('L','ext',blue,'-'),('L','flx',red,'-'),('R','ext',blue,'--'),('R','flx',red,'--')]
    for row,prefix in enumerate(('RG','Hip PF','Hip MN')):
        ax=axes[row]
        for side,action,color,ls in styles:
            name=side+' '+prefix+' '+action;y=trace[mask,nnames.index(name)] if prefix=='Hip MN' else decoded[nm,names.index(name)];tx=t[mask] if prefix=='Hip MN' else nt[nm];ax.plot(tx,y,color=color,ls=ls,lw=1.1,label=side+' '+('E' if action=='ext' else 'F'))
        ax.set_ylabel(('RG activity\n(spike-derived)' if row==0 else 'Hip PF activity\n(spike-derived)' if row==1 else 'Hip MN\n(mV)'),fontsize=10)
    axes[0].legend(ncol=4,loc='upper right',bbox_to_anchor=(1.,1.25),frameon=False,fontsize=10)
    for ax,col,label in zip(axes[3:6],range(3),('Hip flexion\n(deg)','Knee flexion\n(deg)','Ankle DF\n(deg)')):
        for side,k,color,ls in [('L',col,blue,'-'),('R',col+3,red,'--')]:ax.plot(t[mask],z['angles'][mask,k],color=color,ls=ls,lw=1.1,label=side)
        ax.set_ylabel(label,fontsize=10)
    axes[3].legend(ncol=2,loc='upper right',frameon=False,fontsize=10)
    force=z['foot_force_world'][:,:,1]/bw
    axes[6].plot(t[mask],force[mask,0],color=blue,lw=1.,label='L foot');axes[6].plot(t[mask],force[mask,1],color=red,ls='--',lw=1.,label='R foot');axes[6].set_ylabel('Foot normal\nforce (BW)');axes[6].legend(ncol=2,loc='upper right',frameon=False)
    channels=np.column_stack((z['rail_force_world'][:,1],z['guide_force_world'][:,1],z['harness_force'][:,1]));up=np.maximum(channels,0).sum(axis=1)/bw;down=np.maximum(-channels,0).sum(axis=1)/bw
    axes[7].plot(t[mask],up[mask],color='#009E73',lw=1.,label='Upward');axes[7].plot(t[mask],down[mask],color='#CC79A7',ls='--',lw=1.,label='Downward');axes[7].set_ylabel('Support\nreaction (BW)');axes[7].legend(ncol=2,loc='upper right',frameon=False)
    clear=np.minimum(z['foot_clearance'][:,::2],z['foot_clearance'][:,1::2]);axes[8].plot(t[mask],clear[mask,0]*1000,color=blue,lw=1.,label='L');axes[8].plot(t[mask],clear[mask,1]*1000,color=red,ls='--',lw=1.,label='R');axes[8].axhline(0,color='#777777',lw=.7);axes[8].set_ylabel('Whole-foot\nclearance (mm)');axes[8].set_xlabel('Time (s)')
    for side,color in [(0,blue),(1,red)]:
        swung=mask&(clear[:,side]>.001)&(force[:,side]<.025);axes[8].fill_between(t,0,1,where=swung,transform=axes[8].get_xaxis_transform(),color=color,alpha=.07,lw=0)
    for ax in axes:
        ax.grid(axis='y',color='#dddddd',lw=.45);ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=10);ax.set_xlim(start,end)
    fig.suptitle('W2L spiking/graded controller: four left-RG cycles',fontsize=11,y=.995);fig.subplots_adjust(top=.95,bottom=.12,left=.16,right=.985);fig.text(.04,.025,textwrap.fill('RG/PF: filtered population spike activity; MNs: graded voltage. Lower shading: unloaded whole-foot clearance. Passive rails/guides; full-run 10/9 swings, 36.48% BW mean net support.',85),fontsize=10);savefig(fig,'W2L_704_80_four_cycles');plt.close(fig)
    report=dict(source=str(BASE),source_sha256=record['source_sha256'],window_s=[start,end],selection_rule='First5leftRG-F decoded-activity peaks after3s; four complete intervals, fixed minimum0.5speak spacing and prominence25%ofpost3ssignal range. No parameter/outcome selection.',RG_PF_signal='normalized filtered spike population activity',MN_signal='graded membrane voltage mV',joint_signal='actual saved anatomical angles: positive hip/knee flexion and ankle dorsiflexion',clearance='minimum of actual oriented main-foot/toe box clearance above moving floor',support='Positive/negative reactions summed separately over rails,guides,harness; not net cancellation',all_data_hashes={name:digest(BASE/name) for name in ('trajectory.npz','central_spikes.npz','sensory_spikes.npz','SPIKING_COMPARISON.json')},mean_upward_BW=float(up[mask].mean()),mean_downward_BW=float(down[mask].mean()),full_run_assessment=record['assessment'])
    dump(OUT/'FOUR_CYCLE_PROVENANCE.json',report);return report

def validate_exports(ledger,models):
    known={'SN-Ia','SN-II','SN-Ib','SN-heel','SN-toe','PORT-load','IN-V3','IN-C','HC-RG-E','HC-RG-F','IN-InE','IN-InF','HC-PF-E','HC-PF-F','IN-PF','IN-IaIN','RC','MN','MUSCLE','SUB'};checked=[]
    def check(data):
        labels=[n['label'] for n in data['nodes']];assert len(labels)==len(set(labels));assert all(n['type'] in known for n in data['nodes']);assert all(e['from'] in labels and e['to'] in labels for e in data['edges'])
        for n in data['nodes']:
            if 'sub' in n:check(n['sub'])
    for model in models:check(model['data']);checked.append(dict(model=model['name'],nodes=len(model['data']['nodes']),edges=len(model['data']['edges'])))
    full=models[-1]['data'];assert len(full['nodes'])==108 and len(full['edges'])==216
    for e,saved in zip(ledger['chemical_edges'],full['edges']):assert saved['gain']==e['gmax_uS'] and saved['tag']==tag(e)
    import subprocess
    editor=REPO/'Code/MuJoCo_SNS/spinal/connectome_block_editor.html';html=editor.read_text(encoding='utf-8');start=html.index('function specOf(');end=html.index('function serialize()',start);function=html[start:end]
    javascript=function+'''
const input=JSON.parse(require('fs').readFileSync(0,'utf8'));let contexts=0,edges=0;
function check(spec){
 const ids=new Map(spec.nodes.map((n,i)=>[n.label,'check_'+i]));
 const ns=spec.nodes.map(n=>({...n,id:ids.get(n.label)}));
 const es=spec.edges.map(e=>({...e,from:ids.get(e.from),to:ids.get(e.to)}));
 const actual=specOf(ns,es,spec.groups||{},new Set());
 if(actual.nodes.length!==spec.nodes.length||actual.edges.length!==spec.edges.length)throw Error('Lost graph content');
 for(let k=0;k<actual.edges.length;k++)for(const field of ['from','to','sign','gain','tag'])if(actual.edges[k][field]!==spec.edges[k][field])throw Error('Changed '+field);
 contexts++;edges+=actual.edges.length;for(const n of spec.nodes)if(n.sub)check(n.sub);
}
for(const m of input.models)check(m.data);process.stdout.write(JSON.stringify({contexts,edges,exact_edge_fields_preserved:true}));
'''
    result=subprocess.run([r'D:\NodeJS\node-v22.23.2-win-x64\node.exe','-e',javascript],input=json.dumps(dict(models=models)),text=True,capture_output=True,check=True);roundtrip=json.loads(result.stdout)
    return dict(census=ledger['census'],all_source_rows_exact=True,saved_output_gain_edits_exact=True,all_editor_types_supported=True,all_endpoints_valid=True,models=checked,actual_editor_specOf_roundtrip=roundtrip,editor_sha256=digest(editor),GUI_import_performed=False,source_or_browser_storage_modified=False,scope='Actual editor specOf JavaScript serialization executed in Node on every nested exported context, alongside source matrix/census checks. No DOM,browser,editor,or localStorage interaction. Full-field numeric edge tags/gains preserved; extended metadata authoritatively remains in ledger.')

def write_readme():
    text=f'''# W2L spiking/graded controller: editable schematic

Machine: **easteregg2**. Selected saved live run:
`{BASE}`

The controller has four RG and eight PF populations of 16 calcium/adapting LIF cells, and 16 sensory populations of 32 LIF cells (12 homonymous II and four heel/toe inputs). That is **704 spiking and 80 graded membranes**: 108 logical channels and 216 chemical pathways. Motor neurons remain graded. The 28 population spike filters/readouts are computational interfaces, not additional neurons or independently stored cellwise synapse matrices.

## Open the drawing in the block editor

Preserve your current tab first. Use **Templates → Blank** to create a new tab, then **Load JSON from file…** and select:

`{OUT / 'W2L_704_80_hierarchy_connectome.json'}`

Double-click a SUB node to inspect its circuit. Double-click a population to inspect its actual 16 or 32 cells and nonneural input/readout blocks; F4 returns to the parent. Individual schematic files also load directly:

- `{OUT / 'W2L_704_80_left_RG_PF_connectome.json'}`
- `{OUT / 'W2L_704_80_left_knee_pair_connectome.json'}`
- `{OUT / 'W2L_704_80_exact_connectome.json'}`

**W2L_704_80_workspace.json is an archival multi-tab snapshot, not a direct input to the editor's single-tab Load button.** No browser, editor, or localStorage state was changed to produce these files.

## Exact graph and editable figures

`W2L_704_80_graph_ledger.json` preserves the actual compiled matrix, 784 membrane memberships, 44 physical adapters, parameter units, 12 motor outputs and anatomical moment-arm signs, runtime ankle target swaps, gain calibration, and pooled hip force-law settings. Full original source names and IDs remain alongside any display alias. `W2L_704_80_chemical_edges.csv` and every editor edge tag preserve the four full-precision effective synaptic fields: conductance, reversal potential, presynaptic threshold, and saturation. The PDF/SVG synaptic-type table lists 49 effective variants with their original source-named types. Runtime calibration can give several effective variants of one named type.

The overview, left RG/PF circuit, knee motor/reflex pair, and four-cycle result are provided as editable text/vector SVG, PDF, and PNG. Source fonts are at least 10 pt at approximately thesis-column width. The source graph is unchanged; these are schematic exports, not executable controller replacements. The hierarchy uses explicit pathway counts, not summed conductance, for inter-subsystem bundles. The full graph is authoritative; publication detail panels are induced subnetworks with boundary paths retained in the hierarchy/full graph.

## Source-specific details

- One source Renshaw cell has a duplicated L label but actual right-side neighbours. Its display alias explicitly marks the right-side role; the original name/ID and numerical connections are preserved.
- RG/PF traces show filtered population spike activity, not individual soma voltage. Source calcium-bearing Renshaw cells remain graded.
- In this source Hill implementation, Ia is gain times series length, II is gain times parallel length, and Ib is gain times tension. Do not silently substitute a generic Ia=velocity equation. Physical Ia adapters enter Ia-IN cells in this model rather than its separate Ia-labelled cells.
- Heel/toe inputs are geometric contact levels, not GRF or heelstrike impulses. Actual GRF, support, and clearance are recorded diagnostics. This baseline adds no GRF-to-body-balance pathway.
- Source ankle MN ext actuates dorsiflexion; source ankle MN flx actuates plantarflexion. Four explicitly recorded runtime PF-target swaps are retained. Physical route names alone are not anatomical action labels.
- Excitatory/inhibitory arrow categories are source-type shorthand. Actual instantaneous chemical current is g(Erev−Vpost), so its sign depends on postsynaptic voltage.

## Saved-data result

`FOUR_CYCLE_PROVENANCE.json` identifies the fixed 3.842–7.488 s window: four intervals between the first five left-RG-F activity peaks after 3 s. Data are from the actual live neural controller, not the later offline motor-wave replay. Panels show RG/PF spike activity, graded hip MN voltage, anatomical hip/knee/ankle angles, foot normal forces, separately summed upward/downward support, and minimum oriented main-foot/toe clearance. Large impact peaks are retained. The whole 10 s run has 10/9 unloaded, cleared swings and 36.48% BW net support; passive rails/guides remain. This is supported stepping, not independent walking or a periodic-orbit result.

Regenerate on easteregg2:

```powershell
& 'D:\\Anaconda\\envs\\myo\\python.exe' '{Path(__file__).resolve()}'
```

`EXPORT_VERIFICATION.json` records exact census, saved gain-ledger agreement, valid native editor types/endpoints, and full-matrix export checks. Browser import was not performed.
'''
    (OUT/'README.md').write_text(text,encoding='utf-8')

def main():
    OUT.mkdir(parents=True,exist_ok=True);ledger,record,c=extract();dump(OUT/'W2L_704_80_graph_ledger.json',ledger);groups,models=editor_exports(ledger);models=detail_figures(ledger)+models;# full exact graph stays last
    dump(OUT/'W2L_704_80_workspace.json',dict(models=models,active=2));overview_figure(ledger);type_tables(ledger);window=four_cycle_figure(record);check=validate_exports(ledger,models);dump(OUT/'EXPORT_VERIFICATION.json',check)
    with (OUT/'W2L_704_80_chemical_edges.csv').open('w',encoding='utf-8-sig',newline='') as f:
        keys=['index','source','target','source_type_name','type_key','gmax_uS','reversal_mV','threshold_mV','saturation_mV','source_index','target_index'];writer=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore');writer.writeheader();writer.writerows(ledger['chemical_edges'])
    write_readme()
    dump(OUT/'SOURCE_MANIFEST.json',dict(machine='easteregg2',generator=str(Path(__file__)),generator_sha256=digest(__file__),source_sha256=ledger['source_sha256'],verified_run=str(BASE),files={p.name:digest(p) for p in OUT.iterdir() if p.is_file() and p.name!='SOURCE_MANIFEST.json'}));print(json.dumps(dict(output=str(OUT),census=ledger['census'],window=window['window_s'],verified=True),indent=2))

if __name__=='__main__':main()


