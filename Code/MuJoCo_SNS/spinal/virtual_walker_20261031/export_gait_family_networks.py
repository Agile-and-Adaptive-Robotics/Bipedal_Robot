"""Read-only standing-candidate graph extraction and editable scientific figures.

Executes the actual standing runner's nested build function, not the physical
simulation. Strict saved builder/parameter hashes and saved node counts gate it.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import sys,json,ast,inspect,hashlib,collections,copy,csv,subprocess,textwrap
from pathlib import Path
import numpy as np

REPO=Path(r'D:\GitHub\Bipedal_Robot')
HERE=Path(__file__).resolve().parent
FOLDER=REPO/'Code/MuJoCo_SNS/robot_torso/bpa_neural_coordination'
sys.path.insert(0,str(FOLDER))
import gait_standing_family as F
import export_baseline_network as DRAW
PRISTINE_G=copy.deepcopy(F.P.G)
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
OUT=REPO/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families'

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def plain(o):
    if isinstance(o,np.ndarray):return o.tolist()
    if isinstance(o,np.generic):return o.item()
    if isinstance(o,dict):return {str(k):plain(v) for k,v in o.items()}
    if isinstance(o,(list,tuple)):return [plain(v) for v in o]
    if isinstance(o,Path):return str(o)
    if hasattr(o,'params'):return dict(class_name=type(o).__name__,parameters=plain(o.params))
    return o
def dump(p,o):p.write_text(json.dumps(plain(o),indent=2,ensure_ascii=False),encoding='utf-8')

CASES={
 's3k':('other_family_pressure/s3k_bio_prepared_v1','rustuna_s3k_v1/BEST.json'),
 'w2lvar':('normalized_prepared_transfer/w2lvar_bio_graded_v1','group_diagnosis/NORMALIZED_STANDING_CHECKPOINT.json'),
 'syn6':('pinned_preparation/syn6_prepared_bio_nominal_v1','rustuna_syn6_v1/BEST.json')}

def family_class(name):
    if name.startswith(('RG_','ADAP_','InE_','InF_')):return 'RG'
    if name.startswith('PF_'):return 'PF'
    if name.startswith('MN_'):return 'MN'
    if name.startswith(('Ia_','II_','Ib_')):return 'Afferents'
    if name.startswith(('IaIN_','IIX_','IIIN_','IBIN_','IBEXC','RC_')):return 'Reflex INs'
    return 'Descending and commissural'

def extract(family):
    # w2lvar.build deliberately leaves its runner overlay in params.G. Never
    # allow it to contaminate the next family's physical-input metadata.
    F.P.G.clear();F.P.G.update(copy.deepcopy(PRISTINE_G))
    stem,cfgpath=CASES[family];prefix=ROOT/stem
    result=json.loads(prefix.with_suffix('.result.json').read_text())
    assert result['failure'] is None and result['family']==family
    cfg=json.loads((ROOT/cfgpath).read_text())
    settings=cfg['global_parameters'] if family=='w2lvar' else cfg['trial']['params']
    checks={}
    for name in ('runner.py','params.py','build_network.py','build_network_w2lvar.py','build_network_syn6.py'):
        path=F.SPINAL/name;checks[name]=dict(path=str(path),actual=digest(path),saved=result['source_sha256'][name])
        assert checks[name]['actual']==checks[name]['saved'],name+' changed: no winner substitution allowed'
    dependencies=['muscle_map.py']
    if family=='w2lvar':dependencies+=['joint_pf_weights.json']
    if family=='syn6':dependencies+=['synergy_basis.npz','fsa_results/fsa_backsolve.npz','fsa_backsolve.py']
    archive=json.loads((ROOT/'followup_source_snapshot/MANIFEST.json').read_text())
    archived={row['source']:row for row in archive['files']+archive.get('supplemental_direct_dependencies',[])}
    dependency_checks={}
    for name in dependencies:
        path=F.SPINAL/name;row=archived.get(str(path))
        if row is None:
            original=json.loads((ROOT/'source_snapshot/MANIFEST.json').read_text())
            row=next((x for x in original['files'] if x['source']==str(path)),None)
        if row is None:
            if name!='fsa_backsolve.py':raise ValueError('Missing dependency snapshot: '+str(path))
            # This analytical helper was not archived with the run. Record the
            # limitation explicitly; saved weight projections and actual neural
            # response replay below independently gate its mapped output.
            dependency_checks[name]=dict(path=str(path),sha256=digest(path),snapshot=None,scope='Unarchived analytical helper; exact saved bilateral projection and actual saved neural-current response are separately required.')
            continue
        actual=digest(path);assert actual==row['sha256'] and digest(row['snapshot'])==actual
        dependency_checks[name]=dict(path=str(path),sha256=actual,snapshot=row['snapshot'])
    for key,saved_key in [('pf_mn_gain','shared_pf_mn_gain'),('ib_group_gain','ib_group_gain'),('standing_bias_scale','standing_bias_scale')]:
        assert settings[key]==result[saved_key]
    model=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL))
    acts=[model.actuator(i).name for i in range(model.nu)]
    assert len(acts)==92
    os.environ['AARL_NET']=family
    # These are the only two pre-build modifications of --reflex-stand.
    F.P.G['ii_set_pose']=1.;F.P.G['posture_to_mn']=0.
    capture={}
    class Captured(BaseException):pass
    def snapshot(inner,names):
        capture['net']=inspect.currentframe().f_back.f_locals['n']
        capture['inner']=inner;capture['names']=names
        raise Captured()
    tree=ast.parse(Path(F.__file__).read_text())
    outer=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='run')
    node=next(x for x in outer.body if isinstance(x,ast.FunctionDef) and x.name=='build')
    module=ast.Module(body=[node],type_ignores=[]);ast.fix_missing_locations(module)
    env=dict(F.__dict__)
    env.update(oldbuild=F.B.build,capture=capture,gain=settings['pf_mn_gain'],ib_group_gain=settings['ib_group_gain'],
       pathway_group_scales=cfg.get('pathway_group_scales'),backend='sparse',compiled_transform=snapshot,
       fast_motor_io=True,pinned_actuator_settling=True)
    exec(compile(module,str(Path(F.__file__))+'::actual_nested_build','exec'),env)
    try:env['build'](acts,F.P.DT,True)
    except Captured:pass
    else:raise RuntimeError('Did not intercept exact candidate build')
    n=capture['net'];c=capture['inner'].c;names=capture['names']
    assert len(names)==result['neurons'] and not c.spiking
    saved_projection=prefix.parent/(prefix.name+'_projection.json')
    assert json.loads(saved_projection.read_text())==capture['projection'],'Saved bilateral weight projection differs'
    with np.load(prefix.with_suffix('.release_state.npz')) as release:
        assert release['neural_V'].shape==c.V.shape
    source_connections=collections.defaultdict(list)
    for con in n.net.connections:
        pre=n.net.populations[con['source']]['name'];post=n.net.populations[con['destination']]['name']
        source_connections[(pre,post)].append(dict(name=con.get('name'),type=con.get('type'),parameters=plain(con['params'])))
    nodes=[]
    for i,name in enumerate(names):
        pop=next(p for p in n.net.populations if p['name']==name)
        nodes.append(dict(index=i,id=family+':'+str(i),name=name,source_name=name,kind='graded',population=1,
             layer=family_class(name),source_neuron_parameters=plain(pop.get('params',getattr(pop.get('type'),'params',{})))))
    variants={};edges=[]
    for post,pre in zip(*np.nonzero(c.g_max_non)):
        g,er,lo,hi=map(float,(c.g_max_non[post,pre],c.del_e[post,pre],c.e_lo[post,pre],c.e_hi[post,pre]))
        signature=(g,er,lo,hi);variants.setdefault(signature,'T'+str(len(variants)+1).zfill(3));key=variants[signature]
        src=source_connections.get((names[pre],names[post]),[])
        edge=dict(index=len(edges),id=family+':edge:'+str(len(edges)),source_index=int(pre),target_index=int(post),
          source=names[pre],target=names[post],gmax_uS=g,reversal_mV=er,threshold_mV=lo,saturation_mV=hi,
          type_key=key,source_type_name=' / '.join(sorted(set(x.get('name') or 'unnamed SNS connection' for x in src))) or 'bilateral pooled counterpart',
          source_connection_definitions=src,sign_category='exc' if er>0 else 'inh')
        edges.append(edge)
    disabled=[]
    index={name:i for i,name in enumerate(names)}
    for (pre,post),defs in source_connections.items():
        if c.g_max_non[index[post],index[pre]]==0:
            j,i=index[post],index[pre]
            disabled.append(dict(source=pre,target=post,gmax_uS=0.,reversal_mV=float(c.del_e[j,i]),threshold_mV=float(c.e_lo[j,i]),saturation_mV=float(c.e_hi[j,i]),source_connection_definitions=defs))
    input_edges=[]
    for post,port in zip(*np.nonzero(c.input_connectivity)):
        input_edges.append(dict(input=n.inputs[port],target=names[post],weight=float(c.input_connectivity[post,port]),units='nA'))
    flags={k:plain(getattr(n,k)) for k in ('interleg','joint_pf','renshaw','vest','sides') if hasattr(n,k)}
    policy=dict(model=str(F.R.MODEL),model_sha256=digest(F.R.MODEL),actuator_entries=92,motor_neurons=sum(x['layer']=='MN' for x in nodes),
      biologically_force_pruned=sorted(F.R.PRUNE_MUSCLES),biologically_active_routes=86,
      policy='All92 MNs and afferent arcs remain; runner zeros force parameters of six listed routes. This export is biological, not the86-route virtual BPA replacement.',
      motor_conversion='clip(MN voltage / E_HI,0,1); E_HI='+str(F.R.E_HI)+' mV; native activation and contraction dynamics follow.',
      input_policy='--reflex-stand: II setpoint at initial standing length; solved per-muscle POST_ bias with saved scale; supraspinal BAL surrogate inputs zero; posture_to_mn=0.',
      apparatus=result['apparatus'],frame=result['frame'])
    ledger=dict(machine='easteregg2',family=family,scope='Exact reconstructed graded biological standing candidate; not a gait/curriculum or spiking result.',
       saved_result=str(prefix.with_suffix('.result.json')),saved_result_sha256=digest(prefix.with_suffix('.result.json')),
       saved_configuration=str(ROOT/cfgpath),saved_configuration_sha256=digest(ROOT/cfgpath),settings=settings,
       pathway_group_scales=cfg.get('pathway_group_scales'),source_checks=checks,dependency_snapshot_checks=dependency_checks,
       saved_bilateral_projection_exact=True,saved_bilateral_projection_file=str(saved_projection),
       construction='Actual AST-extracted F.run.build; actual two --reflex-stand prebuild flags. Captured after tying, candidate scaling, and64-step exact sparse parity, before any physics.',
       current_standing_runner_sha256=digest(F.__file__),saved_standing_runner_sha256=result['source_sha256']['gait_standing_family.py'],
       standing_runner_scope='Runner file has acquired subsequent optional experiments; selected flags and actual nested build are reconstructed. Core builder/parameter/physical-runner hashes match saved run exactly.',
       census=dict(neurons=len(nodes),inputs=len(n.inputs),nonzero_chemical_pairs=len(edges),zero_conductance_defined_pairs=len(disabled),spiking_neurons=0,
           layer_counts=dict(collections.Counter(x['layer'] for x in nodes))),
       compiled_flags=dict(flags,spiking=bool(c.spiking),gated=bool(c.gated)),global_G=plain(F.P.G),dt_s=float(c.dt),voltage_reference='SNS voltage relative to resting potential (0mV), not absolute−60mV.',
       compiled_matrix_sha256={key:hashlib.sha256(np.asarray(getattr(c,key)).tobytes()).hexdigest() for key in ('g_max_non','del_e','e_lo','e_hi','input_connectivity')},
       nodes=nodes,chemical_edges=edges,disabled_source_pairs=disabled,input_adapters=input_edges,actuation=policy,
       effective_type_variants=[dict(key=key,gmax_uS=sig[0],reversal_mV=sig[1],threshold_mV=sig[2],saturation_mV=sig[3],source_name='Compiled effective chemical type') for sig,key in variants.items()],
       encoder_parameters=dict(AFF=plain(F.P.AFF),MOD=plain(F.P.MOD),PHASE_RESET=plain(F.P.PHASE_RESET)),
       encoder_units_and_signs=dict(length='m',length_rate='m/s',tension='abs(native muscle force),N',current='nA',
          Ia='max(gIa * clip(length_rate/ia_vel_ref,−1,1),0)',
          II='max(i0_ii*(0.3+0.7*RGstance)+gII*clip((length−initial_pose_length)/half_length_range,−1,1),0)',
          Ib='max(0.5*gIb*clip(abs(force)/Fmax,0,1.5),0)',
          gating='Exact gIa/gII/gIb phase/drive/gamma factors are retained in the hashed runner.py and AFF/MOD/G dictionaries. No GRF is silently substituted for source tendon Ib.'),
       backend_parity_max_error_mV=capture['parity_max_mV'],source_backend_parity_max_error_mV=result['backend_parity_max_error_mV'])
    replay_prefix=prefix
    if not replay_prefix.with_suffix('.inputs.npz').exists():
        replay_prefix=ROOT/('other_family_pressure/s3k_bpa_prepared_v1' if family=='s3k' else 'prepared_spiking/syn6_graded_bpa_current_v2')
        replay_result=json.loads(replay_prefix.with_suffix('.result.json').read_text())
        assert replay_result['family']==family and replay_result.get('compiled_transform') is None
        for key in ('shared_pf_mn_gain','ib_group_gain','standing_bias_scale'):
            assert replay_result[key]==result[key]
        assert json.loads((replay_prefix.parent/(replay_prefix.name+'_projection.json')).read_text())==capture['projection']
    input_file=replay_prefix.with_suffix('.inputs.npz')
    if input_file.exists():
        with np.load(replay_prefix.with_suffix('.release_state.npz')) as release:
            for key in ('V','V_last','b_gate','c_gate','b_gate_last','c_gate_last'):
                if 'neural_'+key in release:setattr(c,key,release['neural_'+key].copy())
            first=int(release['runner_neural_steps'])-1
        with np.load(input_file) as saved,np.load(replay_prefix.with_suffix('.apparatus.npz')) as observed:
            assert saved['names'].tolist()==list(n.inputs)
            maximum=0.
            count=min(64,len(observed['neuro']))
            assert count>=32
            for j,u in enumerate(saved['inputs'][first:first+count]):
                capture['inner'].forward(u);maximum=max(maximum,float(np.max(abs(c.V-observed['neuro'][j]))))
        ledger['saved_release_replay']=dict(steps=count,max_error_mV=maximum,input_file=str(input_file),saved_input_offset=first,
            replay_condition='same biological case' if replay_prefix==prefix else 'same graded neural configuration in saved pressure-actuation sibling, not a biological/pressure state-equivalence claim',
            scope='Actual saved release state and actual saved currents; neural only, no mechanics replay or new physical simulation. A physical failure in a pressure sibling is not promoted by neural reconstruction parity.')
        assert maximum<1e-10,'Reconstructed graph does not reproduce saved neural response'
    return ledger

def node_type(name):
    if name.startswith('MN_'):return 'MN'
    if name.startswith('RC_'):return 'RC'
    if name.startswith('RG_'):return 'HC-RG-E' if '_E_' in name else 'HC-RG-F'
    if name.startswith('PF_'):return 'IN-PF' if '_IN_' in name else 'HC-PF-E'
    if name.startswith('Ia_'):return 'SN-Ia'
    if name.startswith('II_'):return 'SN-II'
    if name.startswith('Ib_'):return 'SN-Ib'
    return 'IN-IaIN'

def edge_tag(e):
    return f"{e['type_key']} | G={e['gmax_uS']:.17g}uS; Erev={e['reversal_mV']:.17g}mV; pre-threshold={e['threshold_mV']:.17g}mV; pre-saturation={e['saturation_mV']:.17g}mV | matrix row {e['target_index']}, col {e['source_index']}"

def editor_spec(ledger,names,positions=None):
    selected=set(names);nodes=[]
    for i,n in enumerate(x for x in ledger['nodes'] if x['name'] in selected):
        x,y=positions[n['name']] if positions else (i%8*1.5,-i//8*1.5)
        kind=node_type(n['name'])
        if n['name'].startswith('PF_') and '_IN_' not in n['name']:
            drives=[e['source'] for e in ledger['chemical_edges'] if e['target']==n['name'] and e['source'].startswith('RG_') and e['sign_category']=='exc']
            if any(x.startswith('RG_F_') for x in drives):kind='HC-PF-F'
        nodes.append(dict(label=n['name'],type=kind,x=round(120*(x+2)),y=round(120*(15-y)),gain=1.))
    edges=[dict(**{'from':e['source'],'to':e['target']},sign=e['sign_category'],gain=e['gmax_uS'],tag=edge_tag(e)) for e in ledger['chemical_edges'] if e['source'] in selected and e['target'] in selected]
    return dict(name=ledger['family']+' standing candidate',nodes=nodes,edges=edges,groups={})

def central_positions(ledger):
    names={n['name'] for n in ledger['nodes']};pos={};labels={}
    for name,x,y,label in [('RG_E_l',3.,10.,'RG E'),('RG_F_l',7.,10.,'RG F'),('InE_l',3.,12.,'IN E'),('InF_l',7.,12.,'IN F')]:
        if name in names:pos[name]=(x,y);labels[name]=label
    pf=[n['name'] for n in ledger['nodes'] if n['name'].startswith('PF_') and n['name'].endswith('_l') and '_IN_' not in n['name']]
    for i,name in enumerate(pf):
        x=1.+8.*i/max(1,len(pf)-1);pos[name]=(x,5.);labels[name]=name[3:-2].replace('-','\n').replace('_','\n')
    inter=[n['name'] for n in ledger['nodes'] if n['name'].startswith('PF_IN_') and n['name'].endswith('_l')]
    for i,name in enumerate(inter):
        x=2.+6.*i/max(1,len(inter)-1);pos[name]=(x,7.5);labels[name]='IN\n'+name[6:-2].replace('-','')
    return pos,labels

def motor_positions(ledger):
    names={n['name'] for n in ledger['nodes']};pos={};labels={}
    for muscle,x in [('vas_lat_l',2.5),('bifemsh_l',8.)]:
        for prefix,y,dx,short in [('MN_',3.5,0.,'MN'),('RC_',1.5,0.,'RE'),('Ia_',10.,-1.5,'Ia'),('II_',10.,0.,'II'),('Ib_',10.,1.5,'Ib'),('IaIN_',7.5,-1.5,'Ia IN'),('IIX_',7.5,0.,'IIX'),('IBIN_',7.5,1.5,'Ib IN'),('IIIN_',5.5,1.5,'II IN')]:
            name=prefix+muscle
            if name in names:pos[name]=(x+dx,y);labels[name]=('MN\n'+('VL' if muscle=='vas_lat_l' else 'BFsh')) if prefix=='MN_' else short
    # Candidate-effective PF inputs, and explicit zeroed source PF inputs.
    pf=set(e['source'] for e in ledger['chemical_edges']+ledger['disabled_source_pairs'] if e['target'] in ('MN_vas_lat_l','MN_bifemsh_l') and e['source'].startswith('PF_'))
    for i,name in enumerate(sorted(pf)):
        pos[name]=(.5+10*i/max(1,len(pf)-1),12.5);labels[name]=name[3:-2].replace('_','\n').replace('-','\n')
    return pos,labels

def circuit(ledger,folder,stem,positions,labels,title,caption,show_disabled=False):
    plt=DRAW.setup_plot();from matplotlib.patches import Circle,FancyArrowPatch
    edges=[e for e in ledger['chemical_edges'] if e['source'] in positions and e['target'] in positions]
    routes=DRAW.orthogonal_routes(positions,edges)
    fig,ax=plt.subplots(figsize=(6.5,9.5));ax.set_aspect('equal');ax.axis('off');ax.set_title(title,loc='left',fontsize=11,pad=12)
    placed=[]
    for e,route in routes:
        line,=ax.plot(*np.array(route).T,color='#333333',lw=.8,zorder=1);line.set_gid(e['id'])
        if e['sign_category']=='inh':ax.add_patch(Circle(route[-1],.045,fc='#222222',ec='#222222',zorder=4))
        else:ax.add_patch(FancyArrowPatch(route[-2],route[-1],arrowstyle='-|>',mutation_scale=9,fc='white',ec='#222222',lw=.8,zorder=4))
        points=np.asarray(route);lengths=np.linalg.norm(np.diff(points,axis=0),axis=1);chosen=None
        for segment in np.argsort(lengths)[::-1]:
            if lengths[segment]<.9:continue
            for fraction in (.5,.25,.75):
                p=(1-fraction)*points[segment]+fraction*points[segment+1];bounds=(p[0]-.35,p[0]+.35,p[1]-.16,p[1]+.16)
                if any(not(bounds[1]<b[0] or b[1]<bounds[0] or bounds[3]<b[2] or b[3]<bounds[2]) for b in placed):continue
                if any(abs(p[0]-x)<.7 and abs(p[1]-y)<.55 for x,y in positions.values()):continue
                chosen=(p,bounds);break
            if chosen is not None:break
        if chosen is not None:
            p,bounds=chosen;placed.append(bounds);ax.text(*p,e['type_key'],fontsize=10,ha='center',va='center',bbox=dict(fc='white',ec='none',pad=.4),zorder=4)
    disabled=[]
    if show_disabled:
        disabled=[e for e in ledger['disabled_source_pairs'] if e['source'] in positions and e['target'] in positions and e['source'].startswith('PF_') and e['target'].startswith('MN_')]
        for i,e in enumerate(disabled):
            a,b=positions[e['source']],positions[e['target']];middle=10.9+.12*i
            ax.plot([a[0],a[0],middle,middle,b[0]],[a[1]-.4,11.4,11.4,b[1],b[1]],color='#888888',ls=':',lw=.7,zorder=0)
    for name,(x,y) in positions.items():
        patch=Circle((x,y),.38,fc='#fff3bf' if family_class(name)=='Afferents' else '#87cefa',ec='#222222',lw=1.,zorder=5)
        patch.set_gid('neuron_'+name);ax.add_patch(patch);ax.text(x,y,labels.get(name,name),fontsize=10,ha='center',va='center',zorder=6)
    xy=np.array(list(positions.values()));ax.set_xlim(xy[:,0].min()-.6,max(xy[:,0].max()+.6,12. if disabled else 0.));ax.set_ylim(xy[:,1].min()-.6,xy[:,1].max()+.6)
    fig.subplots_adjust(bottom=.23,top=.93,left=.025,right=.98)
    note=caption+'\n'+str(len(positions))+' actual cells, '+str(len(edges))+' active induced chemical pairs.'
    if disabled:note+=' Dotted PF→MN routes: '+str(len(disabled))+' source-defined pairs with effective Gmax=0.'
    fig.text(.04,.105,'\n'.join(textwrap.fill(t,86) for t in note.splitlines()),fontsize=10,va='bottom')
    fig.text(.04,.025,textwrap.fill('All cells shown are graded. Open triangle: excitatory type; filled circle: inhibitory. Exact fields and omitted boundary inputs remain in the graph ledger. Current = g(Erev − Vpost).',86),fontsize=10,va='bottom')
    for suffix in ('svg','pdf','png'):fig.savefig(folder/(stem+'.'+suffix),dpi=200,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    data=editor_spec(ledger,positions,positions)
    for ee,(e,route) in zip(data['edges'],routes):
        assert ee['tag']==edge_tag(e);ee['pts']=[dict(x=round(120*(x+2)),y=round(120*(15-y))) for x,y in route[1:-1]]
    for e in disabled:
        data['edges'].append(dict(**{'from':e['source'],'to':e['target']},sign='exc' if e['reversal_mV']>0 else 'inh',gain=0.,tag=f"DISABLED source-defined PF→MN: G=0uS; Erev={e['reversal_mV']:.17g}mV; pre-threshold={e['threshold_mV']:.17g}mV; pre-saturation={e['saturation_mV']:.17g}mV. No synaptic current."))
    dump(folder/(stem+'_connectome.json'),data)
    return data

def overview(ledger,folder):
    plt=DRAW.setup_plot();from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
    classes=['RG','PF','MN','Afferents','Reflex INs','Descending and commissural'];nodes=ledger['nodes'];lookup={n['name']:n['layer'] for n in nodes}
    matrix=np.zeros((6,6),int)
    for e in ledger['chemical_edges']:matrix[classes.index(lookup[e['source']]),classes.index(lookup[e['target']])]+=1
    assert matrix.sum()==ledger['census']['nonzero_chemical_pairs']
    xy={'RG':(1.7,6.5),'PF':(5.2,6.5),'MN':(8.7,6.5),'Afferents':(8.7,3.4),'Reflex INs':(5.2,3.4),'Descending and commissural':(1.7,3.4)}
    counts=ledger['census']['layer_counts'];fig,ax=plt.subplots(figsize=(6.5,6.5));ax.set_xlim(0,10.8);ax.set_ylim(0,8.5);ax.axis('off')
    family=ledger['family'];ax.set_title(f'{family}: actual graded standing candidate',loc='left',fontsize=11)
    for cls,(x,y) in xy.items():
        label=('Descending /\ncommissural' if cls.startswith('Descending') else cls)+'\n'+str(counts.get(cls,0))+' cells'
        ax.add_patch(FancyBboxPatch((x-1.4,y-.55),2.8,1.1,boxstyle='round,pad=.04',fc='#fff3bf' if cls=='Afferents' else '#d9effb',ec='#333333'));ax.text(x,y,label,ha='center',va='center',fontsize=10)
    for a,b in [('RG','PF'),('PF','MN'),('Afferents','Reflex INs'),('Reflex INs','MN'),('Afferents','MN'),('Descending and commissural','RG')]:
        count=matrix[classes.index(a),classes.index(b)];x,y=xy[a];xx,yy=xy[b]
        if not count:
            if (a,b)==('PF','MN'):ax.text((x+xx)/2,y-.95,'PF→MN =0',fontsize=10,ha='center',color='#666666')
            continue
        if y==yy:route=[(x+np.sign(xx-x)*1.45,y),(xx-np.sign(xx-x)*1.45,yy)]
        elif x==xx:route=[(x,y+np.sign(yy-y)*.6),(xx,yy-np.sign(yy-y)*.6)]
        else:route=[(x,y+.6),(x,(y+yy)/2),(xx,(y+yy)/2),(xx,yy-.6)]
        ax.plot(*np.array(route).T,color='#555555',lw=.85,zorder=0);ax.add_patch(FancyArrowPatch(route[-2],route[-1],arrowstyle='->',mutation_scale=10,color='#555555'))
    c=ledger['census'];s=ledger['settings'];ax.text(.15,1.45,f"{c['neurons']} graded cells; {c['inputs']} external-current inputs\n{c['nonzero_chemical_pairs']} active chemical pairs; {c['zero_conductance_defined_pairs']} zero-conductance defined pairs\n92 MNs; 86 force-active muscle routes; 6 force-pruned routes",fontsize=10,va='center')
    fig.subplots_adjust(bottom=.23,left=.025,right=.98,top=.92)
    note=f"PF→MN multiplier {s['pf_mn_gain']:.8g}; IBEXC→MN {s['ib_group_gain']:.8g}; standing POST_ bias ×{s['standing_bias_scale']:.8g}."
    if family=='w2lvar':note+=' Antagonist IIIN gains use the saved source-fanin mean-pool normalization.'
    fig.text(.04,.145,textwrap.fill(note,86),fontsize=10);fig.text(.04,.02,textwrap.fill('Selected hierarchy arrows only. All actual connections, full numeric fields, candidate flags, and physical input/output policies are preserved in editable exports. Standing only; no spiking or walking claim.',86),fontsize=10)
    for suffix in ('svg','pdf','png'):fig.savefig(folder/('overview.'+suffix),dpi=200,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    aggregate=[];top=[]
    for cls in classes:
        x,y=xy[cls];selected=[n['name'] for n in nodes if n['layer']==cls]
        top.append(dict(label=cls+' ('+str(len(selected))+')',type='SUB',x=x*130,y=(9-y)*130,sub=editor_spec(ledger,selected)))
    for a in classes:
        for b in classes:
            for sign in ('exc','inh'):
                selected=[e for e in ledger['chemical_edges'] if lookup[e['source']]==a and lookup[e['target']]==b and e['sign_category']==sign]
                if not selected or a==b:continue
                aggregate.append(dict(**{'from':a+' ('+str(counts[a])+')','to':b+' ('+str(counts[b])+')'},sign=sign,gain=1.,tag='AGGREGATED '+str(len(selected))+' actual chemical pairs; drawing gain1 is NOT a physical conductance. Constituent rows: '+','.join(str(e['index']) for e in selected)))
    hierarchy=dict(name=family+' standing hierarchy (aggregated boundary edges)',nodes=top,edges=aggregate,groups={})
    internal=sum(len(node['sub']['edges']) for node in top)
    boundary=sum(len(edge['tag'].split('Constituent rows: ')[1].split(',')) for edge in aggregate)
    assert internal+boundary==ledger['census']['nonzero_chemical_pairs']
    dump(folder/'hierarchy_connectome.json',hierarchy);dump(folder/'layer_connection_counts.json',dict(order=classes,matrix_source_rows_target_columns=matrix.tolist()))
    return hierarchy

def verify(ledger,models):
    html=(REPO/'Code/MuJoCo_SNS/spinal/connectome_block_editor.html').read_text(encoding='utf-8');function=html[html.index('function specOf('):html.index('function serialize()',html.index('function specOf('))]
    js=function+'''const models=JSON.parse(require('fs').readFileSync(0,'utf8'));let contexts=0,edges=0;
    function check(s){const ids=new Map(s.nodes.map((n,i)=>[n.label,'n'+i]));if(ids.size!==s.nodes.length)throw Error('Duplicate node');const ns=s.nodes.map(n=>({...n,id:ids.get(n.label)}));const es=s.edges.map(e=>({...e,from:ids.get(e.from),to:ids.get(e.to)}));if(es.some(e=>!e.from||!e.to))throw Error('Missing endpoint');const a=specOf(ns,es,{},new Set());if(a.nodes.length!==s.nodes.length||a.edges.length!==s.edges.length)throw Error('Graph loss');for(let i=0;i<a.edges.length;i++)for(const k of ['from','to','sign','gain','tag'])if(a.edges[i][k]!==s.edges[i][k])throw Error('Changed '+k);contexts++;edges+=a.edges.length;for(const n of s.nodes)if(n.sub)check(n.sub);}
    models.forEach(check);process.stdout.write(JSON.stringify({contexts,edges,edge_fields_exact:true}));'''
    r=subprocess.run([r'D:\NodeJS\node-v22.23.2-win-x64\node.exe','-e',js],input=json.dumps(models),capture_output=True,text=True,check=True)
    assert len(models[-1]['nodes'])==ledger['census']['neurons'] and len(models[-1]['edges'])==ledger['census']['nonzero_chemical_pairs']
    return dict(census=ledger['census'],editor_specOf_roundtrip=json.loads(r.stdout),source_builder_hashes_match_saved_run=True,compiled_matrix_sha256=ledger['compiled_matrix_sha256'],GUI_imported=False,simulation_run=False,source_or_storage_mutation=False)

def render(family,ledger):
    folder=OUT/family;folder.mkdir(exist_ok=True)
    hierarchy=overview(ledger,folder)
    p,l=central_positions(ledger);core=circuit(ledger,folder,'left_RG_PF',p,l,f'{family}: left RG and PF subnetworks','Actual induced RG/PF circuit; candidate-effective weights. External sensory, descending, and commissural boundary connections remain in the full graph.')
    p,l=motor_positions(ledger);motor=circuit(ledger,folder,'left_knee_reflex',p,l,f'{family}: left knee motor / reflex subsystem','Actual induced vastus-lateralis / biceps-femoris-short-head circuit. PF sources are shown at top; other muscles’ antagonist-IN and group load-sharing boundary paths remain in the complete graph.',True)
    exact=editor_spec(ledger,[n['name'] for n in ledger['nodes']]);dump(folder/'exact_connectome.json',exact)
    dump(folder/'EXPORT_VERIFICATION.json',verify(ledger,[hierarchy,core,motor,exact]))
    with (folder/'chemical_edges.csv').open('w',encoding='utf-8',newline='') as stream:
        keys=['index','source','target','type_key','gmax_uS','reversal_mV','threshold_mV','saturation_mV'];w=csv.DictWriter(stream,fieldnames=keys,extrasaction='ignore');w.writeheader();w.writerows(ledger['chemical_edges'])
    dump(folder/'effective_types.json',ledger['effective_type_variants'])
    from matplotlib.backends.backend_pdf import PdfPages
    visible={(e['from'],e['to']) for data in (core,motor) for e in data['edges']}
    keys={e['type_key'] for e in ledger['chemical_edges'] if (e['source'],e['target']) in visible}
    types=[r for r in ledger['effective_type_variants'] if r['key'] in keys]
    plt=DRAW.setup_plot()
    with PdfPages(folder/'visible_panel_synaptic_fields.pdf') as pdf:
        for page,start in enumerate(range(0,len(types),18)):
            rows=types[start:start+18];fig,ax=plt.subplots(figsize=(6.5,8.5));ax.axis('off');ax.set_title(f'{family}: visible-panel synaptic fields',loc='left',fontsize=11)
            cells=[[r['key'],f"{r['gmax_uS']:.9g}",f"{r['reversal_mV']:.9g}",f"{r['threshold_mV']:.9g}",f"{r['saturation_mV']:.9g}"] for r in rows]
            table=ax.table(cellText=cells,colLabels=['Key','Gmax (µS)','Erev (mV)','Pre-thr (mV)','Pre-sat (mV)'],loc='upper left',colWidths=[.12,.28,.20,.20,.20],cellLoc='left');table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,1.9)
            fig.text(.04,.035,textwrap.fill('Keys refer to actual compiled effective types. Full precision and every endpoint are preserved in chemical_edges.csv, effective_types.json, and editor tags. Source SNS connection names are not unique named AnimatLab types; none are invented here.',85),fontsize=10);pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
    print(json.dumps(dict(family=family,rendered=True,folder=str(folder))),flush=True)

def readme():
    rows=[]
    for family in CASES:
        l=json.loads((OUT/(family+'_graph_ledger.json')).read_text(encoding='utf-8'));c=l['census'];rows.append(f"| {family} | {c['neurons']} | {c['inputs']} | {c['nonzero_chemical_pairs']} | {c['zero_conductance_defined_pairs']} |")
    text='''# Gait2392 standing-controller editable schematics

Machine: **easteregg2**. These are reconstructed saved **graded biological standing candidates**, not generic default builds, walking results, or spiking networks. Persistent-sodium RG membranes and their fixed h-gate kinetics remain graded in these selected cases.

| Family | Neurons | External inputs | Active chemical pairs | Zero-conductance defined pairs |
| --- | ---: | ---: | ---: | ---: |
'''+ '\n'.join(rows)+'''

Each family has an overview, a left RG/PF panel, and a left knee motor/reflex panel in editable SVG, PDF, and PNG. The small panels are induced circuits; omitted boundary connections remain in the complete graph. The historical light-blue, black-outline cells, inhibitory interneurons above E/F cells, Renshaw cells below MNs, orthogonal routes, open excitatory triangles, and filled inhibitory circles follow the verified LH drawing style. Labels and tables use a 10 pt floor.

The s3k and normalized w2lvar standing candidates have PF→MN and IBEXC→MN multipliers exactly zero. Their dotted PF routes are source-defined zero-conductance connections and contribute no current; they are not evidence of phasic PF motor recruitment during standing. The w2lvar antagonist IIIN→MN gains use the saved ten group-specific mean-pool factors, retaining inhibition rather than ablating it. Syn6 retains its saved weak nonzero PF drive. Every source definition, effective Gmax, reversal potential, presynaptic threshold, and saturation is stored in the graph ledger, CSV, and editor tags. SNS connection names are retained as supplied; these builders do not have unique named AnimatLab synapse types, so none are invented. T-keys identify effective numeric field tuples.

All 92 MNs and their afferent arcs remain. The physical runner force-prunes quad_fem_r/l, gem_r/l, and peri_r/l, leaving 86 biological force-active routes. This is separate from the 86-route hypothetical BPA replacement and does not depict CAD pneumatic hardware. Selected biological runs use the corrected initial pose, the recorded support apparatus, and 0.6 s pinned initialization; they are not harness-free balance claims. The source remains legacy z-up with an explicit frame distinction from the y-up CAD robot.

Physical input policy: positive lengthening velocity drives rectified Ia current; length relative to the initial standing pose drives II current with a gated tonic component; absolute source muscle tension normalized by Fmax drives Ib. Currents are nA, lengths m, velocity m/s, tension N. Exact phase/drive/gamma coefficients and normalizations are in each ledger and the hashed runner. No ground-reaction-force signal is silently substituted for tendon Ib. POST_ currents carry the saved gravity-only static-solve biases; the --reflex-stand flag disables the supraspinal BAL surrogate inputs. MN output is clip(V/E_HI,0,1), followed by native activation and muscle contraction dynamics. Voltage is relative to resting potential, not absolute −60 mV.

## Open in the block editor

Editor: `D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal\\connectome_block_editor.html`.

Save your current work, select **Templates → Blank** to create a new tab, then use **Load JSON from file…** on that tab. Open a family's `hierarchy_connectome.json` for nested layers, `left_RG_PF_connectome.json` or `left_knee_reflex_connectome.json` for the compact functional panel, or `exact_connectome.json` for the full active matrix graph. Double-click a SUB block to enter; F4 exits. Hierarchy boundary edges aggregate explicitly listed actual matrix rows; their drawing gain of 1 is not a biological conductance. Individual-node subnetworks and the exact graph preserve actual per-pair conductances. Zero source-defined PF edges are included as zero-gain entries in the motor panel; the exact graph contains only nonzero matrix entries. No browser, localStorage, or existing editor work was changed.

## Provenance and verification

Extraction executes the actual nested build function from `gait_standing_family.run` with the exact selected gains, bilateral tying, source group scales, and the two actual pre-build --reflex-stand flags, then captures the compiled network before physics. It does not run the physical simulation. Builder, parameter, and physical-runner hashes must match the saved candidate; neuron count and release-state dimensions must also match. The standing-runner wrapper has acquired later optional studies, and its current/saved hashes are recorded separately rather than falsely equated. All three reconstructions reproduce 64 actual saved post-release neural steps with **0 mV maximum error**. W2lvar uses its biological case directly; s3k and syn6 use the same graded neural configuration in saved pressure-actuation siblings for this neural-only validation, without equating biological and pressure mechanical states or promoting a failed pressure case. All saved bilateral weight projections match exactly. Muscle classification, joint-weight tables, synergy basis, and phase-profile dependencies are checked against archived snapshots as applicable. Syn6's analytical conductance helper was not archived; its current hash and that limitation are explicit in the ledger, and its effective output is independently constrained by exact saved weight projections and neural-response parity. The editor exports are functional schematics, not executable replacement runtime checkpoints; full neuron parameters and external-input rules remain in the ledgers and original source.

Every hierarchy and exact export is round-tripped through the real block editor's `specOf` JavaScript serializer in Node, without a DOM or localStorage. All endpoints, counts, signs, conductances, and full-field tags are checked. Source matrix hashes and individual checks are recorded in each family's `EXPORT_VERIFICATION.json`. Each family ledger records the exact saved result/configuration paths, their hashes, all current source hashes, flags, input adapters, source neuron parameters, and the 92/86 actuation policy. Raw experiment traces remain outside the repository.

Generator: `D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal\\virtual_walker_20261031\\export_gait_family_networks.py`.

Rebuild with `D:\\Anaconda\\envs\\myo\\python.exe` and that full script path. `--render-only` reuses the exported ledgers; `--family s3k`, `--family w2lvar`, or `--family syn6` limits the export. No production files or manuscripts are modified.
'''
    for family in CASES:
        text+='\n### '+family+'\n\nFolder: `'+str(OUT/family)+'`.\n\n'
        for stem in ('hierarchy_connectome.json','left_RG_PF.svg','left_knee_reflex.svg','visible_panel_synaptic_fields.pdf','exact_connectome.json'):
            text+='- `'+str(OUT/family/stem)+'`\n'
        text+='\nSaved result: `'+str(ROOT/CASES[family][0])+'.result.json`.\n'
    (OUT/'README.md').write_text(text,encoding='utf-8')

def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--family',choices=CASES);p.add_argument('--render-only',action='store_true');a=p.parse_args()
    OUT.mkdir(exist_ok=True)
    if a.family is None and not a.render_only:
        for family in CASES:subprocess.run([sys.executable,str(Path(__file__)),'--family',family],check=True)
        return
    for family in [a.family] if a.family else CASES:
        ledger=json.loads((OUT/(family+'_graph_ledger.json')).read_text(encoding='utf-8')) if a.render_only else extract(family)
        dump(OUT/(family+'_graph_ledger.json'),ledger)
        print(json.dumps(dict(family=family,census=ledger['census'],flags=ledger['compiled_flags'])),flush=True)
        render(family,ledger)
    if all((OUT/(f+'_graph_ledger.json')).exists() for f in CASES):readme()

if __name__=='__main__':main()
