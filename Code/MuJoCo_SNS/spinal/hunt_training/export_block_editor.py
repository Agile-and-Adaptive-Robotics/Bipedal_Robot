"""Export actual compiled edges into Ben's v3.3 editor, without altering his workspace."""
import json,re,hashlib,argparse,collections,os
from pathlib import Path
H=Path(__file__).resolve().parent;S=H.parent;repo=S.parents[2]
p=argparse.ArgumentParser();p.add_argument('--graph',type=Path,default=H/'runs/fourcell_dynamic_20261029_v1/gate_verification/standing.graph.json');args=p.parse_args()
source=args.graph.resolve();graph=json.loads(source.read_text());labels=[n['label'] for n in graph['nodes']]
dest=repo/'Neuromechanical_Models/Mujoco_SNS_models';assert dest.is_dir()
editor=S/'connectome_block_editor.html';text=editor.read_text(encoding='utf-8');types=set(re.findall(r'^  "([A-Za-z0-9+_-]+)":',text,re.M))
def group(name):
    side='R' if name.endswith('_r') else 'L' if name.endswith('_l') else 'Shared'
    if name.startswith(('PT_BAL_','PT_IN_')):return 'Shared task balance'
    if name.startswith(('RG_','InE_','InF_','HUNT_COMM_','HUNT_AF_')):return side+' rhythm and timing'
    if name.startswith('HUNT_BANK_'):return side+' phase filters'
    if name.startswith(('PF_BI_JOINT_HAM_','PF_BI_JOINT_RECT_','PF_IN_BI_PAIR_')):return side+' RF HAM four-cell PF'
    if name.startswith('PF') and 'BI_' in name:return side+' other biarticular PF'
    if name.startswith('PF'):return side+' joint PF'
    if name.startswith('MN_'):return side+' motor neurons'
    if name.startswith(('Ia_','II_','Ib_','HEEL_','TOE','LOAD_')):return side+' sensory neurons'
    if name.startswith(('IaIN_','IIIN_','IIX_','IBIN_','IBEXC_','RC_')):return side+' reflex interneurons'
    return 'Shared descending and posture'
def typ(n):
    if n.startswith('RG_E_'):return 'HC-RG-E'
    if n.startswith('RG_F_'):return 'HC-RG-F'
    if n.startswith('InE_'):return 'IN-InE'
    if n.startswith('InF_'):return 'IN-InF'
    if n.startswith(('PF_IN_','HUNT_BANK_IN_')):return 'IN-PF'
    if n.startswith('PF'):return 'HC-PF-F' if any(x in n for x in ('HIP-F','KNEE-F','HAM')) else 'HC-PF-E'
    if n.startswith('MN_'):return 'MN'
    if n.startswith('Ia_'):return 'SN-Ia'
    if n.startswith('II_'):return 'SN-II'
    if n.startswith('Ib_'):return 'SN-Ib'
    if n.startswith('IaIN_'):return 'IN-IaIN'
    if n.startswith('IIIN_'):return 'IN-IIi'
    if n.startswith('IIX_'):return 'IN-IIe'
    if n.startswith('IBIN_'):return 'IN-IbIN'
    if n.startswith('IBEXC_'):return 'IN-Ib+'
    if n.startswith('RC_'):return 'RC'
    if n.startswith('HEEL_'):return 'SN-heel'
    if n.startswith('TOE'):return 'SN-toe'
    return 'IN-C'
groups=sorted(set(map(group,labels)));gids={g:'g'+str(i) for i,g in enumerate(groups)}
enabled=[e for e in graph['edges'] if max(e['standing'],e['walking'])>1e-9]
def edge(e):return dict(**{'from':labels[e['source']],'to':labels[e['target']]},sign='exc' if e['reversal']>0 else 'inh',gain=max(e['standing'],e['walking']),tag='gmax',pts=[])
nodes=[];serial=collections.Counter();row_y=0;orig={}
for start in range(0,len(groups),4):
    row=groups[start:start+4]
    for col,g in enumerate(row):
        for i,label in enumerate([n for n in labels if group(n)==g]):
            nd=dict(type=typ(label),label=label,x=180+col*2000+(i%7)*250,y=160+row_y+(i//7)*100,grp=gids[g]);nodes.append(nd);orig[label]=nd
    row_y+=max((sum(group(n)==g for n in labels)+6)//7 for g in row)*100+300
flat=dict(nodes=nodes,edges=list(map(edge,enabled)),groups={gids[g]:g for g in groups},scope=graph['scope'],source=str(source),representation='Exact enabled neuron graph; no channel dynamics reconstructed')
overview=dict(nodes=[],edges=[],scope='Module overview. Boundary PORT nodes are display references, not extra neurons. Bundled overview edges count connections, not synaptic conductance. Exact flat graph is authoritative.')
for gi,g in enumerate(groups):
    owned=[n for n in labels if group(n)==g];spec=dict(nodes=[],edges=[],groups={'inside':g,'input':'Boundary inputs','output':'Boundary outputs'})
    for i,label in enumerate(owned):spec['nodes'].append(dict(type=typ(label),label=label,x=550+(i%6)*240,y=180+(i//6)*100,grp='inside'))
    ports={}
    for e in enabled:
        a,b=labels[e['source']],labels[e['target']];ga,gb=group(a),group(b)
        if ga!=g and gb!=g:continue
        en=edge(e)
        if ga!=g:
            key=('input',ga);ports.setdefault(key,'INPUT '+ga);en['from']=ports[key];en['tag']='Boundary source '+a+' -> '+b
        if gb!=g:
            key=('output',gb);ports.setdefault(key,'OUTPUT '+gb);en['to']=ports[key];en['tag']='Boundary target '+a+' -> '+b
        spec['edges'].append(en)
    for direction in ('input','output'):
        for j,((d,g2),label) in enumerate([(k,v) for k,v in ports.items() if k[0]==direction]):spec['nodes'].append(dict(type='PORT-load',label=label,x=100 if direction=='input' else 2200,y=180+j*140,grp=direction))
    overview['nodes'].append(dict(type='SUB',label=g,x=220+(gi%4)*400,y=150+(gi//4)*180,sub=spec))
bundles=collections.Counter((group(labels[e['source']]),group(labels[e['target']]),'exc' if e['reversal']>0 else 'inh') for e in enabled if group(labels[e['source']])!=group(labels[e['target']]))
overview['edges']=[dict(**{'from':a,'to':b},sign=sign,gain=count,tag=f'DISPLAY BUNDLE: {count} neuron edges; gain is count, not conductance',pts=[]) for (a,b,sign),count in bundles.items()]
models=[dict(name='Hunt modules (double-click SUB)',data=overview),dict(name='Hunt exact full graph',data=flat)]
for side in ('r','l'):
    selected={n for n in labels if n.endswith('_'+side) and (n.startswith(('PF_BI_JOINT_HAM_','PF_BI_JOINT_RECT_','PF_IN_BI_PAIR_')) or n in ['PF_'+j+'_'+side for j in ('HIP-E','HIP-F','KNEE-E','KNEE-F')])}
    # Show the four anatomical RF/hamstring motor pools, including pools whose
    # fitted edges are below threshold; other motors remain in the exact graph.
    selected.update('MN_'+family+'_'+side for family in ('rect_fem','semimem','semiten','bifemlh'))
    focus=dict(nodes=[],edges=[edge(e) for e in enabled if labels[e['source']] in selected and labels[e['target']] in selected],groups={'inputs':'Joint PF inputs','paired':'Four-cell RF/HAM PF','motor':'Motor outputs'},scope='Induced subnetwork; boundary afferents and other PF inputs are omitted, not deleted.')
    positions={
        'PF_HIP-E':(180,150),'PF_KNEE-F':(180,260),'PF_HIP-F':(180,440),'PF_KNEE-E':(180,550),
        'PF_BI_JOINT_HAM':(600,190),'PF_BI_JOINT_RECT':(600,510),
        'PF_IN_BI_PAIR_HAM':(820,190),'PF_IN_BI_PAIR_RECT':(820,510),
        'MN_semimem':(1180,150),'MN_semiten':(1180,270),'MN_bifemlh':(1180,390),'MN_rect_fem':(1180,550)}
    for label in labels:
        if label not in selected:continue
        x,y=positions[label[:-2]];kind='motor' if label.startswith('MN_') else 'paired' if 'BI_' in label else 'inputs'
        focus['nodes'].append(dict(type=typ(label),label=label,x=x,y=y,grp=kind))
    for e in focus['edges']:
        if e['from']=='PF_BI_JOINT_HAM_'+side and e['to'].startswith('MN_'):
            route_y=70+40*('semimem','semiten','bifemlh').index(e['to'][3:-2]);e['pts']=[dict(x=940,y=route_y),dict(x=1060,y=route_y)]
        if e['from']=='PF_BI_JOINT_RECT_'+side and e['to'].startswith('MN_'):e['pts']=[dict(x=940,y=640),dict(x=1060,y=640)]
    models.append(dict(name='RF HAM '+('right' if side=='r' else 'left'),data=focus))
def validate(spec):
    ns=spec['nodes'];lab={n['label'] for n in ns};assert len(lab)==len(ns)
    assert all(n['type'] in types for n in ns)
    for e in spec['edges']:assert e['from'] in lab and e['to'] in lab and e['gain']>=0,e
    for n in ns:
        if 'sub' in n:validate(n['sub'])
for model in models:validate(model['data'])
assert len(flat['nodes'])==len(graph['nodes']) and len(flat['edges'])==len(enabled)
stem='Hunt_fourcell_session1';files=[]
for suffix,spec in [('full',flat),('subsystems',overview),('right_RF_HAM',models[2]['data']),('left_RF_HAM',models[3]['data'])]:
    path=dest/(stem+'_'+suffix+'_connectome.json');tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(spec,indent=2),encoding='utf-8');os.replace(tmp,path);files.append(path.name)
workspace=dict(models=models,active=0);content=json.dumps(models).replace('</','<\\/');seed='let models = [{ name: "Ben", data: { nodes: [], edges: [] } }];';assert text.count(seed)==1
text=text.replace(seed,'let models = '+content+';').replace('cme_models_v1','cme_models_hunt_session1_'+hashlib.sha256(content.encode()).hexdigest()[:12])
original_model_spec=re.search(r'function modelSpec\(\) \{[\s\S]*?\n\}',text).group(0)
fixed_model_spec='''function modelSpec() {
  let ns = nodes, es = edges, gn = grpNames, hid = hiddenG;
  for (let i = ctxStack.length - 1; i >= 0; i--) {
    const f = ctxStack[i];
    f.hostNode.sub = specOf(ns, es, gn, hid);
    ns = f.nodes; es = f.edges; gn = f.grpNames; hid = f.hiddenG;
  }
  return specOf(ns, es, gn, hid);
}'''
text=text.replace(original_model_spec,fixed_model_spec,1)
text=text.replace('<title>Connectome Block Editor v3.3','<title>Hunt four-cell — Connectome Block Editor v3.3',1)
out=dest/(stem+'_editor.html');temp=out.with_suffix('.tmp');temp.write_text(text,encoding='utf-8');os.replace(temp,out)
(dest/(stem+'_workspace.json')).write_text(json.dumps(workspace,indent=2),encoding='utf-8')
(dest/(stem+'_README.md')).write_text('''# Session 1: Hunt four-cell drawing in the v3.3 block editor

Machine: easteregg2. Open Hunt_fourcell_session1_editor.html for the preloaded
editor, or use the Hunt Session1 templates registered in the original block
editor under Templates > SNS walker versions. The prior embedded template
library, canonical 110/216 model, and user browser storage are preserved.
The copy uses a separate content-specific
storage key; no TEST_* entries are created.

Both the original editor and standalone copy fix modelSpec write-back so
saving while inside a subsystem retains the outer model.
This is checked with the editor's actual JavaScript serialization functions.

The module overview contains native SUB nodes; double-click to enter, F4 to
return. Its boundary ports are drawing references. Overview gain fields count
bundled edges and are not conductances. The exact full-graph tab preserves
every enabled compiled neuron edge and is authoritative. RF/HAM tabs show
induced subgraphs. Sartorius and gracilis remain separate within the other
biarticular subsystem and in the exact graph.

Conductances are the maximum of the recorded standing and gait-expression
views at one pose, not a guarantee of recruitment through a gait cycle.
Disabled edges below 1e-9 are excluded. These are editor drawings, not runnable
SNS reconstructions: channel equations and external mechanics remain in code.
The supplied source candidate has failed ground stepping; this is not a
walking result. Live localhost UI verification covered module entry and both RF/HAM detail
views after Ben opened the HTTP page. The full flat graph was checked offline.
The saved screenshots capture the tested state before later user editing.
''',encoding='utf-8')
manifest=dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),neurons=len(labels),enabled_edges=len(enabled),modules=len(groups),files=files,preloaded_editor=out.name,all_endpoints_valid=True,all_types_valid=True)
(dest/(stem+'_manifest.json')).write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
