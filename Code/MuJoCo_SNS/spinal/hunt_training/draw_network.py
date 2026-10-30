"""Draw actual compiled connectivity and SNS-Toolbox subnetwork projections."""
import os,sys,json,collections,html
from pathlib import Path
H=Path(__file__).resolve().parent;O=Path(sys.argv[1]) if len(sys.argv)>1 else H/'runs/network_diagram';V=Path(sys.argv[2]) if len(sys.argv)>2 else H/'network_figures'/('candidate' if len(sys.argv)>1 else '');V.mkdir(exist_ok=True)
os.environ['PATH']=str(Path(r'C:\Program Files\Graphviz\bin'))+os.pathsep+str(Path(sys.executable).parent/'Library/bin')+os.pathsep+os.environ.get('PATH','')
import graphviz
from sns_toolbox.networks import Network
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.renderer import render
z=json.loads((O/'standing.graph.json').read_text());names=[n['label'] for n in z['nodes']]
def group(n):
    if n.startswith('HUNT_COMM_'):return 'Candidate commissural inhibition'
    if n.startswith(('RG_','InE_','InF_')):return 'RG + local inhibition'
    if n.startswith('HUNT_AF_'):return 'Calibrated timing afferents'
    if n.startswith('HUNT_BANK_'):return 'Causal phase filters'
    if n.startswith(('PT_BAL_','PT_IN_')):return 'Task balance + relays'
    if n.startswith('PF_HUNT_TONIC'):return 'Tonic feedforward'
    if 'BI_' in n and n.startswith('PF'):return 'Biarticular PF / integrators'
    if n.startswith('PF'):return 'Joint PF + local inhibition'
    if n.startswith('MN_'):return 'Motor neurons'
    if n.startswith(('Ia_','II_','Ib_')):return 'Muscle sensory populations'
    if n.startswith(('IaIN_','IIIN_','IIX_','IBIN_','IBEXC_')):return 'Sensory interneurons'
    if n.startswith('RC_'):return 'Renshaw network'
    if n.startswith(('HEEL_','TOE','LOAD_')):return 'Foot / load populations'
    return 'Descending / posture / balance'
groups=sorted({group(n) for n in names});counts=collections.Counter(map(group,names))
edges=[e for e in z['edges'] if max(e['standing'],e['walking'])>1e-9]
g=graphviz.Digraph('implemented',format='svg');g.attr(rankdir='LR',bgcolor='white',label='biw2lvar-Hunt-PT: enabled compiled connections\nGreen: excitatory; red: inhibitory. Counts aggregate real synapses; enabled does not prove functional recruitment.',labelloc='t',fontsize='17',fontname='Arial')
g.attr('node',shape='box',style='rounded,filled',fillcolor='#EEF3F8',fontname='Arial',fontsize='12')
for i,name in enumerate(groups):g.node(str(i),name+'\n'+str(counts[name])+' neurons')
agg=collections.Counter()
for e in edges:
    a,b=group(names[e['source']]),group(names[e['target']]);sign='exc' if e['reversal']>0 else 'inh'
    if a!=b:agg[(a,b,sign)]+=1
for (a,b,sign),count in agg.items():g.edge(str(groups.index(a)),str(groups.index(b)),label=str(count),color='#228833' if sign=='exc' else '#CC3311',arrowhead='normal' if sign=='exc' else 'tee',fontsize='10')
g.render(str(V/'overview'),cleanup=True);g.format='png';g.render(str(V/'overview'),cleanup=True)
def project(label,selected):
    selected=set(selected);net=Network(name=label)
    if not selected:return dict(neurons=0,edges=0,scope='Not built in this stage snapshot.')
    def display(name):
        if name.startswith('HUNT_BANK_IN_'):return name.replace('HUNT_BANK_IN_','Filter inhibitory relay ')
        if name.startswith('HUNT_BANK_'):return name.replace('HUNT_BANK_','RG phase filter ')
        return name
    for i in sorted(selected):net.add_neuron(NonSpikingNeuron(),name=display(names[i]))
    edgecount=0
    for e in edges:
        if e['source'] in selected and e['target'] in selected:
            net.add_connection(NonSpikingSynapse(max_conductance=max(e['standing'],e['walking']),reversal_potential=e['reversal']),display(names[e['source']]),display(names[e['target']]));edgecount+=1
    render(net,view=False,save=True,filename=str(V/label),img_format='svg')
    render(net,view=False,save=True,filename=str(V/label),img_format='png')
    return dict(neurons=len(selected),edges=edgecount,scope='Actual induced subnetwork; omitted boundary connections are not deleted from the model. Graphics use generic neuron shapes; dynamics remain in the source model.')
panels={}
panels['rhythm']=project('rhythm',[i for i,n in enumerate(names) if n.startswith(('RG_','InE_','InF_','HUNT_AF_','HUNT_COMM_'))])
panels['joint_and_biarticular']=project('joint_and_biarticular',[i for i,n in enumerate(names) if n.endswith('_r') and n.startswith('PF') and any(k in n for k in ('HIP-E','HIP-F','KNEE-E','KNEE-F','ANKLE_Z','BI_JOINT','BI_PAIR'))])
paired=any(n.startswith('PF_IN_BI_PAIR_') for n in names)
if paired:panels['paired_rf_ham']=project('paired_rf_ham',[i for i,n in enumerate(names) if n.endswith('_r') and (n.startswith('PF_IN_BI_PAIR_') or n.startswith(('PF_BI_JOINT_HAM_','PF_BI_JOINT_RECT_')) or n in ('PF_HIP-E_r','PF_HIP-F_r','PF_KNEE-E_r','PF_KNEE-F_r'))])
panels['hip_axis_pf']=project('hip_axis_pf',[i for i,n in enumerate(names) if n.endswith('_r') and ('HIP_X' in n or n.startswith('HUNT_BANK_'))])
panels['other_biarticular']=project('other_biarticular',[i for i,n in enumerate(names) if n.endswith('_r') and n.startswith('PF') and any(k in n for k in ('BI_SAR','BI_GRAC','BI_TFL','HIP-E','HIP-F','KNEE-E','KNEE-F','HIP_X'))])
selected={i for i,n in enumerate(names) if n.startswith(('PT_BAL_','PT_IN_'))}
for source in list(selected):
    targets=sorted([e for e in edges if e['source']==source and names[e['target']].startswith('MN_')],key=lambda e:max(e['standing'],e['walking']),reverse=True)[:2]
    selected.update(e['target'] for e in targets)
panels['balance_example']=project('balance_example',selected)
layout={};row_y=0
for start in range(0,len(groups),4):
    row=groups[start:start+4]
    for col,name in enumerate(row):layout[name]=(col*1700,row_y)
    row_y+=max((counts[name]+9)//10*80 for name in row)+240
serial=collections.Counter();nodes=[]
for i,name in enumerate(names):
    category=group(name);j=serial[category];serial[category]+=1;x,y=layout[category]
    typ='MN' if name.startswith('MN_') else 'HC-RG-E' if name.startswith('RG_E') else 'HC-RG-F' if name.startswith('RG_F') else 'HC-PF-E' if name.startswith('PF') else 'IN-C'
    nodes.append(dict(id='n'+str(i),type=typ,label=name,x=x+(j%10)*150,y=y+(j//10)*80,subnetwork=category))
editor=dict(nodes=nodes,edges=[dict(**{'from':names[e['source']],'to':names[e['target']]},sign='exc' if e['reversal']>0 else 'inh',gain=max(e['standing'],e['walking']),tag='compiled',pts=[],standing_gain=e['standing'],walking_expression_gain=e['walking']) for e in edges],scope=z['scope'],purpose='Visualization export; does not reconstruct full neuronal/channel equations or a runnable SNS model.',drawing_threshold=1e-9)
(O/'biw2lvar_Hunt_PT_connectome.json').write_text(json.dumps(editor))
summary=dict(neurons=len(names),enabled_edges=len(edges),disabled_edges=len(z['edges'])-len(edges),groups=dict(counts),panels=panels,commissural_neurons=[n for n in names if n.startswith(('C1_','V3_','V0','HUNT_COMM_'))],scope=z['scope'])
(V/'NETWORK_CENSUS.json').write_text(json.dumps(summary,indent=2))
reciprocal=any(names[e['source']].startswith('IaIN_') and names[e['target']].startswith('MN_') for e in edges)
intro='These diagrams come from compiled conductances, not literature expectations. '+('A candidate commissural pathway is built. ' if summary['commissural_neurons'] else 'No active commissural subnetwork is built. ')+('Ia reciprocal inhibition is enabled in this candidate. ' if reciprocal else 'Ia reciprocal inhibition outputs to MNs are disabled. ')+'IIX/IIIN/Renshaw outputs to MNs remain disabled. The raw export retains disabled edges; the editor export contains enabled edges only. Walking-expression gating at this recorded pose is not proof of walking.'
page='<html><meta charset="utf-8"><title>biw2lvar-Hunt-PT network</title><style>body{font:16px Arial;margin:24px;color:#142638}img{max-width:100%;border:1px solid #ddd}summary{font-size:20px;padding:15px;cursor:pointer}a{color:#17649c}p{max-width:1000px;line-height:1.5}</style><h1>biw2lvar-Hunt-PT: implemented subnetworks</h1><p>'+intro+'</p><p>1,029 neurons · 424 inputs. Overview: green excites, red inhibits. Detailed panels use the SNS-Toolbox renderer and omit boundary edges for readability.</p><p>Timing afferents 0–5: hip-flexor II, hip-extensor II, hip-flexor Ia, hip-extensor Ia, ankle-extensor Ib, and knee-extensor Ib. HUNT_BANK cells are causal phase filters, not extra biological oscillators. Joint PF pairs retain inhibitory relays. HAM/RECT/GAS use single integrating populations; SAR/GRAC/TFL retain paired engineering modules.</p>'
page=page.replace('1,029 neurons · 424 inputs',str(len(names))+' neurons · '+str(len(z['inputs']))+' inputs')
if paired:page=page.replace('HAM/RECT/GAS use single integrating populations; SAR/GRAC/TFL retain paired engineering modules.','RF and HAM form one four-cell module per leg: two integrating PF cells and two reciprocal inhibitory relays. GAS remains a single integrator. Sartorius, gracilis, and TFL retain separate paired engineering modules.')
page+='<p>The added time-filter cells are engineering features, not neurons specified by Hunt. The PT input computation is a body-level controller outside SNS; it is separate from local muscle-stretch feedback. <a href="'+os.path.relpath(H/'BALANCE_AND_PF_CLARIFICATION.md',V).replace('\\','/')+'">Read the clarification and measured muscle functions</a>.</p>'
for extra,title in [('balance_dataflow.svg','What the balance layer computes'),('biarticular_four_cell_candidate.svg','Four-cell RF/HAM motif (topology schematic)')]:
    if (V/extra).exists():page+='<details><summary>'+title+'</summary><img src="'+extra+'"></details>'
for label in ['overview',*panels]:
    page+='<details '+('open' if label=='overview' else '')+'><summary>'+html.escape(label.replace('_',' '))+'</summary>'
    page+=('<img src="'+label+'.svg">' if label=='overview' or panels[label]['neurons'] else '<p>Not built in this stage snapshot.</p>')+'</details>'
editor_link=os.path.relpath(O/'biw2lvar_Hunt_PT_connectome.json',V).replace('\\','/');raw_link=os.path.relpath(O/'standing.graph.json',V).replace('\\','/')
page+='<p><a href="'+editor_link+'">Import JSON for your connectome editor</a> · <a href="'+raw_link+'">Full compiled edge data, including disabled connections</a></p>'
(V/'index.html').write_text(page,encoding='utf-8');(V/'overview.alt.txt').write_text(intro)
print(json.dumps(summary,indent=2))
