"""Draw saved connections only. No controller imports, neural updates, or physics."""
from pathlib import Path
from collections import defaultdict
import csv,json,re,subprocess,xml.etree.ElementTree as ET
import pymupdf
from PIL import Image, ImageOps, ImageDraw
from audit_saved_graphs import ROOT, OUT, table, digest

FIG=OUT/'drawings';FIG.mkdir(exist_ok=True)
DOT=Path('C:/Users/Ben Bolen/.conda/envs/myo/Library/bin/dot.exe')
BASE=ROOT/'Neuromechanical_Models/Mujoco_SNS_models'
PANELS=[]
def q(x):return json.dumps(str(x))
def label(name):
    s=name.replace('_',' ').replace('bifemsh','BFSH').replace('vas med','vastus medialis')
    words=s.split();lines=[];current=''
    for w in words:
        if len(current)+len(w)>17:lines.append(current);current=w
        else:current=(current+' '+w).strip()
    if current:lines.append(current)
    return '\n'.join(lines)
def layer(n):
    if n.get('layer'):return n['layer']
    s=n['name'].lower()
    if 'rg' in s:return 'Rhythm generator'
    if 'pf' in s:return 'Pattern formation'
    if any(x in s for x in ('heel','toe','contact')):return 'Foot contact'
    if any(x in s for x in (' ia',' ii',' ib')):return 'Sensory / reflex'
    if ' mn' in s:return 'Motor / recurrent'
    return 'Other source cells'
def color(n):
    s=(n.get('layer','')+' '+n['name']).lower()
    if any(x in s for x in ('ia','ib','ii','affer','sensory')):return '#E8D9ED'
    if 'mn' in s or 'motor' in s:return '#F4C5B8'
    if 'rg' in s or 'rhythm' in s:return '#C3DEF0'
    if 'pf' in s or 'pattern' in s:return '#D0E4BF'
    return '#F3E5B5'
def render(slug,title,nodes,edges,note,overview=False,boundary=None):
    boundary=set(boundary or [])
    names={n['name']:n for n in nodes}
    assert len(names)==len(nodes)
    assert all(e['source'] in names and e['target'] in names for e in edges)
    lines=['digraph G {','graph [rankdir=TB, bgcolor="white", pad="0.3", nodesep="0.4", ranksep="0.7", splines=polyline, outputorder=edgesfirst, fontname="Arial", fontsize=16, labelloc=t, label='+q(title)+' ];','node [fontname="Arial", fontsize=12, margin="0.12,0.08", penwidth=1.2];','edge [fontname="Arial", fontsize=10, arrowsize=0.75, penwidth=1.1];']
    for name,n in names.items():
        shape=n.get('shape','box' if overview else 'pentagon' if name in boundary else 'ellipse')
        style='rounded,filled' if overview else 'dashed,filled' if name in boundary else 'filled'
        lines.append(q(name)+' [id='+q(n.get('id',name))+',label='+q(n.get('display',label(name)))+',shape='+shape+',style='+q(style)+',fillcolor='+q(color(n))+'];')
    for i,e in enumerate(edges):
        sign=e.get('sign_category',e.get('sign','unknown'))
        head='odot' if sign=='inh' else 'inv' if sign=='exc' else 'normal' if sign=='signal' else 'none'
        g=e.get('gmax_uS',e.get('gain'))
        txt=str(e['count'])+' edges' if overview else ''
        if not overview and len(edges)<=30 and g is not None:txt=f'{g:.3g}'
        lines.append(q(e['source'])+' -> '+q(e['target'])+' [id='+q(e.get('id',f'edge:{i}'))+',arrowhead='+head+',label='+q(txt)+',style='+q('dashed' if sign=='signal' else 'solid')+',color='+q('#4477AA' if sign=='signal' else '#555555' if sign=='unknown' else '#222222')+'];')
    lines+=['legend [shape=plaintext, fontsize=11, label='+q(note+'\nArrow: saved excitatory category. Open circle: saved inhibitory category.\nNeuron colours: RG blue; PF green; motor salmon; afferent purple.\nDashed pentagons: boundary references, not additional neurons.\nEdge numbers: gmax (uS). Full fields and IDs accompany each panel.')+'];','}']
    dot=FIG/(slug+'.dot');dot.write_text('\n'.join(lines),encoding='utf-8')
    for fmt in ('svg','pdf'):
        subprocess.run([str(DOT),'-T'+fmt,str(dot),'-o',str(FIG/(slug+'.'+fmt))],check=True,capture_output=True)
    doc=pymupdf.open(FIG/(slug+'.pdf'));doc[0].get_pixmap(matrix=pymupdf.Matrix(1.5,1.5)).save(FIG/(slug+'.png'));size=list(doc[0].rect);doc.close()
    table(FIG/(slug+'_edges.csv'),edges)
    (FIG/(slug+'_editable.json')).write_text(json.dumps({'title':title,'nodes':nodes,'edges':edges,'boundary_references':sorted(boundary),'source_scope':note},indent=2),encoding='utf-8')
    PANELS.append({'slug':slug,'title':title,'nodes':len(nodes),'edges':len(edges),'pdf_points':size,'source_scope':note})

def overview(family,nodes,edges,scope):
    groups=defaultdict(list)
    for n in nodes:groups[layer(n)].append(n)
    lookup={n['name']:k for k,v in groups.items() for n in v}
    bundles=defaultdict(list)
    for e in edges:
        a,b=lookup[e['source']],lookup[e['target']]
        if a!=b:bundles[a,b].append(e['id'])
    ns=[{'name':k,'layer':k,'display':k+'\n'+str(len(v))+' neurons'} for k,v in groups.items()]
    es=[{'id':family+':bundle:'+str(i),'source':a,'target':b,'count':len(ids),'underlying_edge_ids':ids,'sign':'unknown'} for i,((a,b),ids) in enumerate(bundles.items())]
    render(family+'_overview',family+' | connection overview',ns,es,scope+'\nModule links bundle actual edges; within-module edges are in expanded panels/tables.',True)

def subset(family,tag,nodes,edges,chosen,scope,incoming=False):
    selected=set(chosen);boundary=set()
    if incoming:
        boundary={e['source'] for e in edges if e['target'] in selected and e['source'] not in selected}
    ids=selected|boundary
    es=[e for e in edges if e['source'] in ids and e['target'] in selected]
    render(family+'_'+tag,family+' | '+tag.replace('_',' '),[n for n in nodes if n['name'] in ids],es,scope,boundary=boundary)

for family in ('s3k','w2lvar','syn6','W2L704_80'):
    p=BASE/('W2L_704_80_baseline/W2L_704_80_graph_ledger.json' if family=='W2L704_80' else 'Gait2392_standing_families/'+family+'_graph_ledger.json')
    d=json.loads(p.read_text());nodes=d['nodes'];edges=d['chemical_edges']
    scope='Saved standing condition' if family!='W2L704_80' else 'Saved W2L hybrid condition; source neuron interfaces shown, not 704 individual somas'
    overview(family,nodes,edges,scope)
    rg=[n['name'] for n in nodes if n.get('layer')=='RG' or ' RG ' in ' '+n['name']+' ']
    subset(family,'bilateral_RG',nodes,edges,rg,scope+' | RG neurons and their actual incoming boundary connections.',True)
    knee=[n['name'] for n in nodes if any(t in n['name'] for t in ('vas_med_l','bifemsh_l','IBEXC_knee_ext_l'))] if family!='W2L704_80' else [n['name'] for n in nodes if n['name'].startswith('L ') and any(t in n['name'] for t in ('Knee MN','Knee ext I','Knee flx I','Knee II'))]
    subset(family,'left_knee_reflex',nodes,edges,knee,scope+' | Selected left knee motor/afferent circuits and incoming boundary references.',True)

# Native rat: endpoints and type parameters come directly from the saved standalone.
rat=ROOT/'tmp/functional_subnetwork_remote/named_conditions/rat134/Rat134_ReceptorGeometry_Standalone.asim'
r=ET.parse(rat).getroot();ns=list(r.iter('Neuron'));names={x.findtext('ID'):x.findtext('Name') for x in ns};types={x.findtext('ID'):x for x in r.iter('SynapseType')};nmap={x.findtext('ID'):x for x in ns}
nodes=[{'id':x.findtext('ID'),'name':x.findtext('Name')} for x in ns];edges=[]
for x in r.iter('Connexion'):
    t=types[x.findtext('SynapseTypeID')];target=x.findtext('TargetID');rev=float(t.findtext('Equil'));rest=float(nmap[target].findtext('RestingPot'))
    edges.append({'id':x.findtext('ID'),'source':names[x.findtext('SourceID')],'target':names[target],'gmax_uS':float(t.findtext('SynAmp')),'reversal_mV':rev,'threshold_mV':float(t.findtext('ThreshV')),'saturation_mV':float(t.findtext('SaturateV')),'type_name':t.findtext('Name'),'sign_category':'exc' if rev>rest else 'inh','sign_scope':'relative to target resting Vm only; operating Vm can change current sign'})
overview('Rat134',nodes,edges,'Curated134 receptor-geometry baseline; refinement12 gains not installed')
for side in ('LH','RH'):
    chosen=[n['name'] for n in nodes if n['name'].startswith(side) and ('RG' in n['name'] or 'Hip PF' in n['name'])]
    subset('Rat134',side+'_RG_hipPF',nodes,edges,chosen,'Curated134 | native source graph; arrows/circles classified at target rest. Boundary afferents retained.',True)

# Direct-ODE mechanism reference, extracted as literal wiring without importing it.
with (OUT/'reduced_Deng14/nodes.csv').open(encoding='utf-8-sig') as f:nodes=[dict(x) for x in csv.DictReader(f)]
with (OUT/'reduced_Deng14/edges.csv').open(encoding='utf-8-sig') as f:raw=list(csv.DictReader(f))
rest={n['name']:float(n['rest_mV']) for n in nodes}
edges=[{**x,'gmax_uS':float(x['gmax_uS']),'sign_category':'exc' if float(x['Erev_mV'])>rest[x['target']] else 'inh'} for x in raw]
render('Deng14_complete','Reduced Deng14 | full 14-neuron circuit',nodes,edges,'Fixed-h direct ODE reference, 350ms h; no physical afferents or whole-body connection.\nArrow/circle classification is relative to resting Vm. ODE gains retained.')

for side in ('left','right'):
    p=BASE/('Hunt_fourcell_session1_'+side+'_RF_HAM_connectome.json');d=json.loads(p.read_text())
    nodes=[{'name':n['label'],'type':n['type']} for n in d['nodes']]
    edges=[{**e,'id':'fourcell:'+side+':'+str(i),'source':e['from'],'target':e['to'],'gmax_uS':e['gain']} for i,e in enumerate(d['edges'])]
    render('RF_HAM_'+side,'RF/HAM four-cell PF | '+side,nodes,edges,'Recorded Session1 editor projection, not a newly asserted current circuit.\nRF/HAM retained separately from sartorius/gracilis; supplied excitation/inhibition categories.')

d=json.loads((BASE/'Gait2392_standing_families/syn6_graph_ledger.json').read_text())
nodes=[dict(n) for n in d['nodes']]
for n in nodes:
    if n['name'].startswith('RG_'):n['display']=label(n['name'])+'\n16 spiking somas\nNaP h retained'
subset('Gait_v4','central_RG_interfaces',nodes,d['chemical_edges'],[n['name'] for n in nodes if n['name'].startswith('RG_')],'TimingZeroSFA | source population interfaces; each RG contains16 spiking somas.\nChemical edges use spike-derived proxy voltage; these are NOT event synapses.\nZero extra SFA; h tau multiplier0.9496883348174533; 0.2ms internal clock.',True)

cfg=json.loads((ROOT/'tmp/functional_subnetwork_remote/named_conditions/cad_normalization/RECRUITMENT_PREREGISTRATION.json').read_text())
nodes=[];edges=[]
for side in ('r','l'):
    mn='MN_vas_med_'+side
    nodes.append({'name':mn})
    for source,g in cfg['absolute_shared_weights_uS'].items():
        pre=source+'_'+side;nodes.append({'name':pre})
        edges.append({'id':'CADpair:'+pre+':'+mn,'source':pre,'target':mn,'gmax_uS':g,'sign':'exc','source_evidence':'cad_pair_unclamped.py FreePair.__init__; RECRUITMENT_PREREGISTRATION.json; saved endpoint names'})
    previous=mn
    for component in ('Inverse force to pressure','Pressure dynamics','vas med BPA'):
        name=component+' '+side;nodes.append({'name':name,'display':label(name),'shape':'box'})
        edges.append({'id':'CADpair:signal:'+side+':'+component,'source':previous,'target':name,'sign':'signal','source_evidence':'Saved assumed pressure-pair configuration and reproduction manifest; engineering signal path, not chemical synapse'})
        previous=name
render('CAD262_pressure_pair','CAD262 | the two neural pressure-controlled routes',nodes,edges,'Support085 mixed condition | only bilateral vas_med pair is expanded here.\nDashed blue arrows: engineering signal conversion. Pair POST current=0.\n52 other routes retain external hybrid allocation and documented force floor.\nMeasured upward support87.3372%BW; guidance retained; not an all-BPA controller.')

(FIG/'DRAWING_INDEX.json').write_text(json.dumps(PANELS,indent=2))
combined=pymupdf.open()
for panel in PANELS:
    with pymupdf.open(FIG/(panel['slug']+'.pdf')) as src:combined.insert_pdf(src)
combined.save(FIG/'Connections.pdf');combined.close()
thumbs=[]
for panel in PANELS:
    im=Image.open(FIG/(panel['slug']+'.png')).convert('RGB');im.thumbnail((660,450))
    tile=Image.new('RGB',(700,500),'white');tile.paste(im,((700-im.width)//2,30));ImageDraw.Draw(tile).text((12,480),panel['title'],fill='black');thumbs.append(tile)
sheet=Image.new('RGB',(1400,500*((len(thumbs)+1)//2)), '#e4e4e4')
for i,t in enumerate(thumbs):sheet.paste(t,((i%2)*700,(i//2)*500))
sheet.save(FIG/'Contact_sheet.png')
print(json.dumps({'panels':len(PANELS),'combined_pdf':str(FIG/'Connections.pdf'),'circuits_run':0}))
