"""Prepare editable hierarchy copies from saved data, without executing circuits."""
from __future__ import annotations
import csv,hashlib,json,re,copy
from pathlib import Path
from collections import defaultdict
from xml.etree import ElementTree as E
from project import normalize,flatten

REPO=Path(__file__).resolve().parents[3]
STANDING=REPO/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families'
NATIVE=REPO/'Documentation/Reports and Papers/Dissertation/Notes/sns_rebuild_session1_20261030/current_W2L_results/native_schematic'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def side(name):return 'Left' if name.endswith('_l') or name.startswith('L ') else 'Right' if name.endswith('_r') or name.startswith('R ') else 'Shared'
def kind(name,layer):
    if layer=='PF' and name.startswith('PF_S'):return 'PF'
    if 'RG'==layer and ('IN' in name or name.startswith('In')):return 'IN-InE' if 'E' in name or 'ext' in name else 'IN-InF'
    if layer in ('RG','PF') and not ('IN' in name):return ('HC-RG-' if layer=='RG' else 'HC-PF-')+('E' if '_E_' in name or name.endswith(' ext') or '-E_' in name else 'F')
    if layer=='PF':return 'IN-PF'
    if name.startswith('RC') or name.endswith(' RE'):return 'RC'
    if layer=='MN':return 'MN'
    if 'heel' in name.lower() or name.startswith('HEEL'):return 'SN-heel'
    if 'toe contact' in name.lower() or name.startswith('TOE_'):return 'SN-toe'
    if name.startswith('II_') or ' II' in name:return 'SN-II'
    if name.startswith('Ib_') or ' Ib' in name:return 'SN-Ib'
    if name.startswith('Ia_') or (' Ia' in name and not name.endswith(' IN')):return 'SN-Ia'
    if name.startswith('V3') or name.endswith(' V3'):return 'IN-V3'
    if 'c1' in name.lower():return 'IN-C'
    if name.startswith('DRIVE'):return 'DRIVE'
    if name=='MLR':return 'MLR'
    return 'IN-IaIN'

def arrange(nodes):
    rows=defaultdict(list)
    for n in nodes:
        t=n['type'];rank=0 if t in ('MLR','DRIVE') else 1 if t.startswith('HC-') or t=='PF' else 2 if t in ('IN-InE','IN-InF','IN-PF') else 4 if t=='MN' else 5 if t=='RC' else 3
        rows[rank].append(n)
    for rank,items in sorted(rows.items()):
        items.sort(key=lambda n:(0 if re.search(r'(^InE|_E_|-E_|ext)',n['label']) else 1,n['label']))
        for i,n in enumerate(items):n['x']=180+(i%6)*210;n['y']=90+rank*150+(i//6)*105

def hierarchy(flat):
    """Two nested levels, with true child paths for every inter-layer edge."""
    nodes=flat['nodes'];layers=defaultdict(list)
    for n in nodes:layers[n['layer']].append(copy.deepcopy(n))
    root=copy.deepcopy(flat);root['nodes']=[];root['edges']=[];paths={};homes={}
    pos={'Descending and commissural':(700,65),'RG':(700,240),'PF':(700,420),'MN':(700,760),'Reflex INs':(1060,590),'Afferents':(330,590)}
    for layer,items in layers.items():
        mid='module:'+hashlib.sha256((flat['title']+layer).encode()).hexdigest()[:12]
        graph={'nodes':[],'edges':[],'presentation':{'style':'schematic'}};subgroups=defaultdict(list)
        for n in items:subgroups[side(n['label'])].append(n)
        for k,(label,group) in enumerate(sorted(subgroups.items())):
            sid=mid+':'+label.lower();arrange(group);sub={'nodes':group,'edges':[],'presentation':{'style':'schematic'}}
            graph['nodes'].append({'id':sid,'type':'SUB','label':label+' · '+str(len(group))+' cells','scope':{'layer':layer,'side':label},'x':240+k*380,'y':180,'sub':sub,'view':{'width':260,'height':72}})
            for n in group:paths[n['id']]=[mid,sid,n['id']];homes[n['id']]=sub
        x,y=pos.get(layer,(700,1000));root['nodes'].append({'id':mid,'type':'SUB','label':layer+'\n'+str(len(items))+' cells','scope':{'layer':layer},'x':x,'y':y,'sub':graph,'view':{'width':265,'height':82}})
    for original in flat['edges']:
        a,b=paths[original['fromId']],paths[original['toId']];e=copy.deepcopy(original)
        if a[:2]==b[:2]:target=homes[original['fromId']];start=2
        elif a[0]==b[0]:target=next(n['sub'] for n in root['nodes'] if n['id']==a[0]);start=1
        else:target=root;start=0
        for key,p in [('from',a),('to',b)]:e[key+'Id']=p[start];e[key+'Path']=p[start+1:]
        target['edges'].append(e)
    root['title']=flat['title']+' — hierarchy';root['presentation']={'style':'schematic','bundle_module_wires':True}
    root=normalize(root)
    ns,es,_=flatten(root)
    assert {n['id'] for n in ns}=={n['id'] for n in flat['nodes']}
    before={e['id']:(e['fromId'],e['toId'],e['synapse']) for e in flat['edges']}
    after={e['id']:(e['source_id'].split('/')[-1],e['target_id'].split('/')[-1],e['synapse']) for e in es}
    assert before==after
    return root

def standing(family):
    path=STANDING/(family+'_graph_ledger.json');d=json.loads(path.read_text(encoding='utf-8'));nodes=[];byname={}
    for raw in d['nodes']:
        p=raw['source_neuron_parameters'];name=raw['name'];model={'kind':'persistent_NaP' if p.get('gated') else 'spiking' if p.get('spiking') else 'graded','capacitance_nF':p['membrane_capacitance']*1000,'conductance_uS':p['membrane_conductance'],'rest_mV':p['resting_potential'],'bias_nA':p['bias']}
        n={'id':raw['id'],'label':name,'type':kind(name,raw['layer']),'layer':raw['layer'],'model':model,'source_record':raw,'x':0,'y':0};nodes.append(n);byname[name]=n
    edges=[]
    for i,raw in enumerate(d['chemical_edges']+d.get('disabled_source_pairs',[])):
        a,b=byname[raw['source']],byname[raw['target']]
        edges.append({'id':raw.get('id',family+':disabled:'+str(i)),'fromId':a['id'],'toId':b['id'],'from':a['label'],'to':b['label'],'sign':'exc' if raw['reversal_mV']>b['model']['rest_mV'] else 'inh','gain':raw['gmax_uS'],'synapse':{k:raw[k] for k in ('gmax_uS','reversal_mV','threshold_mV','saturation_mV')},'source_record':raw,'pts':[]})
    arrange(nodes)
    return normalize({'title':family+' — saved standing condition','clock':{'dt_ms':2},'presentation':{'style':'schematic'},'source':{'file':str(path),'sha256':sha(path),'condition':'Historical standing export. This is not the later selected Hunt controller.'},'nodes':nodes,'edges':edges,'source_context':{k:v for k,v in d.items() if k not in ('nodes','chemical_edges','compiled_edge_matrix')},'editor_notice':'Saved standing condition. Hierarchy wires retain individual synapses. Expand a module to inspect its cells.'})

def native108():
    names=NATIVE/'NEURON_NAMES.csv';ep=NATIVE/'EXECUTED_EDGES.csv'
    source=REPO/'tmp/functional_subnetwork_remote/W2L108_verified/editable/W2L_VirtualWalker_Newer_Standalone.asim'
    raw={n.findtext('ID'):E.tostring(n,encoding='unicode') for n in E.parse(source).getroot().iter('Neuron')} if source.exists() else {}
    nodes=[]
    with names.open(encoding='utf-8-sig',newline='') as stream:records=list(csv.DictReader(stream))
    for r in records:
        name=r['NativeName'];layer='Descending and commissural' if name.endswith((' V3',' c1')) else 'RG' if ' RG ' in name else 'PF' if ' PF ' in name else 'MN' if ' MN ' in name and not name.endswith(' RE') else 'Reflex INs' if name.endswith((' IN',' RE')) else 'Afferents'
        nodes.append({'id':r['ID'],'label':name,'type':kind(name,layer),'layer':layer,'model':{'kind':'native_graded'},'source_record':r,'source_neuron_xml':raw.get(r['ID']),'x':0,'y':0})
    byid={n['id']:n for n in nodes};edges=[]
    with ep.open(encoding='utf-8-sig',newline='') as stream:edge_records=list(csv.DictReader(stream))
    for r in edge_records:
        edges.append({'id':r['source_connection_ID'],'fromId':r['source_ID'],'toId':r['target_ID'],'from':byid[r['source_ID']]['label'],'to':byid[r['target_ID']]['label'],'sign':'exc' if float(r['reversal_mV'])>-60 else 'inh','gain':float(r['gmax_uS']),'synapse':{k:float(r[k]) for k in ('gmax_uS','reversal_mV','threshold_mV','saturation_mV')},'source_record':r,'pts':[]})
    assert len(nodes)==108 and len(edges)==216
    arrange(nodes)
    return normalize({'title':'W2L108 — GENTLE108 saved native execution','nodes':nodes,'edges':edges,'clock':{'dt_ms':2},'presentation':{'style':'schematic'},'source':{'neuron_file':str(names),'neuron_sha256':sha(names),'edge_file':str(ep),'edge_sha256':sha(ep),'native_source':str(source),'native_source_sha256':sha(source) if source.exists() else None,'execution_provenance':json.loads((NATIVE/'SCHEMATIC_PROVENANCE.json').read_text())},'editor_notice':'GENTLE108 saved execution: 108 neurons, 216 synapses. Native dynamics are preserved as source data; the desktop builder does not yet implement them.'})

def reviewed_syn6():
    path=Path('D:/Users/Ben Bolen/NeuroGraphics_20261102/editor_preview/sources/Syn6_FINAL_COMPILED_GRAPH.json')
    expected='697678b3281d7fb3b2a9322c2bc77b0671e7d05beb932fcb4364b8aa5ddaddea'
    if not path.exists() or sha(path)!=expected:raise ValueError('The designated reviewed Syn6 source is missing or changed; no historical condition was substituted.')
    d=json.loads(path.read_text(encoding='utf-8'));old=json.loads((STANDING/'syn6_graph_ledger.json').read_text(encoding='utf-8'));layers={n['name']:n['layer'] for n in old['nodes']};nodes=[]
    for r in d['nodes']:
        name=r['name'];layer=layers[name];p={'kind':'persistent_NaP' if r['ionic_channels'] else 'graded','capacitance_nF':r['C_microF']*1000,'conductance_uS':r['G_uS'],'rest_mV':r['rest_mV'],'bias_nA':r['bias_nA']}
        nodes.append({'id':'syn6-reviewed:'+str(r['index']),'label':name,'type':kind(name,layer),'layer':layer,'model':p,'source_record':r,'x':0,'y':0})
    byname={n['label']:n for n in nodes};edges=[]
    for i,r in enumerate(d['chemical_edges']):
        a,b=byname[r['source']],byname[r['target']]
        edges.append({'id':'syn6-reviewed:edge:'+str(i),'fromId':a['id'],'toId':b['id'],'from':a['label'],'to':b['label'],'sign':'exc' if r['reversal_mV']>b['model']['rest_mV'] else 'inh','gain':r['gmax_uS'],'synapse':{k:r[k] for k in ('gmax_uS','reversal_mV','threshold_mV','saturation_mV')},'source_record':r,'pts':[]})
    assert len(nodes)==794 and len(edges)==3200
    arrange(nodes)
    return normalize({'title':'Syn6 — reviewed final compiled graph','nodes':nodes,'edges':edges,'clock':{'dt_ms':d['dt_s']*1000},'presentation':{'style':'schematic'},'source':{'file':str(path),'sha256':expected,'condition':'Designated reviewed final compiled graph; separate from the older standing export.'},'source_context':{k:v for k,v in d.items() if k not in ('nodes','chemical_edges')},'editor_notice':'Reviewed Syn6: 794 cells, 3,200 synapses. This is separate from the older standing condition. Native ionic dynamics are retained as source records.'})

def prepare(folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);catalog=[]
    for family in ('syn6','s3k','w2lvar','W2L108'):
        flat=reviewed_syn6() if family=='syn6' else standing(family) if family!='W2L108' else native108();tree=hierarchy(flat);stem='syn6_reviewed' if family=='syn6' else family
        for suffix,doc in [('hierarchy',tree),('exact',flat)]:
            path=folder/(stem+'_'+suffix+'.sns.json')
            if not path.exists():path.write_text(json.dumps(doc,ensure_ascii=False,indent=2),encoding='utf-8')
        catalog.append({'family':family,'hierarchy':str(folder/(stem+'_hierarchy.sns.json')),'exact':str(folder/(stem+'_exact.sns.json')),'neurons':len(flat['nodes']),'synapses':len(flat['edges']),'condition':flat['title']})
    (folder/'LIBRARY_INDEX.json').write_text(json.dumps(catalog,indent=2),encoding='utf-8')
    return catalog

if __name__=='__main__':
    import sys
    print(json.dumps(prepare(Path(sys.argv[1])),indent=2))
