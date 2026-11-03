"""Lossless editable connectome documents and atomic file sessions."""
from __future__ import annotations
import copy,hashlib,json,math,os,tempfile,uuid
from pathlib import Path

SCHEMA='sns-connectome-project/1'
PHYSICAL={'MUSCLE','ADAPTER','SPINDLE','TENDON','PORT-load'}
class ProjectError(ValueError):pass
def clone(x):return copy.deepcopy(x)
def fresh(prefix):return prefix+'_'+uuid.uuid4().hex
def fingerprint(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def node_model(kind='graded'):
    return {'kind':kind,'capacitance_nF':30.0,'conductance_uS':1.0,'rest_mV':0.0,'bias_nA':0.0}
def normalize(document):
    """Preserve arbitrary fields; normalize only editable graph identities."""
    doc=clone(document)
    if not isinstance(doc,dict) or not isinstance(doc.get('nodes'),list):
        raise ProjectError('Open an editor project with nodes and edges. Compiled-source ledgers require a source-specific adapter.')
    if 'chemical_edges' in doc and 'edges' not in doc:
        raise ProjectError('This is a compiled-source ledger, not an editor project. Use a source-specific import adapter.')
    def graph(g):
        g.setdefault('edges',g.get('synapses',[]))
        for n in g['nodes']:
            n.setdefault('id',fresh('n'));n.setdefault('label',n.get('name',n['id']));n.setdefault('x',0);n.setdefault('y',0)
            if 'type' not in n:raise ProjectError('A node lacks an editor block type. Source ledgers need an explicit adapter.')
            if n.get('sub') is not None:graph(n['sub'])
            elif n['type'] not in PHYSICAL:n.setdefault('model',{'kind':'unresolved'})
        byid={n['id']:n for n in g['nodes']};bylabel={}
        for n in g['nodes']:bylabel.setdefault(n['label'],[]).append(n)
        for e in g['edges']:
            e.setdefault('id',fresh('e'))
            for side in ('from','to'):
                key=side+'Id'
                if key not in e:
                    if e.get(side) in byid:e[key]=e[side]
                    elif len(bylabel.get(e.get(side),[]))==1:e[key]=bylabel[e[side]][0]['id']
                    else:raise ProjectError(f'Unresolved or ambiguous endpoint: {e.get(side)}')
                if e[key] not in byid:raise ProjectError('Missing endpoint '+str(e[key]))
                e[side]=byid[e[key]]['label']
            e.setdefault('pts',[])
            if e.get('kind')!='physical' and 'synapse' not in e:
                e['synapse']={'gmax_uS':e.get('gain',0.5),'reversal_mV':None,'threshold_mV':None,'saturation_mV':None}
    graph(doc);doc.setdefault('schema',SCHEMA);doc.setdefault('clock',{'dt_ms':2.0})
    validate(doc)
    return doc
def validate(doc):
    def finite(x):
        if isinstance(x,float) and not math.isfinite(x):raise ProjectError('Non-finite parameter')
        if isinstance(x,dict):
            for v in x.values():finite(v)
        elif isinstance(x,list):
            for v in x:finite(v)
    finite(doc)
    def graph(g):
        ids=[n['id'] for n in g['nodes']]
        if len(ids)!=len(set(ids)):raise ProjectError('Duplicate node IDs in a subsystem')
        edges=g.get('edges',[]);eids=[e['id'] for e in edges]
        if len(eids)!=len(set(eids)):raise ProjectError('Duplicate edge IDs')
        by={n['id']:n for n in g['nodes']}
        def endpoint(e,side):
            key=e[side+'Id']
            if key not in by:raise ProjectError('Dangling edge '+e['id'])
            n=by[key]
            for child in e.get(side+'Path',[]):
                if not n.get('sub'):raise ProjectError('Boundary path enters a non-subsystem')
                n=next((x for x in n['sub']['nodes'] if x['id']==child),None)
                if n is None:raise ProjectError('Missing child boundary endpoint '+str(child))
        for e in edges:
            endpoint(e,'from');endpoint(e,'to')
        for n in g['nodes']:
            if n.get('sub') is not None:graph(n['sub'])
    graph(doc);return True
def flatten(doc):
    """Return actual neurons and edges behind nested display modules."""
    doc=normalize(doc);ns=[];es=[];ports=[]
    def graph(g,prefix=()):
        by={n['id']:n for n in g['nodes']}
        def resolve(e,side):
            path=[e[side+'Id'],*e.get(side+'Path',[])];n=by[path[0]]
            for child in path[1:]:n=next(x for x in n['sub']['nodes'] if x['id']==child)
            if n.get('sub') is not None:raise ProjectError('Bind the boundary wire to a child neuron before Build: '+e['id'])
            return '/'.join((*prefix,*path))
        for n in g['nodes']:
            if n.get('sub') is not None:graph(n['sub'],(*prefix,n['id']))
            else:ns.append({**clone(n),'runtime_id':'/'.join((*prefix,n['id']))})
        for e in g.get('edges',[]):es.append({**clone(e),'source_id':resolve(e,'from'),'target_id':resolve(e,'to')})
        for p in g.get('inputs',[]):ports.append({**clone(p),'runtime_node':'/'.join((*prefix,p['node_id'],*p.get('node_path',[])))})
    graph(doc);return ns,es,ports

class FileSession:
    def __init__(self):self.paths={};self.hashes={}
    def open(self,path,tab):
        path=Path(path).resolve();raw=json.loads(path.read_text(encoding='utf-8-sig'));doc=normalize(raw)
        self.paths[tab]=path;self.hashes[tab]=fingerprint(path)
        return doc
    def save(self,doc,tab,path=None):
        doc=normalize(doc);destination=Path(path).resolve() if path else self.paths.get(tab)
        if destination is None:raise ProjectError('Use Save As to choose a project file.')
        if destination==self.paths.get(tab) and destination.exists() and fingerprint(destination)!=self.hashes.get(tab):
            raise ProjectError('This file changed outside the editor. Reopen it or use Save As to preserve both versions.')
        destination.parent.mkdir(parents=True,exist_ok=True)
        data=json.dumps(doc,ensure_ascii=False,indent=2,allow_nan=False).encode('utf-8')
        fd,temp=tempfile.mkstemp(prefix='.'+destination.name+'.',suffix='.tmp',dir=destination.parent)
        try:
            with os.fdopen(fd,'wb') as stream:stream.write(data);stream.flush();os.fsync(stream.fileno())
            os.replace(temp,destination)
        finally:
            if os.path.exists(temp):os.unlink(temp)
        self.paths[tab]=destination;self.hashes[tab]=fingerprint(destination)
        return destination,doc
