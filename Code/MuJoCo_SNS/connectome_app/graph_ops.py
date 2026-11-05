"""Recompose nested graphs without changing canonical node or connection IDs."""
import copy
from project import normalize,flatten

def flat_document(document):
    ns,es,_=flatten(document);out=copy.deepcopy(document)
    out['nodes']=[{k:v for k,v in n.items() if k!='runtime_id'} for n in ns]
    out['edges']=[]
    for edge in es:
        e={k:v for k,v in edge.items() if k not in ('source_id','target_id','fromPath','toPath')}
        e['fromId']=edge['source_id'].split('/')[-1];e['toId']=edge['target_id'].split('/')[-1];out['edges'].append(e)
    return out

def preserve_hierarchy(actual,layout):
    """Use old ownership/layout; all model fields and endpoints come from actual."""
    actual=flat_document(actual);by={n['id']:copy.deepcopy(n) for n in actual['nodes']};paths={};containers={}
    old_edges={e['id']:e for e in flat_document(layout)['edges']}
    out=copy.deepcopy(actual)
    def graph(old,prefix=()):
        new={k:copy.deepcopy(v) for k,v in old.items() if k not in ('nodes','edges','synapses','native')};new['nodes']=[];new['edges']=[];containers[prefix]=new
        for prior in old.get('nodes',[]):
            if prior.get('sub') is not None:
                host=copy.deepcopy(prior);host['sub']=graph(prior['sub'],(*prefix,prior['id']));new['nodes'].append(host)
            elif prior['id'] in by:
                n=by.pop(prior['id'])
                for k in ('x','y','view','grp'):
                    if k in prior:n[k]=copy.deepcopy(prior[k])
                new['nodes'].append(n);paths[n['id']]=(*prefix,n['id'])
        return new
    tree=graph(layout)
    for n in by.values():tree['nodes'].append(n);paths[n['id']]=(n['id'],)
    for raw in actual['edges']:
        a,b=paths[raw['fromId']],paths[raw['toId']];common=[]
        for x,y in zip(a[:-1],b[:-1]):
            if x!=y:break
            common.append(x)
        offset=len(common);e=copy.deepcopy(raw)
        for k in ('pts','fromPort','toPort'):
            if k in old_edges.get(e['id'],{}):e[k]=copy.deepcopy(old_edges[e['id']][k])
        for side,p in [('from',a),('to',b)]:e[side+'Id']=p[offset];e[side+'Path']=list(p[offset+1:])
        containers[tuple(common)]['edges'].append(e)
    out['nodes']=tree['nodes'];out['edges']=tree['edges']
    for k in ('presentation','editor_views'):
        if k in layout:out[k]=copy.deepcopy(layout[k])
    return normalize(out)
