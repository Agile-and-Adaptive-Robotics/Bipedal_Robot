"""Recompose nested graphs without changing canonical node or connection IDs."""
import copy
from project import normalize,flatten

def flat_document(document):
    ns,es,ports=flatten(document);out=copy.deepcopy(document)
    out['nodes']=[{k:v for k,v in n.items() if k!='runtime_id'} for n in ns]
    out['edges']=[]
    for edge in es:
        e={k:v for k,v in edge.items() if k not in ('source_id','target_id','fromPath','toPath')}
        e['fromId']=edge['source_id'].split('/')[-1];e['toId']=edge['target_id'].split('/')[-1];out['edges'].append(e)
    if ports or 'inputs' in out:
        out['inputs']=[]
        for port in ports:
            p={k:copy.deepcopy(v) for k,v in port.items() if k not in ('runtime_node','node_path')};p['node_id']=port['runtime_node'].split('/')[-1];out['inputs'].append(p)
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
    for port in out.get('inputs',[]):
        path=paths[port['node_id']];port['node_id']=path[0];port['node_path']=list(path[1:])
    for k in ('presentation','editor_views'):
        if k in layout:out[k]=copy.deepcopy(layout[k])
    return normalize(out)

def bilateral_overview(document):
    """Expose left/right layers at the overview, retaining exact leaf bindings.

    Only the generated layer/side folders are promoted. A version marker lets
    subsequent user grouping and layout remain intentional on the next open.
    """
    if document.get('presentation',{}).get('bilateral_overview_version')==1:return document
    layers={'RG','PF','MN','Afferents','Reflex INs'};layout=copy.deepcopy(document);promoted=[];changed=False
    for host in layout['nodes']:
        layer=host.get('scope',{}).get('layer',host['label'].split('\n')[0]);children=host.get('sub',{}).get('nodes',[])
        def side(node):return node.get('scope',{}).get('side',node['label'].split(' · ')[0])
        if layer not in layers or host.get('scope',{}).get('side') or not children or not all(n.get('sub') is not None and side(n) in ('Left','Right','Shared') for n in children):
            promoted.append(host);continue
        changed=True
        for child in children:
            label=side(child);child['scope']={'layer':layer,'side':label};count=len(child['sub']['nodes']);child['label']=f'{label} {layer}\n{count} cells'
            child['view']={**child.get('view',{}),'width':265,'height':82}
            child['x']=host.get('x',700)+(-230 if label=='Left' else 230 if label=='Right' else 0);child['y']=host.get('y',240)
            if layer in ('Afferents','Reflex INs'):
                child['x']=100 if label=='Left' else 1300 if label=='Right' else 700;child['y']=430 if layer=='Afferents' else 680
            promoted.append(child)
        old_group=host.get('grp');hidden=old_group and old_group in layout.get('hidden',[])
        for child in children:
            child.setdefault('grp','module:'+child['id']);layout.setdefault('groups',{})[child['grp']]=child['label']
            if hidden:layout.setdefault('hidden',[]).append(child['grp'])
    if not changed:return document
    layout['nodes']=promoted
    # Flatten the original, not the changed layout whose wires still name the
    # old folders. Recompose each wire at its new lowest common container.
    actual=flat_document(document);paths={};containers={}
    def walk(graph,prefix=()):
        containers[prefix]=graph;graph['edges']=[]
        for node in graph['nodes']:
            if node.get('sub') is not None:walk(node['sub'],(*prefix,node['id']))
            else:paths[node['id']]=(*prefix,node['id'])
    walk(layout)
    for raw in actual['edges']:
        a,b=paths[raw['fromId']],paths[raw['toId']];common=[]
        for x,y in zip(a[:-1],b[:-1]):
            if x!=y:break
            common.append(x)
        e=copy.deepcopy(raw);start=len(common)
        for key,path in [('from',a),('to',b)]:e[key+'Id']=path[start];e[key+'Path']=list(path[start+1:])
        containers[tuple(common)]['edges'].append(e)
    if actual.get('inputs'):
        # All explicit input ports retain their actual leaf neuron, too.
        for graph in containers.values():graph.pop('inputs',None)
        layout['inputs']=[]
        for raw in actual['inputs']:
            p=copy.deepcopy(raw);path=paths[p['node_id']];p['node_id']=path[0];p['node_path']=list(path[1:]);layout['inputs'].append(p)
    layout.setdefault('presentation',{})['bilateral_overview_version']=1
    return normalize(layout)
