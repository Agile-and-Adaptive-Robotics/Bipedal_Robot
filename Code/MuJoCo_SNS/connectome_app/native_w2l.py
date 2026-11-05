"""W2L AnimatLab model IO. Static XML translation; never runs a circuit."""
from __future__ import annotations
import base64,copy,hashlib,json,math,re,uuid
from pathlib import Path
from lxml import etree as E
from project import ProjectError,normalize
from graph_ops import flat_document,preserve_hierarchy

SYN={'gmax_uS':'SynAmp','reversal_mV':'Equil','threshold_mV':'ThreshV','saturation_mV':'SaturateV'}
PARAM={'rest_mV':'RestingPot','tau_ms':'TimeConst','size':'Size','tonic_nA':'TonicStimulus','gmax_ca_uS':'GMaxCa','initial_threshold_mV':'InitialThresh','noise_mV':'Noise','ca_activation_midpoint_mV':'CaActivation/MidPoint','ca_activation_slope':'CaActivation/Slope','ca_activation_tau_ms':'CaActivation/TimeConstant','ca_inactivation_midpoint_mV':'CaDeactivation/MidPoint','ca_inactivation_slope':'CaDeactivation/Slope','ca_inactivation_tau_ms':'CaDeactivation/TimeConstant'}
PROJECT_PARAM={'rest_mV':('RestingPotential','milli',1e-3),'tau_ms':('TimeConstant','milli',1e-3),'size':('RelativeSize','None',1.),'tonic_nA':('TonicStimulus','nano',1e-9),'gmax_ca_uS':('MaxCaConductance','micro',1e-6),'initial_threshold_mV':('InitialThreshold','milli',1e-3),'noise_mV':('TonicNoise','milli',1e-3),'ca_activation_midpoint_mV':('CaActivation/MidPoint','milli',1e-3),'ca_activation_slope':('CaActivation/Slope','None',1.),'ca_activation_tau_ms':('CaActivation/TimeConstant','milli',1e-3),'ca_inactivation_midpoint_mV':('CaDeactivation/MidPoint','milli',1e-3),'ca_inactivation_slope':('CaDeactivation/Slope','None',1.),'ca_inactivation_tau_ms':('CaDeactivation/TimeConstant','milli',1e-3)}
PROJECT_SYN={'gmax_uS':('MaxSynapticConductance','micro',1e-6),'reversal_mV':('EquilibriumPotential','milli',1e-3),'threshold_mV':('PreSynapticThreshold','milli',1e-3),'saturation_mV':('PreSynapticSaturationLevel','milli',1e-3)}
PARAM.update(relative_accommodation='RelativeAccom',accommodation_tau_ms='AccomTimeConst',ahp_conductance_uS='AHPAmp',ahp_tau_ms='AHPTimeConst')
PROJECT_PARAM.update(relative_accommodation=('RelativeAccomodation','None',1.),accommodation_tau_ms=('AccomodationTimeConstant','milli',1e-3),ahp_conductance_uS=('AHP_Conductance','micro',1e-6),ahp_tau_ms=('AHP_TimeConstant','milli',1e-3))
enc=lambda b:base64.b64encode(b).decode('ascii')
dec=lambda s:base64.b64decode(s)
sha=lambda b:hashlib.sha256(b).hexdigest()
def parse(raw):return E.fromstring(raw,E.XMLParser(remove_blank_text=False,strip_cdata=False,resolve_entities=False))
def xml(r):return E.tostring(r,encoding='UTF-8',xml_declaration=True)
def number(n,key,default=0.):
    v=n.findtext(key)
    return float(v) if v not in (None,'') else default
def field(n,key,value):
    x=n.find(key)
    if x is None:
        parent=n
        for part in key.split('/'):
            x=parent.find(part)
            if x is None:x=E.SubElement(parent,part)
            parent=x
    text='True' if value is True else 'False' if value is False else format(value,'.17g') if isinstance(value,(float,int)) else str(value)
    if x.text!=text:x.text=text
def scaled(n,key,value,scale='None',factor=1.):
    x=n.find(key)
    if x is None:
        field(n,key,'');x=n.find(key)
    x.attrib.update({'Value':format(float(value),'.17g'),'Scale':scale,'Actual':format(float(value)*factor,'.17g')});x.text=None
def typed_id(value):
    try:uuid.UUID(value)
    except ValueError:raise ProjectError('Native objects require GUID IDs; create the object with the current desktop editor: '+value)
    return value
def reroll(element,new_id):
    for item in element.iter('ID'):item.text=str(uuid.uuid4())
    field(element,'ID',new_id)
def remove(element):element.getparent().remove(element)
def catalog(root):return {x.findtext('ID'):x for x in root.iter('Neuron')}
def runtime_module(root):return next(x for x in root.iter('NeuralModule') if x.find('Neurons') is not None)

def read_model(path):
    path=Path(path).resolve();project=None
    if path.suffix.lower()=='.aproj':
        project=path;pr=parse(path.read_bytes());asim=path.parent/pr.findtext('SimulationFile','')
        if not asim.is_file():asim=path.with_suffix('.asim')
        if not asim.is_file():raise ProjectError('The native project needs its paired .asim export in the same folder.')
    else:
        asim=path;candidates=[path.with_suffix('.aproj'),path.with_name(path.stem.removesuffix('_Standalone')+'.aproj')]
        project=next((p for p in candidates if p.exists()),None)
    raw=asim.read_bytes();r=parse(raw);module=runtime_module(r);real=catalog(r)
    if len(real)<1:raise ProjectError('No native neurons were found.')
    from controller_library import kind,arrange,hierarchy
    ns=[]
    for ident,n in real.items():
        name=n.findtext('Name');layer='Descending and commissural' if name.endswith((' V3',' c1')) else 'RG' if ' RG ' in name else 'PF' if ' PF ' in name else 'MN' if ' MN ' in name and not name.endswith(' RE') else 'Reflex INs' if name.endswith((' IN',' RE')) else 'Afferents'
        values={k:number(n,p) for k,p in PARAM.items()}
        ns.append({'id':ident,'type':kind(name,layer),'label':name,'layer':layer,'model':{'kind':'native_graded','rest_mV':values['rest_mV'],'native_parameters':values,'enabled':n.findtext('Enabled')!='False','burst_init_at_bottom':n.findtext('BurstInitAtBottom')=='True'},'source_record':{'native_id':ident,'native_name':name},'x':0,'y':0})
    ts={x.findtext('ID'):x for x in module.findall('Synapses/NonSpikingSynapses/SynapseType')};edges=[]
    for e in module.iter('Connexion'):
        t=ts[e.findtext('SynapseTypeID')];target=e.findtext('TargetID');p={k:number(t,v) for k,v in SYN.items()}
        edges.append({'id':e.findtext('ID'),'fromId':e.findtext('SourceID'),'toId':target,'from':real[e.findtext('SourceID')].findtext('Name'),'to':real[target].findtext('Name'),'synapse':p,'gain':p['gmax_uS'],'sign':'exc' if p['reversal_mV']>number(real[target],'RestingPot') else 'inh','native_type_id':t.findtext('ID'),'native_type_name':t.findtext('Name'),'pts':[]})
    adapters=[];bodies={}
    for b in r.iter('RigidBody'):
        if b.findtext('ID'):bodies[b.findtext('ID')]={'id':b.findtext('ID'),'name':b.findtext('Name'),'type':b.findtext('Type'),'fields':{c.tag:c.text for c in b if len(c)==0 and c.text is not None}}
    for a in r.iter('Adapter'):
        adapters.append({'id':a.findtext('ID'),'name':a.findtext('Name'),'direction':a.findtext('Type'),'source_id':a.findtext('SourceID'),'target_id':a.findtext('TargetID'),'source_data_type':a.findtext('SourceDataType'),'target_data_type':a.findtext('TargetDataType'),'enabled':a.findtext('Enabled')!='False','gain':{k:number(a.find('Gain'),k) for k in 'ABCD'}})
    native={'format':'AnimatLab-asim','asim_name':asim.name,'asim_b64':enc(raw),'asim_sha256':sha(raw),'project_name':project.name if project else None,'project_b64':enc(project.read_bytes()) if project else None,'project_sha256':sha(project.read_bytes()) if project else None,'primary_suffix':path.suffix.lower(),'adapters':adapters,'bodies':bodies,'source_files':{str(asim):sha(raw),**({str(project):sha(project.read_bytes())} if project else {})}}
    arrange(ns);flat=normalize({'title':path.stem,'nodes':ns,'edges':edges,'clock':{'dt_ms':number(module,'TimeStep')*1000},'native':native,'source':{'file':str(path),'condition':'Native W2L file; no external tuning/runtime overlays were applied.'},'editor_notice':'Native model: Save writes the real .asim and paired .aproj. Save As creates another model copy.'})
    doc=hierarchy(flat)
    sidecar=path.with_suffix('.editor.json')
    if not sidecar.exists() and project:sidecar=project.with_suffix('.editor.json')
    if sidecar.exists():
        old=json.loads(sidecar.read_text(encoding='utf-8'));doc=preserve_hierarchy(flat,old['document'])
    from graph_ops import bilateral_overview
    return bilateral_overview(doc)

def model_neurons(document):
    d=flat_document(document);return d,{n['id']:n for n in d['nodes'] if n.get('model',{}).get('kind')!='physical' and n['type'] not in ('MUSCLE','SPINDLE','TENDON','ADAPTER','PORT-load')}

def build_bytes(document):
    """Prepare both native formats entirely in memory, before any file write."""
    doc,neurons=model_neurons(document);native=doc.get('native',{})
    if native.get('format')!='AnimatLab-asim':raise ProjectError('This project has no native W2L source binding.')
    raw=dec(native['asim_b64']);r=parse(raw);before=E.tostring(r);module=runtime_module(r);old=catalog(r);old_edges={e.findtext('ID'):e for e in module.iter('Connexion')};types={x.findtext('ID'):x for x in module.findall('Synapses/NonSpikingSynapses/SynapseType')}
    template=next(iter(old.values()));newtypes={};deleted=set(old)-set(neurons)
    for ident,n in neurons.items():
        typed_id(ident)
        if ident not in old:
            original=n.get('source_record',{}).get('native_id');x=copy.deepcopy(old.get(original,template));reroll(x,ident);module.find('Neurons').append(x);old[ident]=x
        x=old[ident];p=n.get('model',{});values=p.get('native_parameters')
        if values is None:
            if p.get('kind')!='graded':raise ProjectError('Specify native neuron parameters before saving '+n['label'])
            values={k:number(x,v) for k,v in PARAM.items()};values.update(rest_mV=p['rest_mV'],tau_ms=p['capacitance_nF']/p['conductance_uS'],size=p['conductance_uS'],tonic_nA=p['bias_nA'],gmax_ca_uS=0.)
        if any(not isinstance(v,(float,int)) or not math.isfinite(v) for v in values.values()):raise ProjectError('Native neuron values must be finite numbers.')
        if values['tau_ms']<=0 or values['size']<=0:raise ProjectError('Neuron time constant and size must be positive.')
        field(x,'Name',n['label']);field(x,'Enabled',p.get('enabled',True));field(x,'BurstInitAtBottom',p.get('burst_init_at_bottom',True))
        for k,key in PARAM.items():
            if k in values and number(x,key)!=values[k]:field(x,key,values[k])
    for ident in deleted:remove(old[ident])
    wanted={};conn_template=next(iter(old_edges.values()))
    for e in doc['edges']:
        if e.get('kind')=='physical':continue
        typed_id(e['id']);a,b=e['fromId'],e['toId']
        if a not in neurons or b not in neurons:raise ProjectError('A native synapse needs two neural endpoints.')
        p=e.get('synapse',{})
        if any(not isinstance(p.get(k),(int,float)) or not math.isfinite(p[k]) for k in SYN):raise ProjectError('Complete all four native synapse values before saving.')
        if p['gmax_uS']<0 or p['saturation_mV']<=p['threshold_mV']:raise ProjectError('Invalid native synapse conductance or voltage range.')
        x=old_edges.get(e['id'])
        if x is None:x=copy.deepcopy(conn_template);field(x,'ID',e['id']);conn_template.getparent().append(x)
        tid=x.findtext('SynapseTypeID');base=types.get(tid,next(iter(types.values())))
        if any(number(base,v)!=p[k] for k,v in SYN.items()):
            # Clone for this connection: a shared named type must not change
            # other connections when an individual wire is edited.
            tid=str(uuid.uuid5(uuid.NAMESPACE_URL,'sns-editor:'+e['id']+json.dumps(p,sort_keys=True)))
            if tid not in types:
                t=copy.deepcopy(base);field(t,'ID',tid);field(t,'Name',base.findtext('Name')+' [editor '+e['id'][:8]+']')
                for k,v in SYN.items():field(t,v,p[k])
                base.getparent().append(t);types[tid]=t;newtypes[tid]={'base':base.findtext('ID'),'values':p}
        field(x,'SourceID',a);field(x,'TargetID',b);field(x,'SynapseTypeID',tid)
        if e['id'] not in old_edges or any(number(base,v)!=p[k] for k,v in SYN.items()):field(x,'G',p['gmax_uS'])
        wanted[e['id']]=x
    for ident,x in old_edges.items():
        if ident not in wanted:remove(x)
    adapter_records={a['id']:a for a in native.get('adapters',[])}
    for x in list(r.iter('Adapter')):
        ident=x.findtext('ID');a=adapter_records.get(ident)
        if a is None or a['source_id'] in deleted or a['target_id'] in deleted:remove(x);continue
        for key,fieldname in [('source_id','SourceID'),('target_id','TargetID'),('source_data_type','SourceDataType'),('target_data_type','TargetDataType'),('enabled','Enabled')]:
            value=a[key];oldvalue=x.findtext(fieldname)
            if (value if isinstance(value,str) else str(value))!=oldvalue:field(x,fieldname,value)
        for k,v in a['gain'].items():
            if not math.isfinite(v):raise ProjectError('Adapter polynomial values must be finite.')
            if number(x.find('Gain'),k)!=v:field(x.find('Gain'),k,v)
    for x in list(r.iter('DataColumn')):
        if x.findtext('TargetID') in deleted:remove(x)
    for x in list(r.findall('.//ExternalStimuli/Stimulus')):
        if x.findtext('TargetNodeID') in deleted:remove(x)
    out=raw if E.tostring(r)==before else xml(r)
    project=update_project(dec(native['project_b64']),neurons,doc['edges'],newtypes,adapter_records,deleted,r,parse(raw)) if native.get('project_b64') else None
    verify_native(r)
    return out,project

def verify_native(root):
    nodes=catalog(root);types={x.findtext('ID') for x in root.iter('SynapseType')};ids=[x.text for x in root.iter('ID') if x.text]
    # ID reference elements elsewhere are not object declarations.
    declarations=[x.findtext('ID') for x in root.iter() if x.tag in ('Neuron','Connexion','SynapseType','Adapter')]
    if len(declarations)!=len(set(declarations)):raise ProjectError('Duplicate native object IDs.')
    for e in root.iter('Connexion'):
        if e.findtext('SourceID') not in nodes or e.findtext('TargetID') not in nodes or e.findtext('SynapseTypeID') not in types:raise ProjectError('Dangling native synapse endpoint/type.')
    return True

def update_project(raw,neurons,edges,newtypes,adapters,deleted,asim,baseline):
    """Keep native GUI objects, type triplets, page links, and drawings in sync."""
    root=parse(raw);before=E.tostring(root);allnodes={x.findtext('ID'):x for x in root.iter('Node') if x.findtext('ID')};real={i:x for i,x in allnodes.items() if x.findtext('ClassName','').endswith('Neurons.NonSpiking')};links={x.findtext('ID'):x for x in root.iter('Link') if x.findtext('ClassName','').endswith('Behavior.Synapse')};typebank={x.findtext('ID'):x for x in root.iter('Link') if '.SynapseTypes.' in x.findtext('ClassName','')};bodymap={x.findtext('ID'):x for x in root.iter('RigidBody') if x.findtext('ID')}
    original_neurons=catalog(baseline);original_edges={x.findtext('ID'):x for x in baseline.iter('Connexion')};original_types={x.findtext('ID'):x for x in baseline.iter('SynapseType')};original_adapters={x.findtext('ID'):x for x in baseline.iter('Adapter')}
    topology=bool(deleted or set(neurons)!=set(real) or {e['id'] for e in edges if e.get('kind')!='physical'}!=set(links))
    base_node=next(iter(real.values()));base_link=next(iter(links.values()));off_template=next(x for x in allnodes.values() if x.findtext('ClassName','').endswith('.OffPage'))
    def page(node):return next(x for x in node.iterancestors() if x.find('Nodes') is not None and x.find('Links') is not None)
    for ident,n in neurons.items():
        x=real.get(ident)
        if x is None:
            x=copy.deepcopy(real.get(n.get('source_record',{}).get('native_id'),base_node));reroll(x,ident)
            for k in ('InLinks','OutLinks'):x.find(k).clear()
            base_node.getparent().append(x);real[ident]=x;allnodes[ident]=x
            for gate in ('CaActivation','CaDeactivation'):
                field(x,gate+'/ID',catalog(asim)[ident].findtext(gate+'/ID'))
        p=n['model'];values=p.get('native_parameters')
        if values is None:values={k:number(catalog(asim)[ident],key) for k,key in PARAM.items()}
        old_neuron=original_neurons.get(ident)
        if old_neuron is None or old_neuron.findtext('Name')!=n['label']:field(x,'Text',n['label'])
        if old_neuron is None or (old_neuron.findtext('Enabled')!='False')!=p.get('enabled',True):field(x,'Enabled',p.get('enabled',True))
        if old_neuron is None or (old_neuron.findtext('BurstInitAtBottom')=='True')!=p.get('burst_init_at_bottom',True):field(x,'InitAtBottom',p.get('burst_init_at_bottom',True))
        for k,(key,scale,factor) in PROJECT_PARAM.items():
            old=x.find(key);actual=float(old.get('Actual','nan')) if old is not None else float('nan')
            if k in values and (old_neuron is None or number(old_neuron,PARAM[k])!=values[k]):scaled(x,key,values[k],scale,factor)
    for tid,row in newtypes.items():
        if tid in typebank:continue
        t=copy.deepcopy(typebank[row['base']]);field(t,'ID',tid);field(t,'Name',t.findtext('Name')+' [editor]')
        for k,(key,scale,factor) in PROJECT_SYN.items():scaled(t,key,row['values'][k],scale,factor)
        typebank[row['base']].getparent().append(t);typebank[tid]=t
    def represented(owner,ident):
        if real[ident].getparent() is owner.find('Nodes'):return ident
        existing=next((n for n in owner.find('Nodes') if n.findtext('LinkedNodeID')==ident),None)
        if existing is not None:return existing.findtext('ID')
        proxy=copy.deepcopy(off_template);pid=str(uuid.uuid5(uuid.NAMESPACE_URL,'sns-editor:offpage:'+owner.findtext('ID','root')+ident));reroll(proxy,pid);field(proxy,'LinkedNodeID',ident);field(proxy,'Text',neurons[ident]['label'])
        for k in ('InLinks','OutLinks'):proxy.find(k).clear()
        owner.find('Nodes').append(proxy);allnodes[pid]=proxy;return pid
    wanted={}
    asim_edges={x.findtext('ID'):x for x in asim.iter('Connexion')}
    for e in edges:
        if e.get('kind')=='physical':continue
        ident=e['id'];x=links.get(ident)
        if x is None:x=copy.deepcopy(base_link);field(x,'ID',ident);base_link.getparent().append(x)
        owner=page(x);previous=original_edges.get(ident)
        for side,key,nativekey in [('from','OriginID','SourceID'),('to','DestinationID','TargetID')]:
            if previous is None or previous.findtext(nativekey)!=e[side+'Id']:
                field(x,key,represented(owner,e[side+'Id']));topology=True
        if previous is None or previous.findtext('SynapseTypeID')!=asim_edges[ident].findtext('SynapseTypeID'):field(x,'SynapticTypeID',asim_edges[ident].findtext('SynapseTypeID'))
        old=x.find('SynapticConductance');value=float(old.get('Actual')) if old is not None else float('nan')
        if previous is None or number(original_types[previous.findtext('SynapseTypeID')],'SynAmp')!=e['synapse']['gmax_uS']:scaled(x,'SynapticConductance',e['synapse']['gmax_uS'],'micro',1e-6)
        wanted[ident]=x
    for ident,x in links.items():
        if ident not in wanted:remove(x)
    for ident,x in list(allnodes.items()):
        if ident in deleted or x.findtext('LinkedNodeID') in deleted:remove(x);allnodes.pop(ident,None)
        elif x.findtext('LinkedNodeID') in neurons:
            target=x.findtext('LinkedNodeID');prior=original_neurons.get(target)
            if prior is None or prior.findtext('Name')!=neurons[target]['label']:field(x,'Text',neurons[target]['label'])
    for ident in set(original_adapters)-set(adapters):
        x=allnodes.pop(ident,None)
        if x is not None:
            for link in list(root.iter('Link')):
                if ident in (link.findtext('OriginID'),link.findtext('DestinationID')):remove(link)
            remove(x);topology=True
    for ident,a in adapters.items():
        x=allnodes.get(ident)
        if x is None:continue
        if a['source_id'] in deleted or a['target_id'] in deleted:
            # Remove the adapter and its two GUI wires, not unrelated mechanics.
            for l in list(root.iter('Link')):
                if l.findtext('OriginID')==ident or l.findtext('DestinationID')==ident:remove(l)
            remove(x);allnodes.pop(ident,None);continue
        previous=original_adapters.get(ident)
        if previous is None or (previous.findtext('Enabled')!='False')!=a['enabled']:field(x,'Enabled',a['enabled'])
        for k,v in a['gain'].items():
            old=x.find('Gain/'+k);actual=float(old.get('Actual','nan')) if old is not None else float('nan')
            if previous is None or number(previous.find('Gain'),k)!=v:scaled(x,'Gain/'+k,v)
    # Rebuild only logical synapse registrations; retain physical/adapter links.
    synids=(set(links)-set(wanted))|{i for i,x in wanted.items() if i not in original_edges or x.findtext('OriginID')!=links.get(i,x).findtext('OriginID') or original_edges[i].findtext('SourceID')!=asim_edges[i].findtext('SourceID') or original_edges[i].findtext('TargetID')!=asim_edges[i].findtext('TargetID')}
    for x in allnodes.values():
        for key in ('InLinks','OutLinks'):
            block=x.find(key)
            if block is not None:
                for child in list(block):
                    if child.text in synids:
                        logical=wanted.get(child.text)
                        endpoint='DestinationID' if key=='InLinks' else 'OriginID'
                        if logical is None or logical.findtext(endpoint)!=x.findtext('ID'):block.remove(child)
    for ident,x in wanted.items():
        if ident not in synids:continue
        for key,endpoint in [('OutLinks','OriginID'),('InLinks','DestinationID')]:
            n=allnodes[x.findtext(endpoint)];block=n.find(key)
            if block is None:block=E.SubElement(n,key)
            if not any(c.text==ident for c in block):E.SubElement(block,'ID').text=ident
    surviving_links={x.findtext('ID') for x in root.iter('Link')}
    for x in allnodes.values():
        for key in ('InLinks','OutLinks'):
            block=x.find(key)
            if block is not None:
                for ref in list(block):
                    if ref.text not in surviving_links:block.remove(ref)
    if topology:synchronize_drawings(root)
    for x in root.iter('DiagramXml'):
        if x.text:
            drawing=E.fromstring(x.text.encode());changed=False
            for item in drawing.findall('.//AddFlow/Node'):
                ident=item.findtext('Tag');node=allnodes.get(ident)
                canonical=node.findtext('LinkedNodeID') or ident if node is not None else None
                old=original_neurons.get(canonical);new=neurons.get(canonical)
                if new is not None and (old is None or old.findtext('Name')!=new['label']):field(item,'Text',new['label']);changed=True
            if changed:x.text=E.CDATA(E.tostring(drawing,encoding='unicode'))
    return raw if E.tostring(root)==before else xml(root)

def synchronize_drawings(root):
    """Keep AddFlow interleaved node/link indexes valid after topology edits."""
    for owner in root.iter():
        payload=owner.find('DiagramXml');nodes=owner.find('Nodes');links=owner.find('Links')
        if payload is None or nodes is None or links is None or not payload.text:continue
        full=E.fromstring(payload.text.encode());drawing=full.find('.//AddFlow')
        if drawing is None:continue
        old=E.tostring(full);by={n.findtext('ID'):n for n in nodes};linkmap={l.findtext('ID'):l for l in links}
        entries=[x for x in drawing if x.tag in ('Node','Link')]
        for x in entries:
            tag=x.findtext('Tag')
            if (x.tag=='Node' and tag not in by) or (x.tag=='Link' and tag not in linkmap):drawing.remove(x)
        drawn={x.findtext('Tag') for x in drawing if x.tag=='Node'}
        node_template=next((x for x in drawing if x.tag=='Node'),None);link_template=next((x for x in drawing if x.tag=='Link'),None)
        for ident,n in by.items():
            if ident in drawn:continue
            loc=n.find('Location');size=n.find('Size');x=copy.deepcopy(node_template) if node_template is not None else E.Element('Node');drawing.append(x)
            x.attrib.update({'Left':loc.get('x','100') if loc is not None else '100','Top':loc.get('y','100') if loc is not None else '100','Width':size.get('Width','100') if size is not None else '100','Height':size.get('Height','55') if size is not None else '55'});field(x,'Tag',ident);field(x,'Text',n.findtext('Text',''))
        drawnlinks={x.findtext('Tag') for x in drawing if x.tag=='Link'}
        for ident,l in linkmap.items():
            if ident not in drawnlinks:
                x=copy.deepcopy(link_template) if link_template is not None else E.Element('Link');drawing.append(x);field(x,'Tag',ident)
        indexes={x.findtext('Tag'):i for i,x in enumerate(c for c in drawing if c.tag in ('Node','Link')) if x.tag=='Node'}
        for x in drawing:
            if x.tag!='Link':continue
            link=linkmap[x.findtext('Tag')]
            if link.findtext('OriginID') not in indexes or link.findtext('DestinationID') not in indexes:raise ProjectError('An AddFlow link has no local endpoint representation.')
            x.set('Org',str(indexes[link.findtext('OriginID')]));x.set('Dst',str(indexes[link.findtext('DestinationID')]))
        drawing.set('Nodes',str(sum(x.tag=='Node' for x in drawing)));drawing.set('Links',str(sum(x.tag=='Link' for x in drawing)))
        if E.tostring(full)!=old:payload.text=E.CDATA(E.tostring(full,encoding='unicode'))
