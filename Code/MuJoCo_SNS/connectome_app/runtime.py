"""SNS-Toolbox consumer for saved editable projects. No physics or auto-run."""
from project import ProjectError,normalize,flatten
def build_project(document,compile_backend=True):
    from sns_toolbox.networks import Network
    from sns_toolbox.neurons import NonSpikingNeuron
    from sns_toolbox.connections import NonSpikingSynapse
    doc=normalize(document);nodes,edges,ports=flatten(doc);net=Network(name=doc.get('title','Editable connectome'))
    neural=[];lookup={}
    for n in nodes:
        model=n.get('model',{})
        if model.get('kind')=='physical' or n['type'] in {'MUSCLE','ADAPTER','SPINDLE','TENDON','PORT-load'}:continue
        if model.get('kind')!='graded':raise ProjectError(f"{n['label']}: this first builder supports explicit graded neuron models; '{model.get('kind')}' is preserved but not implemented.")
        required=('capacitance_nF','conductance_uS','rest_mV','bias_nA')
        if any(k not in model for k in required):raise ProjectError('Incomplete neuron parameters: '+n['label'])
        if model['capacitance_nF']<=0 or model['conductance_uS']<=0:raise ProjectError('C and G must be positive: '+n['label'])
        neuron=NonSpikingNeuron(membrane_capacitance=model['capacitance_nF'],membrane_conductance=model['conductance_uS'],resting_potential=model['rest_mV'],bias=model['bias_nA'])
        net.add_neuron(neuron,name=n['runtime_id']);lookup[n['runtime_id']]=n;neural.append(n)
    count=0
    for e in edges:
        if e.get('kind')=='physical':continue
        if e['source_id'] not in lookup or e['target_id'] not in lookup:raise ProjectError('A neural edge terminates on a physical block: '+e['id'])
        p=e.get('synapse',{});keys=('gmax_uS','reversal_mV','threshold_mV','saturation_mV')
        if any(p.get(k) is None for k in keys):raise ProjectError('Enter all synapse parameters before Build: '+e['id'])
        if p['gmax_uS']<0 or p['saturation_mV']<=p['threshold_mV']:raise ProjectError('Invalid conductance/voltage range: '+e['id'])
        net.add_connection(NonSpikingSynapse(max_conductance=p['gmax_uS'],reversal_potential=p['reversal_mV'],e_lo=p['threshold_mV'],e_hi=p['saturation_mV']),e['source_id'],e['target_id']);count+=1
    if not neural:raise ProjectError('The project has no buildable neural cells.')
    if not ports:ports=[{'runtime_node':n['runtime_id'],'name':n['label']} for n in neural if n.get('external_input')]
    for p in ports:
        if p['runtime_node'] not in lookup:raise ProjectError('Unresolved input port '+str(p))
        net.add_input(p['runtime_node'])
    outputs=[n for n in neural if n.get('external_output') or n['type']=='MN']
    for n in outputs:net.add_output(n['runtime_id'])
    dt=float(doc.get('clock',{}).get('dt_ms',2.0))
    if dt<=0:raise ProjectError('Neural dt_ms must be positive.')
    compiled=net.compile(dt=dt,backend='numpy') if compile_backend else None
    return {'network':net,'compiled':compiled,'neurons':len(neural),'synapses':count,'input_names':[p.get('name',p['runtime_node']) for p in ports],'output_names':[n['label'] for n in outputs],'dt_ms':dt,'physical_bindings':'Physical symbols and maps are saved; external muscle/physics binding is not executed by this neural builder.','steps_executed':0}
