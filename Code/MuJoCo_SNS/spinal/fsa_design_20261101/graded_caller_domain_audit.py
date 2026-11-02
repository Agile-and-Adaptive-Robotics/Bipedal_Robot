"""Source-only graded caller-domain supplement; no neural/physics simulation.

Uses saved source graphs, actual static activation targets, and runner formulas.
Does not import MuJoCo, SNS, or a builder, and does not infer domains from traces.
"""
import argparse,array,hashlib,json,math
from pathlib import Path
import xml.etree.ElementTree as ET
HERE=Path(__file__).resolve().parent;SPINAL=HERE.parent;REPO=SPINAL.parents[2]
BASE=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031')
SHARED=BASE/'fsa_design_20261101/spiking_family_review/canonical_central/CALLER_DOMAIN_382.json'
COMP=BASE/'fsa_design/gait_full_composition_v1'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def validate(values,domains,context,expected_context):
    if context!=expected_context:raise ValueError('Changed caller parameters require a new domain ledger')
    if set(values)!=set(domains):raise ValueError('Missing/unused caller ports')
    if any(not math.isfinite(float(x)) or not domains[k][0]<=float(x)<=domains[k][1] for k,x in values.items()):raise ValueError('Caller current outside declared interval; no clipping')
    return dict(values)


def audit(family):
    reference=json.loads(SHARED.read_text());gp=REPO/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families'/(family+'_graph_ledger.json');graph=json.loads(gp.read_text());result_path=Path(graph['saved_result']);saved=json.loads(result_path.read_text());comp=json.loads((COMP/(family+'_COMPOSITION.json')).read_text());anchor=comp['saved_standing']['anchor'];anchor=anchor.get('baseline_saved_validation',anchor)
    files=[SHARED,gp,result_path,SPINAL/'runner.py',SPINAL.parent/'robot_torso/bpa_neural_coordination/gait_standing_family.py',Path(graph['actuation']['model'])]
    frozen={str(p):sha(p) for p in files}
    for item in reference['source_files']:
        if sha(item['path'])!=item['sha256']:raise ValueError('Reference caller source changed')
    if sha(result_path)!=graph['saved_result_sha256'] or sha(graph['actuation']['model'])!=graph['actuation']['model_sha256']:raise ValueError('Saved static/model source does not match graph ledger')
    AFF=anchor['actual_AFF'];MOD=anchor['actual_MOD'];G=graph['global_G'];scale=float(saved['standing_bias_scale']);affscale=float(saved['muscle_afferent_scale'])
    if AFF!=reference['caller_parameter_context']['AFF'] or MOD!=reference['caller_parameter_context']['MOD']:raise ValueError('AFF/MOD differ; do not reuse canonical bounds')
    if affscale!=1 or G.get('gamma_coact',0)!=0 or G.get('vest_prop',0)!=0 or G.get('contact_onset',0)!=0 or G.get('contra_swing',0)!=0:raise ValueError('Unsupported caller modulation or contact-event additions')
    maxima={'Ia':AFF['ia_gain']*(1+MOD['ia']),'II':AFF['i0_ii']+AFF['ii_gain']*(1+MOD['ii']),'Ib':.5*AFF['ib_gain']*1.5}
    act_names=[element.attrib['name'] for element in ET.parse(graph['actuation']['model']).getroot().find('actuator')]
    activation=saved['standing_activation_targets']
    if len(act_names)!=92 or len(activation)!=92:raise ValueError('Expected actual92-name activation mapping')
    act=dict(zip(act_names,activation));same_activation=act==reference['static_calibration']['per_muscle_activation']
    if not same_activation:raise ValueError('Static activation differs; cannot reuse reference mapping')
    inputs={r['input'] for r in graph['input_adapters']};order=[name for name in reference['input_order'] if name in inputs]
    if len(order)!=graph['census']['inputs'] or set(order)!=inputs:raise ValueError('Unmatched caller input names')
    # Reconstruct the actual archived matrix from named adapters; validates the
    # filtered reference port ordering and all weights without executing a net.
    matrix=array.array('d',[0.])*(len(graph['nodes'])*len(order));nodes={r['name']:i for i,r in enumerate(graph['nodes'])};ports={name:i for i,name in enumerate(order)}
    for row in graph['input_adapters']:matrix[nodes[row['target']]*len(order)+ports[row['input']]]=row['weight']
    matrix_hash=hashlib.sha256(matrix.tobytes()).hexdigest()
    if matrix_hash!=graph['compiled_matrix_sha256']['input_connectivity']:raise ValueError('Named adapter reconstruction does not match actual compiled matrix')
    candidate_hash=comp['saved_standing']['anchor'].get('candidate_compiled_sha256',{}).get('input_connectivity',matrix_hash)
    if candidate_hash!=matrix_hash:raise ValueError('New projection changed input mapping; separate adapter derivation required')
    domains={};recipes={}
    for name in order:
        if name.startswith(('Ia_','II_','Ib_')):domains[name]=[0.,maxima[name.split('_')[0]]];recipes[name]='Runner clipped mechanics plus actual AFF/MOD; gamma0,vest_prop0,afferent multiplier1'
        elif name.startswith('POST_'):
            current=6.*act[name[5:]]*scale;domains[name]=[current,current];recipes[name]='6*saved static activation*fixed standing_bias_scale; DRIVE0 gives stand_frac1 and ankletrim1'
        elif name=='DRIVE':domains[name]=[0.,0.];recipes[name]='Actual --stand-eval schedule,not a global descending-current bound'
        elif name=='POSTURE':domains[name]=[1.,1.];recipes[name]='Actual --stand-eval drive_posture returns POSTURE1nA,not5mV'
        elif name.startswith(('BAL_','BILAT_')):domains[name]=[0.,0.];recipes[name]='Actual --reflex-stand caller zeros supraspinal balance ports'
        elif name.startswith(('HEEL_c_','TOE_c_','LOAD_c_')):
            key='heel_rge' if name.startswith('HEEL') else 'toe_rge' if name.startswith('TOE') else 'ib_rge';gain=float(G[key]);domains[name]=[0.,1.5*gain];recipes[name]='Runner min(normal force/(.35/.50/.60*BW),1.5)*actual'+key+'; no onset/contralateral additions in this condition'
        else:raise ValueError('Unresolved caller port: '+name)
    context={'family':family,'standing_bias_scale':scale,'muscle_afferent_scale':affscale,'gamma_coact':G.get('gamma_coact',0),'vest_prop':G.get('vest_prop',0),'heel_rge':G['heel_rge'],'toe_rge':G['toe_rge'],'ib_rge':G['ib_rge'],'contact_onset':G['contact_onset'],'contra_swing':G['contra_swing'],'mode':'--stand-eval --reflex-stand'}
    endpoints=0;outside=0;base={name:bounds[0] for name,bounds in domains.items()}
    for name,bounds in domains.items():
        for x in bounds:
            values={**base,name:x}
            if validate(values,domains,context,context)!=values:raise AssertionError('Endpoint adapter altered current')
        endpoints+=1
        try:validate({**base,name:bounds[1]+1.},domains,context,context)
        except ValueError:outside+=1
    changed={**context,'standing_bias_scale':scale+.01}
    try:validate(base,domains,changed,context);changed_rejected=False
    except ValueError:changed_rejected=True
    if outside!=len(order) or not changed_rejected:raise AssertionError('Caller domain rejection contract failed')
    if frozen!={str(p):sha(p) for p in files}:raise RuntimeError('Source changed during audit')
    return {'schema':'aarl_graded_all_caller_domains_v1','role':'Gait neural pathway composition reviewer','family':family,'units':'nA','input_order':order,'domains_nA':domains,'recipes':recipes,'caller_context':context,'mathematical_muscle_maxima_nA':maxima,'reference_central_ledger':str(SHARED),'source_sha256':frozen,'static_activation_name_mapping_exact_to_reference':same_activation,'actual_compiled_input_matrix_sha256':matrix_hash,'named_input_matrix_reconstruction_exact':True,'new_syn6_projection_input_mapping_unchanged':candidate_hash==matrix_hash,'all_port_endpoints_identity_checked':endpoints,'all_port_outside_upper_rejected':outside,'changed_caller_rejected':changed_rejected,'standing_result_not_transferred_to_new_neural_candidate':True,'scope':'Actual saved standing caller only. Fixed commands/bias and conservative muscle/contact envelopes from source,not successful traces. New expression fixtures and wider commands require their own explicit contract. No simulation executed.','neural_steps':0,'physics_steps':0,'training':False,'mechanics_changed':False,'global_release':False}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    for family in ('s3k','w2lvar','syn6'):
        result=audit(family);target=args.out/(family+'_ALL_CALLER_DOMAINS.json');target.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps({'family':family,'ports':len(result['input_order']),'standing_bias_scale':result['caller_context']['standing_bias_scale'],'heel_max_nA':result['domains_nA'].get('HEEL_c_l'),'endpoints':result['all_port_endpoints_identity_checked'],'simulations':0,'output':str(target)}))
