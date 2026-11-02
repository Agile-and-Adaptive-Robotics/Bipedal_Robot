"""Default-off A/B heel-routing hypotheses; static compiled-current proofs only.

A: ipsi InE excitation + ipsi InF inhibition (selective ipsilateral extension).
B: ipsi InE excitation + contra InF excitation (W2L laterality precedent).
Neither adds direct heel->contra RG-F excitation or changes caller encoding.
No neuronal or body time integration is used by the CLI proof.
"""
from __future__ import annotations
import argparse,copy,hashlib,json,os,sys
from pathlib import Path
if __name__=='__main__':
    for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
HERE=Path(__file__).resolve().parent


def utilities():
    try:from . import graded_profile as G
    except ImportError:import graded_profile as G
    return G


def make_policy(family,hypothesis):
    return {'schema':'aarl_optional_heel_stance_v1','family':family,'hypothesis':hypothesis,
            'conductance_policy':'canonical_right_existing_ipsi_InE',
            'event_contract':'Upstream caller must declare pulse or held contact; this policy neither detects onset nor changes source current.'}


def apply_policy(spinal,policy=None):
    if policy is None:return {'applied':False}
    G=utilities();G.assert_campaign_allowed('optional_profile_implementation')
    p=json.loads(Path(policy).read_text()) if isinstance(policy,(str,Path)) else copy.deepcopy(policy)
    G.exact_keys(p,('schema','family','hypothesis','conductance_policy','event_contract'),'heel policy')
    if p['schema']!='aarl_optional_heel_stance_v1' or p['hypothesis'] not in ('A_ipsilateral_extension','B_contralateral_extension_inhibition') or p['conductance_policy']!='canonical_right_existing_ipsi_InE' or not p['event_contract'].strip():raise ValueError('Unsupported/unbound heel hypothesis')
    actual={'build_network':'s3k','build_network_w2lvar':'w2lvar','build_network_syn6':'syn6'}.get(type(spinal).__module__.split('.')[-1])
    if actual!=p['family']:raise ValueError('Actual source family does not match heel policy')
    c=G.raw_compiled(spinal)
    if spinal.compiled is not c or c.spiking or c.electrical or c.electrical_rectified:raise ValueError('Apply to graded dense chemical SNS before constructing an accelerator')
    names=G.population_names(spinal);ix={n:i for i,n in enumerate(names)}
    required=[f'{kind}_{s}' for kind in ('HEEL','InE','InF','RG_E','RG_F') for s in ('r','l')]
    if any(n not in ix for n in required):raise ValueError('Contact-enabled design fixture required; no missing cells invented')
    for s in ('r','l'):
        for pre,post in ((f'InE_{s}',f'RG_F_{s}'),(f'InF_{s}',f'RG_E_{s}')):
            i,j=ix[post],ix[pre]
            if c.g_max_non[i,j]<=0 or c.del_e[i,j]>=0:raise ValueError('Expected inhibitory RG relay not present')
        for phase in ('E','F'):
            for t in ('r','l'):
                if c.g_max_non[ix[f'RG_{phase}_{t}'],ix[f'HEEL_{s}']]:raise ValueError('Existing direct heel->RG route needs a separate declared composition; no silent stacking')
    i,j=ix['InE_r'],ix['HEEL_r'];g=float(c.g_max_non[i,j]);exc=float(c.del_e[i,j]);lo=float(c.e_lo[i,j]);hi=float(c.e_hi[i,j]);inh=float(c.del_e[ix['RG_E_r'],ix['InF_r']])
    if g<=0 or exc<=0 or hi<=lo or inh>=0:raise ValueError('Missing actual source reference synapse')
    before=G.field_digest(c);state=c.V.copy();old={};changes=[];targets={};before_fields={k:getattr(c,k).copy() for k in ('g_max_non','del_e','e_lo','e_hi')}
    for s in ('r','l'):
        t='l' if s=='r' else 'r';src=ix[f'HEEL_{s}']
        chosen={f'InE_{s}':exc,(f'InF_{s}' if p['hypothesis']=='A_ipsilateral_extension' else f'InF_{t}'):(inh if p['hypothesis']=='A_ipsilateral_extension' else exc)}
        targets[f'HEEL_{s}']=chosen
        for side in ('r','l'):
            for kind in ('InE','InF'):
                name=f'{kind}_{side}';dst=ix[name];old[(dst,src)]={k:float(getattr(c,k)[dst,src]) for k in ('g_max_non','del_e','e_lo','e_hi')}
                c.g_max_non[dst,src]=g if name in chosen else 0.
                if name in chosen:
                    c.del_e[dst,src]=chosen[name];c.e_lo[dst,src]=lo;c.e_hi[dst,src]=hi
                new={k:float(getattr(c,k)[dst,src]) for k in old[(dst,src)]}
                if new!=old[(dst,src)]:changes.append({'source':f'HEEL_{s}','target':name,'before':old[(dst,src)],'after':new})
    if not np.array_equal(c.V,state):raise AssertionError('Dynamic states changed')
    for key,value in before_fields.items():
        if any((int(a),int(b)) not in old for a,b in zip(*np.nonzero(value!=getattr(c,key)))):raise AssertionError('Fields changed outside declared heel->RG-interneuron scope')
    scoped=[{'target_index':i,'source_index':j,**{k:float(getattr(c,k)[i,j]) for k in ('g_max_non','del_e','e_lo','e_hi')}} for i,j in old]
    record={'applied':True,'policy':p,'actual_reference':{'conductance_uS':g,'exc_reversal_mV':exc,'inh_reversal_mV':inh,'window_mV':[lo,hi]},'graph_delta':changes,'scoped_final_fields':scoped,'outside_scope_chemical_fields_bitexact':True,'before_parameter_sha256':before,'after_parameter_sha256':G.field_digest(c),'dynamic_states_changed':False,'targets':targets,'source_current_encoding_unchanged':True,'direct_contra_RG_F_added':False,'other_heel_PF_paths_preserved':True,'scope':'Optional engineering hypothesis,not a unique physiological correction or a walking result'}
    spinal._heel_stance_policy_record=record
    return record


def verify_policy(spinal):
    record=getattr(spinal,'_heel_stance_policy_record',None)
    if record is None:return {'applied':False}
    c=utilities().raw_compiled(spinal)
    for item in record['scoped_final_fields']:
        i,j=item['target_index'],item['source_index']
        if any(float(getattr(c,k)[i,j])!=item[k] for k in ('g_max_non','del_e','e_lo','e_hi')):raise ValueError('Final consumed heel routing differs from the named hypothesis')
    return {'applied':True,'scoped_fields_match':True,'hypothesis':record['policy']['hypothesis'],'sparse_cache_parity_tested':False}


def isolated_static_operator(spinal,source_current_nA=.75,operating_mV=2.5):
    """Exact clipped-conductance steady relay operator from actual SNS arrays.

    Other incoming paths are excluded explicitly; IN operating bias is supplied
    only in this analytical fixture. This is not the network's live equilibrium.
    """
    G=utilities();c=G.raw_compiled(spinal);names=G.population_names(spinal);ix={n:i for i,n in enumerate(names)};rows=[]
    for side in ('r','l'):
        h=ix[f'HEEL_{side}']
        if np.any(c.g_max_non[h]):raise ValueError('Heel encoder has incoming paths; direct isolated input assumption invalid')
        heelv=float(c.V_rest[h]+(source_current_nA+c.i_b[h])/c.g_m[h]);out=[]
        for t in ('r','l'):
            for kind,rg in (('InE','RG_F'),('InF','RG_E')):
                i=ix[f'{kind}_{t}'];q=ix[f'{rg}_{t}'];geff=float(c.g_max_non[i,h])*G.activation_fraction(c,i,h,heelv)
                # Declared engineering bias holds relay at operating_mV without heel.
                v=(float(c.g_m[i])*operating_mV+geff*float(c.del_e[i,h]))/(float(c.g_m[i])+geff)
                dI=float(c.g_max_non[q,i])*(G.activation_fraction(c,q,i,v)-G.activation_fraction(c,q,i,operating_mV))*(float(c.del_e[q,i])-operating_mV)
                out.append({'via':f'{kind}_{t}','target':f'{rg}_{t}','baseline_relay_mV':operating_mV,'heel_high_relay_mV':v,'delta_target_current_nA':dI})
        rows.append({'heel_source':f'HEEL_{side}','input_current_nA':source_current_nA,'heel_encoder_steady_mV':heelv,'two_hop':out})
    return rows


def audit(family):
    import xml.etree.ElementTree as ET
    S=HERE.parent;repo=S.parents[2];sys.path.insert(0,str(S));os.environ['AARL_NET']=family
    import build_network as B,params as P
    graph=json.loads((repo/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families'/(family+'_graph_ledger.json')).read_text());acts=[a.attrib['name'] for a in ET.parse(graph['actuation']['model']).getroot().find('actuator')]
    n=B.build(acts);G=utilities();before=G.field_digest(G.raw_compiled(n));apply_policy(n,None)
    if before!=G.field_digest(G.raw_compiled(n)):raise AssertionError('None/default changed parameters')
    fixture='Raw canonical source constructor; no saved standing/body success inherited'
    if family=='s3k':
        try:apply_policy(n,make_policy(family,'A_ipsilateral_extension'))
        except ValueError as e:missing_rejected=str(e)
        else:raise AssertionError('Missing HEEL was not rejected')
        P.G.update(full_rules=1.,heel_rge=.5);n=B.build(acts);fixture='NEW contacts/full-rules-enabled source constructor: full_rules1,heel_rge.5; not standing baseline'
    else:missing_rejected=None
    baseline=isolated_static_operator(n);cases={}
    for hypothesis in ('A_ipsilateral_extension','B_contralateral_extension_inhibition'):
        unit=copy.deepcopy(n);record=apply_policy(unit,make_policy(family,hypothesis));rows=isolated_static_operator(unit)
        for row in rows:
            s=row['heel_source'][-1];t='l' if s=='r' else 'r';values={r['target']:r['delta_target_current_nA'] for r in row['two_hop']}
            expected={f'RG_F_{s}':-1,f'RG_E_{s}':1} if hypothesis=='A_ipsilateral_extension' else {f'RG_F_{s}':-1,f'RG_E_{t}':-1}
            if any(values[k]*sign<=0 for k,sign in expected.items()) or any(abs(v)>1e-14 for k,v in values.items() if k not in expected):raise AssertionError('Incorrect two-hop functional sign/laterality')
        final_verifier=verify_policy(unit);zero=isolated_static_operator(unit,0.)
        if any(abs(q['delta_target_current_nA'])>1e-14 for r in zero for q in r['two_hop']):raise AssertionError('Zero heel current has an effect')
        bad=make_policy(family,hypothesis);bad['unused_key']=True
        try:apply_policy(unit,bad)
        except ValueError:unused_rejected=True
        else:raise AssertionError('Unused policy key accepted')
        grid=[]
        for op in (0.,2.5,4.):
            sat={}
            for high in (0.,.125,.75,1.5,5.,20.):
                rr=isolated_static_operator(unit,high,op);grid.append({'operating_point_mV':op,'heel_current_nA':high,'rows':rr})
                for row in rr:
                    s=row['heel_source'][-1];t='l' if s=='r' else 'r';values={r['target']:r['delta_target_current_nA'] for r in row['two_hop']};expected={f'RG_F_{s}':-1,f'RG_E_{s}':1} if hypothesis=='A_ipsilateral_extension' else {f'RG_F_{s}':-1,f'RG_E_{t}':-1}
                    if any(values[k]*sign < -1e-14 for k,sign in expected.items()) or any(abs(v)>1e-14 for k,v in values.items() if k not in expected):raise AssertionError('Sign/laterality reversal in static domain grid')
                if high in (5.,20.):sat[high]=[[q['delta_target_current_nA'] for q in r['two_hop']] for r in rr]
            if sat[5.]!=sat[20.]:raise AssertionError('Expected heel synaptic saturation not consumed')
        cases[hypothesis]={'applied_fields':record,'final_scoped_verifier':final_verifier,'static_high':rows,'static_zero':zero,'static_domain_grid':grid,'domain_scope':'0..0.75nA covers W2L saved current; larger currents are engineering saturation probes,not admissible saved-caller inputs. Operating points are declared isolated fixtures.','disinhibition_floor':'A cannot disinhibit an already-silent InF output; at relay0mV its RG-E current change is zero,not positive.','unused_key_rejected':unused_rejected,'dense_operator_consumed':'Actual compiled g_max_non/del_e/e_lo/e_hi; exact clipped SNS conductance formula','sparse_scope':'Apply before accelerator construction; post-wrapper application rejected,not a sparse timestep parity claim'}
    return {'role':'Gait neural pathway composition reviewer','family':family,'fixture':fixture,'neurons':len(G.population_names(n)),'None_parameter_bitexact':True,'default_no_HEEL_enable_rejected':missing_rejected,'baseline_static':baseline,'cases':cases,'event_contract':{'pulse':'High operator applies while encoded heel onset is high; low operator is zero at equilibrium. No onset detector installed.','held':'Same nonzero high operator persists during held contact; it is tonic bias,not a one-off phase reset.','neural_time_constants':'Preserved; dynamic timing/release untested in this static audit.'},'source_current_probe_nA':.75,'probe_origin':'W2L saved caller heel maximum; separate engineering unit fixture for Syn6 whose standing caller gain is zero. Not a new caller contract.','operating_point_mV':2.5,'neuronal_time_steps':0,'physics_time_steps':0,'training':False,'production_wiring_edited':False,'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),S/'build_network.py',S/'build_network_w2lvar.py',S/'build_network_syn6.py')},'direct_contra_RG_F_variant':'Separate future hypothesis; not implemented or stacked'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--family',required=True,choices=('s3k','w2lvar','syn6'));p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=audit(a.family);a.out.mkdir(exist_ok=True,parents=True);f=a.out/(a.family+'_HEEL_AB_STATIC.json');f.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps({'family':a.family,'None_bitexact':r['None_parameter_bitexact'],'A_B_static_signs':True,'neural_steps':0,'physics_steps':0,'output':str(f)}))
