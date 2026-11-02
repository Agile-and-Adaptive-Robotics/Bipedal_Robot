"""Read-only analytical range/compatibility manifests from actual Gait ledgers.

No network construction, optimizer, neural step, physical integration, or NMF.
Missing physical domains stay unresolved. A saved successful trace is not used
as a design-domain endpoint. Numerical probes require separate root approval.
"""
import argparse,json,hashlib,math
from pathlib import Path
from table2_glif import analytic_passive_voltage_bound
HERE=Path(__file__).parent
REPO=HERE.parents[3]
LEDGERS=REPO/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families'
SYNERGY=REPO/'Code/MuJoCo_SNS/spinal/synergy_basis.npz'

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def profile(family,declared_input_domains=None,activation_choices=None):
    path=LEDGERS/(family+'_graph_ledger.json');source_hash=digest(path);z=json.loads(path.read_text());domains=declared_input_domains or {};steps=[];nodes=[]
    outgoing={n['index']:[] for n in z['nodes']};incoming={n['index']:[] for n in z['nodes']};inputs={n['name']:[] for n in z['nodes']}
    for e in z['chemical_edges']:outgoing[e['source_index']].append(e);incoming[e['target_index']].append(e)
    for a in z['input_adapters']:inputs[a['target']].append(a)
    for n in z['nodes']:
        p=n['source_neuron_parameters'];ionic=any(float(v)!=0 for v in p.get('Gion',[]));unresolved=[];lo=hi=0.
        for a in inputs[n['name']]:
            d=domains.get(a['input'])
            if d is None:unresolved.append(a['input']);continue
            if len(d)!=2 or d[0]>d[1]:raise ValueError('Invalid current bounds for '+a['input'])
            vals=[a['weight']*v for v in d];lo+=min(vals);hi+=max(vals)
        gm=p['membrane_conductance'];cm=p['membrane_capacitance'];rest=p['resting_potential'];bounds=analytic_passive_voltage_bound(gm,rest,p['bias'],None if unresolved else [lo,hi],[(e['gmax_uS'],e['reversal_mV']) for e in incoming[n['index']]])
        if ionic:bounds=dict(resolved=False,reason='NaP/Ca ionic state invalidates passive convex equilibrium envelope')
        pole=1-z['dt_s']*gm/cm
        nodes.append(dict(name=n['name'],source_index=n['index'],compatibility='ionic_burster_NOT_pure_GLIF_Table2' if ionic else 'passive_graded_candidate_for_GLIF_function_design',source_parameters=p,source_continuous_leak_tau_s=cm/gm,source_discrete_leak_pole=pole,source_discrete_equivalent_tau_s=-z['dt_s']/math.log(pole) if 0<pole<1 else None,source_numeric_C=cm,effective_C_nF_if_G_is_uS=1000*cm,input_current_nA_bounds=None if unresolved else [lo,hi],unresolved_inputs=unresolved,passive_equilibrium_bound=bounds,outgoing_activation_intervals_mV=sorted({(e['threshold_mV'],e['saturation_mV']) for e in outgoing[n['index']]}),synaptic_fanin=len(incoming[n['index']]),inhibitory_fanin=sum(e['reversal_mV']<rest for e in incoming[n['index']]),fsa_function_target=None,fsa_function_target_status='Must specify quantity/sign/range/operating point; existing Gmax is not itself a functional design target'))
    reasons=['Network-wide Fmax,R,theta0 and physical range/fastest-cycle decisions are pending; no arbitrary paper example substituted.','Per-pathway transient function m has not been selected.','Graded effective tau is recorded, but conductance operating-point calibration and Eq26 verification remain pending.','Eq21 is implemented in analytical module; bias is not yet applied to any model.','Eq23 is implemented; exact/rate endpoint and input-domain checks remain pending.','Eq27/28/29 are implemented; delta,bounds,clock convergence,and reset-synapse tests remain pending.','Function-derived Gmax and simultaneous fan-in/offset design remain pending; inherited matrix values are read-only.']
    for k,reason in enumerate(reasons,1):steps.append(dict(step=k,status='pending_design_or_verification',reason=reason))
    return dict(machine='easteregg2',family=family,source_ledger=str(path),source_ledger_sha256=source_hash,source_candidate=z['saved_configuration'],source_candidate_sha256=z['saved_configuration_sha256'],source_settings=z['settings'],source_census=z['census'],source_matrix_sha256=z['compiled_matrix_sha256'],source_NMF_basis=dict(path=str(SYNERGY),sha256=digest(SYNERGY),policy='READ_ONLY; ZCode owns8-synergy redesign; this module never recomputes or substitutes the K6 basis'),unit_adapter=dict(source_dt_s=z['dt_s'],source_update='dU=dt_s/C_numeric*I_numeric; effective tau_s=C_numeric/G_numeric',physical_equivalent='Cm_nF=1000*C_numeric if G is uS and I is nA',warning='Do not label source numerical C=0.05 as0.05nF when its implemented leak tau is50ms'),activation_choices=activation_choices,table2_steps=steps,nodes=nodes,architecture_scope='Proposed GLIF relay/MN hybrid only. NaP RG remains graded first; full RG spiking requires separate ionic-spiking design and validation.',integration_status='No optional adapter installed in current SNS backend. Actual conductance port exists in standalone GLIF runtime; source backend interception and pathway parity remain unimplemented.',individual_pathways_verified=False,training_allowed=False,curriculum_allowed=False,physical_trials_allowed=False)

def write_profiles(output,domains=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    before=digest(SYNERGY)
    for family in ('s3k','w2lvar','syn6'):
        data=profile(family,domains);(output/(family+'_DESIGN_PROFILE.json')).write_text(json.dumps(data,indent=2),encoding='utf-8')
        assert digest(data['source_ledger'])==data['source_ledger_sha256']
    assert digest(SYNERGY)==before

if __name__=='__main__':
    p=argparse.ArgumentParser(description='Static source-ledger design inventory only; no probes or training');p.add_argument('--output',required=True);p.add_argument('--declared-input-domains');a=p.parse_args()
    write_profiles(a.output,json.loads(Path(a.declared_input_domains).read_text()) if a.declared_input_domains else None)
