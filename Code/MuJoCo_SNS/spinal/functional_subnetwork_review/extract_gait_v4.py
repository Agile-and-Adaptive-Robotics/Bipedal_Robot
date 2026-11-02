"""Expand saved population identities; preserve chemical edges at their real interface."""
import json
from audit_saved_graphs import ROOT, OUT, table, digest
ledger=ROOT/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families/syn6_graph_ledger.json'
candidate=ROOT/'tmp/functional_subnetwork_remote/named_conditions/gait_v4/TIMING_CANDIDATE.json'
d=json.loads(ledger.read_text());cfg=json.loads(candidate.read_text())
assert digest(ledger)==cfg['source_graph_sha256']
central=[n for n in d['nodes'] if n['name'].startswith(('RG_','PF_S'))]
assert len(central)==16 and sum(n['name'].startswith('RG_') for n in central)==4
rows=[]
for n in d['nodes']:
    is_central=n in central
    if not is_central:
        rows.append({**n,'effective_id':n['id'],'effective_kind':'graded retained','soma_index':None})
        continue
    rg=n['name'].startswith('RG_')
    for cell in range(16):
        rows.append({**n,'effective_id':n['id']+f':soma:{cell}','effective_kind':'NaP spiking soma' if rg else 'fixed-threshold spiking soma','soma_index':cell,'threshold_relative_mV':(4. if rg else 1.)*(1+(.02 if rg else .15)*(cell/15.-.5)),'reset_relative_mV':2.5 if rg else 0.,'refractory_s':.001,'spike_history_tau_s':.025,'RG_adaptation_increment_nA':0.,'RG_h_tau_multiplier':cfg['parameters']['rg_h_tau_scale'] if rg else None,'parameter_scope':'source_neuron_parameters are construction fields; runtime compiled state and population transform remain separately versioned'})
assert len(rows)==1034
out=OUT/'Gait_TimingZeroSFA'
table(out/'effective_cell_identities.csv',rows)
table(out/'population_chemical_edges.csv',[{**e,'interface':'Graded conductance gated by spike-derived proxy voltage for central endpoints; NOT a direct pairwise event synapse'} for e in d['chemical_edges']])
steps=[
{'step':1,'operation':'Fmax,R,theta0','finding':'No declared network-wide Table2 Fmax/R/theta0 contract. RG/PF soma thresholds and source voltage window differ. Encoding rate is not model bandwidth.','status':'incomplete'},
{'step':2,'operation':'role-dependent m','finding':'Soma thresholds fixed per cell, m=0 special case. RG NaP h is retained separately; zero adaptation increment does not imply removal of the h gate.','status':'fixed-threshold special case; role design not demonstrated'},
{'step':3,'operation':'tau_theta','finding':'No adaptive-threshold state. 25ms AP-history and scaled NaP h time constant are not tau_theta.','status':'not an adaptive-threshold implementation'},
{'step':4,'operation':'derived bias Eq21','finding':'Bias copied from source compiled neuron. No evidence of recalculation from Gmem*theta0/(2-m).','status':'not demonstrated'},
{'step':5,'operation':'tau_mem and rate range Eq23','finding':'Source C/G plus conductance and NaP dynamics retained. Fixed reset/refractory, nonlinear inverse passage-time decoder. No declared Fmax rate-range proof by Eq23.','status':'not demonstrated'},
{'step':6,'operation':'tau_s/delta Eq27','finding':'25ms exponential spike-history estimator with event-marked conductance and h; no declared delta or equivalent reset-to-peak chemical event law. Do not relabel observer memory as compliant tau_s.','status':'different interface; literal step not demonstrated'},
{'step':7,'operation':'functional Gmax Eq29','finding':'Saved syn6 graded chemical conductances consume spike-derived proxy voltage. They are not newly derived Eq29 spike-event conductances.','status':'graded interface retained; literal step not demonstrated'}]
table(out/'TABLE2_AUDIT.csv',steps)
(out/'INTERFACE_CONTRACT.json').write_text(json.dumps({'configuration_sha256':digest(candidate),'base_graph_sha256':digest(ledger),'cells':{'spiking':256,'graded':778},'central_population_count':16,'source_graph_policy':'Preserve original syn6 graph at population interfaces; no invented all-to-all soma edges.','readout':'Causal rate/event-marked g/h inverse decoder; silent rate <0.1 maps to -4mV. Support weighting uses count=rate*16*(2*.025).','clocks':{'internal_spike_s':.0002,'neural_physics_s':.002},'source_state_warning':'Frozen physical continuation and initial compiled state were not reconstructed or stepped in this static extraction.','design_complete':False},indent=2))
print(json.dumps({'effective_cells':len(rows),'chemical_population_edges':len(d['chemical_edges']),'Table2_complete':False}))
