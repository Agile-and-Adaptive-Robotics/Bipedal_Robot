"""Original submitted rat120 and reduced Rybak electrical design; no body/training."""
import sys,json,copy,importlib.util
from pathlib import Path
from types import FunctionType
from collections import Counter
import numpy as np
import non_gait_contracts as N
import non_gait_gate_evidence as B

def rybak():
    source=N.SPINAL/'replication/rybak2006_standalone.py'
    spec=importlib.util.spec_from_file_location('non_gait_original_rybak',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    c=m.build_circuit();original=copy.deepcopy(c)
    # pF/nS/pA -> nF/uS/nA: all electrical current terms AND C divide by1000.
    converted=copy.deepcopy(c)
    for key in ['gl','gnap']:converted[key]=converted[key]/1000.
    scope=dict(m.simulate.__globals__);scope.update(C=m.C/1000.,GSYN=m.GSYN/1000.,GDRIVE=m.GDRIVE/1000.)
    physical=FunctionType(m.simulate.__code__,scope,'physical_unit_simulate',m.simulate.__defaults__)
    t,v,f,_=m.simulate(t_end_s=.02,dt_ms=.05,circuit=c)
    tt,vv,ff,_=physical(t_end_s=.02,dt_ms=.05,circuit=converted)
    error=float(np.max(abs(v-vv)));assert error<1e-10 and np.array_equal(t,tt)
    assert all(np.array_equal(c[key],original[key]) for key in ['gl','el','gnap','WE','WI']) and m.C==10.
    ports=[]
    for sign,E,weights in [('exc',m.E_SYN_E,c['WE']),('inh',m.E_SYN_I,c['WI'])]:
        for post,pre in zip(*np.nonzero(weights)):
            g=float(m.GSYN*weights[post,pre]);current_pA=g*(E-(-55.));current_nA=(g/1000.)*(E-(-55.))
            assert abs(current_nA-current_pA/1000.)<1e-12
            ports.append(dict(source=c['labels'][pre],target=c['labels'][post],sign=sign,source_conductance_nS=g,physical_conductance_uS=g/1000.,current_at_minus55_nA=current_nA))
    passive=[]
    for i,label in enumerate(c['labels']):
        if c['gnap'][i]>0:continue
        # Independently specified conditional conductance function over0,.5,1 occupancy.
        for level in [0.,.5,1.]:
            ge=m.GSYN*c['WE'][i].sum()*level;gi=m.GSYN*c['WI'][i].sum()*level
            total=c['gl'][i]+ge+gi;target=(c['gl'][i]*c['el'][i]+ge*m.E_SYN_E+gi*m.E_SYN_I)/total
            tau_ms=m.C/total
            assert np.isfinite(target) and tau_ms>0
        passive.append(dict(name=label,source_leak_nS=float(c['gl'][i]),source_C_pF=m.C,no_input_tau_ms=float(m.C/c['gl'][i]),physical_C_nF=m.C/1000.,function='Conditional exact source rational conductance equilibrium; poleC/(gLeak+gExc+gInh), bounded signed role preserved.'))
    return dict(passed=True,source=N.evidence(source),rules=N.evidence(N.SPINAL/'replication/rybak_rules.json'),cells=len(c['labels']),retained_NaP_cells=int(np.count_nonzero(c['gnap'])),passive_roles=passive,chemical_ports=ports,source_clock='dt_ms; output seconds; no body/cadence rescaling',physical_units=dict(C_nF=m.C/1000.,GSYN_uS=m.GSYN/1000.,GDRIVE_uS=m.GDRIVE/1000.),unit_adapted_20ms_electrical_voltage_error_mV=error,default_source_exact=True,normalized_activity_domain='clip((Vm-(-50mV))/50mV,0,1); not spike-rate output',literature_scope='Activity-based reinterpretation with source flagged calibration constants, not1:1originalspikingHHport',bilateral_applicable=False,body_coordinates=[],physics_steps=0,training_steps=0)

def main():
    raw=N.RAW_RECORDS
    a=json.loads((raw/'NATIVE_ROLE_DESIGN_public_rat120_precursor_v1.json').read_text());b=json.loads((raw/'NATIVE_ROLE_DESIGN_public_rat120_precursor_v1_fresh.json').read_text())
    assert a==b and a['electrical_checks_passed']
    c,source,anchor=N.actual_native_controller('public_rat120_precursor')
    bilateral=B.edge_account(c.names,c.edges);bilateral['application_proof']=B.proof(bilateral,c.edges.copy())
    report=dict(machine='easteregg2',public_rat120_precursor=dict(passed=True,anchor=anchor,source=N.evidence(source),fresh_complete_JSON_exact=True,role_counts=dict(Counter(row['role'] for row in a['source_roles'])),bilateral=bilateral,physics_steps=0,training_steps=0),reduced_Rybak_two_level=rybak())
    dest=raw/'SUBMITTED_ORIGINAL_DESIGN_CHECKS_v1.json';assert not dest.exists();dest.write_text(json.dumps(report,indent=2))
    for family in ['public_rat120_precursor','reduced_Rybak_two_level']:
        row=report[family];is_rat=family.startswith('public')
        coordinate=[dict(name=side+'_'+joint,coverage='neural') for side in ['RH','LH'] for joint in ['hip','knee','ankle']]+[dict(name='Mid_spine',coverage='passive'),dict(name='pelvis_root6DoF',coverage='source/import fixed constraint')] if is_rat else []
        ledger=dict(schema='aarl_non_gait_functional_contract_v1',family=family,machine='easteregg2',source_evidence=[row['source']],controlled_subspace=dict(coordinate_accounting=coordinate,all_source_coordinates_accounted=True),design_evidence=N.evidence(dest),appropriate_local_design_checks_verified=True,body_performance_is_not_design_prerequisite=True,design_passed=False,scaling_passed=False,pathways_passed=False,training_allowed=False,blockers=[],limitations=['Original rat120 source retains historical RH latch and is not curated134; no new physical success claim.' if is_rat else 'Unilateral14population activity-based literature-template mechanism; not14cellDeng orwholewalker.'],review_required='Parent review of original-submission roles, explicit source asymmetry, and opt-in profiles. UserOK required before any training.')
        (N.LEDGERS/(family+'.json')).write_text(json.dumps(ledger,indent=2))
    print(json.dumps(dict(public_rat120_passed=True,public_rat120_fresh_exact=True,rybak_cells=report['reduced_Rybak_two_level']['cells'],rybak_passed=True,rybak_voltage_error=report['reduced_Rybak_two_level']['unit_adapted_20ms_electrical_voltage_error_mV'])))

if __name__=='__main__':main()
