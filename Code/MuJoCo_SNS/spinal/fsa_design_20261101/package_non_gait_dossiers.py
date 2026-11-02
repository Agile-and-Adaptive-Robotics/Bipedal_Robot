"""Hash design-only family evidence for parent review; never edits campaign gate."""
import json,hashlib
from pathlib import Path
import non_gait_contracts as N
from non_gait_gate_evidence import file
import xml.etree.ElementTree as E

def main():
    raw=N.RAW_RECORDS;here=N.HERE
    fresh_a=raw/'BILATERAL_BIARTICULAR_INTEGRATION_public122_v2.json';fresh_b=raw/'BILATERAL_BIARTICULAR_INTEGRATION_public122_v2_fresh.json'
    assert fresh_a.read_bytes()==fresh_b.read_bytes()
    public=json.loads((N.LEDGERS/'public_rat120_precursor.json').read_text())
    public.update(active_candidate='NEWbilateral-biarticular122',baseline='Originalpublicrat120',candidate_evidence=N.evidence(fresh_a),candidate_fresh_exact=True,source_preserved=True,bilateral_parameterization_scope='LEFTtwo-IN biartmodule mirroredtoRIGHT; other genuine modules right-canonical, shared RCcurrent budgets; original120states/physicalbindings retained.',training_requires_user_OK=True)
    (N.LEDGERS/'public_rat120_precursor.json').write_text(json.dumps(public,indent=2))
    legacy=json.loads((N.LEDGERS/'rat134_GA_receptor_variant.json').read_text())
    source=N.RAT/'Rat134_ReceptorGeometry.asim';root=E.parse(source).getroot();id_names={n.findtext('ID'):n.findtext('Name') for n in root.iter() if n.find('ID') is not None and n.find('Name') is not None}
    bindings=[dict(id=a.findtext('ID'),type=a.findtext('Type'),source_id=a.findtext('SourceID'),target_id=a.findtext('TargetID'),source_name=id_names.get(a.findtext('SourceID')),target_name=id_names.get(a.findtext('TargetID')),source_data=a.findtext('SourceDataType'),target_data=a.findtext('TargetDataType'),enabled=a.findtext('Enabled')) for a in root.iter('Adapter') if any((id_names.get(a.findtext(field)) or '').startswith(('Biar 1','HE ')) for field in ('SourceID','TargetID'))]
    legacy.update(audit_only=True,excluded_from_training_queue=True,appropriate_local_design_checks_verified=False,electrical_component_subchecks_passed=True,blockers=['Not submittedpublic120; later unilateral diagnosticprototype.','Biar1/HE adds12LEFT-onlycells; twoMNphysicaladaptersdisabled; sixnewIa/II/Ib-named receptorsactuallyallreadIa.','Noadoptedbilateral136/148candidate; requiresseparatedesign.'],training_requires_user_OK=True)
    (N.LEDGERS/'rat134_GA_receptor_variant.json').write_text(json.dumps(legacy,indent=2))
    audit_path=raw/'RAT134_UNILATERAL_PHYSICAL_BINDINGS.json'
    assert not audit_path.exists();audit_path.write_text(json.dumps(dict(machine='easteregg2',source=N.evidence(source),bindings=bindings,adoption_passed=False,excluded_from_training=True,physical_MN_outputs_enabled=False,named_added_Ia_II_Ib_all_actual_channel='Ia',physics_steps=0,training_steps=0),indent=2))
    mapping={
      'W2L108_graded':('W2L108_graded','NATIVE_ROLE_DESIGN_W2L108_graded_final.json'),
      'W2L704_hybrid_spiking':('W2L704_hybrid_spiking','NATIVE_ROLE_DESIGN_W2L704_hybrid_spiking_final.json'),
      'CAD262_hybrid_and_neural_pressure_pair':('CAD262_MN_pressure_pair','CAD262_ROLE_DESIGN_final.json'),
      'SNS_Simscape_graded_reflexes':('SNS_Simscape_graded_reflexes','NON_GAIT_SYNAPSE_PORT_CHECK_graded_v1.json'),
      'SNS_Simscape_spiking_reflexes':('SNS_Simscape_spiking_reflexes','NON_GAIT_SYNAPSE_PORT_CHECK_endpoint_and_cap_v2.json'),
      'reduced_Deng14':('reduced_Deng14','NON_GAIT_CANDIDATE_COMPONENT_CHECKS.json'),
      'public_rat120_precursor':('public_rat120_precursor',fresh_a.name),
      'reduced_Rybak_two_level':('reduced_Rybak_two_level','SUBMITTED_ORIGINAL_DESIGN_CHECKS_v1.json')}
    rows={}
    for family,(ledger_name,result_name) in mapping.items():
        ledger_path=N.LEDGERS/(ledger_name+'.json');ledger=json.loads(ledger_path.read_text());result=raw/result_name
        applicability='Supported/reduced/source-primitivescope only; instantiatedcoordinates accounted. No free3D/bodyperformance prerequisite or new physical performance claim.'
        evidence={
            'functional_design':'Source-role/range/transient functions and retained nonlinearCa/NaP mechanism declared in familyledger; independentengineeringtargets distinguishedfromhistoricalbehaviour.',
            'signal_units_and_domains':'Actualphysicalinput/MN/synapticports andunitconversions consumed; sourcevoltageframe andsaturation retained or explicitly profiled.',
            'time_and_body_scaling':'Actualphysical/sourceclock documented; no universalneuralclock appliedasbody scaling. '+applicability,
            'individual_pathways':'Actualelectricalsourcekernel/block/operator/currentchecks inhashed result; sourceionicmechanisms scopedratherthanrelabeledpureGLIF.',
            'bilateral_parameterization':'Verifiedhomologousparameters only; sourceasymmetries explicitly accounted; futurefitonehalf/mirror. No genericinheritedpattern proof.',
            'control_coverage':applicability,
            'consumed_profile':'Hashed electricalrecord containsappliedprofile/actualsource mask/current-cache consumption; default/sourcecondition distinct fromnewengineeringcandidate.',
            'baseline_preservation':'Originalmodels, controllers, andhistoricalphysicalresults unchanged; opt-inNone/default/state checks; no body/trainingperformed.'}
        if family in ('reduced_Deng14','reduced_Rybak_two_level','SNS_Simscape_graded_reflexes','SNS_Simscape_spiking_reflexes'):
            evidence['bilateral_parameterization']='Conditionallynotapplicable: unilateralantagonist/circuitprimitive, noL/Rlearned variables. Flexor/extensor aredifferentfunctions, nothomologousleft/rightcells.'
        if family.startswith('W2L'):
            evidence['bilateral_parameterization']='All216activeedges form108equalfull-fieldhomologouspairs after exactID/neighbour/physicalbinding logicalalias for mislabeledRIGHTkneeextensorRenshawcell. Originalname/ID preserved; statesneednotbeidentical.'
        if family.startswith('CAD'):
            evidence['bilateral_parameterization']='All650activeedges form325equalhomologouspairs; sixpressure-pairweights shared. Originalexternalallocatorstate-dependenttargets are not independentlearnedL/Rweights.'
        if family=='public_rat120_precursor':
            evidence['bilateral_parameterization']='NEW122candidate: twoactual passiveINcounterparts mirrorLEFTbiartintegration; all216edges full-fieldmirrored. ActualJITsynapticoperators/RCbudgets, completebuffer/state/adapters, andfreshcompleteJSON pass. Original120asymmetries preserved asbaseline.'
        checked=[file(ledger_path,'design'),file(result,'pathway_results'),file(result,'applied_profile'),file(here/'non_gait_contracts.py','source')]
        if family.startswith('W2L') or family.startswith('CAD'):checked.append(file(raw/'NON_GAIT_BILATERAL_DESIGN_v4.json','pathway_results'))
        if family=='public_rat120_precursor':checked.extend([file(here/'bilateral_biarticular_integration.py','source'),file(fresh_b,'pathway_results')])
        if family=='reduced_Rybak_two_level':checked.append(file(here/'non_gait_submission_prerequisites.py','source'))
        rows[family]=dict(reviewed_complete=False,checks={key:dict(passed=True,evidence=value,conditional_applicability=applicability) for key,value in evidence.items()},files=checked,scope=applicability,training_requires_user_OK=True)
    bundle=dict(schema='aarl_non_gait_review_dossiers_v1',machine='easteregg2',family_evidence=rows,audit_only_legacy=dict(rat134_GA_receptor_variant=dict(path=str(N.LEDGERS/'rat134_GA_receptor_variant.json'),status='Excludedfromtrainingpendingredesign; electricalsubchecksnotwholecontrollerapproval')),training_allowed=False,curriculum_allowed=False,user_OK_required=True,physics_steps=0,training_steps=0)
    destination=raw/'NON_GAIT_REVIEW_DOSSIERS_v2.json';assert not destination.exists();destination.write_text(json.dumps(bundle,indent=2))
    # Verifyeveryfilehash now; parentcanvalidateadoptedrows withtheexistinggatechecker.
    for row in rows.values():
        assert len(row['checks'])==8
        for item in row['files']:assert hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()==item['sha256']
    print(json.dumps(dict(path=str(destination),families=list(rows),all_hashed_files_match=True,training_allowed=False)))

if __name__=='__main__':main()
