"""Build-only FSA/scaling applicability audit of actual legacy spiking selector."""
import os,sys,json,hashlib,argparse
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
HERE=Path(__file__).parent;sys.path[:0]=[str(HERE),str(HERE.parent)]
import numpy as np
from graded_profile import assert_campaign_allowed
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main(out):
    assert_campaign_allowed('unit_and_individual_pathway_probes');os.environ['AARL_NET']='spiking'
    import mujoco,build_network as B,build_network_spiking as S,sns_toolbox.backends as backend
    from params import DT,TAU
    model=HERE.parents[3]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
    files=[Path(B.__file__),Path(S.__file__),Path(backend.__file__),HERE.parent/'params.py',HERE.parent/'spiking_calibration.json'];before={str(p):sha(p) for p in files}
    m=mujoco.MjModel.from_xml_path(str(model));n=B.build([m.actuator(i).name for i in range(m.nu)]);c=n.compiled
    assert isinstance(n,S.SpinalNetworkSpiking)
    names=[p['name'] for p in n.net.populations];sp=np.array([bool(p['type'].params.get('spiking',False)) for p in n.net.populations]);pure=sp & (c.theta_increment==0)
    classes={}
    for kind in ('aff','in','rc','rg','pf'):
        neu=S._spk(float(TAU['afferent']),kind);p=neu.params;reset=float(p['reset_potential']);theta=float(p['threshold_initial_value'])-reset;bias=float(p['bias'])+float(p['membrane_conductance'])*(float(p['resting_potential'])-reset)
        classes[kind]=dict(raw_parameters=p,reset_origin_theta0_mV=theta,reset_origin_effective_bias_nA=bias,literal_Table2_m0_bias_nA=float(p['membrane_conductance'])*theta/2,Table2_bias_equal=bool(abs(bias-float(p['membrane_conductance'])*theta/2)<1e-12),scope='Class probe uses supplied afferent membrane tau only to expose class electrical definitions; actual cell taus are separately inventoried below')
    def bounds(a,mask):
        v=np.asarray(a)[mask];return [float(v.min()),float(v.max())] if len(v) else None
    checks={
      '1_Fmax_R_theta0':dict(status='partial_legacy',reason='Absolute thresholds/rest/reset are explicit, but f_ref40Hz and readoutF_FULL50Hz are nominal calibration points, not a declared per-cell maximum or physical input domain.'),
      '2_threshold_m':dict(status='partial_legacy',reason='All source proportional threshold slopes are0. RG/PF additionally use post-spike threshold increments; that is not the Table2 voltage-dependent m adaptation.'),
      '3_tau_theta':dict(status='not_derived',reason='RG reuses graded rg_nap_h; PF uses0.15s. No tau_theta=tau_graded*(1-m/2) derivation. Nonadapting relays have unused0.1s threshold pole.'),
      '4_Ibias':dict(status='not_derived',reason='All spiking cells retain16nA tonic bias. Even after the required reset-origin coordinate conversion, their effective biases differ from G*theta0/2. See electrical class ledger.'),
      '5_tau_mem':dict(status='not_derived',reason='Copies graded membrane tau by class, not R/(Fmax*theta_star). Numeric capacitance uses microfarads with seconds; numericC*1000 givesnF.'),
      '6_tau_syn_delta':dict(status='different_model',reason='5/20ms decay, additive increments capped atGmax (usually8increments), forward-Euler decay. This is not Table2 event-reset-to-Gmax and has no declared delta bound.'),
      '7_Gmax_function':dict(status='legacy_calibration_not_Table2',reason='Per-class CAL multipliers and nominal40/50Hz two-cell bisections. No full sender/receiver operating-domain conductance derivation or fan-in proof.'),
    }
    report=dict(family='production_spiking_mirror',selector='AARL_NET=spiking',actual_class=type(n).__name__,source_files_sha256=before,source_unchanged=before=={str(p):sha(p) for p in files},census=dict(total=len(names),actual_spiking=int(sp.sum()),remaining_graded=int((~sp).sum()),constant_threshold_LIF=int(pure.sum()),post_spike_adapting=int((sp & ~pure).sum()),graded_MN=sum(x.startswith('MN_') for x in names),graded_readout_sinks=sum(x.startswith('RD_') for x in names),chemical_spike_edges=int(np.count_nonzero(c.g_max_spike)),chemical_graded_edges=int(np.count_nonzero(c.g_max_non))),clock=dict(physical_host_s=float(DT),neural_substep_s=float(c.dt),substeps=n.n_sub,actual_elapsed_per_host_s=float(c.dt*n.n_sub),global_clock_knob_present=False,carrier_and_rhythm_share_substep=True,clock_scope='Substepping changes numerical resolution, not a5x biological cadence multiplier'),source_calibration=S.CAL,electrical_class_reset_coordinate_audit=classes,actual_spiking_membrane_tau_s=bounds(c.c_m/c.g_m,sp),actual_spiking_effective_C_nF=bounds(1000*c.c_m,sp),raw_theta0_mV=bounds(c.theta_0,sp),actual_spike_synapse_time_constants_s=sorted(set(float(v) for v in c.tau_syn[c.g_max_spike>0])),readout=dict(raw_cell_indices=n.ridx,runner_indices={k:n.idx[k] for k in n.ridx},sink_pole_s=.05,nominal_full_scale_rate_Hz=S.F_FULL,scope='Sink-only graded taps; raw RG/PF chemical outputs still use actual spikes. The taps are not evidence of Table2 transmitter calibration.'),steps=checks,NaP_scope='This mirror replaces graded NaP RG with adapting LIF; it does not preserve a spiking NaP conductance or the acceptedv4 observer.',opt_in_FSA_compatibility='Pure constant-threshold relay redesign is mathematically possible after declared input/rate ranges. Adaptive RG/PF and spiking receiver conductance budgets require separate design; the graded-receiver soleus input-port proofs do not cover this mirror.',FSA_design_verified=False,training_allowed=False,no_body_or_neural_time_integration=True)
    report['use_inventory']=dict(classification='Legacy reachable API/template; not the current accepted central-spiking standing source',current_selected_anchor=str(Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_dynamic_readout_v4\TIMING_CANDIDATE.json')),selected_anchor_implementation='TimingZeroSFA custom256-cell central wrapper overgraded syn6, not this310-spiking-cell mirror',reachable_via=str(HERE.parent/'build_network.py')+' AARL_NET=spiking selector used by genericrunner',direct_nonhistorical_references=[str(HERE.parent/x) for x in ('calibrate_spiking.py','draw_spiking_schematic.py','make_editor_templates.py')],historical_experiments=str(HERE.parent/'reports_spiking_20261002'),evidence_scope='Read-only source-reference and accepted-anchor audit; no proof that an external user never selects this publicAPI',training_dependency_decision='Root must classify intendedactive targets. This unused-by-current-anchor API should not silently become a prerequisite for anothercontroller, but requires its own design before anyone trainsit.')
    assert report['source_unchanged'];out=Path(out);out.mkdir(parents=True,exist_ok=True);(out/'PRODUCTION_MIRROR_APPLICABILITY.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(census=report['census'],clock=report['clock'],source_unchanged=True,FSA_design_verified=False,use_classification=report['use_inventory']['classification'])),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();main(a.output)
