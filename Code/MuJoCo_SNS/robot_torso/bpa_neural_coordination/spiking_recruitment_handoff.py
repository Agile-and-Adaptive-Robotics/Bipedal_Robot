"""Read-only result figure and explicit CAD-transfer feasibility ledger."""
import json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import spiking_recruitment_experiment as S
D=S.DEST;HERE=Path(__file__).resolve().parent;ROBOT=HERE.parent;CAD=D.parent.parent/'cad_bpa'
def read(p):return json.loads(p.read_text())
def main():
 audit=read(D/'FRESH_REPEAT_AUDIT.json');assert audit['initial']['nominal_initial']['accepted'] and not audit['initial']['strong_initial']['accepted'] and audit['all_fresh_arrays_exact']
 stem=D/'strong_initial/w2lvar_case'
 with np.load(stem.with_suffix('.final_neural_state.npz')) as z:r=z['recruit_logs'];events=z['recruit_lif_events'];totals=z['recruit_lif_total']
 with np.load(stem.with_suffix('.pressure.npz')) as z:p=z['samples'];pt=z['time'];idx=z['names'].tolist().index('pect_r')
 nt=(np.arange(len(r))-300+1)*.002;et=(np.arange(len(events))+1)*.002
 plt.rcParams.update({'font.family':'Arial','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'});fig,axes=plt.subplots(4,1,figsize=(8.3,8),sharex=True,constrained_layout=True)
 axes[0].plot(nt,r[:,0,4],color='#AA3377');axes[0].set_ylabel('Error (0–1)');axes[0].set_title('(A) Right-pectineus measured-tension error input')
 axes[1].step(et,events[:,0],where='post',color='#4477AA');axes[1].set_ylabel('Spikes / 2 ms');axes[1].set_title('(B) Actual events across 128 LIF cells')
 axes[2].plot(nt,r[:,0,6],color='#EE7733');axes[2].set_ylabel('Added MN (nA)');axes[2].set_ylim(0,5.5);axes[2].set_title('(C) Decoded positive conductance current (cap 5 nA)')
 axes[3].plot(pt,p[:,1,idx],color='#4477AA',label='Actual pressure');axes[3].plot(pt,p[:,0,idx],color='#EE7733',ls='--',label='Reference');axes[3].set_ylabel('Pressure (kPa)');axes[3].set_xlabel('Released-physics time (s)');axes[3].set_title('(D) Unchanged pressure-response model');axes[3].legend()
 for ax in axes:ax.set_xlim(17.15,17.85);ax.grid(alpha=.2)
 fig.suptitle('Two-error-relay spiking mirror: right pectineus only\nKnown strong challenge FAILS later at 19.374 s',fontsize=12)
 for ext in ('png','pdf','svg'):fig.savefig(D/('SPIKING_RELAY_STRONG_WINDOW.'+ext),dpi=200)
 plt.close(fig)
 caption=f'Right pectineus only. Only two error relays are spiking (128 LIF each); 890 other modeled neurons remain graded. Event counts represent actual threshold/reset/refractory LIF spikes, summed in2ms bins; this is not a membrane-voltage raster. Full strong trajectory fails unchanged domain guard at19.374s released; displayed interval17.15–17.85s is not evidence of a completed challenge. The8s nominal case passes and repeats every array/state exactly. Total event counts right/left are {totals.tolist()}; left error population is silent in this condition, so no bilateral neural rhythm claim. Virtual BPA geometry and strong lateral/rotational guidance remain.'
 circuit=dict(graded_core_cells=888,graded_valve_cells=2,spiking_cells_per_population=128,spiking_population_names=['pect_r geometric_tension_error','pect_l geometric_tension_error'],population_parameters=read(D/'OFFLINE_CHECKS.json')['policy']|dict(individual_frequency_multiplier_range=[.8,1.2],decoder_gain=read(D/'OFFLINE_CHECKS.json')['decoder_gain']),actual_strong_event_counts_R_L=totals.tolist(),left_population_silent_in_observed_condition=True,parameters_mirrored=True,output_map='Measured delayed length/positive velocity and actual tension -> same dimensionless error -> LIF input drive -> real spike counts -> postsynaptic filter -> same positive conductance current into homonymous pectineus MN; cap5nA. No graded error-neuron shadow.',gmax_uS=.125,reversal_mV=40.,current_formula='clip(.125*clip(decoded,0,1)*(40-V_MN),0,5)nA',old_valve_path_unchanged=True,force_pressure_voltage_projection=False,oracle_gain_input=False,source_snapshot=str(D/'source_snapshot/spiking_recruitment_experiment.py'),full_state_and_nominal_fresh_audit=str(D/'FRESH_REPEAT_AUDIT.json'),strong_case_accepted=False)
 (D/'CIRCUIT.json').write_text(json.dumps(circuit,indent=2),encoding='utf-8')
 (D/'SPIKING_RELAY_STRONG_WINDOW_caption.txt').write_text(caption,encoding='utf-8');(D/'SPIKING_RELAY_STRONG_WINDOW_alt.txt').write_text('Four aligned short traces show right-sided error, population spike counts, bounded added motor current, and pressure. Title explicitly reports later strong-case failure.',encoding='utf-8')
 template=read(HERE/'CAD_BPA_PARAMETER_TEMPLATE.json');names=set(template['actuators']);fixture=read(Path(S.E.registration()['training'][0]['path']));human=set(fixture['actuators']);proposal=read(CAD/'mount_subset/PAIRED_REMOTE_SLEEVE_MOUNT_PROPOSAL.json');subset=read(CAD/'mount_subset/SLEEVE_PAIR_CLEARANCE.json')['adopted_simultaneous_sleeve_subset'];rows=[]
 for name in sorted(names):
  v=template['actuators'][name];m=v['reference_metadata_only'];rows.append(dict(name=name,exact_name_in_human_virtual_bpa=name in human,error_relay_tested_here=name in ('pect_r','pect_l'),simultaneous_candidate_sleeve=name in subset,missing_parameters=[k for k,val in v.items() if val is None],master_label=m.get('master_name'),pressure_model_runnable=not any(val is None for val in v.values())))
 ledger=dict(machine='easteregg2',cad_model=str(ROBOT/'ah001_neural.xml'),cad_actuator_count=len(rows),exact_name_overlap_with_human_virtual_bpa=len(names&human),tested_pectineus_present_on_CAD={n:n in names for n in ('pect_r','pect_l')},candidate_simultaneous_routes=subset,routes=rows,subset_geometry=[r for r in proposal['routes'] if r['name'] in subset],moving_interface=str(HERE/'cad_moving_force.py'),same_pressure_assignment_permitted=False,requirements=['Declare route-specific rest/nonmuscle lengths, series stiffness, parallel count, pressure limits, fill/vent/delay and uncertainty; template remains null, not automatically filled.','Disable all native gains/biases; apply complete moving-route -J.T@T exactly once. Preserve support force buffer.','NN force feedback from actual bank tension; geometry feedback from complete dynamic route length/velocity, not native zero actuator_force/fixed-site moment arms.','Do not map pectineus relays to vastus by name substitution. The action and signed torque must be re-audited and fitted explicitly.','Sleeve pair sample clearance is not mounted hardware identification or continuous swept-volume/fitting/cable validation.'],source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'CAD_BPA_PARAMETER_TEMPLATE.json',ROBOT/'robotbody_routes.json',CAD/'mount_subset/PAIRED_REMOTE_SLEEVE_MOUNT_PROPOSAL.json')})
 (D/'CAD_TRANSFER_FEASIBILITY.json').write_text(json.dumps(ledger,indent=2),encoding='utf-8')
 text=f'''# CAD transfer feasibility after spiking-relay check

Machine: easteregg2. The current CAD model has {len(names)} route actuators; {len(names&human)} have exact names in the tested human virtual-BPA assignment. **Neither pectineus route is present**, so the two tested error relays have no direct CAD recipient. Do not silently replace them with vastus or another muscle.

The source-derived simultaneous sleeve proposal remains vas_med_r/l only. Saved CAD geometry gives a63.43mm minimum mutual surface gap across17sampled poses; bifemsh sleeve envelopes overlap18.07mm and are not an accepted simultaneous pair. The active vas_med paths are labeled Vastus Intermedius in the master ledger; this discrepancy remains explicit.

The all-null54-actuator pressure template is intentionally non-runnable. Missing fields and exact name overlap are listed route-by-route in CAD_TRANSFER_FEASIBILITY.json. Existing geometric placements do not supply effective series stiffness, nonmuscle transmission length, valve response/delay, friction/hysteresis, component mass, or hardware force calibration. No values were invented or copied from the human's assumed sleeve design.

The existing moving-route bridge is {HERE / 'cad_moving_force.py'}. It applies complete moving-route virtual-work forces with all native gains/biases explicitly off, once per step. Neural force sensors must use actual bank tension because native actuator_force is zero; length/velocity sensors must use complete dynamic-route values. Existing finite-difference/power parity is recorded under {CAD / 'CAD_MOVING_FORCE_PARITY.json'}. No new CAD simulation was justified by the present unassigned parameter ledger.

Candidate coordinates: {CAD / 'mount_subset/PAIRED_REMOTE_SLEEVE_MOUNT_PROPOSAL.json'}.
Geometry figure: {CAD / 'mount_subset/CAD_SLEEVE_MOUNT_PROPOSAL.png'}.
'''
 (D/'CAD_TRANSFER_FEASIBILITY.md').write_text(text,encoding='utf-8');print(json.dumps(dict(spiking_counts=totals.tolist(),cad_routes=len(names),same_name_overlap=len(names&human),pectineus_in_cad=ledger['tested_pectineus_present_on_CAD'],figure=str(D/'SPIKING_RELAY_STRONG_WINDOW.png')),indent=2))
if __name__=='__main__':main()
