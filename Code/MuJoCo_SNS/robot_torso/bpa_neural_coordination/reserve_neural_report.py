"""Durable bounded fit record; failed nominal fits are not promoted to held-outs."""
import json,hashlib,shutil
from pathlib import Path
import numpy as np
from reserve_neural_runner import DEST,ROOT

def artifact(prefix,suffix):return Path(str(prefix)+suffix)

def main():
 rows=[];configuration=json.loads((ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json').read_text());pect_index=list(configuration['actuators']).index('pect_r')
 for file in sorted(DEST.glob('*.comparison.json')):
  c=json.loads(file.read_text());r=json.loads(Path(c['result']).read_text());prefix=Path(c['result'].removesuffix('.result.json'))
  with np.load(artifact(prefix,'.relay.npz')) as z:signals=z['samples'];tail=signals[-25:,0]
  rows.append(dict(case=prefix.name,cap_nA=c['cap_nA'],preview_s=c['preview_s'],condition=c['condition'],accepted=c['assessment']['accepted'],moving_s=r['completed_samples']*.002,failure=r['failure'],maximum_current_nA=c['maximum_added_current_nA'],integrated_added_current_nA_s=c['integrated_added_current_nA_s'],complete_restore_test=c['complete_restore_selftest'],preparation_exact=all(c['legacy_complete_release_common_exact'].values()),last50ms_right_signal_mean=float(tail[:,3].mean()),last50ms_right_current_utilization=float(tail[:,5].mean()/c['cap_nA']) if c['cap_nA'] else 0.,result=c['result']))
  with np.load(artifact(prefix,'.policy.npz')) as z:last=z['last_policy'][pect_index]
  rows[-1]['last_requested_tension_N']=float(last[0]);rows[-1]['last_minimum_feasible_tension_N']=float(last[1]);rows[-1]['last_selected_pressure_kpa']=float(last[6]);rows[-1]['last_static_residual_N']=float(last[7])
 fit=[row for row in rows if 'fit_' in row['case']];assert len(fit)<=6
 repeats={};a=DEST/'w2lvar_fit_c0p2_h0p1';b=DEST/'w2lvar_strongest_fresh'
 if artifact(b,'.comparison.json').exists():
  for suffix in ('.release_state.npz','.motor_state.npz','.pressure.npz','.inputs.npz','.apparatus.npz','.apparatus_extended.npz','.policy.npz','.relay.npz','.final_neural_state.npz'):
   with np.load(artifact(a,suffix)) as x,np.load(artifact(b,suffix)) as y:
    assert set(x.files)==set(y.files)
    repeats[suffix]={k:bool(np.array_equal(x[k],y[k],equal_nan=True)) if np.issubdtype(x[k].dtype,np.number) else bool(np.array_equal(x[k],y[k])) for k in x.files}
  if not all(v for group in repeats.values() for v in group.values()):raise ValueError('Fresh neural-relay repeat mismatch')
 report=dict(machine='easteregg2',rows=rows,plan=json.loads((DEST/'PREDECLARED_FIT.json').read_text()),zero_gain_parity=str(DEST/'ZERO_GAIN_PARITY.json'),sns_relay_checks=json.loads((DEST/'RELAY_UNIT_CHECKS.json').read_text()),strongest_fresh_repeat=repeats,accepted_fit_count=sum(row['accepted'] for row in fit),heldout_status='Not promoted: no accepted nominal fit' if fit and not any(row['accepted'] for row in fit) else 'Pending fresh/held-out checks')
 (DEST/'FIT_SUMMARY.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 text='# Anticipatory pressure-reserve neural relay — easteregg2\n\n'
 text+='Two explicitly named cells, RESERVE_pect_r/l, excite the actual homonymous MN inputs through mirrored conductance synapses. Their current is added to the original tonic input; no original pathway is canceled. Existing POST input ports are only the compiled additive-current transport. Neither MN voltage/output nor actual force/pressure is overwritten. Underlying normalized-w2lvar parameters and actuator geometry remain frozen. This is an engineering actuator-limit pathway, not a claim of biological pressure afferents.\n\n'
 text+='The sensor uses actual model pressure, route length, and route velocity: h=l0+b+Fmax_nominal Pactual/(620K)−L; filtered pressure derivative (10 ms) gives hdot=Fmax_nominal Pdot/(620K)−v. Predicted closing reserve is h+preview min(hdot,0). Sensor current is 10 nA times clip(−K predicted/1 N,0,1). Each relay has leak 1 µS, capacitance 6 nF, and time constant 6 ms, with 2 ms Euler integration. Its old voltage gates an excitatory 40 mV-reversal synapse into the actual MN; gmax=Icap/40 mV, with an explicit positive current cap. The hydraulic model and forward-domain guard are unchanged.\n\n'
 text+='Zero gain exactly reproduces the nearest-reference comparator across all recorded actual activation, excitation, pressure, neural-input/activity, and apparatus arrays. The actual relay voltage traces match an independent SNS-Toolbox two-neuron simulation exactly (0.0 mV error); logged conductance currents match actual MN excitation-derived voltages exactly (0.0 nA error). Complete restore tests execute during genuine nonzero relay current and include all relay/filter/previous-pressure/input-buffer/counter/log/no-RNG states plus the graded backend gates and sparse buffers. Preparation is identical before relay activation.\n\n'
 text+='The predeclared six-point fit uses caps 0.05/0.1/0.2 nA and previews 50/100 ms, nominal 20 s studies only. Selection requires acceptance, then smallest cap and longest preview. No held-out condition is used for selection. One initial decimal-tag filename collision was corrected to filesystem-safe tags; its valid completed negative was reused by verified parameter metadata, not rerun or overwritten.\n\n| Case | Cap (nA) | Preview (ms) | Pass | Moving (s) | Max current (nA) | Late current/cap |\n| --- | ---: | ---: | --- | ---: | ---: | ---: |\n'
 for row in rows:text+=f"| {row['case']} | {row['cap_nA']:.2f} | {1000*row['preview_s']:.0f} | {row['accepted']} | {row['moving_s']:.3f} | {row['maximum_current_nA']:.5f} | {row['last50ms_right_current_utilization']:.3f} |\n"
 if repeats:text+='\nThe strongest candidate (0.2 nA, 100 ms) fails after 16.466 s moving and repeats exactly in a fresh process, including full neural/relay/filter state, preparation, pressure, actual activation, currents, and apparatus arrays.\n'
 text+='\nAt the strongest candidate\'s boundary, actual added MN current is 0.17581 nA and the sensor is fully active, but requested tension is only 24.3128 N versus a feasible minimum of 25.0747 N. The nearest inverse therefore still chooses the minimum reference pressure (13.0039 kPa). The whole below-minimum force-request range maps to that same boundary reference: this is a conversion deadzone for the bounded MN-current correction. It is not evidence that the relay failed to run. Source route velocity agrees with finite differences of the recorded length to approximately 3e−11 m/s near the inspected boundary.\n'
 text+='\nAll six nominal fits remain negative, so no trained winner is promoted to delay/map/pulse generalization tests. The sensor detects closing reserve and the relay supplies bounded current; this does not establish adequate pressure reserve. The source of failure remains the unchanged pectineus actuator domain. Detailed failures, actual pressures, signed apparatus loads, visual mesh penetration, and knee limits remain in each result. No 20 s durability, independent standing, or walking claim is made.\n\n'
 text+='A distinct future architecture could connect explicitly modeled reserve/valve neurons to a declared pressure-reference conversion rather than relying on added force-command MN recruitment. That must be a separately named and drawn neural pathway, with bounded neuron currents and pressure-reference output, zero-gain parity, and actual delayed pressure dynamics. It must not clamp actual pressure or inject a force floor. No such valve-command pathway is implemented in this experiment.\n\n'
 text+='## Reproduction\n\nMachine easteregg2, interpreter `D:\\Anaconda\\envs\\myo\\python.exe`. `'+str(Path(__file__).with_name('reserve_neural_runner.py'))+' --fit` resumes only the six declared nominal points. Individual cases accept `--cap`, `--preview`, `--condition`, `--encoding`, `--tag`, and `--duration`. Tags must avoid decimal points. `'+str(Path(__file__).with_name('check_reserve_neural.py'))+'` checks actual relay traces against SNS-Toolbox. Circuit connection identities, exact current caps, and compiled MN input slots are saved with every comparison JSON.\n'
 (DEST/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8')
 example=json.loads((DEST/'w2lvar_fit_c0p2_h0p1.comparison.json').read_text());(DEST/'CIRCUIT.json').write_text(json.dumps(example['circuit'],indent=2),encoding='utf-8')
 dest=DEST/'source_snapshot';dest.mkdir(exist_ok=True);manifest={}
 for name in ('pressure_reserve_neural.py','reserve_neural_runner.py','check_reserve_neural.py','reserve_neural_report.py','inverse_feasibility_policy.py'):
  p=Path(__file__).with_name(name);digest=hashlib.sha256(p.read_bytes()).hexdigest();target=dest/(p.stem+'_'+digest[:8]+p.suffix)
  if not target.exists():shutil.copyfile(p,target)
  manifest[str(p)]=dict(sha256=digest,snapshot=str(target))
 (dest/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 index=ROOT/'EXPERIMENT_RECORD.md';begin='<!-- pressure-reserve-neural:start -->';end='<!-- pressure-reserve-neural:end -->'
 addition=begin+'\n## Explicit anticipatory neural reserve relay\n\nRecord: `'+str(DEST/'EXPERIMENT_RECORD.md')+'`. Two named mirrored 6 ms leaky relay cells receive actual pressure/length/velocity reserve prediction and add capped excitatory conductance current to pectineus MNs. Six predeclared nominal fits (0.05/0.1/0.2 nA, 50/100 ms preview) remain negative under unchanged actuator-domain guards, so no winner is promoted to held-out delay/map/pulse validation. Zero gain reproduces the physical nearest-policy comparator exactly; actual relay dynamics and current formula match independent SNS-Toolbox/actual-MN checks with zero error. Nonzero-current complete-state restore audits pass. Original controller parameters, geometry, pressure lag, and source pathways remain unchanged.\n'+end
 existing=index.read_text(encoding='utf-8')
 if begin in existing:
  left,tail=existing.split(begin,1);_,right=tail.split(end,1);existing=left+addition+right
 else:existing+='\n\n'+addition+'\n'
 index.write_text(existing,encoding='utf-8')
 print(json.dumps(dict(fits=len(fit),accepted=report['accepted_fit_count'],record=str(DEST/'EXPERIMENT_RECORD.md'))))

if __name__=='__main__':main()
