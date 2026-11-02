"""Report matched interface-policy experiments without hiding failed forward guards."""
import json,hashlib,shutil,ast
from pathlib import Path
import numpy as np
from normalized_prepared_transfer import ROOT
from inverse_feasibility_policy import feasible_interval
from cad_interface import RoutedBPABank

def main():
 out=ROOT/'inverse_boundary_policy';rows=[];audits={}
 config=json.loads((ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json').read_text());bank=RoutedBPABank(list(config['actuators']),config,.002)
 for path in sorted(out.glob('*.comparison.json')):
  c=json.loads(path.read_text());r=json.loads(Path(c['result']).read_text());prefix=Path(c['result']).with_suffix('').with_suffix('')
  with np.load(prefix.with_suffix('.policy.npz')) as z:events=z['events'];names=z['route_names'];channels=list(z['channels'])
  row=dict(case=prefix.name,accepted=c['assessment']['accepted'],policy=c['policy'],encoding=c['encoding'],moving_s=r['completed_samples']*.002,wall_s=r['wall_s'],failure=r['failure'],events=len(events),maximum_static_residual_N=float(np.max(abs(events[:,channels.index('static_tension_residual_N')]))) if len(events) else 0.,release_all_exact=all(c['legacy_release_state_exact'].values()),result=c['result'])
  if r['failure'] and 'outside the declared contraction domain' in r['failure']:
   d=ast.literal_eval(r['failure'].split('failing routes: ',1)[1])[0];i=bank.names.index(d['name']);pmin,pmax,tmin,tmax=feasible_interval(d['route_length_m'],bank.parameters[i]);row['dynamic_domain_failure']=dict(route=d['name'],actual_pressure_kpa=d['actual_pressure_kpa'],minimum_static_pressure_at_failure_kpa=pmin,pressure_deficit_kpa=pmin-d['actual_pressure_kpa'],excess_length_mm=d['excess_over_pressure_supported_length_mm'])
  rows.append(row)
 for encoding in ('graded','spike'):
  a=out/f'w2lvar_nearest_{encoding}_20s_v1';b=out/f'w2lvar_nearest_{encoding}_20s_fresh'
  if not b.with_suffix('.comparison.json').exists():continue
  suffixes=['.release_state.npz','.motor_state.npz','.pressure.npz','.inputs.npz','.apparatus.npz','.apparatus_extended.npz','.policy.npz']
  if encoding=='spike':suffixes+=['.spikes.npz','.final_spiking_state.npz']
  result={}
  for suffix in suffixes:
   with np.load(a.with_suffix(suffix)) as x,np.load(b.with_suffix(suffix)) as y:
    assert set(x.files)==set(y.files)
    result[suffix]={k:bool(np.array_equal(x[k],y[k],equal_nan=True)) if np.issubdtype(x[k].dtype,np.number) else bool(np.array_equal(x[k],y[k])) for k in x.files}
  if not all(v for group in result.values() for v in group.values()):raise ValueError('Fresh policy repeat mismatch: '+encoding)
  audits[encoding]=dict(all_equal=True,comparisons=result)
 report=dict(machine='easteregg2',rows=rows,fresh_repeat_audits=audits,unit_gate=str(out/'BOUNDARY_UNIT_CHECKS.json'))
 (out/'COMPARISON.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 text='# Inverse-feasibility boundary policy — easteregg2\n\n'
 text+='The frozen normalized-w2lvar controller and virtual actuator geometry are unchanged. Original zero-sentinel failures remain in `'+str(ROOT/'normalized_prepared_transfer')+'`. Policies activate only after exactly matched complete preparation. This is an interface-policy comparison, not new neural training or calibrated CAD hardware validation.\n\n'
 text+='For route length L, fixed length b, sleeve rest length l0, series stiffness K, and tensile force T, sleeve length is L−b−T/K. Nonnegative contraction requires Tmin,geometry=max(0,K(L−b−l0)). At zero contraction, F=Fmax P/620, giving Pmin=620 Tmin,geometry/Fmax. Feasible equilibrium tension endpoints are computed with the unchanged forward law at Pmin and Pmax, including the original zero-force shortened-route branch. Empty intervals abort. Only directed ULP rounding is used at the lower endpoint; no domain tolerance, arbitrary pressure margin, or actual-state projection is added.\n\n'
 text+='Strict abort reports an infeasible request before changing pressure/history/force. Nearest policy selects the closest attainable static reference, logs tension residuals and original inverse flags, and advances actual pressure with the original lag/delay/history. Bounds use nominal actuator parameters, not held-out force-map oracle information. Actual uncertain forward-force guards remain unchanged.\n\n'
 text+='Unit checks cover four actual failure geometries, 258 design-length regimes, atomic abort, unchanged rejection of historical invalid pressure states, and bit-exact old/new 86-route pressure/history/force evolution over 20 feasible steps.\n\n| Case | Pass | Moving (s) | Saturation events | Maximum static residual (N) |\n| --- | --- | ---: | ---: | ---: |\n'
 for row in rows:text+=f"| {row['case']} | {row['accepted']} | {row['moving_s']:.3f} | {row['events']} | {row['maximum_static_residual_N']:.6f} |\n"
 text+='\n## Results and interpretation\n\n'
 for row in rows:
  if row['failure']:text+='- '+row['case']+': '+row['failure']+'\n'
 text+='\nA feasible static pressure reference does not ensure that delayed actual pressure remains feasible while the route lengthens. The nearest-policy negative must not be described as a solved 20 s standing result. Actual pressure deficits are retained in the JSON, as are all signed apparatus/mesh/knee diagnostics in each original result. Fresh negative repeats test reproducibility without converting failure to acceptance.\n\n'
 text+='## Separate neural-feedback assessment\n\n'
 text+='The static residual is an explicitly engineering actuator-limit signal, not biological Ia, II, or Ib. A future experiment could route positive lower-limit residual through named limit afferents/interneurons to the homonymous MN and negative upper-limit residual through an inhibitory interneuron. Bilateral gains must be tied, additional cells/edges drawn, inputs/voltages logged, and zero gain must recover the frozen network. Direct motor-command edits or force allocation must not be called neural coordination.\n\n'
 text+='Residual-only feedback activates after a request is already infeasible and may be too late for pressure lag. An anticipatory candidate would sense contraction reserve h=l0+b+Fmax Pactual/(620K)−L and its decrease, using measured pressure plus route length/velocity. Any preview horizon, reserve threshold, or gain is a new declared engineering parameter. Pressure/length uncertainty and withheld perturbations must be tested. A pressure-reference barrier is a separate actuator-local control policy; it must not be hidden inside a claimed neural pathway. No such feedback or barrier is implemented here. Existing Ib already reports force but does not explicitly encode inverse-command infeasibility.\n\n'
 text+='## Reproduction\n\nMachine easteregg2; `D:\\Anaconda\\envs\\myo\\python.exe` with `'+str(Path(__file__).with_name('check_inverse_feasibility_policy.py'))+'` for unit gates, then `'+str(Path(__file__).with_name('boundary_policy_runner.py'))+' --batch --duration 20`. Completed cases are reused. Per-policy/encoding runs also accept `--policy`, `--encoding`, `--duration`, and `--repeat`. No production runner/helper edits, geometry resizing, controller retuning, or neural saturation-input additions were made.\n'
 text+='\n## CAD scope\n\nNo CAD subset pressure run is claimed here. The adopted two-vas_med geometric sleeve subset still lacks assigned Bowden housing/cable route, fittings/brackets/load validation, effective transmission stiffness/friction, pressure-response identification, and component masses. A future assumed-transmission demonstration must supply and label a complete configuration, explicitly disable or otherwise account for all other actuators, and keep moving-route work parity. Source ledger: `'+str(ROOT.parent/'cad_bpa/EXPERIMENT_RECORD.md')+'`.\n'
 (out/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8')
 dest=out/'source_snapshot';dest.mkdir(exist_ok=True);manifest={}
 for name in ('inverse_feasibility_policy.py','check_inverse_feasibility_policy.py','boundary_policy_runner.py','boundary_policy_report.py','actuator_interface.py','cad_interface.py','native_bpa_force.py'):
  p=Path(__file__).with_name(name);digest=hashlib.sha256(p.read_bytes()).hexdigest();target=dest/(p.stem+'_'+digest[:8]+p.suffix)
  if not target.exists():shutil.copyfile(p,target)
  manifest[str(p)]=dict(sha256=digest,snapshot=str(target))
 (dest/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps(dict(cases=len(rows),fresh=list(audits),record=str(out/'EXPERIMENT_RECORD.md'))))
 index=ROOT/'EXPERIMENT_RECORD.md';begin='<!-- inverse-boundary-policy:start -->';end='<!-- inverse-boundary-policy:end -->'
 addition=begin+'\n## Opt-in inverse feasibility boundary experiment\n\nRecord: `'+str(out/'EXPERIMENT_RECORD.md')+'`. Six frozen-controller cases compare explicit abort with nearest-feasible static pressure targeting for graded/sensory-spiking variants. Abort stops at first invalid reference (16.274/16.560 s moving); nearest still fails the unchanged forward guard at16.292/16.584 s. Both nearest negative outcomes repeat byte-exact in fresh processes, including full preparation and final spiking state. The old zero-sentinel failures remain preserved. Pressure/geometry/neural parameters and actual pressure lag are unchanged; static clipping is insufficient for dynamic route-length feasibility. Neural saturation/reserve feedback is assessed separately, with no hidden motor-neuron edits or new feedback path implemented. CAD configuration remains incomplete.\n'+end
 existing=index.read_text(encoding='utf-8')
 if begin in existing:
  left,tail=existing.split(begin,1);_,right=tail.split(end,1);existing=left+addition+right
 else:existing+='\n\n'+addition+'\n'
 index.write_text(existing,encoding='utf-8')

if __name__=='__main__':main()
