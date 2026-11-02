"""Frozen two-relay spiking transfer: gated fresh repeats and compact handoff."""
import json,hashlib,subprocess,sys,shutil
from pathlib import Path
import numpy as np
import spiking_recruitment_experiment as S
D=S.DEST;HERE=Path(__file__).resolve().parent
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def j(path):return json.loads(path.read_text())
def audit_sources():
 S.E.registration();pre=j(D/'PHYSICAL_PREFLIGHT.json');assert pre['passed'] and pre['source_sha256']==digest(Path(S.__file__))
 files=[Path(S.__file__),D/'OFFLINE_CHECKS.json',D/'PHYSICAL_PREFLIGHT.json',S.BASE/'pressure/SELECTED.json',S.BASE/'PREREGISTRATION.json',S.BASE/'source_snapshot/MANIFEST.json'];manifest={str(p):digest(p) for p in files};path=D/'FROZEN_SOURCE_MANIFEST.json'
 if path.exists():assert j(path)==manifest,'Frozen spiking or inherited source changed'
 else:
  path.write_text(json.dumps(manifest,indent=2));snap=D/'source_snapshot';snap.mkdir(exist_ok=True);shutil.copyfile(S.__file__,snap/'spiking_recruitment_experiment.py')
 return manifest
def run(name):
 folder=D/name;folder.mkdir(exist_ok=True);summary=folder/'CASE_SUMMARY.json'
 if not summary.exists():
  with (folder/'RUN.log').open('w',encoding='utf-8') as f:p=subprocess.run([sys.executable,str(Path(S.__file__)),'--case',name],stdout=f,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
  if p.returncode:raise RuntimeError(str(folder/'RUN.log'))
 return j(summary)
def compare(a,b):
 checks={}
 for path in sorted(a.glob('*.npz')):
  other=b/path.name;assert other.exists(),str(other)
  with np.load(path) as x,np.load(other) as y:
   assert set(x.files)==set(y.files)
   checks[path.name]={k:bool(np.array_equal(x[k],y[k])) for k in x.files}
 assert checks and all(ok for row in checks.values() for ok in row.values()),'Fresh array/state mismatch'
 return checks
def main():
 audit_sources();initial={name:j(D/name/'CASE_SUMMARY.json') for name in ('nominal_initial','strong_initial')}
 eligible=[n for n,r in initial.items() if r['accepted'] and r['neural_restore']];fresh={n.replace('_initial','_fresh'):run(n.replace('_initial','_fresh')) for n in eligible};assert all(r['accepted'] and r['neural_restore'] for r in fresh.values())
 parity={n.removesuffix('_initial'):compare(D/n,D/n.replace('_initial','_fresh')) for n in eligible};status=dict(initial=initial,fresh=fresh,all_fresh_arrays_exact=bool(parity),parity=parity,scope='Only two newly added error relays spiking; failed initial cases are NOT repeated. No physical tuning or independent validation. Prior graded3/4 result is not a hybrid validation result.')
 (D/'FRESH_REPEAT_AUDIT.json').write_text(json.dumps(status,indent=2))
 with np.load(D/'strong_initial/w2lvar_case.final_neural_state.npz') as z:spikes=z['recruit_lif_total'].tolist()
 strong=initial['strong_initial'];nom=initial['nominal_initial'];off=j(S.BASE/'pressure/off_train_strong/CASE_SUMMARY.json');cal=j(D/'OFFLINE_CHECKS.json')
 record=f'''# Two-error-relay spiking mirror

Machine: **easteregg2**. The frozen spiking mirror passes the nominal 8 s test, with exact fresh-process repetition of every saved array and state. It **fails** the known strong 20 s study at 19.374 s released, before the required 19.4 s. The unchanged pectineus sleeve-domain guard stops the simulation. The failed condition was not repeated or tuned. The graded relay version remains the successful strong-case configuration. No prior validation-set tuning was performed.

## What changed

Only the two recently added pectineus error relays are replaced by 128 genuine reset/refractory LIF cells per side. Original 888 normalized-w2lvar neurons and two valve relays remain graded: 256 spiking and 890 graded cells. This does not make the whole SNS spiking. The decoded population signal supplies the same bounded 5 nA positive MN conductance current. The nominal pressure interface, valve10 architecture, force map, pressure lag, sensor delay/noise, geometry, knee limits, and apparatus remain frozen.

LIF integration uses 0.1 ms substeps, 5 ms membrane decay, 1 ms refractory duration, deterministic initial membrane phases, and a mirrored 0.8–1.2 frequency spread around the requested population mean (maximum mean 200 Hz; individual maximum 240 Hz). Spike-count filtering uses {-0.002/np.log(2/3):.9f} s to match the graded relay's discrete decay. Decoder gain {cal['decoder_gain']:.12f} was fitted only to constant electrical inputs, with no physical outcomes involved. Outputs retain the same previous-sample 2 ms hold. Membrane, refractory, decode filter, drive cache, counts, histories, and explicit absence of RNG are serialized and restored.

The first identical-frequency electrical prototype synchronized and failed its dynamic mapping gate; OFFLINE_IDENTICAL_RATE_FAILED.json preserves that failure. The heterogeneous version has maximum static error {max(r['absolute_error'] for r in cal['static_rows']):.6f} full scale and dynamic RMSE {cal['dynamic_rmse_fraction']:.6f} full scale. Complete 64-step continuation is byte-exact; zero-input quiescence and zero-gain original forward outputs pass. Preparation matches all inherited state fields, allowing only declared resting LIF states. This is not a separate full-duration zero-gain physical rerun.

## Physical results

| Condition | Study / released duration | Result | Fresh arrays and states |
| --- | --- | --- | --- |
| Nominal | 8 / 7.4 s | Pass | All exact |
| Known strong | Requested 20 / 19.4 s; stopped 19.374 s released | Fail, unchanged domain guard | Not repeated |
| Frozen valve-only known strong | Domain failure at {off['moving_s']:.3f} s | Fail | Previously retained graded control |

At failure, pect_r exceeds its pressure-supported maximum route length by 3.98153e-05 mm at 71.5256 kPa. That small excess is preserved, not excused by rounding or a relaxed tolerance. The failure tests the declared actuator domain rather than demonstrating a physical fall. Strong-case upward/downward apparatus forces before failure are {100*strong['upward_bw']:.4f}%/{100*strong['downward_bw']:.5f}% BW; lateral RMS/peak forces are {100*strong['lateral_rig_rms_bw']:.4f}%/{100*strong['lateral_rig_peak_bw']:.4f}% BW. Peak pressure is {strong['actual_pressure_peak_kpa']:.3f} kPa; maximum added MN current is {strong['max_added_MN_current_nA']:.4f} nA. Visual-foot clearance minimum is {1000*strong['mesh_min_m']:.3f} mm. No knee hyperextension gate was relaxed. Actual spike counts are {spikes} across the two populations. Strong lateral and rotational guidance remains; this is not independent standing, walking, calibrated CAD, or a curriculum result.

This is a nominal/known-challenge transfer check, not independent uncertainty validation. The graded study's 3/4 validation and weaker-map failure remain frozen; no spiking validation claim is made. The valve-only control is a previously frozen graded control and does not isolate all possible spike encoders or prove whole-SNS necessity.

## Files and commands

FROZEN_SOURCE_MANIFEST.json records hashes checked before fresh repetitions; source_snapshot contains the exact spiking module. PHYSICAL_PREFLIGHT.json records electrical/current/default checks. FRESH_REPEAT_AUDIT.json contains every array/state comparison. Each case contains complete pressure, current, event-count, mechanical, and signed apparatus evidence. Original independent archives remain separate.

```powershell
& 'D:\\Anaconda\\envs\\myo\\python.exe' '{HERE / 'spiking_recruitment_followup.py'}'
```

The command reuses completed cases and audits hashes. A new explicitly named `nominal_*` or `strong_*` case can be run using spiking_recruitment_experiment.py --case NAME with the same full Python path. No manuscripts or accepted production checkpoint were changed.
'''
 (D/'EXPERIMENT_RECORD.md').write_text(record,encoding='utf-8');print(json.dumps(dict(initial_acceptance={n:r['accepted'] for n,r in initial.items()},fresh_acceptance={n:r['accepted'] for n,r in fresh.items()},fresh_all_arrays_exact=bool(parity),spike_counts=spikes,record=str(D/'EXPERIMENT_RECORD.md')),indent=2))
if __name__=='__main__':main()
