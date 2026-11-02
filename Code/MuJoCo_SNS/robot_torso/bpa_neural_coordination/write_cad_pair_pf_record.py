"""Current-path audit handoff; no candidate or model changes."""
import json,hashlib,shutil
from pathlib import Path
import cad_pair_pf_audit as A
D=A.OUT;HERE=Path(__file__).resolve().parent
def main():
 current=json.loads((D/'CURRENT_AUDIT.json').read_text());stance=json.loads((D/'STANCE_PATH_AUDIT.json').read_text());fit=stance['conditional_PF_and_IBEXC_common_gain'];lines=f'''# CAD vastus-pair neural recruitment: source audit

Machine: **easteregg2**. No motor-neuron clamp has yet been removed in a live run. The requested PF-only fit fails its offline feasibility check because a second, strong existing excitatory pathway remains. The successful mixed-pressure hybrid baseline and all prior graded/spiking evidence are unchanged.

## Actual supported record

A single 2 s instrumentation repeat records every neural input and all 262 before/after voltages, aligned to 2 ms physical ticks. It uses the same existing meaningful 5% hybrid teacher: two assumed pressure BPAs, 52 reference-force routes, and approximately 92% BW apparatus support. This is a teaching signal, not measured EMG. Its normalized pair commands are 0.05 and MN voltages are 0.25 mV. Replaying those literal commands and recorded lengths through the unchanged pressure-history/series-force model reproduces actual pressures and tensions byte-exactly.

## Why PF-only fitting cannot work

At the observed target voltage, PF_E1/E2 contribute 1.550/2.325 nA per MN. The IBEXC_knee_ext pathway contributes approximately 3.55 nA; POSTURE contributes only approximately 0.143 nA. Leak is −0.25 nA, intrinsic bias is zero, and the remaining sensory currents are smaller. The current clamp cancels this strong net excitation. Without the clamp, instantaneous equilibrium is approximately 3.94–4.02 mV, rather than 0.25 mV.

A shared, nonnegative PF-phase fit drives its sensitive phase gains to zero yet leaves {current['fit_RMS_current_residual_nA']:.3f} nA RMS excess. Removing POSTURE as well still leaves approximately 3.33 nA. Those are structural range failures; no PF gain sweep or physical trial is justified by them.

## What IBEXC actually implements

The source adds an ordinary 1.0 µS excitatory RG_E→IBEXC synapse, excitatory Ib inputs, and an excitatory IBEXC→MN output. Although its comment calls this stance-gated, the inputs are additive, not a multiplicative load gate. At these observed states, mean RG input is approximately 3.809 nA, whereas summed actual Ib input is approximately 0.0060–0.0071 nA, only 0.157–0.187% of positive input. A local zero-Ib, RG-only equilibrium remains approximately 3.815 mV. This local calculation is not a closed-loop force-feedback ablation. No source rule has been repaired or renamed globally.

## Conditional next option, not a tested controller

An additional read-only diagnosis scales existing PF and IBEXC output conductances to these two MNs with one shared bilateral factor, preserving all their ratios and leaving POSTURE, leak, intrinsic bias, and other sensory weights unchanged. Its fitted factor is {fit:.12f}; current residual is {stance['conditional_current_RMS_residual_nA']:.6f} nA RMS. An isolated-neuron replay driven by frozen measured presynaptic voltages gives MN ranges {stance['isolated_MN_voltage_range_mV']} mV. The corresponding fixed-geometry pressure replay is domain-valid, peaks at {stance['conditional_pressure_peak_kpa']:.3f} kPa, and has force RMS errors {stance['conditional_force_RMS_error_N']} N versus the teacher trajectory.

This extends beyond the requested PF-only output fit and has not been applied to source weights, a free network, or live mechanics. It is a possible standing-specific recruitment normalization, not a successful controller or a literature rule. The parent has been asked to decide whether to include the existing RG/Ib-integrator output in that fit. A live branch would also have to remove direct pair POST-clamp input and target-based pair initialization, verify literal SNS state/default behavior, and establish a short nominal gate before any feedback or longer-hold tests.

For that future freed pair, 5% would be the offline teaching target, not an imposed activation or force floor. Other 52 routes may retain the explicit original hybrid context, preventing any whole-model autonomous-control claim. No new all-BPA assignment or hardware identification is inferred.

## Files

TEACHER_NEURAL.npz contains the measured electrical record. CURRENT_AUDIT.json/.npz contains current decomposition and PF-only infeasibility. TONIC_AND_PF_FEASIBILITY.json lists every actual named incoming synapse. STANCE_PATH_AUDIT.json and CONDITIONAL_CENTRAL_REPLAY.npz contain the additive RG/Ib attribution and conditional replay. AUDIT_PROTOCOL.json records the original audit scope; source_snapshot preserves the analysis files. Original CAD XML, neural sources, manuscript, and accepted checkpoints are untouched.
'''
 (D/'EXPERIMENT_RECORD.md').write_text(lines,encoding='utf-8');snap=D/'source_snapshot';snap.mkdir(exist_ok=True);manifest={}
 for name in ('cad_pair_pf_audit.py','cad_pair_pf_feasibility.py','cad_pair_stance_path_audit.py','write_cad_pair_pf_record.py'):
  p=HERE/name;dest=snap/name
  if not dest.exists():shutil.copyfile(p,dest)
  manifest[str(p)]=dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),snapshot=str(dest))
 (snap/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));print(str(D/'EXPERIMENT_RECORD.md'))
if __name__=='__main__':main()
