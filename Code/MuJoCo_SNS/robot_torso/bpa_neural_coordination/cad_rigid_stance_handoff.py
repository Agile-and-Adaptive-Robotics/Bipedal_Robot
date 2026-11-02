"""Readable handoff for the frozen full-CAD native static audit."""
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
import cad_rigid_stance as C


def main():
    folder = C.OUT
    first = json.loads((folder / 'STATIC_RESULT.json').read_text())
    refined = json.loads((folder / 'REFINEMENT_RESULT.json').read_text())
    capacity = json.loads((folder / 'FORCE_CONE_DIAGNOSTIC.json').read_text())
    rows = []
    for label, result in [('Initial fit', first), ('Same-condition continuation', refined)]:
        rows.append(f"| {label} | {result['qacc_max']:.5f} | "
                    f"{np.linalg.norm(result['COM_acceleration']):.6f} | "
                    f"{result['ground_N']:.5f} | "
                    f"{1000*result['closure_gaps_m'][0]:.6f} / "
                    f"{1000*result['closure_gaps_m'][1]:.6f} | Fail |")
    ledger = refined['native_accounting']
    record = f'''# Full-CAD mechanical statics: frozen numerical variant

Machine: **easteregg2**. This is a new mechanical-only ideal-54-route preload audit. It is not a neural, pressure-actuated, or hardware-calibrated standing result.

## Frozen condition

The previously derived knee-loop equality response is fixed at time constant 0.004 s, damping ratio 1, constant impedance 0.9984069577958117, 500 maximum solver iterations, and 2 ms timestep. The earlier stronger fixture's 115.07 µm peak still misses the original 100 µm design target; this case uses the separately authorized 0.5 mm bound. No further equality-parameter selection was performed.

Both feet have zero heading and the source knee flexion remains protected. Geometry, mass/inertia (75 kg), damping, actuator capacities, and original force floors are unchanged. The other 52 routes retain lower force bounds of 0.05 times their maximum isometric force; the previously free vastus-medialis pair has zero lower bounds. All 54 forces are mechanical preload variables, not fitted neural synapses. The pelvis apparatus provides exactly 0.85 body weight upward at the declared initial anchor, with no initial horizontal force or applied moment.

## Native validation

| Case | Maximum angular acceleration (rad/s²) | COM acceleration norm (m/s²) | Ground load (N) | R/L loop gaps (mm) | Static decision |
| --- | ---: | ---: | ---: | --- | --- |
{chr(10).join(rows)}

Required ground force is 110.3625 N. Registered gates are maximum acceleration ≤0.05 rad/s², COM acceleration norm ≤0.01 m/s², ground-load error ≤0.1 N, each foot >1 N, each heel/toe channel >0.5 N, loop gap ≤0.5 mm, knee flexion ≥0.02 rad, heading magnitude ≤10°, and mesh minimum height between −4 and 0 mm.

The continuation's four heel/toe channels are {refined['contact_loads_N']} N. Both feet, all contact channels, loop gaps, knees, heading, mesh height, physical constraint rank, and EOM residual pass. Acceleration, COM acceleration, and total ground load fail. Native EOM maximum residual is {ledger['native_EOM_residual_max']:.9g}; the captured forward/inverse generalized-force discrepancy is {ledger['forward_inverse_force_residual_max']:.9g}. Each knee's physical three-hinge Jacobian has rank two, totaling four. The raw full equality Jacobian has rank {ledger['raw_full_rank']} at this slightly deflected pose; no actual rows were removed or replacement constraint forces inserted.

The first Powell fit used 1,500 evaluations. A predeclared finite-difference consistency probe then found relative Jacobian disagreement 8.6657, so the same-condition continuation used the bounded derivative-free branch, also 1,500 Powell evaluations (3,044 inverse calls including derivative probes and the objective's two residual evaluations). Neither reached solver convergence. The original result and both exact preload states are preserved. Failure is not a proof of global stance infeasibility.

## Contact-force capacity diagnosis

At the saved refined pose, an ideal rank-four contact-force allocation is feasible when frictional horizontal forces are allowed (maximum equality residual {capacity['full_friction']['residual_max']:.9g}), but infeasible when contact forces are restricted to vertical. Its minimum sum of absolute horizontal force components is {capacity['minimum_total_absolute_shear_N']:.6f} N. The actual native contacts already produce {capacity['native_total_absolute_shear_N']:.6f} N by the same absolute-component measure. Therefore the diagnosis does **not** establish that friction is absent; the native force distribution is unbalanced.

A diagnostic-only zero-floor mathematical force cone restores vertical-only feasibility (residual {capacity['normal_only_zero_floor_DIAGNOSTIC']['residual_max']:.9g}). This identifies sensitivity to the imposed five-percent-of-MIF force-floor policy. It does not alter that policy in the physical model, demonstrate a controller fix, or equate five-percent force with five-percent physiological activation.

## Decision and preserved scope

**No 5 s hold, fresh physical repeat, or full-robot 1 ms timestep run was launched**, because native static gates failed. There was no rig/controller/BPA geometry retuning. Prior pressure-domain failures and prior neural/pressure results remain separate and unchanged. A free ideal contact-force allocation is not a native static or walking result.

## Reproduction

Python: `D:\\Anaconda\\envs\\myo\\python.exe`.

Source folder: `D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\robot_torso\\bpa_neural_coordination`.

`cad_rigid_stance.py --solve` creates the first bounded fit; `cad_rigid_refine.py` performs the single registered continuation; `cad_rigid_capacity.py` performs saved-pose capacity diagnostics; `cad_rigid_stance_handoff.py` regenerates this handoff only. Existing outputs are the record, not instructions to rerun completed fits. `PRELOAD_STATE.npz` and `REFINED_PRELOAD_STATE.npz` contain complete static qpos/qvel/qacc, all 54 force preloads, initialization variables, and the explicit root anchor. These are static artifacts, with no neural, pressure, or encoder states.

Authoritative results are `STATIC_RESULT.json`, `REFINEMENT_RESULT.json`, `FORCE_CONE_DIAGNOSTIC.json`, and their preregistrations in this folder. Source snapshots and hashes accompany the handoff.
'''
    (folder / 'EXPERIMENT_RECORD.md').write_text(record, encoding='utf-8')
    snapshot = folder / 'source_snapshot'
    snapshot.mkdir(exist_ok=True)
    sources = [Path(__file__), Path(C.__file__)] + [Path(__file__).with_name(name) for name in
               ('cad_rigid_refine.py', 'cad_rigid_capacity.py', 'rigid_fourbar_numerics.py',
                'cad_ideal54_stance.py', 'cad_compliant_stance.py', 'cad_stance_repair.py')]
    hashes = {}
    for source in sources:
        shutil.copy2(source, snapshot / source.name)
        hashes[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
    (folder / 'HANDOFF_SOURCE_HASHES.json').write_text(json.dumps(hashes, indent=2))
    decision = dict(accepted=False, native_static_gate=refined['checks'], hold5s_launched=False,
                    fresh_repeat_launched=False, full_robot_dt1ms_launched=False,
                    reason='Native acceleration, COM acceleration, and total ground load fail.',
                    frozen_numerics=True, original_fixture_100um_target_passed=False,
                    geometry_or_mass_or_forcefloor_or_neural_change=False)
    (folder / 'DECISION.json').write_text(json.dumps(decision, indent=2))
    print(folder / 'EXPERIMENT_RECORD.md')


if __name__ == '__main__':
    main()
