"""Document source activation semantics and the distinct ideal-force fixture."""
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cad_activation_semantics as A


def main():
    folder = A.OUT
    audit = json.loads((folder/'ACTIVATION_SEMANTICS_AUDIT.json').read_text())
    result = json.loads((folder/'STATIC_RESULT.json').read_text())
    previous = json.loads((folder.parent/'REFINEMENT_RESULT.json').read_text())
    com = float(np.linalg.norm(result['COM_acceleration']))
    ground_error = abs(result['ground_N']-result['target_ground_N'])
    ledger = result['native_accounting']
    text = f'''# CAD activation semantics and zero-floor mechanical capacity

Machine: **easteregg2**. This explicitly new condition is mechanical-only. No SNS, BPA pressure, physiological coactivation, or standing success is represented.

## Source-specific meaning of five percent

All 54 actuators in `D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\robot_torso\\ah001_neural.xml` are first-order general/filter actuators with 30 ms time constant, fixed negative peak-force gain, and zero bias. Their normalized internal state a follows da/dt = (u−a)/0.03. Positive tensile-reference magnitude is **T = F0 a**, independent of route length or shortening velocity. There is no passive-fiber or compliant-tendon law in these native reference channels.

Direct native actuation checks at normalized state 0, 0.05, and 0.10, using both the exact closed geometry and the saved refined geometry, reproduce this map with maximum error zero. Five percent normalized state therefore means exactly five percent peak reference tension **within this ideal model**. It is not a verified mapping from physiological five-percent activation, EMG, or BPA pressure. Human-named peak capacities alone do not identify force-length/velocity or passive-law parameters for these custom CAD routes; assigning human Hill factors would add unsupported assumptions.

The existing 52-route floor sums to **{audit['sum_internal_reference_tensions_N']:.1f} N of internal reference tension**, with largest single floor {audit['largest_single_floor_N']:.2f} N. This sum is not an external support force. Imposing this floor is also distinct from tonic neural excitation or maintaining antagonistic neural activity. The original floor policy and all its negative results remain unchanged.

## Explicit new mechanical definition

All 54 force variables are allowed in [0,F0], with zero modeled passive force. This is a mathematical and mechanical force-feasibility fixture, **not a zero-floor SNS case or a five-percent-coactivation controller**. The ideal model does not provide a defensible source-specific physiological five-percent tension bound for the other 52 routes.

The new fit starts from the previous saved refined initialization. Geometry, mass/inertia (75 kg), damping, source actuator capacities, the 85%-BW apparatus force, and the existing grounded zero-heading pose-family bounds are unchanged. Equality numerics remain frozen at 4 ms, constant impedance 0.9984069577958117, 500 maximum iterations, and 2 ms timestep. The fixture's prior 100 µm target miss remains preserved. No sleeve, pressure, geometry, floor policy in an existing run, or neural parameter was changed.

## Actual native checks

| Measure | Previous 52×5%-of-MIF floor | New zero-floor ideal fixture | Required |
| --- | ---: | ---: | ---: |
| Maximum angular acceleration (rad/s²) | {previous['qacc_max']:.6f} | {result['qacc_max']:.6f} | ≤0.05 |
| COM acceleration norm (m/s²) | {np.linalg.norm(previous['COM_acceleration']):.6f} | {com:.6f} | ≤0.01 |
| Actual ground force (N) | {previous['ground_N']:.6f} | {result['ground_N']:.6f} | 110.3625 ±0.1 |
| Right / left loop gap (mm) | {1000*previous['closure_gaps_m'][0]:.6f} / {1000*previous['closure_gaps_m'][1]:.6f} | {1000*result['loop_gaps_m'][0]:.6f} / {1000*result['loop_gaps_m'][1]:.6f} | ≤0.5 each |

The four heel/toe loads are {result['contact_loads_N']} N. Contact channels, loop gaps, knees, zero heading, mesh height, native EOM, and physical rank four pass. Ground-load error is {ground_error:.6f} N. Acceleration, COM acceleration, and ground-load gates fail. EOM maximum residual is {ledger['native_EOM_residual_max']:.9g}; captured forward/inverse generalized-force discrepancy is {ledger['forward_inverse_force_residual_max']:.9g}. No actual constraint rows were removed. The raw full equality rank is {ledger['raw_full_rank']}, with physical knee-hinge rank two per side.

The single bounded Powell case uses {result['nfev']} evaluations and {result['wall_s']:.3f} s of solver wall time, ending at its evaluation budget rather than convergence. It improves this local native equilibrium approximation materially, but does not prove global feasibility or infeasibility. No second fit, 5 s hold, fresh physical repeat, or full-robot timestep run was launched.

## Durable artifacts and entry points

`ACTIVATION_SEMANTICS_AUDIT.json` contains all 54 reference maps at two saved poses and the floor generalized-force vector. `PREREGISTRATION.json` fixes the one new condition, gates, and source hashes. `STATIC_RESULT.json` contains actual native checks and full force/EOM accounting. `PRELOAD_STATE.npz` contains qpos, qvel, actual qacc, all 54 preloads, pose variables, and the explicit root anchor. None contains neural/encoder/pressure state because none was simulated. `STATIC_FLOOR_COMPARISON.png/.pdf/.svg` is a diagnostic figure, not a walking result.

Python: `D:\\Anaconda\\envs\\myo\\python.exe`. Source folder: `D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\robot_torso\\bpa_neural_coordination`. Audit entry: `cad_activation_semantics.py`; single-case entry: `cad_zero_floor_static.py`; documentation/figure entry: `cad_activation_handoff.py`. Completed numerical cases need not be rerun. Source snapshots and hashes accompany this record. Original XML, robot runner, dynamic routes, neural weights, and manuscript files were not edited.
'''
    (folder/'EXPERIMENT_RECORD.md').write_text(text, encoding='utf-8')
    plt.rcParams.update({'font.family':'Arial','font.size':10})
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.35))
    measures = [(previous['qacc_max'], result['qacc_max'], .05, 'Peak angular acceleration', 'rad/s²'),
                (np.linalg.norm(previous['COM_acceleration']), com, .01, 'COM acceleration norm', 'm/s²'),
                (abs(previous['ground_N']-previous['target_ground_N']), ground_error, .1, 'Ground-force error', 'N')]
    for ax, (old, new, gate, title, units) in zip(axes, measures):
        ax.bar([0,1], [old,new], color=['#4477AA','#EE6677'], width=.58)
        ax.axhline(gate, color='#222222', linestyle='--', linewidth=1.2, label='Acceptance limit')
        ax.set_yscale('log'); ax.set_ylim(gate*.5, max(old,new)*4)
        ax.set_xticks([0,1], ['5% force floor', 'Zero-floor fixture'])
        ax.set_title(title, fontsize=10); ax.set_ylabel(units)
        ax.grid(axis='y', alpha=.18); ax.set_axisbelow(True)
        for i,value in enumerate([old,new]): ax.text(i,value*1.12,f'{value:.3g}',ha='center',va='bottom')
        ax.spines[['right','top']].set_visible(False)
    axes[0].legend(frameon=False,fontsize=9,loc='upper right')
    fig.suptitle('CAD static capacity: both cases fail native equilibrium gates', fontsize=12)
    fig.tight_layout(rect=[0,0,1,.93])
    for suffix in ('png','pdf','svg'): fig.savefig(folder/('STATIC_FLOOR_COMPARISON.'+suffix),dpi=220,bbox_inches='tight')
    plt.close(fig)
    (folder/'STATIC_FLOOR_COMPARISON_alt.txt').write_text('Three logarithmic bar charts compare unchanged ideal54 force-floor and zero-floor mechanical cases. Acceleration, COM acceleration, and ground-force error improve with zero floors, but all remain above declared dashed acceptance limits. No neural or pressure actuation was simulated.')
    snapshot = folder/'source_snapshot'; snapshot.mkdir(exist_ok=True)
    files = [Path(__file__),Path(A.__file__),Path(__file__).with_name('cad_zero_floor_static.py')]
    hashes={}
    for source in files:
        shutil.copy2(source,snapshot/source.name);hashes[str(source)]=hashlib.sha256(source.read_bytes()).hexdigest()
    (folder/'HANDOFF_SOURCE_HASHES.json').write_text(json.dumps(hashes,indent=2))
    print(folder/'EXPERIMENT_RECORD.md')


if __name__ == '__main__': main()
