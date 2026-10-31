# Model-specific tendon, damping, and runtime evidence

Read-only evidence audit on easteregg2. No new simulations or manuscript edits.
Repository root for all paths below: `D:\GitHub\Bipedal_Robot`.

## Findings suitable for the results narrative

**W2L supported stepping and its authored muscle physics.** The improved newer
W2L uses the ported AnimatLab `LinearHillMuscle` equations. These already contain
series stiffness, parallel stiffness, and a viscous internal-force state.
They are not MuJoCo's native inelastic-biological-tendon muscle law, a Millard
muscle, or the separate new Thelen-informed prototype. The supported-stepping
results in this package were obtained with these source-specific Hill muscles.
The finite-difference anatomical audit and near-flat supported-stance IK
corrected the reference posture before ankle phase and recruitment tuning;
they did not substitute a new tendon model.

Source: `Code\MuJoCo_SNS\spinal\virtual_walker_20261031\native_kernel.py`
(`hill_step`, `bank_parameters`), and
`Code\MuJoCo_SNS\spinal\physics_fidelity\W2L108_IK_REFERENCE.md`.

**Implemented elastic biological-tendon prototype for the Gait2392 physics
investigation.** A separate 92-unit, Numba/preallocated muscle bank adds an
elastic biological tendon, projected fiber-length state, constant-thickness
pennation, a backward-Euler equilibrium solve, and explicit fiber damping.
It uses the Gait2392 source's Thelen2003 tendon and static force–length
parameters. Across six matched-parameter OpenSim isometric cases, maximum
force difference was **2.1574e-6 N**, and fiber-length difference was
**8.3177e-12 m**. An independent SciPy root trajectory differed by at most
**3.2828e-14 m**. The 92-unit prescribed-motion bank had no bound hits and
maximum equilibrium residual **2.2841e-8 N**. A one-DoF MuJoCo spatial-tendon
coupling remained finite for **4 s**, with residual below **7.21e-10 N**.
Force RMS difference against a 0.5 ms reference decreased from **0.13756 N**
at 2 ms to **0.04632 N** at 1 ms.

Implementation: `Code\MuJoCo_SNS\spinal\physics_fidelity\series_elastic.py`.
Checks: `benchmark.py`, `opensim_static_check.py` in the same directory.
Durable raw evidence:
`Code\MuJoCo_SNS\spinal\physics_fidelity\runs\audit_20261030\RESULTS.json`
and `OPENSIM_STATIC_CHECK.json`.

**Verified damping and contact diagnostics.** A 1 kg, 100 N/m analytic oscillator
check with 10 N·s/m viscous damping gave **0.06441 mm** maximum position error,
monotonically decreasing sampled energy, and a final/initial energy ratio
of **4.1402e-9**. A separate 1 kg, 0.10 m sphere-drop benchmark quantified
the contact tradeoff: baseline peak force/penetration **142.86 N / 8.726 mm**;
the existing less-viscous setting **281.01 N / 3.143 mm**. The audit also found
that the historical contact option called “nonlinear” edited the wrong
impedance index and selected constant impedance; its trace matched an explicit
constant-impedance control exactly. This is a completed implementation audit,
not a calibrated damping coefficient for either walker.

Sources: `physics_fidelity\benchmark.py`, `contact_benchmark.py`,
`runs\audit_20261030\RESULTS.json`, and `CONTACT_RESULTS.json`.

## Runtime claims and their exact workloads

| Workload | Recorded result | Scope |
| --- | --- | --- |
| Original W2L, matched Python and compiled equations | 4 simulated s: **17.20797 s** reference versus **0.969586 s** compiled; **17.7478× speedup**, **4.1255× real time** | **80 neurons, 150 edges**, 24-unit Hill/receptor bank, air configuration, MuJoCo 2.3.7 implicitfast at 2 ms and native 0.2 ms neural substeps; includes diagnostics/post-step forward, excludes loading, JIT/warm-up, saving, and rendering |
| Separate 92-unit compliant biological-tendon bank | 16 simulated s in **3.54658 s**, **4.5114× real time**; matched rigid counterpart **0.100933 s** | Prescribed length/activation, muscle-bank calculation only; compliance costs **35.1378×** the rigid counterpart in this workload |
| Near-flat W2L fixed-force mechanical hold | 2 simulated s in **0.21935 s**, approximately **9.1× real time** | Actual MuJoCo with twelve constant forces, stationary belt, no SNS or dynamic Hill update; a mechanics check, not the controller benchmark |

The first comparison is a matched **end-to-end simulation-loop** measurement,
not a neural-only kernel test. Saved maximum errors were **1.3142e-14** in
qpos, **7.5318e-13 degrees** in reported joint angles, **3.2685e-13 mV** in
neural voltages, and **7.9410e-11 N** in forces. The reference uses Python
versions of the same source equations, not an SNS-Toolbox dense backend or
AnimatLab. It does not establish the same speedup for the newer 108-neuron
ground configuration.

Sources:

- `Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\BENCHMARK_FINAL.json`
  and `benchmark_port.py`.
- `Code\MuJoCo_SNS\spinal\physics_fidelity\runs\audit_20261030\RESULTS.json`.
- `Code\MuJoCo_SNS\spinal\physics_fidelity\runs\canonical_12II\ik_reference\FIXED_FORCE_HOLD.json`.

Use the durable compliant-bank timing above, rather than silently mixing it
with the earlier README's historical 2.47 s / 6.47× real-time timing. The
physics audit was regenerated after scratch outputs disappeared;
`REGENERATED_EVIDENCE.json` records the successful reruns. No new timing was
collected for this evidence note. The earlier 26.94× W2L speedup predates the
additional guide observers; `BENCHMARK_FINAL.json` is the applicable final
matched comparison.

## One scope statement for manuscript use

The elastic-tendon implementation is a tested Thelen-informed prototype with
an approximate force–velocity relation and experimental fiber damping, not a
complete dynamic Thelen/Millard or Hyfydy replication, a calibrated biological
model, or an installed explanation of the W2L stepping improvement. Its
isometric checks establish numerical agreement under matched parameters;
dynamic OpenSim equivalence and whole-walker integration remain separate work.
The contact tests quantify numerical behavior without identifying physiological
damping or validating ground-impact forces.
