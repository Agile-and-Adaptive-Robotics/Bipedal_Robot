# Four-cell RF/hamstring PF: training and physical validation

Session 1, easteregg2, 2026-10-29. Human Gait2392, 92 muscle actuators.
This is separate from CHATGPT_REPORT.md. No dissertation source was changed.

## Outcome

Three fits of the four-cell version followed Hunt's sequence and were physically tested.
No variant passed the complete corrected curriculum. Assisted standing passes;
ground stepping remains unsuccessful. Stage 5 was withheld because verified
alternating load-bearing steps were not obtained. Later-stage runs after an
earlier failure are diagnostics, not curriculum graduation.

The RF/HAM module has two graded PF integrators and two reciprocal inhibitory
relays per leg. It is not another intrinsic NaP oscillator. Sartorius and
gracilis retain separate function-specific modules and independent MN outputs.
The conditional census is 1,025 neurons in deafferented air, 1,027 with the
commissural candidate, and 1,035/424 inputs with the added PT branch.

## Hunt training

A/B were revalidated before fitting C: intrinsic period 2.466 s; isolated
sensory-driven right/left periods 1.236/1.212 s. The coupled expected-afferent
network gives 1.264/1.246 s. These are timing checks, not physical gait results.

C fits joint PF, RF/HAM integration, and MN output sequentially, without old D
weights contributing to the feedforward fit. D then fits Ia/Ib and optional
IaIN refinement with C frozen. All motor targets are inverse-model estimates,
not EMG. The 80% C target does not establish an additive physiological 80/20 split.

The first fit uses expected sensory-driven timing. The second fits C jointly
at that cadence and at the slower, deafferented intrinsic cadence. The third
also fits the native four-cell dynamics, including inhibitory conductances,
then refits MN output and D. Its small preallocated model reproduces the full
network's recorded module dynamics to 4.44e-16 mV before optimization.

| Fit | C+D activation RMSE | Stage-1 phase error | Stage-1 smaller knee ROM (deg) | Stage-2 phase error |
|---|---:|---:|---:|---:|
| Single cadence | 0.12700 | 0.123 | 37.68 | 0.201 |
| Two cadences | 0.13502 | 0.074 | 16.85 | 0.225 |
| Two cadences + native paired fit | 0.13500 | 0.067 | 17.33 | 0.359 |

The table uses nominal, unscaled physical cases. Stage 1 in the first fit
fails coverage/frontal checks. The later fits improve deafferented phase and
frontal behavior, but the left knee remains below the 20-degree ROM gate.
Their afferented air motion still fails phase consistency. Smaller supervised
error alone is not a sufficient controller-selection objective.

A separate single-cadence calibration probe passes the air motion/frontal
screen: 4/4 knee cycles, phase error 0.1437, period 2.166 s, and no sampled
crossing. It changes commissural/reflex and motor gains after the fit, has no
all-synapse robustness certificate, and is not a curriculum winner.

## Actual ground contact and support

Ground diagnostics include static harnesses, gradual lowering by 40–120 mm,
explicit constant-speed towing, lateral restraint reduction, intact stance
output, and timing/motor-reflex ablations. Towing is mechanical assistance;
it is not self-propelled locomotion. None supplies a periodic stride clock.

The strongest loaded diagnostic from the native paired fit transfers
55.81% BW through the feet and 44.05% through the vertical harness.
Both feet remain in contact. This is a load-sharing result, not walking or a
successful support-weaning rung. Removing motor reflexes and timing inputs
also fails to produce sustained alternating ground contacts.

Recorded bilateral RG activity continues bursting on the ground. The initial
paired PF module is biased toward right RF and left HAM, and the feet fail
to unload. Thus continuous stance contact must not be described as proof that
the RG itself stopped. The motor/force transfer and phase coordination remain
unresolved; lower neural replay error has not solved them.

The corrected physical gate now uses once-counted GRFs, 0.05/0.025 BW
hysteresis, at least 80 ms preceding unloading and 100 ms loaded stance,
three contacts per leg, alternating contact phases, appropriate duty, physical
motion, frontal clearance, and pelvis tilt. Foot-body-origin heights are not
sole clearance; full foot-edge clearance/self-collision is not certified.
Legacy duplicated heel/toe contact totals remain diagnostics, not acceptance
evidence. Synthetic checks accept alternating loading and reject synchronous
loading, permanent double stance, and airborne motion. These synthetic cases
are tests of the gate, not successful walker trials.

A fresh replay with the corrected gate is bit-identical in qfull, control,
and recorded neural traces to the earlier loaded diagnostic. The new gate
does not change the model dynamics or rescue the failed result.

## Acceleration and records

The first complete fit takes 158.71 s; the two-cadence fit takes 383.13 s.
The native paired optimization performs 768 module trials in
10.91 s, excluding MuJoCo physics and the subsequent MN/D refit.
Sparse SNS updates, preallocated arrays, cached edge indices, and vectorized
afferents are used throughout. The scalar/vector afferent encoders agree exactly
in the sampled check. A 500-update fixed-conductance neural comparison takes
0.245 s accelerated versus 24.534 s native (100.2x), with maximum voltage disagreement 3.55e-15 mV.
This excludes physics, changing conductances, setup, and logging, and was measured
during concurrent physical work. It is not a 100x full-simulation claim.
Earlier matched full-model measurements remain about 6.7–9.2x end to end.

Raw weights, traces, logs, source snapshots, hashes, and trial configurations
are under the ignored runs/fourcell_* directories. INPUT_PROVENANCE.json and
CURATED_SUMMARY.json accompany the native-pair campaign. Raw convenience
deafferented replay JSON used a nominal-period warmup; the curated summary
also recomputes that score after four intrinsic periods. No raw result was
silently overwritten to change its outcome. Earlier candidates are preserved.

The current viewer is network_figures/fourcell_dynamic/index.html; it includes
the actual four-cell subnetwork. The grounded diagnostic figure uses the
single-cadence full-feedback trial and is explicitly labeled with that scope.
No harness-free walking is demonstrated.

Latest-conductance verification: 0.2530 s sparse versus 14.3064 s native over 500 updates (56.55x); maximum voltage disagreement 4.34e-19 mV. The physical queue had completed. The encoder comparison remains exact. This has the same neural-only scope; timings vary across runs.
