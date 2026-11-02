# Integrated Hunt neural design review

Machine: easteregg2. Role: Gait neural pathway composition reviewer.
The runnable candidate uses the current Hunt-PT architecture with 1,033 neurons,
the actual four-cell RF/hamstring PF circuit per leg, and the same saved weight
files as `aligned_support_ladder_v1/packaged_configuration.json`.
Its numeric `rg_tonic_current` is **1.2**; the older folder name containing
`tonic_0p02` is not the numeric configuration. No earlier air or ground pass is
transferred to these new candidates.

## Applied composition

`integrated_hunt_design.py` builds the actual controller and consumes:

1. The optional graded local-reflex design, including restoration of all 92
   homonymous II paths into the authoritative runtime reflex parameter set.
2. An explicit HUNT_AF input/unit adapter and fast relay design, with the
   timing-strength alternatives below kept separate.
3. Existing IBEXC neurons configured as thresholded phase/load coincidence
   gates. The declared full-stance reference is 0.225 nA at F/Fmax = 1;
   autogenic Ib routes retain their separate function.
4. The real Peterka/Ting-inspired task-balance object, retained for actual model
   binding and tested through its original calibration method.

No NaP intrinsic parameters, PF shape, muscle dimensions, contact geometry,
force law, tendon compliance, damping, timestep, or harness mechanics are changed
by this composition. Sartorius and gracilis remain separate from the RF/hamstring
PF pair. An allocated interneuron is not automatically an active motor pathway;
the saved graph and command tests distinguish these cases.

The requested heel hypothesis B is **not applied**. The saved Hunt controller
already has a direct heel-to-contralateral-RG-F excitatory connection with
conductance 6, while inherited heel-to-InE/InF paths are masked. Installing B
would reactivate masked paths and stack a different contact policy. The existing
choice is preserved, and the actual named contact connections are serialized.

## Two named timing alternatives

| Alternative | E-pool reference current | F-pool reference current | Interpretation |
| --- | ---: | ---: | --- |
| `function_equivalent` | 0.500991937637077 nA | 0.0518116865856527 nA | Preserves the saved unsaturated source input-output strength; primary comparison baseline |
| `canonical_right_role_budget` | 0.5 nA | 0.5 nA | Explicit balanced engineering budgets; separately changed F strength |

These are reference currents at RG voltage 2.5 mV and the declared coactive
feature reference, not global current clamps. Both mirror the selected source
function bilaterally. The second increases the F reference substantially;
equal numerical budgets are not presumed physiologically better or required
for FSA correctness.

Both use a relay capacitance/leak time constant of 3.8962288367 ms at the unchanged
2 ms neural update. The discrete encoder phase allocation is 10 degrees at the
declared 7.155635 Hz engineering band. PF and intrinsic RG clocks are retained.
This is an encoder allocation, not a certification of whole-loop bandwidth.

## Verified checks

For each timing alternative:

- 128 actual accelerated neural steps across commands 0, 0.5, 1, and 0, with
  final authoritative gain, timing-window, coincidence, and cache checks.
- 64 fixed-gain comparisons among zero conductance, connection removal/rebuild,
  and the dense reference: maximum voltage difference **0.0 mV**.
- Actual isolated motor-current subtraction error below 4e-15 nA.
- Exact parameter reconstruction from the saved profile and gain checkpoint,
  followed by an independent fresh-process final-graph comparison.
- Exact `None` behavior and a fresh default constructor after candidate use.
  The actual default `physical_eval.build` constructor also matches the source
  AST with only the new optional hook removed. No whole-body default trajectory
  parity is inferred from this constructor check.

The runtime keeps zero-valued allocated slots that later legitimate command
changes may use. The deletion/rebuild comparison concerns frozen final gains
with finite inputs; it does not delete live fitted parameters from the model.

## Stage-specific protocol handling

The exact serialized design was consumed through actual `physical_eval.build`
setup for stages 1, 2, and 3 without body integration.

| Stage | Constructed neurons | Timing afferents | PT object | Protocol treatment |
| --- | ---: | --- | --- | --- |
| 1 | 1,025 | Disabled | Omitted | Reflex outputs remain zero; active timing/Ib/afferent/PT checks explicitly omitted |
| 2 | 1,033 | Enabled | Retained | Full composition installed; mechanical binding follows normal setup |
| 3 | 1,033 | Enabled | Retained | Full composition installed; source protocol separately suppresses tonic RG current |

Removing PT changes sparse slot indices. Gain checkpoints therefore bind by
exact source/target names, with only the explicitly omitted PT endpoints allowed
in stage 1. They do not assume that indices match the 1,033-neuron graph. No
deafferented motor pathway is silently reactivated to make an electrical check
pass. The profile specification remains fixed while protocol masks are recorded.

## Actual bound-PT preflight

Each alternative was run through the original `physical_eval` model/pose/harness
setup and original `HuntNetwork.bind_mujoco`. The process stopped immediately
after binding; `mujoco.mj_step` was replaced by a rejecting guard. Results:

- Original fixed-pose PT recruitment calibration completed, with bilateral
  recruitment verified and all four directional bases nonzero.
- 372 actual PT-to-MN connections were active after calibration.
- Actual bound-dispatcher motor-current error was 5.995204332975845e-15 nA.
- Both initial knees were −10 degrees; no positive extension was introduced.
- Actual compiled mechanical parameter hashes were saved for comparison.
- **Zero body steps and zero training trials.** This establishes binding and
  electrical consumption, not standing, walking, or stability.

The source condition enables `SwingDiagnostics`, which installs an instance
forward closure. Copying a wrapper dictionary does not rebind that closure;
an initial cloned ON/OFF diagnostic therefore measured the wrong object's
state. The corrected check invokes the real dispatcher with saved/restored
neural state and conductances. The initial diagnostic failure is retained;
no production forward method or physical parameter was changed to resolve it.

## Entry points and artifacts

The machine-ready profiles are:

```text
D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\v2\function_equivalent\PROFILE.json
D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\v2\canonical_right_role_budget\PROFILE.json
```

Each directory contains the exact gain checkpoint, numerical compiled graph,
fresh-process check, integrated review, actual bound graph, and
`STATIC_BIND_PREFLIGHT.json`. Use these explicit v2 subdirectories. The v1
profiles/results remain in the original named directories; all 35 referenced
source files and both evidence sets were frozen byte-for-byte before the v2
change, with a checked 56-file archive and manifest in `integration/v1_frozen`.
Older profiles intentionally fail current source hashes; the archive preserves
their exact source revision. The earlier top-level `integration/PROFILE.json`
is an initial superseded checkpoint.

For the existing evaluator, add only:

```json
{"integrated_fsa_profile": "D:\\Users\\Ben Bolen\\NeuroCampaign_20261102_9h\\integration\\v2\\function_equivalent\\PROFILE.json"}
```

The default-off hook runs **after** the evaluator's existing gain transforms and
before binding. Per-step verification follows actual command-dependent weight
updates; it does not reinitialize fitted gains. Conflicting caller gain scales
are rejected when an explicit checkpoint supplies those same gains, rather than
silently ignoring a proposed training knob.

Programmatic entry points are `apply_to_network(net, profile, protocol_stage=...)`,
`build(acts, dt=.002, options=..., profile=..., protocol_stage=...)`, and
`verify(net, full=..., require_bound_PT=...)`. Explicit learned parameter
checkpoints can be written with `snapshot_parameters`; source changes require
an explicitly regenerated, reviewed profile rather than bypassing hashes.

For a bounded two-parameter Hunt-B study, use the already supported timing API:

```python
p = integrated_hunt_design.load_profile(primary_profile)
p["timing_profile"] = hunt_timing_profile.default_profile(
    mode="custom_budgets", relay_timing="fast_role_budget")
p["timing_profile"]["per_RG_budget_nA"] = {"E": E_budget, "F": F_budget}
```

This derives a fresh ledger from the two shared budgets, preserving source
within-pool relative contributions and bilateral tying. It must not be reported
using the old 0.5/0.5 labels. The initializer, PF/MN parameters, and local II
functions can remain fixed. Training launch and its bounded budget belong to
the authorized campaign supervisor. Current saved-trace evidence below favors
local motor/reflex calibration before changing the already regular RG clocks.

## Validator-only v2 change

The Ib verifier previously subtracted two full 1,033 by 1,033 matrices before
selecting cached edges. V2 selects those same entries first, then subtracts.
All guards remain active. The change removes an 8,536,712-byte temporary per
call; the corresponding selected range vector is 79,232 bytes.

`VALIDATOR_PARITY.json` verifies that this is the only Ib code change, that all
32 serialized numerical graph arrays are identical for both timing alternatives,
and that 128 commanded voltage/NaP steps are exact under old and new verifiers.
Cached-range, cached-conductance, and protected-window corruption are rejected
by both versions. The shared voltage hash is
`8903494995bba62b59a22665f76a94748658bfd0bc05407d18561be0bcdaedd4`.

On the same graph, 100 complete light per-tick checks took 0.587116 seconds
before and 0.270576 seconds after the change, a 2.17-fold validator speedup.
The Ib check alone improved 2.54-fold. These are diagnostic timings, not a
whole-model speedup. V2 also explicitly hashes the actual Hunt targets and
sensor metadata. No neural arithmetic or physical parameter changed.

## Observed short verification and next calibration

The supervisor's recorded source S1/S2 cases pass. The composed S1 result is
identical to source S1. In the recorded S2 cases, both new timing alternatives
remain finite and avoid crossing and hyperextension, but fail the physical
knee-phase criterion: function-equivalent phase error 0.208953, balanced phase
error 0.268619. V2's body repeat preserves the v1 assessment. These are suspended
air verifications, not evidence of ground walking.

Read-only saved-trace analysis over the actual assessment interval, 5–14 seconds,
locates the primary irregularity downstream of the RG/PF clocks. Primary right
RG-E period CV is 0.000560, compared with source 0.000485; right hip-E PF CV is
also 0.000560 versus 0.000485. Right knee-F PF timing remains similarly regular,
while physical knee peak intervals vary substantially. Only the two existing
heel-to-contralateral-RG-F connections provide explicit cross-central routes;
the intrinsic commissural gains are zero and contact is zero in these air
records. That can limit robustness, but its absence is not a demonstrated cause:
the source without intrinsic coupling passes and the new RG clocks stay regular.

`LOCAL_CURRENT_DIAGNOSIS.json` reconstructs fixed local currents using the actual
recorded pre-update full-neuron voltages. Independent recorded Ia/Ib aggregate
currents match within 3.47e-18 nA, and IaIN currents match exactly. In the new
right psoas during RG-E dominance, restored homonymous II supplies 0.88454 nA
versus zero in source, while activation rises from 0.05567 to 0.19074. Right
bifemlh and soleus receive 0.48668 and 0.44680 nA mean II current in the same
neural phase. Psoas reciprocal inhibition during RG-F dominance is −0.48399 nA
versus zero in source. By contrast, the representative vastus-intermedius Ib
group contribution changes from 0.00460 nA to zero. These are measured-trajectory
component replays, not new causal ablations or proof that every II route should
be weakened. RG dominance is not measured physical stance/swing in an air case.

The requested Hunt-D API has therefore been checked without body integration:

```python
p = integrated_hunt_design.load_profile(primary_profile)
p["parameter_checkpoint"] = None
p["local_profile"]["contracts"] = {
    "Ia_homonymous": ia_delta,
    "II_homonymous": ii_delta,
    "Ib_autogenic": -ib_absolute_delta,
    "Ia_reciprocal_pool": -reciprocal_absolute_delta,
}
```

Clearing the checkpoint is essential: an old saved reflex array would otherwise
overwrite a newly derived target. `LOCAL_TARGET_PREFLIGHT.json` verifies the
baseline ±0.5 mV contracts with checkpoint cleared reproduce the exact base
arrays. It also verifies II = 0.02 alone and the shared weak-all targets
(+0.1, +0.02, −0.1, −0.1). Each consumes 368 actual pool targets, persists through
all four command values, and preserves PF, FF, biarticular, and timing buffers
exactly. `LOCAL_TARGET_GAINS.npz` contains the exact named gains. Suggested small
engineering ranges are Ia 0.05–0.5 mV, II 0.005–0.15 mV, and inhibitory magnitudes
0.05–0.5 mV; these are not identified physiological bounds. Zero is excluded by
the signed functional-contract definition. No training was launched by this
reviewer.

## Dissertation-family comparison scope

This implemented Hunt candidate is distinct from the dissertation S3K, W2L, and
Syn6 saved recipes. Their selected source choices remain:

- Full S3K: four phase PF cells per leg, explicit local interneurons, and active
  flexor Ia/II and extensor Ib/II central routes. The pruned version retains
  central Ia/II and local Ib circuitry while removing the declared contact,
  direct central Ib, weak-RG, and interleg components.
- W2L: separate hip PF pair and a shared knee/ankle PF pair, explicit reflex
  interneurons, active normalized central Ia/II/Ib routing, and the actual
  ipsilateral/contralateral heel relay pattern. It does not supply an independent
  ankle PF pair or the newer four-cell RF/hamstring circuit by name alone.
- Syn6: the frozen K6 factors and its own explicit reflex dressing. Its saved
  source recipe lacks direct F-Ia/II central projections and uses the questioned
  positive-both-ipsilateral heel pattern. Optional rules or a K8 inverse fit do
  not constitute an implemented forward repair or walking result.

The paper-to-code comparison and unresolved functional rules remain recorded
separately in `GAIT_DISSERTATION_RECIPE_REVIEW.md` and the parent's literature
review. The Hunt composition does not establish a global all-family design,
physiological completeness, free-3D control, or curriculum pass.
In particular, the Hunt motor mask does not automatically restore every
allocated RC or antagonistic-II motor projection. The local composition covers
the four explicitly declared functional classes; broader biological pathway
choices need their own review and cannot be inferred from the neuron census.
