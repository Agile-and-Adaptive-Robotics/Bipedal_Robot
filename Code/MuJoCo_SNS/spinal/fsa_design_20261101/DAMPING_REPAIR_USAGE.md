# Optional damping-realization candidate

Machine: **easteregg2**. The implementation is
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\damping_repair.py`.
No production source, historical result, default builder, or controller is changed.
No body simulation, Hunt fitting, or curriculum run was performed.

## Why a distinct physics variant is needed

`runner.patch_xml` intends to put damping on each added support spring, but its
regular expression only replaces a damping attribute already written on a joint
tag. Thirteen internal joints inherit damping from the default class, so the
legacy patch leaves their compiled damping at 0.05. Pelvis tags already have
explicit damping, and their patch already works.

The optional repair first writes each named joint's **actual compiled raw
baseline damping** into its existing raw XML tag when absent, then calls the
original `runner.patch_xml`. It preserves the rest of the XML byte structure,
including the exact compiler and joint-axis patterns on which that function
depends. It does not serialize the XML tree or duplicate the patch's numerical
spring/damper tables. It follows original caller overrides and rig scaling.

In the tested default ground patch, five lumbar/hip-rotation **rig** coordinates
change from 0.05 to 25 N m s/rad. Eight ankle/subtalar/MTP/hip-adduction restraint
surrogates change from 0.05 to 3 N m s/rad. These are the existing engineering
specifications being realized, not newly measured anatomical damping values.

## Interfaces and safe application order

`materialize_damping(raw_xml, raw_model, policy=None)` returns XML and a binding
ledger. `patch_xml(raw_xml, raw_model, raw_data, policy=None, **original_kwargs)`
is a pure wrapper around the original patch. `enabled(policy=None)` is a scoped
callback for callers that already use `runner.apply_harness`.

With `None`, materialization is a no-op, the pure wrapper returns XML identical
to the legacy patch, and the context installs no callback. A supplied policy
must equal `default_policy()`; unknown or unused keys are rejected.

For future model setup in an isolated process, with existing import paths set:

```python
from fsa_design_20261101.damping_repair import default_policy, enabled

# raw_model/raw_data must be the original converted model and starting pose.
with enabled(default_policy()) as damping_ledger:
    model = R.apply_harness(raw_model, raw_data, **existing_harness_options)
# The original R.patch_xml callback has now been restored.
```

The callback is process-local and not intended for concurrent threads or nested
repair contexts. It restores the original callback on normal exit and exceptions.
The existing loader creates and deletes its scoped temporary XML in the model's
own directory to preserve relative asset paths. Accepted controllers must not be
described as validated under this candidate until separately evaluated under its
versioned configuration after the campaign gate permits it.

Run the bounded **static audit only** in PowerShell:

```powershell
$env:CONDA_PREFIX='D:\Anaconda\envs\myo'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\damping_repair.py' --output-dir 'D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\damping_repair_v1'
```

Both models receive the same shared-force-law v3 parameter patch. The comparison
checks **396 non-damping compiled fields exactly**, including mass, inertia,
axes, muscle/force arrays, equalities, compiler-derived arrays, options, and
statistics. Only `dof_damping` differs. M and K are exactly equal at all four
reference postures; D is independently recomputed from native force derivatives.
The context-restoration and unknown-policy rejection tests also pass.

## Static mechanics interpretation

The files use four old walking postures, hold all six pelvis coordinates fixed,
and reduce to 17 internal coordinates through the ideal equality projection E:

    M = E.T @ M_native @ E
    r = E.T @ (qfrc_bias - qfrc_passive - qfrc_actuator), at activation=0, v=0
    K = dr/dq
    D = d[E.T @ (qfrc_bias-qfrc_passive-qfrc_actuator)]/dv, qvel_native=E@v

The report separates the frozen-inertia impedance matrix

    A_frozen = [[0, I], [-M^-1 K, -M^-1 D]]

from the actual reduced zero-velocity acceleration tangent. Because these
postures are not equilibria, a0=-M^-1 r is generally nonzero, and its position
derivative includes the additional inertia-derivative term:

    da/dq_j = -M^-1 (K[:,j] + (dM/dq_j) @ a0).

That derivative is independently evaluated by central differences of the entire
reduced acceleration calculation. Spectral damping ratios are reported for
complex modes as -Re(lambda)/abs(lambda); real poles have no asserted oscillatory
damping ratio. All spectra and actual 34x34 A matrices are saved.

At phases 0/25/50/75%, the actual local maximum-real poles are respectively
1.649/0.930/1.351/0.822 per second for legacy damping, and
1.465/0.715/1.084/0.462 for the candidate. They remain positive. Frozen-inertia
poles are negative, illustrating why the simplified calculation cannot establish
stability at a non-equilibrium posture. Neither calculation includes contacts,
soft follower-constraint compliance, floating-pelvis motion, neural feedback, or
activation/fiber dynamics. This is not standing or walking validation.

## Log all apparatus assistance, not only vertical support

`rig_assistance(model, data, raw_baseline_model)` reads q/v and current kinematic
axes without changing state. The caller supplies data whose kinematics correspond
to that timestamp. It enumerates **all 11 declared rig coordinates**: pelvis
tx/ty/tz/tilt/list/rotation, lumbar extension/bending/rotation, and bilateral hip
rotation. The latter five remain rig assistance even at zero vertical BW support.

For each scalar coordinate:

    tau_spring = -k*(q-qref)
    tau_damper = -c*v
    power_spring = tau_spring*v
    power_damper = -c*v*v
    spring_energy = 0.5*k*(q-qref)^2

Translations use N and m; rotations use N m and radians; both powers use W.
The helper records the full implemented forces and their increment over the
explicit raw baseline's spring/damping coefficients, so inherited baseline
damping is not silently relabeled as apparatus. Root translation force is mapped
through actual world axes and projected onto gravity-derived up; normalization
uses moving-body mass, excluding fixed ground. Signed and positive/negative power
are retained rather than concealing assistance through cancellation.

The static formula fixture with only internal rig offsets/velocities has exactly
zero vertical support but nonzero internal torques and -0.023475 W incremental
rig power. This is an artificial logging test, not a measured walker trajectory.

For a later authorized run, preallocate one column per named rig coordinate for
q, v, spring/damper force, incremental force, and signed power at actual physics
timestamps. Preserve native units and actual rig scale/reference values. Integrate
positive and negative work separately; record spring energy and damper work.
Use timestamp-consistent solver/kinematic sampling, and distinguish instantaneous
constitutive power from numerical step work. Do not call an extra full dynamics
update merely to log, since that can alter solver warm-start behavior. Zero
vertical support is not apparatus-free balance, and dissipated apparatus energy
is not metabolic cost.

Evidence:
`D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\damping_repair_v1\`
contains `DAMPING_REPAIR_STATIC.json`, plus separate legacy/repaired complete
mechanics records. Each stores the frozen input hashes. Global campaign release
remains false.
