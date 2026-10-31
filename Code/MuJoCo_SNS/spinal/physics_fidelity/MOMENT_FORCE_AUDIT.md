# Matched-coordinate route and bilateral force-law audit

easteregg2, 30 October 2026. No production model or parent controller changes.
Active raw evidence: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\runs\audit_20261030\`.

## Coordinate adapter and method

`bilateral_moment_arm_audit.py` uses the repaired forward MuJoCo model and
`CoupledArms` to include polynomial equality-follower motion in moment arms.
It compares twelve leg muscle families over all seven ipsilateral coordinates
at neutral and at identical **bilateral anatomical** hip flexion 30°, knee
angle −60°, and ankle angle 10°; other selected coordinates are zero.
The source OpenSim model's `computeMomentArm` supplies the independent reference.

Both poses and moment arms must use the established FK sign adapter:

- Native MuJoCo hip_flexion and hip_adduction have sign **−1** relative to
  OpenSim on both sides. Knee/ankle and the other audited coordinates use +1.
- Anatomical moment arm is `native_moment_arm * coordinate_sign`.
- The old converted model's world frame is x anterior, y left, z superior.
  OpenSim/ISB is x anterior, y superior, z right. This audit explicitly adapts
  the older baseline; it does not redefine the anatomical frame.

Ignoring these signs produced sixteen apparent hip-sign disagreements; all
disappear with the adapter. They are not sixteen reversed muscle routes.
MuJoCo analytical effective moment arms agree with tendon-length finite
differences to at most 1.29e−11 m/rad in the two tested poses.

## Main route findings

There are **no significant sign disagreements** after coordinate adaptation.
The largest L/R route discrepancy among all twelve families and seven axes
is **0.0958 mm**, semimembranosus about hip rotation at the swing pose.
BFSH/vasti/soleus/tibialis-anterior main-action routes are essentially bilateral.
This is a two-pose audit, not full-ROM certification.

| Route | Neutral MuJoCo / OpenSim (mm) | Swing MuJoCo / OpenSim (mm) |
| --- | --- | --- |
| BFSH knee | −21.134 / −20.361 | −33.660 / −33.113 |
| Vastus medialis knee | 51.848 / 46.301 | 38.792 / 35.538 |
| Vastus intermedius knee | 51.334 / 48.719 | 40.312 / 37.338 |
| Vastus lateralis knee | 51.088 / 42.782 | 38.748 / 35.266 |
| Soleus ankle | See full JSON | −42.0628 / −42.0626 |
| Tibialis anterior ankle | See full JSON | 43.2213 / 43.2208 |

There are real cross-engine differences, for example vastus lateralis differs
by 8.31 mm (19.4%) at neutral. These are predominantly shared bilateral
differences, not a unilateral route reversal that directly explains the left
extra pulse. Small nonzero near-zero moment arms should not be called failures.

## Route length versus fitted normalization

| Quantity | Right | Left |
| --- | --- | --- |
| BFSH Lmt, neutral | 0.24598156 m | 0.24598156 m |
| BFSH Lmt, swing | 0.21431849 m | 0.21431849 m |
| BFSH actuator lengthrange | [0.191264, 0.249268] m | Same |
| BFSH inferred fitted fiber optimum | 0.08538614 m | 0.09818306 m |
| BFSH inferred fitted tendon offset | 0.13204308 m | 0.12454311 m |
| Soleus Lmt, neutral | 0.29139796 m | 0.29139796 m |
| Soleus Lmt, swing | 0.29904149 m | 0.29904149 m |
| Soleus actuator lengthrange | [0.229942, 0.326094] m | Same |
| Soleus inferred fitted fiber optimum | 0.09315610 m | 0.09043347 m |
| Soleus inferred fitted tendon offset | 0.22901044 m | 0.22606735 m |

The inferred values follow MuJoCo's operating-range transformation and are
surrogate normalization parameters, not measured biological fiber/tendon lengths.
BFSH's source OpenSim Lmt is 0.24625076 m neutral and 0.21421546 m swing.
The small route-length error cannot explain a large L/R difference: the
converted L/R lengths and lengthranges themselves match here.

## Why maximum force alone concealed the problem

Independent full-activation-only curve fits can trade active and passive force.
For BFSH at neutral, direct MuJoCo zero-velocity measurements give:

| Component | Right | Left |
| --- | --- | --- |
| Active force increment at activation 1 | 520.62 N | 605.54 N |
| Passive tension | 215.57 N | 127.07 N |
| Total tension at activation 1 | 736.19 N | 732.60 N |
| Total tension at activation 0.05 | 241.60 N | 157.34 N |

Thus the full-activation totals look similar while low-drive forces differ by
84.26 N. At swing, BFSH active capacities are 668.81/672.06 N and passive
tension is zero. The asymmetry is pose- and activation-dependent; it is not
captured by one Fmax ratio. This is a plausible contributor to unequal motor
responses with tied neural weights, but it is not a causal gait experiment.

The old cached OpenSim forces were checked against fresh equilibria for BFSH
and soleus: maximum discrepancies are about 1.14e−13/4.55e−13 N. A provisional
concern about activation transients in the MuJoCo conversion evaluator was also
tested: its preparatory 5 ms step takes activation 0→1 with this model's actual
activation parameters, so subsequent sampled steps use activation 1. **The
transient-contamination hypothesis is not supported for this configuration.**
Do not present it as the cause of the original conversion differences.

## Joint bilateral surrogate diagnostic

`shared_force_surrogate_audit.py` fits one four-parameter set per muscle family
to both legs using fresh OpenSim static equilibrium at activations 0.05, 0.3,
and 1. It shares the right lengthrange, fits normalized operating-range ends,
passive maximum, and force scale, and retains the native MuJoCo force law.
Eight functional training poses vary hip, knee, and ankle; two intermediate
poses are held out. The source targets include compliant-tendon equilibrium.
Fitting all six families takes about one second after source-data collection.

| Muscle | Baseline held-out RMSE | Copy-right held-out RMSE | Joint-fit held-out RMSE |
| --- | --- | --- | --- |
| BFSH | 10.37 N | 11.82 N | 3.07 N |
| Vastus medialis | 37.65 N | 41.87 N | 7.90 N |
| Vastus intermedius | 35.00 N | 40.29 N | 9.10 N |
| Vastus lateralis | 50.96 N | 50.54 N | 12.24 N |
| Soleus | 81.79 N | 78.95 N | 29.27 N |
| Tibialis anterior | 16.66 N | 19.84 N | 3.53 N |

The joint fit reduces held-out source-force error for all six families and
makes matched-pose bilateral force differences smaller than 0.0083 N. Five
families agree to numerical precision; vastus intermedius retains a small
geometry-dependent residual (0.00822 N versus 180.23 N baseline). Copying the
right parameters alone makes held-out error worse for four families. Prefer
the jointly fitted variant for an explicit diagnostic, rather than describing
arbitrary right-to-left copying as improved OpenSim fidelity.

The BFSH upper operating endpoint hits its allowed lower bound, and soleus's
lower endpoint hits its lower bound; extrapolation must be checked. These
results do not validate the full RoM, dynamic force–velocity behavior, extreme
joint angles, or a whole-body controller. Only six families are included.
No claim of compliant-tendon dynamics follows from static surrogate agreement.

## Opt-in parameter file

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\shared_force_law_functional_v1.json`

Schema `aarl_shared_surrogate_force_law_v1`: top-level `muscles`, keyed by
family; each contains `actuator_names`, `gainprm`, `biasprm`, `dynprm`,
`lengthrange`, and validation metrics. Assign the four parameter arrays to both
listed actuators. Activation-dynamics parameters already match bilaterally;
the copied dynamics do not establish any new dynamic validation.
Keep actuator types, control bounds, routes, masses, damping, timestep, and
integrator unchanged. Apply before caching Fmax/normalization or solving
standing activations, and record the JSON hash. No default loader is changed
by these audit scripts. `export_shared_force_law.py` refuses an export if
held-out error worsens or matched-pose bilateral differences exceed 0.01 N.

Run the matched-pose audit, shared-surrogate audit, and exporter with
`D:\Anaconda\envs\myo\python.exe` from `D:\GitHub\Bipedal_Robot`.

Raw JSON files: `BILATERAL_MOMENT_ARMS.json`, `SHARED_FORCE_SURROGATE.json`,
`CONVERSION_CACHE_STATIC_CHECK.json`, and `CONVERSION_ACTIVATION_CHECK.json`.

## Hip/gastrocnemius extension: functional v2

`shared_force_surrogate_audit.py --hip-v2` evaluates psoas, iliacus, three
gluteus-maximus divisions, and medial/lateral gastrocnemius. The original grid
had only two training hip angles; that is insufficient coverage for a confident
four-parameter hip-law fit. The extension uses ten functional training poses
spanning hip −15° to 60°, knee 0° to −90°, and ankle −20° to 20°, plus four
held-out intermediate combinations. Activations remain 0.05/0.3/1, and fresh
source OpenSim equilibrium supplies the target tensions. There are no failed
source equilibria. Seven fits take approximately 0.79 s after data collection.

| Family | Baseline held-out RMSE | Joint held-out RMSE | Largest L/R residual | v2 decision |
| --- | --- | --- | --- | --- |
| Psoas | 13.14 N | 3.92 N | Numerical precision | Include |
| Iliacus | 31.72 N | 2.79 N | 2.9245 N | Withhold under strict 0.01 N symmetry gate |
| Gluteus maximus 1 | 32.60 N | 4.96 N | Numerical precision | Include |
| Gluteus maximus 2 | 23.36 N | 5.22 N | Numerical precision | Include |
| Gluteus maximus 3 | 9.52 N | 1.97 N | Numerical precision | Include |
| Medial gastrocnemius | 88.02 N | 48.07 N | Numerical precision | Include |
| Lateral gastrocnemius | 14.31 N | 6.49 N | 0.1976 N | Withhold under strict 0.01 N symmetry gate |

All seven improve held-out static force error. The two withheld families do
not have failed fits: small residual L/R route lengths prevent near-exact force
parity with a common law. Iliacus differs by at most 0.1896 mm in Lmt and retains
2.9245 N mismatch, versus 159.09 N baseline (source F0=1073 N). Lateral
gastrocnemius differs by 0.02060 mm and retains 0.1976 N mismatch, versus
13.07 N baseline (F0=683 N). These are small residual mechanical differences,
not evidence of reversed routes or complete failure. The 0.01 N inclusion
gate is a deliberately strict diagnostic convention, not a biological standard.

NEW file: `shared_force_law_functional_v2.json`. It preserves every v1 muscle
parameter exactly and appends the five included families: eleven families,
twenty-two actuators total. It keeps the same format/schema and adds explicit
extension training/held-out poses and rejection reasons. `export_shared_force_law_v2.py`
verifies that v1 remains byte-unchanged. Raw fits are in
`runs/audit_20261030/SHARED_FORCE_SURROGATE_HIP_V2.json` under this directory.

v1 SHA256: `67f0cad8e69d81f1042fac30e1ebe023e93cd6d0a924897ef8bf1a2e19d9bb5a`.
v2 SHA256: `a288dd88ceb3888e69b13f40bcbd140189f6d7aad1e774d45499413bf65271a7`.

Neither v2's fitting nor its export runs or certifies locomotion, dynamic force
responses, full-ROM extrapolation, or compliant-tendon dynamics. Ground/early-stage
evaluations remain separate parent-session work with the explicit variant.

## Biarticular and low-drive extension: functional v3

NEW file: `shared_force_law_functional_v3.json`; v1/v2 remain byte-unchanged.
v3 preserves every v2 law and adds iliacus, lateral gastrocnemius,
semimembranosus, semitendinosus, biceps femoris long head, rectus femoris, and
TFL: eighteen families, thirty-six actuators. The new hip/knee profile samples
twelve training poses and five held-out combinations spanning hip −15° to 60°,
knee −90° to 0°, and hip adduction −10° to 10°. Adduction variation is included
to give TFL a meaningful additional functional dimension; hip rotation remains
zero and its effects are not certified. Activations are 0.01/0.05/0.15/0.3/1.

| Family | All-drive held-out RMSE, baseline → shared | Low-drive held-out RMSE (0.01/0.05), baseline → shared |
| --- | --- | --- |
| Semimembranosus | 44.85 → 25.80 N | 6.30 → 5.53 N |
| Semitendinosus | 9.03 → 1.71 N | 7.02 → 0.84 N |
| Biceps femoris long head | 23.03 → 11.13 N | 12.19 → 5.14 N |
| Rectus femoris | 259.69 → 22.78 N | 61.78 → 11.56 N |
| TFL | 8.07 → 7.31 N | 9.68 → 8.01 N |

The rectus-femoris force discrepancy is especially important: the forward
repair reroutes it through a patellar tracking point but does not refit its
old Step-3 force parameters afterward. The v3 fit evaluates the repaired route.
This is an alignment issue worth correcting, but its causal effect on stepping
still requires a controlled model comparison.

### Parameter symmetry and source force-difference checks

The strict v2 requirement of delivered-force equality within 0.01 N is relaxed
in v3. **Shared parameters do not imply equal forces when geometry or pose
differs.** At matched anatomical poses, v3 compares the L/R force difference
with the source OpenSim L/R difference. Inclusion requires improved held-out
source-force RMSE and no worsening of the maximum source force-difference
error. An explicitly declared diagnostic residual bound is 1% of source Fmax,
not a physiological standard.

This admits iliacus and lateral gastrocnemius using their already improved
fits. Their residual errors are 2.9245 N and 0.1976 N respectively, reflecting
small route-length differences. For the new biarticular profile, maximum
source force-difference errors drop from baseline to shared as follows:
semimembranosus 130.47→0.8629 N; semitendinosus 36.57→0.00731 N; biceps long
head 128.35 N→numerical precision; rectus femoris 571.57→0.00819 N; and TFL
18.74 N→numerical precision. No actuator route or body geometry is altered.

### Passive hamstring torque: important causal qualification

`v3_passive_component_check.py` separately compares actual source passive
components at activation 0.01 with the original, v2, and v3 rigid surrogate
biases at matched neutral, bent-standing, and swing poses.

At neutral, the original four-hamstring passive knee torque is −4.56/−2.69 N·m
right/left versus source −0.070 N·m each. At bent standing it is −4.02/−2.19
versus source −0.078 N·m. **v2 already removes essentially all this original
artifact through its BFSH law:** its summed passive torque is zero at both
poses. v3 adds only approximately −1.2e−5 N·m at bent standing, and all
variants have zero hamstring passive torque at the audited swing pose.
Therefore the hypothesis that *remaining* hamstring passive bias explains
the parent's current v2 ground failure is not supported by these three poses.
The new laws improve other pose/activation errors, but do not demonstrate an
additional several-N·m standing hamstring-bias removal. A broader dynamic
trajectory may visit other poses and remains to be checked.

Source compliance makes fiber state and passive force activation-dependent;
the rigid MuJoCo passive bias is activation-independent. Agreement of total
tension at low drive does not establish exact component-level equivalence.

### Evidence durability and provenance

Active evidence now lives under the local ignored `runs/audit_20261030/`,
because earlier repo `tmp/` outputs disappeared for an unknown reason.
`REGENERATED_EVIDENCE.json` records eight successfully reproduced checks and
confirms that all three curated profile hashes remained unchanged. Regenerated
raw files are new evidence, not claimed byte-identical historical recovery.

Raw v3 files: `SHARED_FORCE_SURROGATE_BIARTICULAR_V3.json`,
`BIARTICULAR_V3_LOW_DRIVE_DIAGNOSIS.json`, and `V3_PASSIVE_COMPONENT_CHECK.json`.
The shared v3 parameter file SHA256 is
`9dc7d7ba444a961e06fd6efd15210ea324a84ee70992b251994cb7d7ee7165cf`.
The existing loader schema remains unchanged. Neither fitting nor export
certifies a curriculum pass, full-ROM/dynamic fidelity, or walking causality.
