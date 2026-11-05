# Physical kinematic response objective

**Scope correction:** the v1 grounded-reference gate described below is retained
as history. The current default is `stage_specific_kinematics_v2`. During physical
air stages, database gain/phase/amplitude/cadence differences are **diagnostic**,
not grounds for rejection or optimization targets. Air ranking uses measured
bilateral phase error, period variability, and cycle variability, subject to the
existing motion and safety gates. Ground-stage reference matching remains a
separate target. See the condition-specific note at the end.

Session 1; easteregg2; 2026-11-05. No manuscript or physics changes.

## Implementation

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\controller_profiles_20261104\kinematic_response.py`
scores measured hip, knee, and ankle coordinates on **both** sides against IK
reference cycles. `body_probe.py` now writes `body.kinematics.json` after each
completed physical evaluation. `fit_body_feedback.py` selects the new objective
by default and freezes its source/reference hashes in the study manifest.
`body_fit_adapter.py` consumes that physical score when the manifest selects
`objective_profile=database_kinematics_v1`. Old manifests retain their old
objective; no historical results were rewritten.

This corrects the latest follow-up's reliance on neural replay and broad physical
motion gates. The repository's older `kine_ref.py` already contained bilateral
mean, amplitude, phase, and periodicity comparisons; the newer offline/body-fit
path had not carried the complete comparison through to its optimization
objective. Existing calibration block C/D scores remain separate diagnostics.

## Definition

For each of the six joint coordinates, the metric reports:

- Actual stride period and error against that side's reference period.
- Cycle excursion gain, `ROM_sim / ROM_ref`.
- Mean coordinate bias in degrees, retaining the OpenSim zero/sign convention.
- Uncentered waveform RMSE, normalized by reference excursion with a 5-degree
  denominator floor.
- Complex first-, second-, and third-harmonic ratios, including gain in linear
  units and dB and signed phase lag in degrees.
- Cycle-to-cycle RMS variability, reported separately.

For the Fourier convention `q(phi)=mean+Re(C_k exp(i*2*pi*k*phi))`, the
descriptive harmonic ratio is `H_k=C_sim,k/C_ref,k`; gain is `abs(H_k)` and
lag is `-arg(H_k)`. Positive lag means simulation lags the reference. A reference
harmonic smaller than 1 degree is omitted from ratios. Near-zero simulated
harmonics do not receive a favorable phase score.

One shared right-knee maximum-flexion event anchors all six coordinates in air.
There is **no independent per-joint phase alignment**. Stored database left/right
cycles start at their respective loading onsets, so the reference left cycle is
first shifted by its recorded left/right loading phase. These air-stage phase
measurements are relative to the knee anchor, not absolute heel-strike phase.
Ground-stage event anchoring should use measured loading onsets instead.

Cycle normalization permits shape comparison across cadences, but actual cadence
is retained and penalized independently. No dynamic time warping is performed.

The scalar physical ranking loss is:

`J = 2*mean(abs(log(T_sim/T_ref)))`

`  + mean(normalized waveform RMSE)`

`  + 0.5*mean(abs(log(ROM_sim/ROM_ref)))`

`  + 0.5*mean(abs(mean_bias)/reference_ROM)`

`  + 0.5*mean(abs(fundamental_lag_degrees)/180)`

`  + 0.5*mean(worst_cycle_absolute_error/reference_ROM)`

The phase term now uses each cycle's absolute wrapped phase error, rather than
only the phase of the averaged waveform. The maximum-error term covers every
retained cycle and was added with the statistics revision below. Historical
v1/v2 scalar values remain in their original reports.

The reference-ROM denominator floor is also used for normalized mean bias.
These weights are explicit engineering choices, not literature-derived constants.
They can be examined through the separately saved terms. Normalizing keeps large
hip excursions from dominating merely because their numerical scale is larger.

Existing finite/full-duration, protected-knee, frontal-clearance-origin, and
principal-motion gates remain. Either a failed gate or more than 15% period error
adds 1000. Invalid/frozen/fast traces without at least two complete usable cycles
receive 10000. No waveform-matching cutoff has been validated yet: `accepted`
means those movement/safety/cadence gates passed, **not** that the waveform is
fully validated. Full mesh clearance, measured load, and ground alternation
remain separate requirements for ground-stage progression. This scorer is
currently connected to air-stage probes, not a replacement for those requirements.

## Saved-run findings

The primary reference is the existing `subject01_walk1` IK/GRF pair, with periods
1.2334/1.2133 s. It is fixed before scoring; neither controller selects the
easiest reference.

| Measurement | s3k repaired-readout stage 2 | syn6 ankle candidate stage 2 |
| --- | ---: | ---: |
| Measured R/L period (s) | 1.386 / 1.446 | 2.004 / 2.008 |
| Hip excursion gain R/L | 2.879 / 2.943 | 0.587 / 0.599 |
| Knee excursion gain R/L | 0.354 / 0.249 | 0.290 / 0.322 |
| Ankle mean bias R/L (deg) | -33.14 / -32.49 | -25.51 / -26.21 |
| Shape/phase/cadence loss before failure penalty | 2.30755 | 2.35239 |
| Database-response gate | Fail | Fail |

The negative ankle bias indicates more plantarflexion relative to this reference.
Cycle-mean gains can be smaller than peak-to-peak ranges over the entire record
because cycle averaging exposes phase drift and variability.

Both saved trajectories were also compared with all **43** local database NPZ
references. Every comparison is retained in a compact `DATABASE_COMPARISON.json`.
This is fixed-command coverage, **not** proof of speed adaptation or a held-out
validation campaign. Walking and running conditions should be selected explicitly
for later condition-specific fitting/testing. Do not report the minimum loss
over unrelated references as successful tracking.

## Verification

Analytical synthetic checks recover perfect-match loss 0, doubled amplitude gain
2, and delayed phase 45 degrees. They also reject phase reversal as a perfect
match, penalize coordinate offsets, and penalize half-cadence motion even when
its normalized waveform matches exactly.

A fresh physical run at
`D:\Users\Ben Bolen\ControllerTraining_20261104\kinematic_metric_integration_20261105`
automatically generated `body.kinematics.json`. Its joint positions/velocities,
recorded neural voltages, foot forces, and harness wrench are array-exact matches
to the earlier syn6 candidate. Only postprocessing/scoring changed. No new
parameter optimization was launched in this scoring pass.

## Gain/phase plots versus Bode plots

The provided harmonic plots compare **simulation output to reference output**.
They are not neural-input-to-joint-output transfer functions and are not labeled
as Bode plots. Their horizontal axis is stride harmonic order, not an imposed
perturbation frequency.

A true input-output frequency-response experiment needs a declared input
(for example, small command-current modulation), its measured time series, and
the resulting joint-response perturbation. The unperturbed baseline, amplitude
sensitivity, transients, and excitation quality must also be checked. A periodic
walking controller additionally needs phase-aware interpretation. The separate
standing inverted-pendulum model is a useful simpler setting for an initial
small-signal input/output frequency sweep. No such sweep is claimed here.

## Files for the defense builder

Portable metrics, cycle overlays, and gain/phase PNG/PDFs:

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104\kinematic_response_20261105`

Each `s3k`/`syn6` subfolder contains `KINEMATIC_RESPONSE.json`,
`DATABASE_COMPARISON.json`, `IK_OVERLAY.png/.pdf`, and
`HARMONIC_GAIN_PHASE.png/.pdf`.

Historical source trajectories:

`D:\Users\Ben Bolen\ControllerTraining_20261104\s3k_repaired_body_20261104\expression_0.2\stage2`

`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_ankle_selected_repeat_20261104\stage2`

The source module's command line accepts `--run`, `--output`, an optional explicit
`--reference` NPZ, and an optional `--database` directory. With no `--reference`,
it uses the fixed original calibration IK/GRF pair. `--self-check` runs the
analytical metric checks without physics.

## Condition-specific interpretation: correction after Ben's clarification

Hunt, Szczecinski, and Quinn (2017), section 2.4.1, deliberately tune the
endogenous CPG period to twice the intended stepping period. Section 2.4.2
then tunes sensory entrainment using afferents calculated from intended walking
kinematics. Source:
https://pmc.ncbi.nlm.nih.gov/articles/PMC5378996/#s4
This is an endogenous-oscillator design setting, not a general empirical law
that all afferented suspended legs move with twice the grounded period.

Our 2.466 s intrinsic calibration target was explicitly set from 2*1.233 s.
Reproducing it is not independent evidence for that ratio. Deafferented
intrinsic activity, afferented unloaded stepping, and loaded ground stepping
must be reported separately. Their relation remains an experimental question.

The current air objective is
`5*bilateral_phase_error + period_CV + normalized_cycle_variability`, with the
existing safety/motion rejection penalties. Grounded IK differences are retained
but omitted from that objective. The proposed ground objective above retains
cadence, amplitude, bias, and phase matching. Neither establishes a universal
normal air-step trajectory. In particular, the earlier 26-degree ankle bias
quantifies an air-versus-ground difference; by itself it does not establish
pathological air stepping.

The same fixed saved syn6 candidate passes the air gate under v2 while still
showing its database differences. s3k remains an air-motion failure because its
hip excursions violate the pre-existing motion criterion. The original v1
comparison files remain untouched; condition-correct reports are under
`kinematic_response_20261105\air_scope_v2` in the portable result folder.
No training was launched under the short-lived v1 objective.

A useful follow-up is a paired, fixed-parameter series: intrinsic rhythm,
afferented air stepping, then supported ground stepping, with multiple commands
and repeated initial conditions. Measure period ratios, relative joint phases,
amplitude changes, and whether air features predict successful ground stepping.
Fit any air-to-ground relation on one set of conditions and test it on held-out
conditions; do not assume a fixed factor for all three preparations.

## Full goodness-of-fit statistics

The current `amplitude_frequency_cycle_phase_maxerror_v1` reporting revision
generates `RESULTS.md` automatically beside every new kinematic report. Its
primary tables contain total peak-to-peak amplitude, minima/maxima, per-cycle
amplitude mean/sample SD and percentage error; period and frequency mean/sample
SD and percentage error; per-cycle fundamental phase circular mean/SD and worst
absolute phase error; and worst absolute waveform error over **all** retained
cycles, with per-cycle maximum-error mean/SD. RMSE and mean absolute error are
secondary diagnostics. Undefined low-amplitude phase is identified rather than
assigned zero error.

The ground-reference ranking loss also includes normalized worst-cycle maximum
error. Air ranking remains independent of ground-reference mismatch. Exact
sample counts accompany the tables. These are descriptive statistics from one
run, not independent-subject significance tests. A spike check verifies that
one 100-degree error in one cycle remains a 100-degree maximum error even when
another cycle is perfect; a circular test checks the +179/-179-degree boundary.

Updated tables:

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104\kinematic_response_20261105\air_scope_v2\s3k\RESULTS.md`

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104\kinematic_response_20261105\air_scope_v2\syn6\RESULTS.md`
