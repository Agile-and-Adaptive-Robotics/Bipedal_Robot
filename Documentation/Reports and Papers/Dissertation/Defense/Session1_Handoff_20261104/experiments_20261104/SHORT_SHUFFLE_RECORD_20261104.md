# Short, slow W2L walking and sensory-gate diagnostics

Machine: **easteregg2**, 4 November 2026. Newer Li–Deng-derived W2L108.
Ben requested short strides, slower cadence, longer double support, and a
treadmill-direction check. No manuscript, CAD, muscle/tendon law, or damping
parameter was changed. Optional controller/readout variants are explicit.

## Verified partial result

The best slow candidate completed 16 s, with 7/7 alternating unloaded-clear
swings and approximately 2.107/2.120 s limb periods. Upward rail/guide assistance
averages 10.105% BW over 2–16 s. Both feet fall below 0.05 BW load for 0.714%
of that window. This load statistic is not a geometric airborne measurement.
The entire saved numeric trajectory repeats byte-exactly in a new process.

The target short shuffle is **not achieved**: fore-aft foot excursion relative
to the root remains 0.351/0.417 m; whole-foot peak clearance is 75/107 mm.
Loaded double support is 3.571%, and overlap in bouts lasting at least 40 ms
occupies only 1.929% of the window. The new gate requires >=8% sustained overlap,
<=0.20 m excursion, low load-loss time, >=5 qualified swings per side, slow
cycles, no hyperextension, and limited upward/downward apparatus assistance.
No candidate passed all criteria; the criteria were not relaxed after the runs.

Selected shared parameters: belt speed 0.05290055285 m/s; neural time scale
2.101662273; extensor RG tonic current 0.0989108632 nA; hip PF→MN scale
1.037563074. Existing PF→MN scale 1.115222745 and toe→MN conductance
0.0488441541 µS are retained. These are settings, not measured biological values.

## Direction and mechanical checks

The original candidate's toes extend +0.0842 m in world X from the feet.
Its actual belt velocity is approximately −0.1001 m/s: correct rearward motion
relative to the model's forward +X. The slow candidate has a similarly negative
measured belt velocity at its lower command. New previews include a green
forward arrow and a blue arrow derived from saved belt displacement and the
actual world slide axis. Arrow length does not encode velocity magnitude.

Imported forward kinematics was compared with 120 native AnimatLab samples of
feet, toes, and hip-wing positions. With translation removed and unrecorded root
rotation fitted, the current joint signs give mean positional RMS 0.808 mm;
reversing signs gives 137 mm. This is a geometric check, not dynamics parity.
Right-hip muscle names are functionally reversed, but the native MN-output
adapters already compensate; tracing only names would incorrectly suggest a
wiring repair. Ankle MN-output adapters are crossed as well. Effective torque
directions support the existing knee-PF→ankle-MN phase correction.
Muscle Euler timestep/time-constant ratios are 0.165–0.397, below the linear
instability threshold; this check did not implicate muscle integration instability.

## Completed bounded tests

- 48 Rustuna timing/recruitment trials: 351.37 s, no full shuffle pass.
- 48 output-adapter calibration trials: 358.11 s, no full pass. These optionally
  scale MN-output voltage ranges, not muscle force or strength parameters.
- 10 soft hip-excursion/PF-relay trials, including a baseline: no full pass.
- 10 extensor-burst-duration/PF-relay trials: no full pass.
- 10 transient, load-confirmed heel-reset trials: no full pass.
- Four causal heel controls, including baseline: no full pass.

These 130 trials are separate from the earlier 45 pelvic-overlap comparisons.
The new slow candidate also has one fresh-process repeat. Failed outputs remain
in their owned external folders; raw NPZ/model/log files were not put in the repo.

## Does the contact pathway need an interneuron?

The current heel→contralateral flexor-RG pathway is excitatory. A functional
sensory-gating subnetwork could combine contact, acceptance of load, and phase
permission, then generate one transient reset rather than sustained flexor
excitation. Descending modulation could adjust its gain. A relay alone controls
gain; onset memory/adaptation requires additional state or cells. The tested
software gate is not yet a constructed neuronal subnetwork.

However, the causal controls do not establish excessive heel reset as the main
problem. Actually deleting the two contra-flexor reset edges worsens load loss
and overlap; reducing their gains also does not produce the target shuffle.
Ipsilateral extensor/PF projections were retained. The deletion removes the
edges, rather than relying on zero conductance. Thus broad inhibition should
not be adopted as a demonstrated fix.

Rybak et al. (2024), equations 14–15, model sensory-gain suppression by alpha
as SF_i = k_i E / (1 + k_PSi alpha_ipsi). Gamma supplies extensor drive in that
model; it is not the same sensory-suppression parameter. Applying this mechanism
to our heel→contralateral RG-F branch would be an explicit engineering extension,
not a replication of that paper's hip-length/load pathways. A proposed test is
local gating alone, descending gain modulation alone, both, and unchanged
baseline, with all other weights held fixed. The foot input remains excitatory;
descending inhibition targets its gate, rather than introducing inhibitory heel
feedback to both RG antagonists.

Source: https://elifesciences.org/articles/98841 (Table 1 and equations 14–15).
No descending-gated neuronal circuit has been declared trained or deployed by
these software-gate tests. All trials in this batch have completed.

## Full paths on easteregg2

Root: `D:\Users\Ben Bolen\NeuroCampaign_20261031\`.

- `W2L_short_slow_shuffle_20261104`: LAUNCH, MANIFEST, trials, SUMMARY,
  BELT_DIRECTION, FRESH_VERIFICATION, and source_snapshot.
- Exact slow model: `W2L_short_slow_shuffle_20261104\fresh_repeat\W2L_VirtualWalker_Newer\repeat\model.xml`.
- Same recording: that folder's trajectory.npz.
- Arrow GIF/still: that folder's SUPPORTED_STEPPING_REPLAY.gif/.png
  (the still is SUPPORTED_STEPPING_STILL.png).
- `W2L_shuffle_readout_20261104`: second Rustuna study and source snapshot.
- `W2L_shuffle_PF_limits_20261104`: soft sensory limits.
- `W2L_shuffle_extensor_duty_20261104`: duty-profile comparisons.
- `W2L_shuffle_heel_handover_20261104`: event-gate comparisons.
- `W2L_shuffle_heel_causal_20261104`: actual deletion/reduced-gain controls.
- `W2L_native_FK_sign_audit_20261104\RESULT.json`: native geometry comparison.

Portable media/record: the existing
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104`
folder, as LiDeng_W2L108_slow_direction.gif/.png and this record.
