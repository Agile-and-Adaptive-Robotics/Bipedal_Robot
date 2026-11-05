# W2L108 pelvic posture and double-support tests

Machine: **easteregg2**, 4 November 2026. Ben observed goose stepping and proposed
anterior pelvic tilt and more double support. Forty-five bounded 10 s comparisons
were run, followed by two independent repeats. Raw trials remain outside the repo.

The observation is supported by the data: the previous winner has only 1.51%
double support after 1 s, with both feet above 0.05 BW normal force. Its overlap
bouts last only 2–18 ms. Both feet are below that load threshold for 18.16% of the
window. Swing knee flexion averages 10.91° over the existing qualified swing
intervals; 62.1% of those samples are below 10°. Passing the unloaded-clear
alternation gate did not establish a satisfactory walking waveform.

## Results and tradeoffs

| Condition | Mean anterior pelvic tilt | Double support | Upward assistance | Swing knee flexion |
|---|---:|---:|---:|---:|
| Previous winner | −7.39° | 1.51% | 29.52% BW | 10.91° |
| Posture PF, kp=1, target 0° | −3.88° | 2.80% | 21.13% BW | 9.95° |
| Posture PF, kp=1.4, target 0° | +4.16° | 10.04% | 32.63% BW | 10.03° |

All percentages use the same 1–10 s window. Positive tilt is anterior, measured
from the actual root orientation, not from the command. The two posture cases
repeat every recorded numeric array byte-exactly in new processes and retain
positive anatomical knee flexion. They are **tradeoff examples, not newly
accepted best gaits**. None of the tested conditions meets the predeclared
combined improvement gate without excess assistance or other shortcomings.

The anterior/overlap case sways from −15.59° to +27.50°, so the favorable mean
does not establish stable 4° posture regulation. Its loaded-overlap bouts include
70–174 ms intervals, but many remain brief impact events. Restricting the
descriptive calculation to bouts lasting at least 40 ms gives 4.82% of the
window, versus zero in the baseline. This extra statistic was added after the
tests as a diagnostic, not substituted into the saved acceptance criteria.
Knee motion remains too straight in swing. The reduced-support case preserves
8/7 swings and full qualified alternation; the anterior case has 6/6 swings and
90.9% qualified alternation. Neither is silently substituted for the prior winner.

Initial 2°/4° anterior root rotations with compensating hip coordinates did not
maintain anterior tilt. Shared extensor RG tonic inputs, contact→PF gain,
V3 coupling, modest knee-flexor gain, and a slower neural clock did not provide
the required combined improvement. Lowering the rails by 5–30 mm increased
support and posterior lean; those explicit apparatus variants also failed.

## Posture pathway actually tested

Simulated IMU-like root tilt and angular rate → explicit PD task command →
measured stance gate → existing antagonistic hip PF flexor/extensor cells →
existing MN/muscle pathways. Gains are shared across legs; each side's stance
gate uses its own measured ground force. There is no direct posture→MN input,
root torque, added mechanical posture spring, or new postural neuron population.
This is a diagnostic use of the existing PF pair, not a completed neural task
estimator or the full proposed Gait2392 postural subnetwork. Commands and sensed
states are recorded in trajectory.npz/pelvis_pf_posture. A zero option is the
unchanged baseline. Additional angular-rate feedback is a controller change;
joint/muscle damping, tendon law, and muscle physics were not changed.

The first posture ladder uses kd=0.15 nA/(rad/s), current cap 0.6 nA, and starts
at 1 s. The follow-up tests higher derivative gains and a 0.3 nA cap. These
reduced some effects but did not produce an accepted combined improvement.
The remaining work is stable pelvic regulation and stance/swing waveform
shaping, with sustained overlap, flight fraction, and knee flexion in the cost
and acceptance tests, rather than optimizing rail support alone.

## Full paths

All roots below are on easteregg2 under
`D:\Users\Ben Bolen\NeuroCampaign_20261031\`:

- `W2L_pose_overlap_20261104`: initial tilt, extensor drive, hip-flexor comparisons.
- `W2L_overlap_feedback_20261104`: contact/PF, commissural, knee, and clock tests.
- `W2L_pelvis_PF_posture_20261104`: initial posture-command ladder and low-support repeat.
- `W2L_pelvis_PF_refine_20261104`: intermediate posture gains and anterior/overlap repeat.
- `W2L_pelvis_PF_derivative_20261104`: additional velocity-feedback/current-limit tests.
- `W2L_rail_height_overlap_20261104`: explicitly changed rail heights.

Each root contains MANIFEST.json, STATUS.json, COMPLETE.json, and trials.
The two repeated roots additionally contain FRESH_REPEAT.json and
`fresh_repeat\W2L_VirtualWalker_Newer\repeat\model.xml`, trajectory.npz,
SUPPORTED_STEPPING_REPLAY.gif, and SUPPORTED_STEPPING_STILL.png.
For the rail ladder, per-case configurations specify the actual rail-height
changes; the inherited generic manifest phrase 'fixed rail/guide' does not
describe those explicitly changed-height variants.

Portable diagnostic GIFs and this record are in
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104`.
Use LiDeng_W2L108_pelvis_overlap_trial.gif for the anterior/overlap tradeoff and
LiDeng_W2L108_lower_support_posture_trial.gif for the lower-support alternative.
Recorded 5–8.2 s previews are subsets of the trials. Review whole-trial traces
and the sway/assistance figures before interpreting a preview.
