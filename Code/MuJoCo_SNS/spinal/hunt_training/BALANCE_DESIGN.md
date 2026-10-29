# Stage 3 balance: source interpretation and next controller tests

The experimental pathway is now implemented in balance.py and tested through
physical_eval.py. This is a Peterka/Ting-inspired engineering adaptation, not
a replication or validated physiological balance model. Distributed harness
restrictions remain assistance. Results and backend comparisons are recorded
separately in BALANCE_RESULTS.md.

## Implemented pathway

Four directional leaky SNS populations (anterior, posterior, left, and right)
and four inhibitory relay neurons project to native muscle MNs. Signals combine
delayed CoM displacement and velocity; filtered acceleration is independently
optional. Default delay is 100 ms, velocity/acceleration filter 40 ms, and task
and relay membrane time constants 20 ms. The implemented cue mixture combines
true simulated CoM error and an orientation-based inverted-pendulum proxy with
an adjustable, fixed weight. It does not implement adaptive sensory reweighting,
vision, or biological sensory transduction.

Recruitment is calibrated once per starting pose using the native model's
finite activation-to-CoM-acceleration response with ground contact. Four bounded
nonnegative activation bases approximate opposite horizontal accelerations
while discouraging vertical cross-effects. Their weights become fixed SNS
conductances; opposing recruitment reaches MNs through inhibitory relays.
There is no runtime inverse-dynamics allocator or MN voltage cancellation.

Replacement and additive modes are tested explicitly. The additive mode retains
the inherited posture, trunk, and balance pathways. An identical-topology
disabled-task control isolates the added pathway. The inherited controller
passes its initial standing screen after correcting the initial stance; that
improvement cannot be attributed solely to the new pathway.

Tuning penalizes CoM position error, torso tilt, measured horizontal harness
reaction, insufficient foot loading, and activation effort. The acceptance gate
is separate and unchanged. Low restraint force alone is not accepted as balance.

## Standing versus stepping

For quiet standing, regulate center-of-mass (CoM) displacement and velocity
relative to the current foot support area. Keep bilateral ground loading,
bounded joint/torso motion, and low measured harness reaction as separate gates.
Being upright inside a stiff harness is not independent balance. Low sway alone
can indicate a restrictive harness rather than a successful neural controller.

For walking, a forward CoM trajectory and controlled redirection through ground
reaction forces and foot placement are appropriate. Hip flexion contributes to
swing advance, but it is not the sole mechanism that prevents falling. Do not
train stage 3 by rewarding forward lean or continuous stepping. A recovery step
is a separately scored balance-recovery task.

## Primary sources and what they justify

- Peterka (2002), Sensorimotor Integration in Human Postural Control,
  https://doi.org/10.1152/jn.2002.88.3.1097. Human visual/support-surface rotations
  support a feedback model with sensory reweighting, stiffness, damping, and
  delay. This motivates explicit, separately perturbable orientation cues;
  it does not prescribe a universal gain set for this SNS or a literal circuit.
- Welch and Ting (2008), A Feedback Model Reproduces Muscle Activity During
  Human Postural Responses to Support-Surface Translations,
  https://doi.org/10.1152/jn.01110.2007. Delayed CoM displacement, velocity, and
  acceleration can describe temporal muscle recruitment. This motivates testing
  a task-level error pathway, rather than joint-angle feedback alone.
- Welch and Ting (2009), https://doi.org/10.1152/jn.90775.2008. Initial and later
  components of human postural EMG scale differently with perturbation
  acceleration and velocity. Use finite perturbations of distinct magnitude
  and duration; do not infer robust feedback from one quiet-stance trial.
- Ting and Macpherson (2005), A Limited Set of Muscle Synergies for Force
  Control During a Postural Task, https://doi.org/10.1152/jn.00681.2004.
  Cat postural responses support directionally organized muscle synergies.
  This motivates a recruitment basis fitted from force-producing functions;
  it does not establish that gait PF populations and postural synergies are
  identical or that a particular number transfers to the human model.
- Hof et al. (2007), Control of lateral balance in walking,
  https://pubmed.ncbi.nlm.nih.gov/16740390/. Velocity-sensitive extrapolated CoM
  and foot placement are useful for dynamic lateral balance. The inverted-
  pendulum approximation needs care under strong harness forces and distributed
  restrictions; do not present its unmodified margin as proof of stability.

## Candidate implementation and verification

1. Log CoM state, foot contact points/forces, pelvis orientation, and measured
   harness wrench. Plot lateral restraint force against stance/swing timing.
   Include out-of-plane hip/lumbar restraint torques in the assistance audit.
2. Build separate anterior/posterior and right/left task-error channels with
   adjustable delayed CoM position/velocity feedback; test filtered acceleration
   separately because numerical differentiation and feedback delay can amplify
   noise. Report the chosen delay and filter, not only the gain.
3. Fit recruitment against actual muscle moment arms and ground-force effects,
   with antagonist inhibition through interneurons. Distinguish this standing
   task-level pathway from the timed RG/PF locomotor output. Retain all 3D axes.
4. Compare inherited feedback, added CoM feedback, and pathway ablations at
   identical harness settings. Hold perturbation trials out of gain fitting.
5. Reduce vertical support and lateral restraint separately after passing
   ground-load/posture gates. Use known finite forward/backward and lateral
   pulses, then reproducible band-limited disturbances. A disturbance trained
   and tested at the same instant alone is not a robustness test.

The virtual pelvis wrench is derived from springs/dampers and generalized
joint torques; it is not yet a physical one-attachment six-axis transducer.
External applied perturbations and explicit support forces are logged separately.
