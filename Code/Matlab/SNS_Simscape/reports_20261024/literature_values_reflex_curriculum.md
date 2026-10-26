# Literature values for the SNS demo suite + curriculum — first pass (2026-10-24)

Ben's order (2026-10-24): ground the demos in actual studies, run the
SNS-toolbox and AnimatLab tutorials and match their responses, and define
the five-part curriculum. Provenance convention: **L** = literature value
sourced from the web today (URL given; full-text verification pending —
do not cite in the dissertation until verified against the paper).

## 1. Stretch reflex / load compensation (beer-cup demo target)

- **L** Nichols & Houk 1976, "Improvement in linearity and regulation of
  stiffness that results from actions of stretch reflex", J Neurophysiol
  (https://pubmed.ncbi.nlm.nih.gov/1249597): reflex action COMPENSATES for
  variations in intrinsic muscle stiffness — the reflex is a stiffness
  regulator, not just added force.
- **L** Sinkjaer et al. 1988, human ankle dorsiflexors
  (https://pubmed.ncbi.nlm.nih.gov/3171659): reflex contribution to total
  stiffness begins ~50 ms after stretch onset, PEAKS 150–300 ms after.
- **L** Mirbagheri et al. 2000
  (https://pubmed.ncbi.nlm.nih.gov/11156307): reflex stiffness up to ~55%
  of total human ankle stiffness; typical literature range ~30–60%
  depending on activation and position.

Our beer-cup numbers (dynamically reproduced 2026-10-24): sag reflex-OFF
7.85° vs ON 2.00° under the same pour → apparent stiffness ratio ≈ 3.9×,
i.e. the reflex supplies ≈74% of total steady stiffness. That is somewhat
above the human literature range (30–60%) — plausibly because the demo's
biceps series compliance (0.02 m/rad) is deliberately soft "so there is
real work for the reflex." First-order conclusion: same order of
magnitude, same mechanism (finite-gain servo); to claim more, tune the
compliance to a human musculotendon value and re-derive the fraction.
Latency comparison: demo reflex path ≈ SN τ 16.7 ms + MN τ 10 ms + BPA
activation τ 150 ms → force response ~150-200 ms, inside Sinkjaer's
150–300 ms functional-stiffness peak window (but the neural part alone,
~27 ms, matches the short-latency M1 response).

## 2. Air stepping / infant stepping (curriculum stages 1–2 target)

- **L** Thelen 1986, treadmill-elicited stepping in 7-month-olds
  (https://pubmed.ncbi.nlm.nih.gov/3802974): ALL infants showed immediate
  alternating stepping when supported over a treadmill — pattern present,
  expression gated by biomechanics.
- **L** Thelen et al. 1984, newborn stepping and physical growth
  (https://pzacad.pitzer.edu/~dmoore/psych152_S21/21_1984_thelen-et-al_relati.pdf):
  adding small leg weights suppresses stepping; water immersion (reducing
  effective weight) restores it → strength-to-weight ratio gates stepping.
  DIRECT design guidance: air stepping = unloaded-limb expression of the
  same CPG the harness curriculum later loads.
- **L** Teulier et al. 2012/2015
  (https://pmc.ncbi.nlm.nih.gov/articles/PMC3424088): only 70–84% of
  healthy newborns respond with step-like patterns; high variability.
  Expect variability across seeds/initial states, not a single universal
  response.

## 3. Harness-supported locomotion (curriculum stages 3a–4 target)

- **L** Clinical Lokomat/BWSTT protocols: initial body-weight support
  typically 30–50%, progressively reduced (Hocoma Lokomat
  https://www.hocoma.com/us/solutions/lokomat; clinical descriptions
  https://www.umms.org/shore/health-services/rehabilitation/outpatient-services/physical-therapy/lokomat-robotic-gait;
  Westlake & Pohl 2009 https://pubmed.ncbi.nlm.nih.gov/19523207). The
  devices guide the legs in the sagittal plane — matching Ben's note that
  literature harnesses often constrain to sagittal-plane locomotion.
  Curriculum harness design should copy this: %BWS reduced in stages,
  sagittal constraint, no torques injected by the harness itself.

## 4. Tutorial targets to match

- SNS-toolbox tutorials 1–9 (repo: `Code\MuJoCo_SNS\spinal\sns_tutorials`;
  fresh re-runs of 1/2/4/8 in progress 2026-10-24). Existing anchor: the
  units test already matches the toolbox numpy backend on a tutorial-2-
  style 2-neuron circuit at 6.07e-04 mV; extend to the other tutorial
  circuits.
- **L** AnimatLab canonical tutorial: the official **Human Stretch
  Reflex** example (video tutorial,
  https://www.animatlab.com/Help/Tutorials/Examples/Human-Stretch-Reflex)
  and Szczecinski's "Central control components of a 'simple' stretch
  reflex" (ModelDB, https://modeldb.science/modellist/150700?all_simu=true);
  the founding AnimatLab paper is Cofer et al. 2010
  (https://www.sciencedirect.com/science/article/abs/pii/S0165027010000087),
  whose demo is arm flexion with myotatic + withdrawal reflexes. Plan:
  obtain the tutorial project, run headless on EB475WS4
  (AnimatSimulator.exe), extract traces, build the Simulink SNS twin,
  match responses.
