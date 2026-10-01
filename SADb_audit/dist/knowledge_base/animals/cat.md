# Animal: Cat

219 papers in the corpus.

- **Rahmati 2025** — [Role of forelimb morphology in muscle sensorimotor functions during locomotion in the cat](https://doi.org/10.1113/JP287448)  
  - animals: Cat  
  - Measured 46 cat forelimb muscles plus walking mechanics and EMG to compute force- and length-dependent afferent patterns: moment arm, PCSA, and fascicle length strongly shape muscle forces and proprioceptive signals, with forelimb morphology contributing mainly to lateral stability and turning control rather than propulsion — with direct implications for cervical afferent mapping in neuromechanical models.
- **Klishko 2025** — [Effects of spinal transection and locomotor speed on muscle synergies of the cat hindlimb](https://doi.org/10.1113/JP288089)  
  - animals: Cat  
  - Fifteen hindlimb muscle EMGs in cats across 0.4–1.0 m/s before and after low-thoracic transection factor into five synergies in every condition; synergy number, composition, and activation patterns survive transection and speed changes, supporting spinal-level synergy control and grounding a pattern-formation network organization for the two-level CPG.
- **Yassine 2025** — [Speed-dependent locomotor adjustments following staggered thoracic lateral hemisections in adult cats](https://doi.org/10.1152/jn.00331.2025)  
  - animals: Cat  
  - Staggered thoracic lateral hemisections in cats shift the ipsilesional hindlimb from stance/extensor dominance toward swing/flexor while the contralateral limb compensates with prolonged stance and reduced swing, and forelimbs take extra steps independent of speed or lesion side — lesion signatures that reorganized sensorimotor network models must reproduce.  
  - *Robot/sim:* Model staggered hemisection as unilateral loss of descending drive plus partial crossed-connection removal in a quadruped CPG; test whether ipsilesional stance shortening, contralateral stance prolongation, and extra forelimb steps emerge without retuning, and how they scale with speed.
- **Shinohara 2025** — [Mechanisms of adaptive interlimb coordination to sudden ground loss: a neuromusculoskeletal modeling study](https://doi.org/10.1101/2025.11.11.687930)  
  - animals: Cat  
  - A neuromusculoskeletal cat model with two-level half-center CPGs per hindlimb reproduces the adaptive interlimb coordination observed when spinalized cats suddenly lose ground support — without any parameter re-optimization — and nullcline analysis attributes the robustness to afferent feedback controlling the transitions between fast and slow neuronal dynamics.  
  - *Robot/sim:* Implement paired two-level half-center CPGs whose fast-slow transitions are controlled by afferent feedback, tuned on flat ground only; drop the foot into a hole without re-optimization and check that adaptive interlimb coordination emerges, then analyze the transition nullclines.
- **Rybak 2025** — [Operation of spinal sensorimotor circuits controlling phase durations during tied-belt and split-belt locomoti](https://doi.org/10.7554/elife.103504)  
  - animals: Cat  
  - Combining a computational model of cat locomotion with experiments shows that after thoracic lateral hemisection the contralesional ('intact') side remains governed mainly by supraspinal drives while the ipsilesional side is dominated by somatosensory feedback; simulations correctly predicted changes in cycle and phase durations during tied-belt and split-belt locomotion, and placing the ipsilesional hindlimb on the slow belt substantially reduced the effects of hemisection.  
  - *Robot/sim:* Simulate hemisection in a two-level CPG model by removing descending drive to one side; test whether ipsilesional stepping is sustained by somatosensory feedback alone, and reproduce the reduction of hemisection effects when the ipsilesional limb steps on the slow belt during split-belt locomotion.
- **Mari 2024** — [Changes in intra- and interlimb reflexes from hindlimb cutaneous afferents after staggered thoracic lateral he](https://doi.org/10.1113/jp286151)  
  - animals: Cat · afferents: Cutaneous  
  - In cats receiving staggered thoracic lateral hemisections on opposite sides of the cord, short-latency cutaneous reflex responses from superficial peroneal stimulation were largely preserved in homonymous and crossed hindlimb muscles, while mid- and long-latency homonymous and crossed responses in both hindlimbs, and forelimb responses, occurred less frequently; whenever responses were present, all latency classes retained their phase-dependent modulation. Loss of longer-latency cutaneous transmission therefore tracks impaired balance and interlimb coordination after spinal cord injury, making  
  - *Robot/sim:* Implement short-, mid-, and long-latency cutaneous reflex arcs with phase-dependent gating in a quadruped model; ablating the long-latency inter-enlargement transmission should degrade interlimb coordination while sparing short-latency responses.
- **Mari 2024** — [Changes in intra- and interlimb reflexes from forelimb cutaneous afferents after staggered thoracic lateral he](https://doi.org/10.1113/jp286808)  
  - animals: Cat · afferents: Cutaneous  
  - In cats during quadrupedal locomotion, staggered thoracic lateral hemisections largely spare short-, mid- and long-latency homonymous and crossed forelimb reflexes from superficial radial nerve stimulation together with their phase modulation, but significantly reduce or abolish mid- and long-latency homolateral and diagonal responses in hindlimb muscles. This shows a considerable loss of cutaneous reflex transmission from cervical to lumbar levels after incomplete spinal cord injury despite preserved phase modulation, likely impairing coordinated responses to external perturbations across the  
  - *Robot/sim:* Model phase-gated cutaneous reflex pathways with interlimb homolateral and diagonal projections onto all four limbs; ablate the cervico-lumbar links to reproduce the loss of hindlimb responses and weakened interlimb coordination after incomplete spinal cord injury.
- **Rybak 2024** — [Operation regimes of spinal circuits controlling locomotion and the role of supraspinal drives and sensory fee](https://doi.org/10.7554/elife.98841)  
  - animals: Cat, Mammals  
  - Computational model reproducing intact and spinal-transected cats in tied- and split-belt locomotion, showing the spinal network changes operating regime with speed: a non-oscillating state-machine at slow speeds (<0.4 m/s) whose phase transitions require sensory feedback related to limb extension (removing it prevents slow-speed oscillation), a flexor-driven oscillatory regime at intermediate speeds, and a classical half-center regime at high speeds; after transection only the state-machine regime remains.  
  - *Robot/sim:* Build the tri-regime spinal controller (state-machine with sensory phase transitions at slow speed, flexor-driven, half-center) in a neuromechanical simulation; ablate the limb-extension feedback term to confirm loss of slow-speed oscillation and lock the spinalized configuration to state-machine behavior.
- **Harnie 2024** — [Forelimb movements contribute to hindlimb cutaneous reflexes during locomotion in cats](https://doi.org/10.1152/jn.00104.2024)  
  - animals: Cat · afferents: Cutaneous  
  - During quadrupedal locomotion in cats, forelimb movement modulates hindlimb cutaneous reflexes evoked by superficial peroneal stimulation, particularly the occurrence of long-latency responses — evidence that interlimb (forelimb-hindlimb) signals gate spinal cutaneous reflex pathways during walking.  
  - *Robot/sim:* Add forelimb-state-dependent gain modulation to hindlimb cutaneous reflex pathways in a quadruped neuromechanical model; ablating the modulation should remove the long-latency reflex gating present during quadrupedal walking.
- **Mari et al. 2023** — [A sensory signal related to left-right symmetry modulates intra- and interlimb cutaneous reflexes during locom](https://doi.org/10.3389/fnsys.2023.1199079)  
  - animals: Cat · pathways: Cutaneous stance modification  
  - A sensory signal related to left-right symmetry modulates intra- and interlimb cutaneous reflexes during locomotion in intact cats, extending the split-belt symmetry findings to intact walking.
- **Li et al. 2023** — [Identified interneurons contributing to locomotion in mammals](https://doi.org/10.1016/b978-0-12-819260-3.00009-3)  
  - animals: Mammals, Mice, Cat  
  - Chapter cataloguing identified interneuron classes contributing to mammalian locomotion and their synaptic relations; useful companion to Sengupta & Bagnall 2023.
- **Audet 2023** — [Spinal Sensorimotor Circuits Play a Prominent Role in Hindlimb Locomotor Recovery after Staggered Thoracic Lat](https://doi.org/10.1523/ENEURO.0191-23.2023)  
  - animals: Cat  
  - Staggered thoracic lateral hemisections in cats: hindlimb locomotion recovers spontaneously but forelimb–hindlimb coordination degrades into weaker 2:1 patterns, left–right stance/swing asymmetries appear then reverse after the second hemisection, and posture plus interlimb coordination remain impaired — lumbar sensorimotor circuits carry hindlimb recovery but not supraspinal-dependent coordination.
- **Audet 2022** — [Control of fore- and hindlimb movements and their coordination during quadrupedal locomotion across speeds in ](https://doi.org/10.1089/neu.2022.0042)  
  - animals: Cat  
  - Quadrupedal fore-hindlimb coordination across speeds in adult spinal cats: the isolated spinal cord replays speed-dependent coordination, including gallop-like patterns.
- **Kohler 2022** — [Diversified physiological sensory input connectivity questions the existence of distinct classes of spinal int](https://doi.org/10.1016/j.isci.2022.104083)  
  - animals: Cat  
  - The spinal cord is engaged in all forms of motor performance but its functions are
far from understood. Because network connectivity defines function, we
explored the connectivity of muscular, tendon, and tactile sensory inputs among
a wide population of spinal interneurons in the lower cervical segments. Using
low noise intracellular whole cell recordings in the decerebrated, non-anesthetized cat in vivo, we could define mono-, di-, and trisynaptic inputs as well as
the weights of each input. Whereas each neuron had a highly specific input, and
each indirect input could moreover be explained 
- **Harnie 2022** — [State- and Condition-Dependent Modulation of the Hindlimb Locomotor Pattern in Intact and Spinal Cats Across S](https://doi.org/10.3389/fnsys.2022.814028)  
  - animals: Cat  
  - Comparing quadrupedal and hindlimb-only locomotion in the same intact cats before and after spinal transection separates state-dependent from condition-dependent changes in the hindlimb pattern: the spinal state produces convergence of stance and swing durations at high speed, improper ankle-hip coordination, and altered flexor burst timing, whereas the hindlimb-only condition itself shifts paw placement, yield magnitude, and some burst durations - so hindlimb-only intact locomotion, not quadrupedal, is the proper baseline for spinal-state comparisons.  
  - *Robot/sim:* When modeling spinal-state gait changes, validate against a hindlimb-only intact baseline rather than quadrupedal locomotion; implement intact versus spinal (supraspinal drive removed) controller variants and separate state-dependent effects (stance/swing convergence at speed, ankle-hip coordination) from condition-dependent ones.
- **Kim 2022** — [Contribution of Afferent Feedback to Adaptive Hindlimb Walking in Cats: A Neuromusculoskeletal Modeling Study](https://doi.org/10.3389/fbioe.2022.825149)  
  - animals: Cat  
  - A cat hindlimb neuromusculoskeletal model coupled to a half-center CPG reproduces normal walking and the adaptive movement changes when the foot steps into a hole; dynamical-systems and nullcline analysis of the coupled CPG-musculoskeleton-environment loop shows how afferent feedback mediates adaptive locomotion in ways that immobilized fictive preparations cannot reveal.  
  - *Robot/sim:* Reproduce the architecture — half-center CPG coupled to a cat hindlimb musculoskeletal model with limb afferent feedback — and ablate each afferent channel during hole-step perturbations to map which feedback is necessary for the adaptive response.
- **Merlet 2021** — [Cutaneous inputs from perineal region facilitate spinal locomotor activity and modulate cutaneous reflexes fro](https://doi.org/10.1002/jnr.24791)  
  - animals: Cat, Mammals · afferents: Cutaneous  
  - In spinal cats, mechanical perineal stimulation triggers and intensifies rhythmic hindlimb activity - shortening cycle and burst durations and raising flexor and extensor amplitudes - while simultaneously decreasing short-latency ipsilateral and contralateral cutaneous reflexes across joints and limbs, indicating that perineal facilitation of locomotion and weight support acts by increasing the excitability of CPG circuitry through state-dependent interneuronal modulation rather than by increasing foot cutaneous afferent excitation.  
  - *Robot/sim:* Model perineal input as a tonic excitability gain on CPG neurons paired with a divisive reduction of cutaneous reflex pathway gains; ablating the reflex-gain reduction tests whether unmodulated cutaneous feedback destabilizes the evoked rhythm.
- **Klishko 2021** — [Common and distinct muscle synergies during level and slope locomotion in the cat.](https://doi.org/10.1152/jn.00310.2020)  
  - animals: Cat  
  - Cat downslope walking's atypical EMG (silent one-joint hip extensors, stance-related flexor bursts) shares the majority of its burst groups and muscle synergies with level and upslope walking, and slope changes swing/stance phase durations but not cycle duration — consistent with one shared CPG whose output is reshaped by somatosensory and supraspinal inputs rather than task-specific circuit reconfiguration.  
  - *Robot/sim:* Drive a fixed muscle-synergy output layer from a single CPG whose phase-duration scaling reproduces slope-dependent EMG; verify the model predicts the downslope silent hip extensors and stance flexor bursts without changing synergy vectors.
- **Zholudeva 2021** — [Spinal Interneurons as Gatekeepers to Neuroplasticity after Injury or Disease.](https://doi.org/10.1523/jneurosci.1654-20.2020)  
  - animals: Mice, Cat, Rat, Human, Vertebrates  
  - Review positioning spinal interneurons — a heterogeneous population that modulates motor, sensory, and autonomic function — as key components of plasticity and recovery after spinal cord injury, and arguing that treatments should be optimized for how they engage interneuron circuits rather than viewed as acting around them.
- **Akay and Murray 2021** — [Relative Contribution of Proprioceptive and Vestibular Sensory Systems to Locomotion: Opportunities for Discov](https://doi.org/10.3390/ijms22031467)  
  - animals: Cat · pathways: Ia or II stance to swing  
  - Hypothesis: The Weight of the role of segmental feedback is less important at slower speeds but increases at higher speeds, where as in teh vestibular system has the opposite effect

Mouse model is key tool to undertadning terrestial locmotion (Pg 1)

Vestivular feedback being critical at lower speed and somatosensory feedback is necessary at higher velocites. (Pg1)


CPG works with both sensory feedback from the leg or supraspinal center to generate locomotor patter that is able to deal with changes in terrain (Pg2)

The foot is on the ground during the stance phase and moves in the opposite 
- **Merlet 2021** — [Inhibition and Facilitation of the Spinal Locomotor Central Pattern Generator and Reflex Circuits by Somatosen](https://doi.org/10.3389/fnins.2021.720542)  
  - animals: Cat, Mammals · afferents: Cutaneous  
  - Review showing that in mammals with complete spinal cord injury, somatosensory input from the lumbar region powerfully inhibits hindlimb locomotion while perineal input facilitates it, and the two regions oppositely regulate cutaneous reflexes from the foot (lumbar input increases reflex gain, perineal input decreases it); spinal cord injury can also cause loss of functional specificity, with somatosensory feedback abnormally co-activating functions such as locomotion and micturition.  
  - *Robot/sim:* Add region-specific somatosensory gating to a spinal locomotor model: lumbar input inhibits and perineal input facilitates the CPG, with opposite modulation of cutaneous reflex gain; simulating complete spinal-cord injury then tests emergence of maladaptive co-activation of locomotion with other spinal functions.
- **Domínguez-Rodríguez 2020** — [Candidate Interneurons Mediating the Resetting of the Locomotor Rhythm by Extensor Group I Afferents in the Ca](https://doi.org/10.1016/j.neuroscience.2020.09.017)  
  - animals: Cat · afferents: Ib · pathways: Ib stance to swing; Ib swing to stance; Ib excitatory  
  - Extensor group I afferent stimulation resets fictive locomotion to extension - prolonging ongoing extension and terminating ongoing flexion - through a polysynaptic excitation of extensor motoneurons (latencies near 3.5-4.0 ms, compatible with three interposed interneurons) that replaces the classical Ib non-reciprocal inhibition during locomotion. Extension-phase interneurons receiving short-latency group I excitation satisfy the criteria for the pathway that resets the rhythm and may belong to the rhythm-generating layer of the CPG.  
  - *Robot/sim:* Implement extensor group I load afferents as inputs to rhythm-layer interneurons that both prolong the extensor half-center and reset flexion to extension; ablating them should remove load-dependent stance prolongation and rhythm resetting.
- **Harnie 2020** — [The Spinal Control of Backward Locomotion.](https://doi.org/10.1523/jneurosci.0816-20.2020)  
  - animals: Cat  
  - After complete spinal transection, cats can generate backward walking - given tonic perineal somatosensory facilitation - using speed-modulation strategies, muscle activation patterns, and five muscle synergies shared with forward locomotion, including split-belt stepping; the spinal locomotor network is therefore shared between directions, with limb sensory feedback controlling direction and backward locomotion requiring greater spinal excitability.  
  - *Robot/sim:* Implement one CPG whose phase-sequence direction is selected by limb afferent signals, with a raised excitability threshold for the backward mode, and test whether a fixed synergy basis serves both directions.
- **Latash 2020** — [On the organization of the locomotor CPG: insights from split-belt locomotion and mathematical modeling](https://doi.org/10.1101/2020.07.17.205351)  
  - animals: Cat  
  - Review of organization locomotor cpg: insights. Computational model study. We suggest that the operation of RGs is state -dependent, so that with an increase of external excitation the rhythmogen esis changes from “flexor-driven” oscillations to a “classical half-center” mechanism.
- **Latash 2020** — [On the Organization of the Locomotor CPG: Insights From Split-Belt Locomotion and Mathematical Modeling](https://doi.org/10.3389/fnins.2020.598888)  
  - animals: Cat  
  - Uses split-belt experiments in intact and spinal cats together with mathematical modeling of CPG circuits to support a state-dependent organization: with increasing excitatory input or locomotor speed, the rhythmogenic mechanism of each limb's rhythm generator changes from flexor-driven rhythmicity toward a classical half-center mechanism, with specific commissural interactions between left and right rhythm generators accounting for stance and swing duration changes under symmetric and asymmetric conditions.  
  - *Robot/sim:* Build a two-level CPG with a separate rhythm generator per limb and commissural interactions; reproduce stance and swing duration changes under symmetric and asymmetric (split-belt) drives and verify that the flexor-driven-to-half-center regime switch emerges with increasing drive.
- **Merlet 2020** — [Cutaneous inputs from perineal region facilitates and modulates spinal locomotor activity and reduces cutaneou](https://doi.org/10.1101/2020.07.29.226530)  
  - animals: Cat · afferents: Cutaneous  
  - In spinal cats, mechanical stimulation of the perineal region triggers and re-times rhythmic locomotor activity — shortening cycle and burst durations and increasing flexor and extensor burst amplitudes — while simultaneously decreasing short-latency ipsilateral and contralateral cutaneous reflexes from the foot across joints and limbs. The locomotor facilitation therefore reflects excitation of CPG circuitry with state-dependent reflex suppression by spinal interneurons, not an increase in cutaneous reflex gain.  
  - *Robot/sim:* Implement a tonic cutaneous (perineal-analog) input that raises CPG drive while scaling down short-latency reflex gains; the model should show shorter cycles, larger bursts, and reduced cutaneous reflex responses, testing state-dependent reflex gating layered on a CPG.
- **Akay 2020** — [Sensory Feedback Control of Locomotor Pattern Generation in Cats and Mice.](https://doi.org/10.1016/j.neuroscience.2020.05.008)  
  - animals: Cat, Mice  
  - Review framing sensory feedback control of locomotion across the two dominant models: cat experiments established the broad picture of phase- and task-dependent reflex modulation during walking, while mouse molecular genetics now enables population-specific deletion and modulation to assign those roles to identified afferent and interneuron classes.  
  - *Robot/sim:* Use the cat-derived sensory-control map to choose which afferent pathways to implement in a simulated walker, then emulate mouse population deletions as targeted pathway ablations to predict functional losses.
- **Fujiki 2019** — [Phase-Dependent Response to Afferent Stimulation During Fictive Locomotion: A Computational Modeling Study.](https://doi.org/10.3389/fnins.2019.01288)  
  - animals: Cat  
  - A half-center CPG model with a slowly inactivating persistent sodium current reproduces the phase-dependent effects of brief afferent stimulation on fictive locomotion in cats, in which stimulation can shorten or prolong the current locomotor phase and reset the cycle depending on stimulation phase and nerve. Dynamic-systems analysis of the model identifies the mechanisms of locomotor rhythm resetting, supporting phase-dependent afferent access to the CPG as the substrate for sensory correction of the locomotor rhythm.  
  - *Robot/sim:* Implement the persistent-sodium half-center CPG and apply brief phase-tagged stimuli to reproduce the cat phase-response behavior (phase-dependent shortening, prolongation, and resetting); use the model's dynamic-systems structure to analyze the resetting mechanisms.
- **Desrochers 2019** — [Spinal control of muscle synergies for adult mammalian locomotion](https://doi.org/10.1113/JP277018)  
  - animals: Cat  
  - Key points
The control of locomotion is thought to be generated by activating groups of muscles that perform similar actions, which are termed muscle synergies.
Here, we investigated if muscle synergies are controlled at the level of the spinal cord.
We did this by comparing muscle activity in the legs of cats during stepping on a treadmill before and after a complete spinal transection that abolishes commands from the brain.
We show that muscle synergies were maintained following spinal transection, validating the concept that muscle synergies for locomotion are primarily controlled by circui
- **Duysens and Forner-Cordero 2019** — [A controller perspective on biological gait control: Reflexes and central pattern generators](https://doi.org/10.1016/j.arcontrol.2019.04.004)  
  - animals: Human, Cat  
  - Reviews the neural control of gait from a control-engineering perspective and argues the biological CPG is asymmetric: the rhythm generator drives the flexor half in feedforward fashion while the extensor half operates in feedback mode, gated by external inputs such as ground contact — an architecture with direct implications for how biped robots should couple CPGs to contact sensing.
- **Habu 2019** — [Three-dimensional walking of a simulated muscle-driven quadruped robot with neuromorphic two-level central pat](https://doi.org/10.1177/1729881419885288)  
  - animals: Cat  
  - Computational model study. For each leg, we use a two-level central pattern generator consisting of a rhythm generation part to produce basic rhythms and a pattern formation part to synergistically activate a different set of muscles in each of the four sequential phases (swing, touchdown, stance, and liftoff).
- **Hurteau 2018** — [Intralimb and Interlimb Cutaneous Reflexes during Locomotion in the Intact Cat](https://doi.org/10.1523/JNEUROSCI.3288-17.2018)  
  - animals: Cat · pathways: Cutaneous flexor excitation; Cutaneous stance modification  
  - The skin contains receptors that, when activated, send inputs to spinal circuits, signaling a perturbation. Rapid responses, or
reflexes,inmuscles ofthe contacted limb and opposite homologous limb helpmaintain balance andforward progression.Here, we
investigated reflexes during quadrupedal locomotion in the cat by electrically stimulating cutaneous nerves in each of the four
limbs. Functionally, responses appear to modify the trajectory or stabilize the movement of the stimulated limb while modifying
the support phase ofthe other limbs. Reflexes between limbs aremediated byfast-conducting path
- **Nichols 2018** — [Distributed force feedback in the Spinal Cord and the regulation of limb mechanics](https://doi.org/10.1152/jn.00216.2017)  
  - animals: Cat, Human · pathways: Ib disynaptic inhibition; Ib disynaptic excitation; Ia monosynaptic  
  - Review Update Paper; This paper show able to give insight to the inhibitory and excitatory force feedback during locomotion
- **Côté 2018** — [Spinal Control of Locomotion: Individual Neurons, Their Circuits and Functions.](https://doi.org/10.3389/fphys.2018.00784)  
  - animals: Cat, Mice, Rat, Human, Mammals  
  - Review of spinal interneurones carrying muscle-stretch and load afferent information across mammals: rest-state physiology is conserved from rodents to humans, but locomotor activation profiles differ markedly between species, plausibly reflecting species differences in afferent distribution and interneuronal interactions; targeted neuromodulation of these circuits is the emerging rehabilitation strategy after spinal cord injury.
- **Macefield 2018** — [Functional properties of human muscle spindles](https://doi.org/10.1152/jn.00071.2018)  
  - animals: Cat, Human  
  - Review of functional properties human muscle spindles. First published April 18, 2018; doi:10.1152/ jn.00071.2018.—Muscle spindles are ubiquitous encapsulated mechanoreceptors found in most mammalian muscles.
- **Proske and Gandevia 2018** — [Kinesthetic Senses](https://doi.org/10.1002/cphy.c170036)  
  - animals: Human, Cat  
  - The kinesthetic senses are the senses of position and movement of the body, senses we are awareof only on introspection. A method used to study kinesthesia is muscle vibration, which engagesafferents of muscle spindles to trigger illusions of movement and changed position. When vibratingelbow ﬂexors, it generates sensations of forearm extension, when vibrating extensors, sensationsof forearm ﬂexion. Vibrating the elbow joint produces no illusion. Vibrating ﬂexors and extensorstogether at the same frequency also produces no illusion, because what is perceived is the signaldifference between ant
- **Hurteau and Frigon 2018** — [A Spinal Mechanism Related to Left–Right Symmetry Reduces Cutaneous Reflex Modulation Independently of Speed D](https://doi.org/10.1523/JNEUROSCI.1082-18.2018)  
  - animals: Cat · pathways: Cutaneous stance modification  
  - A spinal mechanism related to left-right symmetry reduces cutaneous reflex modulation independently of speed during split-belt locomotion — symmetry detection is a distinct spinal computation, not a speed by-product.
- **Yakovenko et al. 2018** — [Analytical CPG model driven by limb velocity input generates accurate temporal locomotor dynamics.](https://doi.org/10.7287/peerj.preprints.26734v1)  
  - animals: Cat  
  - Analytical CPG model driven by limb velocity input generates accurate temporal locomotor dynamics across speeds — a compact account of speed-dependent phase timing.
- **Hurteau 2017** — [Nonlinear Modulation of Cutaneous Reflexes with Increasing Speed of Locomotion in Spinal Cats](https://doi.org/10.1523/JNEUROSCI.3042-16.2017)  
  - animals: Cat  
  - When walking, receptors located in the skin respond to mechanical pressure and send signals to the CNS to correct the trajectory
of the limb and to reinforce weight support. These signals produce different responses, or reflexes, if they occur when the foot is
contacting the ground or in the air. This is known as phase-dependent modulation of reflexes. However, when walking at faster
speeds, we do not know if and how these reflexes are changed. In the present study, we show that reflexes from the skin are
modulated with speed and that this is controlled at the level of the spinal cord. This mo
- **Vincent 2017** — [Muscle proprioceptors in adult rat: mechanosensory signaling and synapse distribution in spinal cord](https://doi.org/10.1152/jn.00497.2017)  
  - animals: Rat, Cat · afferents: Ia, II, Ib  
  - In adult rats, functionally identified triceps surae proprioceptors (Ia, II, Ib) project and distribute provisional synapses in the spinal cord much as in the cat, but rat Ib afferents fire robustly during passive muscle stretch and Ia afferents show exaggerated dynamic responses even after locomotor scaling is accounted for — mechanosensory coding is species-adapted, so cat-derived afferent transfer functions cannot be transferred to other species without correction.  
  - *Robot/sim:* Instantiate afferent encoders with rat-specific gains (exaggerated Ia dynamic response, stretch-sensitive Ib) rather than cat transfer functions in a rat-scale neuromuscular leg model; test controller sensitivity to the species differences.
- **Kuczynski 2017** — [Lack of adaptation during prolonged split-belt locomotion in the intact and spinal cat.](https://doi.org/10.1113/jp274518)  
  - animals: Cat, Human  
  - Neither intact nor spinal-transected cats adapt to 10 min of split-belt locomotion: step-length and double-support asymmetries persist with no after-effect on returning to tied belts (symmetry is restored immediately), and spinal cats show no EMG modulation while intact cats raise extensor EMG throughout - unlike humans - suggesting that restoring left-right symmetry is not required for balance in quadrupedal gait and that split-belt adaptation is not a cat spinal plasticity.  
  - *Robot/sim:* Run a spinal-only CPG plus reflex quadruped model on split-belt input - it should reproduce the cat phenotype (persistent asymmetry, immediate re-symmetrization, no after-effect), whereas adding a supraspinal adaptive layer is required for the human-like after-effect.
- **Bondy 2016** — [Control of Cat Walking and Paw-Shake by a Multifunctional Central Pattern Generator](https://doi.org/10.1007/978-1-4939-3267-2_12)  
  - animals: Cat, Mammals  
  - A single half-center CPG model with two coexisting activity regimes — fast ~10 Hz paw-shake and slow ~2 Hz walking — drives a neuromechanical cat-hindlimb model to produce either behavior when regime-appropriate spinal synaptic weights are selected; paw skin afferent input is proposed as the trigger for selecting both the CPG regime and the spinal circuitry.
- **Danner et al. 2016** — [Central control of interlimb coordination and speed-dependent gait expression in quadrupeds.](https://doi.org/10.1113/jp272787)  
  - animals: Cat  
  - Computational demonstration that one spinal circuit with commissural and intralimb coupling reorganizes across speeds to express walk, trot, and gallop — gait expression emerges from speed-dependent coordination, not separate CPGs.
- **Frigon 2016** — [Left–right coordination from simple to extreme conditions during split-belt locomotion in the chronic spinal a](https://doi.org/10.1113/JP272740)  
  - animals: Cat  
  - These people make cat zombies
- **Frigon 2015** — [Modulation of forelimb and hindlimb muscle activity during quadrupedal tied-belt and split-belt locomotion in ](https://doi.org/10.1016/j.neuroscience.2014.12.084)  
  - animals: Cat  
  - Cats were chronically implanted for EMG, which was
obtained from six muscles: biceps brachii, triceps brachii,
flexor carpi ulnaris, sartorius, vastus lateralis and medial
gastrocnemius. During tied-belt locomotion, cats stepped
from 0.4 to 1.0 m/s in 0.1 m/s increments whereas during
split-belt locomotion, cats stepped with left–right speed
differences of 0.1 to 0.4 m/s in 0.1 m/s increments. During
tied-belt locomotion, EMG burst durations and mean EMG
amplitudes of all muscles respectively decreased and
increased with increasing speed. During split-belt locomotion, there was a clear differe
- **Hunt 2015** — [A biologically based neural system coordinates the joints and legs of a tetrapod](https://doi.org/10.1088/1748-3190/10/5/055004)  
  - animals: Rat, Cat  
  - A biologically based neural controller built from cat hindlimb pathway data coordinates sagittal-plane trotting in a 14-joint planar rat model actuated by antagonistic Hill muscle pairs; changing the strength of a single inter-leg connection suffices to account for the phase-timing differences observed between individual trotting rats, identifying inter-leg coupling strength as a dominant determinant of interlimb coordination.  
  - *Robot/sim:* Implement cat-derived hind-leg and inter-leg coordination networks on a 14-joint antagonist-Hill-muscle tetrapod model, vary the strength of the single inter-leg connection to reproduce inter-individual trot-timing differences, and ablate inter-leg connections to probe gait breakdown.
- **Dambreville 2015** — [The spinal control of locomotion and step-to-step variability in left-right symmetry from slow to moderate spe](https://doi.org/10.1152/jn.00419.2015)  
  - animals: Cat  
  - During
experiments, EMG and bilateral video recordings were made during
treadmill locomotion from 0.1 to 0.4 m/s in 0.05 m/s increments.
Cycle and stance durations significantly decreased with increasing
speed, whereas swing duration remained unaffected. Extensor burst
duration significantly decreased with increasing speed, whereas sartorius burst duration remained unchanged. Stride length, step length,
and the relative distance of the paw at stance offset significantly
increased with increasing speed, whereas the relative distance at
stance onset and both the temporal and spatial phasing betw
- **Hurteau 2015** — [Effect of stimulating the lumbar skin caudal to a complete spinal cord injury on hindlimb locomotion](https://doi.org/10.1152/jn.00739.2014)  
  - animals: Cat  
  - Mechanical stimulation of the lumbar skin above a complete spinal transection disrupts spinal locomotion in walking cats — pinching stops stepping and abolishes weight support, focal clip stimulation reduces flexor and extensor activity and shifts paw placement caudally, with strongest effects from L3-L5 — showing cutaneous input modulates the excitability of the spinal rhythm-generating segments themselves.
- **Frigon 2014** — [Speed-dependent modulation of phase variations on a step-by-step basis and its impact on the consistency of in](https://doi.org/10.1152/jn.00524.2013)  
  - animals: Cat  
  - In intact walking cats, step-by-step phase variations are speed dependent: at 0.3 m/s variations are stance/extensor dominated and interlimb coordination is least consistent, while at faster speeds stance and swing covary in roughly equal proportion; interlimb phasing shifts with speed for homolateral and diagonal but not homologous limb pairs, and hindlimb phase variation linearly predicts the consistency of interlimb coordination.  
  - *Robot/sim:* Implement a CPG whose stance/swing duration allocation shifts from stance-dominated to covarying with speed, and verify it reproduces the linear relation between hindlimb phase variation and interlimb-coordination consistency across speeds.
- **D’Angelo 2014** — [Modulation of phase durations, phase variations, and temporal coordination of the four limbs during quadrupeda](https://doi.org/10.1152/jn.00160.2014)  
  - animals: Cat  
  - During experiments, EMG and bilateral video recordings were made during
treadmill locomotion from 0.1 to 0.4 m/s in 0.05 m/s increments.
Cycle and stance durations significantly decreased with increasing
speed, whereas swing duration remained unaffected. Extensor burst
duration significantly decreased with increasing speed, whereas sartorius burst duration remained unchanged. Stride length, step length,
and the relative distance of the paw at stance offset significantly
increased with increasing speed, whereas the relative distance at
stance onset and both the temporal and spatial phasing betw
- **Shevtsova 2013** — [Two-Level Model of Mammalian Locomotor CPG](https://doi.org/10.1007/978-1-4614-7320-6_49-2)  
  - animals: Cat  
  - Chapter presentation of the two-level model of the mammalian locomotor CPG: a half-center rhythm generator (RG) drives a separate pattern-formation network (PF) that determines the coordinated motoneuron activity pattern, an architecture inferred from deletion analysis and afferent stimulation during brainstem-evoked fictive locomotion in the decerebrate, immobilized cat.
- **Frigon 2013** — [The Cat Model of Spinal Cord Injury](https://doi.org/10.1007/978-1-62703-197-4_8)  
  - animals: Cat  
  - Synthesis of landmark cat spinal cord injury studies: cats recover robust hindlimb locomotion after partial or complete spinalization, which established the spinal-network basis of locomotor recovery, the adaptive changes taking place after injury, and the treadmill-training principles now applied to rehabilitating spinal cord-injured humans.  
  - *Robot/sim:* Simulate controllers with complete loss of supraspinal drive (spinal-transection analog) and identify which local feedback rules - load gating, perineal-like excitability boosts - suffice to recover stepping.
- **Wosnitza 2013** — [Inter-leg coordination in the control of walking speed in Drosophila](https://doi.org/10.1242/jeb.078139)  
  - animals: Cat, Insects  
  - The transition from one gait to another is discontinuous and it can be shown that quadrupeds select the energetically optimal gait at a given speed (Hoyt and Taylor, 1981).
- **Frigon et al. 2013** — [Split-belt walking alters the relationship between locomotor phases and cycle duration across speeds in intact](https://doi.org/10.1523/jneurosci.3931-12.2013)  
  - animals: Cat  
  - Split-belt walking alters the phase-cycle duration relationship similarly in intact and chronic spinalized cats; the spinal cord alone performs substantial speed-dependent processing of afferent timing information.
- **Krouchev and Drew 2013** — [Motor cortical regulation of sparse synergies provides a framework for the flexible control of precision walki](https://doi.org/10.3389/fncom.2013.00083)  
  - animals: Cat  
  - Motor cortex regulates sparse synergies during precision walking; cortical control selects and modulates a small set of muscle synergies rather than individual muscles, a framework for goal-dependent gait modification.
- **McKay 2012** — [Optimization of Muscle Activity for Task-Level Goals Predicts Complex Changes in Limb Forces across Biomechani](https://doi.org/10.1371/journal.pcbi.1002465)  
  - animals: Cat, Human  
  - Computational model study. In an unrestrained balance task in cats, we demonstrate that achieving task-level constraints center of mass forces and moments while minimizing control effort predicts detailed patterns of muscle activity and ground reaction forces in an anatomically-realistic musculoskeletal model.
- **Hatz 2012** — [Control of ankle extensor muscle activity in walking cats](https://doi.org/10.1152/jn.00944.2011)  
  - animals: Cat, Mammals · pathways: Ib excitatory; type II excitatory  
  - In conscious cats with an isolated medial gastrocnemius, stance-phase and slope-dependent modulation of ankle-extensor activity is well described by constant central drive with constant proprioceptive gains: Ib force feedback is the primary modulator, group II adds a small tonic contribution, and Ia feedback contributes none — terrain compensation without changing central commands or nervous-system gains.
- **Markin 2012** — [Motoneuronal and muscle synergies involved in cat hindlimb control during fictive and real locomotion: a compa](https://doi.org/10.1152/jn.00865.2011)  
  - animals: Cat  
  - Markin and Rybak model. Cat EMG data. Focus on synergies.
- **Harrison 2012** — [Forelimb muscle activity during equine locomotion](https://doi.org/10.1242/jeb.065441)  
  - animals: Cat, Human  
  - Computational model study. The patterns of net muscular torques developed over one gait cycle have been calculated for a wide variety of species using measurements of joint kinematics and ground reaction forces (Clayton et al., 2000; Colborne et al., 1997; Dogan et al., 1991; Fowler et al., 1993;...
- **Gossard 2011** — [Chapter 2--the spinal generation of phases and cycle duration.](https://doi.org/10.1016/b978-0-444-53825-3.00007-3)  
  - animals: Cat  
  - In decerebrate cats, fictive locomotion without sensory feedback shows a built-in asymmetry - cycle period varies predominantly with extensor phase duration - and ankle dorsiflexion greatly prolongs the extension phase during fictive locomotion but not fictive scratching, evidence that locomotor and scratch rhythms rely on distinct spinal interneuronal components.  
  - *Robot/sim:* Build a rhythm generator whose period is dominated by extensor-phase duration with separate locomotor and scratch components, plus ankle-dorsiflexion afferent input that prolongs extension only in locomotor mode; test mode-specific ablation.
- **Frigon 2010** — [Effects of ankle and hip muscle afferent inputs on rhythm generation during fictive locomotion.](https://doi.org/10.1152/jn.01028.2009)  
  - animals: Cat · afferents: Ia, Ib, II  
  - During spontaneous fictive locomotion in decerebrate cats, group I afferents from the plantaris (ankle extensor) nerve always reset the locomotor rhythm, prolonging extension when stimulated in extension and terminating flexion to initiate extension when stimulated in flexion, whereas hip muscle nerve stimulation (rectus femoris, caudal gluteal, sartorius) reset the rhythm only in restricted epochs, with sartorius, particularly at group II strength, shortening the flexion phase. The epoch-specific access of hip afferents implies the rhythm generator operates with several subdivisions to determ  
  - *Robot/sim:* Implement a rhythm generator with phase subdivisions such that extensor group I input prolongs extension and terminates flexion globally while hip afferent input acts only within specific phase windows; ablating the subdivision should make hip afferents reset the rhythm globally.
- **Brownstone and Bui 2010** — [Spinal interneurons providing input to the final common path during locomotion](https://doi.org/10.1016/b978-0-444-53613-6.00006-x)  
  - animals: Mammals, Cat  
  - Spinal interneurons providing input to the final common path during locomotion; framework relating interneuron diversity to the assembly of locomotor commands.
- **Jankowska and Edgley 2010** — [Functional subdivision of feline spinal interneurons in reflex pathways from group Ib and II muscle afferents;](https://doi.org/10.1111/j.1460-9568.2010.07354.x)  
  - animals: Cat, Mammals · afferents: Ib, II  
  - Reassesses the subdivision of adult mammalian spinal intermediate-zone interneurons and finds no compelling reason to separate those with group Ib input from those co-excited by group II afferents: distributed input patterns, projection targets (ipsilateral, contralateral, bilateral), excitatory versus inhibitory identity, and task-dependent reflex changes during locomotion are all consistent with shared premotor interneurons integrating group I and II signals - prompting the proposed renaming to group I/II interneurons.  
  - *Robot/sim:* Replace separate Ib and II reflex channels in a locomotion model with a single premotor interneuron population receiving convergent group I and II input, and test whether the merged population reproduces the task-dependent reflex modulation and reversals seen during locomotion.
- **Donelan 2009** — [Force Regulation of Ankle Extensor Muscle Activity in Freely Walking Cats](https://doi.org/10.1152/jn.90918.2008)  
  - animals: Cat · afferents: Ib, Ia, II · pathways: Ib excitatory  
  - In freely walking cats with the medial gastrocnemius functionally isolated, changes in ankle extensor activity across walking conditions correlate strongly with Ib (force) afferent activity rather than with Ia or group II activity, and force feedback contributes about 30% of total muscle activity during level walking. The pathway gain from muscle force to motoneuron depolarization is length-independent while loop gain rises with length through the intrinsic force-length property, and force feedback's contribution increases during upslope and decreases during downslope walking, providing a simp  
  - *Robot/sim:* Implement autogenic Ib force feedback on the ankle extensor with length-dependent loop gain (via the force-length property) supplying roughly 30% of stance activation; ablate it during slope changes to test terrain compensation.
- **Karayannidou 2009** — [Maintenance of Lateral Stability During Standing and Walking in the Cat](https://doi.org/10.1152/jn.90934.2008)  
  - animals: Cat  
  - Lateral perturbations in the cat recruit context-specific postural mechanisms: during standing, corrective redistribution of muscle activity between symmetrical limbs (contralateral hindlimb extensor excitation with ipsilateral inhibition); during walking, a corrective lateral step whose direction depends on push timing relative to the step cycle and whose EMG signature is marked hip abductor and adductor reprogramming — reconfiguring the base of support rather than redistributing activity.  
  - *Robot/sim:* Implement phase-dependent lateral-perturbation correction in a legged model — standing mode redistributes extensor activity across limbs, walking mode triggers a corrective lateral step via step-cycle-gated hip abductor/adductor modulation — and test recovery against push timing.
- **Chiel 2009** — [The Brain in Its Body: Motor Control and Sensing in a Biomechanical Context](https://doi.org/10.1523/jneurosci.3338-09.2009)  
  - animals: Cat, Human, Insects, Lamprey, Rat, Salamander · pathways: Biomechanically mediated preflexive feedback  
  - Reviews molluscan feeding, postural control in cats and humans, locomotion simulations in lamprey, insect, cat and salamander, and rat vibrissal sensing to argue that adaptive behavior emerges from nervous-system-body-environment interaction: control is shared between nervous system and periphery, neural activity organizes degrees of freedom into biomechanically meaningful subsets, mechanics alone can play crucial roles in enforcing gait patterns, and the mechanics of sensors is crucial for their function.  
  - *Robot/sim:* Embed morphologically realistic muscle and sensor mechanics so that body dynamics contribute to gait enforcement (preflexes); progressively remove neural correction loops and quantify the locomotor stability retained by mechanics alone.
- **Pérez 2009** — [An Intersegmental Neuronal Architecture for Spinal Wave Propagation under Deletions](https://doi.org/10.1523/jneurosci.1737-09.2009)  
  - animals: Cat  
  - During cat scratching, traveling electrical waves along the spinal cord persist through extensor-burst deletions that leave the cycle unaltered, whereas deletions that perturb the cycle coincide with loss of the traveling wave; numerical simulations of an asymmetric two-layer CPG distributed longitudinally along the cord reproduce both the sinusoidal waves and both deletion classes, supporting a longitudinal chain organization of CPG networks.  
  - *Robot/sim:* Implement a longitudinally distributed asymmetric two-layer CPG and perturb individual unit oscillators to reproduce the deletion dichotomy (cycle-preserving versus cycle-perturbing deletions); use it to calibrate intersegmental coupling in multi-segment walkers.
- **Ross and Nichols 2009** — [Heterogenic Feedback Between Hindlimb Extensors in the Spontaneously Locomoting Premammillary Cat](https://doi.org/10.1152/jn.90338.2008)  
  - animals: Cat · pathways: Ib disynaptic inhibition  
  - During spontaneous treadmill stepping in the premammillary decerebrate cat, force-dependent heterogenic inhibition between hindlimb extensors persists (quadriceps onto gastrocnemius, gastrocnemius onto plantaris/FHL) but distal-onto-proximal inhibition is weaker than during the crossed-extension reflex, yielding a proximal-to-distal gradient of Ib inhibition that supports interjoint coordination and limb stability.
- **Frigon and Gossard 2009** — [Asymmetric control of cycle period by the spinal locomotor rhythm generator in the adult cat.](https://doi.org/10.1113/jphysiol.2009.176669)  
  - animals: Cat · pathways: Fictive locomotion without sensory feedback  
  - In decerebrate and spinal cats, fictive locomotor cycle period varied predominantly with the extension phase while the flexion phase stayed relatively invariant; this walking-speed asymmetry (stance varies, swing constant) persisted without phasic sensory feedback, supraspinal structures, pharmacology, or sustained stimulation, establishing it as an intrinsic property of spinal CPG organization.  
  - *Robot/sim:* Build a half-center CPG whose cycle period is modulated primarily through extension (stance) duration with swing duration near-invariant; verify the asymmetry survives removal of all simulated sensory feedback, then add feedback to test modulation.
- **Nichols and Ross 2009** — [The Implications of Force Feedback for the λ Model](https://doi.org/10.1007/978-0-387-77064-2_36)  
  - animals: Cat, Mammals · pathways: Ib excitatory; Ib inhibition  
  - In the λ-model framework, autogenic length feedback compensates muscle nonlinearities while positive force feedback — during level cat stepping largely restricted to gastrocnemius — reinforces the stiff ankle–knee linkage; heterogenic inhibitory force feedback spans different joints and axes, so it coordinates interjoint action and shifts activation thresholds, making threshold a feedback-dependent quantity rather than a pure descending control variable.
- **Barriere 2008** — [Prominent Role of the Spinal Central Pattern Generator in the Recovery of Locomotion after Partial Spinal Cord](https://doi.org/10.1523/jneurosci.5692-07.2008)  
  - animals: Cat  
  - A dual-lesion paradigm in cats shows the spinal CPG undergoes plastic changes after partial spinal cord injury that are shaped by locomotor training: cats trained on the treadmill after a partial thoracic lesion expressed bilateral hindlimb locomotion within hours of a subsequent complete transection, whereas untrained cats walked asymmetrically with the limb on the partially lesioned side recovering first. Re-expression of the hindlimb locomotor pattern after partial SCI therefore arises mostly from intrinsic changes below the lesion in the CPG and afferent inputs rather than from remnant des  
  - *Robot/sim:* Model a spinal CPG whose synapses adapt during repetitive rhythmic 'training' activation, then remove descending drive completely; the trained network should re-express bilateral rhythm rapidly while the untrained network recovers slowly and asymmetrically, reproducing the lesion hierarchy.
- **McCrea and Rybak 2008** — [Organization of mammalian locomotor rhythm and pattern generation](https://doi.org/10.1016/j.brainresrev.2007.08.006)  
  - animals: Cat, Mammals  
  - Central pattern generators (CPGs) located in the spinal cord produce the coordinated activation of flexor and extensor motoneurons during locomotion. Previously proposed architectures for the spinal locomotor CPG have included the classical half-center oscillator and the unit burst generator (UBG) comprised of multiple coupled oscillators. We have recently proposed another organization in which a two-level CPG has a common rhythm generator (RG) that controls the operation of the pattern formation (PF) circuitry responsible for motoneuron activation. These architectures are discussed in relatio
- **Pearson 2008** — [Role of sensory feedback in the control of stance duration in walking cats](https://doi.org/10.1016/j.brainresrev.2007.06.014)  
  - animals: Cat · pathways: Ib stance to swing; Ia or II stance to swing  
  - Review Paper
- **Scrivens 2008** — [A robotic device for understanding neuromechanical interactions during standing balance control.](https://doi.org/10.1088/1748-3182/3/2/026002)  
  - animals: Cat, Human  
  - Computational model study. Here we demonstrate that independent variations in either stance width or delayed neural feedback gains can have profound and often surprisingly detrimental effects on the postural stability of the system.
- **Jankowska 2008** — [Spinal interneuronal networks in the cat: elementary components.](https://doi.org/10.1016/j.brainresrev.2007.06.022)  
  - animals: Cat  
  - Establishes that feline commissural interneuron networks operate as elementary building blocks incorporated into larger networks: they receive mono- and disynaptic reticulospinal and vestibulospinal input plus mono- to oligosynaptic muscle-afferent input, and reach motoneurons mono- or disynaptically, so one population can serve reflex, postural, locomotor, and voluntary coordination.  
  - *Robot/sim:* Model left-right coordination as commissural interneuron blocks receiving mono-/disynaptic descending and muscle-afferent input; ablate individual blocks to predict which interlimb coordination components fail during walking.
- **Deliagina 2008** — [Spinal and supraspinal postural networks](https://doi.org/10.1016/j.brainresrev.2007.06.017)  
  - animals: Cat, Human, Lamprey  
  - In the lamprey, the postural control system is driven by vestibular input.
- **Welch 2008** — [A Feedback Model Reproduces Muscle Activity During Human Postural Responses to Support-Surface Translations](https://doi.org/10.1152/jn.01110.2007.)  
  - animals: Cat, Human  
  - Computational model study. We investigated whether a simple feedback law could explain temporal patterns of muscle activation in response to support-surface translations in human subjects.
- **McKay 2008** — [Functional muscle synergies constrain force production during postural tasks.](https://doi.org/10.1016/j.jbiomech.2007.09.012)  
  - animals: Cat  
  - Computational model study. We recently demonstrated that a set of ﬁve functional muscle synergies were sufﬁcient to characterize both hindlimb muscle activity and active forces during automatic postural responses in cats standing at multiple postural conﬁgurations.
- **Maas 2007** — [The effects of self-reinnervation of cat medial and lateral gastrocnemius muscles on hindlimb kinematics in sl](https://doi.org/10.1007/s00221-007-0938-8)  
  - animals: Cat · afferents: Ia, II, Cutaneous  
  - Self-reinnervation of cat gastrocnemius muscles (motor function recovers, proprioceptive feedback permanently absent) leaves level- and upslope-walking kinematics recovered within 14-19 weeks but produces permanent ankle and interjoint deficits in downslope walking - indicating MG/LG proprioceptive feedback is specifically required for regulating ankle extensors on downslope, with compensation by other sensory sources (e.g., cutaneous) or altered central drive elsewhere.  
  - *Robot/sim:* Ablate length and force feedback from ankle extensors in a slope-walking simulator: level and upslope gaits should tolerate the loss via central-drive/cutaneous compensation while downslope develops persistent ankle yield - a test of feedback redundancy.
- **Guevremont 2007** — [Physiologically based controller for generating overground locomotion using functional electrical stimulation.](https://doi.org/10.1152/jn.01177.2006)  
  - animals: Cat  
  - In spinal cats stepping via functional electrical stimulation, an intrinsically timed controller achieves overground stepping more easily (lower sensitivity to initial stimulation parameters) but cannot compensate walkway resistance or muscle fatigue, while a sensory-driven controller that switches between unloaded flexion and loaded extension phases adapts to loading — motivating a combined controller that runs on intrinsic timing but resets phase from sensory signals.  
  - *Robot/sim:* Implement a fixed-timing CPG with load-based sensory phase resetting (combined controller) in a neuromuscular biped; ablating the sensory reset should cause failures under fatigue-like actuator force decay and changing walkway resistance, reproducing the timed controller's limitation.
- **Procházka et al. 2007** — [Predictive and reactive tuning of the locomotor CPG](https://doi.org/10.1093/icb/icm065)  
  - animals: Cat  
  - The locomotor CPG is tuned both predictively (state-dependent adjustment to expected limb mechanics) and reactively (perturbation-driven corrections); proprioceptive feedback does not merely correct errors but continuously calibrates the pattern.
- **Lockhart 2007** — [Optimal sensorimotor transformations for balance](https://doi.org/10.1038/nn1986)  
  - animals: Cat, Human  
  - Computational model study. Optimal sensorimotor transformations for balance Daniel B Lockhart 1 & Lena H Ting 2 Here we have identiﬁed a sensorimotor transformation that is used by a mammalian nervous system to produce a multijoint motor behavior.
- **McCrea and Rybak 2007** — [Modeling the mammalian locomotor CPG: insights from mistakes and perturbations](https://doi.org/10.1016/S0079-6123(06)65015-2)  
  - animals: Cat  
  - A computational model of the mammalian spinal cord circuitry incorporating a two-level central pattern generator (CPG) with separate half-center rhythm generator (RG) and pattern formation (PF) networks.
The model consists of interacting populations of interneurons and motoneurons described in the Hodgkin-Huxley style. Locomotor rhythm generation is based on a combination of intrinsic (persistent sodium current dependent) properties of excitatory RG neurons and reciprocal inhibition between the two half-centers comprising the RG. The two-level architecture of the CPG was suggested from an anal
- **Torres-oviedo 2007** — [Muscle Synergies Characterizing Human Postural Responses](https://doi.org/10.1152/jn.01360.2006.)  
  - animals: Cat, Human  
  - These results suggest that muscle synergies represent a general neural strategy underlying muscle coordination in postural tasks.
- **Hultborn and Nielsen 2007** — [Spinal control of locomotion--from cat to man.](https://doi.org/10.1111/j.1748-1716.2006.01651.x)  
  - animals: Cat, Human, Vertebrates  
  - Review establishing that spinal networks generate the basic locomotor rhythm across vertebrates including man, with limb sensory feedback essential for effective locomotion: sensory regulation reaches motoneurons via reflex pathways that bypass the rhythm generators and also acts on the locomotor networks themselves, controlling phase timing, shaping muscle-activity patterns, adding excitatory drive, and driving long-term adaptation - the basis for treadmill-training rehabilitation after spinal cord injury.  
  - *Robot/sim:* Implement both sensory routes in a locomotion model - direct reflex pathways to motoneurons plus afferent input to the rhythm generator - and ablate each separately to test their predicted contributions to phase timing, pattern shaping, and net excitatory drive.
- **Rossignol 2006** — [Plasticity of connections underlying locomotor recovery after central and/or peripheral lesions in the adult m](https://doi.org/10.1098/rstb.2006.1889)  
  - animals: Cat, Mice, Rat, Human, Mammals  
  - Review concluding that locomotor recovery after spinal lesions in adult mammals is partly due to plasticity within existing spinal locomotor networks: locomotor training changes the excitability of simple reflex pathways and more complex circuitry, adaptation to lesions entails changes at both spinal and supraspinal levels, and the cat-derived framework extends to rat, mouse, and human spinal pattern generation.  
  - *Robot/sim:* Model lesion and recovery by removing pathways and re-tuning residual connection weights through training-like afferent input, testing which plasticity rules restore stepping when descending drive is lost.
- **Perreault and Raastad 2006** — [Contribution of morphology and membrane resistance to integration of fast synaptic signals in two thalamic cel](https://doi.org/10.1113/jphysiol.2006.113043)  
  - animals: Cat  
  - Morphology and membrane resistance determine integration of fast synaptic signals in two thalamic cell types; dendritic structure sets the window for synaptic summation.
- **Mileusnic 2006** — [Mathematical Models of Proprioceptors. I. Control and Transduction in the Muscle Spindle](https://doi.org/10.1152/jn.00868.2005)  
  - animals: Cat  
  - Computational model study. In the case of simultaneous static and dynamic fusimotor efferent stimulation, we demonstrated the importance of including the experimentally observed effect of partial occlusion.
- **Gregor 2006** — [Mechanics of slope walking in the cat: quantification of muscle load, length change, and ankle extensor EMG pa](https://doi.org/10.1152/jn.01300.2004)  
  - animals: Cat · afferents: Cutaneous  
  - During slope walking in the cat, stance-phase mechanics invert between conditions: downslope walking increases ankle and knee extensor muscle-tendon stretch and peak stretch velocities while decreasing paw-pad forces and ankle and hip extensor moments, with the opposite upslope. Because these mechanical variables are what muscle length, muscle force, and paw-pad cutaneous feedback encode, slope-dependent EMG synergy changes plausibly arise from altered motion-dependent afferent input rather than CPG reorganization.  
  - *Robot/sim:* Drive a walking controller's extensor burst regulation with simulated length, force, and paw-pad feedback and replay the measured slope mechanics; test whether the same CPG produces the observed slope-dependent EMG changes without retuning.
- **Rybak 2006** — [Modelling spinal circuitry involved in locomotor pattern generation: insights from the effects of afferent sti](https://doi.org/10.1113/jphysiol.2006.118711)  
  - animals: Cat · pathways: Ia reciprocal inhibition  
  - Computational model of 2 level CPG
- **Frigon and Rossignol 2006** — [Functional plasticity following spinal cord lesions.](https://doi.org/10.1016/s0079-6123(06)57016-5)  
  - animals: Cat  
  - After spinal cord injury, reflex pathways caudal to the lesion are initially depressed by low motoneuron excitability and then recover — sometimes to exaggeration (spasticity) — and in spinal cats step training normalizes transmission in simple reflex pathways, suggesting that the modified afferent inflow must itself be normalized for a stable locomotor rhythm to be re-expressed.
- **Rybak 2006** — [Modelling spinal circuitry involved in locomotor pattern generation: insights from deletions during fictive lo](https://doi.org/10.1113/jphysiol.2006.118703)  
  - animals: Cat  
  - Rybak and McCrea. Mammalian CPG model. Fictive locomotor activity in decerebrated cats.
- **Hultborn 2006** — [Spinal reflexes, mechanisms and concepts: from Eccles to Lundberg and beyond.](https://doi.org/10.1016/j.pneurobio.2006.04.001)  
  - animals: Cat · afferents: Ia, II, Ib, Flexor reflex afferents  
  - Historical review of Eccles' intracellular recordings defining five classical spinal reflex systems - recurrent inhibition of motoneurons via motor axon collaterals and Renshaw cells, pathways from muscle spindles and Golgi tendon organs, presynaptic inhibition, and the flexor reflex - and of Lundberg's demonstration that spinal interneurons converge segmental sensory afferents with descending pathways, the keystone for modern hypotheses of spinal motor control later tested in behaving animals and non-invasively in humans.  
  - *Robot/sim:* Implement the classical circuit library - Renshaw recurrent inhibition, spindle and tendon-organ reflex pathways, presynaptic inhibition, flexor reflex - as separable modules in a spinal network model, ablating each to quantify its contribution to motor output.
- **McVea 2005** — [A Role for Hip Position in Initiating the Swing-to-Stance Transition in Walking Cats](https://doi.org/10.1152/jn.00511.2005)  
  - animals: Cat · pathways: Ia or II swing to stance  
  - This investigation obtained data that support the hypothesis that afferent signals associated with hip flexion play a role in initiating the swing-to-stance transition of the hind legs in walking cats
- **Langlet 2005** — [Mid-Lumbar Segments Are Needed for the Expression of Locomotion in Chronic Spinal Cats](https://doi.org/10.1152/jn.00909.2004)  
  - animals: Cat  
  - In chronically spinalized cats, a second transection at caudal L3 or L4 permanently abolishes treadmill locomotion even after weeks of training, while lesions at L2 or rostral L3 spare it; fast paw shakes and clonidine-induced hindlimb hyperextension persist below the lesion, so mid-lumbar segments are necessary specifically for locomotor rhythm expression, not for other rhythmic motor patterns or caudal motoneuron function.  
  - *Robot/sim:* Model a spatially distributed spinal network where only mid-lumbar-equivalent modules generate the locomotor pattern while caudal segments retain reflex and rhythmic capacity; ablating the mid-lumbar module should eliminate gait but leave paw-shake-like rhythmicity intact.
- **Stecina et al. 2005** — [Parallel reflex pathways from flexor muscle afferents evoking resetting and flexion enhancement during fictive](https://doi.org/10.1113/jphysiol.2005.095505)  
  - animals: Cat · pathways: Ia or II stance to swing  
  - Parallel reflex pathways from flexor muscle afferents evoke resetting and flexion enhancement during fictive locomotion and scratch; group I/II flexor afferents access the rhythm generator through multiple routes.
- **Lafrenière-Roula et al. 2005** — [Deletions of Rhythmic Motoneuron Activity During Fictive Locomotion and Scratch Provide Clues to the Organizat](https://doi.org/10.1152/jn.00216.2005)  
  - animals: Cat · pathways: Fictive locomotion without sensory feedback  
  - Deletions of rhythmic motoneuron activity during fictive locomotion and scratch provide key clues to CPG organization — evidence for the two-level rhythm/pattern architecture.
- **Nielsen et al. 2005** — [Organization of common synaptic drive to motoneurones during fictive locomotion in the spinal cat](https://doi.org/10.1113/jphysiol.2005.091744)  
  - animals: Cat  
  - Common drive to motoneurones during fictive locomotion is organized in low-frequency coherent components that group muscles into functional synergies; locomotor drive is shared rather than muscle-specific.
- **Ekeberg and Pearson 2005** — [Computer simulation of stepping in the hind legs of the cat: an examination of mechanisms regulating the stanc](https://doi.org/10.1152/jn.00065.2005)  
  - animals: Cat · pathways: Ib stance to swing  
  - A three-dimensional simulation of cat hind-leg stepping shows that stance termination governed by ankle extensor force (load) signals — alone or combined with hip position — yields stable stepping and correct interleg timing, whereas the hip position signal alone fails; mutual inhibition between controllers restores stability but not correct timing. Coordination depends critically on load-sensitive signals from each leg, with mechanical linkages mediated by these signals playing a significant role in establishing the alternating gait.  
  - *Robot/sim:* Reproduce directly: gate the stance-to-swing transition on thresholded ankle-extensor force (with and without hip-angle gating) in a hind-leg/biped model, ablate each channel, and verify the load channel is necessary for stable alternation and correct contralateral timing.
- **Angel 2005** — [Candidate interneurones mediating group I disynaptic EPSPs in extensor motoneurones during fictive locomotion ](https://doi.org/10.1113/jphysiol.2004.076034)  
  - animals: Cat · afferents: Ia, Ib · pathways: Ib disynaptic excitation  
  - Identifies a candidate population of excitatory interneurons in intermediate laminae of mid-to-caudal L7 that are activated by extensor group I afferents at monosynaptic latency specifically during the extensor phase of MLR-evoked fictive locomotion, project to extensor motor nuclei, and burst rhythmically in extension even without afferent stimulation — a previously unknown population that could mediate the group I disynaptic excitation of extensor motoneurones while also supplying central extensor drive.  
  - *Robot/sim:* Implement stance-gated group I afferent excitation of extensor motoneurones through a disynaptic interneuron layer whose cells also carry central extensor drive; ablating either the afferent gate or the interneuron population should remove load-dependent extensor reinforcement.
- **Yakovenko 2005** — [Control of locomotor cycle durations.](https://doi.org/10.1152/jn.00991.2004)  
  - animals: Cat  
  - In MLR-evoked fictive locomotion in cats, spontaneous cycle-duration variation is carried more often by flexor than extensor phase durations (22 of 31 experiments), so the CPG is not inherently extensor- or flexor-biased; a half-center model with two parameters per timing element — background drive ('bias') and sensitivity ('gain') — fits the phase-duration plots, background drive determines which half-center is 'dominant', and normal-cat data are fitted equally well, implying sensory input and central drive combine to set locomotor phase durations.  
  - *Robot/sim:* Implement a half-center oscillator with independent bias and gain on each timing element; reproduce the phase-duration scaling and switch the dominant phase by redistributing background drive — a direct mechanism for duty-factor control in CPG walkers.
- **Büschges 2005** — [Sensory control and organization of neural networks mediating coordination of multisegmental organs for locomo](https://doi.org/10.1152/jn.00615.2004)  
  - animals: Stick Insect, Cat, Lamprey  
  - Review comparing stick insect, cat, and lamprey: locomotor patterns arise from central pattern generating networks, local sensory feedback about movements and forces in the locomotor organs, and coordinating signals from neighboring segments or appendages; the central network controlling a multi-segmental organ comprises multiple segmental CPGs matching the organ's structure, and common schemes of sensory feedback operate across walking species.  
  - *Robot/sim:* Structure the controller as one CPG module per joint or segment with local force and movement feedback plus inter-module coordinating connections; ablate the local feedback to test the multi-CPG organization hypothesis.
- **Quevedo 2005** — [Intracellular analysis of reflex pathways underlying the stumbling corrective reaction during fictive locomoti](https://doi.org/10.1152/jn.00176.2005)  
  - animals: Cat, Human · afferents: Cutaneous · pathways: Cutaneous flexor excitation  
  - Intracellular recordings in decerebrate cats during fictive locomotion show the stumbling corrective reaction is built from di- and trisynaptic cutaneous excitation of knee-flexor and ankle-extensor motoneurons, with locomotor-phase-dependent reconfiguration of inhibition (increased onto ankle flexors, suppressed onto extensors) and motoneuron membrane potential sculpting whether short-latency EPSPs actually recruit firing.  
  - *Robot/sim:* Implement phase-gated cutaneous reflex pathways: di-/trisynaptic excitation of knee flexors and ankle extensors with reciprocal inhibitory gating (suppressed onto extensors, increased onto flexors during flexion) plus motoneuron membrane-potential sculpting, and test foot-dorsum contact responses during swing.
- **Quevedo 2005** — [Stumbling corrective reaction during fictive locomotion in the cat.](https://doi.org/10.1152/jn.00175.2005)  
  - animals: Cat · pathways: Cutaneous flexor excitation  
  - Stimulation of the cutaneous superficial peroneal nerve during the flexion phase of fictive locomotion in decerebrate cats reproduces the full stumbling corrective reaction — brief flexor excitation, flexor inhibition, recruitment of knee and ankle extensors, then post-stimulation flexor excitation with a prolonged flexion phase — showing the corrective synergy is assembled by spinal locomotor circuitry from short-latency cutaneous reflexes plus actions of the rhythm-generating network.
- **Donelan and Pearson 2004** — [Contribution of Force Feedback to Ankle Extensor Activity in Decerebrate Walking Cats](https://doi.org/10.1152/jn.00325.2004)  
  - animals: Cat · afferents: Ib, Ia, II · pathways: Ib excitatory  
  - Measured the loop gain of homonymous positive force feedback from ankle extensor Golgi tendon organs in decerebrate walking cats: about 0.2 at short muscle lengths (roughly 20 percent of total activity and force) rising to 0.5 at long lengths (about 50 percent), with the length dependence arising from the intrinsic force-length property while the force-to-motoneuron gain stays length-independent - establishing positive Ib feedback as a substantial contributor to stance extensor drive, not a minor modulator.  
  - *Robot/sim:* Implement homonymous positive force feedback on stance extensors with loop gain 0.2-0.5, scaled by the muscle force-length property; ablating it should cut extensor activity by a length-dependent 20-50 percent.
- **Yamaguchi 2004** — [The central pattern generator for forelimb locomotion in the cat.](https://doi.org/10.1016/s0079-6123(03)43011-2)  
  - animals: Cat  
  - In decerebrate cats with fictive forelimb locomotion evoked by repetitive cervical lateral funiculus stimulation, the shortest descending pathway to cervical motoneurons via the CPG is disynaptic, and the intercalated interneurons are proposed to form mutually exciting, reverberating circuits that constitute the CPG itself. Disynaptic, trisynaptic, and polysynaptic PSPs in motoneurons are phase-related to the rhythm with reciprocal modulation between extensor and flexor responses, indicating interactions between the mediating pathways.  
  - *Robot/sim:* Implement the forelimb CPG as a reverberating, mutually excitatory interneuron loop driven disynaptically from a descending command; ablating the mutual excitation should collapse the rhythm and its phase-locked PSP modulation.
- **Duysens et al. 2004** — [Sensory Influences on Interlimb Coordination During Gait](https://doi.org/10.1007/978-1-4419-9056-3_1)  
  - animals: Cat, Human  
  - Review of sensory influences on interlimb coordination during gait across species and preparations; organizes cutaneous and proprioceptive contributions to interlimb phase coupling.
- **Pearson 2004** — [Generating the walking gait: role of sensory feedback.](https://doi.org/10.1016/s0079-6123(03)43012-4)  
  - animals: Cat  
  - Synthesizes cat walking evidence that feedback from muscle proprioceptors establishes the timing of major phase transitions, contributes to burst production, generates some features of the motor pattern, and is required for adaptive modification after alterations in leg mechanics; argues that afferent signals likely reorganize the functioning of central networks, making 'afferent modulation of a hard-wired CPG' too simplistic a framework.  
  - *Robot/sim:* Model phase-transition timing and burst production as afferent-controlled processes rather than fixed CPG timing; ablate individual proprioceptive channels and quantify the loss of adaptive motor-pattern modification to altered leg mechanics.
- **Donelan and Pearson 2004** — [Contribution of sensory feedback to ongoing ankle extensor activity during the stance phase of walking](https://doi.org/10.1139/y04-043)  
  - animals: Human, Cat · afferents: II, Ib · pathways: Ib excitatory; type II excitatory  
  - Quantitative review establishing that load-related sensory feedback contributes up to 60% of ongoing ankle extensor activity during the stance phase of walking, with secondary spindle endings (human) and Golgi tendon organs (human and cat) the likely receptors — autogenic positive feedback reinforcing stance extensor force. Argues that resolving which receptor groups set extensor magnitude across locomotor tasks requires network simulations coupled to forward-dynamic musculoskeletal models, since experimental strategies alone cannot dissociate the distributed contributions.  
  - *Robot/sim:* Implement autogenic positive force feedback (Ib) and length/velocity feedback (II) onto stance-gated ankle extensor motoneurons; ablate each channel in turn and test whether unloading reduces extensor activity by up to 60% as reported.
- **Ting 2004** — [Ratio of Shear to Load Ground-Reaction Force May Underlie the Directional Tuning of the Automatic Postural Res](https://doi.org/10.1152/jn.00773.2003)  
  - animals: Cat  
  - This study sought to identify the sensory signals that encode perturbation direction rapidly enough to shape the directional tuning of the automatic postural response.
- **Ivashko 2003** — [Modeling the spinal cord neural circuitry controlling cat hindlimb movement during locomotion](https://doi.org/10.1016/S0925-2312(02)00832-9)  
  - animals: Cat · pathways: Ib disynaptic excitation; Ib disynaptic excitation; Ib stance to swing; Ib swing to stance; Ia stance to swing; Ia swing to stance; Ia monosynaptic; Ia monosynaptic excitation; Mechanosensory monosynaptic excitation; Ia disynaptic inhibition  
  - Abstract
A computational model of the spinal cord neural circuitry that controls locomotor movements of simulated cat hindlimbs. The neural circuitry includes two central pattern generators integrated with reflex circuits. All neurons were modeled in the Hodgkin–Huxley style. The musculoskeletal system includes two three-joint hindlimbs and the trunk. Each
hindlimb is actuated by nine one- and two-joint muscles (a Hill-type model). Our simulations allow us to suggest a specific network architecture in the spinal cord and a pattern of feedback connectivities (from Ia and Ib fibers and touch s
- **Geyer 2003** — [Positive force feedback in bouncing gaits?](https://doi.org/10.1098/rspb.2003.2454)  
  - animals: Cat · pathways: Ib stance to swing  
  - Bring more context on how the transition from gaits to running are formed and give models examples of how operations roll out.
- **Poppele and Bosco 2003** — [Sophisticated spinal contributions to motor control](https://doi.org/10.1016/S0166-2236(03)00073-0)  
  - animals: Vertebrates, Frog, Turtle, Cat  
  - Review paper.
There is an extensive propriospinal network of reciprocal excitatory and inhibitory connections active during locomotion that may activate MNs far from the site of pattern generator.

A key issue in motor control is how sensory
inputs direct and inform motor output, – that is, the
sensorimotor process. Other major issues involve the
actual control of the motor apparatus. In general, there
are at least three basic requirements for motor control:
the transformations that map information from sensory
to motor coordinates, the specification of individual
muscle activations to achieve
- **Bouyer and Rossignol 2003** — [Contribution of Cutaneous Inputs From the Hindpaw to the Control of Locomotion. II. Spinal Cats](https://doi.org/10.1152/jn.00497.2003)  
  - animals: Cat, Mammals · pathways: Cutaneous stance modification  
  - Sequential hindpaw cutaneous denervation in spinal cats shows cutaneous inputs are necessary for plantar foot placement and weight bearing during spinal locomotion — fully denervated cats never recovered either despite 35–71 days of treadmill training — while partial denervations reveal substantial adaptive capacity of the spinal cord.
- **Wilmink and Nichols 2003** — [Distribution of Heterogenic Reflexes Among the Quadriceps and Triceps Surae Muscles of the Cat Hind Limb](https://doi.org/10.1152/jn.00833.2002)  
  - animals: Cat  
  - Heterogenic proprioceptive reflexes among cat quadriceps and triceps surae are organized by articulation rather than motor-unit composition: excitatory length feedback strongly links the uniarticular vastus muscles (and vastus to soleus) to regulate joint stiffness, while force-related inhibition is absent between vastus muscles but strong and bidirectional between vastus and rectus femoris and between triceps surae and quadriceps - regulating interjoint coupling and, together with length feedback, the endpoint mechanical properties.  
  - *Robot/sim:* Implement heterogenic feedback by articulation: excitatory length feedback among uniarticular synergists for joint stiffness, bidirectional force inhibition between muscles spanning different joints for interjoint coupling; ablate the interjoint force inhibition to test endpoint property regulation.
- **Gagnon 2003** — [Contribution of Cutaneous Inputs From the Hindpaw to the Control of Locomotion. I. Intact Cats](https://doi.org/10.1152/jn.00496.2003)  
  - animals: Cat · afferents: Cutaneous  
  - Bilateral hindpaw cutaneous denervation at the ankle in intact cats barely disrupts level treadmill walking but severely impairs ladder and incline walking early on; flexor EMG remains permanently elevated and mediolateral ground reaction forces increase by 200 percent, showing cutaneous feedback matters most for demanding locomotor contexts and that active adaptive mechanisms compensate without fully restoring the control pattern. The companion paper shows these same inputs become critical for foot placement after spinalization.  
  - *Robot/sim:* Add foot/ankle cutaneous afferent channels to a walking simulation and ablate them: the model should reproduce near-normal level walking with degraded ladder/incline behavior, persistently elevated flexor drive, and increased mediolateral ground reaction forces.
- **Lam and Pearson 2002** — [The Role of Proprioceptive Feedback in the Regulation and Adaptation of Locomotor Activity](https://doi.org/10.1007/978-1-4615-0713-0_40)  
  - animals: Cat · pathways: Ia or II swing to stance  
  - muscle spindle ends swing
- **Duysens 2002** — [A walking robot called human: lessons to be learned from neural control of locomotion.](https://doi.org/10.1016/s0021-9290(01)00187-7)  
  - animals: Cat, Human  
  - Distills cat and human locomotion into design principles for walking robots: control at three levels (actuator/motoneuron, whole-limb flexion-extension oscillators with mutual inhibition, interlimb coordination), with the most essential feedback at the limb level, where activation of the extensor part of the limb oscillator must be triggered by feedback signalling onset of loading via limb load sensors, and flexor activation must require unloading below a threshold plus a hip position within the normal end-of-stance range.  
  - *Robot/sim:* Implement the three-level architecture with limb-level decision rules: extensor oscillator activation on limb loading (load-sensor threshold) and flexor initiation gated on unloading below threshold together with hip position in the end-of-stance range; ablate each rule and quantify gait robustness loss.
- **Yakovenko 2002** — [Spatiotemporal activation of lumbosacral motoneurons in the locomotor step cycle](https://doi.org/10.1152/jn.00479.2001)  
  - animals: Cat  
  - A 3-D digital reconstruction of 27 cat hindlimb motoneuron pools, modulated by compiled EMG profiles, reveals a rostrocaudal oscillation of spinal activity across the step cycle: the caudal third of the lumbosacral enlargement dominates during stance, activation shifts abruptly rostral during swing, and a transient caudal focus at the stance-swing transition coincides with retractor muscles (gracilis, posterior biceps, semimembranosus, semitendinosus) that clear the foot from the ground.  
  - *Robot/sim:* Distribute a simulated cat-scale MN pool model rostrocaudally per these digitized coordinates and phase-modulate the pools with EMG-derived profiles; test whether the rostrocaudal oscillation and the caudal stance-swing transition focus emerge, and how they degrade when sensory inputs are ablated.
- **Perreault 2002** — [Motoneurons have different membrane resistance during fictive scratching and weight support.](https://doi.org/10.1523/jneurosci.22-18-08259.2002)  
  - animals: Cat  
  - Motoneuron membrane resistance differs between fictive scratching and weight support; motor state changes motoneuron conductance and therefore the gain of synaptic inputs — the same afferent signal has state-dependent effects.
- **Burke 2001** — [Patterns of locomotor drive to motoneurons and last-order interneurons: clues to the structure of the CPG.](https://doi.org/10.1152/jn.2001.86.1.447)  
  - animals: Cat · afferents: Cutaneous  
  - During fictive locomotion in decerebrate cats, low-threshold cutaneous reflex pathways to the flexor digitorum longus motor pool are differentially controlled at the level of last-order excitatory interneurons, and the modulation pattern is preserved whether FDL fires in early flexion or during extension. The motor pool and reflex-modulation patterns reveal distinct early versus late flexion components, and pattern formation can be separated from rhythm generation — evidence the CPG embodies partially distinct neural organizations for these two functions.  
  - *Robot/sim:* Implement a CPG with separable rhythm-generator and pattern-formation modules plus phase-gated cutaneous reflex gains following the early- versus late-flexion modulation; ablate the pattern module to test rhythm persistence with altered motoneuron scheduling.
- **Lam and Pearson 2001** — [Proprioceptive modulation of hip flexor activity during the swing phase of locomotion in decerebrate cats.](https://doi.org/10.1152/jn.2001.86.3.1321)  
  - animals: Cat · pathways: Ia or II swing to stance  
  - muscle spindle ends swing
- **Quevedo 2000** — [Group I disynaptic excitation of cat hindlimb flexor and bifunctional motoneurones during fictive locomotion.](https://doi.org/10.1111/j.1469-7793.2000.t01-1-00549.x)  
  - animals: Cat · afferents: Ia, Ib · pathways: Ib disynaptic excitation  
  - During fictive locomotion in decerebrate cats, group I stimulation evokes locomotor-dependent disynaptic excitation (single interneurone, mean latency 1.64 ms) of most flexor (89 percent) and many bifunctional (64 percent) motoneurones, largest during flexion and from homonymous nerves, evoked by both tendon-organ (Ib) and muscle-spindle (Ia) afferents and separate from the extensor group I pathway - indicating distinct flexor and extensor excitatory interneurone groups that reinforce ongoing locomotor activity throughout the limb.  
  - *Robot/sim:* Add disynaptic group I excitation of flexor and bifunctional motoneurone pools, phase-weighted to peak during flexion and driven by homonymous afferents, as load-reinforcement of ongoing activity; ablate it to quantify the loss of flexor reinforcement.
- **Pang and Yang 2000** — [The initiation of the swing phase in human infant stepping: importance of hip position and leg loading](https://doi.org/10.1111/j.1469-7793.2000.00389.x)  
  - animals: Human, Cat  
  - In supported stepping of human infants, hip flexion combined with high limb load prolongs stance and delays swing, whereas hip extension with low load shortens stance and advances swing, remarkably similar to reduced cat preparations. Hip position and load show an inverse relationship at the time of swing initiation, indicating the two sensory factors combine to regulate the stance-to-swing transition, consistent with similar brainstem and spinal walking circuitry in infants and cats.  
  - *Robot/sim:* Implement swing-initiation gating as a combined hip-extension and low-load condition with an inverse interaction between hip angle and limb load; ablate either input to reproduce prolonged stance and delayed swing.
- **Marcoux and Rossignol 2000** — [Initiating or Blocking Locomotion in Spinal Cats by Applying Noradrenergic Drugs to Restricted Lumbar Spinal S](https://doi.org/10.1523/jneurosci.20-22-08577.2000)  
  - animals: Cat  
  - Focal alpha2-noradrenergic action at restricted lumbar levels both starts and stops spinal walking: topical clonidine over L3-L4 or L5-L7, or microinjection at L3-L5, induces treadmill locomotion in acutely spinalized cats, whereas yohimbine at L3-L5 (not L6) or a transection at L3-L4 blocks it — mid-lumbar segments are a sufficient pharmacological trigger and a necessary substrate for the rhythm.  
  - *Robot/sim:* Build a spatially segmented spinal CPG in which only mid-lumbar-equivalent modules can be tonically activated to launch the full hindlimb pattern; ablating those modules while leaving caudal motoneuron pools functional should abolish gait but spare other rhythmic outputs.
- **Gosgnach 2000** — [Depression of group Ia monosynaptic EPSPs in cat hindlimb motoneurones during fictive locomotion.](https://doi.org/10.1111/j.1469-7793.2000.00639.x)  
  - animals: Cat · pathways: Ia presynaptic inhibition  
  - During brainstem-evoked fictive locomotion in decerebrate cats, monosynaptic group Ia EPSPs are tonically depressed to about two-thirds of control in most motoneurones, with group I field potentials similarly reduced and only weak correlation with decreased motoneurone input resistance — evidence that presynaptic inhibition of Ia transmission underlies the tonic depression of stretch reflexes during locomotion.
- **Abelew 2000** — [Local Loss of Proprioception Results in Disruption of Interjoint Coordination During Locomotion in the Cat](https://doi.org/10.1152/jn.2000.84.5.2709)  
  - animals: Cat, Human  
  - These results indicate an important role for the stretch reﬂex and stiffness regulation during locomotion.
- **Dietz and Duysens 2000** — [Significance of load receptor input during locomotion: a review](https://doi.org/10.1016/s0966-6362(99)00052-1)  
  - animals: Human, Cat · pathways: Ib stance to swing  
  - Review making the case that extensor load receptors are central to locomotor control: in the cat, Golgi tendon organ input switches function during walking from Ib inhibition to extensor facilitation proportional to load, and in humans leg extensor activation during stance scales with body weight — one load-regulatory mechanism spanning quadrupedal and bipedal gait.
- **Barbeau 1999** — [Tapping into spinal circuits to restore motor function.](https://doi.org/10.1016/s0165-0173(99)00008-9)  
  - animals: Cat, Human, Vertebrates  
  - Multidisciplinary review arguing that electrically activating spinal interneuronal circuits (reflex and pattern-generating) could coordinate many muscles at once for neuroprostheses, far more tractable than stimulating each muscle individually; draws on phase-adaptable cat hindlimb reflexes, chick embryo rhythmogenic network development, the cat spinal locomotor pattern generator, and locomotor training in incomplete SCI patients to identify candidate circuits and the control problems engineers must solve.  
  - *Robot/sim:* Model the control layer for FNS-style stimulation as a spinal circuit hierarchy (phase-modulated reflexes plus CPG) driving a musculoskeletal limb, instead of per-muscle open-loop stimulation.
- **Prochazka 1999** — Quantifying proprioception  
  - animals: Cat, Human, Mammals  
  - Quantifies proprioceptive afferent discharge during natural movement: spindle and GTO firing ranges, gain changes with movement, and the case that proprioceptive feedback operates with state-dependent gain.
- **Perreault 1999** — [Depression of muscle and cutaneous afferent-evoked monosynaptic field potentials during fictive locomotion in ](https://doi.org/10.1111/j.1469-7793.1999.00691.x)  
  - animals: Cat · afferents: Ia, Ib, II, Cutaneous · pathways: Ia presynaptic inhibition  
  - During MLR-evoked fictive locomotion in decerebrate cats, monosynaptic field potentials evoked by group I, group II, and cutaneous afferents are tonically depressed (intermediate-lamina group II potentials down to a mean of 49% of control), with smaller phase-dependent cyclic modulation superimposed. The depression begins with tonic MLR drive before rhythm onset, indicating that the locomotor state reduces synaptic transmission from primary afferents onto first-order spinal interneurons — afferent influx to the locomotor circuitry is gated down at its earliest synapse during locomotion.  
  - *Robot/sim:* Implement locomotor-state gating on afferent input channels: scale group I/II/cutaneous afferent synaptic gains by a tonic depression factor (~0.8 for dorsal/group I, ~0.5 for intermediate group II) plus a small phase-dependent component, and verify reflex amplitudes are suppressed during locomotion while phase-cycling as observed.
- **Zehr and Stein 1999** — [What functions do reflexes serve during human locomotion](https://doi.org/10.1016/s0301-0082(98)00081-1)  
  - animals: Lamprey, Cat, Human · afferents: Cutaneous · pathways: Cutaneous flexor excitation; Ia monosynaptic excitation  
  - Review establishing that reflexes during locomotion are task-, phase-, and context-dependent, with a division of labor: cutaneous reflexes alter swing-limb trajectory to avoid stumbling, stretch reflexes stabilize limb trajectory and assist force production during stance, and load receptor reflexes support body weight and influence step-cycle timing - functions dynamically reassigned across the cycle and clinically exploitable after neurotrauma.  
  - *Robot/sim:* Implement phase- and task-gated cutaneous (swing trajectory), Ia stretch (stance stabilization and force), and load-receptor (weight support, cycle timing) pathways in a walker; ablate each to reproduce the predicted functional losses.
- **Mori 1999** — [Stimulation of a Restricted Region in the Midline Cerebellar White Matter Evokes Coordinated Quadrupedal Locom](https://doi.org/10.1152/jn.1999.82.1.290)  
  - animals: Cat  
  - Identifies a cerebellar locomotor region in the midline white matter (hook bundle of Russell, crossed fastigiofugal fibers) whose low-threshold microstimulation evokes well-coordinated quadrupedal locomotion in decerebrate cats; it summates with subthreshold mesencephalic locomotor region stimulation and remains effective after MLR lesions, establishing the fastigial nucleus as an independent supraspinal site that triggers both brainstem and spinal locomotor subprograms.  
  - *Robot/sim:* Implement two independent supraspinal drive sites (fastigial/CLR-analog and MLR-analog) converging on the spinal CPG with subthreshold summation; ablate each to test whether the other alone sustains locomotion, and test cycle-time shortening under combined suprathreshold drive.
- **Hiebert and Pearson 1999** — [Contribution of Sensory Feedback to the Generation of Extensor Activity During Walking in the Decerebrate Cat](https://doi.org/10.1152/jn.1999.81.2.758)  
  - animals: Cat  
  - Letting the foot step into a hole reduces knee and ankle extensor activity to roughly 70% of normal in decerebrate cats, and selectively resisting ankle extension restores it; with segmental dorsal-root sections mapped accordingly, more than half of stance-phase extensor drive is proprioceptive — a continuous load-scaled contribution that automatically matches burst intensity to external demands.  
  - *Robot/sim:* Implement autogenic load-dependent extensor burst reinforcement during stance in a neuromechanical walker; removing ground support should cut extensor activity to about 70%, and resisting extension should restore it, directly testing the afferent share of stance drive.
- **Schomburg 1998** — [Flexor reflex afferents reset the step cycle during fictive locomotion in the cat](https://doi.org/10.1007/s002210050522)  
  - animals: Cat · afferents: Flexor reflex afferents, II, III/IV, Cutaneous  
  - Brief trains to flexor reflex afferents (FRA: joint, cutaneous, group II-III muscle afferents) reset the L-dopa-induced fictive locomotor rhythm in high-spinal cats - interrupting extensor activity and initiating flexion when delivered in extension, or prolonging flexion in late flexion - demonstrating that FRA interneurones are constituent elements of the rhythm-generating network organized along half-centre lines.  
  - *Robot/sim:* Implement FRA input as a phase-dependent reset of the half-centre oscillator (extensor interruption plus flexion trigger, flexion prolongation late in flexion) and verify that brief afferent trains shorten or lengthen single cycles after which the original rhythm period resumes.
- **Van de Crommert 1998** — [Neural control of locomotion: sensory control of the central pattern generator and its relation to treadmill t](https://doi.org/10.1016/s0966-6362(98)00010-1)  
  - animals: Cat, Human  
  - Review synthesizing that locomotor recovery through treadmill training in spinalized cats and spinal-cord-injured patients rests on a spinal central pattern generator whose activation and regulation depend on locomotor-related afferent input: adequate sensory input during training can drive and calibrate the spinal locomotor circuitry, a principle for designing gait rehabilitation programs.  
  - *Robot/sim:* Implement a CPG controller whose activation and phase regulation are gated by locomotor-related afferent input, and ablate that input to show that training-like recovery of rhythmic output fails under tonic drive alone.
- **Carlson-Kuhta 1998** — [Forms of forward quadrupedal locomotion. II. A comparison of posture, hindlimb kinematics, and motor patterns ](https://doi.org/10.1152/jn.1998.79.4.1687)  
  - animals: Cat  
  - Per the abstract provided: across steep downslope grades in walking cats, stance-phase yield increases at the ankle, swing knee and ankle flexion decreases with ankle extension lowering the paw to contact, hip extensors fall silent while hip flexors brake the rate of hip extension, and ankle extensor bursts truncate around paw contact - the motor pattern reorganizes to counteract external braking forces rather than simply rescaling level gait.  
  - *Robot/sim:* Gate hip extensor drive off and use hip flexors as stance brakes on downslopes, with ankle extensor bursts centered on paw contact, testing whether a fixed CPG with slope-dependent gating reproduces the graded downslope kinematics.
- **Duysens and Van de Crommert 1998** — [Neural control of locomotion; Part 1: The central pattern generator from cats to humans](https://doi.org/10.1016/s0966-6362(97)00042-8)  
  - animals: Cat, Human  
  - Traces the spinal CPG concept from cats with complete spinal cord transection that recover locomotor function to evidence for a human spinal locomotor CPG drawn from treadmill-training recovery in incomplete spinal cord injury, framing locomotor rehabilitation as exploitation of spinal rhythm-generating circuitry conserved from cats to humans.  
  - *Robot/sim:* Implement a spinal CPG block whose basic rhythmic activation survives removal of step-specific sensory feedback, mirroring transection-cat recovery, and study training-like adaptation regimes of the circuit in simulation.
- **Procházka and Gorassini 1998** — [Ensemble firing of muscle afferents recorded during normal locomotion in cats](https://doi.org/10.1111/j.1469-7793.1998.293bu.x)  
  - animals: Cat · afferents: Ia, II, Ib  
  - Chronic recordings during normal cat stepping show that ensemble spindle afferent firing is largely predictable from muscle length and velocity alone (with significant EMG-linked fusimotor action in triceps surae), and that the ensemble of triceps surae tendon organ afferents accurately encodes whole-muscle force. To a first approximation, large muscle afferents in the cat hindlimb signal muscle velocity (Ia), length (II), and force (Ib) at locomotor speeds and amplitudes.  
  - *Robot/sim:* Instantiate afferent feedback as ensemble codes (Ia proportional to muscle velocity, II to length, Ib ensemble averaging to whole-muscle force) rather than single-unit transfer functions in a locomotion controller, and verify feedback control survives realistic firing variance and fusimotor modulation.
- **Kakuda 1998** — [Dynamic response of human muscle spindle afferents to stretch during voluntary contraction](https://doi.org/10.1111/j.1469-7793.1998.621bb.x)  
  - animals: Cat, Human  
  - This study was conducted to investigate the basic pattern of dynamic and static fusimotor actions on human muscle spindles.
- **Procházka and Gorassini 1998** — [Models of ensemble firing of muscle spindle afferents recorded during normal locomotion in cats](https://doi.org/10.1111/j.1469-7793.1998.277bu.x)  
  - animals: Cat · afferents: Ia  
  - Mathematical models of muscle spindle primary (presumed Ia) firing, validated against chronically recorded hamstring afferents over 132 cat step cycles, are dominated by muscle velocity — firing rate scales approximately with the square root of velocity (power-law exponents 0.5-0.6) — with only small EMG-linked fusimotor components; inverse models recover muscle length accurately from firing profiles, implying the CNS could decode muscle length simply from Ia ensemble firing.  
  - *Robot/sim:* Use a square-root-of-velocity-dominant spindle encoder with a small EMG-linked fusimotor term as the Ia block in reflex models, and run it in reverse to estimate muscle length from afferent firing for state estimation in a neurally controlled walker.
- **Rossignol 1998** — [Pharmacological Activation and Modulation of the Central Pattern Generator for Locomotion in the Cat](https://doi.org/10.1111/j.1749-6632.1998.tb09061.x)  
  - animals: Cat  
  - Box 6128, Station Centre-Ville, Montreal (Qc) H3C 3J7, Canada ABSTRACT: Pharmacological agents have been shown to be capable of inducing a pat- tern of rhythmic activity recorded in muscle nerves or motoneurons of paralyzed spinal cats that closely resembles the locomotor pattern seen in intact...
- **Stephens and Yang 1996** — [Short latency, non-reciprocal group I inhibition is reduced during the stance phase of walking in humans.](https://doi.org/10.1016/s0006-8993(96)00977-8)  
  - animals: Human, Cat · afferents: Ib · pathways: Ib disynaptic inhibition  
  - In intact humans, short-latency non-reciprocal group I inhibition from the medial gastrocnemius nerve onto the conditioned soleus H-reflex, disynaptic and present at rest in most subjects, is significantly reduced during treadmill walking, with some subjects showing significant excitation, mirroring the reduction of resting Ib inhibition toward excitation seen in walking cats. The reduction may be partially accounted for by activation of the triceps surae itself.  
  - *Robot/sim:* Implement disynaptic Ib inhibition onto ankle extensors with a locomotor-phase-dependent gain reduction (partial reversal toward excitation during gait); ablate the phase switch to test its effect on stance extensor activity.
- **Angel 1996** — [Group I extensor afferents evoke disynaptic EPSPs in cat hindlimb extensor motorneurones during fictive locomo](https://doi.org/10.1113/jphysiol.1996.sp021538)  
  - animals: Cat · afferents: Ia, Ib · pathways: Ib disynaptic excitation  
  - During the extension (stance-equivalent) phase of MLR-evoked fictive locomotion in decerebrate cats, group I afferents from ankle extensors evoke disynaptic EPSPs (mean central latency 1.55 ms, single interposed interneuron) in hip, knee, and ankle extensor motoneurons — also elicited by selective group Ia activation via muscle stretch — with hip extensor motoneurons receiving excitation from both homonymous and ankle extensor nerves. The effect arises from cyclic disinhibition of excitatory interneurons rather than motoneuronal voltage-dependent conductances or flexion-phase presynaptic inhib  
  - *Robot/sim:* Implement a disynaptic excitatory pathway from extensor group I afferents to extensor motoneurons, gated ON during extension via interneuron disinhibition (OFF during flexion); ablate it and test for loss of stance force reinforcement.
- **Hiebert 1996** — [Contribution of Hind Limb Flexor Muscle Afferents to the Timing of Phase Transitions in the Cat Step Cycle](https://doi.org/10.1152/jn.1996.75.3.1126)  
  - animals: Cat · afferents: Ia, II · pathways: Ia or II stance to swing  
  - In walking decerebrate cats, flexor muscle lengthening during stance terminates extensor activity and resets the rhythm to ipsilateral flexion with contralateral extension: EDL and iliopsoas act through group Ia spindle afferents and tibialis anterior through group II afferents, providing a flexor-length gate on the stance-to-swing transition.  
  - *Robot/sim:* Implement flexor-length feedback (Ia from iliopsoas and EDL, II from tibialis anterior) that inhibits the extensor half-center to terminate stance; ablation should delay swing onset and prolong stance.
- **Whelan 1996** — [CONTROL OF LOCOMOTION IN THE DECEREBRATE CAT](https://doi.org/10.1016/0301-0082(96)00028-7)  
  - animals: Cat · pathways: Ib stance to swing; Ia or II stance to swing  
  - In the decerebrate cat, locomotion is initiated by the mesencephalic locomotor region acting through the medial medullary reticular formation and the ventrolateral funiculus, and the phase transitions are set by identified afferents: group I Golgi tendon organ afferents prolong stance, while length- and velocity-sensitive afferents from extensor muscles signal leg extension to permit swing — unloaded and extended is the swing-permitting state.
- **McCrea 1995** — [Disynaptic group I excitation of synergist ankle extensor motoneurones during fictive locomotion in the cat.](https://doi.org/10.1113/jphysiol.1995.sp020897)  
  - animals: Cat · afferents: Ia, Ib · pathways: Ib disynaptic inhibition; Ib disynaptic excitation  
  - In decerebrate cats, plantaris nerve group I stimulation evokes short-latency disynaptic IPSPs in medial gastrocnemius motoneurons at rest but EPSPs during the extensor phase of MLR-evoked fictive locomotion, also elicited by selective group Ia activation via Achilles tendon stretch, indicating an excitatory group Ia and Ib feedback system that reinforces ongoing extensor activity during stance. In clonidine-treated spinal preparations the resting inhibition disappears without appearing as excitation, suggesting locomotion inhibits the inhibitory interneurons operating at rest.  
  - *Robot/sim:* Implement phase-switched group I pathways onto synergist ankle extensor motoneurons, disynaptic inhibition at rest and disynaptic excitation during the extensor phase; ablate the switch to test reinforcement of extensor activity during stance.
- **Pearson 1995** — [Proprioceptive regulation of locomotion.](https://doi.org/10.1016/0959-4388(95)80107-3)  
  - animals: Cat, Insects, Arthropods · pathways: Ib stance to swing; Ia or II stance to swing  
  - Review of proprioceptive regulation across walking systems: during locomotion, Golgi tendon organ feedback from extensor muscles reverses from inhibition to excitation, maintaining stance while the extensors are loaded and helping time the stance-to-swing transition, while primary and secondary spindle afferents also influence the timing of the rhythm — establishing phase-dependent reflex reversal as a general feature of cat and arthropod walking.
- **Guertin 1995** — [Ankle extensor group I afferents excite extensors throughout the hindlimb during fictive locomotion in the cat](https://doi.org/10.1113/jphysiol.1995.sp020871)  
  - animals: Cat · pathways: Type 1 swing to stance; type 1 stance to swing  
  - 1. The effects of stimulating hindlimb extensor nerves (100-200 ms trains, 100 Hz, < or = 2 times threshold) during the flexor and extensor phases of the locomotor step cycle were analysed in the decerebrate, paralysed cat during fictive locomotion evoked by stimulation of the mesencephalic locomotor region. 2. Stimulation during extension of either the medial gastrocnemius (MG), lateral gastrocnemius-soleus (LGS) or plantaris (Pl) nerves was equally effective in increasing the duration and amplitude of electroneurogram (ENG) activity recorded in ipsilateral ankle, knee and hip extensor nerves
- **Whelan 1995** — [Stimulation of the group I extensor afferents prolongs the stance phase in walking cats.](https://doi.org/10.1007/bf00241961)  
  - animals: Cat · afferents: Ia, Ib  
  - Group I stimulation of extensor nerves (LG-Sol, plantaris, VL/VI) in walking decerebrate cats prolongs stance and delays flexor burst onset, with spatial summation across ankle and knee extensors and oligosynaptic excitation of heteronymous extensor EMG; stimulation early in flexion abruptly terminates swing and reinstates stance - supporting late-stance hindlimb unloading as a necessary condition for swing initiation.  
  - *Robot/sim:* Implement the stance-to-swing gate as extensor group I feedback (with spatial summation across ankle and knee extensors) that prevents flexor burst onset until unloading; simulate late-stance unloading to trigger swing and premature unloading to precipitate it.
- **Perreault 1995** — Effects of stimulation of hindlimb flexor group II afferents during fictive locomotion in the cat  
  - animals: Cat · pathways: Ia or II stance to swing; type II excitatory  
  - Excitatory input to flexor type II afferent.
- **Kriellaars 1994** — [Mechanical entrainment of fictive locomotion in the decerebrate cat.](https://doi.org/10.1152/jn.1994.71.6.2074)  
  - animals: Cat · pathways: Ia or II stance to swing  
  - Sinusoidal hip movements as small as 5-20 degrees entrain fictive locomotion in decerebrate cats via low-threshold stretch-sensitive afferents from intrinsic hip muscles (joint capsular afferents unnecessary), with extensor bursts locked to imposed hip flexion as the extensors are stretched and loaded - a positive-feedback mechanism by which the rhythm generator prolongs extension. Subharmonic entrainment and frequency-dependent phase shifts mark the locomotor rhythm generator as a nonlinear oscillator whose rhythm-generating interneurons receive highly convergent afferent input from muscles s  
  - *Robot/sim:* Implement hip stretch afferents as direct inputs to the extensor half-center of the rhythm generator that prolong stance and entrain the cycle to imposed hip motion; ablating them should abolish mechanical entrainment and stance prolongation.
- **Gossard 1994** — [Transmission in a locomotor-related group Ib pathway from hindlimb extensor muscles in the cat](https://doi.org/10.1007/bf00228410)  
  - animals: Cat · pathways: Ib stance to swing  
  - Phasic stimulations of group I afferents from ankle and knee extensor muscles may entrain the intrinsic locomotor rhythm.
 
The intrinsic locomotor rhythm then acts on motoneurons through the spinal rhythm generators , and concluded that the major part of these effects come from Golgi tendon organ Ib afferents

Injection of nialamide and L-DOPA evoked long lasting reflexes upon stimulations of high threshold afferents before spontaneous fictive locomotion commenced.

Interneurons must therefore be located in the L7 and S1 spinal segments.
- **Perreault 1994** — [Microstimulation of the medullary reticular formation during fictive locomotion.](https://doi.org/10.1152/jn.1994.71.1.229)  
  - animals: Cat  
  - Microstimulation of the medullary reticular formation during fictive locomotion evokes site-dependent modulation of extensor activity — reticulospinal effects are state-dependent.
- **Beloozerova and Sirota 1993** — [The role of the motor cortex in the control of accuracy of locomotor movements in the cat.](https://doi.org/10.1113/jphysiol.1993.sp019498)  
  - animals: Cat  
  - Motor cortex controls the accuracy of locomotor movements: cortical neurons encode precise foot placement, and cortical inactivation impairs accuracy without abolishing walking.
- **Perreault 1993** — [Activity of medullary reticulospinal neurons during fictive locomotion.](https://doi.org/10.1152/jn.1993.69.6.2232)  
  - animals: Cat · pathways: Fictive locomotion without sensory feedback  
  - In paralyzed, high-decerebrate cats, most medullary reticulospinal neurons projecting to the lumbar cord are phasically modulated during fictive locomotion - ENG-related, rhythm-related, or locomotion-onset tonic subpopulations, intermingled, with forelimb-related cells dorsal to hindlimb-related ones - demonstrating that reticulospinal phase modulation persists without phasic peripheral afferent input and matches the complexity seen in intact animals.  
  - *Robot/sim:* Add a descending reticulospinal layer whose units are phasically modulated by the CPG even when afferent feedback is disabled; ablate it to test its contribution to interlimb coordination in the biped model.
- **Buford and Smith 1993** — [Adaptive control for backward quadrupedal walking. III. Stumbling corrective reactions and cutaneous reflex se](https://doi.org/10.1152/jn.1993.70.3.1102)  
  - animals: Cat · pathways: Cutaneous flexor excitation  
  - In cats walking forward or backward on a treadmill, paw taps that obstruct swing evoke direction-appropriate stumbling corrective reactions — draw the limb away from the obstacle, lift it over, then reposition for stance — while stance-phase taps leave the ongoing step largely unchanged but exaggerate the next swing; mechanical taps produce far richer responses than electrical pulses, so cutaneous reflex sensitivity is structured by phase, gait direction, and stimulus modality.
- **Pearson and Collins 1993** — [Reversal of the influence of group Ib afferents from plantaris on activity in medial gastrocnemius muscle duri](https://doi.org/10.1152/jn.1993.70.3.1009)  
  - animals: Cat · afferents: Ib, Ia · pathways: Ib excitatory; Ib inhibition  
  - In clonidine-treated acute and chronic spinal cats, group I stimulation of the plantaris nerve entrained the locomotor rhythm and, during locomotor activity, group Ib afferents from plantaris exerted an excitatory action on medial gastrocnemius bursts (30-50 ms latency) through the extensor half-center of the rhythm generator — reversing the inhibitory effect the same stimuli produced on tonic activity without locomotion; selective Ia activation by muscle vibration neither entrained the rhythm nor augmented the bursts.  
  - *Robot/sim:* Implement state-dependent Ib feedback in a half-center CPG: autogenic inhibition at rest switching to excitation of the extensor half-center during locomotion, with rhythm entrainment by group I stimulation but not by Ia-specific input; ablate the reversal to test loss of load-dependent stance burst augmentation.
- **Pearson 1992** — [Entrainment of the locomotor rhythm by group Ib afferents from ankle extensor muscles in spinal cats.](https://doi.org/10.1007/bf00230939)  
  - animals: Cat · pathways: Ib stance to swing  
  - Rhythmic contractions of ankle extensor muscles or stimulation of their group I afferents entrain the spinal locomotor rhythm in spinal cats, with flexor bursts timed to follow release of extensor stretch by about 200 ms — direct evidence that Ib input during stance inhibits flexor-burst generation and promotes extensor activity, so a decline in Ib activity near the end of stance helps time the stance-to-swing transition.
- **LaBella 1992** — [Low-threshold, short-latency cutaneous reflexes during fictive locomotion in the "semi-chronic" spinal cat.](https://doi.org/10.1007/bf00231657)  
  - animals: Cat · afferents: Cutaneous · pathways: Cutaneous stance modification; Cutaneous flexor excitation  
  - During fictive locomotion in semi-chronic spinal cats, low-threshold cutaneous nerve stimulation evokes predominantly short-latency (5-15 ms) excitatory reflexes in extensor and flexor motor nerves that are phase-modulated to peak during each nerve's active period, with occasional longer-latency responses modulated out of phase; the stereotypic modulation across muscles and animals implies the locomotor CPG exerts stereotypic control over cutaneous reflex interneurons.  
  - *Robot/sim:* Add phase-gated cutaneous reflex pathways to motoneuron pools with gain peaking during each pool's own burst, plus an out-of-phase longer-latency component, and verify the stereotypic phase-dependent modulation the CPG exerts over cutaneous reflex transmission.
- **Pearson and Rossignol 1991** — [Fictive motor patterns in chronic spinal cats](https://doi.org/10.1152/jn.1991.66.6.1874)  
  - animals: Cat · afferents: Cutaneous · pathways: Fictive locomotion without sensory feedback  
  - Chronic spinal cats generate three distinctly different fictive motor patterns, locomotion, paw shake (approximately 8 Hz rhythmic activity), and slow rhythmic leg flexion, depending only on the site and mode of stimulation (perineal skin, water jet to the paw, paw squeeze), with clonidine facilitating evocation. Passive leg movement from extension to flexion progressively shortened flexor bursts, lengthened the cycle period, and made the pattern harder to evoke, and trained late-spinal animals expressed more complex, position-dependent patterns than untrained early-spinal ones, demonstrating   
  - *Robot/sim:* Build a pattern generator supporting multiple coexistent fictive modes (stepping, paw shake, rhythmic flexion) switched by input site, and modulate flexor burst termination by a position-related afferent signal to reproduce the effects of passively moving the legs from extension to flexion.
- **Pratt 1991** — [Functionally complex muscles of the cat hindlimb. IV. Intramuscular distribution of movement command signals a](https://doi.org/10.1007/bf00229407)  
  - animals: Cat · afferents: Cutaneous · pathways: Cutaneous flexor excitation  
  - In awake walking cats, four of five broad bifunctional thigh muscles contain intramuscular subregions that are differentially recruited across locomotion, scratching, and paw shaking, and low-threshold cutaneous afferents subdivide the same muscles, but with muscle-subregion-specific and often complex responses that do not simply follow the background locomotor synergies. Cutaneous reflex premotoneuronal circuitry is therefore not synonymous with the locomotor CPG even though the two systems interact powerfully, implying specialized convergence onto task-dependent muscle subunits.  
  - *Robot/sim:* Model bifunctional muscles as multiple subregion pools with task-dependent recruitment and subregion-specific cutaneous reflex gains; removing the reflex subdivision tests whether locomotor synergies alone reproduce the observed activation patterns.
- **Pratt 1991** — [Functionally complex muscles of the cat hindlimb. I. Patterns of activation across sartorius.](https://doi.org/10.1007/bf00229406)  
  - animals: Cat  
  - During sixteen horizontal platform-translation directions in the standing cat, the three neuromuscular compartments of biceps femoris are activated non-homogeneously, with petal-shaped EMG tuning curves rotating progressively counterclockwise from anterior (hip-extensor-like) through middle to posterior (knee-flexor-like) compartments; the active region moves continuously across anterior and middle compartments with a discontinuity at the middle-posterior border, and functional units do not match the anatomical compartments exactly.  
  - *Robot/sim:* Represent multi-compartment muscles with direction-tuned regional activations (rotated tuning curves per compartment) instead of one activation per muscle, and test improvements in directional postural-response control.
- **Drew 1991** — [Functional organization within the medullary reticular formation of the intact unanesthetized cat. III. Micros](https://doi.org/10.1152/jn.1991.66.3.919)  
  - animals: Cat  
  - Microstimulation of the cat medullary reticular formation during locomotion evokes phase-dependent responses incorporated into the ongoing pattern — flexor excitation in each limb during that muscle's natural activity period, inhibition of ipsilateral extensors, mixed responses in contralateral vastus lateralis — with response latencies shortest during a muscle's active period, showing that reticulospinal transmission is gated by the locomotor network itself.  
  - *Robot/sim:* Implement descending reticulospinal-like drive onto a CPG with phase-dependent pathway gating (gain high during each muscle's active phase) and reproduce the shortening of response latency with activity state as a test of state-dependent transmission.
- **Orsal et al. 1990** — [Interlimb coordination during fictive locomotion in the thalamic cat.](https://doi.org/10.1007/bf00228795)  
  - animals: Cat  
  - Interlimb coordination during fictive locomotion in the thalamic cat persists without phasic sensory feedback, demonstrating central interlimb coupling that sensory input can then modulate.
- **Gossard 1990** — [Phase-dependent modulation of primary afferent depolarization in single cutaneous primary afferents evoked by ](https://doi.org/10.1016/0006-8993(90)90334-8)  
  - animals: Cat · afferents: Cutaneous  
  - Cutaneous primary afferent terminals in the cat hindlimb are rhythmically depolarized during fictive locomotion (L-PAD), and primary afferent depolarization evoked by peripheral nerve stimulation is phase-modulated by the CPG - minimal during flexion, maximal during extension - so presynaptic inhibitory gating of cutaneous reflex transmission is locomotor-phase dependent. Centrally generated and peripherally evoked presynaptic mechanisms are in part separately controlled, since evoked PAD peaked in extension regardless of L-PAD amplitude.  
  - *Robot/sim:* Implement locomotor-phase-dependent presynaptic gain on cutaneous feedback channels (lowest during flexion, highest during extension) in a reflex-gated walker, then ablate the gating to quantify its contribution to reflex timing and gait stability.
- **Buford and Smith 1990** — [Adaptive control for backward quadrupedal walking. II. Hindlimb muscle synergies](https://doi.org/10.1152/jn.1990.64.3.756)  
  - animals: Cat  
  - Backward walking in cats reuses the forward reciprocal synergies — flexors in swing, extensors in stance — but retunes temporal parameters and amplitudes (gastrocnemius active from midswing with ramping stance EMG, near-absent knee and ankle yield, early-ending iliopsoas with prolonged semitendinosus), indicating one direction-general circuit with task-specific reweighting.  
  - *Robot/sim:* Retune only burst timing and amplitude parameters of a single CPG-synergy controller to reproduce backward-walking EMG benchmarks (LG from midswing, ramping stance EMG, absent yield), testing whether one synergy set covers both directions with temporal reweighting alone.
- **Shefchyk 1990** — [Activity of interneurons within the L4 spinal segment of the cat during brainstem-evoked fictive locomotion](https://doi.org/10.1007/bf00228156)  
  - animals: Cat · afferents: II  
  - During MLR-evoked fictive locomotion in the decerebrate cat, midlumbar (L4) interneurons receiving group II input from quadriceps, sartorius, and pretibial flexors — and projecting to L7 motor nuclei — fired predominantly with ipsilateral flexor bursts and were less responsive to peripheral input during extension: a phase-gated group II pathway positioned to shape flexor activity.  
  - *Robot/sim:* Implement a midlumbar interneuron layer receiving group II afferent input from quadriceps, sartorius, and pretibial flexor pathways and projecting to motor nuclei, with gain reduced during extension; test its contribution to flexor burst generation in a CPG model.
- **Gossard 1989** — [Intra-axonal recordings of cutaneous primary afferents during fictive locomotion in the cat](https://doi.org/10.1152/jn.1989.62.5.1177)  
  - animals: Cat · afferents: Cutaneous  
  - Intra-axonal recordings during fictive locomotion in decorticated cats show that all cutaneous primary afferents develop locomotor-rhythmic membrane potential fluctuations with two depolarization waves per cycle (maximal during flexor activity in the majority of units), demonstrating that the CPG phasically controls the efficacy of transmission in cutaneous pathways at a presynaptic level as part of the locomotor program.  
  - *Robot/sim:* Implement phase-dependent presynaptic gating of cutaneous feedback channels: scale cutaneous synaptic gains by a locomotor-phase gate with maximal presynaptic inhibition phased to flexor activity, and compare reflex behavior against an ungated model.
- **Barbeau and Rossignol 1987** — [Recovery of locomotion after chronic spinalization in the adult cat](https://doi.org/10.1016/0006-8993(87)91442-9)  
  - animals: Cat  
  - Adult cats spinalized at T13 recover plantar weight-bearing treadmill walking after weeks to months of interactive training, with cycle duration, stance/swing timing, joint coupling, and hindlimb EMG resembling intact cats at comparable speeds — establishing the chronic spinal cat as the standard model for testing how training, drugs, and other treatments influence locomotor recovery below the lesion.
- **Pratt and Jordan 1987** — [Ia inhibitory interneurons and Renshaw cells as contributors to the spinal mechanisms of fictive locomotion.](https://doi.org/10.1152/jn.1987.57.1.56)  
  - animals: Cat · afferents: Ia · pathways: Ia reciprocal inhibition  
  - During MLR-evoked fictive locomotion in decerebrate cats, Renshaw cells discharge rhythmically in phase with the motoneuron pools that excite them, with extensor Renshaw activity peaking at the end of extension coincident with extensor motoneuron hyperpolarization; Renshaw inhibition of quadriceps Ia inhibitory interneurons is not sufficient to disrupt the phasing of reciprocal inhibition, though extensor Renshaw cells may contribute to the late-extension decline in IaIN discharge. Ia reciprocal IPSPs in motoneurons vary significantly in amplitude across the fictive step cycle.  
  - *Robot/sim:* Implement Renshaw cells phase-locked to motoneuron pools with recurrent inhibition onto Ia inhibitory interneurons; ablate the Renshaw-to-IaIN connection to test the late-extension decline of reciprocal inhibition.
- **Noga 1987** — [The role of Renshaw cells in locomotion: antagonism of their excitation from motor axon collaterals with intra](https://doi.org/10.1007/bf00236206)  
  - animals: Cat  
  - Blocking nicotinic transmission from motor axon collaterals with mecamylamine during fictive locomotion abolished recurrent IPSPs and rhythmic Renshaw firing yet left locomotor drive potentials, burst timing, and Ia inhibitory interneuron rhythmicity intact, with motoneuron and Ia-inhibitory-interneuron bursts carrying more spikes at higher frequency - establishing that Renshaw cells are not part of the locomotor CPG and do not time bursts, but limit firing rates within them.  
  - *Robot/sim:* Implement Renshaw cells as recurrent negative feedback acting only on motoneuron and Ia-inhibitory-interneuron firing rates; ablating them in simulation should raise within-burst firing while leaving burst timing and the locomotor rhythm unchanged.
- **Mori 1987** — [Integration of posture and locomotion in acute decerebrate cats and in awake, freely moving cats](https://doi.org/10.1016/0301-0082(87)90010-4)  
  - animals: Cat  
  - Integration of posture and locomotion in decerebrate cats: postural set and locomotion share circuitry, with tonically controlled parallel descending systems.
- **Edgley and Jankowska 1987** — [An interneuronal relay for group I and II muscle afferents in the midlumbar segments of the cat spinal cord.](https://doi.org/10.1113/jphysiol.1987.sp016676)  
  - animals: Cat · afferents: Ia, Ib, II, Cutaneous · pathways: type II excitatory; II inhibitory  
  - Midlumbar (L4) interneurons with dominant group II input - also driven by group I, cutaneous, joint, and interosseous afferents - project to hindlimb motor nuclei and are activated monosynaptically by group I and II afferents, making them a disynaptic relay that adds 0.4-2.5 ms latency to length- and force-afferent actions on motoneurons relative to caudal lumbar pathways.  
  - *Robot/sim:* Insert a midlumbar relay interneuron layer (adding 0.4-2.5 ms latency) between group I/II afferents and motoneuron pools; ablate it to test the functional cost of losing the longer-latency afferent actions during perturbed stepping.
- **Conway 1987** — [Proprioceptive input resets central locomotor rhythm in the spinal cat.](https://doi.org/10.1007/bf00249807)  
  - animals: Cat · pathways: Ib stance to swing  
  - Extensor group I afferents have access to central rhythm generators
 >Important in reflex regulation of stepping.


Rythm generators came mainly from Golgi Tendon Organ Ib afferents for group I

Increased load of limb extensor during the stance phase enhance and prolong extensor activity while simultaneously delaying the transition to the swing phase of the step cycle.
- **Armstrong 1986** — [Supraspinal contributions to the initiation and control of locomotion in the cat.](https://doi.org/10.1016/0301-0082(86)90021-3)  
  - animals: Cat
- **Lovely 1986** — [Effects of training on the recovery of full-weight-bearing stepping in the adult spinal cat.](https://doi.org/10.1016/0014-4886(86)90094-4)  
  - animals: Cat  
  - Most adult spinal cats (14 of 16) recover full-weight-bearing treadmill stepping within a month of transection when cued by tail pinch or crimping, and daily treadmill training emphasizing complete weight bearing raises maximum treadmill speed about 2.6-fold over untrained cats (0.619 versus 0.240 m/s) - task-specific training markedly improves spinal stepping capacity below the lesion.  
  - *Robot/sim:* Model a below-lesion spinal stepping network with training-dependent plasticity driven by weight-bearing activity; disabling the adaptation should reproduce the untrained-cat speed plateau.
- **Abraham 1985** — [The distal hindlimb musculature of the cat. Cutaneous reflexes during locomotion.](https://doi.org/10.1007/bf00235875)  
  - animals: Cat · afferents: Cutaneous · pathways: Cutaneous flexor excitation  
  - Brief electrical stimulation of cutaneous nerves during cat treadmill locomotion evoked phase-gated short-latency reflexes — inhibition of extensors and excitation of flexors, expressed only during locomotor phases in which the respective motoneuron pools were active — while longer-latency components were gated similarly but not identically, indicating partial independence from the basic locomotor drive on the pool.  
  - *Robot/sim:* Implement phase-gated cutaneous reflex pathways — flexor excitation and extensor inhibition enabled only during each pool's active locomotor phase, with separately gated late components — and test perturbation responses across the step cycle.
- **Fleshman 1984** — [Peripheral and Central Control of Flexor Digitorum Longus and Flexor Hallucis Longus Motoneurons: The Synaptic](https://doi.org/10.1007/bf00235825)  
  - animals: Cat · afferents: Ia, Cutaneous · pathways: Ia monosynaptic excitation; Fictive locomotion without sensory feedback  
  - Cat FDL and FHL have identical mechanical actions yet are used differently in locomotion, and the differentiation is central: during fictive locomotion FHL motoneurons are coactive with ankle extensors in extension while FDL fires a brief burst at flexion onset; both pools share monosynaptic Ia excitation from both muscles, and only some cutaneous circuits differ (disynaptic EPSPs from saphenous and superficial peroneal nerves in FDL only; forelimb cutaneous IPSPs chiefly in FHL) - reflex organization alone cannot explain the functional divergence.  
  - *Robot/sim:* Drive two synergist pools at one joint with distinct CPG burst targets (extensor-coactive versus flexion-onset) while sharing monosynaptic Ia excitation across both; ablate the differentiated central drive to show common reflex organization cannot produce the functional divergence.
- **Grillner and Zangger 1984** — [The effect of dorsal root transection on the efferent motor pattern in the cat's hindlimb during locomotion.](https://doi.org/10.1111/j.1748-1716.1984.tb07400.x)  
  - animals: Cat · pathways: Fictive locomotion without sensory feedback  
  - After transecting all dorsal roots from one hindlimb, mesencephalic walking cats retain the limb's complex motor pattern, including the double-burst knee flexors and appropriately timed toe dorsiflexor, although with greater variability and occasional breakdown, showing the central network generates the full pattern without phasic afferent input while afferents stabilize and fine-tune it.  
  - *Robot/sim:* Ablate all sensory feedback channels in a CPG walker and test whether the joint-specific burst structure (double-burst flexors) persists with increased variability, a lesion test separating central pattern generation from afferent stabilization.
- **Andersson and Grillner 1983** — [Peripheral control of the cat's step cycle. II. Entrainment of the central pattern generators for locomotion b](https://doi.org/10.1111/j.1748-1716.1983.tb07267.x)  
  - animals: Cat, Mammals  
  - In low-spinal curarized cats with pharmacologically evoked fictive locomotion, sinusoidal hip movements entrain the burst pattern in strict 1:1 coordination over roughly ±5–70% of the resting burst period, with flexion-phase movements driving flexor efferents and extension driving extensors — direct evidence that hip-joint afferent feedback can pace the spinal locomotor CPG.
- **Jankowska and McCrea 1983** — [Shared reflex pathways from Ib tendon organ afferents and Ia muscle spindle afferents in the cat.](https://doi.org/10.1113/jphysiol.1983.sp014663)  
  - animals: Cat · afferents: Ia, Ib  
  - Ia spindle and Ib tendon-organ afferents act through shared premotoneuronal interneurons: costimulation evokes PSPs exceeding the arithmetic sum of the separate effects, in both excitatory and inhibitory pathways at disynaptic or trisynaptic latency - so fusimotor-set spindle input can modulate tendon-organ reflex action, and Ia and Ib afferents operate as one common feedback system rather than parallel channels.  
  - *Robot/sim:* Route Ia and Ib afferents through a shared interneuron layer so spindle input gain-modulates tendon-organ reflex action; compare against independent parallel channels to reproduce the beyond-linear-summation facilitation.
- **Akazawa 1982** — [Modulation of stretch reflexes during locomotion in the mesencephalic cat](https://doi.org/10.1113/jphysiol.1982.sp014319)  
  - animals: Cat · afferents: Ia · pathways: Ia monosynaptic  
  - The soleus stretch reflex in walking decerebrate cats is deeply modulated across the step cycle, peaking at or before the locomotor extensor EMG peak, with the modulation attributable to postsynaptic or presynaptic changes at the alpha-motoneuron rather than fusimotor or afferent-volley changes; reflex gain parallels the cycle rise in intrinsic muscle stiffness, supporting a load-compensation role in stance.  
  - *Robot/sim:* Give a walking model phase-modulated Ia stretch-reflex gain peaking with stance extensor activity and pair it with cycle-varying intrinsic muscle stiffness; ablate the gain modulation (constant gain) to test load-compensation benefit.
- **Andersson and Grillner 1981** — [Peripheral control of the cat's step cycle I. Phase dependent effects of ramp-movements of the hip during "fic](https://doi.org/10.1111/j.1748-1716.1981.tb06867.x)  
  - animals: Cat  
  - Imposed ramp movements of the hip during fictive locomotion in acute spinal curarized cats reset the locomotor rhythm with a phase-dependent sign, prolonging the cycle when applied early (during flexor activity) and markedly shortening it later, with direction-sensitive positive feedback on flexor bursts at the end of flexion and stronger flexor responses at more extended hip positions, demonstrating potent phase-specific peripheral control of the central pattern generator.  
  - *Robot/sim:* Implement a phase-dependent hip-angle reset input to the rhythm generator: early-cycle ramps prolong and late-cycle ramps shorten the cycle, with end-of-flexion directional sensitivity and position-dependent gain; ablating this input should remove the cycle-timing modulation.
- **Jankowska 1981** — [Common interneurones in reflex pathways from group 1a and 1b afferents of ankle extensors in the cat.](https://doi.org/10.1113/jphysiol.1981.sp013556)  
  - animals: Cat · pathways: Ib disynaptic inhibition; Ib disynaptic excitation  
  - Most lamina V-VI interneurones take convergent input from both Ia spindle and Ib tendon-organ afferents of ankle and toe extensors — about 64% shared pathways showing co-excitation, co-inhibition, or opposite actions — so reflex circuitry is not strictly partitioned between length and force feedback, and interneurones projecting to motor nuclei show the same convergence.
- **Perret and Cabelguen 1980** — [Main characteristics of the hindlimb locomotor cycle in the decorticate cat with special reference to bifuncti](https://doi.org/10.1016/0006-8993(80)90207-3)  
  - animals: Cat · afferents: Ia, Flexor reflex afferents  
  - In decorticate cats, all studied hindlimb muscles show alpha-gamma coactivation, pure flexors and extensors alternate simply, and bifunctional pluriarticular muscles receive both flexor and extensor commands, dividing the locomotor cycle into a flexion phase, an extension phase, and two transition phases. The relative weight of these excitations depends on interactions between central commands and peripheral inflow, especially afferents acting according to the flexor reflex pattern, providing a scheme whereby an initially simple rhythmic command plus afferent gating produces the biomechanicall  
  - *Robot/sim:* Implement bifunctional motoneuron pools receiving weighted flexor and extensor CPG commands modulated by flexor-reflex-pattern afferent input; the model should generate transition-phase complexity and task-adapted output from a simple alternating rhythm.
- **Forssberg 1980** — [The locomotion of the low spinal cat. II. Interlimb coordination](https://doi.org/10.1111/j.1748-1716.1980.tb06534.x)  
  - animals: Cat  
  - Chronic spinal kittens on a split-belt treadmill maintain a common hindlimb rhythm across two- to threefold belt speed differences by prolonging the flexion or first extension phase of the fast-belt limb and shortening the flexion phase of the slow-belt limb, with extreme differences producing 2, 3, or even 4 fast-limb steps per slow-limb cycle; support-phase duration mainly follows each limb's own belt speed. All phases of the step cycle are thus modifiable, with several spinal mechanisms coordinating the limbs.  
  - *Robot/sim:* Implement spinal interlimb coordination with per-limb phase modification under asymmetric timing commands (a split-belt analog), including shared rhythm up to two- to threefold speed differences and integer multiple stepping beyond that.
- **Duysens and Pearson 1980** — [Inhibition of flexor burst generation by loading ankle extensor muscles in walking cats](https://doi.org/10.1016/0006-8993(80)90206-1)  
  - animals: Cat · pathways: Ib stance to swing  
  - It has been found that decreased load detected by Golgi organs encourages swing muscles to activate
- **Wand 1980** — [Neuromuscular responses to gait perturbations in freely moving cats.](https://doi.org/10.1007/bf00237937)  
  - animals: Cat · afferents: Cutaneous · pathways: Cutaneous flexor excitation  
  - Impeding the forward swing of a cat's foot evokes rapid elevating reactions whose mechanically elicited EMG sequences are more complex and longer lasting than electrically elicited ones; anesthetizing the foot dorsum abolishes ankle-extensor and knee-flexor responses and reduces ankle-flexor responses to simple stretch reflexes — evidence for cutaneous mediation of the placing-like reaction during stepping, with points of difference from chronic spinal kittens that challenge a purely spinal account.  
  - *Robot/sim:* Implement a swing-phase obstacle-contact reflex: a foot-dorsum cutaneous-analog sensor triggers a phased flexor activation sequence to elevate the limb over the obstacle; ablating the sensor should abolish the elevating strategy, matching the anesthesia result.
- **McCrea 1980** — [Renshaw cell activity and recurrent effects on motoneurons during fictive locomotion.](https://doi.org/10.1152/jn.1980.44.3.475)  
  - animals: Cat  
  - During fictive locomotion in pre- and post-mammillary cats, Renshaw cells fire rhythmic bursts tied to single phases of the step cycle, and recurrent inhibitory IPSPs persist in motoneurons through all phases with amplitudes inversely correlated with motoneuron membrane potential — evidence that recurrent inhibition is phase-structured within the locomotor pattern rather than a static gain setting.
- **Duysens and Loeb 1980** — [Modulation of ipsi- and contralateral reflex responses in unrestrained walking cats](https://doi.org/10.1152/jn.1980.44.5.1024)  
  - animals: Cat · afferents: Cutaneous · pathways: Cutaneous flexor excitation; Cutaneous stance modification  
  - In unrestrained walking cats, cutaneous stimulation evokes phase-modulated responses: two excitatory peaks (about 10 and 25 ms) in flexors that grow toward the end of stance, early inhibition then late excitation (P3, about 35 ms) in extensors when stimuli fall in early stance, and crossed responses in the contralateral limb at 20-25 ms. Concludes the walking cat shows modulation of transmission in a flexor-excitatory/extensor-inhibitory pathway - likely by the flexor part of the spinal locomotor oscillator - rather than strict reflex reversal, plus specialized flexor-inhibitory and extensor-e  
  - *Robot/sim:* Implement cutaneous reflex pathways with locomotor-phase gating: flexor excitation whose gain peaks at late stance, extensor inhibition during early stance, and crossed extensor/flexor responses; ablating the phase gating tests the oscillator-modulation hypothesis versus fixed reflex reversal.
- **Forssberg 1980** — [The locomotion of the low spinal cat. I. Coordination within a hindlimb.](https://doi.org/10.1111/j.1748-1716.1980.tb06533.x)  
  - animals: Cat  
  - The low spinal cat walks with coordinated hindlimb bursts: spinal locomotion is possible below transection and is shaped by afferent input — foundation for spinal locomotion training concepts.
- **Grillner and Zangger 1979** — [On the central generation of locomotion in the low spinal cat](https://doi.org/10.1007/BF00235671)  
  - animals: Cat · pathways: Total afferent inhibition  
  - Cited by Geyer and Herr, 2010.
Total afferent inhibition.

A central network of neurones in the spinal cord has been shown to produce a rhythmic motor output similar to locomotion after suppression of all afferent inflow. The experiments were performed mainly in acute spinal cats (th. 12), which had received DOPA i.v. and the monoamine oxidase inhibitor Nialamide. In some preparations all dorsal roots supplying the spinal cord were transected, in others phasic afferent activity was suppressed by curarization. The activity was recorded as neurograms from nerve filaments or as electromyograms.


- **Prochazka 1979** — [Muscle spindle discharge in normal and obstructed movements.](https://doi.org/10.1113/jphysiol.1979.sp012645)  
  - animals: Cat · afferents: Ia, II  
  - During voluntary movements in the cat, both primary and secondary spindle endings fall silent during active shortening (the faster the shortening, the deeper the silence) and resume firing when shortening is unexpectedly arrested - with increases too gradual to indicate strong alpha-gamma coactivation; above about 0.2 resting lengths per second discharge is dominated by length and velocity variations, with fusimotor action possibly predominating only below that.  
  - *Robot/sim:* Model Ia and II afferents as velocity-dominant encoders that go nearly silent during rapid active shortening and resume discharge on unexpected arrest, with only modest gamma bias - use as the spindle feedback front end in a biped controller.
- **Forssberg 1979** — [Stumbling corrective reaction: a phase-dependent compensatory reaction during locomotion](https://doi.org/10.1152/jn.1979.42.4.936)  
  - animals: Cat · afferents: Cutaneous, Nociceptive · pathways: Cutaneous flexor excitation; Cutaneous stance modification  
  - The classic characterization of the stumbling corrective reaction: dorsum contact during swing evokes short-latency flexor activation lifting the paw over the obstacle via early (~10 ms) and late (~25 ms) pathways to knee flexors, whereas the same stimulus during support inhibits-then-excites extensors and boosts the next swing's flexion, phase-dependent reflex gating organized by the spinal locomotor generator, with painful input instead evoking withdrawal throughout the cycle.  
  - *Robot/sim:* Implement a phase-switched stumble reflex: paw-dorsum tactile input drives short-latency flexor bursts during swing for obstacle clearance and stance-phase extensor modulation with enhanced next-swing flexion; ablate the reflex to quantify its contribution to fall avoidance.
- **Duysens and Stein 1978** — [Reflexes induced by nerve stimulation in walking cats with implanted cuff electrodes.](https://doi.org/10.1007/bf00239728)  
  - animals: Cat · afferents: Cutaneous · pathways: Cutaneous stance modification; Cutaneous flexor excitation  
  - In walking cats, low-threshold cutaneous nerve input powerfully modulates the step cycle: posterior tibial stimulation just before ankle-extensor EMG onset prolongs ipsilateral flexion and contralateral extension, sural and peroneal stimulation activates ankle extensors, and identical gastrocnemius-soleus muscle-nerve stimulation has no effect - the reflexes are cutaneously mediated and phase-dependent.  
  - *Robot/sim:* Implement phase-gated cutaneous reflex channels - pre-extensor-onset flexion prolongation, sural/peroneal-type extensor activation, and stance yielding - in a walking model and ablate them to test perturbation recovery.
- **Grillner and Rossignol 1978** — [On the initiation of the swing phase of locomotion in chronic spinal cats.](https://doi.org/10.1016/0006-8993(78)90973-3)  
  - animals: Cat  
  - In chronic spinal cats walking on a treadmill, a limb held by the paw lifts off and rejoins walking once the hip is brought far enough back, at a hip angle very close to that of normal swing initiation, and lift-off also tends to occur at particular contralateral cycle phases. Establishes hip position and contralateral step-cycle phase as two key factors determining swing initiation, i.e., the stance-to-swing transition is regulated by proprioceptive and interlimb signals rather than a fixed internal clock.  
  - *Robot/sim:* Implement swing initiation in a walking model as a hip-extension position threshold gated by contralateral cycle phase; ablating the hip-position input should delay or prevent swing onset, reproducing the held-paw behavior.
- **Prochazka 1977** — [Ia afferent activity during a variety of voluntary movements in the cat](https://doi.org/10.1113/jphysiol.1977.sp011864)  
  - animals: Cat · afferents: Ia, Cutaneous  
  - Chronic dorsal-root recordings of spindle primary afferents in freely moving cats showed that fusimotor drive to ankle extensors operates mainly in the extension phases of stepping, that unexpected loading (back thrust) evokes a potent fusimotor response, and that during falls Ia discharge decreases despite rising extensor EMG - clear evidence for independent alpha and gamma control - while clonus-like paw-shaking belongs to the normal motor repertoire.  
  - *Robot/sim:* Model fusimotor action as a phase-gated gain on spindle feedback - raised in extension, boosted by load transients, suppressed during falls - and test whether state-dependent afferent gain improves load compensation in a walking controller.
- **Forssberg 1977** — [Phasic gain control of reflexes from the dorsum of the paw during spinal locomotion.](https://doi.org/10.1016/0006-8993(77)90710-7)  
  - animals: Cat · afferents: Cutaneous · pathways: Cutaneous flexor excitation; Cutaneous stance modification  
  - In chronic spinal cats walking on a treadmill, tactile stimulation of the paw dorsum during swing evokes a short-latency flexion response with concomitant crossed extension, while the same stimulus during stance increases ipsilateral extension, a phase-dependent reflex reversal that compensates for unpredicted obstacles by lifting the paw in swing and reinforcing support in stance. The responses are well adapted to ongoing locomotion and leave interlimb coordination intact, except when delivered as the foot approaches the ground after flexion, where the alternating pattern is disturbed.  
  - *Robot/sim:* Implement a phase-switched paw-dorsum cutaneous reflex, flexion plus crossed extension in swing and ipsilateral extension reinforcement in stance, in a walking model; ablating it tests its contribution to obstacle compensation.
- **Duysens 1977** — [Reflex control of locomotion as revealed by stimulation of cutaneous afferents in spontaneously walking premam](https://doi.org/10.1152/jn.1977.40.4.737)  
  - animals: Cat · pathways: Cutaneous stance modification; Cutaneous flexor excitation  
  - Low-intensity stimulation of cutaneous afferents in spontaneously walking premammillary cats prolongs the extensor burst when delivered during stance but shortens swing flexion and advances extensor onset, while high-intensity stimulation prolongs flexion — cutaneous reflex effects on the locomotor pattern are phase dependent, with large fibres inhibiting and small fibres exciting the flexion-generating circuitry.
- **Hulliger 1977** — [Effects of combining static and dynamic fusimotor stimulation on the response of the muscle spindle primary en](https://doi.org/10.1113/jphysiol.1977.sp011840)  
  - animals: Cat · afferents: Ia  
  - In cat soleus muscle spindle primary endings under sinusoidal stretch, combined static and dynamic fusimotor stimulation is dominated by static action at small amplitudes (up to about 50 micrometers) with the dynamic contribution growing progressively with stretch amplitude; at response peaks the two actions sum (dynamic stronger), whereas in troughs static action occludes the weaker dynamic action, and phase differences between conditions remain below 20 degrees. The findings characterize how fusimotor set shapes Ia encoding of stretch.  
  - *Robot/sim:* Implement spindle-primary afferent encoding with separate static and dynamic fusimotor drives reproducing the amplitude-dependent summation and occlusion behavior, and modulate fusimotor set during locomotion to shape Ia feedback.
- **Schomburg 1977** — [Phase-dependent transmission in the excitatory propriospinal reflex pathway from forelimb afferents to lumbar ](https://doi.org/10.1016/0304-3940(77)90187-2)  
  - animals: Cat  
  - During fictive locomotion in high spinal paralyzed cats, forelimb nerve stimulation evokes EPSPs in hindlimb extensor motoneurons only during extension and in flexor motoneurons only during flexion, establishing that transmission in the descending excitatory propriospinal reflex pathway is cyclically gated at the lumbar level by the locomotor rhythm.  
  - *Robot/sim:* Add forelimb-to-hindlimb excitatory propriospinal coupling to a quadruped CPG, gated by lumbar phase (excitation reaching extensor motoneurons during extension and flexor motoneurons during flexion); removing the phase gate should scramble interlimb coordination.
- **Procházka 1976** — [Discharges of single hindlimb afferents in the freely moving cat.](https://doi.org/10.1152/jn.1976.39.5.1090)  
  - animals: Cat · afferents: Ia, Ib  
  - Chronic single-fiber recordings from L7 dorsal roots in unrestrained walking cats: spindle Ia primaries of ankle extensors fire fastest during the phase in which they are passively stretched, with fusimotor drive during active contraction insufficient to fully overcome the unloading effect of rapid shortening, and Ib tendon organs of toe extensors discharge mainly during stance with some swing-phase activity. Cycle-to-cycle firing variability was far higher during active contraction, and brisk muscle stretches evoked rapid (disynaptic or trisynaptic) reflex arcs, arguing that mesencephalic-pre  
  - *Robot/sim:* Implement spindle Ia models in which gamma drive biases firing but rapid shortening unloads the spindle, and Ib firing gated by stance loading; simulated afferent phase profiles can then be validated against these recorded discharge patterns during walking.
- **Duysens and Pearson 1976** — [The role of cutaneous afferents from the distal hindlimb in the regulation of the step cycle of thalamic cats](https://doi.org/10.1007/bf00235013)  
  - animals: Cat · pathways: Cutaneous stance modification  
  - Cutaneous afferents from the distal hindlimb regulate the step cycle in thalamic cats: paw loading suppresses flexion during stance while stimulation can trigger and reset the cycle — load-related cutaneous feedback contributes to swing-stance transitions.
- **Grillner and Zangger 1975** — [How detailed is the central pattern generation for locomotion](https://doi.org/10.1016/0006-8993(75)90401-1)  
  - animals: Cat · pathways: Fictive locomotion without sensory feedback  
  - Deafferenting one or both hindlimbs of mesencephalic treadmill-walking cats leaves the detailed EMG activation pattern intact — brief muscle-specific bursts within both flexion and extension phases persist (more variable after bilateral deafferentation) — showing that the central program does not simply alternate flexors and extensors but sequentially starts and terminates each muscle at the correct instant, while afferents serve to handle external perturbations rather than to time the pattern.
- **Feldman and Orlovsky 1975** — [Activity of interneurons mediating reciprocal 1a inhibition during locomotion](https://doi.org/10.1016/0006-8993(75)90974-9)  
  - animals: Cat · pathways: Ia reciprocal inhibition  
  - The data obtained show that during locomotion there are at least two sources of inhibition of motoneurons of antagonistic muscles: (i) the activity of 1a afferents of the active muscle, mediating by corresponding interneurons; and (ii) signals coming through the same interneurons from the central mechanisms generating stepping movements.
- **Goslow 1973** — [The cat step cycle: Hind limb joint angles and muscle lengths during unrestrained locomotion](https://doi.org/10.1002/jmor.1051410102)  
  - animals: Cat · afferents: Ia, II, Ib  
  - The reference cat hindlimb kinematics dataset: from walking to galloping the F/E1/E2/E3 phase sequence is preserved while E3 compresses most and F expands, knee and ankle move in near unity, and stance extensors undergo a single stretch-shorten cycle, implying muscle spindles and tendon organs are driven primarily by lengthening and isometric contractions rather than passive stretch, with synchronous activation of Ia, group II, and Ib endings as muscles become active.  
  - *Robot/sim:* Use the F/E1/E2/E3 phase-duration scaling and knee-ankle unity as kinematic targets across speeds, and drive afferent models from lengthening/isometric contraction phases rather than absolute length to reproduce natural spindle and tendon-organ activation timing.
- **Hultborn 1971** — [Recurrent inhibition from motor axon collaterals of transmission in the Ia inhibitory pathway to motoneurones](https://doi.org/10.1113/jphysiol.1971.sp009487)  
  - animals: Cat · pathways: Ia inhibitory  
  - Ia Inhibitoru
- **Engberg and Lundberg 1969** — [An electromyographic analysis of muscular activity in the hindlimb of the cat during unrestrained locomotion.](https://doi.org/10.1111/j.1748-1716.1969.tb04415.x)  
  - animals: Cat · afferents: Ia  
  - Classic unrestrained-cat EMG study: hindlimb extensor activity is rather uniform across muscles while flexor activity is individualized by functional group, and the precise timing of extensor EMG onset argues against a reflex origin from limb receptors, supporting centrally programmed alternating extensor-flexor activation with possible Ia-mediated reflex regulation superimposed.  
  - *Robot/sim:* Program alternating extensor-flexor activation centrally with superimposed Ia feedback in a walker; ablate the central timing to test whether Ia reflexes alone can set extensor onset timing.
- **Jankowska et al. 1967** — [The Effect of DOPA on the Spinal Cord 5. Reciprocal organization of pathways transmitting excitatory action to](https://doi.org/10.1111/j.1748-1716.1967.tb03636.x)  
  - animals: Cat · pathways: Ia reciprocal inhibition  
  - DOPA evokes a reciprocal organization of excitatory pathways to flexor and extensor alpha motoneurones — the pharmacological foundation of the half-center model of fictive locomotion.
- **Jankowska 1967** — [The Effect of DOPA on the Spinal Cord 6. Half‐centre organization of interneurones transmitting effects from t](https://doi.org/10.1111/j.1748-1716.1967.tb03637.x)  
  - animals: Cat · afferents: Flexor reflex afferents, Ia · pathways: Ia presynaptic inhibition  
  - Classic DOPA study revealing the half-center organization of interneurons transmitting late flexor-reflex-afferent effects in spinal cats: reciprocally organized populations activated from ipsilateral versus contralateral FRA project to flexor versus extensor motoneurons, while a third type depolarizes Ia afferent terminals, adding presynaptic inhibition to the long-latency alternating discharge.  
  - *Robot/sim:* Implement reciprocally organized FRA-driven half-centers plus a presynaptic-inhibition channel onto Ia terminals; ablate the Ia presynaptic inhibition to test how afferent gating shapes the DOPA-like alternating rhythm.
- **Houk 1967** — [Responses of Golgi tendon organs to active contractions of the soleus muscle of the cat.](https://doi.org/10.1152/jn.1967.30.3.466)  
  - animals: Cat  
  - RESPONSES OF GOLGI TENDON ORGANS TO ACTIVE CONTRACTIONS OF THE SOLEUS MUSCLE OF THE CAT1 JAMES HOUK2 AND ELWOOD HENNEMAN Department of Physiolog-y, Harvard Medical School, Boston, Massachusetts (Received for publication July 29, 1966) WHEN A MUSCLE CONTRACTS it develops forces which are applied directly...
- **Shik 1966** — Control of walking and running by means of electric stimulation of the midbrain  
  - animals: Cat  
  - Electrical stimulation of the mesencephalic locomotor region elicits walking and running graded by stimulus intensity — the founding demonstration of the MLR and descending gait control.
- **Matthews 1964** — [Muscle Spindles and Their Motor Control](https://doi.org/10.1152/physrev.1964.44.2.219)  
  - animals: Cat  
  - Computational model study. The existence of two distinct types of afferent nerve ending, however, has been widely suspected since Ruflini described the histologically distinct primary and secondary endings, while very recent histological work strongly suggests the existence of two distinct types of intrafusal muscle fiber with independent motor...
- **Eccles and Lundberg 1958** — [Integrative pattern of Ia synaptic actions on motoneurones of hip and knee muscles.](https://doi.org/10.1113/jphysiol.1958.sp006101)  
  - animals: Cat · pathways: Ia monosynaptic excitation  
  - Classic integrative map of Ia synaptic actions on hip and knee motoneurones: Ia excitation diverges broadly to synergists and across joints, fractionated into flexor and extensor recipient pools.
- **Eccles and Lundberg 1957** — [The convergence of monosynaptic excitatory afferents on to many different species of alpha motoneurones](https://doi.org/10.1113/jphysiol.1957.sp005794)  
  - animals: Cat  
  - Intracellular recording from more than 400 cat spinal motoneurones defined the receptor field of group Ia monosynaptic excitation: a motoneurone receives EPSPs not only from its own muscle (homonymous) but from a characteristic set of synergists (heteronymous), establishing the convergent organization of the Ia monosynaptic pathway.
- **Brown 1911** — [The intrinsic factors in the act of progression in the mammal](https://doi.org/10.1098/rspb.1911.0077)  
  - animals: Mammals, Cat  
  - Rhythmic stepping persists in the low-spinal cat after complete deafferentation of the hindlimb muscles, and Graham Brown concluded the cycle is generated by paired antagonistic spinal centres whose balance sways through fatigue and rebound — the origin of the half-center concept — with proprioceptive input playing a regulative, not causative, role in grading steps to the terrain.
- **Sherrington 1910** — [Flexion-reflex of the limb, crossed extension-reflex, and reflex stepping and standing.](https://doi.org/10.1113/jphysiol.1910.sp001362)  
  - animals: Cat, Dog · pathways: Cutaneous flexor excitation  
  - Defined the flexion-reflex as a protective type-reflex of the whole limb and its companion crossed extension-reflex in the spinal cat and dog, and showed that reflex stepping fractionates them by phase: the extensor phase of one limb is reinforced by a crossed extension-reflex driven by the opposite limb's flexion, with stimulus-dependent reversal (Umkehr) marking the turning points of the step.
- **Smith ** — [Forms of Forward Quadrupedal Locomotion. III. A Comparison of Posture, Hindlimb Kinematics, and Motor Patterns](https://doi.org/10.1152/jn.1998.79.4.1702)  
  - animals: Cat  
  - Across downslope grades (25-100%) in walking cats, stance-phase yield increases at the ankle, swing-phase knee and ankle flexion decreases while ankle extension lowers the paw to contact, hip extensors fall silent and hip flexors instead brake the rate of hip extension, and ankle extensor bursts truncate around paw contact - a systematic reorganization of the motor pattern to counteract external braking forces rather than a rescaled level-gait pattern.  
  - *Robot/sim:* Gate hip extensor stance drive off and engage hip flexors as mid-stance brakes on downslope terrain, with ankle extensor bursts truncated around paw contact, and test whether a fixed CPG plus slope-dependent gating reproduces the measured downslope kinematics.
