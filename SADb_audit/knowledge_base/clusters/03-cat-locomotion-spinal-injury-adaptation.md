# Cluster 4: Cat locomotion & spinal-injury adaptation

94 papers · years 1910–2026 · 90 curated

## Top authors

- Frigon (7)
- Forssberg (4)
- Rossignol (4)
- Grillner and Zangger (3)
- Yakovenko (3)
- Hurteau (3)
- Harnie (3)
- Merlet (3)
- Pratt (2)
- Markin (2)

## Pathways represented

- Cutaneous flexor excitation (8)
- Cutaneous stance modification (8)
- Fictive locomotion without sensory feedback (5)
- Total afferent inhibition (1)
- Ia monosynaptic excitation (1)

## Key papers (by citations)

- **Sherrington 1910** — [Flexion-reflex of the limb, crossed extension-reflex, and reflex stepping and standing.](https://doi.org/10.1113/jphysiol.1910.sp001362)  
  - animals: Cat, Dog · pathways: Cutaneous flexor excitation  
  - Defined the flexion-reflex as a protective type-reflex of the whole limb and its companion crossed extension-reflex in the spinal cat and dog, and showed that reflex stepping fractionates them by phase: the extensor phase of one limb is reinforced by a crossed extension-reflex driven by the opposite limb's flexion, with stimulus-dependent reversal (Umkehr) marking the turning points of the step.
- **Grillner and Zangger 1979** — [On the central generation of locomotion in the low spinal cat](https://doi.org/10.1007/BF00235671)  
  - animals: Cat · pathways: Total afferent inhibition  
  - Cited by Geyer and Herr, 2010.
Total afferent inhibition.

A central network of neurones in the spinal cord has been shown to produce a rhythmic motor output similar to locomotion after suppression of all afferent inflow. The experiments were performed mainly in acute spinal cats (th. 12), which had received DOPA i.v. and the monoamine oxidase inhibitor Nialamide. In some preparations all dorsal roots supplying the spinal cord were transected, in others phasic afferent activity was suppressed by curarization. The activity was recorded as neurograms from nerve filaments or as electromyograms.


- **Barbeau and Rossignol 1987** — [Recovery of locomotion after chronic spinalization in the adult cat](https://doi.org/10.1016/0006-8993(87)91442-9)  
  - animals: Cat  
  - Adult cats spinalized at T13 recover plantar weight-bearing treadmill walking after weeks to months of interactive training, with cycle duration, stance/swing timing, joint coupling, and hindlimb EMG resembling intact cats at comparable speeds — establishing the chronic spinal cat as the standard model for testing how training, drugs, and other treatments influence locomotor recovery below the lesion.
- **Goslow 1973** — [The cat step cycle: Hind limb joint angles and muscle lengths during unrestrained locomotion](https://doi.org/10.1002/jmor.1051410102)  
  - animals: Cat · afferents: Ia, II, Ib  
  - The reference cat hindlimb kinematics dataset: from walking to galloping the F/E1/E2/E3 phase sequence is preserved while E3 compresses most and F expands, knee and ankle move in near unity, and stance extensors undergo a single stretch-shorten cycle, implying muscle spindles and tendon organs are driven primarily by lengthening and isometric contractions rather than passive stretch, with synchronous activation of Ia, group II, and Ib endings as muscles become active.  
  - *Robot/sim:* Use the F/E1/E2/E3 phase-duration scaling and knee-ankle unity as kinematic targets across speeds, and drive afferent models from lengthening/isometric contraction phases rather than absolute length to reproduce natural spindle and tendon-organ activation timing.
- **Lovely 1986** — [Effects of training on the recovery of full-weight-bearing stepping in the adult spinal cat.](https://doi.org/10.1016/0014-4886(86)90094-4)  
  - animals: Cat  
  - Most adult spinal cats (14 of 16) recover full-weight-bearing treadmill stepping within a month of transection when cued by tail pinch or crimping, and daily treadmill training emphasizing complete weight bearing raises maximum treadmill speed about 2.6-fold over untrained cats (0.619 versus 0.240 m/s) - task-specific training markedly improves spinal stepping capacity below the lesion.  
  - *Robot/sim:* Model a below-lesion spinal stepping network with training-dependent plasticity driven by weight-bearing activity; disabling the adaptation should reproduce the untrained-cat speed plateau.
- **Engberg and Lundberg 1969** — [An electromyographic analysis of muscular activity in the hindlimb of the cat during unrestrained locomotion.](https://doi.org/10.1111/j.1748-1716.1969.tb04415.x)  
  - animals: Cat · afferents: Ia  
  - Classic unrestrained-cat EMG study: hindlimb extensor activity is rather uniform across muscles while flexor activity is individualized by functional group, and the precise timing of extensor EMG onset argues against a reflex origin from limb receptors, supporting centrally programmed alternating extensor-flexor activation with possible Ia-mediated reflex regulation superimposed.  
  - *Robot/sim:* Program alternating extensor-flexor activation centrally with superimposed Ia feedback in a walker; ablate the central timing to test whether Ia reflexes alone can set extensor onset timing.
- **Forssberg 1979** — [Stumbling corrective reaction: a phase-dependent compensatory reaction during locomotion](https://doi.org/10.1152/jn.1979.42.4.936)  
  - animals: Cat · afferents: Cutaneous, Nociceptive · pathways: Cutaneous flexor excitation; Cutaneous stance modification  
  - The classic characterization of the stumbling corrective reaction: dorsum contact during swing evokes short-latency flexor activation lifting the paw over the obstacle via early (~10 ms) and late (~25 ms) pathways to knee flexors, whereas the same stimulus during support inhibits-then-excites extensors and boosts the next swing's flexion, phase-dependent reflex gating organized by the spinal locomotor generator, with painful input instead evoking withdrawal throughout the cycle.  
  - *Robot/sim:* Implement a phase-switched stumble reflex: paw-dorsum tactile input drives short-latency flexor bursts during swing for obstacle clearance and stance-phase extensor modulation with enhanced next-swing flexion; ablate the reflex to quantify its contribution to fall avoidance.
- **Nilsson 1985** — [Changes in leg movements and muscle activity with speed of locomotion and mode of progression in humans](https://doi.org/10.1111/j.1748-1716.1985.tb07612.x)  
  - animals: Human  
  - Human treadmill dataset (walking 0.4-3.0 m/s, running 1-9 m/s) showing speed adaptation proceeds by increasing both frequency and amplitude of leg movements, with mode-specific EMG re-timing: rectus femoris shifts from knee-extension to hip-flexion emphasis with speed, gastrocnemius lateralis pre-activates before foot contact in running but not walking, and tibialis anterior peaks before rather than after touchdown in running - the basic stride structure matches other animals, suggesting shared neural control.  
  - *Robot/sim:* Use these speed- and mode-dependent EMG phase relationships (rectus femoris role switch, gastrocnemius pre-contact onset, tibialis anterior timing) as quantitative validation targets for a CPG plus reflex biped model across walking and running speeds.
- **Rossignol 2011** — [Neural Control of Stereotypic Limb Movements](https://doi.org/10.1002/cphy.cp120105)  
  - afferents: Cutaneous  
  - Comprehensive chapter on the neural control of stereotypic limb movements, spanning locomotor kinematics and EMG across speeds, slopes, and directions; interlimb coordination; spinal, complete, and partial lesion preparations and pharmacology; central pattern generation and locomotion interneurons; and proprioceptive and cutaneous afferent control with mechanisms of reflex modulation and supraspinal initiation.
- **Forssberg 1980** — [The locomotion of the low spinal cat. II. Interlimb coordination](https://doi.org/10.1111/j.1748-1716.1980.tb06534.x)  
  - animals: Cat  
  - Chronic spinal kittens on a split-belt treadmill maintain a common hindlimb rhythm across two- to threefold belt speed differences by prolonging the flexion or first extension phase of the fast-belt limb and shortening the flexion phase of the slow-belt limb, with extreme differences producing 2, 3, or even 4 fast-limb steps per slow-limb cycle; support-phase duration mainly follows each limb's own belt speed. All phases of the step cycle are thus modifiable, with several spinal mechanisms coordinating the limbs.  
  - *Robot/sim:* Implement spinal interlimb coordination with per-limb phase modification under asymmetric timing commands (a split-belt analog), including shared rhythm up to two- to threefold speed differences and integer multiple stepping beyond that.
- **Forssberg 1977** — [Phasic gain control of reflexes from the dorsum of the paw during spinal locomotion.](https://doi.org/10.1016/0006-8993(77)90710-7)  
  - animals: Cat · afferents: Cutaneous · pathways: Cutaneous flexor excitation; Cutaneous stance modification  
  - In chronic spinal cats walking on a treadmill, tactile stimulation of the paw dorsum during swing evokes a short-latency flexion response with concomitant crossed extension, while the same stimulus during stance increases ipsilateral extension, a phase-dependent reflex reversal that compensates for unpredicted obstacles by lifting the paw in swing and reinforcing support in stance. The responses are well adapted to ongoing locomotion and leave interlimb coordination intact, except when delivered as the foot approaches the ground after flexion, where the alternating pattern is disturbed.  
  - *Robot/sim:* Implement a phase-switched paw-dorsum cutaneous reflex, flexion plus crossed extension in swing and ipsilateral extension reinforcement in stance, in a walking model; ablating it tests its contribution to obstacle compensation.
- **Forssberg 1980** — [The locomotion of the low spinal cat. I. Coordination within a hindlimb.](https://doi.org/10.1111/j.1748-1716.1980.tb06533.x)  
  - animals: Cat  
  - The low spinal cat walks with coordinated hindlimb bursts: spinal locomotion is possible below transection and is shaped by afferent input — foundation for spinal locomotion training concepts.
- **Grillner and Zangger 1975** — [How detailed is the central pattern generation for locomotion](https://doi.org/10.1016/0006-8993(75)90401-1)  
  - animals: Cat · pathways: Fictive locomotion without sensory feedback  
  - Deafferenting one or both hindlimbs of mesencephalic treadmill-walking cats leaves the detailed EMG activation pattern intact — brief muscle-specific bursts within both flexion and extension phases persist (more variable after bilateral deafferentation) — showing that the central program does not simply alternate flexors and extensors but sequentially starts and terminates each muscle at the correct instant, while afferents serve to handle external perturbations rather than to time the pattern.
- **Perret and Cabelguen 1980** — [Main characteristics of the hindlimb locomotor cycle in the decorticate cat with special reference to bifuncti](https://doi.org/10.1016/0006-8993(80)90207-3)  
  - animals: Cat · afferents: Ia, Flexor reflex afferents  
  - In decorticate cats, all studied hindlimb muscles show alpha-gamma coactivation, pure flexors and extensors alternate simply, and bifunctional pluriarticular muscles receive both flexor and extensor commands, dividing the locomotor cycle into a flexion phase, an extension phase, and two transition phases. The relative weight of these excitations depends on interactions between central commands and peripheral inflow, especially afferents acting according to the flexor reflex pattern, providing a scheme whereby an initially simple rhythmic command plus afferent gating produces the biomechanicall  
  - *Robot/sim:* Implement bifunctional motoneuron pools receiving weighted flexor and extensor CPG commands modulated by flexor-reflex-pattern afferent input; the model should generate transition-phase complexity and task-adapted output from a simple alternating rhythm.
- **Dietz 1994** — [HUMAN NEURONAL INTERLIMB COORDINATION DURING SPLIT-BELT LOCOMOTION](https://doi.org/10.1007/bf00227344)  
  - animals: Human  
  - Human split-belt adaptation occurs within 10-20 stride cycles through reorganization of the stride cycle (support shortens and swing lengthens on the fast leg), with ipsilateral gastrocnemius activity scaling almost linearly with belt speed under proprioceptive control while the contralateral tibialis anterior is centrally modulated, revealing an ipsilateral-proprioceptive / contralateral-central coupling between the support phase of one leg and the swing phase of the other.  
  - *Robot/sim:* Implement split-belt adaptation in a biped model as ipsilateral proprioceptive extensor feedback that scales stance muscle activity with belt speed, plus a central contralateral coupling from support to opposite-limb swing; ablating the proprioceptive branch should abolish automatic speed matching.
- **Pearson and Rossignol 1991** — [Fictive motor patterns in chronic spinal cats](https://doi.org/10.1152/jn.1991.66.6.1874)  
  - animals: Cat · afferents: Cutaneous · pathways: Fictive locomotion without sensory feedback  
  - Chronic spinal cats generate three distinctly different fictive motor patterns, locomotion, paw shake (approximately 8 Hz rhythmic activity), and slow rhythmic leg flexion, depending only on the site and mode of stimulation (perineal skin, water jet to the paw, paw squeeze), with clonidine facilitating evocation. Passive leg movement from extension to flexion progressively shortened flexor bursts, lengthened the cycle period, and made the pattern harder to evoke, and trained late-spinal animals expressed more complex, position-dependent patterns than untrained early-spinal ones, demonstrating   
  - *Robot/sim:* Build a pattern generator supporting multiple coexistent fictive modes (stepping, paw shake, rhythmic flexion) switched by input site, and modulate flexor burst termination by a position-related afferent signal to reproduce the effects of passively moving the legs from extension to flexion.
- **Beloozerova and Sirota 1993** — [The role of the motor cortex in the control of accuracy of locomotor movements in the cat.](https://doi.org/10.1113/jphysiol.1993.sp019498)  
  - animals: Cat  
  - Motor cortex controls the accuracy of locomotor movements: cortical neurons encode precise foot placement, and cortical inactivation impairs accuracy without abolishing walking.
- **Armstrong 1986** — [Supraspinal contributions to the initiation and control of locomotion in the cat.](https://doi.org/10.1016/0301-0082(86)90021-3)  
  - animals: Cat
- **Barriere 2008** — [Prominent Role of the Spinal Central Pattern Generator in the Recovery of Locomotion after Partial Spinal Cord](https://doi.org/10.1523/jneurosci.5692-07.2008)  
  - animals: Cat  
  - A dual-lesion paradigm in cats shows the spinal CPG undergoes plastic changes after partial spinal cord injury that are shaped by locomotor training: cats trained on the treadmill after a partial thoracic lesion expressed bilateral hindlimb locomotion within hours of a subsequent complete transection, whereas untrained cats walked asymmetrically with the limb on the partially lesioned side recovering first. Re-expression of the hindlimb locomotor pattern after partial SCI therefore arises mostly from intrinsic changes below the lesion in the CPG and afferent inputs rather than from remnant des  
  - *Robot/sim:* Model a spinal CPG whose synapses adapt during repetitive rhythmic 'training' activation, then remove descending drive completely; the trained network should re-express bilateral rhythm rapidly while the untrained network recovers slowly and asymmetrically, reproducing the lesion hierarchy.
- **Duysens and Loeb 1980** — [Modulation of ipsi- and contralateral reflex responses in unrestrained walking cats](https://doi.org/10.1152/jn.1980.44.5.1024)  
  - animals: Cat · afferents: Cutaneous · pathways: Cutaneous flexor excitation; Cutaneous stance modification  
  - In unrestrained walking cats, cutaneous stimulation evokes phase-modulated responses: two excitatory peaks (about 10 and 25 ms) in flexors that grow toward the end of stance, early inhibition then late excitation (P3, about 35 ms) in extensors when stimuli fall in early stance, and crossed responses in the contralateral limb at 20-25 ms. Concludes the walking cat shows modulation of transmission in a flexor-excitatory/extensor-inhibitory pathway - likely by the flexor part of the spinal locomotor oscillator - rather than strict reflex reversal, plus specialized flexor-inhibitory and extensor-e  
  - *Robot/sim:* Implement cutaneous reflex pathways with locomotor-phase gating: flexor excitation whose gain peaks at late stance, extensor inhibition during early stance, and crossed extensor/flexor responses; ablating the phase gating tests the oscillator-modulation hypothesis versus fixed reflex reversal.
- **Orlovsky 1972** — [Activity of vestibulospinal neurons during locomotion.](https://doi.org/10.1016/0006-8993(72)90007-8)
- **Gagnon 2003** — [Contribution of Cutaneous Inputs From the Hindpaw to the Control of Locomotion. I. Intact Cats](https://doi.org/10.1152/jn.00496.2003)  
  - animals: Cat · afferents: Cutaneous  
  - Bilateral hindpaw cutaneous denervation at the ankle in intact cats barely disrupts level treadmill walking but severely impairs ladder and incline walking early on; flexor EMG remains permanently elevated and mediolateral ground reaction forces increase by 200 percent, showing cutaneous feedback matters most for demanding locomotor contexts and that active adaptive mechanisms compensate without fully restoring the control pattern. The companion paper shows these same inputs become critical for foot placement after spinalization.  
  - *Robot/sim:* Add foot/ankle cutaneous afferent channels to a walking simulation and ablate them: the model should reproduce near-normal level walking with degraded ladder/incline behavior, persistently elevated flexor drive, and increased mediolateral ground reaction forces.
- **Grillner and Zangger 1984** — [The effect of dorsal root transection on the efferent motor pattern in the cat's hindlimb during locomotion.](https://doi.org/10.1111/j.1748-1716.1984.tb07400.x)  
  - animals: Cat · pathways: Fictive locomotion without sensory feedback  
  - After transecting all dorsal roots from one hindlimb, mesencephalic walking cats retain the limb's complex motor pattern, including the double-burst knee flexors and appropriately timed toe dorsiflexor, although with greater variability and occasional breakdown, showing the central network generates the full pattern without phasic afferent input while afferents stabilize and fine-tune it.  
  - *Robot/sim:* Ablate all sensory feedback channels in a CPG walker and test whether the joint-specific burst structure (double-burst flexors) persists with increased variability, a lesion test separating central pattern generation from afferent stabilization.
- **Carlson-Kuhta 1998** — [Forms of forward quadrupedal locomotion. II. A comparison of posture, hindlimb kinematics, and motor patterns ](https://doi.org/10.1152/jn.1998.79.4.1687)  
  - animals: Cat  
  - Per the abstract provided: across steep downslope grades in walking cats, stance-phase yield increases at the ankle, swing knee and ankle flexion decreases with ankle extension lowering the paw to contact, hip extensors fall silent while hip flexors brake the rate of hip extension, and ankle extensor bursts truncate around paw contact - the motor pattern reorganizes to counteract external braking forces rather than simply rescaling level gait.  
  - *Robot/sim:* Gate hip extensor drive off and use hip flexors as stance brakes on downslopes, with ankle extensor bursts centered on paw contact, testing whether a fixed CPG with slope-dependent gating reproduces the graded downslope kinematics.
- **Frigon 2017** — [The neural control of interlimb coordination during mammalian locomotion.](https://doi.org/10.1152/jn.00978.2016)  
  - animals: Mammals, Human  
  - Review of interlimb coordination in mammalian locomotion: coordination between forelimb and hindlimb spinal networks is achieved by mechanisms intrinsic to the spinal cord, somatosensory feedback from the limbs, and supraspinal pathways; incomplete spinal cord injury disrupts this coordination, but lesion-based inference is confounded by compensatory strategies, redundant control, and plasticity in remaining circuits.  
  - *Robot/sim:* Build a locomotion model with separate interlimb-coupling modules (intraspinal coupling, limb somatosensory feedback, supraspinal modulation) and ablate them individually to reproduce the coordination disruptions seen after incomplete spinal lesions.
