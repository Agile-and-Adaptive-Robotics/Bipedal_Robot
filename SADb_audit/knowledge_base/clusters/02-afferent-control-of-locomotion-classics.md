# Cluster 3: Afferent control of locomotion (classics)

141 papers · years 1957–2025 · 130 curated

## Top authors

- Duysens (7)
- Pearson (5)
- Eccles and Lundberg (3)
- Prochazka (3)
- Quevedo (3)
- Perreault (3)
- Jankowska (2)
- Duysens and Pearson (2)
- Procházka and Prochazka (2)
- Andersson and Grillner (2)

## Pathways represented

- Ib excitatory (12)
- Ib stance to swing (10)
- Cutaneous flexor excitation (8)
- Ib disynaptic excitation (7)
- Ia monosynaptic excitation (5)
- Ia or II stance to swing (5)
- Ib disynaptic inhibition (5)
- Cutaneous stance modification (4)
- Ib inhibition (4)
- type II excitatory (4)
- Ia or II swing to stance (3)
- Ia monosynaptic (2)

## Key papers (by citations)

- **Rossignol 2006** — [Dynamic Sensorimotor Interactions in Locomotion](https://doi.org/10.1152/physrev.00028.2005)  
  - afferents: Cutaneous · pathways: Cutaneous stance modification  
  - Definitive synthesis of locomotion as a dynamic interaction between a genetically determined spinal CPG and sensory feedback: extensor proprioceptive input adjusts stance timing and amplitude to locomotor speed while being silenced in the opposite phase, and skin afferents predominantly correct limb and foot placement during stance on uneven terrain. Transmission in locomotor pathways is modulated in a state- and phase-dependent manner through shared presynaptic, interneuronal, and motoneuronal mechanisms, keeping the central program and feedback congruous.  
  - *Robot/sim:* Build the controller as CPG plus phase-gated feedback: extensor proprioceptive channels active only during stance (silenced in swing) adjusting timing and amplitude to speed, and cutaneous placement-correction reflexes active during stance — then test gait adaptation to speed and uneven terrain.
- **Capaday et al. 1986** — [Amplitude modulation of the soleus H-reflex in the human during walking and standing](https://doi.org/10.1523/jneurosci.06-05-01308.1986)  
  - animals: Human · pathways: Ia monosynaptic excitation  
  - The soleus H-reflex is strongly modulated over the step cycle during walking versus standing; central gating of the monosynaptic Ia pathway is phase- and task-dependent.
- **Jankowska 1992** — [INTERNEURONAL RELAY IN SPINAL PATHWAYS FROM PROPRIOCEPTORS](https://doi.org/10.1016/0301-0082(92)90024-9)  
  - animals: Mammals  
  - review pathway and sensory
- **Eccles and Lundberg 1957** — [The convergence of monosynaptic excitatory afferents on to many different species of alpha motoneurones](https://doi.org/10.1113/jphysiol.1957.sp005794)  
  - animals: Cat  
  - Intracellular recording from more than 400 cat spinal motoneurones defined the receptor field of group Ia monosynaptic excitation: a motoneurone receives EPSPs not only from its own muscle (homonymous) but from a characteristic set of synergists (heteronymous), establishing the convergent organization of the Ia monosynaptic pathway.
- **Forssberg 1975** — [Phase dependent reflex reversal during walking in chronic spinal cats](https://doi.org/10.1016/0006-8993(75)91013-6)
- **Harkema 1997** — [Human Lumbosacral Spinal Cord Interprets Loading During Stepping](https://doi.org/10.1152/jn.1997.77.2.797)  
  - animals: Human  
  - In spinal-cord-injured humans during assisted stepping, soleus, gastrocnemius, and tibialis anterior EMG amplitude scales directly with peak limb load regardless of supraspinal input, and load predicts phasic EMG better than muscle-tendon stretch or stretch velocity - the lumbosacral cord uses limb-loading cues, not just stretch reflexes, to grade efferent output for stepping.  
  - *Robot/sim:* Implement stance-gated limb-load feedback that scales extensor and flexor motoneuron pool gain; ablating it should mimic reduced load-bearing and depress stepping EMG, reproducing body-weight-support effects.
- **Jami 1992** — [Golgi tendon organs in mammalian skeletal muscle: functional properties and central actions](https://doi.org/10.1152/physrev.1992.72.3.623)
- **Grillner and Rossignol 1978** — [On the initiation of the swing phase of locomotion in chronic spinal cats.](https://doi.org/10.1016/0006-8993(78)90973-3)  
  - animals: Cat  
  - In chronic spinal cats walking on a treadmill, a limb held by the paw lifts off and rejoins walking once the hip is brought far enough back, at a hip angle very close to that of normal swing initiation, and lift-off also tends to occur at particular contralateral cycle phases. Establishes hip position and contralateral step-cycle phase as two key factors determining swing initiation, i.e., the stance-to-swing transition is regulated by proprioceptive and interlimb signals rather than a fixed internal clock.  
  - *Robot/sim:* Implement swing initiation in a walking model as a hip-extension position threshold gated by contralateral cycle phase; ablating the hip-position input should delay or prevent swing onset, reproducing the held-paw behavior.
- **Eccles and Lundberg 1957** — [Synaptic actions on motoneurones caused by impulses in Golgi tendon organ afferents.](https://doi.org/10.1113/jphysiol.1957.sp005849)
- **Duysens 2000** — [Load-regulating mechanisms in gait and posture: Comparative aspects.](https://doi.org/10.1152/physrev.2000.80.1.83)  
  - pathways: Ib inhibition  
  - Type Ib feedback has an inhibitory effect on muscles in non-moving animals
- **Duysens and Pearson 1980** — [Inhibition of flexor burst generation by loading ankle extensor muscles in walking cats](https://doi.org/10.1016/0006-8993(80)90206-1)  
  - animals: Cat · pathways: Ib stance to swing  
  - It has been found that decreased load detected by Golgi organs encourages swing muscles to activate
- **Houk 1967** — [Responses of Golgi tendon organs to active contractions of the soleus muscle of the cat.](https://doi.org/10.1152/jn.1967.30.3.466)  
  - animals: Cat  
  - RESPONSES OF GOLGI TENDON ORGANS TO ACTIVE CONTRACTIONS OF THE SOLEUS MUSCLE OF THE CAT1 JAMES HOUK2 AND ELWOOD HENNEMAN Department of Physiolog-y, Harvard Medical School, Boston, Massachusetts (Received for publication July 29, 1966) WHEN A MUSCLE CONTRACTS it develops forces which are applied directly...
- **Procházka and Prochazka 1989** — [Sensorimotor gain control: a basic strategy of motor systems?](https://doi.org/10.1016/0301-0082(89)90004-x)
- **Zehr and Stein 1999** — [What functions do reflexes serve during human locomotion](https://doi.org/10.1016/s0301-0082(98)00081-1)  
  - animals: Lamprey, Cat, Human · afferents: Cutaneous · pathways: Cutaneous flexor excitation; Ia monosynaptic excitation  
  - Review establishing that reflexes during locomotion are task-, phase-, and context-dependent, with a division of labor: cutaneous reflexes alter swing-limb trajectory to avoid stumbling, stretch reflexes stabilize limb trajectory and assist force production during stance, and load receptor reflexes support body weight and influence step-cycle timing - functions dynamically reassigned across the cycle and clinically exploitable after neurotrauma.  
  - *Robot/sim:* Implement phase- and task-gated cutaneous (swing trajectory), Ia stretch (stance stabilization and force), and load-receptor (weight support, cycle timing) pathways in a walker; ablate each to reproduce the predicted functional losses.
- **Conway 1987** — [Proprioceptive input resets central locomotor rhythm in the spinal cat.](https://doi.org/10.1007/bf00249807)  
  - animals: Cat · pathways: Ib stance to swing  
  - Extensor group I afferents have access to central rhythm generators
 >Important in reflex regulation of stepping.


Rythm generators came mainly from Golgi Tendon Organ Ib afferents for group I

Increased load of limb extensor during the stance phase enhance and prolong extensor activity while simultaneously delaying the transition to the swing phase of the step cycle.
- **Dietz 2002** — [Proprioception and locomotor disorders](https://doi.org/10.1038/nrn939)  
  - animals: Human  
  - Review paper
Key Points
Locomotion in mammals depends on central pattern generators — networks of spinal interneurons that can produce rhythmic outputs independently of any modulatory input. However, the spinal activity pattern is influenced by inputs from peripheral afferents, brainstem nuclei and cortical motor centres. Central pattern generators must select the appropriate inputs at each stage of movement and according to external conditions.

Afferent inputs that influence gait include the short-latency stretch reflex, which is mediated by excitatory monosynaptic connections between sensor
- **Pearson 1995** — [Proprioceptive regulation of locomotion.](https://doi.org/10.1016/0959-4388(95)80107-3)  
  - animals: Cat, Insects, Arthropods · pathways: Ib stance to swing; Ia or II stance to swing  
  - Review of proprioceptive regulation across walking systems: during locomotion, Golgi tendon organ feedback from extensor muscles reverses from inhibition to excitation, maintaining stance while the extensors are loaded and helping time the stance-to-swing transition, while primary and secondary spindle afferents also influence the timing of the rhythm — establishing phase-dependent reflex reversal as a general feature of cat and arthropod walking.
- **Zehr and Duysens 2004** — [Regulation of Arm and Leg Movement during Human Locomotion](https://doi.org/10.1177/1073858404264680)  
  - animals: Human  
  - Synthesis that both arms and legs during human locomotion are regulated by central pattern generators, with sensory feedback regulating CPG activity and assisting interlimb coordination; coupling strength is stronger between the legs than between the arms, but all four limbs are similarly governed by CPG activity and reflex control.  
  - *Robot/sim:* Implement four limb oscillators (a quadrupedal-style CPG) with weaker arm-leg than leg-leg coupling plus afferent phase modulation of each oscillator; test whether arm-swing entrainment matches human interlimb coordination data.
- **Eccles and Lundberg 1958** — [Integrative pattern of Ia synaptic actions on motoneurones of hip and knee muscles.](https://doi.org/10.1113/jphysiol.1958.sp006101)  
  - animals: Cat · pathways: Ia monosynaptic excitation  
  - Classic integrative map of Ia synaptic actions on hip and knee motoneurones: Ia excitation diverges broadly to synergists and across joints, fractionated into flexor and extensor recipient pools.
- **Pearson and Collins 1993** — [Reversal of the influence of group Ib afferents from plantaris on activity in medial gastrocnemius muscle duri](https://doi.org/10.1152/jn.1993.70.3.1009)  
  - animals: Cat · afferents: Ib, Ia · pathways: Ib excitatory; Ib inhibition  
  - In clonidine-treated acute and chronic spinal cats, group I stimulation of the plantaris nerve entrained the locomotor rhythm and, during locomotor activity, group Ib afferents from plantaris exerted an excitatory action on medial gastrocnemius bursts (30-50 ms latency) through the extensor half-center of the rhythm generator — reversing the inhibitory effect the same stimuli produced on tonic activity without locomotion; selective Ia activation by muscle vibration neither entrained the rhythm nor augmented the bursts.  
  - *Robot/sim:* Implement state-dependent Ib feedback in a half-center CPG: autogenic inhibition at rest switching to excitation of the extensor half-center during locomotion, with rhythm entrainment by group I stimulation but not by Ia-specific input; ablate the reversal to test loss of load-dependent stance burst augmentation.
- **Pearson 2004** — [Generating the walking gait: role of sensory feedback.](https://doi.org/10.1016/s0079-6123(03)43012-4)  
  - animals: Cat  
  - Synthesizes cat walking evidence that feedback from muscle proprioceptors establishes the timing of major phase transitions, contributes to burst production, generates some features of the motor pattern, and is required for adaptive modification after alterations in leg mechanics; argues that afferent signals likely reorganize the functioning of central networks, making 'afferent modulation of a hard-wired CPG' too simplistic a framework.  
  - *Robot/sim:* Model phase-transition timing and burst production as afferent-controlled processes rather than fixed CPG timing; ablate individual proprioceptive channels and quantify the loss of adaptive motor-pattern modification to altered leg mechanics.
- **Dietz and Harkema 2004** — [Locomotor activity in spinal cord-injured persons](https://doi.org/10.1152/japplphysiol.00942.2003)  
  - animals: Human  
  - After a spinal cord injury (SCI) of the cat or rat, neuronal centers below the level of lesion exhibit plasticity that can be exploited by specific training paradigms. In individuals with complete or incomplete SCI, human spinal locomotor centers can be activated and modulated by locomotor training (facilitating stepping movements of the legs using body weight support on a treadmill to provide appropriate sensory cues). Individuals with incomplete SCI benefit from locomotor training such that they improve their ability to walk over ground. Load- or hip joint-related afferent input seems to be 
- **Dietz and Duysens 2000** — [Significance of load receptor input during locomotion: a review](https://doi.org/10.1016/s0966-6362(99)00052-1)  
  - animals: Human, Cat · pathways: Ib stance to swing  
  - Review making the case that extensor load receptors are central to locomotor control: in the cat, Golgi tendon organ input switches function during walking from Ib inhibition to extensor facilitation proportional to load, and in humans leg extensor activation during stance scales with body weight — one load-regulatory mechanism spanning quadrupedal and bipedal gait.
- **Bosco and Poppele 2001** — [Proprioception from a spinocerebellar perspective.](https://doi.org/10.1152/physrev.2001.81.2.539)  
  - Develops the dorsal spinocerebellar tract system as a model showing that spinal proprioceptive processing produces a global representation of whole-limb parameters rather than a muscle-by-muscle or joint-by-joint code, with the functional organization originating in the biomechanical linkages of the limb and a distributed spinal processing network.  
  - *Robot/sim:* Replace muscle-by-muscle proprioceptive feedback in a locomotor controller with a spinocerebellar-like whole-limb kinematic encoding stage, and test whether the global representation improves interjoint coordination and balance relative to a joint-by-joint representation.
- **Windhorst 2007** — [Muscle proprioceptive feedback and spinal networks](https://doi.org/10.1016/j.brainresbull.2007.03.010)  
  - afferents: Ia, II, Ib  
  - Review assigning computational roles to segmental muscle spindle, Golgi tendon organ, and Renshaw cell feedback: generation of anti-gravity thrust during stance, locomotor phase timing, linearization of muscle nonlinearities, lever-arm and fatigue compensation, synergy formation, and selection of perturbation responses, a functional catalog of what spinal feedback contributes to motor control.  
  - *Robot/sim:* Implement spindle (length/velocity), tendon-organ (force), and Renshaw recurrent feedback in a neuromuscular locomotion model and ablate each pathway to test the review's functional catalog: phase timing, fatigue compensation, linearization, and synergy stabilization.
