# Feedback pathway: Ib stance to swing

17 papers in the corpus.

- **Domínguez-Rodríguez 2020** — [Candidate Interneurons Mediating the Resetting of the Locomotor Rhythm by Extensor Group I Afferents in the Ca](https://doi.org/10.1016/j.neuroscience.2020.09.017)  
  - animals: Cat · afferents: Ib · pathways: Ib stance to swing; Ib swing to stance; Ib excitatory  
  - Extensor group I afferent stimulation resets fictive locomotion to extension - prolonging ongoing extension and terminating ongoing flexion - through a polysynaptic excitation of extensor motoneurons (latencies near 3.5-4.0 ms, compatible with three interposed interneurons) that replaces the classical Ib non-reciprocal inhibition during locomotion. Extension-phase interneurons receiving short-latency group I excitation satisfy the criteria for the pathway that resets the rhythm and may belong to the rhythm-generating layer of the CPG.  
  - *Robot/sim:* Implement extensor group I load afferents as inputs to rhythm-layer interneurons that both prolong the extensor half-center and reset flexion to extension; ablating them should remove load-dependent stance prolongation and rhythm resetting.
- **Kiehn 2016** — [Decoding the organization of spinal circuits that control locomotion](https://doi.org/10.1038/nrn.2016.9)  
  - pathways: Ia reciprocal inhibition; Ia or II stance to swing; Ib stance to swing  
  - Review paper.
Abstract | Unravelling the functional operation of neuronal networks and linking cellular activity
to specific behavioural outcomes are among the biggest challenges in neuroscience. In this
broad field of research, substantial progress has been made in studies of the spinal networks
that control locomotion. Through united efforts using electrophysiological and molecular
genetic network approaches and behavioural studies in phylogenetically diverse experimental
models, the organization of locomotor networks has begun to be decoded. The emergent themes
from this research are that t
- **Duysens 2013** — [The flexion synergy, mother of all synergies and father of new models of gait](https://doi.org/10.3389/fncom.2013.00014)  
  - animals: Human · pathways: Ib stance to swing  
  - The flexion synergy of Sherrington's flexor reflex, the withdrawal-reflex modules of Schouenborg, and the neonatal locomotor modules of Dominici largely overlap, and the end-of-stance facilitation of the flexion synergy — broadly afferent-facilitated but load-suppressed — points to a flexor burst generator as the core of an asymmetric CPG model whose afferent gating has already been implemented in walking bipedal robots.
- **Pearson 2008** — [Role of sensory feedback in the control of stance duration in walking cats](https://doi.org/10.1016/j.brainresrev.2007.06.014)  
  - animals: Cat · pathways: Ib stance to swing; Ia or II stance to swing  
  - Review Paper
- **Ekeberg and Pearson 2005** — [Computer simulation of stepping in the hind legs of the cat: an examination of mechanisms regulating the stanc](https://doi.org/10.1152/jn.00065.2005)  
  - animals: Cat · pathways: Ib stance to swing  
  - A three-dimensional simulation of cat hind-leg stepping shows that stance termination governed by ankle extensor force (load) signals — alone or combined with hip position — yields stable stepping and correct interleg timing, whereas the hip position signal alone fails; mutual inhibition between controllers restores stability but not correct timing. Coordination depends critically on load-sensitive signals from each leg, with mechanical linkages mediated by these signals playing a significant role in establishing the alternating gait.  
  - *Robot/sim:* Reproduce directly: gate the stance-to-swing transition on thresholded ankle-extensor force (with and without hip-angle gating) in a hind-leg/biped model, ablate each channel, and verify the load channel is necessary for stable alternation and correct contralateral timing.
- **Ivashko 2003** — [Modeling the spinal cord neural circuitry controlling cat hindlimb movement during locomotion](https://doi.org/10.1016/S0925-2312(02)00832-9)  
  - animals: Cat · pathways: Ib disynaptic excitation; Ib disynaptic excitation; Ib stance to swing; Ib swing to stance; Ia stance to swing; Ia swing to stance; Ia monosynaptic; Ia monosynaptic excitation; Mechanosensory monosynaptic excitation; Ia disynaptic inhibition  
  - Abstract
A computational model of the spinal cord neural circuitry that controls locomotor movements of simulated cat hindlimbs. The neural circuitry includes two central pattern generators integrated with reflex circuits. All neurons were modeled in the Hodgkin–Huxley style. The musculoskeletal system includes two three-joint hindlimbs and the trunk. Each
hindlimb is actuated by nine one- and two-joint muscles (a Hill-type model). Our simulations allow us to suggest a specific network architecture in the spinal cord and a pattern of feedback connectivities (from Ia and Ib fibers and touch s
- **Geyer 2003** — [Positive force feedback in bouncing gaits?](https://doi.org/10.1098/rspb.2003.2454)  
  - animals: Cat · pathways: Ib stance to swing  
  - Bring more context on how the transition from gaits to running are formed and give models examples of how operations roll out.
- **Dietz 2003** — [Spinal Cord Pattern Generators for Locomotion](https://doi.org/10.1016/s1388-2457(03)00120-2)  
  - animals: Human · pathways: Ib stance to swing  
  - CPG generates rhythm an shapes the pattern of bursts of motoneurons(Pg.2)

Cats with transected spinal cord and with cut dorsal roots still produce rhythmic alternating contraction in ankle flexors and extensors (Pg.2)

Basic mechanism underlying locomotion, no fundamental differences seems to exist between bipeds and quadrupeds.(Pg.2)

Differences between cats and primates may be due to increased importance of the Cortico spinal tract in primates. Additionally the gait of primates relies more on Supraspinal drive.(Pg.2)

Spinal Circuitry for locomotion might be suppressed by Supraspinal input
- **Dietz and Duysens 2000** — [Significance of load receptor input during locomotion: a review](https://doi.org/10.1016/s0966-6362(99)00052-1)  
  - animals: Human, Cat · pathways: Ib stance to swing  
  - Review making the case that extensor load receptors are central to locomotor control: in the cat, Golgi tendon organ input switches function during walking from Ib inhibition to extensor facilitation proportional to load, and in humans leg extensor activation during stance scales with body weight — one load-regulatory mechanism spanning quadrupedal and bipedal gait.
- **Pearson 1998** — [Enhancement and Resetting of Locomotor Activity by Muscle Afferentsa](https://doi.org/10.1111/j.1749-6632.1998.tb09050.x)  
  - animals: Mammals · afferents: Ib · pathways: Ib stance to swing; Ib disynaptic excitation; Ib excitatory  
  - Feedback from ankle extensor group I afferents sets two stance-critical features of the mammalian locomotor pattern: Ib Golgi tendon organ input excites the extensor half-center to prolong extensor bursts, so swing is not initiated while the leg is loaded, and extensor burst magnitude is enhanced through disynaptic and polysynaptic pathways opened only during locomotion. Reflex gains must be calibrated to limb biomechanics, and weakening ankle extensors evokes compensatory recalibration of proprioceptive influences on the pattern.  
  - *Robot/sim:* Implement stance-gated Ib excitation of the extensor half-center plus locomotion-gated disynaptic enhancement of extensor motoneuron gain; ablating the load gate should trigger swing under load and abolish load-proportional burst scaling.
- **Prochazka 1997** — [Positive Force Feedback Control of Muscles](https://doi.org/10.1152/jn.1997.77.6.3226)  
  - animals: Human · pathways: Ib stance to swing; Ib excitatory  
  - Establish some basic properties of positive force feedback in relation to load compensation, stability intrinsic muscle properties an interaction with displacement feedback
Positive force feedback was an effective means of generating load compensation
Suggestion a useful role of Ib-mediated positive force feedback from extensor muscles during gait.
- **Whelan 1996** — [CONTROL OF LOCOMOTION IN THE DECEREBRATE CAT](https://doi.org/10.1016/0301-0082(96)00028-7)  
  - animals: Cat · pathways: Ib stance to swing; Ia or II stance to swing  
  - In the decerebrate cat, locomotion is initiated by the mesencephalic locomotor region acting through the medial medullary reticular formation and the ventrolateral funiculus, and the phase transitions are set by identified afferents: group I Golgi tendon organ afferents prolong stance, while length- and velocity-sensitive afferents from extensor muscles signal leg extension to permit swing — unloaded and extended is the swing-permitting state.
- **Pearson 1995** — [Proprioceptive regulation of locomotion.](https://doi.org/10.1016/0959-4388(95)80107-3)  
  - animals: Cat, Insects, Arthropods · pathways: Ib stance to swing; Ia or II stance to swing  
  - Review of proprioceptive regulation across walking systems: during locomotion, Golgi tendon organ feedback from extensor muscles reverses from inhibition to excitation, maintaining stance while the extensors are loaded and helping time the stance-to-swing transition, while primary and secondary spindle afferents also influence the timing of the rhythm — establishing phase-dependent reflex reversal as a general feature of cat and arthropod walking.
- **Gossard 1994** — [Transmission in a locomotor-related group Ib pathway from hindlimb extensor muscles in the cat](https://doi.org/10.1007/bf00228410)  
  - animals: Cat · pathways: Ib stance to swing  
  - Phasic stimulations of group I afferents from ankle and knee extensor muscles may entrain the intrinsic locomotor rhythm.
 
The intrinsic locomotor rhythm then acts on motoneurons through the spinal rhythm generators , and concluded that the major part of these effects come from Golgi tendon organ Ib afferents

Injection of nialamide and L-DOPA evoked long lasting reflexes upon stimulations of high threshold afferents before spontaneous fictive locomotion commenced.

Interneurons must therefore be located in the L7 and S1 spinal segments.
- **Pearson 1992** — [Entrainment of the locomotor rhythm by group Ib afferents from ankle extensor muscles in spinal cats.](https://doi.org/10.1007/bf00230939)  
  - animals: Cat · pathways: Ib stance to swing  
  - Rhythmic contractions of ankle extensor muscles or stimulation of their group I afferents entrain the spinal locomotor rhythm in spinal cats, with flexor bursts timed to follow release of extensor stretch by about 200 ms — direct evidence that Ib input during stance inhibits flexor-burst generation and promotes extensor activity, so a decline in Ib activity near the end of stance helps time the stance-to-swing transition.
- **Conway 1987** — [Proprioceptive input resets central locomotor rhythm in the spinal cat.](https://doi.org/10.1007/bf00249807)  
  - animals: Cat · pathways: Ib stance to swing  
  - Extensor group I afferents have access to central rhythm generators
 >Important in reflex regulation of stepping.


Rythm generators came mainly from Golgi Tendon Organ Ib afferents for group I

Increased load of limb extensor during the stance phase enhance and prolong extensor activity while simultaneously delaying the transition to the swing phase of the step cycle.
- **Duysens and Pearson 1980** — [Inhibition of flexor burst generation by loading ankle extensor muscles in walking cats](https://doi.org/10.1016/0006-8993(80)90206-1)  
  - animals: Cat · pathways: Ib stance to swing  
  - It has been found that decreased load detected by Golgi organs encourages swing muscles to activate
