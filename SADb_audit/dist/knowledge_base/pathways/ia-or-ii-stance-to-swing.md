# Feedback pathway: Ia or II stance to swing

10 papers in the corpus.

- **Akay and Murray 2021** — [Relative Contribution of Proprioceptive and Vestibular Sensory Systems to Locomotion: Opportunities for Discov](https://doi.org/10.3390/ijms22031467)  
  - animals: Cat · pathways: Ia or II stance to swing  
  - Hypothesis: The Weight of the role of segmental feedback is less important at slower speeds but increases at higher speeds, where as in teh vestibular system has the opposite effect

Mouse model is key tool to undertadning terrestial locmotion (Pg 1)

Vestivular feedback being critical at lower speed and somatosensory feedback is necessary at higher velocites. (Pg1)


CPG works with both sensory feedback from the leg or supraspinal center to generate locomotor patter that is able to deal with changes in terrain (Pg2)

The foot is on the ground during the stance phase and moves in the opposite 
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
- **Akay 2014** — [Degradation of mouse locomotor pattern in the absence of proprioceptive sensory feedback](https://doi.org/10.1073/pnas.1419045111)  
  - animals: Mice · pathways: Ia or II swing to stance; Ia or II stance to swing; Ib swing to stance  
  - Ia/ II stance to swing

Ia/II and Ib swing to stance
- **Pearson 2008** — [Role of sensory feedback in the control of stance duration in walking cats](https://doi.org/10.1016/j.brainresrev.2007.06.014)  
  - animals: Cat · pathways: Ib stance to swing; Ia or II stance to swing  
  - Review Paper
- **Stecina et al. 2005** — [Parallel reflex pathways from flexor muscle afferents evoking resetting and flexion enhancement during fictive](https://doi.org/10.1113/jphysiol.2005.095505)  
  - animals: Cat · pathways: Ia or II stance to swing  
  - Parallel reflex pathways from flexor muscle afferents evoke resetting and flexion enhancement during fictive locomotion and scratch; group I/II flexor afferents access the rhythm generator through multiple routes.
- **Hiebert 1996** — [Contribution of Hind Limb Flexor Muscle Afferents to the Timing of Phase Transitions in the Cat Step Cycle](https://doi.org/10.1152/jn.1996.75.3.1126)  
  - animals: Cat · afferents: Ia, II · pathways: Ia or II stance to swing  
  - In walking decerebrate cats, flexor muscle lengthening during stance terminates extensor activity and resets the rhythm to ipsilateral flexion with contralateral extension: EDL and iliopsoas act through group Ia spindle afferents and tibialis anterior through group II afferents, providing a flexor-length gate on the stance-to-swing transition.  
  - *Robot/sim:* Implement flexor-length feedback (Ia from iliopsoas and EDL, II from tibialis anterior) that inhibits the extensor half-center to terminate stance; ablation should delay swing onset and prolong stance.
- **Whelan 1996** — [CONTROL OF LOCOMOTION IN THE DECEREBRATE CAT](https://doi.org/10.1016/0301-0082(96)00028-7)  
  - animals: Cat · pathways: Ib stance to swing; Ia or II stance to swing  
  - In the decerebrate cat, locomotion is initiated by the mesencephalic locomotor region acting through the medial medullary reticular formation and the ventrolateral funiculus, and the phase transitions are set by identified afferents: group I Golgi tendon organ afferents prolong stance, while length- and velocity-sensitive afferents from extensor muscles signal leg extension to permit swing — unloaded and extended is the swing-permitting state.
- **Pearson 1995** — [Proprioceptive regulation of locomotion.](https://doi.org/10.1016/0959-4388(95)80107-3)  
  - animals: Cat, Insects, Arthropods · pathways: Ib stance to swing; Ia or II stance to swing  
  - Review of proprioceptive regulation across walking systems: during locomotion, Golgi tendon organ feedback from extensor muscles reverses from inhibition to excitation, maintaining stance while the extensors are loaded and helping time the stance-to-swing transition, while primary and secondary spindle afferents also influence the timing of the rhythm — establishing phase-dependent reflex reversal as a general feature of cat and arthropod walking.
- **Perreault 1995** — Effects of stimulation of hindlimb flexor group II afferents during fictive locomotion in the cat  
  - animals: Cat · pathways: Ia or II stance to swing; type II excitatory  
  - Excitatory input to flexor type II afferent.
- **Kriellaars 1994** — [Mechanical entrainment of fictive locomotion in the decerebrate cat.](https://doi.org/10.1152/jn.1994.71.6.2074)  
  - animals: Cat · pathways: Ia or II stance to swing  
  - Sinusoidal hip movements as small as 5-20 degrees entrain fictive locomotion in decerebrate cats via low-threshold stretch-sensitive afferents from intrinsic hip muscles (joint capsular afferents unnecessary), with extensor bursts locked to imposed hip flexion as the extensors are stretched and loaded - a positive-feedback mechanism by which the rhythm generator prolongs extension. Subharmonic entrainment and frequency-dependent phase shifts mark the locomotor rhythm generator as a nonlinear oscillator whose rhythm-generating interneurons receive highly convergent afferent input from muscles s  
  - *Robot/sim:* Implement hip stretch afferents as direct inputs to the extensor half-center of the rhythm generator that prolong stance and entrain the cycle to imposed hip motion; ablating them should abolish mechanical entrainment and stance prolongation.
