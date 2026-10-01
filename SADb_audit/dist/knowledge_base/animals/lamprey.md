# Animal: Lamprey

16 papers in the corpus.

- **Grillner and El Manira 2020** — [Current Principles of Motor Control, with Special Reference to Vertebrate Locomotion](https://doi.org/10.1152/physrev.00015.2019)  
  - animals: Vertebrates, Mammals, Lamprey  
  - Frames vertebrate locomotor control as spinal CPG microcircuits setting muscle timing with sensory compensation for perturbations, brainstem command systems setting CPG activity level and locomotor speed, basal ganglia selecting and initiating/stopping motor programs, with postural and steering systems integrated around the propulsive core.  
  - *Robot/sim:* Architect a controller with distinct CPG modules, brainstem-like drive/gain modules for speed, and a basal ganglia-like action-selection layer that starts and stops locomotion; ablate each level to reproduce the framework's predicted functional dissociations.
- **Guertin 2013** — [Central Pattern Generator for Locomotion: Anatomical, Physiological, and Pathophysiological Considerations](https://doi.org/10.3389/fneur.2012.00183)  
  - animals: Vertebrates, Lamprey, Human · pathways: Fictive locomotion without sensory feedback  
  - Century-spanning review concluding that walking, flying, and swimming are largely controlled by spinal central pattern generator networks, demonstrated across vertebrate species from lamprey to human, which can self-produce basic rhythmic coordinated movements even in the absence of descending or peripheral inputs; CPG plasticity is then linked to Restless Legs Syndrome, Periodic Leg Movement, Alternating Leg Muscle Activation, and Uner Tan Syndrome.  
  - *Robot/sim:* Model excitability or plasticity changes in CPG elements and test whether they generate spontaneous locomotor-like bursts resembling the reviewed disorders, identifying which parameter changes reproduce each pathology.
- **Kiehn 2011** — [Development and functional organization of spinal locomotor circuits](https://doi.org/10.1016/j.conb.2010.09.004)  
  - animals: Vertebrates, Mice, Zebrafish, Lamprey, Mammals  
  - This review outlines and compares recent
advances that have revealed the developmental and functional
organization of these fundamental spinal motor networks in
limbed and non-limbed animals. The comparison will highlight
common principles and divergence in the organization of the
spinal locomotor network structure in these different species as
well as point to unresolved issues regarding the assembly and
functioning of these networks.
- **Chiel 2009** — [The Brain in Its Body: Motor Control and Sensing in a Biomechanical Context](https://doi.org/10.1523/jneurosci.3338-09.2009)  
  - animals: Cat, Human, Insects, Lamprey, Rat, Salamander · pathways: Biomechanically mediated preflexive feedback  
  - Reviews molluscan feeding, postural control in cats and humans, locomotion simulations in lamprey, insect, cat and salamander, and rat vibrissal sensing to argue that adaptive behavior emerges from nervous-system-body-environment interaction: control is shared between nervous system and periphery, neural activity organizes degrees of freedom into biomechanically meaningful subsets, mechanics alone can play crucial roles in enforcing gait patterns, and the mechanics of sensors is crucial for their function.  
  - *Robot/sim:* Embed morphologically realistic muscle and sensor mechanics so that body dynamics contribute to gait enforcement (preflexes); progressively remove neural correction loops and quantify the locomotor stability retained by mechanics alone.
- **Dubuc 2008** — [Initiation of locomotion in lampreys.](https://doi.org/10.1016/j.brainresrev.2007.07.016)  
  - animals: Lamprey  
  - Review of how locomotion is initiated in lampreys: escape swimming is triggered by a simple cutaneous pathway through reticulospinal command cells whose intrinsic membrane properties convert brief sensory input into long-lasting excitation, while goal-directed locomotion is channeled through the mesencephalic locomotor region, which drives reticulospinal cells bilaterally in graded fashion via glutamatergic and cholinergic transmission.
- **Deliagina 2008** — [Spinal and supraspinal postural networks](https://doi.org/10.1016/j.brainresrev.2007.06.017)  
  - animals: Cat, Human, Lamprey  
  - In the lamprey, the postural control system is driven by vestibular input.
- **Büschges 2005** — [Sensory control and organization of neural networks mediating coordination of multisegmental organs for locomo](https://doi.org/10.1152/jn.00615.2004)  
  - animals: Stick Insect, Cat, Lamprey  
  - Review comparing stick insect, cat, and lamprey: locomotor patterns arise from central pattern generating networks, local sensory feedback about movements and forces in the locomotor organs, and coordinating signals from neighboring segments or appendages; the central network controlling a multi-segmental organ comprises multiple segmental CPGs matching the organ's structure, and common schemes of sensory feedback operate across walking species.  
  - *Robot/sim:* Structure the controller as one CPG module per joint or segment with local force and movement feedback plus inter-module coordinating connections; ablate the local feedback to test the multi-CPG organization hypothesis.
- **Ijspeert 2001** — [A connectionist central pattern generator for the aquatic and terrestrial gaits of a simulated salamander](https://doi.org/10.1007/s004220000211)  
  - animals: Salamander, Lamprey  
  - A connectionist CPG model of leaky-integrator neurons tuned by a genetic algorithm reproduces the salamander's aquatic traveling-wave swimming and terrestrial standing-wave trotting by coupling a lamprey-like body CPG with a limb CPG. Tonic, non-oscillating excitation applied to the single network modulates speed, direction, and gait type, demonstrating that one distributed circuit with two coupled subnetworks can underlie both locomotor programs.  
  - *Robot/sim:* Implement the two-level body-plus-limb CPG with a genetic-algorithm parameter search; sweep the tonic drive to verify gait, speed, and direction modulation, and ablate the limb-to-body CPG coupling to isolate the swimming program.
- **Leeuwen 1999** — [Simulations of neuromuscular control in lamprey swimming](https://doi.org/10.1098/rstb.1999.0441)  
  - animals: Lamprey  
  - Cited by Geyer and Herr.
The neuronal generation of vertebrate locomotion has been extensively studied in the lamprey. Models at different levels of abstraction are being used to describe this system, from abstract nonlinear oscillators to interconnected model neurons comprising multiple compartments and a Hodgkin–Huxley representation of the most relevant ion channels. To study the role of sensory feedback by simulation, it eventually also becomes necessary to incorporate the mechanical movements in the models. By using simplifying models of muscle activation, body mechanics, counteracting wa
- **Ijspeert 1999** — [Evolution and Development of a Central Pattern Generator for the Swimming of a Lamprey](https://doi.org/10.1162/106454699568773)  
  - animals: Lamprey  
  - Computational model study. Inspired by the central pattern generators found in animals, we develop neural controllers that can produce the patterns of oscillations necessary for the swimming of a simulated lamprey.
- **Zehr and Stein 1999** — [What functions do reflexes serve during human locomotion](https://doi.org/10.1016/s0301-0082(98)00081-1)  
  - animals: Lamprey, Cat, Human · afferents: Cutaneous · pathways: Cutaneous flexor excitation; Ia monosynaptic excitation  
  - Review establishing that reflexes during locomotion are task-, phase-, and context-dependent, with a division of labor: cutaneous reflexes alter swing-limb trajectory to avoid stumbling, stretch reflexes stabilize limb trajectory and assist force production during stance, and load receptor reflexes support body weight and influence step-cycle timing - functions dynamically reassigned across the cycle and clinically exploitable after neurotrauma.  
  - *Robot/sim:* Implement phase- and task-gated cutaneous (swing trajectory), Ia stretch (stance stabilization and force), and load-receptor (weight support, cycle timing) pathways in a walker; ablate each to reproduce the predicted functional losses.
- **Grillner 1998** — [Intrinsic function of a neuronal network — a vertebrate central pattern generator](https://doi.org/10.1016/s0165-0173(98)00002-2)  
  - animals: Lamprey  
  - Computational model study. The role of different subtypes of Ca 2q channels, Ca2q dependent Kq channels and voltage dependent NMDA channels at the neuronal and network level is in focus as well as the effects of different metabotropic, aminergic and peptidergic modulators that target these ion channels.
- **Grillner 1995** — [Neural networks that co-ordinate locomotion and body orientation in lamprey](https://doi.org/10.1016/0166-2236(95)80008-p)  
  - animals: Lamprey  
  - Lamprey networks coordinating locomotion and body orientation; the lamprey as the reference preparation for cellular-level CPG and reticulospinal control.
- **Rossignol and Dubuc 1994** — [Spinal pattern generation](https://doi.org/10.1016/0959-4388(94)90139-2)  
  - animals: Rat, Lamprey, Vertebrates  
  - Review framing spinal pattern generation around three foci: transmitter actions on rhythms in reduced preparations, activity-dependent changes in membrane properties of the generating circuits, and CPG-afferent interactions, emphasizing that new reflex responses and membrane properties emerge only during rhythmic activity.  
  - *Robot/sim:* Implement CPG circuits whose reflex pathway gains and signs are gated by rhythmic state, and test whether reflex responses absent in quiescent simulation appear when the CPG is active, mirroring state-dependent reflex modification.
- **Cohen 1992** — [Modelling of intersegmental coordination in the lamprey central pattern generator for locomotion](https://doi.org/10.1016/0166-2236(92)90006-t)  
  - animals: Lamprey  
  - Review of the lamprey intersegmental coordination program combining mathematical analysis, biological experimentation, and computer simulation: swimming holds approximately one body wavelength at all speeds because ascending intersegmental coupling dominates descending coupling and thereby sets the intersegmental phase lag, with long-range coupling adding significant functional influence.  
  - *Robot/sim:* Build a chain of segmental oscillators with asymmetric coupling (dominant ascending, including long-range connections) and verify constant body-wavelength phase relationships across swimming frequencies; ablating ascending coupling should break the phase lag.
- **Grillner 1985** — [Neurobiological Bases of Rhythmic Motor Acts in Vertebrates](https://doi.org/10.1126/science.3975635)  
  - animals: Vertebrates, Mammals, Lamprey  
  - Classic formulation of the principles underlying nervous control of innate rhythmic motor acts in vertebrates, establishing the in vitro lamprey spinal cord preparation in which the locomotor motor pattern can be elicited in isolated cord sections for detailed mechanism analysis, and proposing that parts of the locomotor network are reused in other motor acts, including learned ones.  
  - *Robot/sim:* Implement modular CPG subnetworks shared across multiple motor programs in a simulated swimmer or walker and verify that reuse preserves timing and stability while reducing controller complexity.
