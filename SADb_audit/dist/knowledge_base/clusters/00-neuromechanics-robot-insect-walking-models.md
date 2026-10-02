# Cluster 1: Neuromechanics & robot/insect walking models

205 papers · years 1972–2026 · 181 curated

## Top authors

- Szczecinski (6)
- Büschges (4)
- Beer (4)
- Hunt (4)
- Young (4)
- Schmitt (3)
- Hilts (3)
- Deng (3)
- Ijspeert (2)
- Burkitt (2)

## Pathways represented

- Biomechanically mediated preflexive feedback (10)
- Chordotonal organ multi-synaptic excitation (5)
- Chordotonal organ multi-synaptic inhibition (5)
- Ib stance to swing (4)
- trochanteral campaniform sensilla load signals to adjust MN magnitude (3)
- Ib disynaptic excitation (2)
- Ia or II stance to swing (1)
- Ib swing to stance (1)
- Ia stance to swing (1)
- Ia swing to stance (1)
- Ia monosynaptic (1)
- Ia monosynaptic excitation (1)

## Key papers (by citations)

- **Wilson 1972** — [Excitatory and Inhibitory Interactions in Localized Populations of Model Neurons](https://doi.org/10.1016/s0006-3495(72)86068-5)  
  - The Wilson-Cowan coupled nonlinear equations for localized excitatory and inhibitory neuron populations exhibit simple and multiple hysteresis, multiple stable states, and limit-cycle oscillations whose frequency is a monotonic function of stimulus intensity, with a formal result linking oscillation under one stimulus class to multistability under another - the foundational dynamical substrate for neural rhythm-generation models.  
  - *Robot/sim:* Use coupled excitatory/inhibitory population units as CPG building blocks, exploiting their monotonic frequency-versus-stimulus relation for drive-dependent cycle frequency and their hysteresis and multistability for gait-state switching and regime transitions.
- **Horak 2006** — [Postural orientation and equilibrium: what do we need to know about neural control of balance to prevent falls](https://doi.org/10.1093/ageing/afl077)  
  - animals: Human  
  - Frames postural control as two interacting goals - orientation (active alignment of trunk and head to gravity, support surface, and visual surround) and equilibrium (movement strategies stabilizing the center of body mass) - achieved by dynamic sensorimotor processes in which the weighting of somatosensory, vestibular, and visual inputs depends on task goals and context; anticipatory postural adjustments precede voluntary limb movement, and damage to different underlying systems yields context-specific instabilities.  
  - *Robot/sim:* Implement a balance controller with task- and context-dependent weighting of somatosensory, vestibular, and visual channels plus anticipatory postural adjustments before limb movements; ablating single channels should reproduce distinct, context-specific instability signatures.
- **Ijspeert 2008** — [2008 Special Issue: Central pattern generators for locomotion control in animals and robots: A review](https://doi.org/10.1016/j.neunet.2008.03.014)  
  - animals: Vertebrates  
  - The canonical review of locomotor CPGs spanning neurobiology and robotics: CPG circuits turn simple low-dimensional inputs into coordinated high-dimensional rhythmic output, CPG models — neural networks or coupled oscillators — are a standard robot-control approach, and robots themselves serve as scientific tools for testing biological hypotheses about CPG function.
- **Marder and Bucher 2001** — [Central pattern generators and the control of rhythmic movements](https://doi.org/10.1016/s0960-9822(01)00581-4)  
  - animals: Vertebrates, Insects, Arthropods  
  - Review of central pattern generators: rhythmic motor output can be produced with no sensory or descending timing cues, and neuromodulation reconfigures the same circuit to generate multiple patterns — the framework in which one spinal network can produce many gaits, with general principles drawn from invertebrate circuits and vertebrate spinal cord and brainstem.
- **Burkitt 2006** — [A Review of the Integrate-and-fire Neuron Model: I. Homogeneous Synaptic Input](https://doi.org/10.1007/s00422-006-0068-6)  
  - Review of review integrate-and-fire neuron model: i.. Computational model study. The integrate-and-ﬁre neuron model has become established as a canonical model for the description of spiking neurons be- cause it is capable of being analyzed mathematically while at the same time being sufﬁciently complex to capture many of the essential features of neural processing.
- **Taga 1991** — [Self-organized control of bipedal locomotion by neural oscillators in unpredictable environment](https://doi.org/10.1007/bf00198086)  
  - Proposed that stable and flexible locomotion emerges as a global limit cycle produced by global entrainment between coupled neural oscillators, the musculoskeletal system, and the environment, rather than by tracking explicit trajectories; the simulated biped withstood mechanical perturbations and environmental changes, and switched between walking and running (with hysteresis) by changing a single nonspecific parameter.  
  - *Robot/sim:* Implement coupled neural oscillators with mutual inhibition, entrained with a musculoskeletal biped through sensorimotor feedback loops; verify perturbation resistance, environmental adaptation, and walk-run transitions with hysteresis under a single-parameter drive change, replicating the paper's claims.
- **Chiel and Beer 1997** — [The brain has a body: adaptive behavior emerges from interactions of nervous system, body and environment.](https://doi.org/10.1016/s0166-2236(97)01149-1)  
  - The brain has a body: adaptive behavior emerges from nervous system–body–environment interaction, so neural analysis alone cannot explain motor control.
- **Geyer and Herr 2010** — [A Muscle-Reflex Model That Encodes Principles of Legged Mechanics Produces Human Walking Dynamics and Muscle A](https://doi.org/10.1109/TNSRE.2010.2047592)  
  - animals: Human  
  - While neuroscientists identify increasingly complex neural circuits that control animal and human gait, biomechanists ﬁnd that locomotion requires little control if principles of legged mechanics are heeded that shape and exploit the dynamics of legged systems. Here, we show that muscle reﬂexes could be vital to link these two observations. We develop a model of human locomotion that is controlled by muscle reﬂexes which encode principles of legged mechanics. Equipped with this reﬂex control, we ﬁnd this model to stabilize into a walking gait from its dynamic interplay with the ground, reprodu
- **Pearson and Kg 1993** — [Common Principles of Motor Control in Vertebrates and Invertebrates](https://doi.org/10.1146/annurev.ne.16.030193.001405)  
  - The abstract attached to this record describes the brain's default mode network rather than motor control, an apparent database mismatch, so no content can be grounded from it and the content fields are left empty; the record needs its abstract re-sourced from the DOI before curation.
- **Taga 1995** — [A model of the neuro-musculo-skeletal system for human locomotion](https://doi.org/10.1007/bf00204048)  
  - animals: Human  
  - A neuromusculoskeletal model in which human gait emerges as a stable limit cycle through global entrainment among neural oscillators, the musculoskeletal body, and the environment: walking is robust to perturbations and loads, graded in speed by a single tonic drive parameter, and entrainable by rhythmic input — the CPG and the body need not be programmed separately.
- **Vogels 2005** — [NEURAL NETWORK DYNAMICS](https://doi.org/10.1146/annurev.neuro.28.061604.135637)  
  - Review of neural network dynamics. Computational model study. We also review propagation of stimulus-driven activity through sponta- neously active networks.
- **Delp 1995** — [A graphics-based software system to develop and analyze models of musculoskeletal structures](https://doi.org/10.1016/0010-4825(95)98882-e)  
  - Computational model study. -We have created a graphics-based software system that enables users to develop and analyze musculoskeletal models without programming.
- **Orlovsky 1999** — [Neuronal Control of LocomotionFrom Mollusc to Man](https://doi.org/10.1093/acprof:oso/9780198524052.001.0001)  
  - animals: Human  
  - Comparative monograph surveying the neural mechanisms of locomotion across evolutionarily diverse species, from leech swimming to human running, synthesizing CPG organization, descending control, and sensory regulation into unifying principles of how nervous systems generate locomotion; a canonical cross-species reference for locomotor neurobiology.  
  - *Robot/sim:* Use its comparative survey of CPG organization and descending control across species to select conserved motifs (half-center rhythm cores, brainstem initiation and speed control) as architecture choices for robotic locomotor controllers.
- **Cruse 1990** — [What mechanisms coordinate leg movement in walking arthropods?](https://doi.org/10.1016/0166-2236(90)90057-H)  
  - animals: Arthropods  
  - Comparing the walk of a six-legged robot with the walk of an insect, the immense differences are immediately obvious. The walking of an
animal is much more versatile, and seems to be more
effective and elegant. Thus it is useful to consider the
corresponding biological mechanisms in order to apply
these or similar mechanisms to the control of walking
legs in machines. Until recently, little information on
the biological control mechanisms has been available;
this paper summarizes recent developments.
- **Grillner 1995** — [Neural networks that co-ordinate locomotion and body orientation in lamprey](https://doi.org/10.1016/0166-2236(95)80008-p)  
  - animals: Lamprey  
  - Lamprey networks coordinating locomotion and body orientation; the lamprey as the reference preparation for cellular-level CPG and reticulospinal control.
- **Gerstner 2000** — [Population Dynamics of Spiking Neurons: Fast Transients, Asynchronous States, and Locking](https://doi.org/10.1162/089976600300015899)  
  - It is analytically sho wn that transien ts from a state of incoheren t / ring can b e immediate/.
- **Ijspeert 2001** — [A connectionist central pattern generator for the aquatic and terrestrial gaits of a simulated salamander](https://doi.org/10.1007/s004220000211)  
  - animals: Salamander, Lamprey  
  - A connectionist CPG model of leaky-integrator neurons tuned by a genetic algorithm reproduces the salamander's aquatic traveling-wave swimming and terrestrial standing-wave trotting by coupling a lamprey-like body CPG with a limb CPG. Tonic, non-oscillating excitation applied to the single network modulates speed, direction, and gait type, demonstrating that one distributed circuit with two coupled subnetworks can underlie both locomotor programs.  
  - *Robot/sim:* Implement the two-level body-plus-limb CPG with a genetic-algorithm parameter search; sweep the tonic drive to verify gait, speed, and direction modulation, and ablate the limb-to-body CPG coupling to isolate the swimming program.
- **Nishikawa et al. 2007** — [Neuromechanics: an integrative approach for understanding motor control.](https://doi.org/10.1093/icb/icm024)  
  - Manifesto for neuromechanics: understanding motor control requires the integrated nervous-system-plus-body as the unit of analysis; the intellectual case with examples across locomotion and feeding.
- **Golowasch 2002** — [Failure of Averaging in the Construction of a Conductance-Based Neuron Model](https://doi.org/10.1152/jn.00412.2001)  
  - A conductance-based model built from the mean maximal conductances of a population of equivalent one-spike bursters fires three action potentials per burst instead of one, because the population occupies a highly concave region of parameter space that does not contain its mean - demonstrating that averaging parameters across multiple preparations can fail to characterize a system whose behavior depends on interactions among highly variable components. A core caution for parameterizing neural models from pooled experimental data.  
  - *Robot/sim:* When fitting neural or CPG models to multiple experimental recordings, select parameter sets by matching behavior rather than by averaging parameters, and always verify that the mean-parameter model reproduces the target dynamics - an explicit modeling-hygiene rule for database-driven model construction.
- **Burkitt 2006** — [A review of the integrate-and-fire neuron model: II. Inhomogeneous synaptic input and network properties](https://doi.org/10.1007/s00422-006-0082-8)  
  - Review of review integrate-and-fire neuron model: ii.. Computational model study. The integrate-and-ﬁre neuron model has become established as a canonical model for the description of spiking neurons be- cause it is capable of being analyzed mathematically while at the same time being sufﬁciently complex to capture many of the essential features of neural processing.
- **Tuthill and Wilson 2016** — [Mechanosensation and Adaptive Motor Control in Insects](https://doi.org/10.1016/j.cub.2016.06.070)  
  - animals: Insects  
  - The ability of animals to flexibly navigate through complex environments depends on the integration of sensory information with motor commands. The sensory modality most tightly linked to motor control is mechanosensation. Adaptive motor control depends critically on an animal’s ability to respond to mechanical forces
generated both within and outside the body. The compact neural circuits of insects provide appealing systems to investigate how mechanical cues guide locomotion in rugged environments. Here, we review our current understanding of mechanosensation in insects and its role in adapti
- **Yu 2014** — [A Survey on CPG-Inspired Control Models and System Implementation](https://doi.org/10.1109/tnnls.2013.2280596)  
  - Survey of two decades of CPG-inspired robot locomotion control across abstraction levels: CPGs modeled as coupled neuron groups generate rhythmic signals without sensory feedback but require sensory feedback to shape those signals; reviews the relative advantages of the models, the main design, optimization, and implementation issues, and trends for multi-DOF articulated control.  
  - *Robot/sim:* Use as a comparative menu of CPG abstractions for a simulated walker: implement several abstraction levels on a common plant and evaluate the survey's design criteria (ease of design, optimization, implementation) alongside sensory-shaping variants.
- **Farina 2014** — [The effective neural drive to muscles is the common synaptic input to motor neurons](https://doi.org/10.1113/jphysiol.2014.273581)  
  - Computational model study. It is shown theoretically that for frequencies smaller than the average discharge rates of the motor neurons, the pool of motor neurons determines a pure ampliﬁcation of the frequency components common to all motor neurons, so that the common input is transmitted almost undistorted and...
- **Mendes 2013** — [Quantification of gait parameters in freely walking wild type and sensory deprived Drosophila melanogaster](https://doi.org/10.7554/elife.00231)  
  - animals: Insects  
  - High-speed optical tracking of freely walking Drosophila quantified wild-type gait, tarsal positioning, and intersegmental and left-right coordination; genetic inactivation of leg sensory neurons (blocking proprioceptive feedback) left interleg coordination and the tripod gait intact but degraded step precision - central coupling suffices for interleg coordination while leg proprioception tunes the accuracy of individual steps.  
  - *Robot/sim:* In a hexapod model, ablate leg-local proprioceptive feedback while keeping central interleg coupling - the phenotype to reproduce is preserved tripod coordination with degraded footfall precision.
- **Hyun 2014** — [High speed trot-running: Implementation of a hierarchical controller using proprioceptive impedance control on](https://doi.org/10.1177/0278364914532150)  
  - Computational model study. The developed controller enables high-speed running of up to 6 m/s (Froude number of Fr ≈ 7.34) incorporating proprioceptive feedback and programmable virtual leg compliance of the MIT Cheetah.
