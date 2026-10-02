# Animal: Arthropods

11 papers in the corpus.

- **Nourse 2019** — [Analyzing the Interplay Between Local CPG Activity and Sensory Signals for Inter-leg Coordination in Drosophil](https://doi.org/10.1007/978-3-030-24741-6_34)  
  - animals: Insects, Arthropods  
  - Tibia-amputated Drosophila legs show a speed-dependent stump oscillation — the wide range of possible periods at low walking speed collapses to a minimum as speed increases — which noisy leg CPGs stabilized by intra-leg load feedback and inter-leg coordinating signals can explain; the measured data anchor a simplified neuromuscular model of inter-leg coordination.
- **Szczecinski and Quinn 2018** — [Leg-local neural mechanisms for searching and learning enhance robotic locomotion](https://doi.org/10.1007/s00422-017-0726-x)  
  - animals: Arthropods  
  - Two leg-local arthropod mechanisms were ported to MantisBot: negative feedback of leg depressor force keeps each stance leg carrying an appropriate share of body weight, and absence of load triggers a ground-search mode - implemented with leg-local memory and command neurons so the robot transitions between searching and stepping while mimicking animal data from the literature.  
  - *Robot/sim:* Implement leg-local negative force feedback onto stance depressor activation plus a load-gated search mode with per-leg memory in a hexapod controller; ablate the force feedback to show body-weight distribution across stance legs fails and searching is never triggered.
- **Edwards and Prilutsky 2017** — [Sensory Feedback in the Control of Posture and Locomotion](https://doi.org/10.1002/9781118873397.ch9)  
  - animals: Arthropods, Vertebrates  
  - Book chapter synthesizing sensory control of posture and locomotion in arthropods and vertebrates through classical feedback-control theory, arguing that experimental paradigms, reduced animal preparations, and neuromechanical modeling must be combined to understand postural control, which in both groups maintains a body-configuration set point (a set of joint angles keeping the body above the feet) during standing, locomotion, and other movements, and resists sudden perturbations.  
  - *Robot/sim:* Build posture controllers for simulated walkers per the chapter's framework - body-configuration set points defended by sensory feedback loops - and ablate loops individually to compare predicted versus resulting postural deficits in arthropod- and vertebrate-style morphologies.
- **Ekeberg 2004** — [Dynamic simulation of insect walking.](https://doi.org/10.1016/j.asd.2004.05.002)  
  - animals: Stick Insect, Insects, Arthropods · pathways: Chordotonal organ multi-synaptic excitation; Chordotonal organ multi-synaptic inhibition  
  - A 3-D biomechanical stick-insect leg driven by a reduced neural controller containing only experimentally established chordotonal-organ reflex mechanisms reproduces the full middle-leg step cycle in both restricted and unrestrained simulations; front-leg stepping works with the same mechanisms, whereas hind-leg control requires reorganized, possibly sign-reversed, chordotonal influence on levator–depressor timing.
- **Koditschek 2004** — [Mechanical aspects of legged locomotion control](https://doi.org/10.1016/j.asd.2004.06.003)  
  - animals: Arthropods · pathways: Biomechanically mediated preflexive feedback  
  - Review enlisting neurophysiology, biomechanics, control engineering, and nonlinear dynamics to explain effective locomotion as the integration of muscular, skeletal, and neural mechanics, using rapid arthropod terrestrial locomotion as the data-rich model system; the outlined control hypotheses and their mathematical underpinnings inspired the design of the hexapedal robot RHex.  
  - *Robot/sim:* Implement the neuromechanical hypotheses in a RHex-style hexapedal simulation with compliant legs, then ablate mechanical (compliance, damping) versus neural (feedback gain) elements to quantify each hypothesis's contribution to stable rapid locomotion.
- **Ritzmann 2004** — [Convergent evolution and locomotion through complex terrain by insects, vertebrates and robots](https://doi.org/10.1016/j.asd.2004.05.001)  
  - animals: Insects, Vertebrates, Arthropods  
  - Argues from convergent evolution that insects and vertebrates independently arrived at similar neural control properties and mechanical schemes for legged locomotion, highlighting leg specialization, body flexion, and complex head-mounted sensing as critical for complex-terrain agility - features most robots of the era lacked. The design-transfer principle is selection rather than copying: pick the properties critical for the target behavior when building a robot.  
  - *Robot/sim:* Build legged-robot models that incrementally add leg specialization, body flexion, and head-region sensing over a homogeneous-leg rigid-body baseline, and compare complex-terrain agility to test which convergent properties carry the performance benefit.
- **Marder and Bucher 2001** — [Central pattern generators and the control of rhythmic movements](https://doi.org/10.1016/s0960-9822(01)00581-4)  
  - animals: Vertebrates, Insects, Arthropods  
  - Review of central pattern generators: rhythmic motor output can be produced with no sensory or descending timing cues, and neuromodulation reconfigures the same circuit to generate multiple patterns — the framework in which one spinal network can produce many gaits, with general principles drawn from invertebrate circuits and vertebrate spinal cord and brainstem.
- **Pearson 1995** — [Proprioceptive regulation of locomotion.](https://doi.org/10.1016/0959-4388(95)80107-3)  
  - animals: Cat, Insects, Arthropods · pathways: Ib stance to swing; Ia or II stance to swing  
  - Review of proprioceptive regulation across walking systems: during locomotion, Golgi tendon organ feedback from extensor muscles reverses from inhibition to excitation, maintaining stance while the extensors are loaded and helping time the stance-to-swing transition, while primary and secondary spindle afferents also influence the timing of the rhythm — establishing phase-dependent reflex reversal as a general feature of cat and arthropod walking.
- **Rossi-Durand 1993** — [Peripheral proprioceptive modulation in crayfish walking leg by serotonin](https://doi.org/10.1016/0006-8993(93)91131-b)  
  - animals: Arthropods · afferents: Mechanosensory  
  - Serotonin acts directly and peripherally on crayfish coxo-basal chordotonal organ afferent terminals, dose-dependently boosting primary afferent discharge (maximal near 1 uM, inhibitory at high concentrations in some units), providing a peripheral neuromodulatory route for gain-setting proprioceptive feedback during locomotion without altering central circuitry.  
  - *Robot/sim:* Add a state-dependent gain multiplier on proprioceptive afferent channels (a serotonin-analog neuromodulation knob ranging from facilitation to suppression) and test how reflex gain and locomotor pattern organization change as it is swept.
- **Elson 1992** — [Identified proprioceptive afferents and motor rhythm entrainment in the crayfish walking system](https://doi.org/10.1152/jn.1992.67.3.530)  
  - animals: Arthropods · afferents: Mechanosensory  
  - In the crayfish walking leg, the thoraco-coxal muscle receptor organ's two graded nonspiking afferents, phasic T and tonic S, mediate opposing feedback pathways that alternately excite opposite burst phases, so proprioceptor stretch entrains the central rhythm through complete, probabilistic resetting of burst transitions within one or two stimulus cycles.  
  - *Robot/sim:* Build a two-phase oscillator resettable by graded phasic and tonic afferent channels with opposing effects on burst transitions; test entrainment ranges, coordination ratios (2:1, 1:2, 1:3), and resetting probability as in the crayfish preparation.
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
