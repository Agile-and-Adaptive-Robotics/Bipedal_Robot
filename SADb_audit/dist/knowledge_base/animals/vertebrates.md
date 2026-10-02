# Animal: Vertebrates

32 papers in the corpus.

- **Ijspeert and Daley 2023** — [Integration of feedforward and feedback control in the neuromechanics of vertebrate locomotion: a review of ex](https://doi.org/10.1242/jeb.245784)  
  - animals: Vertebrates  
  - Review arguing that vertebrate locomotion emerges from distributed control, spinal CPGs as coupled-oscillator feedforward controllers plus multiple reflex feedback loops, activated and modulated by descending pathways, whose relative contributions vary among animals with body size, intrinsic mechanical stability, time to reach locomotor maturity, and speed. It hypothesizes that distal joints rely more on feedback control than proximal joints, and positions robots and neuromechanical simulations as complements to experiments through what-if scenario testing.  
  - *Robot/sim:* Build hybrid CPG-plus-reflex architectures with joint-dependent feedback weighting (stronger feedback at distal joints) and use what-if simulations to test the proposed division of labor across body size, intrinsic mechanical stability, and speed.
- **Dubuc 2023** — [Locomotor pattern generation and descending control: a historical perspective](https://doi.org/10.1152/jn.00204.2023)  
  - animals: Vertebrates  
  - Historical review of vertebrate locomotor control organizing a century of work into three interacting systems: spinal central pattern generators produce the basic rhythm, brainstem descending inputs initiate, maintain, and stop locomotion while controlling speed and direction, and sensory inputs adapt the locomotor program to environmental conditions.  
  - *Robot/sim:* Structure a locomotor controller as a spinal CPG plus brainstem-analog descending commands (start, stop, speed, direction) plus sensory adaptation loops, and test that removing each layer reproduces its predicted functional loss.
- **Wilson and Sweeney 2023** — [Spinal cords: Symphonies of interneurons across species](https://doi.org/10.3389/fncir.2023.1146449)  
  - animals: Vertebrates  
  - Review of spinal interneuron composition across vertebrates: fish run on two basic classes (ipsilateral excitatory and commissural inhibitory neurons) plus ipsilateral inhibition for escape swimming, and with the evolution of limbs these three types proliferated into molecularly, anatomically, and functionally distinct subpopulations — movement elaboration is mirrored by interneuron specialization.
- **Zholudeva 2021** — [Spinal Interneurons as Gatekeepers to Neuroplasticity after Injury or Disease.](https://doi.org/10.1523/jneurosci.1654-20.2020)  
  - animals: Mice, Cat, Rat, Human, Vertebrates  
  - Review positioning spinal interneurons — a heterogeneous population that modulates motor, sensory, and autonomic function — as key components of plasticity and recovery after spinal cord injury, and arguing that treatments should be optimized for how they engage interneuron circuits rather than viewed as acting around them.
- **Grillner and Kozlov 2021** — [The CPGs for Limbed Locomotion-Facts and Fiction.](https://doi.org/10.3390/ijms22115882)  
  - animals: Mammals, Vertebrates  
  - Review separating fact from fiction in limbed-locomotion CPG claims; takes a hard line on evidence quality for interneuronal and computational assertions.
- **Grillner and El Manira 2020** — [Current Principles of Motor Control, with Special Reference to Vertebrate Locomotion](https://doi.org/10.1152/physrev.00015.2019)  
  - animals: Vertebrates, Mammals, Lamprey  
  - Frames vertebrate locomotor control as spinal CPG microcircuits setting muscle timing with sensory compensation for perturbations, brainstem command systems setting CPG activity level and locomotor speed, basal ganglia selecting and initiating/stopping motor programs, with postural and steering systems integrated around the propulsive core.  
  - *Robot/sim:* Architect a controller with distinct CPG modules, brainstem-like drive/gain modules for speed, and a basal ganglia-like action-selection layer that starts and stops locomotion; ablate each level to reproduce the framework's predicted functional dissociations.
- **Steuer and Guertin 2019** — [Central pattern generators in the brainstem and spinal cord: an overview of basic principles, similarities and](https://doi.org/10.1515/revneuro-2017-0102)  
  - animals: Vertebrates  
  - Comparative overview of brainstem and spinal CPGs (deglutition, mastication, respiration, defecation, micturition, ejaculation, locomotion) showing that organization, function, and cellular properties are generally well-preserved phylogenetically across in vivo models, in vitro preparations, and primates, with the locomotor network the most fully characterized — supporting transfer of locomotor CPG design principles across species and preparations.  
  - *Robot/sim:* Use the review's summary of conserved CPG organization (modular rhythmogenic units, descending trigger and steering inputs) as the architectural template for a simulated network controlling multiple rhythmic behaviors, not only locomotion.
- **Ting and Chiel 2017** — [Muscle, Biomechanics, and Implications for Neural Control](https://doi.org/10.1002/9781118873397.ch12)  
  - animals: Vertebrates, Invertebrates  
  - Argues that muscle multifunctionality versus specialization determines the relative division of labor between neural and muscular control in vertebrates and invertebrates: the transformation from neural activation to muscle force depends on intrinsic neuromuscular properties, body structure, environmental forces, and behavioral context, so muscles are not single-function actuators and neural control complexity must be evaluated within that biomechanical context.  
  - *Robot/sim:* Endow simulated muscles with realistic intrinsic and context-dependent properties, then measure how much neural control complexity is required as muscle multifunctionality varies - directly testing the chapter's neural-versus-muscular complexity trade-off.
- **Edwards and Prilutsky 2017** — [Sensory Feedback in the Control of Posture and Locomotion](https://doi.org/10.1002/9781118873397.ch9)  
  - animals: Arthropods, Vertebrates  
  - Book chapter synthesizing sensory control of posture and locomotion in arthropods and vertebrates through classical feedback-control theory, arguing that experimental paradigms, reduced animal preparations, and neuromechanical modeling must be combined to understand postural control, which in both groups maintains a body-configuration set point (a set of joint angles keeping the body above the feet) during standing, locomotion, and other movements, and resists sudden perturbations.  
  - *Robot/sim:* Build posture controllers for simulated walkers per the chapter's framework - body-configuration set points defended by sensory feedback loops - and ablate loops individually to compare predicted versus resulting postural deficits in arthropod- and vertebrate-style morphologies.
- **Luca 2017** — [Bioinspired morphing wings for extended flight envelope and roll control of small drones](https://doi.org/10.1098/rsfs.2016.0092)  
  - animals: Vertebrates, Human  
  - We show that a fully deployed configuration enhances manoeuvrability while a folded configuration offers low drag at high speeds and is beneficial in strong headwinds.
- **Mohamed 2016** — [Development and Flight Testing of a Turbulence Mitigation System for Micro Air Vehicles](https://doi.org/10.1002/rob.21626)  
  - animals: Vertebrates  
  - Development and Flight Testing of a Turbulence Mitigation System for Micro Air Vehicles •••••••••••••••••••••••••••••••••••• A.
- **Guertin 2013** — [Central Pattern Generator for Locomotion: Anatomical, Physiological, and Pathophysiological Considerations](https://doi.org/10.3389/fneur.2012.00183)  
  - animals: Vertebrates, Lamprey, Human · pathways: Fictive locomotion without sensory feedback  
  - Century-spanning review concluding that walking, flying, and swimming are largely controlled by spinal central pattern generator networks, demonstrated across vertebrate species from lamprey to human, which can self-produce basic rhythmic coordinated movements even in the absence of descending or peripheral inputs; CPG plasticity is then linked to Restless Legs Syndrome, Periodic Leg Movement, Alternating Leg Muscle Activation, and Uner Tan Syndrome.  
  - *Robot/sim:* Model excitability or plasticity changes in CPG elements and test whether they generate spontaneous locomotor-like bursts resembling the reviewed disorders, identifying which parameter changes reproduce each pathology.
- **Fouad 2012** — [Comparative Locomotor Systems](https://doi.org/10.1002/9781118133880.hop203007)  
  - animals: Vertebrates  
  - Comparative synthesis arguing that basic principles of locomotor network organization are conserved between vertebrates and invertebrates, permitting mechanistic study at multiple physiological levels; sensory signals and neuromodulatory input reconfigure network properties, and these networks reorganize in response to nervous-system injury.  
  - *Robot/sim:* Implement CPG networks whose properties are reconfigurable by simulated neuromodulatory and sensory inputs, then reproduce the reported reorganization after simulated injury (partial deafferentation or transection) as a benchmark for recovery mechanisms.
- **Grillner 2011** — [Control of Locomotion in Bipeds, Tetrapods, and Fish](https://doi.org/10.1002/cphy.cp010226)  
  - animals: Vertebrates · afferents: Ia, II, III/IV, Cutaneous  
  - Comprehensive synthesis of locomotor control across bipeds, tetrapods, and fish: spinal unit-burst-generator organization, interlimb coordination, load- and hip-position-dependent reflex regulation of the spinal CPG, tonic gating of reflexes from skin and group II/III afferents, and brainstem and cerebellar control - the canonical reference architecture for locomotor CPG models.  
  - *Robot/sim:* Use the chapter's architecture - unit burst generators per limb/joint, interlimb coordination pathways, and load and hip-position feedback gating the rhythm - as the skeleton for a locomotion controller, with the reviewed reflex gating (skin, group II/III) implemented as state-dependent pathway gains.
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
- **Harris-Warrick 2010** — [General Principles of Rhythmogenesis in Central Pattern Generator Networks](https://doi.org/10.1016/b978-0-444-53613-6.00014-9)  
  - animals: Vertebrates  
  - Four general principles of CPG rhythmogenesis: rhythmogenic ionic currents underlie every CPG regardless of network-pacemaker versus endogenous-kernel organization; fast synaptic transmission can evoke slow currents that alter cycle frequency; rhythmogenesis is multiply and redundantly implemented; and glial cells may participate - with CPGs able to drive rhythmic behavior without sensory feedback while afferents modulate cycle frequency and phase timing.  
  - *Robot/sim:* Build CPG neurons with coexisting rhythmogenic mechanisms (intrinsic currents, network reciprocity, slow synaptic currents) and ablate them one at a time - reproducing why single-mechanism lesions typically fail to abolish rhythm.
- **Grant 2010** — [Design and analysis of biomimetic joints for morphing of micro air vehicles](https://doi.org/10.1088/1748-3182/5/4/045007)  
  - animals: Vertebrates  
  - Review of design analysis biomimetic joints morphing. This paper presents an overview of designs that incorporate morphing to enhance their ﬂight characteristics.
- **Goulding 2009** — [Circuits controlling vertebrate locomotion: moving in a new direction](https://doi.org/10.1038/nrn2608)  
  - animals: Vertebrates  
  - Review arguing that spinal central pattern generator networks are an experimentally tractable model system for understanding how moderately complex neuronal ensembles generate specific motor behaviors, and that novel molecular-genetic tools plus advances in knowledge of spinal cord development bring a comprehensive circuit-level account of locomotor organization within reach.  
  - *Robot/sim:* Build CPG models whose components correspond to genetically identified interneuron classes and delete classes one at a time to compare simulated phenotypes with mutant locomotion.
- **Grillner and Jessell 2009** — [Measured motion: searching for simplicity in spinal locomotor networks.](https://doi.org/10.1016/j.conb.2009.10.011)  
  - animals: Vertebrates  
  - Review of the rules governing spinal locomotor network assembly and function: variations in wiring-diagram organization across vertebrates adapt motor programs to environmental demands; interneuron membrane properties and synaptic interactions underlie modulation of motor circuits and encoded motor behaviors; and molecular-genetic approaches now permit mapping and manipulating interneuron connectivity to link perturbations to network function and motor behavior.  
  - *Robot/sim:* Encode the review's interneuron classes as separable functional modules (rhythm-, pattern-, and output-level) in a locomotor controller and run in silico perturbations mirroring genetic ablations to predict behavioral deficits.
- **Ijspeert 2008** — [2008 Special Issue: Central pattern generators for locomotion control in animals and robots: A review](https://doi.org/10.1016/j.neunet.2008.03.014)  
  - animals: Vertebrates  
  - The canonical review of locomotor CPGs spanning neurobiology and robotics: CPG circuits turn simple low-dimensional inputs into coordinated high-dimensional rhythmic output, CPG models — neural networks or coupled oscillators — are a standard robot-control approach, and robots themselves serve as scientific tools for testing biological hypotheses about CPG function.
- **Lentink 2007** — [How swifts control their glide performance with morphing wings](https://doi.org/10.1038/nature05733)  
  - animals: Vertebrates  
  - Computational model study. Here we describe the aerodynamic and structural performance of actual swift wings, as measured in a wind tunnel, and on this basis build a semi- empirical glide model.
- **Hultborn and Nielsen 2007** — [Spinal control of locomotion--from cat to man.](https://doi.org/10.1111/j.1748-1716.2006.01651.x)  
  - animals: Cat, Human, Vertebrates  
  - Review establishing that spinal networks generate the basic locomotor rhythm across vertebrates including man, with limb sensory feedback essential for effective locomotion: sensory regulation reaches motoneurons via reflex pathways that bypass the rhythm generators and also acts on the locomotor networks themselves, controlling phase timing, shaping muscle-activity patterns, adding excitatory drive, and driving long-term adaptation - the basis for treadmill-training rehabilitation after spinal cord injury.  
  - *Robot/sim:* Implement both sensory routes in a locomotion model - direct reflex pathways to motoneurons plus afferent input to the rhythm generator - and ablate each separately to test their predicted contributions to phase timing, pattern shaping, and net excitatory drive.
- **Ritzmann 2004** — [Convergent evolution and locomotion through complex terrain by insects, vertebrates and robots](https://doi.org/10.1016/j.asd.2004.05.001)  
  - animals: Insects, Vertebrates, Arthropods  
  - Argues from convergent evolution that insects and vertebrates independently arrived at similar neural control properties and mechanical schemes for legged locomotion, highlighting leg specialization, body flexion, and complex head-mounted sensing as critical for complex-terrain agility - features most robots of the era lacked. The design-transfer principle is selection rather than copying: pick the properties critical for the target behavior when building a robot.  
  - *Robot/sim:* Build legged-robot models that incrementally add leg specialization, body flexion, and head-region sensing over a homogeneous-leg rigid-body baseline, and compare complex-terrain agility to test which convergent properties carry the performance benefit.
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
- **Marder and Bucher 2001** — [Central pattern generators and the control of rhythmic movements](https://doi.org/10.1016/s0960-9822(01)00581-4)  
  - animals: Vertebrates, Insects, Arthropods  
  - Review of central pattern generators: rhythmic motor output can be produced with no sensory or descending timing cues, and neuromodulation reconfigures the same circuit to generate multiple patterns — the framework in which one spinal network can produce many gaits, with general principles drawn from invertebrate circuits and vertebrate spinal cord and brainstem.
- **Jessell 2000** — [Neuronal specification in the spinal cord: inductive signals and transcriptional codes](https://doi.org/10.1038/35049541)  
  - animals: Vertebrates  
  - Reviews the developmental logic of spinal neuronal specification: inductive signals and transcriptional codes assign class identity to progenitors at defined positions in the neural tube, assembling neural circuits with the selective connectivity that underlies the mature organism's behavioral repertoire.
- **Barbeau 1999** — [Tapping into spinal circuits to restore motor function.](https://doi.org/10.1016/s0165-0173(99)00008-9)  
  - animals: Cat, Human, Vertebrates  
  - Multidisciplinary review arguing that electrically activating spinal interneuronal circuits (reflex and pattern-generating) could coordinate many muscles at once for neuroprostheses, far more tractable than stimulating each muscle individually; draws on phase-adaptable cat hindlimb reflexes, chick embryo rhythmogenic network development, the cat spinal locomotor pattern generator, and locomotor training in incomplete SCI patients to identify candidate circuits and the control problems engineers must solve.  
  - *Robot/sim:* Model the control layer for FNS-style stimulation as a spinal circuit hierarchy (phase-modulated reflexes plus CPG) driving a musculoskeletal limb, instead of per-muscle open-loop stimulation.
- **Tresch 1999** — [The construction of movement by the spinal cord](https://doi.org/10.1038/5721)  
  - animals: Vertebrates, Frog  
  - We used a computational analysis to identify the basic elements with which the vertebrate spinal
cord constructs one complex behavior. This analysis extracted a small set of muscle synergies from
the range of muscle activations generated by cutaneous stimulation of the frog hindlimb. The
flexible combination of these synergies was able to account for the large number of different motor
patterns produced by different animals. These results therefore demonstrate one strategy used by
the vertebrate nervous system to produce movement in a computationally simple manner.
- **Ekeberg 1998** — [A neuro-mechanical model of legged locomotion: single leg control](https://doi.org/10.1007/s004220050468)  
  - animals: Vertebrates  
  - Neuro-mechanical single-leg model uniting CPG and peripheral feedback: a neural phase generator gates fast sensory feedback pathways according to the current step phase and filters inconsistent afferent input, while sensory signals both act directly on the motor system and entrain the phase generator; with elastic actuators the model steps stably over a large velocity range.  
  - *Robot/sim:* Implement a neural phase generator that gates fast feedback per step phase, with elastic actuators and sensory entrainment of the generator; ablate the gating so feedback is always on, and test for over-control and instability.
- **Rossignol and Dubuc 1994** — [Spinal pattern generation](https://doi.org/10.1016/0959-4388(94)90139-2)  
  - animals: Rat, Lamprey, Vertebrates  
  - Review framing spinal pattern generation around three foci: transmitter actions on rhythms in reduced preparations, activity-dependent changes in membrane properties of the generating circuits, and CPG-afferent interactions, emphasizing that new reflex responses and membrane properties emerge only during rhythmic activity.  
  - *Robot/sim:* Implement CPG circuits whose reflex pathway gains and signs are gated by rhythmic state, and test whether reflex responses absent in quiescent simulation appear when the CPG is active, mirroring state-dependent reflex modification.
- **Rome 1993** — [How Fish Power Swimming](https://doi.org/10.1126/science.8332898)  
  - animals: Vertebrates  
  - Isolated red muscle bundles driven through the in vivo length changes and stimulation patterns of steady swimming showed that most swimming power comes from posterior musculature while the anterior mainly transmits force - the opposite of the prevailing assumption - and that contractile properties are regionally tuned along the body to enhance power production.  
  - *Robot/sim:* In undulatory swimmer models, assign posterior-dominant power production with rostral force transmission and regionally tuned contractile properties instead of uniform actuation, and quantify the effect on swimming performance.
- **Grillner 1985** — [Neurobiological Bases of Rhythmic Motor Acts in Vertebrates](https://doi.org/10.1126/science.3975635)  
  - animals: Vertebrates, Mammals, Lamprey  
  - Classic formulation of the principles underlying nervous control of innate rhythmic motor acts in vertebrates, establishing the in vitro lamprey spinal cord preparation in which the locomotor motor pattern can be elicited in isolated cord sections for detailed mechanism analysis, and proposing that parts of the locomotor network are reused in other motor acts, including learned ones.  
  - *Robot/sim:* Implement modular CPG subnetworks shared across multiple motor programs in a simulated swimmer or walker and verify that reuse preserves timing and stability while reducing controller complexity.
