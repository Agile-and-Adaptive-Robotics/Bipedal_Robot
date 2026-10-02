# Animal: Turtle

4 papers in the corpus.

- **Stein 2018** — [Central pattern generators in the turtle spinal cord: selection among the forms of motor behaviors.](https://doi.org/10.1152/jn.00602.2017)  
  - animals: Turtle · afferents: Cutaneous  
  - The turtle spinal cord contains separate CPGs for three forms of scratch, two forms of swim, and one flexion reflex, with multisecond motor memories, interneurons partially shared across behaviors (enabling motor pattern blends) alongside behavior-specific cells, and deletion analysis supporting unit burst generators organized per degree-of-freedom direction (hip-flexor, hip-extensor, knee-flexor, etc.) - an organizational complexity the classic flexor/extensor half-center scheme cannot capture.  
  - *Robot/sim:* Implement the CPG as unit burst generators per joint direction under a selection layer combining shared and behavior-specific interneurons, and reproduce motor pattern blends and deletions as signatures of that organization.
- **Stein 2016** — [Modular organization of the multipartite central pattern generator for turtle rostral scratch: knee-related in](https://doi.org/10.1152/jn.00871.2015)  
  - animals: Turtle  
  - Single-unit interneuron recordings during knee-related deletions in turtle rostral scratch show that knee-extensor and knee-flexor interneurons occupy modules separate from the hip-extensor and hip-flexor modules, establishing at the single-unit level that the limb CPG is multipartite, organized per joint and function rather than as unitary flexor and extensor halves.  
  - *Robot/sim:* Implement a multipartite CPG with separate hip and knee flexor/extensor modules per joint; reproduce deletions by transiently silencing individual modules and verify that other modules continue bursting, then compare against a classical bipartite half-center controller.
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
- **Robertson and Stein 1988** — [Synaptic control of hindlimb motoneurones during three forms of the fictive scratch reflex in the turtle.](https://doi.org/10.1113/jphysiol.1988.sp017281)  
  - animals: Turtle  
  - Intracellular recordings for three distinct forms of turtle scratch reflex in response to tactile stimulation of body surface.
