# Feedback pathway: Ia monosynaptic excitation

7 papers in the corpus.

- **Mayer et al. 2018** — [Role of muscle spindle feedback in regulating muscle activity strength during walking at different speed in mi](https://doi.org/10.1152/jn.00250.2018)  
  - animals: Mice · pathways: Ia monosynaptic excitation  
  - Muscle spindle feedback regulates the strength of muscle activity across walking speeds in mice — Ia feedback scales EMG amplitude, not only timing, so spindle loss degrades force output as speed increases.
- **Procházka and Prochazka 2011** — [Proprioceptive Feedback and Movement Regulation](https://doi.org/10.1002/cphy.cp120103)  
  - animals: Human · afferents: Ia, II, Ib, Cutaneous · pathways: Ia monosynaptic excitation; Ib excitatory  
  - Comprehensive treatment of proprioceptive feedback in movement regulation spanning muscle spindle, tendon organ, joint, ligament, skin, and invertebrate proprioceptors, their firing during active movement (including task-related fusimotor set), and a control-theoretic framing: stretch reflexes as proportional control, positive force feedback via tendon organs, finite-state logic of postural constraints, and inherent feedback from muscle properties.  
  - *Robot/sim:* Build the control hierarchy in simulation: proportional Ia feedback, positive Ib force feedback, finite-state (phase-gated) reflex switching, and muscle-property-mediated preflexes; ablate each layer to quantify its contribution to stable movement.
- **Ivashko 2003** — [Modeling the spinal cord neural circuitry controlling cat hindlimb movement during locomotion](https://doi.org/10.1016/S0925-2312(02)00832-9)  
  - animals: Cat · pathways: Ib disynaptic excitation; Ib stance to swing; Ib swing to stance; Ia stance to swing; Ia swing to stance; Ia monosynaptic; Ia monosynaptic excitation; Mechanosensory monosynaptic excitation; Ia disynaptic inhibition  
  - Abstract
A computational model of the spinal cord neural circuitry that controls locomotor movements of simulated cat hindlimbs. The neural circuitry includes two central pattern generators integrated with reflex circuits. All neurons were modeled in the Hodgkin–Huxley style. The musculoskeletal system includes two three-joint hindlimbs and the trunk. Each
hindlimb is actuated by nine one- and two-joint muscles (a Hill-type model). Our simulations allow us to suggest a specific network architecture in the spinal cord and a pattern of feedback connectivities (from Ia and Ib fibers and touch s
- **Zehr and Stein 1999** — [What functions do reflexes serve during human locomotion](https://doi.org/10.1016/s0301-0082(98)00081-1)  
  - animals: Lamprey, Cat, Human · afferents: Cutaneous · pathways: Cutaneous flexor excitation; Ia monosynaptic excitation  
  - Review establishing that reflexes during locomotion are task-, phase-, and context-dependent, with a division of labor: cutaneous reflexes alter swing-limb trajectory to avoid stumbling, stretch reflexes stabilize limb trajectory and assist force production during stance, and load receptor reflexes support body weight and influence step-cycle timing - functions dynamically reassigned across the cycle and clinically exploitable after neurotrauma.  
  - *Robot/sim:* Implement phase- and task-gated cutaneous (swing trajectory), Ia stretch (stance stabilization and force), and load-receptor (weight support, cycle timing) pathways in a walker; ablate each to reproduce the predicted functional losses.
- **Capaday et al. 1986** — [Amplitude modulation of the soleus H-reflex in the human during walking and standing](https://doi.org/10.1523/jneurosci.06-05-01308.1986)  
  - animals: Human · pathways: Ia monosynaptic excitation  
  - The soleus H-reflex is strongly modulated over the step cycle during walking versus standing; central gating of the monosynaptic Ia pathway is phase- and task-dependent.
- **Fleshman 1984** — [Peripheral and Central Control of Flexor Digitorum Longus and Flexor Hallucis Longus Motoneurons: The Synaptic](https://doi.org/10.1007/bf00235825)  
  - animals: Cat · afferents: Ia, Cutaneous · pathways: Ia monosynaptic excitation; Fictive locomotion without sensory feedback  
  - Cat FDL and FHL have identical mechanical actions yet are used differently in locomotion, and the differentiation is central: during fictive locomotion FHL motoneurons are coactive with ankle extensors in extension while FDL fires a brief burst at flexion onset; both pools share monosynaptic Ia excitation from both muscles, and only some cutaneous circuits differ (disynaptic EPSPs from saphenous and superficial peroneal nerves in FDL only; forelimb cutaneous IPSPs chiefly in FHL) - reflex organization alone cannot explain the functional divergence.  
  - *Robot/sim:* Drive two synergist pools at one joint with distinct CPG burst targets (extensor-coactive versus flexion-onset) while sharing monosynaptic Ia excitation across both; ablate the differentiated central drive to show common reflex organization cannot produce the functional divergence.
- **Eccles and Lundberg 1958** — [Integrative pattern of Ia synaptic actions on motoneurones of hip and knee muscles.](https://doi.org/10.1113/jphysiol.1958.sp006101)  
  - animals: Cat · pathways: Ia monosynaptic excitation  
  - Classic integrative map of Ia synaptic actions on hip and knee motoneurones: Ia excitation diverges broadly to synergists and across joints, fractionated into flexor and extensor recipient pools.
