# Feedback pathway: Ia stance to swing

2 papers in the corpus.

- **Ivashko 2003** — [Modeling the spinal cord neural circuitry controlling cat hindlimb movement during locomotion](https://doi.org/10.1016/S0925-2312(02)00832-9)  
  - animals: Cat · pathways: Ib disynaptic excitation; Ib stance to swing; Ib swing to stance; Ia stance to swing; Ia swing to stance; Ia monosynaptic; Ia monosynaptic excitation; Mechanosensory monosynaptic excitation; Ia disynaptic inhibition  
  - Abstract
A computational model of the spinal cord neural circuitry that controls locomotor movements of simulated cat hindlimbs. The neural circuitry includes two central pattern generators integrated with reflex circuits. All neurons were modeled in the Hodgkin–Huxley style. The musculoskeletal system includes two three-joint hindlimbs and the trunk. Each
hindlimb is actuated by nine one- and two-joint muscles (a Hill-type model). Our simulations allow us to suggest a specific network architecture in the spinal cord and a pattern of feedback connectivities (from Ia and Ib fibers and touch s
- **Verschueren 2003** — [Vibration-Induced Changes in EMG During Human Locomotion](https://doi.org/10.1152/jn.00863.2002)  
  - animals: Human · afferents: Ia · pathways: Ia stance to swing  
  - Continuous tendon vibration during blindfolded walking enhanced the EMG of quadriceps femoris and biceps femoris (mainly in stance) and quadriceps vibration advanced the onset of tibialis anterior within the gait cycle, while vibration of ankle and hip muscles had no significant amplitude effect - implicating Ia afferent input in stance-phase muscle activation and in a limited role in triggering phase transitions (stance to swing) during human locomotion.  
  - *Robot/sim:* Implement Ia input as an amplitude contribution to stance-phase muscle activation plus a weak phase-transition trigger term driven by knee-extensor Ia activity; ablating the timing term should remove the advanced dorsiflexor onset while leaving the stance EMG enhancement intact.
