# Feedback pathway: Ib disynaptic excitation

6 papers in the corpus.

- **Nichols 2018** — [Distributed force feedback in the Spinal Cord and the regulation of limb mechanics](https://doi.org/10.1152/jn.00216.2017)  
  - animals: Cat, Human · pathways: Ib disynaptic inhibition; Ia monosynaptic; Ib disynaptic excitation  
  - Review Update Paper; This paper show able to give insight to the inhibitory and excitatory force feedback during locomotion
- **Ivashko 2003** — [Modeling the spinal cord neural circuitry controlling cat hindlimb movement during locomotion](https://doi.org/10.1016/S0925-2312(02)00832-9)  
  - animals: Cat · pathways: Ib disynaptic excitation; Ib stance to swing; Ib swing to stance; Ia stance to swing; Ia swing to stance; Ia monosynaptic; Ia monosynaptic excitation; Mechanosensory monosynaptic excitation; Ia disynaptic inhibition  
  - Abstract
A computational model of the spinal cord neural circuitry that controls locomotor movements of simulated cat hindlimbs. The neural circuitry includes two central pattern generators integrated with reflex circuits. All neurons were modeled in the Hodgkin–Huxley style. The musculoskeletal system includes two three-joint hindlimbs and the trunk. Each
hindlimb is actuated by nine one- and two-joint muscles (a Hill-type model). Our simulations allow us to suggest a specific network architecture in the spinal cord and a pattern of feedback connectivities (from Ia and Ib fibers and touch s
- **Quevedo 2000** — [Group I disynaptic excitation of cat hindlimb flexor and bifunctional motoneurones during fictive locomotion.](https://doi.org/10.1111/j.1469-7793.2000.t01-1-00549.x)  
  - animals: Cat · afferents: Ia, Ib · pathways: Ib disynaptic excitation  
  - During fictive locomotion in decerebrate cats, group I stimulation evokes locomotor-dependent disynaptic excitation (single interneurone, mean latency 1.64 ms) of most flexor (89 percent) and many bifunctional (64 percent) motoneurones, largest during flexion and from homonymous nerves, evoked by both tendon-organ (Ib) and muscle-spindle (Ia) afferents and separate from the extensor group I pathway - indicating distinct flexor and extensor excitatory interneurone groups that reinforce ongoing locomotor activity throughout the limb.  
  - *Robot/sim:* Add disynaptic group I excitation of flexor and bifunctional motoneurone pools, phase-weighted to peak during flexion and driven by homonymous afferents, as load-reinforcement of ongoing activity; ablate it to quantify the loss of flexor reinforcement.
- **Pearson 1998** — [Enhancement and Resetting of Locomotor Activity by Muscle Afferentsa](https://doi.org/10.1111/j.1749-6632.1998.tb09050.x)  
  - animals: Mammals · afferents: Ib · pathways: Ib stance to swing; Ib disynaptic excitation; Ib excitatory  
  - Feedback from ankle extensor group I afferents sets two stance-critical features of the mammalian locomotor pattern: Ib Golgi tendon organ input excites the extensor half-center to prolong extensor bursts, so swing is not initiated while the leg is loaded, and extensor burst magnitude is enhanced through disynaptic and polysynaptic pathways opened only during locomotion. Reflex gains must be calibrated to limb biomechanics, and weakening ankle extensors evokes compensatory recalibration of proprioceptive influences on the pattern.  
  - *Robot/sim:* Implement stance-gated Ib excitation of the extensor half-center plus locomotion-gated disynaptic enhancement of extensor motoneuron gain; ablating the load gate should trigger swing under load and abolish load-proportional burst scaling.
- **McCrea 1995** — [Disynaptic group I excitation of synergist ankle extensor motoneurones during fictive locomotion in the cat.](https://doi.org/10.1113/jphysiol.1995.sp020897)  
  - animals: Cat · afferents: Ia, Ib · pathways: Ib disynaptic inhibition; Ib disynaptic excitation  
  - In decerebrate cats, plantaris nerve group I stimulation evokes short-latency disynaptic IPSPs in medial gastrocnemius motoneurons at rest but EPSPs during the extensor phase of MLR-evoked fictive locomotion, also elicited by selective group Ia activation via Achilles tendon stretch, indicating an excitatory group Ia and Ib feedback system that reinforces ongoing extensor activity during stance. In clonidine-treated spinal preparations the resting inhibition disappears without appearing as excitation, suggesting locomotion inhibits the inhibitory interneurons operating at rest.  
  - *Robot/sim:* Implement phase-switched group I pathways onto synergist ankle extensor motoneurons, disynaptic inhibition at rest and disynaptic excitation during the extensor phase; ablate the switch to test reinforcement of extensor activity during stance.
- **Jankowska 1981** — [Common interneurones in reflex pathways from group 1a and 1b afferents of ankle extensors in the cat.](https://doi.org/10.1113/jphysiol.1981.sp013556)  
  - animals: Cat · pathways: Ib disynaptic inhibition; Ib disynaptic excitation  
  - Most lamina V-VI interneurones take convergent input from both Ia spindle and Ib tendon-organ afferents of ankle and toe extensors — about 64% shared pathways showing co-excitation, co-inhibition, or opposite actions — so reflex circuitry is not strictly partitioned between length and force feedback, and interneurones projecting to motor nuclei show the same convergence.
