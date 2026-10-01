# Feedback pathway: Ib disynaptic inhibition

5 papers in the corpus.

- **Nichols 2018** — [Distributed force feedback in the Spinal Cord and the regulation of limb mechanics](https://doi.org/10.1152/jn.00216.2017)  
  - animals: Cat, Human · pathways: Ib disynaptic inhibition; Ib disynaptic excitation; Ia monosynaptic  
  - Review Update Paper; This paper show able to give insight to the inhibitory and excitatory force feedback during locomotion
- **Ross and Nichols 2009** — [Heterogenic Feedback Between Hindlimb Extensors in the Spontaneously Locomoting Premammillary Cat](https://doi.org/10.1152/jn.90338.2008)  
  - animals: Cat · pathways: Ib disynaptic inhibition  
  - During spontaneous treadmill stepping in the premammillary decerebrate cat, force-dependent heterogenic inhibition between hindlimb extensors persists (quadriceps onto gastrocnemius, gastrocnemius onto plantaris/FHL) but distal-onto-proximal inhibition is weaker than during the crossed-extension reflex, yielding a proximal-to-distal gradient of Ib inhibition that supports interjoint coordination and limb stability.
- **Stephens and Yang 1996** — [Short latency, non-reciprocal group I inhibition is reduced during the stance phase of walking in humans.](https://doi.org/10.1016/s0006-8993(96)00977-8)  
  - animals: Human, Cat · afferents: Ib · pathways: Ib disynaptic inhibition  
  - In intact humans, short-latency non-reciprocal group I inhibition from the medial gastrocnemius nerve onto the conditioned soleus H-reflex, disynaptic and present at rest in most subjects, is significantly reduced during treadmill walking, with some subjects showing significant excitation, mirroring the reduction of resting Ib inhibition toward excitation seen in walking cats. The reduction may be partially accounted for by activation of the triceps surae itself.  
  - *Robot/sim:* Implement disynaptic Ib inhibition onto ankle extensors with a locomotor-phase-dependent gain reduction (partial reversal toward excitation during gait); ablate the phase switch to test its effect on stance extensor activity.
- **McCrea 1995** — [Disynaptic group I excitation of synergist ankle extensor motoneurones during fictive locomotion in the cat.](https://doi.org/10.1113/jphysiol.1995.sp020897)  
  - animals: Cat · afferents: Ia, Ib · pathways: Ib disynaptic inhibition; Ib disynaptic excitation  
  - In decerebrate cats, plantaris nerve group I stimulation evokes short-latency disynaptic IPSPs in medial gastrocnemius motoneurons at rest but EPSPs during the extensor phase of MLR-evoked fictive locomotion, also elicited by selective group Ia activation via Achilles tendon stretch, indicating an excitatory group Ia and Ib feedback system that reinforces ongoing extensor activity during stance. In clonidine-treated spinal preparations the resting inhibition disappears without appearing as excitation, suggesting locomotion inhibits the inhibitory interneurons operating at rest.  
  - *Robot/sim:* Implement phase-switched group I pathways onto synergist ankle extensor motoneurons, disynaptic inhibition at rest and disynaptic excitation during the extensor phase; ablate the switch to test reinforcement of extensor activity during stance.
- **Jankowska 1981** — [Common interneurones in reflex pathways from group 1a and 1b afferents of ankle extensors in the cat.](https://doi.org/10.1113/jphysiol.1981.sp013556)  
  - animals: Cat · pathways: Ib disynaptic inhibition; Ib disynaptic excitation  
  - Most lamina V-VI interneurones take convergent input from both Ia spindle and Ib tendon-organ afferents of ankle and toe extensors — about 64% shared pathways showing co-excitation, co-inhibition, or opposite actions — so reflex circuitry is not strictly partitioned between length and force feedback, and interneurones projecting to motor nuclei show the same convergence.
