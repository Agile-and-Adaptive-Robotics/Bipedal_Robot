# Afferent type: III/IV

2 papers in the corpus.

- **Grillner 2011** — [Control of Locomotion in Bipeds, Tetrapods, and Fish](https://doi.org/10.1002/cphy.cp010226)  
  - animals: Vertebrates · afferents: Ia, II, III/IV, Cutaneous  
  - Comprehensive synthesis of locomotor control across bipeds, tetrapods, and fish: spinal unit-burst-generator organization, interlimb coordination, load- and hip-position-dependent reflex regulation of the spinal CPG, tonic gating of reflexes from skin and group II/III afferents, and brainstem and cerebellar control - the canonical reference architecture for locomotor CPG models.  
  - *Robot/sim:* Use the chapter's architecture - unit burst generators per limb/joint, interlimb coordination pathways, and load and hip-position feedback gating the rhythm - as the skeleton for a locomotion controller, with the reviewed reflex gating (skin, group II/III) implemented as state-dependent pathway gains.
- **Schomburg 1998** — [Flexor reflex afferents reset the step cycle during fictive locomotion in the cat](https://doi.org/10.1007/s002210050522)  
  - animals: Cat · afferents: Flexor reflex afferents, II, III/IV, Cutaneous  
  - Brief trains to flexor reflex afferents (FRA: joint, cutaneous, group II-III muscle afferents) reset the L-dopa-induced fictive locomotor rhythm in high-spinal cats - interrupting extensor activity and initiating flexion when delivered in extension, or prolonging flexion in late flexion - demonstrating that FRA interneurones are constituent elements of the rhythm-generating network organized along half-centre lines.  
  - *Robot/sim:* Implement FRA input as a phase-dependent reset of the half-centre oscillator (extensor interruption plus flexion trigger, flexion prolongation late in flexion) and verify that brief afferent trains shorten or lengthen single cycles after which the original rhythm period resumes.
