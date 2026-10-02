# Cluster 6: Proprioceptor & neuron models

49 papers · years 1952–2023 · 38 curated

## Top authors

- Izhikevich (2)
- Crowe (2)
- Schiefer (2)
- Hodgkin (1)
- Matthews (1)
- Rall (1)
- Han (1)
- Lv (1)
- Nagumo (1)
- Jordan (1)

## Pathways represented

- Biomechanically mediated preflexive feedback (1)

## Key papers (by citations)

- **Hodgkin 1952** — [A quantitative description of membrane current and its application to conduction and excitation in nerve](https://doi.org/10.1113/jphysiol.1952.sp004764)  
  - Title-level grounding only, because the supplied abstract is publisher boilerplate with no scientific content: this is the foundational quantitative description of membrane current and its application to conduction and excitation in nerve, i.e. the origin of conductance-based membrane modeling that later neuron-level CPG models build on. Re-fetch the true abstract before using this record for mechanism claims.
- **Izhikevich 2003** — [Simple model of spiking neurons](https://doi.org/10.1109/tnn.2003.820440)  
  - animals: Human  
  - Computational model study. As we develop such large-scale brain models consisting of spiking neurons, we must find compromises between two seemingly mutually exclusive requirements: The model for a single neuron must be: 1) computationally simple, yet 2) capable of producing rich firing patterns exhibited by real biological neurons.
- **Izhikevich 2006** — [Dynamical Systems in Neuroscience: The Geometry of Excitability and Bursting](https://doi.org/10.7551/mitpress/2526.001.0001)  
  - Standard reference connecting electrophysiology to nonlinear dynamical systems theory: neuronal information processing depends on dynamical as well as electrophysiological properties, with excitability, bursting, and synchronization treated geometrically (phase planes, bifurcations) starting from one- and two-dimensional Hodgkin-Huxley-type models. Provides the analysis toolkit for classifying the operating regime and burst-termination mechanisms of CPG neuron models.  
  - *Robot/sim:* Use its phase-plane and bifurcation methods to classify the operating regime and burst-termination mechanism of CPG half-center and persistent-sodium neuron models before porting them across simulators, where bifurcation regime — not parameter identity — determines transfer.
- **Matthews 1964** — [Muscle Spindles and Their Motor Control](https://doi.org/10.1152/physrev.1964.44.2.219)  
  - animals: Cat  
  - Computational model study. The existence of two distinct types of afferent nerve ending, however, has been widely suspected since Ruflini described the histologically distinct primary and secondary endings, while very recent histological work strongly suggests the existence of two distinct types of intrafusal muscle fiber with independent motor...
- **Rall 1962** — [Electrophysiology of a Dendritic Neuron Model](https://doi.org/10.1016/s0006-3495(62)86953-7)
- **Han 1995** — [Dephasing and Bursting in Coupled Neural Oscillators](https://doi.org/10.1103/physrevlett.75.3190)  
  - Diffusive voltage coupling between nonlinear neural oscillators can dephase rather than synchronize them, producing a bursting-like regime in which oscillators synchronize at small oscillation amplitude and desynchronize once amplitude grows large — a mechanism by which coupling structure alone generates antiphase and amplitude switching in oscillator pairs.  
  - *Robot/sim:* Couple half-center or inter-leg CPG oscillators diffusively in the voltage variable and sweep coupling strength and amplitude to check for dephasing and amplitude-switching bursts; relevant whenever naive diffusive coupling is assumed to be synchronizing.
- **Crowe 1964** — [The effects of stimulation of static and dynamic fusimotor fibres on the response to stretching of the primary](https://doi.org/10.1113/jphysiol.1964.sp007476)
- **Lv 2016** — [Multiple modes of electrical activities in a new neuron model under electromagnetic radiation](https://doi.org/10.1016/j.neucom.2016.05.004)
- **Nagumo 1972** — [On a response characteristic of a mathematical neuron model](https://doi.org/10.1007/bf00290514)
- **Jordan 1998** — [Initiation of locomotion in mammals.](https://doi.org/10.1111/j.1749-6632.1998.tb09040.x)  
  - animals: Mammals  
  - Reviews brainstem initiation of locomotion: c-fos-labeled active cells during locomotion localize the mesencephalic locomotor region across the periaqueductal gray, cuneiform nucleus, pedunculopontine nucleus, and locus coeruleus, with different subsets of locomotor nuclei apparently activated in different behavioral contexts (exploratory, appetitive, defensive), implying that lesion studies must be interpreted relative to the nuclei active in the specific locomotor task tested.  
  - *Robot/sim:* Implement a supraspinal initiation module with context-selective drives (exploratory, appetitive, defensive) converging on the spinal locomotor network, mimicking the brainstem command structure; ablate submodules to mirror the lesion-study interpretation caveats.
- **Gerstner 2009** — [How Good Are Neuron Models?](https://doi.org/10.1126/science.1181936)  
  - Reports on a competition in which modelers predicted neuronal activity and compares which neuron model performed best - a benchmark of single-neuron model classes against experimental data.
- **Schiefer 2016** — [Sensory feedback by peripheral nerve stimulation improves task performance in individuals with upper limb loss](https://doi.org/10.1088/1741-2560/13/1/016001)  
  - animals: Human  
  - 13 016001 View the article online for updates and enhancements.
- **Tsumoto 2006** — [Bifurcations in Morris–Lecar neuron model](https://doi.org/10.1016/j.neucom.2005.03.006)
- **Mileusnic 2006** — [Mathematical Models of Proprioceptors. I. Control and Transduction in the Muscle Spindle](https://doi.org/10.1152/jn.00868.2005)  
  - animals: Cat  
  - Computational model study. In the case of simultaneous static and dynamic fusimotor efferent stimulation, we demonstrated the importance of including the experimentally observed effect of partial occlusion.
- **Edin 1990** — [Dynamic response of human muscle spindle afferents to stretch](https://doi.org/10.1152/jn.1990.63.6.1297)  
  - animals: Human · afferents: Ia, II, Ib  
  - In human radial-nerve recordings from finger extensor muscles, three discrete and statistically pairwise-independent response markers, the initial burst at stretch onset, the deceleration response at the start of hold, and prompt silencing during imposed shortening, discriminate muscle spindle primaries from secondaries (the dynamic index alone is a poor discriminator with a unimodal distribution), while Golgi tendon organ afferents give negligible stretch responses; the battery supports a probability-based classification of human muscle afferents.  
  - *Robot/sim:* Implement Ia and II spindle models that reproduce the discrete dynamic markers (initial burst, deceleration response, silence during shortening) and validate the model's response statistics against these human single-unit data.
- **Crowe 1964** — [Further studies of static and dynamic fusimotor fibres](https://doi.org/10.1113/jphysiol.1964.sp007477)
- **Macefield 2018** — [Functional properties of human muscle spindles](https://doi.org/10.1152/jn.00071.2018)  
  - animals: Cat, Human  
  - Review of functional properties human muscle spindles. First published April 18, 2018; doi:10.1152/ jn.00071.2018.—Muscle spindles are ubiquitous encapsulated mechanoreceptors found in most mammalian muscles.
- **Kobayashi 2009** — [Made-to-order spiking neuron model equipped with a multi-timescale adaptive threshold](https://doi.org/10.3389/neuro.10.009.2009)  
  - A spiking neuron model combining a nonresetting leaky integrator with a multi-timescale adaptive threshold reproduces and predicts diverse spike responses — regular spiking, intrinsic bursting, fast spiking, and chattering — with only three adaptive threshold parameters, expressing firing characteristics as a continuum rather than discrete classes at low computational cost.  
  - *Robot/sim:* Use the multi-timescale adaptive-threshold IF unit as a low-cost spiking neuron class in CPG networks; sweep its three parameters across regular-spiking to intrinsic-bursting regimes and measure how heterogeneous unit classes influence network rhythm properties.
- **Van Geit 2008** — [Automated neuron model optimization techniques: a review](https://doi.org/10.1007/s00422-008-0257-6)  
  - Review of automated parameter optimization for computational neuron models, distinguishing feature-based, point-by-point voltage-trace, and multi-objective error functions, and detailing search algorithms including brute force, simulated annealing, genetic algorithms, evolution strategies, differential evolution, and particle-swarm optimization; the Neurofitter package combines a phase-plane trajectory density fitness function with several of these searches, replacing intractable hand tuning as model complexity grows.  
  - *Robot/sim:* Adopt automated multi-objective optimization (feature-based or trace-based fitness with evolutionary or annealing search) for tuning CPG and reflex-model parameters instead of hand tuning, per the reviewed error-function and search-algorithm taxonomy.
- **Procházka and Gorassini 1998** — [Models of ensemble firing of muscle spindle afferents recorded during normal locomotion in cats](https://doi.org/10.1111/j.1469-7793.1998.277bu.x)  
  - animals: Cat · afferents: Ia  
  - Mathematical models of muscle spindle primary (presumed Ia) firing, validated against chronically recorded hamstring afferents over 132 cat step cycles, are dominated by muscle velocity — firing rate scales approximately with the square root of velocity (power-law exponents 0.5-0.6) — with only small EMG-linked fusimotor components; inverse models recover muscle length accurately from firing profiles, implying the CNS could decode muscle length simply from Ia ensemble firing.  
  - *Robot/sim:* Use a square-root-of-velocity-dominant spindle encoder with a small EMG-linked fusimotor term as the Ia block in reflex models, and run it in reverse to estimate muscle length from afferent firing for state estimation in a neurally controlled walker.
- **Hulliger 1977** — [Effects of combining static and dynamic fusimotor stimulation on the response of the muscle spindle primary en](https://doi.org/10.1113/jphysiol.1977.sp011840)  
  - animals: Cat · afferents: Ia  
  - In cat soleus muscle spindle primary endings under sinusoidal stretch, combined static and dynamic fusimotor stimulation is dominated by static action at small amplitudes (up to about 50 micrometers) with the dynamic contribution growing progressively with stretch amplitude; at response peaks the two actions sum (dynamic stronger), whereas in troughs static action occludes the weaker dynamic action, and phase differences between conditions remain below 20 degrees. The findings characterize how fusimotor set shapes Ia encoding of stretch.  
  - *Robot/sim:* Implement spindle-primary afferent encoding with separate static and dynamic fusimotor drives reproducing the amplitude-dependent summation and occlusion behavior, and modulate fusimotor set during locomotion to shape Ia feedback.
- **Achuthan 2009** — [Phase-Resetting Curves Determine Synchronization, Phase Locking, and Clustering in Networks of Neural Oscillat](https://doi.org/10.1523/jneurosci.0426-09.2009)  
  - Iterated-map analysis based solely on phase-resetting curves predicts synchronization, phase-locking, and clustering in neural oscillator networks without weak-coupling assumptions: synchrony stability depends on the PRC slope at spike initiation, and between-cluster splay timing can enforce synchrony on subclusters incapable of synchronizing alone.  
  - *Robot/sim:* Use measured phase-resetting curves of coupled oscillators (half-centers or limb CPGs) with the iterated-map stability criteria to predict interlimb phase locking and clustering, selecting coupling strengths that avoid synchrony loss.
- **Sakaguchi 2006** — [Instability of synchronized motion in nonlocally coupled neural oscillators](https://doi.org/10.1103/physreve.73.031907)  
  - Shows that the synchronized state of nonlocally coupled Hodgkin-Huxley neurons with excitatory and inhibitory synaptic coupling is linearly unstable, and that its destabilization yields chimera states, wavy states, clustering, and spatiotemporal chaos, providing a mechanism by which spatially structured, partially synchronized activity can emerge in distributed oscillator networks without parameter heterogeneity.  
  - *Robot/sim:* Build a nonlocally coupled oscillator network with mixed excitatory/inhibitory coupling and test the stability of the synchronized state; parameter sweeps should reproduce chimera, clustering, wavy, and chaotic regimes, a caution for multi-limb CPG designs that assume synchrony.
- **Yakovenko 2004** — [Contribution of stretch reflexes to locomotor control: a modeling study](https://doi.org/10.1007/s00422-003-0449-z)  
  - afferents: Ia, Ib · pathways: Biomechanically mediated preflexive feedback  
  - A two-legged neuromuscular model shows that intrinsic muscle stiffness driven by a stereotyped rhythmic pattern sustains surprisingly stable gait without sensory control, that spindle Ia and tendon organ Ib stretch-reflex feedback (35 ms latency, 30% of the EMG profile, active only when the receptor-bearing muscle contracts) is crucial for load compensation only when central drive is low, and that finite-state control greatly extends adaptive capability.  
  - *Robot/sim:* Add state-gated autogenic Ia/Ib stretch feedback (fixed latency, bounded contribution to motoneuron activation) to a feedforward-driven Hill-muscle biped and ablate it at high versus low central drive to reproduce the state-dependent stabilizing contribution.
- **Bal 1988** — [The pyloric central pattern generator in Crustacea: a set of conditional neuronal oscillators](https://doi.org/10.1007/bf00604049)  
  - This was shown by isolating in situ pyloric neurons from all other elements in the pyloric network.
