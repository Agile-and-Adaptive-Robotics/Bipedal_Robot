# Feedback pathway: type II excitatory

5 papers in the corpus.

- **Hatz 2012** — [Control of ankle extensor muscle activity in walking cats](https://doi.org/10.1152/jn.00944.2011)  
  - animals: Cat, Mammals · pathways: Ib excitatory; type II excitatory  
  - In conscious cats with an isolated medial gastrocnemius, stance-phase and slope-dependent modulation of ankle-extensor activity is well described by constant central drive with constant proprioceptive gains: Ib force feedback is the primary modulator, group II adds a small tonic contribution, and Ia feedback contributes none — terrain compensation without changing central commands or nervous-system gains.
- **Mazzaro 2006** — [Afferent-mediated modulation of the soleus muscle activity during the stance phase of human walking](https://doi.org/10.1007/s00221-006-0451-5)  
  - animals: Human · afferents: Ia, II, Cutaneous · pathways: type II excitatory  
  - In human walking, slow small-amplitude ankle dorsiflexion perturbations during stance modulate ongoing soleus activity through group II afferent feedback (response reduced by tizanidine depression of group II pathways) with at least a partial group Ia contribution; blocking sensory feedback from the foot had no effect, and body-load changes shifted baseline activity without changing response amplitude. Group II and possibly load-sensitive afferents govern stance-phase extensor amplitude, whereas foot cutaneous and intrinsic proprioceptive afferents do not contribute.  
  - *Robot/sim:* Implement ankle stretch-mediated group II feedback (plus partial Ia) onto stance-phase soleus/extensor motoneurons; reproduce the perturbation-response amplitude, its persistence under load change, and its reduction when group II gain is removed.
- **Donelan and Pearson 2004** — [Contribution of sensory feedback to ongoing ankle extensor activity during the stance phase of walking](https://doi.org/10.1139/y04-043)  
  - animals: Human, Cat · afferents: II, Ib · pathways: Ib excitatory; type II excitatory  
  - Quantitative review establishing that load-related sensory feedback contributes up to 60% of ongoing ankle extensor activity during the stance phase of walking, with secondary spindle endings (human) and Golgi tendon organs (human and cat) the likely receptors — autogenic positive feedback reinforcing stance extensor force. Argues that resolving which receptor groups set extensor magnitude across locomotor tasks requires network simulations coupled to forward-dynamic musculoskeletal models, since experimental strategies alone cannot dissociate the distributed contributions.  
  - *Robot/sim:* Implement autogenic positive force feedback (Ib) and length/velocity feedback (II) onto stance-gated ankle extensor motoneurons; ablate each channel in turn and test whether unloading reduces extensor activity by up to 60% as reported.
- **Perreault 1995** — Effects of stimulation of hindlimb flexor group II afferents during fictive locomotion in the cat  
  - animals: Cat · pathways: Ia or II stance to swing; type II excitatory  
  - Excitatory input to flexor type II afferent.
- **Edgley and Jankowska 1987** — [An interneuronal relay for group I and II muscle afferents in the midlumbar segments of the cat spinal cord.](https://doi.org/10.1113/jphysiol.1987.sp016676)  
  - animals: Cat · afferents: Ia, Ib, II, Cutaneous · pathways: type II excitatory; II inhibitory  
  - Midlumbar (L4) interneurons with dominant group II input - also driven by group I, cutaneous, joint, and interosseous afferents - project to hindlimb motor nuclei and are activated monosynaptically by group I and II afferents, making them a disynaptic relay that adds 0.4-2.5 ms latency to length- and force-afferent actions on motoneurons relative to caudal lumbar pathways.  
  - *Robot/sim:* Insert a midlumbar relay interneuron layer (adding 0.4-2.5 ms latency) between group I/II afferents and motoneuron pools; ablate it to test the functional cost of losing the longer-latency afferent actions during perturbed stepping.
