# Feedback pathway: Ia presynaptic inhibition

3 papers in the corpus.

- **Gosgnach 2000** — [Depression of group Ia monosynaptic EPSPs in cat hindlimb motoneurones during fictive locomotion.](https://doi.org/10.1111/j.1469-7793.2000.00639.x)  
  - animals: Cat · pathways: Ia presynaptic inhibition  
  - During brainstem-evoked fictive locomotion in decerebrate cats, monosynaptic group Ia EPSPs are tonically depressed to about two-thirds of control in most motoneurones, with group I field potentials similarly reduced and only weak correlation with decreased motoneurone input resistance — evidence that presynaptic inhibition of Ia transmission underlies the tonic depression of stretch reflexes during locomotion.
- **Perreault 1999** — [Depression of muscle and cutaneous afferent-evoked monosynaptic field potentials during fictive locomotion in ](https://doi.org/10.1111/j.1469-7793.1999.00691.x)  
  - animals: Cat · afferents: Ia, Ib, II, Cutaneous · pathways: Ia presynaptic inhibition  
  - During MLR-evoked fictive locomotion in decerebrate cats, monosynaptic field potentials evoked by group I, group II, and cutaneous afferents are tonically depressed (intermediate-lamina group II potentials down to a mean of 49% of control), with smaller phase-dependent cyclic modulation superimposed. The depression begins with tonic MLR drive before rhythm onset, indicating that the locomotor state reduces synaptic transmission from primary afferents onto first-order spinal interneurons — afferent influx to the locomotor circuitry is gated down at its earliest synapse during locomotion.  
  - *Robot/sim:* Implement locomotor-state gating on afferent input channels: scale group I/II/cutaneous afferent synaptic gains by a tonic depression factor (~0.8 for dorsal/group I, ~0.5 for intermediate group II) plus a small phase-dependent component, and verify reflex amplitudes are suppressed during locomotion while phase-cycling as observed.
- **Jankowska 1967** — [The Effect of DOPA on the Spinal Cord 6. Half‐centre organization of interneurones transmitting effects from t](https://doi.org/10.1111/j.1748-1716.1967.tb03637.x)  
  - animals: Cat · afferents: Flexor reflex afferents, Ia · pathways: Ia presynaptic inhibition  
  - Classic DOPA study revealing the half-center organization of interneurons transmitting late flexor-reflex-afferent effects in spinal cats: reciprocally organized populations activated from ipsilateral versus contralateral FRA project to flexor versus extensor motoneurons, while a third type depolarizes Ia afferent terminals, adding presynaptic inhibition to the long-latency alternating discharge.  
  - *Robot/sim:* Implement reciprocally organized FRA-driven half-centers plus a presynaptic-inhibition channel onto Ia terminals; ablate the Ia presynaptic inhibition to test how afferent gating shapes the DOPA-like alternating rhythm.
