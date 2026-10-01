# Feedback pathway: Ib inhibition

4 papers in the corpus.

- **Nichols and Ross 2009** — [The Implications of Force Feedback for the λ Model](https://doi.org/10.1007/978-0-387-77064-2_36)  
  - animals: Cat, Mammals · pathways: Ib excitatory; Ib inhibition  
  - In the λ-model framework, autogenic length feedback compensates muscle nonlinearities while positive force feedback — during level cat stepping largely restricted to gastrocnemius — reinforces the stiff ankle–knee linkage; heterogenic inhibitory force feedback spans different joints and axes, so it coordinates interjoint action and shifts activation thresholds, making threshold a feedback-dependent quantity rather than a pure descending control variable.
- **Duysens 2000** — [Load-regulating mechanisms in gait and posture: Comparative aspects.](https://doi.org/10.1152/physrev.2000.80.1.83)  
  - pathways: Ib inhibition  
  - Type Ib feedback has an inhibitory effect on muscles in non-moving animals
- **Procházka 1997** — [Implications of Positive Feedback in the Control of Movement](https://doi.org/10.1152/jn.1997.77.6.3237)  
  - animals: Mammals · afferents: Ib · pathways: Ib inhibition; Ib excitatory  
  - Argues with analytic neuromuscular models that tendon-organ afferents switch from negative force feedback in static posture to positive force feedback in locomotion, and that positive force feedback is stable and effective because muscle length-tension properties provide automatic gain control - stability is further supported by concomitant negative displacement feedback and, unexpectedly, by delays in the positive feedback pathway - making positive force feedback a load-compensation mechanism complementing negative displacement and velocity feedback.  
  - *Robot/sim:* Implement stance-phase positive autogenic force feedback (Ib excitatory) on extensors with length-tension automatic gain control plus negative displacement feedback; demonstrate stability at strong gains and with realistic delays, then remove the gain control to show instability.
- **Pearson and Collins 1993** — [Reversal of the influence of group Ib afferents from plantaris on activity in medial gastrocnemius muscle duri](https://doi.org/10.1152/jn.1993.70.3.1009)  
  - animals: Cat · afferents: Ib, Ia · pathways: Ib excitatory; Ib inhibition  
  - In clonidine-treated acute and chronic spinal cats, group I stimulation of the plantaris nerve entrained the locomotor rhythm and, during locomotor activity, group Ib afferents from plantaris exerted an excitatory action on medial gastrocnemius bursts (30-50 ms latency) through the extensor half-center of the rhythm generator — reversing the inhibitory effect the same stimuli produced on tonic activity without locomotion; selective Ia activation by muscle vibration neither entrained the rhythm nor augmented the bursts.  
  - *Robot/sim:* Implement state-dependent Ib feedback in a half-center CPG: autogenic inhibition at rest switching to excitation of the extensor half-center during locomotion, with rhythm entrainment by group I stimulation but not by Ia-specific input; ablate the reversal to test loss of load-dependent stance burst augmentation.
