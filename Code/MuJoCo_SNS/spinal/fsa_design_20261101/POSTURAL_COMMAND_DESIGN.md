# Descending posture commands and the existing balance pathway

Session 1, easteregg2. Design explanation in response to Ben's antagonist-subnetwork diagram. No new body trial or training.

The proposed organization is appropriate: standing commands can feed an extensor/flexor subnetwork while local Ia, II, Ib, and Renshaw circuits retain their functions. During stepping, postural commands can coexist with rhythm/PF drive. A tonic posture reference, sensory corrections, and coactivation are distinct signals.

The current implementation only partly expresses this organization:

- `Code\MuJoCo_SNS\spinal\build_network.py` contains POSTURE-to-MN connections, a POSTURE-to-RG-extensor connection, descending BAL populations, and per-MN POST inputs. Its optional coact parameter recruits antagonist posture groups. This is an existing postural-command mechanism, although it does not duplicate the complete PF-to-interneuron organization in Ben's diagram.
- `Code\MuJoCo_SNS\spinal\hunt_training\balance.py` adds four directional PT_BAL populations and four inhibitory relays. Delayed CoM displacement/velocity and optional acceleration drive them; fixed, pose-calibrated recruitment connects them to MNs. This is task-level feedback. It does not create a complete joint-by-joint descending posture layer, and its desired support reference is captured at initialization.
- The CAD controller has an engineering posture/force allocation path. That must be distinguished from a posture command expressed entirely through the local neural subnetwork.

I implemented the PT pathway as an added task-feedback module and retained the existing standing bias. That left the explicit descending postural subnetwork incompletely developed. It is not a physiological reason to omit it.

For an explicit extension, the joint subnetwork would receive separate nonnegative posture commands p_E and p_F, alongside locomotor PF activity r_E and r_F. For example, descending motor input has the declared form I_desc,E = w_post,E p_E + w_PF,E r_E, with separately designed projections to the appropriate inhibitory interneurons. The motor cells still integrate actual afferent and recurrent synaptic currents; the command is not a motor-voltage clamp.

The postural reference specifies desired alignment or joint position. Local muscle-length/velocity/load feedback and task-level balance error provide corrections. Common antagonist drive can adjust coactivation. Mapping task error to these commands must state its signs, units, gains, saturation, and delay; position offset alone is insufficient to describe dynamic balance. The existing PT recruitment is an engineering design, not a literal Peterka/Ting circuit replication.

Before adopting the extension: draw its actual projections, distinguish posture from PF and sensory inputs, derive operating-point gains with shunting, check mirrored anatomical functions, and test each pathway and combined command. Positive load feedback must retain its declared phase dependence. This memo specifies the design option; it does not claim that the full replacement has been built or validated. Any training/curriculum still awaits Ben's OK.
