# Descending posture commands and antagonist reflex circuits

Machine: **easteregg2**. Source/design audit only; no body run, training, manuscript change, or physiological equivalence claim.

## Actual CAD262 routes

Actual compiled census, not just builder capability: **0 RC, 0 IaIN, 0 IIIN, 0 IIX, and 0 IBIN populations; 8 IBEXC populations with36 active outgoing edges; 50 active POSTURE→MN edges and54 direct POST input ports**. Therefore the original CAD262 cannot be depicted as already containing the complete antagonist IaIN/IIIN/Renshaw circuit. Exact current-source topology is recorded in `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical\CAD_DESCENDING_TOPOLOGY_v1.json`.

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\run_neural_robot.py:68` implements `motor_target_current`. It calculates the instantaneous summed synaptic current onto each MN, cancels that current and the tonic/leak background, then imposes a voltage target proportional to the requested activation. The consumed default is **5 ms**, although its comment says 30 ms. The downstream actuator filter is separately 30 ms. Lines294–302 inject these computed currents into `POST_<muscle>` and supply `POSTURE=0.15`.

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\build_network.py:708` binds each `POST_<muscle>` input **directly to its MN's external current**. It is an input-port label, not an extra postsynaptic neuron. The distinct `POSTURE` population has its own group-weighted excitatory MN connections and extensor-RG connection. `DRIVE` supplies the locomotor network. Thus the original CAD262 path is:

**State/reference → external posture/force allocator → desired recruitment → cancelling current calculation → POST port → MN → actuator interface**, alongside the active RG/PF and sensory network.

It is not simply a physiological descending postural drive replacing PF. While the allocator/cancelling path operates, afferent/PF currents can alter network state, but their instantaneous net current onto each clamped MN is explicitly cancelled. Original hybrid results must retain that qualification.

The separately accepted mixed pressure-pair condition removes that cancelling POST current **only from the two vas_med MNs**. Their existing PF_E1, PF_E2, and IBEXC inputs then determine their recruitment through six shared weights; the other52 routes retain the external hybrid context. IBEXC's additive RG contribution is not a pure force-dependent load gate. Ideal-tension reference, pressure reference, actual pressure, and actual tension remain distinct signals.

## HuntPT routes are different

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\runner.py:1729` supplies the per-muscle standing input `6*act_stand*sf` to the same direct `POST_<muscle>` MN ports. This is a descending bias derived from the standing solution, with explicit standing/walking scaling. It does **not** invoke CAD's instantaneous synaptic-current cancellation.

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\controller.py:179` separately blends inherited POST/BAL→MN conductances by `(1-gate)`; PF expression is gated by the walking command. Existing geometry-dependent recruitment and configured learned reflex weights are separate branches. Lines206–217 supply physical timing afferents and sign-separated joint-position/velocity errors through actual `FUNC_<axis>_<side>_<channel>` inputs. These are engineering joint/task signals, not additional Ia/II receptors by assertion.

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\balance.py:12` constructs explicit `PT_BAL_*` source cells and inhibitory `PT_IN_*` relays. Their delayed/filtered COM position, velocity, and acceleration command is mapped to MN excitation/inhibition. This is a task-level postural pathway. Local muscle-spindle loops remain distinct; the PT state estimate is not relabeled a myostatic afferent.

## Antagonist Ia/II/Ib/Renshaw drawing: what exists and what does not

| Drawing element | Actual code route and qualification |
| --- | --- |
| Descending input to agonist/antagonist MN | Direct POST current ports exist. A physiological command amplitude, coactivation law, and independent effects on interneurons must be declared rather than inferred from CAD's cancellation. |
| Homonymous Ia excitation | Source builder has Ia→MN excitation. Actual conductance/activity depends on the selected model and stage; a shared diagram is not proof every candidate enables it. |
| Reciprocal Ia inhibition | Builder lines868–888 implement Ia→IaIN→antagonist MN when enabled; PF gates IaIN. Replacing that PF input with a postural-command input would be an explicit new connection design, not an existing equivalent pathway. Direct Ia inhibition is an alternate source mode. |
| GroupII feedback | Actual source includes II→central pathways, and optional IIIN/IIX routes. A generic single positive or negative II arrow is insufficient: source target, intermediate relay, and enabled gain determine the net effect. Hunt dynamic MN masking does not automatically preserve every inherited IIIN/IIX output. |
| GroupIb/load feedback | Source Ib/IBIN/IBEXC routes exist in applicable families. Actual force measurement and additive RG recruitment must be separated. A named IBEXC neuron is not by itself proof of pure load feedback. |
| Renshaw recurrent inhibition | Builder lines848–864 conditionally construct MN→RC excitation, RC→MN inhibition, and cross-pool RC inhibition. RC→IaIN inhibition is also conditional. A model without active RC populations/outputs does not acquire this subsystem merely because the reusable builder contains it. |
| Gamma/postural spindle control | Source afferent scaling knobs are not automatically a separate physiological gamma-MN circuit with intrafusal mechanics. Such a circuit is not established by the present CAD cancelling input or PT pathway. |

A meaningful new **postural-drive antagonist reflex demonstration** would explicitly disable locomotor PF outputs, specify descending drive to both MN pools and any chosen interneurons, enable and verify the selected Ia/II/Ib/RC paths, and test their individual signed currents and dynamic effects before a body trial. No such replacement is claimed implemented by this source audit. Existing electrical FSA checks establish the declared source blocks/functions; they do not silently fill missing physiological pathways.

## Bilateral source accounting

The exact ID/neighbourhood audit proves one W2L cell named `L Knee MN ext RE` is the **right** knee extensor Renshaw counterpart: its incoming right extensor-MN connection, outgoing right extensor-MN connection, right flexor-Renshaw relations, and upstream physical right-knee muscle binding agree. An explicit **logical audit/profile alias** resolves all216 active edges into108 equal L/R pairs. Native source names, IDs, and connections are preserved.

The curated rat134 has a real ankle PF-flexor→TA MN gain difference, **right4.522 versus left1.516**, and two paired edges with differing presynaptic windows/reversal fields. It is not blanket symmetric. The source remains frozen; candidate shared fitting must use explicit right-canonical parameters only for verified homologues, with unmatched or functionally different paths separately accounted. Detailed source fields and remaining exceptions are in the hashed bilateral evidence; a label alone is not anatomical or functional proof.
