# Chapter 6: controller-transfer evidence gaps

Session 1, **easteregg2**, 31 October 2026. Read-only audit of the current ProofFinal Chapter 6, Appendix D, applied SNS revision notes, and saved simulation evidence. No manuscript, bibliography, figure inclusion, or Overleaf source was changed.

## The missing demonstration

The chapter already explains the intended progression: a published controller is transferred to different mechanics, then evaluated on Gait2392 and a BPA-actuated robot. What is still needed is a **traceable transfer of the recovered rat controller's eight muscle-output groups per leg to specified Gait2392 or CAD actuators, followed by a recorded closed-loop trial**. The new 108-neuron W2L results supply useful evidence for controller construction, accelerated execution, and hierarchical stimulus experiments. They do not fill that particular transfer step because their mechanical model is the W2L virtual walker.

Current line references below are from:

- [Chapter 6 source](<D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/60-sns.tex>)
- [Appendix D source](<D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex>)
- [Applied revision record](<D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/sns_rebuild_session1_20261030/APPLIED_20261031.md>)

## What the text and records already establish

| Topic | Current text | Evidence and remaining distinction |
|---|---|---|
| Transfer is the chapter's purpose | Chapter 6 lines 4, 44–48, and Table 6.1 lines 59–66 | The framing already matches Ben's goal. Retain this structure and add a concrete transfer example beneath it. |
| Original and newer W2L lineage | Chapter 6 lines 8–12, 27–40 | The original Li-mechanics/Deng-derived controller and the later bilateral RG/contact model are identified. The newer 108-neuron model belongs to this later lineage, not the public rat source. |
| Cross-program neural implementation | Chapter 6 lines 88–101; Appendix D lines 54–58 | Existing NumPy/Simulink comparisons concern specified neural cores. They are separate from moving the recovered rat controller to human or CAD mechanics. |
| Gait2392 mechanics and alternative controllers | Chapter 6 lines 114–142; Appendix D lines 131–176 | Conversion repairs, force checks, network-organization comparisons, and an IK overlay are already present. The `w2lvar` organization is not, by its name, proof of a faithful transfer of the recovered eight-group rat controller. |
| Earlier layer-stimulus evidence | Chapter 6 lines 249–255; Appendix D lines 79–112 | The recorded experiment uses a **14-neuron reduced Deng-derived circuit**, suspended unilateral W2L mechanics, and no afferent return. Appendix D line 81 explicitly identifies those mechanics. Preserve its results as that experiment. |
| Full Gait2392 layer-stimulus series | Appendix D line 176 | The text explicitly says this series was not completed. New W2L stimulus figures must not silently replace that statement with a Gait2392 result. |
| CAD route integration | Chapter 6 lines 267–273; Appendix D lines 304–335 | The author-confirmed SNS-Simscape/BPA standing example and the separate MuJoCo CAD reference-tension experiments have different control/actuator definitions. Appendix D line 331 records MN compensation that cancels PF/reflex input to follow a force allocator. A rat-controller transfer must give the transferred circuit effective actuator authority rather than retain that cancellation unnoticed. |

## New evidence available now

**Newer W2L, 108 neurons and 216 chemical edges.** Afferented suspended MuJoCo runs use the same saved controller and source hash `01ece6edcbcf17c274410b23f1ce7603eab288aadb946d6b548c5102ed5c1154`. The baseline and three separate −10 nA, 100 ms pulses to left RG-extensor, hip-PF-extensor, and hip-MN-extensor neurons have bit-identical pre-pulse neural trajectories. Each recorded injected current integrates to −1.000 nA·s. Figures show both legs, actual current, RG/PF/MN voltages, hip/knee/ankle angles, and unperturbed overlays over a few cycles.

The PF pulse changes local PF/MN activity and joint motion with a much smaller RG change than direct RG stimulation. Descriptive peak left-RG differences over onset–7.8 s are 12.751, 0.193, and 0.0086 mV for RG, PF, and MN pulses, respectively. These measurements help demonstrate the experimental tools; peak differences alone do not classify persistent PF memory or a phase reset. The modest MN-pulse effect is a result at the tested circuit state, not evidence that the MN interface is ineffective.

Full paths, plots, parameters, hashes, and metrics are in [the 108-neuron stimulus-graph index](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/virtual_walker_20261031/CANONICAL108_STIMULUS_GRAPHS.md>). Those figures use virtual-walker mechanics, not rat, Gait2392, or CAD mechanics.

**Recovered public biarticular rat source.** The imported author project has 120 neurons, 16 muscles, and BFA, IP, VA, SO, TA, BFP, RF, and GA identities. Its source-faithful native baseline produces left-side cycling while the right RG remains nearly constant. The MuJoCo mechanical import compiles, but the existing W2L-specific accelerated controller initially rejects its bindings. This source is a public 2020 precursor; numerical differences from the 2023 paper are recorded. Source recovery and model loading are useful deliverables, distinct from a completed coupled rat-to-Gait2392 transfer.

Source scope and native/MuJoCo results: [rat README](<D:/GitHub/Bipedal_Robot/Neuromechanical_Models/Deng_Biarticular_Rat/README.md>) and [separate experiment record](<D:/GitHub/Bipedal_Robot/Neuromechanical_Models/Deng_Biarticular_Rat/EXPERIMENT_RECORD.md>). The transfer agent is developing the actual eight-group mapping; this audit does not invent substitutions from muscle labels.

## Selective placement in the chapter

1. **Table 6.1 / W2L subsection:** retain the table order. After the W2L implementation account around Chapter 6 lines 88–95, add a short, separately identified paragraph about the newer bilateral W2L stimulus experiments. Use **one representative PF-pulse figure**, with actual baseline overlays and a common short time axis. That one figure already includes an unperturbed comparison, so a second full-page baseline figure adds little.
2. **Appendix D, after the current layer-stimulus section (lines 79–112):** add a clearly separate subsection for the **108-neuron bilateral W2L** experiment. A compact RG-versus-MN comparison, plus a small protocol/response table, can hold the complementary results. Retain the existing reduced-circuit result and identify its different circuit/mechanics once. All four full-resolution plots can remain in the repository and defense materials.
3. **Gait2392 subsection, after Chapter 6 line 134:** place the eventual **rat-controller-to-Gait2392 result here**, with one identifiable muscle-group map and a short neural/mechanical response. Until the new run exists, this is a precise evidence gap, not an additional success claim. Do not change the meaning of the earlier production-controller IK overlay at lines 136–142.
4. **Literature-circuit subsection, lines 245–255:** distinguish the reduced Deng-derived demonstration from recovery of the public biarticular rat source. Cite the 2023 paper for its published findings; identify the recovered source version separately. Do not label the 108-neuron W2L implementation as the original rat network.
5. **CAD closing example, Chapter 6 line 267:** retain the brief tools-integration segue. If the transfer reaches CAD mechanics, name its output mapping, actuator law, and support condition. Use the matching recording, rather than substituting the earlier hybrid-allocator GIF.

This would add one main-chapter response figure and, at most, one appendix comparison figure plus a small table. It follows the advisor's request for selected demonstrations and Ben's preference for aligned neural/mechanical traces, rather than a separate plot for every interneuron.

## Deliverable within a few hours versus measurements still required

| Deliverable | Bounded evidence required |
|---|---|
| Complete the existing W2L demonstration now | Already available: source-specific controller identity, baseline, matched RG/PF/MN interventions, actual currents, few-cycle figures, and reproducible raw data. These can strengthen the tools result immediately. |
| Document a rat-to-human transfer configuration | Inventory the actual source circuit and eight outputs per side; map each to named target actuators and its physical joint actions, including biarticular targets. Document unmatched muscles, aggregation, retained source parameters, sensory normalization, and explicit engineering changes. The mapping agent owns this work. |
| First bounded Gait2392 trial | Use the identified circuit with direct neural activation/actuator coupling, explicit sagittal restrictions or harness support, valid contact initialization, and no knee hyperextension. Record RG/PF/MN activity, group activations, joint angles, root height, support reactions, and contact/unloading measures. A finite suspended trial can demonstrate coupling even if the next ground condition needs training. |
| Minimal transfer comparison | Compare the selected target trial with its source or neural-only baseline under the same documented input protocol. Repeat at least one saved target condition independently before reporting a numerical result. |
| Later performance claims | Alternating loaded steps, progressive support reduction, dynamic actuator agreement, sensitivity, and IK agreement require their corresponding measurements. They cannot be inferred from an imported model or from the W2L response figures. |

The strongest short-deadline chapter improvement is therefore **a clearly identified W2L hierarchy experiment now, followed by the first measured, explicitly mapped rat-to-Gait2392 coupling result if the current transfer run completes**. The manuscript already supplies the motivation; the new evidence should make the transfer path concrete.
