# Editor to SNS modeler: D4 verification and replacement decision

## Latest author framing: selective tools chapter

Ben has now clarified the presentation priority. Start with background framing and the controller-transfer knowledge gap, then let Table6.1 guide the developed tools. Li-Deng/W2L is a compact validation/verification example: expectation-matched stimulus/deletion, plausible cycle timing, and a selected supported-ground trial with assistance stated. Do not turn the chapter into a catalogue of walker trials or optimization histories. A biarticular derivative with contact feedback is secondary and should appear only if it adds a distinct demonstrated capability.

The main payoff is alternative neural controllers on detailed Gait2392 mechanics, a neural controller coupled to the SOLIDWORKS robot, assembled graphical models, compiled execution, higher-fidelity muscle options, and BPA actuator dynamics. Ben mentioned a BPA hysteresis ramp: identify the implemented behavior precisely rather than relabeling a lag/saturation block as measured hysteresis. Choose one compact verification per capability; move detailed sweeps, diagnostics, and repeated controls to the appendix or supporting reports. Preserve the table order and finish with integrated robot/BPA supported posture and the next research questions.

Ben explicitly asked the dissertation editor to work with the SNS modeler on easteregg2. Submission is imminent. This is a request for modeler verification and a reviewable proposal, not permission for new tuning or replacing dissertation prose.

Current local Figure D4 is fig:layer_stimulus_validation, included as deng_neuromechanical_stimulus.pdf in 95-AppendixD.tex. Ben rejected its quality. It shows RG, hip PF, and hip MN voltages, highly oscillatory hip/knee/ankle traces, and +50 nA PF current. It was plotted here from Archive/flow_and_validation_20261028/stimulus_long_controls/{baseline.npz,PF_plus50_2cycles.npz,stimulus_results.json}. Generator maps V columns (0,1),(4,5),(12,13), plots q in degrees, and subtracts start_s from time. Caption says reduced fourteen-neuron Deng-derived circuit, suspended unilateral mechanics, no afferent return, 3.876 s stimulus as two baseline periods, and baseline period1.938 s. A copy of the exact currently displayed PDF and source context accompany this request.

Please verify:
1. Exact circuit/mechanics identities, variable mappings, angle conventions, units, and timebase. Is apparent roughly0.97 s RG periodicity versus claimed1.938 s baseline an error, or a valid distinction? Do not alter reported values by assumption.
2. Whether the many mechanical oscillations per neural cycle reflect the model behavior, a mapping/timebase mistake, or a presentation issue. Do not call them a clean neuromechanical demonstration without evidence.
3. Best EXISTING recorded alternative matching Ben's requested Deng-style short aligned RG/PF/MN, hip/knee/ankle and stimulus/baseline display. Keep original shared-RG AnimatLab, later instrumented AnimatLab, canonical108 MuJoCo, reduced14-neuron test, older96/80 W2L, human Gait2392, and CAD/BPA results distinct.
4. The newer canonical108 PF plot is visually promising, but both hips plateau near+20deg. Read-only checks found anatomical bounds about-18 to+20deg after runner offsets (-10deg left,+20deg right); raw maxima20.3185deg/20.7111deg, so this appears limit saturation. Verify intended native limits and whether this is acceptable to show. Do not retune or widen limits without Ben's instruction.

Provide one concrete recommendation, exact run/source paths and hashes, a proposed accurate caption, and a preview if an existing trace can be rendered honestly. Retain raw evidence; no silent deletions, simulation reruns, commits, or Overleaf uploads. Do not edit working60/95 sources while the editor and Ben are submitting. A neural-only version is only an option for Ben, not an approved removal of joint evidence.

Write your reply here: Documentation/Reports and Papers/Dissertation/Notes/sns_rebuild_session1_20261030/EDITOR_SNS_COORDINATION_REPLY_20261031.md
The editor will read it over SSH and handle any author-approved dissertation edits and compile. Please also identify what in your v5 build is newer than EB475WS4 sources so we avoid overwriting either party's work.
