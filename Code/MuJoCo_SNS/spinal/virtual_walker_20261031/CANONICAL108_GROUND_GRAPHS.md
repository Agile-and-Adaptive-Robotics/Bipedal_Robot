# Newer 108-neuron W2L: guided-ground figures

Machine: **easteregg2**. Actual saved fast MuJoCo runs; no new simulation was performed to generate these plots. Each run retains 108 neurons and 216 chemical connection entries, with the specified contact-route intervention.

Each figure shows 4–8 s on a common time axis, covering several neural cycles: heel/toe box signals, RG, hip PF, hip MN, hip/knee/ankle motion, measured vertical foot force, minimum whole-foot/toe clearance, and measured rail/guide/harness vertical assistance. The geometric sensor signal is separate from the actual foot reaction force. Support uses the recorded robot mass times 9.81 m/s²; instantaneous reactions can exceed body weight during acceleration.

| Condition | Full-path PNG |
|---|---|
| Intact baseline | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_baseline_v1\canonical108_ground_neural_load_clearance_4_8.png` |
| All contact output connections disabled | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_all_contact_outputs_off_v1\canonical108_ground_neural_load_clearance_4_8.png` |
| Contact-to-MN disabled | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_contact_MN_off_v1\canonical108_ground_neural_load_clearance_4_8.png` |
| Contact-to-PF disabled | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_contact_PF_off_v1\canonical108_ground_neural_load_clearance_4_8.png` |
| Contact-to-RG disabled | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_contact_RG_off_v1\canonical108_ground_neural_load_clearance_4_8.png` |
| Near-flat initial pose, authored PF phase | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_nearflat_source_phase_v1\canonical108_ground_neural_load_clearance_4_8.png` |
| Near-flat initial pose, four mirrored PF-to-ankle MN targets crossed | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_nearflat_ankle_crossed_v1\canonical108_ground_neural_load_clearance_4_8.png` |

Each PNG has vector `.pdf` and `.svg` counterparts and a `_PROVENANCE.json`. `RESULT.json`, `ASSESSMENT.json`, and `trajectory.npz` in the same run folder specify the actual controller, intervention, measured performance, and recorded signals. Trace appearance is not substituted for the acceptance criteria.

The baseline's heel signals are zero in this window despite toe loading. Saved-state forward kinematics and CoP show toe-supported posture, not an elevated heel sensor mistakenly interpreted as sole contact. The authored ankle MN labels also differ from their verified mechanical actions. See `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\CONTACT_ROUTING_AUDIT.md` and its linked JSON evidence before interpreting a contact-ablation difference.

The later **near-flat** conditions use the separately recorded physical initial pose rather than the original toe-down pose. In the crossed-target condition, heel and toe registration alternate and both legs repeatedly unload and clear the floor; its assessment reports 9/8 qualified swings and 100% alternation with 42.24% mean body-weight assistance. The figure also retains sharp impact-force and support peaks rather than hiding them. These are results of that new configuration, not a relabeling of the earlier baseline or a contact-gain optimization result.

Plot generator: `D:\GitHub\Bipedal_Robot\tmp\w2l_lh_style_20261031\plot_canonical108_ground.py`. No canonical GUI or physics parameter was changed by this plotting pass.

## Supported-stepping replay and slide still

The crossed-target condition above now has a replay GIF and a still selected at an unloaded right-foot swing. Both were rendered from the saved 108-neuron trajectory on **easteregg2**, using MuJoCo 2.3.7 and Pillow; no new trial was run.

- GIF: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_nearflat_ankle_crossed_v1\SUPPORTED_STEPPING_REPLAY.gif`
- Still: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_nearflat_ankle_crossed_v1\SUPPORTED_STEPPING_STILL.png`
- Source hashes, rendering times, and captions: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_nearflat_ankle_crossed_v1\RENDER_PROVENANCE.json`

The GIF spans 3.00–7.24 s (54 frames, 80 ms per frame). The still is the actual 4.84 s frame; its minimum right-foot clearance is 106.44 mm, with right-foot vertical load below 2.5% BW and opposite-foot load above 5% BW. The scene and readable captions were visually inspected. The model uses the recorded IK initial pose, four mirrored PF-to-ankle-MN target swaps, Hill actuators, embedded passive guides with 30 mm clearance, and source neural time constants multiplied by two (runtime scale one). Mean net vertical assistance across the trial is 42.24% BW. These files depict supported stepping in this apparatus.

Renderer: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\render_trial.py`. Camera and visual materials change only the replay appearance; source mechanics and the saved trajectory remain unchanged. Optional title, guidance, detail, and still-time arguments now make these annotations explicit and reusable.

## Lower hip recruitment: separate fresh repeat

The later hip-drive adjustment retains the same physical pose, guides, source clock, nonhip neural gains, and 108-neuron controller. Only the four mirrored hip PF-to-MN gains change: extensor recruitment is 60% and flexor recruitment 30% of the original fitted configuration. Effective conductances are 2.1538912554 µS (extensor) and 1.0769456277 µS (flexor), bilaterally. The additional flexor factor is 0.5 after the common hip factor of 0.6. This is a subsequent lower-effort candidate; the original 39.24% BW publication result and its assets are unchanged.

- Actual fresh replay GIF: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_hip0p6_flex0p5_fresh_v1\SUPPORTED_STEPPING_REPLAY.gif`
- Actual fresh still: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_hip0p6_flex0p5_fresh_v1\SUPPORTED_STEPPING_STILL.png`
- Saved-data before/after hip angles and muscle torques: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_hip0p6_flex0p5_fresh_v1\HIP_ANGLE_TORQUE_COMPARISON.png` (PDF/SVG alongside).
- Exact metric definitions, source hashes, and before/after values: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_hip0p6_flex0p5_fresh_v1\HIP_COMPARISON_PROVENANCE.json`.

Replay: 3.00–7.24 s, 54 frames; selected right-swing still: 4.84 s with 93.94 mm whole-foot clearance. The fresh repeat's numeric trajectory arrays match the preceding trial exactly. Its assessment retains 9/8 qualified swings and 100% alternation, with 40.50% BW mean net assistance. Post-1-s actual hip-tension RMS decreases from approximately 403 to 277 N, and positive muscle work from approximately 683 to 413 J. Hip torque RMS decreases from 35.07/59.68 to 19.50/42.77 N m (left/right). Upper-limit dwell remains 34.5%/39.9%, and right-side acceleration is not uniformly smoother. The GIF is therefore labeled **lower hip recruitment**, rather than claiming that the plateau or all abrupt motion is resolved. The still and comparison plot were visually inspected; this pass performed rendering and saved-array analysis only.
