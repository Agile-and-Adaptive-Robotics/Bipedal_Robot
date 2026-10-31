# Newer W2L: actual 108-neuron MuJoCo stimulus graphs

**Machine: easteregg2.** Actual saved fast MuJoCo runs, not AnimatLab output or an earlier 96-neuron result. All four runs use the same 108-neuron/216-edge canonical controller.

Each figure shows both legs over a few cycles: experimental current, RG, hip PF, hip MN, and hip/knee/ankle angles. Gray overlays in the stimulus figures are the unperturbed baseline, including neural traces. Extensor is blue/solid, flexor is orange/dashed; gray baseline uses the corresponding line styles. The shaded 100 ms window marks the applied current. Voltages are in mV, current in nA, and angles in degrees.

## Full-path figures

Baseline: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_air_baseline_v2\canonical108_hip_neural_kinematics_4_8.png`. Window 4–8 s.

**RG pulse**, target `L RG ext`: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_air_RG_minus10nA_100ms_v1\canonical108_hip_neural_kinematics_3.4_7.8.png`.

**PF pulse**, target `L Hip PF ext`: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_air_PF_minus10nA_100ms_v1\canonical108_hip_neural_kinematics_3.4_7.8.png`.

**MN pulse**, target `L Hip MN ext`: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_air_MN_minus10nA_100ms_v1\canonical108_hip_neural_kinematics_3.4_7.8.png`.

All figures also have `.pdf` and `.svg` companions and a `_PROVENANCE.json` with trajectory and source hashes.

## Common test parameters

Afferented suspended air stepping, 12 s simulated duration, 2 ms physical step and 0.2 ms neural step. The canonical source already contains the doubled neural time constants; runtime neural time scaling is 1. Contacts are disabled in air mode and Root is fixed above the ground. Each test injects **−10 nA for 100 ms beginning at 4.534 s** into one left-side extensor neuron. Actual integrated test current is −1.000 nA·s. Stimulus figures show 3.4–7.8 s, covering a cycle before and several after the pulse.

The baseline and pulse traces match exactly before injection. These tests examine neural and mechanical responses at different hierarchy levels. They do not establish persistent PF memory, free balance, or ground walking. The degree of response depends on the target neuron's state when the pulse arrives.

Detailed descriptive peak-response metrics: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\canonical_12II\STIMULUS_GRAPH_METRICS.json`. No phase-reset criterion is inferred from peak differences alone.

Generator: `D:\GitHub\Bipedal_Robot\tmp\w2l_lh_style_20261031\plot_canonical108.py`. No further simulations were run by the plotting agent.
