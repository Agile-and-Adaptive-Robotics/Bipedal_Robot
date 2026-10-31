W2L comparison and recording projects

104 native AnimatSimulator runs completed: eight ten-second baselines and 96 perturbation trials. All recorded samples are finite. No original walker project, dissertation result, or controller gain was changed.

OPEN IN ANIMATLAB
recording_projects/legacy_treadmill/legacy_treadmill.aproj: original controller, frozen Root, moving platform.
recording_projects/legacy_drop/legacy_drop.aproj: original controller, free Root, stationary platform.
recording_projects/legacy_partial50/legacy_partial50.aproj: original controller, free Root, stationary platform, 205.4031 N upward force at Root.
Equivalent bilateral_* folders use the saved 96-neuron Walker_2_Layer_CPG_BilateralRG controller, with four Type II adapters, four contact-count adapters, and separate L/R rhythm generators. The physical Root subtree and initial pose match the legacy model. Neither controller was retuned.

Select the Sim Window tab and press Play. Each project runs for 10.1 s. Stop/reset to restore the initial pose before recording. Record the actual viewport with PowerPoint Screen Recording or another recorder; ASIM files contain simulation settings, not rendered video. AnimatLab camera controls: left-drag rotates; right-drag zooms. Avoid Ctrl+Shift+right-drag, which applies a mouse-spring force. Official instructions: https://animatlab.com/Help/Documentation/Simulator

CONDITIONS
air: Root frozen; both platform and default Ground structure moved below the legs. All four contact sensors remain zero for both controllers.
treadmill: Root frozen; 4 m slab travels backward in world X. Drive parameter +0.10 begins at 0.5 s. Measured travel is 0.928 m over the driven interval, approximately 0.098 m/s. This is a translating slab, not an endless belt. Contact is intermittent.
drop: Root free, no upward support, stationary platform.
partial50: Root free; constant positive-Y support equal to 50% of the model's declared 41.919 kg weight at gravity 9.8 m/s^2. This force supplies no balance constraint or rotational stabilization.

OBSERVATIONS
Both controllers fall in drop and partial50. For the legacy controller, Root first crosses below 0.5 m at 1.103 s without support and 1.790 s with partial support. Corresponding bilateral-controller times are 0.999 s and 2.077 s. These are operational fall markers, not gait-quality scores. The bilateral partially supported model subsequently moves below belt height; the lower fixed anchor remains part of the scene.

The measured L RG ext period is 0.443 s in legacy air and treadmill trials, 0.493 s in bilateral air, and 0.4645 s in bilateral treadmill. The comparison shows a change in the saved newer controller's rhythm with mechanical contact; it does not isolate the contribution of each added feedback pathway. No unsupported walking is demonstrated.

STIMULI
Deng 2019 protocol transfer: -10 nA to L Hip PF ext or L RG ext at 2.0-2.1 s, separately, in all four conditions and both controllers (16 trials).
Source: Neuromechanical Model of Rat Hindlimb Walking with Two-Layer CPGs, Biomimetics 4(1):21, https://www.mdpi.com/2313-7673/4/1/21

Deng 2022 protocol transfer: +/-0.5 nA to L RG flx IN for 5% of the measured baseline period, at ten phases 0.05 through 0.95, in air and treadmill conditions for both controllers (80 trials). Phase zero is the L RG ext adaptive-midrange upward crossing. phase_response.csv reports first- and third-cycle shifts relative to the unperturbed trial, positive for an advance. These W2L trials do not reproduce the original paper's rat model, frequency retuning, synaptic-strength sweep, or full experimental design.
Source: Biomechanical and Sensory Feedback Regularize the Behavior of Different Locomotor Central Pattern Generators, Biomimetics 7(4):226, https://www.mdpi.com/2313-7673/7/4/226

FILES
output/W2L_Native_AnimatLab_Results.pdf: 27-page review with native captures and all baseline MN/joint and perturbation figures.
figures/: PNG and vector PDF motor-neuron voltage, joint-angle, support, and perturbation plots. Joint angles retain AnimatLab's native convention.
screenshots/: actual AnimatLab captures. legacy_drop_after_release.png shows t=5.448 s; legacy_partial50_after_release.png shows t=2.226 s. The viewport is native 3D rendering.
legacy_drop_clean.png is the clean slide-19 capture at t=2.141 s. legacy_frozen_treadmill_running.png shows t=3.518 s.
legacy/ and bilateral/: each case preserves trial.asim, protocol.json, runner.log, and losslessly compressed signals.npz with channel labels.
baseline_summary.json, phase_plan.json, phase_response.csv, COMPLETE.json, VERIFICATION.json: numerical summary, definitions, and checks.
protocol_sources/: exact input snapshots corresponding to COMPLETE.json hashes.

REPRODUCTION
From the repository root, use the scripts in Neuromechanical_Models/W2L_Treadmill_FrozenRoot/tools/: build_treadmill.py, run_comparison.py, plot_comparison.py, and build_recording_projects.py. On easteregg2 use D:\Anaconda\envs\myo\python.exe. run_comparison.py overwrites this campaign's generated cases. The original source folders remain unchanged.

An initial diagnostic pass exposed an invisible default Ground plane at y=0. That pass was discarded. The reported campaign moves Ground to y=-2 m, and the air platform below the legs, before rerunning every case. No pre-correction result is used here. Generated case data, figures, screenshots, and recording-project copies are ignored by Git; do not add a large ZIP bundle.
