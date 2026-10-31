# W2L Chapter 6 review package

**Staged proposal only.** Machine: **easteregg2**. No manuscript, active chapter figure directory, bibliography, or Overleaf project was changed. This package supplements the present Chapter 6 narrative with a compact test of the newer bilateral W2L controller; it does not replace the earlier matched-input NumPy/Simulink rhythm-core experiment.

## Review and insertion files

Common folder: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results`.

| File | Purpose |
|---|---|
| `build_review\REVIEW.pdf` | Four-page review: short draft and results table; compact stimulus figure; tuned ground figure; working model still. |
| `INSERT.tex` | Proposed Methods/Results/Discussion subsection. |
| `RESULTS_TABLE.tex` | Pose correction, phase correction, and parameter-fitting results; scoring intervals are explicit. |
| `FIGURES.tex` | Captions and labels for the two proposed data figures. |
| `W2L108_left_RG_PF_MN_stimuli.pdf` | Main figure: three left-leg stimulation columns, seven rows, 6.5 inches wide, 10 pt minimum text. PNG and editable SVG alongside. |
| `W2L108_tuned_ground_motion_load.pdf` | Ground evidence: bilateral hip, knee, and ankle motion, actual foot force, whole-foot clearance, and measured net support. PNG and editable SVG alongside. |
| `W2L108_tuned_ground_still.png` | Actual fitted-case replay frame at 5.96 s, showing right-foot swing; minimum right-foot clearance 113.02 mm. |
| `NUMERICAL_RESULTS.json` | Unrounded table values and exact source run folders. |
| `PROVENANCE.json` | Source and trajectory hashes, figure dimensions, font checks, trial ledger, and actual fitted parameters. |

Proposed placement: within `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\60-sns.tex`, under `\subsection{Li--Deng/W2L Walker and Controller}` and before the Gait2392 subsection. Keep the current shared-input rhythm-core verification as a separate result. The compact stimulus figure and results table are the first priorities for the main chapter. Use the ground motion/load figure to substantiate unloaded stepping; use the still in slides or as a model view, rather than treating a still as quantitative evidence. The existing full bilateral neural plots belong in Appendix D if space permits.

If selected for inclusion, copy only the selected figure files to `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\60-sns` and insert the reviewed text/float blocks into the existing chapter. No such copy or insertion has been performed. The staged wrapper `REVIEW.tex` is not a replacement chapter. Its table uses the dissertation's existing `tabularx`/`booktabs` conventions, and its prose uses existing `siunitx` commands. No new citations or bibliography changes are required for these original experimental results.

## Which controller and settings produced the results

The saved neural baseline has 108 graded neurons, 216 chemical connection entries, 12 independent homonymous type-II relays, four retained phase-related type-II channels, and 12 Hill-type muscle actuators. Both source runtime files have SHA-256 `01ece6edcbcf17c274410b23f1ce7603eab288aadb946d6b548c5102ed5c1154`. This is the newer W2L virtual walker, not the public rat controller or a Gait2392 transfer experiment.

The air-stimulus figure uses the actual `canonical108_air_baseline_v2` and its matched RG/PF/MN pulse runs. Each pulse targets the named left extensor neuron: `L RG ext`, `L Hip PF ext`, or `L Hip MN ext`; −10 nA for 100 ms beginning at 4.534 s. The graph shows the recorded substep-averaged applied current. The source neural time constants are already doubled, while the runtime clock multiplier is one. Neural step: 0.2 ms; mechanical step: 2 ms. Voltages and angles are recorded without display rescaling. Gray traces show the unperturbed matched baseline; solid/dashed gray neural traces retain the extensor/flexor identities. Data figures use Tol vibrant blue/orange with redundant solid/dashed styles, Arial, and a 10 pt text floor; each figure has an adjacent alt-text file.

The ground progression makes **two explicit changes before fitting**: select a near-flat physical initial pose by IK, and exchange four mirrored knee-PF-to-ankle-MN targets after verifying the actual actuator actions. Those corrections are distinct from Rustuna fitting. Rails are lowered 20 mm from the earlier source geometry, and the embedded passive guide apparatus has 30 mm clearance. No continuous harness force or torque is applied. The fixed 108-neuron source GUI baseline is therefore not, by itself, the complete tuned ground recipe: the saved run configuration additionally consumes the IK initial pose, four target exchanges, and fitted gains. Do not relabel an unchanged GUI replay as the optimized condition.

The editable source project is `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer.aproj`; its matching runtime is `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer_Standalone.asim`. The trained recipe is versioned separately as `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\RUSTUNA108_GROUND_WINNER.json`. The GUI source has not silently absorbed this runtime winner. Existing older 96-neuron exports and older model families are not inputs to this package.

The earlier replay `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_nearflat_ankle_crossed_v1\SUPPORTED_STEPPING_REPLAY.gif` depicts the **unfitted** pose/phase-corrected case with **42.24% BW** mean net assistance. The package still and the replay in `canonical108_ground_rustuna_winner_fresh_v1` depict the **two-gain fitted** case with **39.24% BW**. Both are real saved trials, and their captions identify the difference.

To repeat the fitted experiment on **easteregg2**, run this PowerShell command with a new, unique tag so that the saved publication case is preserved:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py' --newer --config 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\RUSTUNA108_GROUND_WINNER.json' --duration 10 --tag publication_repeat_UNIQUE
```

To regenerate only these publication figures from existing saved arrays, run `& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\make_figures.py'`. This performs no integration or fitting.

The bounded search varied only a common PF-to-MN recruitment multiplier and the mirrored toe-to-ankle-MN conductance. The selected values are **1.305863499088313** and **0.031839851437587975 µS**, respectively. All 24 recorded trials passed the physical acceptance gate. The objective combines post-1-s upward support with foot-force and relative-belt-speed penalties; it is **not a metabolic cost**. The table deliberately reports whole-trial **net** assistance instead of substituting the optimization objective or its post-transient support component. Minimum fitted knee flexion is 2.50084° left and 2.29370° right.

Qualified swing: the entire foot/toe group clears the belt by more than 1 mm for at least 80 ms, reaches at least 10 mm peak clearance, and carries less than 2.5% BW; the opposite foot bears more than 5% BW for at least half that interval. The gate requires at least three such swings per leg and at least 80% alternating transitions. Swing scoring starts at 1 s; table assistance averages the full 10 s. The ground figure shows 4–8 s without filtering away impact peaks. An isolated high force peak remains visible and is a target for further physics/controller refinement.

## Exact evidence and Appendix D plots

- Source root: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031`.
- Source run root: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer`.
- Full bilateral air plots and exact folders: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\CANONICAL108_STIMULUS_GRAPHS.md`.
- Ground plots and exact folders: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\CANONICAL108_GROUND_GRAPHS.md`.
- Fitted case: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_rustuna_winner_fresh_v1`. This folder contains `trajectory.npz`, `RESULT.json`, `ASSESSMENT.json`, `SUPPORTED_STEPPING_REPLAY.gif`, `SUPPORTED_STEPPING_STILL.png`, and `RENDER_PROVENANCE.json`.
- Consumed fitted configuration: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\RUSTUNA108_GROUND_WINNER.json`.
- Fitting ledger: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\rustuna108_ground_home_v1\TRIALS.jsonl`, `BEST.json`, and `SUMMARY.json`.
- Contact routes, verified mechanical actions, heel geometry, and limit provenance: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\CONTACT_ROUTING_AUDIT.md`.

## Checks and scope

Plots were generated from the actual saved 108-neuron arrays, with source-hash and census checks; no new simulation was run for this package. The figure generator asserts a 10 pt minimum text size at the stated 6.5-inch width. Replay captions distinguish source neural clock ×2 from runtime clock ×1 and identify the embedded physical guides. The still and both plots were visually inspected. The review PDF was compiled locally and every page rendered for inspection; the final log has no undefined references, missing figures, oversized floats, or overfull boxes. No full-dissertation build is claimed because its working source was not changed.

The package reports layer-specific responses and supported stepping in this particular W2L configuration. It does not claim completed Hunt training, a completed five-stage curriculum, a motor-neuron voltage fit, rat-to-Gait2392 transfer, or a calibrated BPA result. No new broad simulation-speed claim is introduced.

Generators: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\make_figures.py` and `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\render_trial.py`. Python: `D:\Anaconda\envs\myo\python.exe`; Matplotlib 3.7.5; MuJoCo 2.3.7; local review compiler: `D:\MiKTeX\miktex\bin\x64\pdflatex.exe`. The prose insertion is approximately 430 whitespace-delimited words including its heading and LaTeX commands.
