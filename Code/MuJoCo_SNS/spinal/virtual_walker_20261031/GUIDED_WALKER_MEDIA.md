# Guided W2L stepping media — Session 1

Machine: **easteregg2**. Created 31 October 2026.

## Files for dissertation and defense editors

| Full path | Role |
|---|---|
| `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\guides_long_v1_t4_c0.03_h0_p5\SUPPORTED_STEPPING_REPLAY.gif` | Saved-motion animation, 4–9.5 s of the 20 s trial. |
| `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\guides_long_v1_t4_c0.03_h0_p5\SUPPORTED_STEPPING_STILL.png` | Still from the same animation, t = 6.72 s. |
| `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\guides_long_v1_t4_c0.03_h0_p5\trajectory.npz` | Full saved neural, muscle, mechanical, contact, and support trajectory. |

Program: MuJoCo **2.3.7**, with source-specific IntegrateFireSim equations and
Hill muscles, accelerated using Numba. Python executable:
`D:\Anaconda\envs\myo\python.exe`.
Runner and renderer:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py`
and
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\render_trial.py`.

## Configuration and test

Original W2L circuit, **80 neurons, 150 chemical connections, and 12 Hill muscle
actuators**. Neural time constants scaled by **4**; physics step **2 ms**, neural
integration step **0.2 ms**. Original horizontal support rails lowered **20 mm**;
passive guide carriage with **30 mm lateral/fore-aft clearance**. Continuous
sagittal harness gain **0**. Phase-triggered **5 N·m yaw perturbation for 100 ms**.
Treadmill motor explicitly enabled through 20 s. Anatomical frame: x anterior,
y superior, and z right.

Test: cadence/guide-clearance comparison with measured ground contact and external
support. Over 20 s, this condition produced **11/11 qualified swings**, **100%
alternating transitions**, median RG period **1.760 s**, and approximately
**28.2% double support** after startup. Mean belt-relative strides were about
**146/141 mm**. Mean measured vertical assistance was **85.4% body weight**.
The knees remained flexed. This is the virtual W2L model with Hill muscles;
the CAD robot and BPA-control experiments have separate records.

Suggested caption: “Rail-supported W2L treadmill stepping with passive guides,
using the original neural circuit at a 1.76 s cycle period. Guide clearance was
30 mm, and measured mean vertical assistance was 85.4% body weight.”

Full settings and assessment are in PARAMETERS.json, RESULT.json, and
GUIDE_ASSESSMENT.json beside the media. The comparison summary is
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\GUIDE_LONG_COMPARISON_V1.json`.
Raw run folders are gitignored and remain local on easteregg2.

A matching short results plot is
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\guides_long_v1_t4_c0.03_h0_p5\SHORT_NEURAL_MECHANICAL.png`.
It shows RG, knee PF, knee MN, bilateral hip/knee/ankle motion, and foot vertical
forces over the same 4–9.5 s interval (approximately three cycles).

## Native AnimatLab companion

`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_GuidedT2_20261031\W2L_VirtualWalker_GuidedT2.aproj`
is the editable native companion with **clock scale 2 and 10 mm guides**. Its
configuration differs from the clock-4, 30 mm MuJoCo animation above. Matching
standalone exports and the native experiment record are provided in the model
folders. No claim of identical contact dynamics between engines is implied.
