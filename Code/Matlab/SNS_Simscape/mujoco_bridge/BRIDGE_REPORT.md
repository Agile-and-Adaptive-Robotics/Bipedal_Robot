# MuJoCo ↔ Simulink Bridge — Install & Prove-It Test (EB475WS4)

Task: prove (or disprove) that the official mathworks-robotics/mujoco-simulink-blockset
can drive MUSCLE actuators of our converted gait2392 model and read back muscle
length/velocity/force — before anyone builds on it.

- Date: 2026-09-12
- Machine: EB475WS4 (12 cores, 32 GB, MATLAB R2025a @ D:\Program Files\MATLAB\R2025a, MinGW-w64 gcc 8.1.0 via MATLAB support package, no VS, winget available, no cmake)
- Model: D:\Github\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml (timestep 0.005 s, 92 muscle actuators, NO sensor section)
- Rule: never modify cvt3.xml; everything new goes in this folder.

## Progress log

- [x] Created working folder + this report.
