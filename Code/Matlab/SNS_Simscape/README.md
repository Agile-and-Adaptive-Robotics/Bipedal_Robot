# SNS-Simscape

Use Simulink to construct neural controllers from the SNS block library. The examples connect neurons and synapses to actuator and joint models. Mechanical examples also require Simscape and the relevant physical-model libraries.

In MATLAB R2025a, add this directory and its `demos` and `simscape_sources` subdirectories to the path. Open `SNS_Library_R2025a.slx` to inspect the available blocks or `KneeReflexDemo_R2025a.slx` to inspect the reflex example. The unsuffixed library is supplied for R2025b.

The optional [MuJoCo bridge](mujoco_bridge/) includes the Simulink blockset source and its license. The Windows build procedure in `mujoco_bridge/matlab/laptop_build_bridge.m` uses MinGW and MuJoCo 2.3.7.
