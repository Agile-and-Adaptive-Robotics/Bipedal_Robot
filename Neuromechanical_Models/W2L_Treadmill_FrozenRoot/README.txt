W2L treadmill recording copy

Open W2L_Treadmill_FrozenRoot.aproj in AnimatLab. The original Walker_2_Layer_CPG folder is unchanged. This copy retains the original neural controller and frozen Root.

Use the Simulation Controller Play button. The drive starts at 0.5 seconds; stop/reset after the ten-second run to return the platform to its starting position. Record the actual 3D viewport with a screen recorder. An ASIM saves simulation settings, not rendered video.

The blue/checkered platform is a 4 m moving slab carried by a motorized prismatic joint, not an endless belt. The fixed anchor is below the surface so it cannot act as a second stationary floor. The checker texture makes surface travel visible.

Native AnimatSimulator test: completed 10.101 s with exit code 0. The belt moved backward in world X from -3.454 m to -4.38192 m, approximately 0.928 m over the driven interval (about 0.098 m/s). Root Y remained 0.99298 m throughout. The belt settled vertically by about 0.035 mm. Both toe contact sensors registered contact; heel sensors remained zero. Contact is intermittent, not stable stance or unsupported walking. TEST_RESULTS.json contains the measurements.

The GUI drive parameter is +0.10. The joint coordinate has the opposite sign to world-X travel. Use measured displacement/time to assess belt speed. Original neural and physical time steps are retained. The invisible default Ground structure is at y=-2 m so it cannot collide with the belt or the feet.

GUI-open validation: project loads without an error dialog. The standalone ASIM is instrumented for position/contact measurements; GUI chart forms from the original project are included. W2L_Treadmill_FrozenRoot.asim is the associated simulation; W2L_Treadmill_FrozenRoot_Standalone.asim is the tested self-contained runner.

The controller and dynamics are not tuned for successful treadmill walking. This setup is for recording supported leg motion and intermittent contact. Keep Root Freeze=True. Raise the platform only in small increments if adjusting contact, and check for initial overlap before recording.

Run the scripts in tools/ from the repository root to reproduce the copy and measurements. Generated text traces are ignored. No original model files were edited.

The matched legacy/newer-controller campaign, including air, treadmill, free-drop, partial support, Deng-style stimuli, native screenshots, motor-neuron voltage, and joint results, is in the sibling W2L_Comparison_20261030 folder. Its README.txt explains the support definition and the limits of the literature-protocol transfer.
