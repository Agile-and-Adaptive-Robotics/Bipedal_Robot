# Native Gait2392 model views

`opensim_human.png` and `opensim_robot.png` are captures of the stock OpenSim 4.6 Viewer in Desktop mode. `NativeScene.java` uses OpenSim's `ModelVisualizationJson` on the default model state. The human source is `mjc/gait2392_simbody/gait2392_simbody.osim`. Its saved display flags hid 88 of 92 muscles. The scene export makes all paths visible in memory, retaining each muscle's original appearance and geometry. The robot source is `gait2392_robotbody.osim`, which carries the modified right-leg paths. Both names are relative to `Solid_Models/OpenSim/Gait2392_Robotbody`.

In the OpenSim Viewer, fit the model in the viewport, expand **OpenSim Environment**, and toggle **Floor** off. In **Snapshot**, select **Make background Transparent** and capture the model. The two OpenSim PNGs retain the alpha channel. The page layout crops to the visible model.

`native_mujoco.py` uses the MuJoCo 2.3.7 OpenGL renderer on the matching `_cvt3.xml` models. The human frame starts at the model's saved default posture. The robot frame uses all 71 joint coordinates returned by the Simulink blockset at 25 ms in the motor-neuron trial. The cameras include the original MuJoCo lighting, sky, floor, model meshes, and tendon paths. Each native model supplies its own muscle-path display properties.

`robot_interface.pdf` is the native Simulink print of the working bridge model. A nonspiking neuron in `SNS_Library_R2025a` drives the actuator map. The MuJoCo block returns a sensor bus from which a Bus Selector takes the modeled muscle force. The MuJoCo block outline and name use black ink for print. The current inputs were 0 and 2 nA in 50 ms trials. The full record, `robot_bridge_readback.json.gz`, includes actuator command, length, force, and all joint positions. The renderer checks both the native force replay and lengths evaluated at the reported poses. The two records from the original and diagram-export runs are byte-identical.

To repeat the Simulink experiment, add this component folder to the MATLAB path and run `test_robot_bridge(repo, outputDirectory)` on a workstation with the compiled MuJoCo 2.3.7 blockset. The output includes the native diagram PDF, logged signals, instrumented XML, and runnable Simulink model. The XML edits resolve mesh paths and add read-only sensors.

To regenerate the two plates on easteregg2, run this from the repository root. `render_models.py` uses PyMuPDF and Pillow in the base environment to crop and lay out native images and the vector Simulink export. It invokes the MuJoCo environment for each native model frame.

```powershell
$figures = 'Documentation/Reports and Papers/Dissertation/Figures/60-sns'
$components = Join-Path $figures 'Components/model_comparisons'
& 'D:/Anaconda/python.exe' (Join-Path $components 'render_models.py') --repo . --mujoco-python 'D:/Anaconda/envs/myo/python.exe' --output $figures
```
