# Single-joint reflex validation

The biological knee tests whether a circuit built in SNSGui excites the knee extensor, inhibits the flexor, and extends the knee. The pneumatic arm tests whether sensory feedback limits elbow sag when a load increases. Both use nonspiking, six-neuron circuits with six signed synapses.

## Run the saved models

From the repository root in MATLAB with Simulink, run the following commands. R2025a Update 1 was used for the recorded trials.

```matlab
repo = pwd;
demo = fullfile(repo,"Code","Matlab","SNS_Simscape","demos");
addpath(demo,fileparts(demo));
destination = tempname;
mkdir(destination);
run_gui_knee_trials(fullfile(destination,"knee"));
run_gui_arm_trials(fullfile(destination,"arm"));
verify_gui_thelen_muscles(destination);
disp(destination)
```

Use Python with NumPy, SciPy, and Matplotlib for the numerical check and PDF. Replace `<destination>` with the directory printed by MATLAB.

```text
python Code/Matlab/SNS_Simscape/demos/check_gui_reflex_results.py "<destination>"
python Code/Matlab/SNS_Simscape/demos/plot_gui_knee_reflex.py biological "<destination>/knee" "<destination>/figures" Code/Matlab/SNS_Simscape/results/reflex_validation/actuator_icon_geometry.json
```

## Rebuild the circuits

`biological_knee_circuit.json` and `pneumatic_arm_circuit.json` are the SNSGui exports. Open either file in the editor to view the saved neurons and synaptic properties. From the MATLAB setup above, the following commands reconstruct the two coupled models using the same **Build Simulink model** callback as the button. These commands replace the corresponding generated demonstrations in the demos directory. Each builder checks the six neuron masks and the source, destination, conductance, and reversal potential of each emitted synapse before attaching sensors, muscles, and mechanics.

```matlab
sns_build_biological_knee(fullfile(destination,"knee_gui"));
sns_build_gui_arm(fullfile(destination,"arm_gui"));
```

## Biological knee

The two right-side muscles, vastus medialis and biceps femoris short head, come from `Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim`. The extraction script uses OpenSim to sample path lengths and knee moment arms from 0° to 95° of flexion, with the remaining coordinates at the model defaults. `demos/gait2392_reflex_geometry.json` contains the native source parameters, geometry, reference states, source SHA-256, and OpenSim version. To regenerate it with Python and the OpenSim 4.4 bindings, run

```text
python Code/Matlab/SNS_Simscape/demos/extract_gait2392_reflex_geometry.py Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim Code/Matlab/SNS_Simscape/demos/gait2392_reflex_geometry.json
```

Thelen fiber force, force–velocity, elastic tendon force, and fixed-height pennation are calculated in the saved Simulink model. The motor command is the membrane voltage mapped from −45 to −35 mV and bounded between zero and one. Activation rises with time constant $0.010(0.5+1.5a)$ s and falls with $0.040/(0.5+1.5a)$ s. The source minimum activation is 0.01. The muscle equations agreed with OpenSim for 18 static and 36 moving-fiber states. [OpenSim muscle implementation](https://github.com/opensim-org/opensim-core/blob/4.4/OpenSim/Actuators/Thelen2003Muscle.cpp) and [activation implementation](https://github.com/opensim-org/opensim-core/blob/4.4/OpenSim/Actuators/MuscleFirstOrderActivationDynamicModel.cpp).

The knee starts at 45° flexion with both muscles at 12% activation. A constant 3.278 N m external flexion moment balances the muscle moment. The reduced joint has inertia 0.12 kg m², damping 0.6 N m s/rad, and a spring of 0.5 N m/rad referenced to the starting angle. At 20 s a 2 nA, 30 ms pulse enters the extensor Ia neuron downstream of the spindle. Relative to the matched unstimulated trial, the knee extended 0.93605° at 167 ms. The half-current trial reached 0.46766°. Removing the six sensory synapses while matching resting motor activation reduced the pulse response below $1.1\times10^{-9}$ degrees. Halving the maximum solver step preserved peak angle and timing. A separate 1 N m flexion pulse for 25 ms lengthened the extensor fibers and increased spindle current. The mechanical pulse also recruited force feedback. Comparing this response with human tendon-tap latency would require matching the mechanical stimulus and modeling neural conduction delays.

## Pneumatic arm and knee

The original arm and its GUI-built counterpart both use the two characterized 20 mm BPA blocks. The continuous ramp adds 0.0625 kg/s and drives the gravity moment directly. The applied mass reaches 0.500 kg at 8 s and 0.625 kg at 10 s. At 10 s, feedback limited sag to 1.99654° compared with 7.84641° when the six synapses were disabled. The sampled elbow angle, cup mass, two motor-neuron voltages, two activations, two forces, and flexor pressure matched between the original and GUI-generated arms in both conditions.

The separate pneumatic knee uses the Chapter 3 pressure–contraction surface on two 20 mm sleeves, each with measured resting and unloaded 620 kPa lengths of 450 and 335 mm. The reduced joint assumes two opposing 35 mm moment arms. In its nonspiking trial, the matched 2 nA pulse produced 0.56129° additional extension at 67 ms. Its spiking derivative replaces the four sensory neurons and their synapses while retaining nonspiking motor neurons. Its unstimulated angle already oscillates, so the paired waveform is phase dependent. The empirical force comparison is recorded in `pneumatic_force_checks.json`. The native pressure blocks reproduce the bounded Chapter 3 force surface to numerical precision. The original unbounded fit has 67.57 N RMSE on the 450 mm specimen and the bounded tension-only block has 65.11 N RMSE.

```matlab
run_knee_pressure_trials(fullfile(destination,"pressure"),"all");
verify_knee_pressure_bench(repo,fullfile(destination,"pressure_bench"));
```

```text
python Code/Matlab/SNS_Simscape/demos/plot_gui_knee_reflex.py pneumatic "<destination>/pressure" "<destination>/figures" Code/Matlab/SNS_Simscape/results/reflex_validation/actuator_icon_geometry.json
```

The actuator drawings in both plots are the coordinates returned by `sns_icon_geom('BIO')` and `sns_icon_geom('BPA')`. The figure builder uses newly drawn vector paths and the recorded simulation signals.

## Standard palette and property editor

The stock `BioMuscle` block uses a rigid tendon and a signed rate of change of musculotendon length. Positive velocity denotes lengthening. Its force–velocity approximation divides fiber speed by optimal fiber length. At 50% activation, optimal fiber length, and a maximum isometric force of 1,000 N, speeds of −0.02, 0, +0.02, and +0.04 m/s produced 453.70, 500.00, 626.42, and 672.93 N. The former eccentric sign produced 156.34 and 0 N at the two positive speeds. A check at +0.25 m/s also verified that both calculated branches stay finite. [OpenSim velocity convention](https://github.com/opensim-org/opensim-core/blob/4.4/OpenSim/Simulation/Model/Muscle.h#L243-L247)

The nonspiking neuron and synapse each have a ramp onset and a voltage span in mV. A 10 mV span beginning at −45 mV gives normalized outputs of 0, 0.5, and 1 at −45, −40, and −35 mV. The property editor regression changes both values on a neuron and a synapse, saves and reloads the circuit, and checks the generated block masks. Run both checks in MATLAB from the repository root.

```matlab
palette = fullfile(pwd,"Code","Matlab","SNS_Simscape");
addpath(palette,fullfile(palette,"gui"));
sns_test_actuators
test_sns_gui_voltage_edits
```

`palette_checks.json` records both shipped libraries, a fresh build, and the source audit of 67 saved Simulink files. The reported biological knee uses its own elastic-tendon muscles and the reported pneumatic arm uses `BPA_20mm`. All 13 biological knee and pneumatic arm signal records matched the pre-correction traces sample for sample.
