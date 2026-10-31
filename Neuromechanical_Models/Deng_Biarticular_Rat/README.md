# Deng biarticular rat import

Machine: **easteregg2**. One imported model, with its original source preserved.

Open this project in AnimatLab 2:

`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\Biarticular_net.aproj`

The local application is `D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatLab2.exe`. Use File → Open and select the project above. Its 21 neural diagram pages are retained, and the six `.aform` chart dependencies are included. There are no external mesh dependencies in this model. The unavailable grass texture was cleared. The GUI itself has not been exercised in this import pass; native simulator loading and execution have been verified.

## Source and scope

Author source: [Rat-SNS-design-and-tuning / Biarticular_net](https://github.com/vipamol/Rat-SNS-design-and-tuning/tree/main/Biarticular_net). This is a public **2020 precursor**, with 120 neurons and 16 physical muscles. It includes BFP, RF, and GA biarticular routes. It is **not verified as the final 2023 paper model**: several published 2023 gains differ, and the public source is not bilaterally symmetric.

`Source\Original_Biarticular_net.aproj` and `Biarticular_net_Standalone.asim` preserve the exact public files. `SOURCE_PROVENANCE.json` records their Git blob and SHA256 hashes. The local project and `Biarticular_net.asim` companion have local paths, a finite stop time, and recording adaptations. No neural connections, muscle parameters, or source joint limits were changed.

The original pelvis is frozen. This is a suspended baseline, not harness-supported ground walking. Source display distance units say Decimeters, but serialized body lengths match the project's SI `Actual` values. Multiplying them by 0.1 would be incorrect; serialized masses are grams.

## Completed native baseline

AnimatSimulator completed 5.01005 seconds. The added chart recorded 5,001 finite samples over 0–5 seconds. The left RG produces approximately 0.485-second cycles after the initial transient. The right RG has no detected peaks and is effectively constant. These are original-source observations; bilateral alternation and ground walking have not passed.

- `Native_baseline_neural_joint.png`: three left-leg cycles, with RG, hip PF, hip MN, and hip/knee/ankle coordinates, all from the same native run.
- `Native_biarticular_activity.png`: BFP, RF, and GA motor-neuron activity and the measured root y coordinate.
- `Native_bilateral_RG_diagnostic.png`: matched voltage axes reveal the right-side latch.
- `NATIVE_BASELINE_RESULTS.json`: cycle diagnostics, coordinate ranges, source scope, and raw trace hash.
- `runs\native_source_smoke\source_baseline.txt`: actual native 39-channel trace.

The joint plots retain native source signs and zero references. They are not anatomically calibrated angles against extension-neutral markers. Source limits were retained; a new knee anatomical calibration and extension gate are required before training or claiming no hyperextension.

## Reproduce

From PowerShell on easteregg2:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\prepare_native_smoke.py'
```

Run the finite, instrumented source runtime:

```powershell
Start-Process -FilePath 'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe' -ArgumentList '"D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\runs\native_source_smoke\Biarticular_source_smoke.asim"' -WorkingDirectory 'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin' -WindowStyle Hidden -Wait
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\analyze_native.py'
```

This overwrites only the designated baseline outputs. Opening/running the GUI project may regenerate its project-named companion; rerun preparation to restore the instrumented runtime.

## MuJoCo import

`Deng_Biarticular_Rat.xml` compiles in MuJoCo 2.3.7 with 16 actuators and the source frozen pelvis preserved. The source-specific import reports no unsupported body or joint types. Initial body positions agree with its source parser within 9.84e-16 m, and initial route lengths within 7.64e-17 m. This verifies the importer's internal transformations, not independent full-trajectory equivalence.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\import_mujoco.py'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\check_accelerated_port.py'
```

The mechanical XML is not yet a validated coupled neural simulation. The existing accelerated W2L controller rejects the rat source at its W2L-specific PF-to-knee-MN binding check. That rejection is saved in `ACCELERATED_PORT_PRECHECK.json`; no substitute neural circuit was used. Next work is a source-specific binding, native-trace parity testing, source-knee calibration, and a timestep convergence test before training. The prospective 2 ms MuJoCo step differs from the source's 50 microseconds, and contact models differ. No rat end-to-end speed-up has yet been measured.

Raw runs and build intermediates are ignored by Git. Curated plots, status summaries, source files, and scripts are retained. Do not label these plots as final-2023 replication, robot standing, or ground stepping results.
