# Gait2392: Hunt four-step training and five-stage validation

Machine: easteregg2. Repository: `D:\GitHub\Bipedal_Robot`.
Python: `D:\Anaconda\envs\myo\python.exe`.

This is the working experiment record, separate from `CHATGPT_REPORT.md`.
Raw arrays, trial configurations, logs, and databases belong in the ignored
`runs/` directory. Existing production winners and the completed exploratory
campaign remain separate.

## Current sequence

1. Audit full IK sources and model constraints for MTP, lumbar x/y/z, hip
   x/y, and calcaneal/subtalar x-related motion before selecting PF targets.
2. Prepare model-consistent kinematics, forces, expected afferents, and MN targets.
3. Fit undriven intrinsic rhythm, expected-afferent transition timing,
   PF-to-MN peaks, and Ia/Ib refinement as four separate functional steps.
4. Run deafferented air, afferented air, standing, air-pattern, and ground
   stages. Failed stage gates remain failures; later diagnostic tests do not
   certify progression.

## Initial IK finding

The 43 cached gait references contain only bilateral hip flexion, knee, and
ankle curves plus summary statistics. Their intake code deliberately discards
the other coordinates. These caches cannot answer whether non-sagittal,
lumbar, or MTP patterns exist. Full original MOT/STO files and model locking
must be inspected; omission from a cache is not evidence of absent movement.

Kinematic independence can justify testing a separate temporal control pattern.
It cannot by itself identify a biological PF population or prove that one
dedicated oscillator per coordinate is necessary. The Gait2392 subtalar joint
is oblique; its coordinate and anatomical-x component will be distinguished.
