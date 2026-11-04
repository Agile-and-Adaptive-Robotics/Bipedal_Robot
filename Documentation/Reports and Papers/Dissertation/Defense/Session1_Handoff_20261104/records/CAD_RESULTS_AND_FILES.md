# New CAD stepping GIF and full paths

Machine: **easteregg2**. Selected: **77.5% measured BW support**,28 s;7R/6L qualifying clear swings in16–28 s;minimum5.10degrees knee flexion. Fresh/reference arrays exactly match. Local checks: 5/8 completed, complete.

| Asset | Full path |
|---|---|
| GIF | [cad3dof_supported.gif](<D:/Users/Ben Bolen/CAD3DOF_20261102/contact_handover_20261103/media/clearance77p5_damped/cad3dof_supported.gif>) |
| Print still | [cad3dof_supported_still.png](<D:/Users/Ben Bolen/CAD3DOF_20261102/contact_handover_20261103/media/clearance77p5_damped/cad3dof_supported_still.png>) |
| Neural and mechanics PDF | [NEURAL_AND_MECHANICS.pdf](<D:/Users/Ben Bolen/CAD3DOF_20261102/contact_handover_20261103/media/clearance77p5_damped/NEURAL_AND_MECHANICS.pdf>) |
| Model | [toe_base_passive_nominal.xml](<D:/Users/Ben Bolen/CAD3DOF_20261102/anatomical_revision/models/toe_base_passive_nominal.xml>) |
| Routes | [routes_toe_base.json](<D:/Users/Ben Bolen/CAD3DOF_20261102/anatomical_revision/models/routes_toe_base.json>) |
| Selected options and source hashes | [clearance77p5_damped.json](<D:/Users/Ben Bolen/CAD3DOF_20261102/contact_handover_20261103/runs/clearance77p5_damped.json>) |
| Trace | [clearance77p5_damped.npz](<D:/Users/Ben Bolen/CAD3DOF_20261102/contact_handover_20261103/runs/clearance77p5_damped.npz>) |
| Executed controller | [control_clearance.py](<D:/Users/Ben Bolen/CAD3DOF_20261102/contact_handover_20261103/source/control_clearance.py>) |
| Frozen runner | [a54c0a7d0f23a497a07fa3b067937b5fd14b0a8b19a0ad3a140aa4db70266ae4.py](<D:/Users/Ben Bolen/CAD3DOF_20261102/source_snapshots/a54c0a7d0f23a497a07fa3b067937b5fd14b0a8b19a0ad3a140aa4db70266ae4.py>) |
| Experiment record | [EXPERIMENT_RECORD.md](<D:/Users/Ben Bolen/CAD3DOF_20261102/contact_handover_20261103/EXPERIMENT_RECORD.md>) |
| Final metrics | [FINAL_CHECKS.json](<D:/Users/Ben Bolen/CAD3DOF_20261102/contact_handover_20261103/FINAL_CHECKS.json>) |
| Repeat and local verification | [FINAL_VERIFICATION.json](<D:/Users/Ben Bolen/CAD3DOF_20261102/contact_handover_20261103/FINAL_VERIFICATION.json>) |
| GIF provenance | [asset_provenance.json](<D:/Users/Ben Bolen/CAD3DOF_20261102/contact_handover_20261103/media/clearance77p5_damped/asset_provenance.json>) |

To repeat, copy source/final_repeat_cases.json to a new filename, change its case names to new tags, and run source/run_cases.py with the copied spec. The driver preserves existing result JSONs. All consumed parameters are in the selected JSON; source and model hashes identify the exact implementation. The controller is hybrid SNS plus posture/clearance allocation; no external foot lift is used.
