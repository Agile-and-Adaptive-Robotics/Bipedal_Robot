# Rat import experiment record

Session 1, easteregg2, 2026-10-31. Separate from CHATGPT_REPORT.md and the BPA actuator experiments.

## Bounded source import

Imported the author's public Biarticular_net project and standalone runtime into this folder. Verified original Git blobs and SHA256 hashes; retained byte-exact originals. Supplied all six `.aform` dependencies. Preserved 21 XML neural diagram pages. Repaired project paths, removed an unavailable visual texture, bounded the simulation to 5.01 seconds, disabled replay recording, and added existing neural/joint/root sensors. No neural or muscle parameters were changed.

Native simulator: `D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe`. The initial attempt used an invalid `BodyPositionY` sensor name and failed loading; corrected to the supported `WorldPositionY`. The failed log remains in the ignored run folder. The corrected run completed at 5.01005 seconds and produced finite 0–5-second traces. Root y remains 0 because the author's pelvis is frozen.

After 1 second, left extensor and flexor RG traces each have eight detected peaks, median period 0.485 seconds. Right RG traces have zero peaks, with ranges below 0.0003 mV. The source baseline therefore does not establish bilateral alternation. Short-window plots cover 3.345–4.8 seconds, approximately three left cycles. Full trace and range diagnostics remain available in `NATIVE_BASELINE_RESULTS.json` and the ignored raw folder.

The initial MuJoCo mechanical conversion compiled with 16 actuators. Explicitly removed the W2L importer's free root to preserve the rat source's frozen pelvis. No body/joint type was rejected. Exact parser-relative initial position/route checks passed, but independently validated native/MuJoCo dynamics have not yet been established. Contact law and timestep differences are declared.

The existing accelerated W2L controller was tested without edits and rejected the rat source at its model-specific shared bilateral PF→knee-flexor MN lookup. That failure is recorded rather than bypassed. Coupled source-specific neural/Hill validation and speed benchmarking remain required. Source joint limits were preserved; anatomical extension-neutral knee calibration and a hard no-hyperextension training gate remain required before subsequent training.

The public source is a 2020 precursor with gains and asymmetries that differ from the 2023 paper. It must not be reported as final-2023 replication. No training sweep, harness release, ground-walking claim, or full-CAD BPA simulation was performed in this import experiment.

## Bounded II-pathway literature and wiring audit

Added `II_PATHWAY_SOURCE_AUDIT.md` and a read-only adapter/edge extractor. Primary Hunt 2017 text explicitly distinguishes central antagonist timing inhibition from local motor feedback. The public rat source contains hip-II inhibition onto antagonist RG/PF branches and LH TA-II excitation onto its own MN. Official AnimatLab source encodes Ia from series length and II from parallel length; these labels cannot substitute for tracing actual mappings and phase. The note includes a local one-DOF feedback equivalence calculation, distinguishes inferred disinhibition from demonstrated wiring, and leaves Lutz/Micklem-specific signs unverified. No neural edits or numerical runs were performed for this audit.
