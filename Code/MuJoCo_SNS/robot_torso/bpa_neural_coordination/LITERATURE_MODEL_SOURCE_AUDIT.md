# Deng model-source and ankle-synergy audit

Read-only source investigation on **easteregg2**, 2026-10-31. No model, muscle, or neural connection was changed. Inspection scripts/cache are in this folder's gitignored `runs\deng_source_audit_20261031\`.

## Ankle action: do not infer it from “EXT” alone

The original ModelDB 2019 project's resolved neural graph does contain **LH_K&A PF E → LH_AnkleZ MN E**, and the corresponding F → F connection. The E motor neuron reaches a NodeToPhysicalAdapter and the physical muscle `LH_AnkleZ-Ext` (body ID `ab2d627c-ba03-44b5-9d25-88ee77acddc3`). Thus, matching the labels reproduces the original neural mapping; it does not establish that a retargeted W2L route has the correct anatomical torque.

The independent geometry agent's read-only rat-project probe reports opposite signed native-axis moment arms for the original ankle antagonists: approximately **+4.9 mm for Ext and −10.5 mm for Flx** on both sides. Its parser has not been verified against native anatomical markers for this original 3D rat pose, so this is **not yet a verified dorsiflexion/plantarflexion classification**. See `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\DENG2019_TORQUE_PROBE.json`. A sign of a joint coordinate, a muscle's name, and its stance function must be checked separately.

[Deng 2019, Sections 3.2 and 4.2](https://pmc.ncbi.nlm.nih.gov/articles/PMC6477610/) motivates the shared PF through knee/ankle kinematics. Importantly, its Figure 13 discussion acknowledges an **ankle stance-phase mismatch with animal data**. The published synergy is therefore not evidence that every similarly named knee/ankle channel is a physiologically correct stance-extensor synergy after mechanical adaptation. Keep the source-faithful baseline; a phase swap can be tested only as an explicitly labeled adaptation, with physical torque and contact consequences measured.

## What the 2019 archive actually contains

Local archive: `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\dev\deng2019_modeldb_264844.zip`.

SHA256: `6b393118c453979255956d9c6eb03940250c1d0c6d04b69a89ded330449ae64d`.

The archive contains 36 entries: six main/backup **.aproj** files, charts, MATLAB comparison files, joint recordings, and a dissertation PDF. Its main project is `K_A wirh RG  march 21/K_A synergy.aproj`, with hip, knee, and ankle antagonist pairs, **12 Hill muscles total**. Four projects have 12 muscles; the two `K_A synergy I.aproj` backups have 16 because they also include left-hip X/Y antagonist pairs. Those extra hip muscles are not BFP/RF/GA biarticular groups. No named BFP/RF/GA successor or standalone ASIM is present in this archive. The complete per-project inventory is in `ARCHIVE_INVENTORY.json`.

Primary distribution: [ModelDB 264844](https://modeldb.science/264844) and its [public repository](https://github.com/ModelDBRepository/264844). The archive is a 2019-model source, not a hidden release of the 2023 model.

## Actual public biarticular AnimatLab sources found

The older published URL `github.com/vipamol/Rat-Waking` redirects to [vipamol/Rat-Walking](https://github.com/vipamol/Rat-Walking). It includes `Biarticular.aproj` and `Biarticular_Standalone.asim`, but its [README](https://github.com/vipamol/Rat-Walking/blob/master/README.md) describes incomplete biomechanical/passive-tension work and points to a separate neural project.

That neural successor is public: **[vipamol/Rat-SNS-design-and-tuning](https://github.com/vipamol/Rat-SNS-design-and-tuning)**. Its [Biarticular_net folder](https://github.com/vipamol/Rat-SNS-design-and-tuning/tree/main/Biarticular_net) provides:

| Actual public file | Git blob SHA |
|---|---|
| [Biarticular_net.aproj](https://github.com/vipamol/Rat-SNS-design-and-tuning/blob/main/Biarticular_net/Biarticular_net.aproj) | `e0b3123513594407f6f09f1e87773f59d75a5dfb` |
| [Biarticular_net_Standalone.asim](https://github.com/vipamol/Rat-SNS-design-and-tuning/blob/main/Biarticular_net/Biarticular_net_Standalone.asim) | `afdc829ced0ead90ebf7f0552da6edc6ee504219` |

Both files were fetched through the GitHub connector and cached without running them at:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\runs\deng_source_audit_20261031\`

`PUBLIC_SOURCE_PROVENANCE.json` verifies both cached files against the exact Git blob SHAs above. The fetched model files are gitignored inspection copies, not changes to an active simulation.

The standalone file has **16 physical muscles (eight per leg), 120 neurons**, and separate hip, knee, and ankle PF groups. Muscle identities include BFA, IP, VA, SO, TA, BFP, RF, and GA. It explicitly labels **SO MN EXT** and **TA MN FLX** and connects ankle PF EXT to SO MN EXT and ankle PF FLX to TA MN FLX. This is concrete public source evidence for the later model's named ankle pools, not a validation of all its physical torque directions.

There is no need to invent an extra ankle muscle in the requested monoarticular W2L baseline. The separate biarticular source really does contain GA; that is an author-provided route and belongs to a separately identified biarticular model.

## Public precursor versus exact 2023 reproduction

The retrieved public files have a GitHub modification date of **2020-11-23**. Their settings do **not** match the supplied 2023 manuscript's Table 2. For example:

| Connection | Public predecessor | Supplied 2023 Table 2 |
|---|---|---|
| Hip/Knee → BFP integration | common “Biar connect” type, 0.5/0.5 µS | 0.967/0.365 µS |
| PF → SO MN | 4.054 µS | 8 µS |
| PF → TA MN | left 1.516 µS; right 4.522 µS | 2 µS |

Some public connections cross side labels, and direct joint-PF-to-biarticular-MN connections coexist with intermediate nodes. Do not assume this older public project is a clean final implementation of the later circuit.

The [2023 publisher page](https://link.springer.com/chapter/10.1007/978-3-031-39504-8_2) does not expose a supplementary model download in the inspected page. No dedicated model-repository URL was found in the supplied author manuscript. The actual public precursor is available; **an exact final-2023 parameterized model release has not been located or verified**.

The supplied manuscript at `D:\GitHub\_zcode_backups\lit_20261026\kaiyu2023.docx` identifies five monoarticular and three biarticular groups. Its Figure 1 describes separate joint PF networks and integration for BFP/RF; GA combines knee-flexor and ankle-extensor signals. These are source-specific requirements for a 2023-inspired reconstruction, not permission to silently change the monoarticular baseline.

## Recommended next checks

1. Preserve both requested W2L sources and their existing signed actuator mappings as the replication baseline.
2. Verify the original rat and adapted W2L anatomical ankle directions with explicit shank/foot markers or a native export. Compare positive tensile torque and virtual work independently of MN names.
3. If testing an ankle phase swap, save it as a distinct adaptation and compare the same neural drive, harness, contact protocol, and knee-extension guards. Record both successful unloading and any degradation of stance support.
4. Use the public biarticular project/body as a documented precursor, then reconcile its graph and parameters against the supplied 2023 manuscript before calling it a reproduction. Do not silently copy the old asymmetries or bypass intermediate integration.

All source findings were sent to the two virtual-walker agents. No literature claim establishes current W2L ground-walking success, and no new simulation outcome was generated in this audit.
