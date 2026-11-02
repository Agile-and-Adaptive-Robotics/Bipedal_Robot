# Gait2392 standing-controller editable schematics

Machine: **easteregg2**. These are reconstructed saved **graded biological standing candidates**, not generic default builds, walking results, or spiking networks. Persistent-sodium RG membranes and their fixed h-gate kinetics remain graded in these selected cases.

| Family | Neurons | External inputs | Active chemical pairs | Zero-conductance defined pairs |
| --- | ---: | ---: | ---: | ---: |
| s3k | 410 | 376 | 836 | 350 |
| w2lvar | 888 | 382 | 7570 | 296 |
| syn6 | 794 | 382 | 2966 | 134 |

Each family has an overview, a left RG/PF panel, and a left knee motor/reflex panel in editable SVG, PDF, and PNG. The small panels are induced circuits; omitted boundary connections remain in the complete graph. The historical light-blue, black-outline cells, inhibitory interneurons above E/F cells, Renshaw cells below MNs, orthogonal routes, open excitatory triangles, and filled inhibitory circles follow the verified LH drawing style. Labels and tables use a 10 pt floor.

The s3k and normalized w2lvar standing candidates have PF→MN and IBEXC→MN multipliers exactly zero. Their dotted PF routes are source-defined zero-conductance connections and contribute no current; they are not evidence of phasic PF motor recruitment during standing. The w2lvar antagonist IIIN→MN gains use the saved ten group-specific mean-pool factors, retaining inhibition rather than ablating it. Syn6 retains its saved weak nonzero PF drive. Every source definition, effective Gmax, reversal potential, presynaptic threshold, and saturation is stored in the graph ledger, CSV, and editor tags. SNS connection names are retained as supplied; these builders do not have unique named AnimatLab synapse types, so none are invented. T-keys identify effective numeric field tuples.

All 92 MNs and their afferent arcs remain. The physical runner force-prunes quad_fem_r/l, gem_r/l, and peri_r/l, leaving 86 biological force-active routes. This is separate from the 86-route hypothetical BPA replacement and does not depict CAD pneumatic hardware. Selected biological runs use the corrected initial pose, the recorded support apparatus, and 0.6 s pinned initialization; they are not harness-free balance claims. The source remains legacy z-up with an explicit frame distinction from the y-up CAD robot.

Physical input policy: positive lengthening velocity drives rectified Ia current; length relative to the initial standing pose drives II current with a gated tonic component; absolute source muscle tension normalized by Fmax drives Ib. Currents are nA, lengths m, velocity m/s, tension N. Exact phase/drive/gamma coefficients and normalizations are in each ledger and the hashed runner. No ground-reaction-force signal is silently substituted for tendon Ib. POST_ currents carry the saved gravity-only static-solve biases; the --reflex-stand flag disables the supraspinal BAL surrogate inputs. MN output is clip(V/E_HI,0,1), followed by native activation and muscle contraction dynamics. Voltage is relative to resting potential, not absolute −60 mV.

## Open in the block editor

Editor: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_block_editor.html`.

Save your current work, select **Templates → Blank** to create a new tab, then use **Load JSON from file…** on that tab. Open a family's `hierarchy_connectome.json` for nested layers, `left_RG_PF_connectome.json` or `left_knee_reflex_connectome.json` for the compact functional panel, or `exact_connectome.json` for the full active matrix graph. Double-click a SUB block to enter; F4 exits. Hierarchy boundary edges aggregate explicitly listed actual matrix rows; their drawing gain of 1 is not a biological conductance. Individual-node subnetworks and the exact graph preserve actual per-pair conductances. Zero source-defined PF edges are included as zero-gain entries in the motor panel; the exact graph contains only nonzero matrix entries. No browser, localStorage, or existing editor work was changed.

## Provenance and verification

Extraction executes the actual nested build function from `gait_standing_family.run` with the exact selected gains, bilateral tying, source group scales, and the two actual pre-build --reflex-stand flags, then captures the compiled network before physics. It does not run the physical simulation. Builder, parameter, and physical-runner hashes must match the saved candidate; neuron count and release-state dimensions must also match. The standing-runner wrapper has acquired later optional studies, and its current/saved hashes are recorded separately rather than falsely equated. All three reconstructions reproduce 64 actual saved post-release neural steps with **0 mV maximum error**. W2lvar uses its biological case directly; s3k and syn6 use the same graded neural configuration in saved pressure-actuation siblings for this neural-only validation, without equating biological and pressure mechanical states or promoting a failed pressure case. All saved bilateral weight projections match exactly. Muscle classification, joint-weight tables, synergy basis, and phase-profile dependencies are checked against archived snapshots as applicable. Syn6's analytical conductance helper was not archived; its current hash and that limitation are explicit in the ledger, and its effective output is independently constrained by exact saved weight projections and neural-response parity. The editor exports are functional schematics, not executable replacement runtime checkpoints; full neuron parameters and external-input rules remain in the ledgers and original source.

Every hierarchy and exact export is round-tripped through the real block editor's `specOf` JavaScript serializer in Node, without a DOM or localStorage. All endpoints, counts, signs, conductances, and full-field tags are checked. Source matrix hashes and individual checks are recorded in each family's `EXPORT_VERIFICATION.json`. Each family ledger records the exact saved result/configuration paths, their hashes, all current source hashes, flags, input adapters, source neuron parameters, and the 92/86 actuation policy. Raw experiment traces remain outside the repository.

Generator: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\export_gait_family_networks.py`.

Rebuild with `D:\Anaconda\envs\myo\python.exe` and that full script path. `--render-only` reuses the exported ledgers; `--family s3k`, `--family w2lvar`, or `--family syn6` limits the export. No production files or manuscripts are modified.

### s3k

Folder: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\s3k`.

- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\s3k\hierarchy_connectome.json`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\s3k\left_RG_PF.svg`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\s3k\left_knee_reflex.svg`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\s3k\visible_panel_synaptic_fields.pdf`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\s3k\exact_connectome.json`

Saved result: `D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\other_family_pressure\s3k_bio_prepared_v1.result.json`.

### w2lvar

Folder: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\w2lvar`.

- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\w2lvar\hierarchy_connectome.json`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\w2lvar\left_RG_PF.svg`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\w2lvar\left_knee_reflex.svg`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\w2lvar\visible_panel_synaptic_fields.pdf`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\w2lvar\exact_connectome.json`

Saved result: `D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\normalized_prepared_transfer\w2lvar_bio_graded_v1.result.json`.

### syn6

Folder: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\syn6`.

- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\syn6\hierarchy_connectome.json`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\syn6\left_RG_PF.svg`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\syn6\left_knee_reflex.svg`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\syn6\visible_panel_synaptic_fields.pdf`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\syn6\exact_connectome.json`

Saved result: `D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\pinned_preparation\syn6_prepared_bio_nominal_v1.result.json`.
