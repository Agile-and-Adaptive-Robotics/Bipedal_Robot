# Rat-to-Gait2392 named-muscle transfer audit

Machine: easteregg2. Read-only inventory, no controller transfer or retraining performed.

Status: rat-to-Gait/CAD transfer is paused at the user's instruction. Current training work concerns the newer108 W2L ground model. The correspondence below remains a recorded preparation, not an executed transfer.

## Actual target files

- Human OpenSim: `D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_simbody.osim`
- Existing human MuJoCo: `D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml`
- Robot-body OpenSim: `D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_robotbody.osim`
- Generated repaired robot OpenSim: `D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_robot.osim`

These four inspected files each contain 92 named muscles/actuators. That is a file-specific inventory, not proof of 92 valid CAD/BPA routes or the user's intended 27-per-side reduced route set. Do not conflate the human simbody file with robotbody, a repaired generated model, or a different reduced assembly.

Source: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\Biarticular_net_Standalone.asim`. All 16 physical muscle-output adapters were traced to their actual MN names. Each preserves its source voltage through a polynomial adapter with linear coefficient 1 and zero offset. This is not yet a voltage-to-human-excitation or pressure conversion.

## Bilateral correspondence after the user's branch clarification

Suffix `{s}` means `_l` for LH and `_r` for RH. The table identifies names, not validated cross-species anatomical homology or torque equivalence.

| Rat MN channel | Actual source muscle name | Assigned Gait2392 names | Disposition |
| --- | --- | --- | --- |
| RF MN FLX | RectusFemoris | `rect_fem{s}` | Clear named muscle correspondence; retain its biarticular route |
| SO MN EXT | Soleus | `soleus{s}` | Clear named muscle correspondence |
| TA MN FLX | TibialisAnterior | `tib_ant{s}` | Clear named muscle correspondence |
| IP MN FLX | Illiopsoas, source spelling | `iliacus{s}`, `psoas{s}` | Explicit shared-drive group fan-out; record force capacity and feedback pooling |
| VA MN EXT | Vastii | `vas_med{s}`, `vas_int{s}`, `vas_lat{s}` | Explicit vasti-group fan-out, not an unnoticed all-vasti substitution |
| GA MN FLX | Gastrocnemius | `med_gas{s}`, `lat_gas{s}` | Explicit gastrocnemius-group fan-out |
| BFP MN EXT | BicepsFemorisPosterior | `bifemlh{s}` | User-directed biarticular branch correspondence |
| BFA MN EXT | BicepsFemorisAnterior | `bifemsh{s}` | User-directed uniarticular branch correspondence; source/target joint function differs |

The user clarified that `bifemsh` is the uniarticular branch and `bifemlh` is the biarticular branch. The working assignment therefore follows that distinction, while retaining this mechanical difference: rat BFA is used as a monoarticular hip extensor, whereas human `bifemsh` is a knee flexor. This is an explicit transfer experiment, not a claim that rat anterior biceps and human short head have equivalent anatomy, moment arms, or sensory timing. BFP/long-head correspondence retains a broadly hip-extensor/knee-flexor role, but still requires moment-arm validation. GA's source MN label FLX is not permission to command ankle dorsiflexors: muscle identity and its actual multi-joint action take priority over that label.

This working map innervates **12 target muscles per side, 24 of 92 total**, from eight source MN channels per side. The remaining 68 receive zero neural control, with passive mechanics retained. No decisions above authorize gluteals, sartorius, gracilis, adductors, peroneals, or lumbar muscles. RF/SO/TA alone formed the earlier six-actuator exact-name alternative, now superseded by the explicit eight-channel assignment. A sagittal harness is an explicit apparatus, not neural innervation of omitted muscles. Group feedback pooling must be recorded so three vasti do not accidentally triple an afferent signal merely because the source used one group.

## Prerequisites before a transfer claim

The imported 120-neuron public rat precursor's left RG oscillates, but its right RG latches in the actual native baseline. It is not yet a verified bilateral controller. First reconcile the source with its intended bilateral circuit, calibrate neural recruitment and sensor mappings, and record a successful baseline. Then establish the explicit target list, a voltage-to-excitation mapping for native Hill muscles or a pressure interface for BPAs, corresponding afferent scalings, sagittal/harness constraints, and the no-hyperextension gate. Preserve the source biarticular integration rather than replacing it with independent joint commands.

The inventory is reproducible with `audit_gait_muscle_names.py`; full target names, routes, adapters, and file hashes are in the ignored `runs\RAT_GAIT_NAME_AUDIT.json`. No additional model project copies were created. The proposed synthetic mass-change study was not executed and is not the requested transfer experiment.
