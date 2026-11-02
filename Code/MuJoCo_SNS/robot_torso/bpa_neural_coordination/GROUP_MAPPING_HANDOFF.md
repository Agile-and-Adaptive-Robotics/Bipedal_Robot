# w2lvar standing mapping — easteregg2

The source expansion gives each MN one homonymous II excitation pathway but **1–9 antagonist-muscle II inhibitory pathways**, each initially at full weight. Matched motor-current and joint-torque audits explain why the unmodified standing mapping suppresses useful recruitment. Both the original all-gain runs and frozen standing fits are tabulated by homologous muscle; these are source-model calculations, not biological recordings.

A distinct opt-in variant divides IIIN-to-MN conductance by the actual target-group fan-in. **All IIIN connections remain nonzero.** Production builders and their original weights are unchanged. The global standing configuration retains PF→MN=0, aggregate IBEXC→MN=0, and tonic-bias scale=5/6. This is a standing-mode experiment, not an unchanged literature controller, a new walking result, or the biw2lvar-Hunt-PT curriculum controller.

Results: a fresh numerical repeat is bit-exact; the 20 s study (including 0.6 s pinned neural initialization) passes guided standing with 3.05° maximum tilt, 42.54 mm sway, and 97.70% mean body-weight ground loading. Mean upward apparatus force is 2.18% BW and downward force is 0.0095% BW. Existing lateral and rotational guidance remains. Knees stay below −1.75°, and visual foot penetration reaches 4.17/3.68 mm. Eight seeded ±1% mapping-gain/tonic-bias checks pass; this is not an exhaustive eleven-dimensional corner study.

The earlier two-group Rustuna calibration failed sway and is preserved as a negative result. Its `BEST.json` is **not** this passing normalized mapping. No pressure-actuator transfer is claimed for this candidate.

## Evidence

- Checkpoint: `D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\group_diagnosis\NORMALIZED_STANDING_CHECKPOINT.json`
- Detailed record: `D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\group_diagnosis\EXPERIMENT_RECORD.md`
- Per-muscle currents: `D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\group_diagnosis\HOMOLOGOUS_MOTOR_CURRENTS.csv`
- Current alignment: `D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\group_diagnosis\PER_MOTOR_CURRENT_AUDIT.json`
- Common torque demand/recruitment: `D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\group_diagnosis\COMMON_JOINT_TORQUE_REFERENCE.json`

The legacy runner NPZ `act` field stores selected **excitation commands**, not native activation. New opt-in `.motor_state.npz` files contain all 92 actual native activations, excitations, and engine force-cache values. Native activation replay agrees exactly with these measured states. No old files or manuscript captions were silently rewritten.

The role audit flags per_tert_r/l as opposite to the source ankle_pf classification at this pose; gemellus also opposes its hip-extensor reference with a much smaller arm. These remain flagged, not globally reclassified. Full-RoM functional classification needs a separate investigation.

## Repeat the exact mapping

Use a new tag. The strict loader consumes every saved group/global gain and rejects changed controller or dependency sources; a changed source needs explicit revalidation.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\run_saved_group_candidate.py' --checkpoint 'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\group_diagnosis\NORMALIZED_STANDING_CHECKPOINT.json' --out 'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\group_diagnosis\manual_repeats' --tag 'NEW_UNIQUE_TAG' --duration '20'
```

Add `--validate-only` to check saved parameters/source identities without starting physics. This loader is needed because the simpler three-parameter loader does not carry the ten mapping gains.
