# Same-run verification capture and few-cycle figures

Machine **easteregg2**; Python `D:\Anaconda\envs\myo\python.exe`. Raw outputs:
`D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures`.

This pipeline records actual neural/mechanical signals without changing a controller or advancing physics. New body trials remain coordinated by the campaign owner. It does not launch calibration, training, or a body trial itself.

## Configuration for the next coordinated verification trial

Add to the actual `hunt_training\physical_eval.py` configuration:

```json
{
  "verification_capture": "D:\\Users\\Ben Bolen\\NeuroCampaign_20261102_9h\\figures\\CASE_NAME",
  "verification_max_samples": 20000,
  "verification_voltage_label": "graded model membrane depolarization (mV; rest-relative)"
}
```

Existing 8 s stage-3 standing and 16 s stage-4 stepping protocols are retained. `runner --time` is not used to shorten these runs because it adds schedule duration rather than specifying the total horizon. A default 20,000-row capture covers these protocols and raises if exceeded, rather than silently truncating.

Default selected cells per side: RG E/F; RF and HAM PF outputs and their two inhibitory relays; iliacus/bifemlh, vastus intermedius/biceps femoris short head, and soleus/tibialis anterior MNs. Exact `network.idx` identities are saved. Missing requested names are explicitly reported; `verification_cells` can supply an exact list. Actual activation states and all controls are saved separately, not labelled as voltages. A caller must supply an explicit different voltage label for a spiking-rate readout.

No pulse is inferred from neuronal activity. Actual stimulus logging is supported by `AlignedCapture.sample(stimulus_nA=...)` but remainsNaN when no actual stimulus vector is supplied.

## Capture point and frames

`runner.py` exposes `VERIFICATION_POST_PHYSICS=None` by default. When enabled through `physical_eval.py`, the observer runs once after the accepted outer physical interval, including adaptive substepping. It never observes rejected adaptive trial steps. The previous observer is restored in `finally`. Held-pose warmup intervals are counted and excluded from the physical trace.

The archive retains logical interval-end/state time, neural-output time, neural-input/interval-start time, and native `data.time`. This matters because preparation can advance the neural schedule while holding mechanics. Contact force is the result of the last solver interval, labelled by interval start/end; root body positions/contact frames retain the solver evaluation phase. There is **no extra `mj_forward`, `mj_step`, or neural step** in capture.

Legacy human world coordinates `(x anterior,y left,z superior)` are explicitly mapped to anatomical `(x anterior,y superior,z right)` by `[x,z,−y]`. The underlying model and saved coordinates are not rotated. New y-up models have an identity anatomical adapter. Recorded joint signs remain source signs; displayed Gait knees use negative native knee angle as positive flexion, with any positive native angle reported as hyperextension.

Foot contact is counted once, using actual contact-frame force transformed to world and then anatomical axes. Hind/forefoot regions use contact position projected from calcaneus toward toes; their superior component is saved. They are virtual regions, not calibrated anatomical heel/toe pads. The old runner's doubled heel/toe proxy is not reused.

The pelvis restraint wrench projects the solver's actual passive generalized root forces onto root axes. It is an apparatus wrench only insofar as those coordinates' passive forces are the known spring/damper restraints; the metadata preserves that qualification. Explicit weight support is recorded separately as the actual configured force inN. **BW normalization uses pelvis subtree mass**, excluding separate fixed ground/apparatus bodies; historical metrics retain their old denominators.

## Outputs from the same trial

- `aligned_trace.npz`: selected actual neural voltage/readout, exact identities, fullqpos/qvel, actualactivation/ctrl, independent regional/full-foot forces, root restraint wrench, explicit support, all clocks, and stimulus telemetry when supplied.
- `model.mjb`: compiled runtime model, including numerical patches; no reconstruction from a different source model during replay.
- `ALIGNED_TRACE_METADATA.json`: source/configuration provenance, compiled-model SHA256/MuJoCo version, frame/timing/sensor definitions, human mass/BW, and missing channels.

## Plot and replay command (no physics integration)

```powershell
& 'D:\Anaconda\envs\myo\python.exe' `
  'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\verification_figures.py' `
  --source 'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures\CASE_NAME' `
  --output 'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures\CASE_NAME\plots' `
  --label 'Exact case identity from the verification result' --cycles 3 --render --fps 20
```

Use `--start 5 --end 9` for an explicit common time window. Without a window, the script selects approximately three recorded RG cycles when peaks can be identified; otherwise it labels the estimate unavailable and shows the last3s. It does not force phase alignment or warp a trace. Keep short labels; identity is retained in the report.

Output plates: `ALIGNED_NEURAL.png/.pdf` (bilateralRG, four-cellRF/HAMPF, representative hip/knee/ankleMNs), `ALIGNED_MECHANICS.png/.pdf` (hip/knee/ankle angles, regional/full-foot loading, restraint force and support), and combined `ALIGNED_VERIFICATION.pdf`. `FIGURE_REPORT.json` records window, missing channels, actual joint limits/observed hyperextension, sourcehash, and timing interpretation. No curriculum acceptance is inferred from a plot.

Replay produces `RECORDED_REPLAY.gif`, `RECORDED_STILL.png`, and `REPLAY_PROVENANCE.json`. It loads the exact hashedMJB, restores saved states, and only calls `mj_forward` for rendering; **zero `mj_step` calls**. It refuses a missing/mismatched compiled model or incompatible MuJoCo version. No source-derived replacement model, corrected pose, or additional support is introduced. State clocks/frame transforms are retained.

Use `--still-only --start 5.999 --end 6.001` to render only the nearest saved state to 6 s, without plots or GIF generation. The exact selected time remains in `REPLAY_PROVENANCE.json`.

## Completed pipeline check, distinguished from new results

The plotting adapter was run on the **existing historical** `fourcell_hunt_20261029_v1\ground_isolation\full_feedback` recording, with6–9.6s selected. Its new plots/report are under `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures\pipeline_check_historical`. Both plates were inspected. The old file lacks iliacus/soleus/tibialis-anterior MNs and regional contact forces; those gaps remain marked rather than fabricated. Its logical pre/post-step timestamps are explicitly flagged as unsuitable for substep latency inference. This check is not a new integrated-controller result.

All capture/plot/runner files pass Python syntax checks. MuJoCo 2.3.7 API introspection verified MJB saving and subtree-mass availability; its renderer lacks `.close()`, so the replay uses version-compatible disposal.

## Actual source-baseline capture verification

The campaign owner's completed stage-2 source-baseline trial now validates the post-physics capture and exact-model replay end-to-end. Outputs are `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures\source_baseline_s2\plots`; source capture review is `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures\source_baseline_s2\CAPTURE_REVIEW.json`.

All 24 requested cells were present and finite across 7,700 captured physical intervals. Neural-output and state timestamps agree exactly; neural input precedes them by 2 ms. Native physics time differs from logical time by the correctly excluded 0.6 s/300 held-pose preparation intervals. Regional superior forces sum exactly to the once-counted whole-foot superior force. The archived MJB hash is exact. These are telemetry/replay checks, separate from the campaign owner's stage acceptance.

The selected figure/replay window is 7.645–13 s, approximately three actual cycles. It excludes the scheduled return to standing, which the automatic last-cycle choice would otherwise include. Use `--end 13 --cycles 3` for this protocol. The source-baseline label is retained; these files do not illustrate the newly composed FSA candidate. Ground contact is zero in this air trial. The disabled floor visual is hidden only in the rendered scene so it does not occlude the feet; the caption and replay provenance state this. Saved model, poses, and contacts are unchanged.

Neural and mechanics plates use a native 6.5-inch figure width with 10-point labels. Both plates and the still were inspected. No separate body trial was launched by the figure worker.

## Optional actual sparse-current telemetry

For new trials only, set **`verification_sparse_currents: true`** together with the selected **`backend: "sparse"`** and an enabled `verification_capture` folder. This option is off by default. Existing captures and figures remain valid when the optional fields are absent. The observer passes the configured backend into `sample`; it does not change the input mask, neural parameters, physics, or stepping.

For each selected neuron, the optional archive contains actual applied external current (`SparsePairwiseForward.i_app`), chemical current from the last actual conductance sums and `V_last`, leak, bias, ionic current, preceding voltage, actual membrane conductance/capacitance, resting voltage, update factor, and update residual. The fields are prefixed `neural_`; `neural_current_status` and `neural_current_available` identify validity per sample. The actual update is

`V_new = V_last + time_factor × (I_leak + I_bias + I_chemical + I_applied + I_ionic)`.

Applied current includes direct mapped POST inputs when present. It is **not** labelled as a stimulus: the independently logged `stimulus_nA` retains its own missing/present meaning. Actual MuJoCo `data.act` remains separate from `data.ctrl`; the legacy main NPZ's potentially command-valued `act` is not used to infer activation or current.

Availability requires the actual sparse forward binding, the actual one-hot applied-input buffer, a `V_last` alias established by a real sparse step, finite fields, and update closure within 1e−11 mV. Native-forward overrides can retain stale sparse buffers; native, unsupported, pre-first-step, or inconsistent samples remain **NaN/missing**, never synthetic zero. Ionic current is an algebraic readout of the actual already-updated gates; no gates or neurons are advanced. Selected passive MNs legitimately have zero ionic current.

The actual `SwingDiagnostics` recorder is also supported. Its exact nested function code and captured network, compiled instance, and native sparse delegate are verified; an arbitrary wrapper or a diagnostic wrapper around native SNS remains unavailable. Buffer-alias and numerical closure checks still apply after this recognition. The two-step electrical fixture at `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures\current_telemetry\diagnostic_wrapper_proof\STATIC_CURRENT_CLOSURE.json` verifies both direct and actual diagnostic-wrapped sparse calls. All selected update residuals are zero, both observations preserve complete neural/buffer hashes, and the wrapped-native negative control is rejected. Observation adds zero steps. Earlier captures that rejected the diagnostic wrapper remain missing-current; they are not retrospectively filled.

Pre-change capture source and its hash are preserved in `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures\current_telemetry\prechange`. The actual 1,033-cell Hunt static proof is `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures\current_telemetry\static_proof\STATIC_CURRENT_CLOSURE.json`: all 24 selected updates close exactly, including 0.7/0.3 nA direct POST currents and nonzero chemical current. Complete wrapper/native array hashes agree before/after observation. Default-off, native/stale, non-one-hot, and intentionally inconsistent-current controls pass. This is one declared electrical step, with the unbound balance updater excluded; no body integration or training occurred.
