# Hunt timing: source saturation and a separate canonical candidate

Machine: **easteregg2**. Existing helper extended: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\hunt_timing_profile.py`.

The old four-cell condition is `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\fourcell_dynamic_20261029_v1\curriculum\stage2.json`. It reconstructs1027neurons, not the HuntPT condition's1033. Its omitted timing_scale defaults to1. A raw feature1.5 therefore gives7.5nA and a passive relay voltage7.5mV, above its5mV outgoing synaptic saturation. The previous **unsaturated** equivalent mode correctly rejects that domain. Its fixed-window/gain-scaling contract remains unchanged and passes a new HuntPT regression using the saved0.01-scale profile.

## Exact saturated-source representation

The explicit `saturated_source_equivalent` mode changes the relay representation while retaining every original Gmax. With feature clipped0..1.5, normalized relay input is `5*feature/1.5`, and `r=a_new/a_old=2/3`. Relay states scale by r, and each existing outgoing threshold/saturation window scales by the same r. Thus:

`clip((r*V-r*elo)/(r*ehi-r*elo),0,1)` equals the original conductance fraction, including saturation for source features≥1. The new encoding spans0..5mV, while the transformed outgoing upper threshold is3.333333333mV. The actual dense e_lo/e_hi and sparse elo/erange caches are updated and verified.

The primary source-preserving variant retains20ms relay C/G. Its old-versus-new native-current grid, including features0/0.25/0.5/0.75/1/1.25/1.5 and RG Vm−2.5/0/2.5/5/8mV, agrees within2.23e-14nA. A120-step source-timing waveform agrees to1.78e-15mV after rescaling and5.33e-15uS in delivered conductance. The optional fast variant retains this static function but uses the separately declared3.896ms relay budget; its transient behavior intentionally differs.

The old four-cell folder has no timing_weights.json and loads separate step_B/right and step_B_left gains. Exact equivalence therefore **retains historical asymmetry**: full-feature reference currents at RG2.5mV are64.2183nA left versus8.54152nA right for RG-E, and1.67087 versus0.910972nA for RG-F. These are electrical endpoints, not measured standing currents. The equivalent representation does not repair those gains or certify the old controller.

## Explicitly new canonical-right role budget

`canonical_right_role_budget` is a different candidate, not an old-function correction. It keeps each RIGHT E/F pool's relative source conductance pattern, derives a0.5nA full-feature incremental budget per RG pool at2.5mV, and copies the resulting exact weights and windows to the homologous left side. It never averages left and right patterns. Six relays per side span0..5mV; existing source eligibility, targets, excitatory reversal, and5mV interface must match. Both selected candidates use the explicit fast relay option.

Actual construction, source live C.step consumption, command-policy persistence, prepared-state mapping, sparse/dense windows, per-feature one-hot currents, and maximum coactive input currents pass on both the old four-cell and HuntPT circuits. Canonical one-hot L/R differences are at most1.12e-14nA; full-pool current errors are at most1.16e-14nA. Future fitting must remain one-half plus mirror. No body, Hunt fit, or curriculum has run.

## Full profile paths

- Source-preserving four-cell: `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_timing\fourcell_saturated_source_v4\TIMING_PROFILE.json`
- Independent fast timing option: `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_timing\fourcell_saturated_fast_v4\TIMING_PROFILE.json`
- New canonical four-cell candidate: `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_timing\fourcell_canonical_right_v4\TIMING_PROFILE.json`
- New canonical HuntPT candidate: `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_timing\HuntPT_canonical_right_v4\TIMING_PROFILE.json`
- Preserved0.01-scale source-equivalence regression: `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_timing\HuntPT_unsaturated_regression_v4\RESULTS.json`

Each of the four new variants has a fresh_process repeat. All result, applied-ledger, and profile JSON triplets are byte-identical across processes. `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_timing\V4_PROFILE_COMPARISON.json` records exact paths/hashes, and the shared exact helper snapshot is SOURCE_hunt_timing_profile_v4.py in that folder. Earlier histories remain unchanged.

The first fresh saturated check exposed an audit issue: the sparse backend allocates last_voltage with np.empty, then overwrites it before every read. Before it aliases authoritative V_last, it is uninitialized scratch, not a prepared neuronal state. The helper now initializes selected scratch entries from causal V rather than range-checking uninitialized memory. Authoritative V/V_last domain guards still reject invalid states before mutation; no voltage projection or tolerance relaxation was introduced.

These are electrical design/component checks. Whole-body verification, Hunt training, and curriculum flags remain false. Existing source models and the parent's composition verifier were not edited.
