# Latest author clarification

The coupled SolidWorks-derived robot assembly is confirmed by Ben and is now explicitly included in P01. Its saved artifact is temporarily unavailable on the laptop. The table gives representative examples for the toolset, not exhaustive validation of every tool. Rybak remains listed as Python. See AUTHOR_CONFIRMATION.md. Earlier evidence distinctions below remain useful for attributing numerical results.

---

# Tools created and representative demonstrations — P01

## Editorial decision

Use **functional verification and proof-of-concept demonstrations**. “Functional verification” means checking the implemented function; “proof-of-concept demonstration” means using it in a representative model or experiment. Avoid “fundamental proof,” which implies a stronger mathematical or scientific result. P01 is a contribution inventory, not a list of unfinished work. Ben rejected the former incomplete-integration framing; its text is archived here, not recommended.

The paragraph and capability table are a paired targeted change. Replace the live paragraph beginning “I developed several software tools” through its capability-table sentence with `P01_replacement.tex`. Replace only the table environment bearing `tab:sns_capability` with `P01_companion_table.tex`. Preserve that label and existing review threads. If Ben’s current paragraph differs, reconcile its wording before applying. Do not upload an entire chapter. No Overleaf edits were made.

## Evidence-supported inventory

| Contribution | What can be claimed | Small representative demonstration and record |
|---|---|---|
| Custom Li–Deng/W2L model and controller implementations | Connor–Ben’s original shared-RG AnimatLab combination, Ben’s later per-side RG/feedback extension, and MuJoCo/SNS-Toolbox implementation; matched neural components in Simulink | Contact-driven core: 1.000 s periods, r≥0.99996, matching contact-removal latch. `spinal/campaigns/20260930/w2l_trace_compare.json`, `sns_simscape.md`. The Simulink result names the neural core, not a complete mechanical twin. |
| SNS-Simscape library and custom GUI | Editable sensory, neural, synaptic, activation, and biological/BPA force blocks; palette, parameter editing, and network import/export | `gui/SNSGui.m`, README GUI smoke record; actuator-reference checks; beer reflex ON/OFF 2.00°/7.85°. These are direct tool-function demonstrations. |
| Gait2392/SNS controller/model implementations | Repaired OpenSim-to-MuJoCo human reference, alternative s3k, w2lvar, and syn6 neural organizations, and custom routing/actuation models | Single-actuator sign/path tests; selected supported human stepping and recorded reference comparisons. Build census and saved parameters distinguish each neural variant. Do not imply identical outcomes for all controllers. |
| BPA-force integration | Custom Festo force computation coupled to mechanical routes; this is not functionality created by MuJoCo or SNS-Toolbox | `spinal/campaigns/20260930/bpa_run/bpa_gate.log`: all 12 W2L muscles use BPA force, bilateral air excursions and supported standing, peak 372 kPa within 620 kPa ceiling. This is distinct from the 92-entry Gait2392 human model and robotbody route. |
| Bidirectional SNS-Simscape/MuJoCo integration | Pairing the new neural library with the existing MathWorks bridge | `mujoco_bridge/BRIDGE_REPORT.md`, tests a/b/c and sensor readback; matched integration data. Credit the bridge authors. Bridge coupling does not convert a Simscape mechanical model into MuJoCo. |
| CAD/musculoskeletal import helpers | Imported mechanisms and model-building/import workflow in MATLAB | Saved SolidWorks XML/URDF imports, inventory logs, and OpenSim-derived Simscape model. Ben confirms that the imported robot assembly is coupled to a neural controller and BPAs using his force model. Its saved model is on the temporarily unavailable laptop. Include this author-confirmed contribution without inventing a file location or numerical result. |
| Automated curriculum, optimization, and experimental controls | Five-stage variant curriculum, bounded Optuna searches, support scheduling, neural stimuli, pathway interventions, timed pushes, body-parameter changes, and output recording | `_curriculum_w2lvar.py`, `_curriculum_syn6.py`, saved stage trials and fresh replays; corrected W2L fixed-controller interventions. Optuna is reused; the staged objectives, drivers, coupling, and protocols are the project’s work. |
| Literature-derived circuit implementations | Reduced Rybak and Shevtsova neural models, Deng-derived circuits, Li-derived adaptation, and Di Russo-derived controller coupled to Gait2392 | Rybak `verify_rybak2006.log`: 13 checks including resetting/non-resetting deletions. Shevtsova has a reduced SNS implementation and recorded neural checks. Di Russo has stage-1 fit and completed stage-2 controller evaluations. List the actual demonstrated function, not successful reproduction of every paper’s locomotion. |
| Selected spiking neural counterparts | Spiking library blocks, selected demos, production-network mirror, and archived AnimatLab conversions | `sns_build_spiking.m`, `sns_units_test_spiking.m`, `build_network_spiking.py`, Appendix B’s implementation/calibration and paired comparisons. “Spiking mirrors of everything above” is broader than the identified records. |

## The five-stage sequence

The variant drivers explicitly implement (1) deafferented air stepping, (2) afferented air stepping, (3) standing balance, (4) walking-pattern tuning without ground contact, and (5) walking with ground contact. The stock `_curriculum.py` has a different historical fifth stage (PF-layer contact variant). Name the driver/sequence used; do not merge their stage labels. “Implemented curriculum” is the tool contribution and does not require claiming every candidate completed every stage.

## Wording of the gap

The project lacked a reusable way to assemble alternative neural controllers with a many-actuator human reference and BPA-actuated models, specify support and sensory conditions, and automate controlled evaluations. The contribution is the implemented models, block-programming environment, custom force integration, staged tuning, and experimental controls. Do not claim Python is inherently difficult or no comparable platform exists. The demonstrated operations explain the gap and its resolution without a subjective language comparison.

Use “many-actuator musculoskeletal model” or “muscle-redundant model.” The latter denotes more muscle actuators than joint degrees of freedom; “over-actuated” may suggest independently controllable actuator directions and needs a more specific control definition.

## Claim attribution and precise placements

- Rybak’s verified implementation located here is Python, not a demonstrated Simulink Rybak model. Put it under the literature-model toolset until the claimed Simulink filename is supplied.
- Shevtsova’s standalone/SNS lumbar circuit is not evidence of a completed Shevtsova–Gait2392 coupling. Shinohara rules/specification are not a runnable Shinohara–Gait2392 model. Di Russo–Gait2392 coupling exists and can be listed as an implemented model.
- A Simulink controller communicating with a MuJoCo model through a bridge does not establish a CAD→Simscape→MuJoCo mechanical export route. List CAD import helpers and neural/physics integration separately unless the exporter/result is identified.
- The new ablation report is a contact+force ablation with Ia/II retained, not full deafferentation. It also retains stale zero-sensitivity interpretations and calls a W2L model BPA without an actuator manifest; do not import those claims into P01 without matching the actual model. Its saved metrics can be evaluated separately.

These are claim-boundary and attribution notes for the editors. The proposed P01 prose and table lead with delivered tools and their demonstrations, not outstanding tasks. No new computation was launched.
