# Uniarticular W2L subnetwork rebuild

Session 1 subagent; easteregg2; 31 October 2026.

## Model changes

Reorganized BOTH requested projects into the historical uniarticular LH-style hierarchy:

- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_20261031\W2L_VirtualWalker.aproj`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer.aproj`

The template is `Biped_2xCPG_wSubs.aproj` from git **23b0e785**, the April 2021 rename carrying the June 2020 LH organization. Earlier commits were inspected and preserved in the local build archive. This precedes the recent partial biarticular additions. The historical LH groups are RG, HipZ, and knee/ankle synergy; synergy contains a four-cell PF page, knee MN/sensory circuitry, and ankle MN/sensory circuitry. The current newer model's contact and commissural circuitry receives its own pages.

This is an actual project reorganization, not only an external diagram. Existing functional objects are moved into pages and cross-page connections use genuine AnimatLab OffPage references. No biological neuron, functional synapse, actuator, or physical support was added. Current model parameters and connectivity were preserved; the old model supplied organization rather than replacing their numerical circuits or mechanics.

| Model | Actual neurons | Existing functional links | Diagram pages | OffPage references |
|---|---:|---:|---:|---:|
| Original | 80 | 222 | 14 | 32 |
| Newer | 96 | 292 | 18 | 70 |

The full diagram-object counts, 140 and 168, also include adapters, physical-body references, and stretch-receptor references; they are not neuron counts. All added OffPage objects reference existing objects and do not introduce simulated neurons.

The original has a shared left-labeled RG driving both sides; its RH organization therefore does not invent a separate right RG. The newer model has bilateral RGs and additional coordination/contact objects. The models retain their differences.

## Preservation and validation

Exact pre-edit project/standalone copies, hashes, and historical sources are preserved under:

`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_20261031\build\subnetwork_rebuild_20261031\`

The builder validates every nested drawing XML, object ID uniqueness, drawn object tags/text, node/link counts, and AddFlow's interleaved zero-based endpoint indices. Resolving each OffPage reference recovers the exact original **project** edge endpoints and synapse type. Original neuronal parameter payloads and physical-body payloads remain identical. Both `.asim` and `_Standalone.asim` remain **byte-identical** to their preserved copies because page organization changes only the project GUI layer. This preserves each supplied source separately; it does not prove the pre-existing project and runtime exports represented the same circuit. A subsequent audit found that the newer pair did not, as detailed below.

One pre-existing labeling error is documented rather than silently changed: node `a3828afc-f8fa-4497-a647-354322ee2bd6` is labeled `L Knee MN ext RE` but has five R-side neighbors and no L-side neighbors. It was placed in the RH knee page; the original label is retained consistently in runtime/project files and telemetry. The other similarly labeled node belongs to LH. Renaming it is a separate label-only cleanup; no parameter change was made here.

GUI-open validation could not run: CUA's separate-instance launch failed with `GetCursorPos failed: Access is denied (0x80070005)`. The existing user window is titled `AnimatLab Pro Biped_2xCPG_wSubs Project *` and has unsaved edits; it was not manipulated, closed, saved, or replaced. XML/drawing validation is not claimed as a substitute for a real GUI-open check.

## Native simulation checks

Executed the existing native `AnimatSimulator.exe` tests for both current standalone files. The default restricted launch first failed with native DLL relocation error 623; a normal elevated execution of the same installed simulator succeeded. No model workaround or physics change was used.

| Model | Exit | Samples | Recorded time | Root Y minimum | Root Y final | Root X displacement | Maximum numerical difference from prior recording |
|---|---:|---:|---|---:|---:|---:|---:|
| Original | 0 | 5,001 | 0–10 s | 0.960805 m | 0.963609 m | −1.006330 m | 0.0 |
| Newer | 0 | 5,001 | 0–10 s | 0.960777 m | 0.965053 m | −1.113340 m | 0.0 |

Both recordings are finite and bit-identical to the preserved numerical recordings. Original runtime was approximately 11.1 s including native startup and logging; this is an AnimatLab check, not a claimed accelerated MuJoCo benchmark.

World **Y is superior** in these models. Root `root_y` position telemetry was already present in the supplied projects' verification chart and has been preserved. A loss of body height is therefore directly observable in headless runs. Both models remain supported by their passive hip-wing/rail geometry; root height near 0.96 m does **not** prove free standing, free walking, or a completed training curriculum. Root X movement on the moving floor does not by itself establish forward locomotion. Contact counts/positions do not constitute measured harness or rail force.

This task did not modify joint limits or claim that old recordings satisfy the later no-hyperextension gate. New training must enforce and measure that requirement separately.

**Source-unit caveat from the parallel MuJoCo conversion audit:** the newer model's four ContactCount-to-neuron polynomial adapters use `C=20` in the native runtime representation, which the conversion audit identifies as 20 A/contact rather than the intended nA scale (native neuron currents are converted to nA with a factor of 1e9). This layout-only pass preserves that source numerically. The MuJoCo agent is creating a separately named, opt-in SI-corrected variant. Finite native execution and parity do not establish that the inherited feedback units are correct. Do not silently replace the supplied source with the corrected variant or combine their outcomes.

### Newer GUI/runtime discrepancy and separate editable companion

The newer project has **96 neurons and 204 synapses**, whereas its supplied standalone has **96 neurons and 208 synapses**. Ten native-only endpoints connect the two toe-contact neurons to their ipsilateral RG/PF/hip-knee MN extensor populations. Six project-only endpoints connect heel/toe to ipsilateral ankle MN extensors and heel to contralateral RG flexors. Several synaptic type choices differ. The project also sets left RG-extensor tonic current to 2 nA, versus 0 in the supplied runtime. Its contact polynomial coefficients are C=1, D=0, versus native C=20, D=−4e−10. These are inherited discrepancies, not changes made by the page reorganization.

The parent requested a separate companion using the **named editable project as graph authority**. The first tested companion uses an explicit **20 nA/contact, −0.4 nA offset engineering choice**, not a magnitude-preserving conversion of the GUI's C=1,D=0. Its existing filename remains:

`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_SIContact_20261031\`

See that folder's EXPERIMENT_RECORD.md, BUILD_SYNC.json, EXPORT_VERIFICATION.json, and NATIVE_SMOKE.json. Its 204-edge graph and used synaptic parameters match the editable companion; it is not the native-faithful 208-edge diagnostic that the MuJoCo agent ran earlier. The companion completes native simulation with bounded neural voltages but **falls**; no gait success is claimed. A separate GUI1nA companion preserves GUI coefficient magnitudes under a declared nA-unit assumption (C=1 nA/contact,D=0). Neither coefficient choice is claimed as Ben's intent or an identified optimum. Both retain the editable source's left/right RG-extensor tonic drive asymmetry of 2/0 nA.

## Review graphics and provenance

Each model folder now has `subnetwork_review/` containing:

- `LH_contact_sheet.png`, a review sheet from the actual saved page objects, locations, and links;
- individual `page_*.svg` exports for inspecting every saved page at readable zoom.

These are reconstructed reviews of the saved project, not screenshots from a successfully opened AnimatLab GUI. Amber dashed boxes are OffPage references; blue subsystem boxes are navigation objects, not neurons. Inspect the individual SVGs for readable labels and connections; the overview intentionally summarizes many pages.

`SUBNETWORK_REBUILD.json` records exact source/output project hashes and structural checks. `SUBNETWORK_VERIFICATION.json` records finite output, baseline parity, root-height results, and page-file paths. Existing BUILD/PRESERVATION records describe the earlier physical-build pass; this record supersedes their statement that drawing text remains unchanged.

## Commands on easteregg2

Replay this frozen subnetwork rebuild from the preserved pre-edit copies:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_20261031\tools\rebuild_subnetwork_style.py'
```

Native checks, run from `D:\GitHub\Bipedal_Robot`:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_20261031\tools\run_test.py'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_20261031\tools\run_test.py' --newer
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_20261031\tools\render_subnetwork_review.py'
```

`inspect_history_subnets.py` captures initial model files only if the corresponding archive file does not already exist. The subnetwork builder deliberately replays the frozen archived source; it must not be presented as a generic updater for later geometry edits. The older geometry/physical-walker builder recreates flat source pages if rerun; reapply the subnetwork builder only when its preserved source is the intended physical version.

All work here is separate from the BPA neural-coordination experiment and the Gait2392 training record. MuJoCo conversion/training is handled by the sibling agent using the unchanged standalone mechanics/circuit. No unsupported walking result, exact literature replication, or biarticular training result is claimed by this reorganization.
