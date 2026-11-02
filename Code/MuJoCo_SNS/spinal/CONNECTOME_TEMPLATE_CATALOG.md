# Curated editor template menu

Machine: easteregg2. Scope: visibility/labels only; no neural graph, trained parameter, simulator, or browser-storage change.

Open `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_block_editor.html` (or the existing localhost editor). The default menu now contains **Blank** plus five explicitly labelled references: Ben's circuit rules, six-component NMF data, and Rybak/Shevtsova/Shinohara rule references. These are not presented as executed current controllers.

**Legacy templates** is off by default. Enabling it reveals all35 retained menu entries, including v9/v10 winners, curriculum snapshots, pre-schematic flat duplicates, pruning experiments, source-model mines, older spiking mirrors, and earlier Hunt exports. The graph data remain embedded/in the original sidecar. Disabling Legacy hides those choices again without loading/deleting a graph or changing a tab.

Current W2L108/704, current Gaitspiking, rat134, CAD262, and saved current production-family graphs are **not automatically embedded**. Import their actual exported graph and source/parameter manifest. An old September family schematic must not be relabelled as a newer executed winner. The catalog records this availability limitation explicitly.

Authoritative visibility policy: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\editor_template_catalog.py`.
Generated small metadata: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_template_catalog.json`, also embedded inside the HTML for file-mode use. Full original graph library: `connectome_templates.json`.

To refresh metadata **without building any neural model**:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\make_editor_templates.py' --catalog-only
```

The ordinary full generator also reapplies the same curated policy after generating its historical library. It was not run in this pass. A future graph import needs a new explicit metadata entry; unknown templates default to Legacy.

Verification: both embedded JavaScript blocks pass Node syntax checks. An isolated fake-DOM test of the actual menu function yields6default choices and35withLegacy, then returns exactly to6. The tested function contains no localStorage access, persist, or loadIntoCanvas call. No user browser page was reloaded or manipulated. Existing localStorage/user drawings remain untouched. Sidecar graph SHA256 remains `208f11b514fc2fb0f67cc0d2218cf9c883d9c4d3d8fb5a9d6c583a6b5ed3e90a`.

**Training gate:** design/scaling and individual-pathway review for all active families precede any training calculations. This menu cleanup performs no simulation, calibration, or optimization.
