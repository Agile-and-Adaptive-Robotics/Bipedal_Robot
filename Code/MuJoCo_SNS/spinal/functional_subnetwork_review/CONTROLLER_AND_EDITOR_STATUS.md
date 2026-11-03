# Schematic coverage and connectome-editor update

All-controller schematic coverage is **incomplete**. The last combined PDF is a set of selected saved-condition views and representative MN examples, not a finished current drawing set for every family.

| Family | Current drawing coverage |
| --- | --- |
| s3k | Saved standing RG/four-PF architecture sheet; no dedicated expanded MN/map sheet. |
| w2lvar | Saved standing RG/PF feedback plus representative shared-PF and biarticular MN/physical-map sheets. Its PF→MN routes have effective zero conductance. This is not the selected later Hunt controller. |
| Syn6 | Earlier PDF uses the preserved historical standing ledger. The current reviewed composition was retrieved and hash-verified; the editor preview has a focused current heel→InE→RG-F view. The full new 22-muscle feedback/RG/PF/MN drawing remains pending. |
| biw2lvar-Hunt-PT | Selected whole-body identity/caller/profile paths obtained from the newer results handoff; no complete selected-controller overview in the new drawing style yet. |
| RF/HAM four-cell | Earlier PDF uses an old Session1 projection. Current executed ledger retrieved; the editor preview has the exact selected left module, 12 cells/12 synapses, with the full bilateral 24-cell/24-edge source ledger retained. |
| W2L108 graded | Native-checkpoint RG/PF overview and paired hip MN circuit. Full native physical adapters are not yet redrawn in the new style. |
| W2L704 hybrid | Saved logical/population-interface sheet; expanded spiking/readout and MN physical-map sheets pending. |
| Rat134 | Curated receptor-geometry central feedback sheet with separate hip, knee, and ankle PF modules; expanded MN/maps pending. |
| CAD262 | Earlier pressure-pair detail exists; matching new-style controller/subnetwork set pending. |
| Gait TimingZeroSFA | Named saved central-population interface sheet; expanded MN/maps pending. Do not replace it with settings from an earlier Gait variant. |
| SNS-Simscape reflex demonstrations | Source evidence located; matching new-style graded/spiking example sheets pending. |
| Reduced Deng14 | Earlier fixed-h 14-cell/18-edge drawing exists; matching polished new-style sheet pending. |

The newer handoff at the top of FUNCTIONAL_SUBNETWORK_OTHER_MACHINE_PROMPT.md supersedes the earlier Syn6 and RF/HAM selections. The reviewed Syn6 graph SHA256 is `697678b3281d7fb3b2a9322c2bc77b0671e7d05beb932fcb4364b8aa5ddaddea`. Its heel→InF branch is absent. Earlier historical drawings must not be relabelled as that controller. The selected Hunt result remains a separately named model/controller.

## Working editor preview

Open `D:\Users\Ben Bolen\NeuroGraphics_20261102\editor_preview\Connectome_Schematic_Preview.html`. Keep schematic_editor.js beside it; the example library is embedded for local-file use. The agent's local preview server also serves this folder at `http://127.0.0.1:8766/Connectome_Schematic_Preview.html` while that process is running.

This is an isolated styling/layout prototype, not a replacement for connectome_block_editor.html. It does not access browser storage or modify the existing editor, controllers, or user drawings. It edits glyph labels/sizes, node positions, wire bends, and source/destination port positions. Source graph records are inspected read-only. It supports zoom/pan, undo/redo, schematic-view JSON import/export, and SVG export. Unexported edits are session-only; switch views after exporting work to keep.

The four editable examples are the two historical w2lvar full MN/map sheets, the current selected left RF/HAM module, and the current reviewed Syn6 heel/RG pathway. The source graph is retained separately from the selected drawing view, including every original field. A bundle's member records remain individual synapses. Full source identity and SHA256 are carried with the document. Source-record round trips and source immutability under layout edits are checked for all four examples. These are data/graphics checks, not neural executions.

Implemented visual conventions: yellow PF cells, blue MN/reflex cells, circular salmon AD/maps, vector spindle/muscle/tendon symbols, physical arrows distinct from synapses, destination-mounted open triangles with their base on the object, inhibitory circles, distributed explicit ports, rounded bends, semicircular wire bridges, and visible dashed disabled routes. Node movement keeps endpoints attached through normalized port coordinates.

## Integrating this into the existing editor

1. **Separate model data from presentation.** Keep stable source IDs, full node/synapse/input/actuation records, and source/configuration identity in the model. Store positions, display labels, shapes, module frames, ports, bends, visibility, and annotations in named views. Changing a diagram must not silently change connectivity or parameters.
2. **Share one renderer between editor and exports.** Port the preview's shape-aware attachment, terminal symbols, rounded paths, and crossing bridges. Add native physical types and distinguish sensor/current/activation signals from chemical synapses. Preserve 45-degree routes when editing rather than force everything onto a right-angle grid.
3. **Use view bundles and reusable subnetworks.** Provide paired MN, shared-PF muscles, biarticular convergence, joint PF, four-cell RF/HAM, and sensory-adapter views. A display bundle carries explicit member IDs and counts; it never replaces several muscles with one neural pool. Inserting a motif must bind its ports to actual cells and preserve variant-specific signs/weights.
4. **Repair import/export before adding current controllers to menus.** The existing editor's specOf/materialize/reload functions copy a restricted set of fields, regenerate UI IDs, and omit arbitrary source metadata. Extend the schema so all Erev, gmax, threshold/saturation, units, clocks, receptor/map parameters, enabled status, and provenance survive. Do not reinterpret a compiled-graph file as an editor document without an adapter and parity checks.
5. **Preserve actual boundary endpoints through subsystem operations.** Existing pack/unpack retargets boundary edges to a SUB glyph and later removes edges incident on that glyph. Record child endpoint/port identity separately, then test pack→edit→unpack, copy/paste, undo, save/reload, and export against the original graph. Display-only grouping should be available without rewriting any graph endpoint.
6. **Keep the curated reference/current distinction.** Current controller imports require their exact executed source and manifest. The new styles are reusable references, not an excuse to make old menu snapshots look like newer winners. Add controllers only when their source adapter preserves all fields and their layout has been reviewed.

The preview demonstrates the renderer and the data/view separation. The desktop milestone below implements file editing and subsystem preservation separately from that preview. Source-specific compiled-schema adapters and complete controller-template coverage remain outstanding. No blanket migration or existing browser-storage rewrite was performed.

## Desktop editor milestone

`Code\MuJoCo_SNS\connectome_app\Start_Connectome_Editor.cmd` launches a Python/Qt
desktop host containing the complete original editor and the schematic renderer.
It adds native Open, Save, Save As, Save Subnetwork As, explicit electrical
parameter editing, and an SNS-Toolbox neural builder. Existing browser storage
and the original editor file are separate. Full node/edge/source fields and
stable identities survive file operations, nested pack/unpack, and subtree
copy/paste. Boundary wires retain child endpoint paths instead of being discarded.

The first editable fixture is the historical w2lvar paired MN extraction:
18 graded neurons, 30 chemical synapses, 14 physical symbols, and 20 physical
paths. Physical maps remain saved external bindings. Build compiles the edited
graded neural graph with zero integration steps; unsupported gated/spiking
models stop with an explicit message. This is not all-controller import coverage
or Simulink generation. See `Code\MuJoCo_SNS\connectome_app\README.md` for launch,
file locations, supported scope, and verification commands.
