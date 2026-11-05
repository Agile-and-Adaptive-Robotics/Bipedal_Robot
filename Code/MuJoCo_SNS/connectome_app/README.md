# Connectome Editor desktop

Launch `Start_Connectome_Editor.cmd`. The desktop app embeds the complete original
connectome editor, adds native project files, and uses the schematic symbols in
the same editable canvas. It does not require the preview web server.

On EB475WS4, projects default to
`D:\Users\Ben Bolen\Connectome_Editor_Projects`. The original HTML editor and
its browser session are separate. The first launch creates `Paired_MN.sns.json`
as an editable copy of the included example. Open accepts existing editor JSON
and `.sns.json` projects, plus native AnimatLab `.aproj` / `.asim` models.
Compiled audit ledgers require explicit adapters.

## Editing and files

To open the prepared controllers, relaunch the desktop app and select **File →
Open controller**. The menu distinguishes reviewed Syn6, the saved s3k and
w2lvar standing conditions, and the editable native W2L108 copy. These are separate conditions;
the menu does not relabel a standing archive as a later Hunt result.

Their working copies are in
`D:\Users\Ben Bolen\Connectome_Editor_Projects\Controllers`. You can also use
**Open / Ctrl+O** and select a `_hierarchy.sns.json` file there. Double-click a
layer module, then Left or Right, to reach its cells. **F4** returns one level.
The companion `_exact.sns.json` files contain the same individual connections
in a flat view. The library retains zero-conductance definitions and source
records. A wire marked **×N** displays N individual synapses; selecting it lets
you choose an individual member and edit its parameters.

A template starts without a file location. Its first **Save** therefore opens
**Save As**. Choose a name in the project folder above. Subsequent **Ctrl+S**
updates that project file. Opening a controller from the new menu opens its
already saved working copy, so ordinary Save has a destination. Use **File →
Open project folder** whenever you want to find your files. Saving a project
writes the editor JSON for JSON projects. For a native W2L108 model, Save updates
its actual `.asim` and paired `.aproj` files. Python controller source write-back
is not yet supported for the other families.

### Native W2L108

**File → Open controller → W2L108 — editable native model** opens
`D:\Users\Ben Bolen\Connectome_Editor_Projects\Models\W2L108\W2L108.aproj`.
The model package includes `W2L108.asim`, its texture/chart assets, and
`W2L108.editor.json`, which stores the hierarchy and drawing layout. Keep these
files together. **Save As** defaults to a new native model name in the same
folder and copies its companion files. It then continues editing that copy.

This working copy comes from the verified native `.aproj` / `.asim` pair and
preserves that pair's saved values. The separately recorded GENTLE108 execution
graph remains in `Controllers`; its external tuning overlays are not applied
to this native file entry.

Native neuron settings, individual chemical synapse parameters, neural additions,
deletions, and reconnections are written to both native formats. Native adapter
polynomial coefficients are editable in their saved SI units. Changing one
synapse clones its shared AnimatLab type when necessary, so other connections
keep their parameters. The native model supplies the effective `SynAmp` value;
the ignored per-connection `G` field is not used as the import authority.

- **Save / Ctrl+S:** write the current project, including edits inside a module.
- **Save As / Ctrl+Shift+S:** choose a new file and continue editing that copy.
- **Save Subnetwork As:** save the current module as a separate project. Its
  external boundary connections are recorded under `interface.boundary_bindings`;
  they are not silently converted into standalone input currents.
- **Open / Ctrl+O:** open a project in a new tab.
- **Export SVG:** save the currently displayed drawing.
- The title and file path mark unsaved changes. Closing the window offers
  Save, Discard, or Cancel. A file changed by another program is not overwritten
  by ordinary Save; reopen it or use Save As.

The original palette, drawing, layer tree, selection, group movement, nested
modules, copy/paste, tabs, undo/redo, and keyboard commands remain. The **Help →
Commands** menu lists shortcuts. Click a neuron or synapse for its parameters.
Electrical values have explicit units; saved source records remain inspectable.
Select a wire to show its bend handles. Source and destination ports may start
at different positions along a box. Rounded bends and semicircular crossing
bridges distinguish routing from junctions. Excitatory triangles have their
broad edge against the receiving object; inhibition uses filled circles.

The selection rectangle includes individual wire control points as well as
neurons. Drag any selected control point or neuron to move that selection.
The blue square handles at a selected wire's ends slide along the current
neuron's perimeter. They preserve the source and destination neuron identities.
Delete removes selected bends; Undo restores their positions.

The first overview shows **Left RG** and **Right RG** as separate blocks, with
their afferent wires attached to their actual receiving neurons. Left/right PF
and MN blocks appear beneath them. Existing generated layer folders are upgraded
on open; the underlying cells and synapses stay the same. **Ctrl+Shift+G** unpacks
the selected block one level, retaining its external child bindings and moving
internal wire bends along with the newly exposed blocks.

Enter **Left MN / Right MN** to open a paired
motor view with the saved Ia, II, Ib, Renshaw, and reflex-interneuron wiring.
The two selectors choose the muscles/cells to display. Gray blocks show the
actual off-page neurons connected to those circuits; their parameters refer to
the same model objects. **Show off-page connections** opens this kind of view
for other layers too. The checkbox controls visibility, not connectivity.

Native W2L108 views include their actual physical maps, muscles, and stretch
receptors. Saved standing exports also show their recorded physical signal
interfaces. Sources which contain only neural/current-port definitions retain
that scope. No unrecorded physical encoder or reflex connection is added.
Grouping and wire layout in a motor view are saved as presentation data while
the underlying neuron and synapse edits update the canonical model. **F4** or
**Layers** returns to the layer hierarchy.

## First supported model

The paired motor-loop example has 18 graded neurons, 30 chemical synapses,
14 physical symbols, and 20 physical paths. It includes Ia, II, Ib, and Renshaw
connections, spindles, muscles, and adapter maps. It is extracted from the saved
historical w2lvar standing ledger identified in the project's `source` field;
it is not a designation of a later Hunt winner. The fixture's complete source
records and conversion notes are preserved.

**Build Neural Model / F5** passes the current edited graph and parameters to
SNS-Toolbox, compiles the neural kernel, and reports its census, inputs, outputs,
and clock. It executes zero integration steps. Muscle/adapter links are stored
but are not executed by this neural-only builder. Persistent-NaP and spiking
models remain saved with their fields, but this first builder explicitly stops
on unsupported models rather than substituting graded dynamics. Legacy editor
files without complete electrical parameters can be drawn and saved; they need
those parameters before Build. All-controller imports and Simulink generation
are later work, not completed features of this milestone.

## Runtime and verification

EB475WS4 uses the separate runtime
`D:\Users\Ben Bolen\ConnectomeEditorRuntime`, created from the existing `myo`
Python with `--system-site-packages`. Its Qt installation does not change the
scientific environment. Elsewhere, install `requirements.txt` in a suitable
environment and set `CONNECTOME_EDITOR_PYTHON` to its Python executable, or launch
`python app.py [project.sns.json]` directly.

Run `python -m unittest discover -s tests -p "test_*.py"`,
`node tests/test_editor.cjs`, and `node tests/test_interaction.cjs`.
`python tests/desktop_smoke.py` exercises the generated Qt editor using a
temporary native copy and the machine-local controller library. These verify save/reopen, preservation of
source and parameter fields, file conflict protection, nested boundary endpoints,
copy ID remapping, and consumption of edited conductances by SNS-Toolbox.
