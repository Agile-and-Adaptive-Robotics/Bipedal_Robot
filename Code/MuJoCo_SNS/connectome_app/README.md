# Connectome Editor desktop

Launch `Start_Connectome_Editor.cmd`. The desktop app embeds the complete original
connectome editor, adds native project files, and uses the schematic symbols in
the same editable canvas. It does not require the preview web server.

On EB475WS4, projects default to
`D:\Users\Ben Bolen\Connectome_Editor_Projects`. The original HTML editor and
its browser session are separate. The first launch creates `Paired_MN.sns.json`
as an editable copy of the included example. Open accepts existing editor JSON
and native `.sns.json` projects; compiled audit ledgers require explicit adapters.

## Editing and files

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

Run `python -m unittest discover -s tests -p "test_*.py"` and
`node tests/test_editor.cjs`. These verify atomic save/reopen, preservation of
source and parameter fields, file conflict protection, nested boundary endpoints,
copy ID remapping, and consumption of edited conductances by SNS-Toolbox.
