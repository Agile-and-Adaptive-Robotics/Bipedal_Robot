# Desktop milestone verification

## Controller library update

The controller library adds native-file menu entries and full graph hierarchy
copies. Static checks preserve all cells, individual synaptic parameters, source
records, zero-conductance entries, and nested endpoint paths: reviewed Syn6
794/3,200; saved standing s3k 410/1,186; saved standing w2lvar 888/7,866; and
W2L108 GENTLE108 108/216. The reviewed Syn6 source is the designated final
compiled graph with SHA256
`697678b3281d7fb3b2a9322c2bc77b0671e7d05beb932fcb4364b8aa5ddaddea`.
It remains separate from the older standing Syn6 files.

Seven Python checks pass. The frontend checks also pass for all four libraries:
2,292 / 1,162 / 2,728 / 164 top-level actual connections are represented by
15 / 12 / 20 / 14 view wires, respectively. Bundling changes no serialized model
field. Selecting a view wire exposes its individual member connections.
The new app loads W2L108 and exposes all bundle counts in its accessibility
tree, with no JavaScript errors logged. Desktop preview activation is limited
by the Windows automation environment; this is separate from loading the file.
The user's older editor window and its unsaved edits were left open.

No neural or body integration was executed for this update. Existing original
model sources were read as saved data. Unsupported native/ionic model dynamics
remain preserved source records; this does not add a whole-controller simulator.

## Original first milestone

Verified on EB475WS4 with the separate ConnectomeEditorRuntime environment.

- Native Save As created
  `D:\Users\Ben Bolen\Connectome_Editor_Projects\MN_desktop_roundtrip.sns.json`.
- The desktop app reopened that copy with its complete graph and provenance.
- A synapse was selected in the desktop canvas; maximum conductance on
  `w2lvar:edge:7250` was changed from `0.6000000238418579` to `0.625 µS`.
  Ctrl+S saved it while the numeric field remained focused.
- A separate read reopened the saved file, checked source/interface records,
  node and edge identities, chemical endpoints, and physical maps, then compiled
  the edited graph with SNS-Toolbox. The changed conductance reaches the built
  network. Census: 18 neurons, 30 synapses, eight external current inputs, two MN
  outputs, 2 ms clock, and zero integration steps.
- Six Python tests passed: complete file round trip and Save As isolation,
  external-edit protection, invalid-source rejection, parameter consumption,
  nested boundary/input resolution, and explicit unsupported-model rejection.
- The frontend check passed for nested grouping, save/reload, ungrouping,
  source-field preservation, subtree ID remapping, explicit input bindings, and
  saving the current value of a focused parameter field.

The original `spinal/connectome_block_editor.html` has no changes from this
milestone. No controller or body simulation was executed. The final native
Build-button visual check was interrupted by Windows access failure
(`GetCursorPos failed: Access is denied. (0x80070005)`); the separate saved-file
builder verification above passed. The desktop window may retain the prior
generated page until the next launch. Relaunch to load the latest editor fixes.

Remaining scope: complete controller-specific import adapters, gated/spiking
builders, execution of physical maps against a body model, and Simulink
generation. This milestone does not designate new accepted controller results.
