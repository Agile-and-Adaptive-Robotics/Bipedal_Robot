# Desktop milestone verification

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
