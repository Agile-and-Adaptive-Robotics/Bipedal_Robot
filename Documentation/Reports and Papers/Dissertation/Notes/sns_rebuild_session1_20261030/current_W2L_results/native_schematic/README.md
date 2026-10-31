# Current W2L native-style publication schematics

**One main figure**, with optional expanded left-leg detail. The remaining views and tables are inspection companions, not a proposal to add seven W2L figures to Chapter 6. No source project, simulation result, manuscript, or active figure folder was changed.

Common folder: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\native_schematic`.

## Main and optional appendix assets

- Main landscape vector PDF: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\native_schematic\W2L108_MAIN_SUBNETWORKS.pdf`.
- Matching editable SVG and preview: `W2L108_MAIN_SUBNETWORKS.svg` and `W2L108_MAIN_SUBNETWORKS.png` in the same folder.
- Optional expanded native left-hip vector PDF: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\native_schematic\W2L108_LEFT_SUBNETWORKS.pdf`; SVG/PNG alongside. Full left knee/ankle detail is `W2L108_NATIVE_LEFT_KNEE_ANKLE.pdf` if required as a companion.
- Proposed caption: `CAPTION.tex`.

The main figure shows two actual four-cell RG subnetworks, the inherited hip PF/MN core with RG input stubs, all 12 individually identified homonymous II relays with visible arrows to their named MN targets, the commissural page, and the left contact page. Each added cross-page stub is backed by its actual executed connection ID. The four retained phase-related II channels and their RG/PF targets are explicitly listed beneath the figure. The right contact counterpart remains in the companion atlas and complete connection table. No biological edge is invented.

The **final main figure is 9 × 5.5 inches**, reflowed to fit the ETD sideways text block plus caption, with **10 pt minimum text at that actual print size**. `MAIN_LAYOUT.json` records this final geometry. Optional expanded views retain their 11 × 8.5 inch canvases and 12.5 pt labels; the named-type audit PDF is a separate full-size table with 10 pt text. Short physiological labels replace long GUI text; saved names and IDs remain in the tables. Native neuron centers and internal explicit route points are retained within each selected panel. Neural/off-page glyphs are enlarged 40% around their unchanged centers to fit print typography; physical sensor and adapter glyphs retain their native sizes. The clearly abstracted RG input and II-to-MN stubs have short display routes and actual executed endpoints. Panels remain separate subnetwork extracts, not a claim that this is the complete 108-neuron sheet.

## Exact source and executed configuration

Source project: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer.aproj`.

Executed graph: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\GENTLE108_EFFECTIVE_CONNECTOME.json`. This is the current lower-hip-recruitment configuration, with 108 neurons and 216 chemical connections. The effective hip PF-to-MN conductances are 2.1538912554 µS (extensor) and 1.0769456277 µS (flexor), bilaterally. The original 39.24% BW fitted trial is not rewritten by this display.

`EXECUTED_EDGES.csv` contains all 216 effective chemical edges, their full source and target IDs/names, exact reversal potential, maximum conductance, presynaptic threshold, saturation, saved type name, and target-change flag. `EXECUTED_NAMED_TYPES.csv` and `.json` group the exact executed settings by saved named type; the same saved type can appear more than once when the runtime applies different recruitment factors. The PDF table is a human-readable companion with rounded display values; CSV/JSON retain full precision. `NEURON_NAMES.csv` provides all 108 original names and physiological publication aliases.

The named type's saved `MaxSynapticConductance` is distinguished from an individual native connection's `G` override. `EXECUTED_EDGES.csv` includes both the executed `gmax_uS` and the unscaled `reference_connection_G_uS`, so that a pre-existing connection-specific value is not mistaken for a runtime fitting change.

The ankle execution overlay preserves the original connection IDs and routing control points while moving only their destination to the actual executed MN. Purple marks these four exchanges; the main overview's hip core has no ankle target overlay. The right and left exchanges are both explicit in `EXECUTED_EDGES.csv` and `SCHEMATIC_PROVENANCE.json` even when only the left ankle page is displayed. E/F are named neural pools; actual ankle mechanical actions must be checked against the actuator map rather than inferred from the labels alone.

Native saved inhibitory circles and excitatory forks are retained. Where the AddFlow drawing omits an `ArrowDst` element, the renderer uses the saved logical link's `ArrowDestination` style. The final audit found **zero** chemical heads requiring completion; the new II-to-MN links already have the appropriate native fallback. The green source node named `6` is labeled `Toe sensor` after verifying its body ID and its ContactCount adapter to the toe-contact neuron. Its original name and ID are retained in the mapping.

## Inspection and reproduction

`SCHEMATIC_PROVENANCE.json` records source hashes, panel node/edge selections, publication aliases, route-point preservation, completed arrowheads, and executed target overlays. `W2L108_NATIVE_SCHEMATIC_ATLAS.pdf` collects all inspection views; it is a companion rather than the main-chapter figure.

`TYPE_II_PHYSICAL_AND_NEURAL_MAP.csv` records all 16 type-II channels: 12 new homonymous relays and four retained phase-related channels, with native physical source names/IDs, the measured data type, and all actual neural targets. The final audit confirms all four executed ankle target changes, zero invented chemical heads, and unchanged source-project bytes.

Generator: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\make_native_schematics.py`. Run with `D:\Anaconda\envs\myo\python.exe` on **easteregg2**. It reads the saved native project and executed graph and writes only these staged assets. No simulation, fitting, native project editing, or physical-model change occurs.

The original Figure 6.1B raster has not been replaced or declared to match an older editable project. Its vector recreation remains contingent on a verified correspondence with that historical shared-RG source.
