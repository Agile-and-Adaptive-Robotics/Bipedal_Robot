# W2L spiking/graded controller: editable schematic

Machine: **easteregg2**. Selected saved live run:
`D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\W2L_VirtualWalker_Newer\pace_0p2_supported10_v1`

The controller has four RG and eight PF populations of 16 calcium/adapting LIF cells, and 16 sensory populations of 32 LIF cells (12 homonymous II and four heel/toe inputs). That is **704 spiking and 80 graded membranes**: 108 logical channels and 216 chemical pathways. Motor neurons remain graded. The 28 population spike filters/readouts are computational interfaces, not additional neurons or independently stored cellwise synapse matrices.

## Open the drawing in the block editor

Preserve your current tab first. Use **Templates → Blank** to create a new tab, then **Load JSON from file…** and select:

`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\W2L_704_80_baseline\W2L_704_80_hierarchy_connectome.json`

Double-click a SUB node to inspect its circuit. Double-click a population to inspect its actual 16 or 32 cells and nonneural input/readout blocks; F4 returns to the parent. Individual schematic files also load directly:

- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\W2L_704_80_baseline\W2L_704_80_left_RG_PF_connectome.json`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\W2L_704_80_baseline\W2L_704_80_left_knee_pair_connectome.json`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\W2L_704_80_baseline\W2L_704_80_exact_connectome.json`

**W2L_704_80_workspace.json is an archival multi-tab snapshot, not a direct input to the editor's single-tab Load button.** No browser, editor, or localStorage state was changed to produce these files.

## Exact graph and editable figures

`W2L_704_80_graph_ledger.json` preserves the actual compiled matrix, 784 membrane memberships, 44 physical adapters, parameter units, 12 motor outputs and anatomical moment-arm signs, runtime ankle target swaps, gain calibration, and pooled hip force-law settings. Full original source names and IDs remain alongside any display alias. `W2L_704_80_chemical_edges.csv` and every editor edge tag preserve the four full-precision effective synaptic fields: conductance, reversal potential, presynaptic threshold, and saturation. The PDF/SVG synaptic-type table lists 49 effective variants with their original source-named types. Runtime calibration can give several effective variants of one named type.

The overview, left RG/PF circuit, knee motor/reflex pair, and four-cycle result are provided as editable text/vector SVG, PDF, and PNG. Source fonts are at least 10 pt at approximately thesis-column width. The source graph is unchanged; these are schematic exports, not executable controller replacements. The hierarchy uses explicit pathway counts, not summed conductance, for inter-subsystem bundles. The full graph is authoritative; publication detail panels are induced subnetworks with boundary paths retained in the hierarchy/full graph.

## Source-specific details

- One source Renshaw cell has a duplicated L label but actual right-side neighbours. Its display alias explicitly marks the right-side role; the original name/ID and numerical connections are preserved.
- RG/PF traces show filtered population spike activity, not individual soma voltage. Source calcium-bearing Renshaw cells remain graded.
- In this source Hill implementation, Ia is gain times series length, II is gain times parallel length, and Ib is gain times tension. Do not silently substitute a generic Ia=velocity equation. Physical Ia adapters enter Ia-IN cells in this model rather than its separate Ia-labelled cells.
- Heel/toe inputs are geometric contact levels, not GRF or heelstrike impulses. Actual GRF, support, and clearance are recorded diagnostics. This baseline adds no GRF-to-body-balance pathway.
- Source ankle MN ext actuates dorsiflexion; source ankle MN flx actuates plantarflexion. Four explicitly recorded runtime PF-target swaps are retained. Physical route names alone are not anatomical action labels.
- Excitatory/inhibitory arrow categories are source-type shorthand. Actual instantaneous chemical current is g(Erev−Vpost), so its sign depends on postsynaptic voltage.

## Saved-data result

`FOUR_CYCLE_PROVENANCE.json` identifies the fixed 3.842–7.488 s window: four intervals between the first five left-RG-F activity peaks after 3 s. Data are from the actual live neural controller, not the later offline motor-wave replay. Panels show RG/PF spike activity, graded hip MN voltage, anatomical hip/knee/ankle angles, foot normal forces, separately summed upward/downward support, and minimum oriented main-foot/toe clearance. Large impact peaks are retained. The whole 10 s run has 10/9 unloaded, cleared swings and 36.48% BW net support; passive rails/guides remain. This is supported stepping, not independent walking or a periodic-orbit result.

Regenerate on easteregg2:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\export_baseline_network.py'
```

`EXPORT_VERIFICATION.json` records exact census, saved gain-ledger agreement, valid native editor types/endpoints, and full-matrix export checks. Browser import was not performed.
