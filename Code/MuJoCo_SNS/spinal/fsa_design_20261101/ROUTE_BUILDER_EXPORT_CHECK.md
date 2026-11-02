# Exact route-builder exports and central wiring

Session 1, 2026-11-01. Read-only verification on **EB475WS4**, hostname checked through SSH, followed by comparison with easteregg2. No remote files changed.

The three files are `D:\Github\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Gait2392_standing_families\{family}_graph_ledger.json` on EB475WS4. All three are byte-identical to the corresponding easteregg2 exports. These are saved standing-condition ledgers, not standalone model files or deafferented air-run exports.

| Export | Neurons / nonzero chemical edges | SHA256 |
| --- | --- | --- |
| syn6 | 794 / 2966 | 0843e32e4e5cf16de974720620c0101e05bf1424018c7815d21b8b1182d6e5b0 |
| s3k | 410 / 836 | 0e232bd04001e342a774d3de08ba5d0ed378907673bc2c5ccfa81daae22ae743 |
| w2lvar | 888 / 7570 | 49b946c627222e030aa7bc741d7de20faeb804646d791b51548e6d3d704d4e2f |

## Direct afferent-to-central edges in the actual exports

This count selects physical Ia_/II_/Ib_ source populations and RG_, InE_, InF_, or PF_ targets. It does not count local motor reflexes or substitute a pooled relay for a direct edge.

| Export | Ia | II | Ib |
| --- | ---: | ---: | ---: |
| s3k | 0 | 0 | 0 |
| syn6 | 0 | 0 | 0 |
| w2lvar | 152 | 152 | 216 |

Local sensory/reflex circuitry remains in s3k and syn6. These numbers do not mean that the entire controller is deafferented.

In all three exports, global_G has full_rules, ia_f_central, ii_f_central, ii_e_central, and ib_e_central equal to zero. W2L uses separate variant-central gains, so those global zeros do not describe its active central routes.

The full s3k builder implements conditional central branches at build_network.py around lines 800–851. The zero gains omit those branches at construction. Calling this a full source builder therefore does not guarantee the full desired topology is enabled. The syn6 builder has an optional extensor-Ib central block around lines 744–780; it lacks corresponding general Ia/II central implementations in that block. This is more than a drawing filter or air-run mask.

## Heel wiring

The syn6 export contains, on each side, HEEL→InE and HEEL→InF, both gmax=0.5 microSiemens and reversal=+8 mV. build_network_syn6.py lines 568–571 hardwire the two excitatory projections. In the same builder, InE inhibits RG_F and InF inhibits RG_E. Thus the two heel paths can recruit inhibition of both ipsilateral RG half-centers. They do not implement a selective extensor-favoring reset by themselves. A common inhibitory reset would require its own declared pulse protocol and functional validation; that was not established by the earlier numerical checks.

The selected syn6 caller has heel_rge=toe_rge=ib_rge=0. Those zero input-current scales silence the contact graph in the saved standing test. They do not correct the graph's intended function or justify enabling it unchanged.

W2L's heel organization differs: ipsilateral InE excitation and contralateral InF excitation. Trace actual targets: the former inhibits the ipsilateral flexor RG; the latter inhibits the contralateral extensor RG. A source comment describing the latter as suppressing the partner flexor is inconsistent with that actual wiring.

## Why earlier checks were insufficient

The earlier checks established conductance/current equations, units, local gains, bounds, cache consumption, and accelerated/reference agreement. They did not establish that every expected physiological pathway was present or that the full signed chain produced the intended phase response. An inventory saying zero qualifying routes must remain a design gap, not a functional pass.

Historical exports and results remain unchanged. A corrected heel policy and an explicit typed central-feedback design are separate new conditions requiring individual signed-path and combined-function checks. No training, body simulation, tendon-compliance change, or damping change is authorized by this audit.

## Fourth file: exact W2L108 AnimatLab export

The supplied EB475WS4 file is `D:\Github\Bipedal_Robot\tmp\functional_subnetwork_remote\W2L108_verified\editable\W2L_VirtualWalker_Newer_Standalone.asim`, SHA256 `f3d1fd84edeec21e2b033495bfc81bd016525d029108341f35be1c5f28c11bc4`. An unchanged read-only audit copy is preserved at `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\remote_w2l108\EB475WS4_drawing.asim` on easteregg2.

Its byte hash differs from the canonical easteregg2 source (`01ece6edcbcf17c274410b23f1ce7603eab288aadb946d6b548c5102ed5c1154`), but complete XML comparison finds exactly one leaf difference: ProjectPath. All 108 neuron IDs and every neuron field, all 216 Connexion IDs and fields, and all 47 synapse-type IDs and fields are identical. Physical adapters and all other XML leaves also match. Thus the drawing file is the same source model relocated into the editable project directory, not a different neural controller. Runtime tuning profiles remain separately specified; they are not retroactively encoded in this raw source export.
