# Drawing revision against Ben's 34 comments

Current drawing pack: `RG_PF_Feedback_v3.pdf` (four pages: grouped W2L108, grouped w2lvar, syn6 RG/PF, and a destination-parameter explanation). Editable SVGs and exact source-connection JSON accompany each page. No controller or native model was changed; no circuit was imported or run.

| Comments | Disposition |
| --- | --- |
| 1,4,5,10,11,26 | Internal RG/PF motif wires removed. Receiving INs remain visible in their own modules. No internal motif leaves the module boundary. |
| 2,6,7,17,18,19,21 | Source boxes have separate outgoing ports, larger rounded bends, and reordered tracks. No forced midpoint fan-out or small artificial up/down jog is used. Input bundle labels identify their source category and count. |
| 3 | Rounded source, neuron, and module boxes. |
| 8 | Crossing detector adds semicircular bridges and a white clearance under the overpassing stroke. Bridges are not junctions. |
| 9,16,31,34 | Removed all absent-route boxes and the blank s3k absent-feedback page. The missing source projections remain documented outside the figures rather than replaced with invented edges. |
| 12,13 | Opposite c1, heel, and V3 sources placed beside RG, rather than below the entire stack. |
| 14 | Actual `R c1 → L RG flx` source edge now included in W2L108. |
| 15 | Verified source difference: the selected W2L108 saved standalone has heel→extensor HC, not heel→extensor IN. Its V3→RG extensor IN is drawn. W2lvar explicitly shows Ib→InE and Ia/II→InF. No unsupported W2L108 heel→IN edge was inserted. |
| 20 | Removed duplicate HCs and the incorrect “off-page” labels. Each destination is drawn once. |
| 22,28,32 | Grouped extensor/flexor layout restored. RG E is left, RG F is right. W2lvar's source-named `PF_KNEE` pair is labelled K&A, consistent with its actual shared knee/ankle builder. |
| 23,24 | Verified source discrepancy: ankle sensory cells exist in W2L108, but their saved outputs go to local motor/reflex circuits, not RG/PF. No invented central ankle path was drawn. W2lvar bundles explicitly include ankle channels. The selected saved source still differs from Ben's expected W2L ankle-feedback diagram; this is not reported as a completed circuit correction. |
| 25 | 45-degree transitions are used for the opposite-V3 routes; long connections otherwise use simple orthogonal tracks. Arbitrary shallow angles removed. |
| 27 | Toe chain is explicitly `Toe contact → Toe relay (TOEDF) → PF destination`, with no ambiguous detached relay glyph. |
| 29 | Added `II_tib_ant_destination_parameters.csv` and the explanatory PDF page. Afferent-neuron identity is separate from each outgoing synapse's parameters. The four central edges share g≈0.0278947uS/Erev8mV in this saved condition; local relay inputs are1uS, excitatory relay→ownMN is0.4uS, and inhibitory relay→antagonistMNs is≈0.133333uS/Erev−5mV. Exact left/right homologues are checked per row. Side suffixes are omitted from visual labels and retained in source tables. |
| 30 | This is present in the saved syn6 source, not merely a plotting origin mistake: `HEEL_l → InE_l` and `HEEL_l → InF_l`. The disputed latter branch is red and identified in the caption. No saved connection was removed or silently corrected. |
| 33 | Syn6's actual RG→S1…S6 projections and PF interneuron outputs now appear together, instead of isolated RG/PF fragments. Its PF channels are six NMF channels, not a renamed hip/K&A pair. |

Source hashes and represented edge fields are in DRAWING_MANIFEST.json and each page's connection JSON. The unchanged source mismatch for W2L108 heel/ankle paths and the disputed syn6 heel branch remain explicit findings, rather than assertions that all circuit concerns have been corrected.
