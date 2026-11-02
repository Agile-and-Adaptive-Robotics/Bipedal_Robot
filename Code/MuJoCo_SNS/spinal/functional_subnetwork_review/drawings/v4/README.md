# Grouped circuit schematics

Open **RG_PF_Schematics.pdf**. Each drawing has an editable SVG, a PNG preview, its exact saved-source connection records, and a geometry check beside it.

1. W2L108 graded central feedback, including opposite-side c1 and V3.
2. W2lvar complete central-feedback overview.
3. W2lvar Ia/II view, with the same module locations.
4. W2lvar Ib/contact view, with the same module locations.
5. Syn6 RG and six PF channels shown together.
6. S3k's actual descending/RG/PF architecture, replacing the previous blank absent-feedback page.
7. Rat134 with separate hip, knee, and ankle PF modules.
8. RF/HAM four-cell PF and its joint inputs/motor outputs, from the saved Session1 left-side projection.
9. W2L704 hybrid, with actual population sizes on the displayed logical interfaces.
10. Gait2392 TimingZeroSFA central-population interfaces; RG/PF-S populations contain16somas, while displayed INs remain graded.
11. Destination-by-destination synapse parameter explanation, including checked left/right values.

Rounded boxes, distributed source ports, semicircular crossing bridges, and triangles with their broad base against the receiving object follow Ben's annotations. A crossing is not a junction. Internal RG motifs are omitted. The RF/HAM page expands its specific four-cell motif. No duplicate destination HCs or absent-route boxes are used. Gray group links on W2L/Rat bundle the two actual extensor/flexor RG-to-PF projections; the individual source edges remain in their connection JSON.

The drawings show selected saved conditions, not newly modified circuits. In the saved W2L108, ankle afferents project to local motor/reflex cells and heel projects to central extensor HCs; the expected additional central ankle/heel→IN routes remain a source discrepancy, not invented arrows. Syn6's disputed saved heel→flexor-IN connection is marked red, and its saved physical contact-input gains remain zero. The RF/HAM projection retains its historical left-side weights and topology; it is not presented as a current tied bilateral winner. These source issues are not claimed resolved by drawing repairs.

`VERIFICATION.json` records saved-edge parity and geometry checks. No circuit, calibration, training, optimization, physics, native model load, or live application was run or changed. Python imports in the drawing pipeline load drawing helpers only; controller source snapshots are never imported or executed.
