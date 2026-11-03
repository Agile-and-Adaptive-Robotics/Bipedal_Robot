# v21 biology and neuroscience review

Independent review from the biology/neuroscience committee perspective, not an endorsement or statement by Radhika Reddy. Reviewed 2 November 2026. No scientific sources, controller settings, or simulations were changed.

## Second-render recheck: current disposition

Individually reinspected final **64, 65, 66, 68, 69, 70, 71, 77, 125, 128, and 129**, and added individual visual inspection of inherited questions **100–111**. The earlier findings below are retained as the first-pass audit; the dispositions here supersede them.

**Scientific corrections pass.** Final 68's native editable schematic has the same four cells and twelve incident connections as `FIGURE_EDGE_AUDIT.json` and `SYNAPSE_FIELDS_FULL_PRECISION.csv`: four excitatory hip/knee inputs, two integrator-to-own-relay excitatory edges, two relay-to-opposing-integrator inhibitory edges, and four excitatory motor-neuron targets. The two inhibitory terminals correspond to the source's −5 mV reversal; the other ten correspond to +8 mV. The schematic preserves the intended endpoint meanings. Final 69 now plots the actual saved `MN_bifemlh_r` response, so the HAM-to-motor-output example follows the circuit just explained. Its generator reads the archived aligned voltage array and retains the bilateral hip-MN row scale; it does not normalize the amplitude or integrate the model.

**Remaining required fixes in the second render:**

- **64:** “Commands” still wraps as “Command” plus a dangling `s`. Its other label is now the clean “Feedback.” Use a shorter label or a wider box.
- **68:** “Semimembranosus” similarly wraps onto a lone `s`. Fit that word deliberately without changing its spelling. The rebuilt diagram is otherwise much clearer and no longer distorted or clipped.
- **66, 70, 71, and 125:** correct proportions and external legends have made the plot labels too small for the planned auditorium. In particular, 70/125's legend and tick labels are approximately 13 pixels in the 1600-pixel-wide render, smaller than the source note. Their plot should be regenerated for the actual narrow/tall panel and with larger fonts. The 66/71 labels are approximately 16 pixels and should also increase. Enlarge plot text while preserving data, axes, colors, and line styles; enlarging or distorting the finished image is not the solution.
- **106:** the Left/Right legend overlaps the first pathway-removal bars. Move it outside the axes or into an empty region that does not cover data.

**Resolved or acceptable after recheck:**

- **65:** four method labels and the citation render without clipping; method attribution remains correct.
- **66:** the new plot restores Baseline/Stimulus identities, numerical time ticks, and the onset-aligned shaded interval. Only the projection text size remains.
- **68:** the circuit uses readable muscle names and defines RF/PF/HAM. Exact connections are checked as described above. Only the one long-word wrap remains.
- **69:** the inappropriate VI/BFSH example is replaced by the actual BFLH motor-neuron trace. All three panels are legible and undistorted, although larger headings would be welcome.
- **70/125:** the legend no longer covers the force peaks. Assistance values and scope are unchanged and correct. Plot readability remains inadequate at its current size.
- **71:** the PF ambiguity is removed: the header says pattern formation, and the legend says plantarflexor/dorsiflexor. The plot is undistorted and its legend is external. It still reports clearance counts rather than plotting geometric clearance directly.
- **77:** no visible collision remains between the source note and “Thank you.” The supported-result scope is intact. Reuse of older recap plot styling remains a nonblocking consistency recommendation.
- **128/129:** the historical motor/sensory crops now have natural proportions. Their labels and cropped-connection context are usable for questions. The source note still identifies historical w2lvar standing data, preventing confusion with ground9.

**Inherited question-slide coverage (100–111):**

| Final slide | Visual assessment |
|---|---|
| 100 | Readable overall. Minor inherited crowding/overlap around the lower summing-node label; fine signal wires will be difficult to follow on a projector. |
| 101 | No clipping. Long code-style internal title and small variable-name axis labels remain; the cleaner implementation comparison on 126 is easier to present. |
| 102 | Text is large and clear. Historical identity is explicit. |
| 103 | Table is readable and fits cleanly. |
| 104 | Replay-count graphic is readable and fits. |
| 105 | Local-perturbation graphic is readable and fits; support and local scope remain visible. |
| 106 | Legend covers the first bars; fix noted above. |
| 107 | Full source plate fits, but small source labels are low-resolution and unsuitable for a detailed projected walkthrough. Keep as an overview or supply a close-up if discussing individual edges. |
| 108 | Four panels, stimulus identities, time axes, and phase-shift statement are legible. No clipping. |
| 109 | The central four-panel result is too small for detailed label/legend reading; the slide works only as a recap. Use a full-size source slide for quantitative questions. |
| 110 | Scope labels correctly distinguish imported rat-model organization from the separate reduced-circuit test. Small circuit subpanels cannot support a detailed projected explanation at this size. |
| 111 | Large two-column text fits without overlap. The numerical-versus-pathway distinction is clear. |

Anchor hashes for this reviewed second render (SHA-256): `slide_064.png` = `85495C48A22BC9F5C6E96D9D374BA6415BE6F8C9126580198374F946D6137971`; `slide_068.png` = `4B4094BAF290E5E05079D8ADE436FBE0BB0DAB585D71EFC667A0A6132CE785E9`; `slide_070.png` = `EC30F550468BA4DAD930911125798B677532FD2C8EFD291E3ADBBB201F00DD15`. Later rerenders need a new visual check before calling these issues resolved.

## Coverage and disposition

Read `neural_content.mjs`, `build/neural_source_plan.md`, the staged campaign results/committee records, selected-controller provenance, support decomposition, the paired runtime comparison, the RF/HAM incident-edge source, the contact-strength caption, and the newly delivered EB475WS4 schematic provenance. Inspected the native PowerPoint PNG renders individually for final slides **57–77 and 124–129**. Inspected the RF/HAM, contact-comparison, and neural-hierarchy source graphics separately.

This is the **first candidate render review**. The numerical and controller-identity audit passes within the stated scope, but the visual/circuit-interpretation issues below require correction. Animation playback was not inspected by this reviewer; the inspected images are native static renders. Slides outside the listed range are not covered by this review.

## Corrections required before delivery

1. **Final 68 / source 65: preserve the diagram proportions and remove clipped material.** The complete RF/HAM SVG is stretched horizontally into a wide slide rectangle. Text and node shapes are visibly distorted, and the bottom legend is clipped on both sides. It also retains a long internal heading, E1–E12 identifiers, code-style muscle names, and separate SAR/GRAC details. For the main talk, use a faithful landscape rearrangement or a focused diagram containing the actual two-joint inputs, four-cell module, and its output targets. Keep the complete incident-edge plate available in questions/source assets. Define **rectus femoris (RF)** at first use and use readable muscle names. Preserve all scientific edge meanings; do not replace the circuit with an invented generic schematic.

2. **Final 64 / source 61: repair the two arrow labels.** “Commands” and “Sensory measurements” wrap onto dangling single-letter `s` lines. Widen the boxes, reposition the labels, or use “Motor output” and “Sensory input.” This is a visible native-render defect, not a subjective content preference.

3. **Final 66 / source 63: make the pulse plot interpretable.** The cropped panel shows the dashed Baseline key, but loses the blue Stimulus key. The source PNG itself has no readable numerical time ticks. The slide adds a time-axis title without restoring the numerical scale. Use the complete genuine source panel or reconstruct its axes from the source record. Show both line identities and a numerical time scale, with the 3.876 s pulse interval unambiguous. The source-circuit identity, +50 nA amplitude, and reduced-circuit scope in the notes are otherwise correct.

4. **Final 69 / source 66: avoid implying an incorrect specific serial pathway.** The preceding diagram's RF/HAM module targets rectus femoris, semimembranosus, semitendinosus, and biceps femoris long head. The rightmost trace panel instead shows vastus intermedius and biceps femoris short head. They are legitimate controller outputs, but are not the downstream muscles just drawn. Either show the actual RF/HAM target motor-neuron traces or explicitly frame the three panels as examples of activity at different levels of the whole controller. The existing broad headline can be read as a serial RF/HAM-to-VI/BFSH explanation.

5. **Final 71 / source 68: preserve the contact plot's aspect ratio and disambiguate PF.** The plot is visibly stretched horizontally. “PF” means pattern formation in its heading, while “PF” means plantarflexor in its legend. Spell out the muscle actions (plantarflexor/dorsiflexor), and move the legend away from the traces. These changes do not require changing colors or data.

## Further presentation improvements

- **Final 70 and 125:** the large in-axis legend obscures foot-force peaks. An external legend would allow the audience to see the alternating loading. The plotted and tabulated values are consistent with the source record.
- **Final 71:** the counts establish the saved clearance-gate result, but the plotted vertical force and pattern-formation voltage do not directly display geometric clearance. A clearance trace or a short corresponding motion comparison would strengthen the assertion-evidence pairing. Do not describe reduced foot force alone as proof of clearance.
- **Final 77:** “Thank you” occupies the source-note row. Remove it or move it into clear space. The final result remains appropriately scoped as supported Gait2392 control.
- **Final 75:** “Hold support settings fixed” is more precise than “Hold support fixed.” The applied support settings can be controlled while the measured restraint contribution changes with motion.
- **Final 127:** “loaded-arm model” would make the simulation status explicit in the headline. Its current caption correctly identifies SNS-Simscape component verification.
- **Final 72 and 77:** the recap reuses the older force-chart legend and diameter typography. Consider reusing the revised main-talk force plot so the recap does not restore the superseded busy presentation. Mechanical claims were not independently re-audited in this neuroscience review.

## Scientific checks that passed

### Narrative and attribution

The main sequence now introduces the suspended/released transfer problem, then a workflow gap, then the model library, block components, CAD mechanics, interface verification, functional calibration, literature circuits, and performance. New controller results follow those tools. This addresses the previous abrupt jump into results.

Hunt's four-part method is attributed to the **INRC2 Meeting, Intel Labs, 15 February 2019**, with the full presentation title in the notes: *Real-time control of dynamic walking in a dog-like robot*. The method is distinguished from the physical air/standing/ground curriculum. Notes explain that intrinsic/timing settings were reused and checked in the new campaign, rather than claiming every step was newly fitted.

### Controller identity and physiology

The selected new result is the human Gait2392 graded **biw2lvar-Hunt-PT** controller with a four-cell RF/HAM module per leg and separately enabled autogenic Renshaw repair. It is not W2L108, historical Syn6, or the separate reduced Deng stimulus experiment. The source retains 1,033 actual cells and 92 muscle outputs; extra logging names are aliases, not extra biological neurons.

The source RF/HAM diagram matches the recorded local graph: hip-flexor and knee-extensor pattern-formation input to RF; hip-extensor and knee-flexor input to HAM; integrators excite their own relays; the relays inhibit the opposing integrator. Its actual muscle-output targets support the correction requested for final 69 above. Voltage traces are labeled model mV, not EMG or recorded human physiology. Engineering model behavior is not presented as a direct validation of human circuitry.

### Support and movement

The chosen S4 record contains **5 right / 4 left qualifying unloaded, whole-foot-mesh-clear swings during 5–13 s**. The displayed GIF covers a shorter 7.663–13.000 s interval; the counts correctly refer to the assessment window in the notes. The result is guided stepping in place, not unsupported walking or established forward progression.

The reported support is correct: **48.848356% BW net vertical assistance**, comprising **40.051043% prescribed lift plus 8.797312% signed vertical restraint**, with **115.560892 N horizontal restraint RMS**. The denominator is the actual pelvis-subtree weight of 768.756726 N, including the modeled follower frames. The main slide reports 48.8% and 115.6 N; the question slide supplies the decomposition. Later support reduction is presented as future work.

### Contact feedback comparison

The quarter-strength condition passes 5/4 qualifying swings, while the original contact-to-pattern-formation strength passes only 2/2. The source identifies toe inhibition of the ankle dorsiflexion pattern-formation cell. The statement is explicitly limited to this test, not generalized as a biological optimum. The next lower-lift test failed; no minimum-support or harness-free-walking claim is made. The planned stance/swing comparisons are clearly separated from completed results.

### Runtime and validation boundaries

The table correctly separates the paired **complete-process 199.623573/29.284280 s = 6.817×** improvement from **neural-call 179.079112/9.779516 s = 18.312×** improvement. The paired source trial is distinct from ground9, and the notes do not imply faster-than-real-time fully instrumented ground execution. The identical recorded-trajectory statement is supported by the source parity audit.

The short MuJoCo–Simulink readback verifies the interface, not walking. The W2L NumPy/Simulink comparison verifies the driven neural-core response, not the selected detailed-human controller. The loaded-arm example verifies a component response. These distinctions are retained in the notes and source captions.

### Newly delivered EB475WS4 diagrams

Final 128–129 are correctly labeled as **historical w2lvar standing-export** close-ups. They are not substituted for ground9's selected controller. AD is defined as an adapter/map, RF as rectus femoris, and MN as motor neuron. The source audit identifies the gray/dashed zero-effective PF-to-MN paths and distinguishes physical signal maps from neural connections. The schematics preserve Ia and II identities; their spindle symbol is an encoded length/rate source, not proof of a complete receptor model.

## Evidence locations

- `assets/neural/sources/RESULTS_AND_FILES.md`
- `assets/neural/sources/COMMITTEE_QA.md`
- `assets/neural/ground9_rc_s4/SUPPORT_DECOMPOSITION_5_TO_13.json`
- `assets/neural/ground9_rc_s4/PRESENTATION_NEURAL_PROVENANCE.json`
- `assets/neural/ground9_rc_s4/REPLAY_PROVENANCE.json`
- `assets/neural/benchmark/BACKEND_AND_VALIDATOR_COMPARISON.md`
- `assets/neural/rfham_selected/RFHAM_FOURCELL_SELECTED.svg`
- `assets/neural/slide_panels/NEURAL_RIGHT_HIERARCHY.svg`
- `assets/neural/slide_panels/CONTACT_PF_AND_FOOT_LOADING.svg`
- `assets/neural/eb475ws4_motor_loops/`
- Campaign originals: `D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/`
- EB originals: `D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/functional_subnetwork_review/drawings/full_motor_loops/`

The earlier `build/neural_source_plan.md` statement that the EB drawings had not been found is superseded by the actual included assets and final 128–129. Update that planning note if it will be delivered as a handoff.

## Final-render closure

Individually viewed final slides **64, 66, 68, 70, 71, 77, 106, 125, and 128** from the fully rendered final candidate. **All required findings from this review are closed. No delivery-blocking biological, controller-identity, or visual issue remains within the reviewed scope.** This closure supersedes the earlier open dispositions above.

- **64:** the label is now “Motor output,” with a deliberate word break and no dangling letter. Feedback and the interface comparison remain clear.
- **66:** both line identities, numerical time, stimulus interval, and readable plot labels are present. The source-architecture and reduced-circuit distinction remains intact.
- **68:** “Semimembranosus” fits correctly, as do the other muscle names. The four cells and twelve ledger-verified connections are preserved, with clear excitatory and inhibitory terminals.
- **70/125:** the revised panels use larger labels, fewer readable axis ticks, and external legends. Neither legend obscures data. Net support, its decomposition, and horizontal guidance remain correctly reported.
- **71:** the enlarged comparison has full muscle-action names, external legends, clear row labels and units, and preserved quarter/original-strength identities. The stated swing counts remain visible and correctly scoped.
- **77:** conclusion content fits without source-note or footer overlap; the neural result remains explicitly supported.
- **106:** the legend is outside the plotting area and no longer covers pathway-removal bars.
- **128:** the natural-aspect historical crop is legible, its adapter/motor definitions are visible, and it remains labeled historical w2lvar standing export.

The parent reports byte identity for unchanged slides; this closure does not claim new visual inspection of those unchanged files. Previous bounded coverage remains **57–77, 100–111, and 124–129**, with the specified corrected slides rechecked. Acceptance concerns the native static renders and scientific records; no live animation-playback inspection is claimed.

Final reviewed render hashes (SHA-256):

| Slide | SHA-256 |
|---|---|
| 64 | `CB39A9AEED421BDC30F92332F4F3B640381BF1FB4113FF4688346BBF13CF7741` |
| 66 | `E031C664E24A00D2E65D3AD98177310E31ABDB172981E7F6F140FEA022845CAB` |
| 68 | `07201E2833D999AAE52854C6A98CCB57EBAD7132B30D0EAEDAFE6B443BCBAE61` |
| 70 | `18E89538755398F464061D805D48802D7B2D9367811285CA3FCF567A969FA142` |
| 71 | `BBB6DD90780001867521C23584433E633837EC5993E95F4F51C9CDC1DFF1AE0E` |
| 77 | `8944243CE7AF5E25A186E18BBE1BCB0A0D096EFBCA29DF5A781E383A35CB039C` |
| 106 | `D59D2407E290D68EBABCE573C752E1858B6F90ABDEEBE3B4DE237DDFBF88C4A0` |
| 125 | `9CEFF145C11BA709ABCFE84F3D702820F79D35894E149B11C122D388B44B7B40` |
| 128 | `B397D6F6CD28E826EDDF4A06949C22EB84AC1661E7270A33DE0DA770204CD91F` |
