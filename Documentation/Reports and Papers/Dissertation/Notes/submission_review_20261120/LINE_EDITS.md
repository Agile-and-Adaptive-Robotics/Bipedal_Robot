# Dissertation wording and figure changes, November 20

## Author’s ethics wording and the search for “correct”

The live and local dissertation both have 327 pages. The project compiled in Overleaf with zero errors and zero warnings. This section records the latest changes. The earlier pass is documented below.

- **Chapter 2, section 2.10, printed pages 38–40.** Line 270 applies the same ethical standard to other animals as to healthy human volunteers. It defines the Three Rs and identifies robotic experiments as replacement. Removed the wording about the “justification” for animal experiments and testing in simulation “before corresponding animal experiments are considered.” Line 272 restores the author’s exact “without repeating these invasive preparations in living animals.” Line 286 retains “A synthetic nervous system therefore replaces future experiments.” The historical account, figure, and references remain in place. [Chapter 2](../../Final/chapters/15-background.tex)
- **Figure 1.1, page 10.** “Corrections, routes” became “Design factors and routes (Ch. 4).” Updated the figure, its SVG, generator, and local text description. The generator starts in embedded Arial. [Roadmap](../../Figures/10-introduction/roadmap.pdf)
- **Figure 7.1, page 167.** “Corrected torque model” became “Torque model with design factors (Ch. 4).” Updated its generator, PDF, and description. The Chapter 4 extensor figure description now also calls its identified parameters design factors. The artwork in both revised diagrams is identical outside the text that changed when checked at twice its native resolution. [Mechanical design workflow](../../Figures/70-overall/design_workflow.pdf)
- **Search results.** The check covered the active LaTeX sources, captions, tables, figure descriptions, and searchable text of the 327-page PDF, including vector labels and references. An English OCR pass over 53 distinct embedded images returned no further matches. It searched case-insensitively for the stems *correct* and *prepar*, including compounds and prefixed forms. The printed plural *corrections* occurs zero times. The five printed instances of the singular *correction* refer to editing literature summaries (page 132), the response to a perturbation and visual feedback in a diagram (both page 177), the recorded encoder-angle adjustment (page 230), and the stance response (page 306). The only plural *preparations* is the author’s sentence on page 39. The singular *preparation* occurs in authorship, publication status, and model or simulation setup. A source label for the active-length equation retains its internal identifier.
- **Appendix A, line 205.** The implementation sets the compliance displacement to zero when the relative strain evaluated before deflection is at least one. The earlier wording said “exceeds unity.” Corrected the inclusive boundary and widened the implementation reference to lines 426–441. The equations and fitted parameter values are preserved. [Evaluator, lines 426–428](../../../../../Testing_Data/2022_02_Festo/minimizeFlxPin.m)
- **Review scope.** Section 2.10 is reserved to the author. The historical advisor-style checker excludes it from stylistic comparisons and recognizes the author’s approved wording separately. [Review guidance](../advisor_feedback_history/README.md)

### Exact source replacements in this pass

- **Chapter 2, line 270.** `on healthy human volunteers. Tonic stimulation` → `on healthy human volunteers. I apply the same ethical standard to other animals. Tonic stimulation`.
- **Chapter 2, line 270.** `The same experiments in animals also carry real moral cost, and their justification is governed by the Three Rs framework: replacement of animal subjects with non-animal methods, reduction of the numbers required, and refinement of the remaining procedures \citep{russell_principles_1959, tannenbaum_russell_2015}. Robotic and computational models allow proposed neural interventions to be tested under controlled conditions before corresponding animal experiments are considered.` → `The Three Rs of animal research are replacement of animal subjects, reduction in the number used, and refinement of procedures \citep{russell_principles_1959, tannenbaum_russell_2015}. This work pursues replacement by using a robot or simulation as the experimental subject.`.
- **Chapter 2, line 270.** `Repeating a model does not create independent biological subjects, and its relevance to animal experiments depends on validation against biological observations.` → `Repeated runs describe one model. Their biological interpretation depends on how well the model reproduces existing observations.`.
- **Chapter 2, line 272.** `without repeating that surgery on another living animal.` → `without repeating these invasive preparations in living animals.`.
- **Chapter 4, line 95.** `The constant term $\chi_0$ corrects the unloaded length discrepancy` → `The design factor $\chi_0$ represents the length offset at zero load`.
- **Chapter 4, line 99.** `The BPA length corrected for this constant offset is` → `Subtracting this offset gives the unloaded BPA length`.
- **Chapter 4, line 245.** `Since $\delta$ is the scalar compliance displacement along the BPA line of action, then the compliance-corrected active BPA length used in the force balance is` → `The axial displacement $\delta$ of the brackets and tendon shortens the active BPA length in the force balance to`.
- **Chapter 7, line 51.** `I then used the corrected model to size two routes` → `I used the torque model with these design factors to size two routes`.
- **Appendix A, line 151.** `(zero in the uncorrected model)` → `(zero in the baseline model)`.
- **Appendix A, line 184.** `The compliance-corrected muscle length entering Equation~(\ref{eq:Lm}) is computed as` → `Before the attachments deflect, Equation~(\ref{eq:Lm}) gives the muscle length as`.
- **Appendix A, line 205.** `when the uncorrected relative strain already exceeds unity (\texttt{minimizeFlxPin.m} lines 434--441)` → `when the relative strain before attachment deflection is at least one (\texttt{minimizeFlxPin.m} lines 426--441)`.

The [source ZIP](../../Bolen_Dissertation.zip) passed its integrity check. Eight of its 165 files changed. The other 157 were verified byte-for-byte against the previous package. The [full PDF](../../Bolen_Dissertation.pdf) contains zero overfull boxes, unresolved references, or citation or font warnings. PDF SHA-256 `4e0257f068beba867a6bfe4ba3dd0b5c0e02455098c091530075f0f7cedda0fb`. ZIP SHA-256 `f29130a87c99ea6f2aeff801def9d68ef05871aefdd860dc69751bc24a1d840b`.

---

> The records below describe the preceding release.


## What the preceding revision changed

The earlier wording can be checked in the [source comparison](SOURCE_CHANGES.diff). Line numbers refer to the LaTeX source and page numbers to the printed manuscript.

- **Chapter 1, line 12.** “An artificial muscles” became “an artificial muscle.” “Interacts with limb, and better predict” became “interacts with the limb and to better predict,” supplying the article and putting “to develop” and “to better predict” in parallel. The word “and” was already present. “In (Chapter 4)” became “in Chapter 4.” Both versions use the same cross-reference, `\ref{ch:torque}`.
- **Chapter 1, line 14.** The added conjunction was “and” between “there is not a complete understanding of this” and “no comprehensive view has been written or collected from scientific papers.” A conjunction is a word that joins words, phrases, or clauses. Here it joins two independent clauses.
- **Chapter 2, line 272.** The last sentence in the submitted draft also appears verbatim in the preceding release. The current edit changes only the final clause of that sentence. Its old and new wording appears below.
- **Chapter 2, line 284.** The sole edit in the preceding revision was “Stuart reproduce” to “Stuart reproduces.” The current revision changes the description of the drawing in the last sentence. It retains the preceding account of the classroom film.
- **Chapter 2, line 286.** The previous opening meant that animal experiments supplied some of the evidence researchers now use to replace animal experiments. The revised paragraph says so directly, identifies the other kinds of work collected in SADb, and gives an example of a sensory intervention before the author’s chosen ending.

## Figures and printed pages

| Location | Revision and verification |
| --- | --- |
| Chapter 4, Figure 4.5, page 87 (page 107 in the PDF viewer) | Moved the second bracket’s axial projection closer to its force vector, rotated its displacement toward the force direction, and rotated the first bracket’s displacement farther from its projection. Each projected length is calculated from the dot product of displacement and force direction. Dotted guides show the projection. Displacements are schematic. [Editable figure source](../../Figures/30-torque/Components/bracket_projection/draw_projection.py) |
| Chapter 4, line 454, Table 4.4, page 98 | TeX reported two underfull horizontal boxes at the end of the table. Full justification stretched spaces in “Negative reward on the worst-angle margin” and the cam-regularizer row. The badness values were 2050 and 10000. Left-aligning the Term and Meaning columns on source line 445 removed both messages. All cells, weights, and column widths are preserved. |
| Chapter 4, Figure 4.9, page 106 | The accepted plot matches the previous source ZIP byte for byte. |
| Chapter 6, line 172, page 152 | The Gait2392 stepping study selects Optuna’s TPE sampler. The staged rhythm and sensory-timing trainer uses Rustuna. The clause defining TPE and the three bibliography entries appear below. [Gait2392 stepping search](../../../../../Code/MuJoCo_SNS/spinal/optuna_walk.py) · [Staged trainer](../../../../../Code/MuJoCo_SNS/spinal/controller_profiles_20261104/train.py) |
| Chapter 7, line 8, page 166 | The mechanical design workflow in Figure 7.1 (page 167) starts from a human torque target. The neuromechanical modeling workflow in Figure 7.5 (page 180) compares each altered trial with a run using the original controller and route. The editor already contained the author’s “so that a researcher,” which is also in this local release. |
| Chapter 7, Figure 7.2(C), page 171 | Removed the full sentence about sectioning a specimen from the caption. In the main text on page 169, the proposed model fills the space between the opposed aramid helices with hyperelastic rubber. Festo identifies aramid and rubber in the product sheet. The patent describes a rubber or elastomer spacer. The thickness of the interlayer in the tested specimens remains unmeasured. [Festo product sheet, p. 2](https://www.festo.com/media/catalog/202851_documentation.pdf) · [Festo patent, Figs. 2–4 and claims 2–5](https://patents.google.com/patent/US6349746B1/en) |

## Exact textual changes in this revision

Each entry quotes the changed words on the specified source line. The complete paragraph is given where it was replaced as a unit. The engineering term “design factors” applies to the identified length offset, stiffness, and wrapping factors. Appendix C still calls the separate encoder-angle correction a correction.

### Chapter 1

[LaTeX source](../../Final/chapters/10-introduction.tex)

- **Line 8.** `were examined in biological preparations` → `were tested in experiments on the salamander spinal cord`
- **Line 14.** `across species, preparations, and experimental protocols` → `across species, behaviors, and experimental conditions`
- **Line 34.** `identify and validate corrections for` → `identify and validate design factors for`
- **Line 50.** `Corrections for` → `Design factors for`

### Chapter 2

[LaTeX source](../../Final/chapters/15-background.tex)

- **Line 158.** `This preparation separates` → `Removing sensory input separates`
- **Line 189.** `depends on the preparation and step phase` → `depends on the experimental conditions and the phase of the step`
- **Line 272.** `The invasive preparations that motivate replacement are not always apparent in circuit schematics of decerebrate-cat locomotion, which emphasize the neural circuits under study but often leave the surgical preparation out of view.` → `Circuit schematics of decerebrate-cat locomotion emphasize the neural circuits under study and often omit the vivisection that made the experiments possible.`
- **Line 272.** `shows the physical preparation behind the circuit diagrams` → `shows the cat on the treadmill after that surgery`
- **Line 272.** `without repeating these invasive preparations in living animals` → `without repeating that surgery on another living animal`
- **Line 284.** `The treadmill preparations that established locomotor neurophysiology are likewise preserved in the historical literature: \citet{stuart_reflections_2007} reproduces the arrangement of the high-decerebrate cat treadmill experiments, with electromyographic recording from the hindlimb muscles, from the original Russian-language reports` → `\citet{stuart_reflections_2007} reproduces a drawing from the original Russian-language report showing how researchers recorded hindlimb muscle activity during treadmill stepping in cats after brainstem transection`
- **Line 286.** `The honest accounting is murkier than the Three Rs framing alone suggests, because the models a robot implements are themselves products of animal experimentation. The half-center organization, the two-level rhythm and pattern architecture, and the length, force, and load signaling of the afferents were each established in preparations like the ones described above \citep{brown_intrinsic_1911, severin_1967, rybak_modelling_2006}, and the databases of afferent properties that parameterize sensory feedback are compiled from animal recordings \citep{bolen_sensory_2023}. A synthetic nervous system therefore replaces future experiments; it does not replace, and cannot repay, the animal work on which its architecture is built. Recognizing this inheritance sharpens rather than undermines the case for robotic subjects: the knowledge purchased at that cost is precisely what makes a robotic replacement scientifically meaningful, and the robot is the instrument that keeps that knowledge productive without adding to the ledger.` → `Some of what we know about the spinal control of walking came from vivisection. Experiments in cats revealed that spinal circuits can generate rhythmic activity and showed how muscle receptors respond as the legs move \citep{brown_intrinsic_1911, severin_1967}. This evidence informs models of rhythm generation, muscle activation, and sensory feedback \citep{rybak_modelling_2006}. SADb brings studies of this kind together with research on human movement, computational models, and reviews \citep{bolen_sensory_2023}. With a validated model, a researcher can disable a sensory pathway and repeat the same locomotor experiment without operating on another animal. A synthetic nervous system therefore replaces future experiments.`

### Chapter 4

[LaTeX source](../../Final/chapters/30-torque.tex)

- **Line 13.** `develop corrections for` → `introduce design factors for`
- **Line 13.** `evaluate the correction terms` → `evaluate the design factors`
- **Line 17.** `the length, compliance, and wrapping corrections` → `the design factors for length offset, compliance, and wrapping loss`
- **Line 93.** `the corrections that were developed` → `the design factors developed`
- **Line 95.** `The correction terms` → `The design factors`
- **Line 253.** `a specified set of correction parameters` → `specified values of the design factors`
- **Line 284.** `The correction is phenomenological; it does not resolve wall stresses or the deformation of the braided actuator.` → `This empirical model represents the loss of usable length. Resolving wall stresses and braid deformation would require a structural model.`
- **Line 314.** `Identification of the Torque-Model Correction Terms` → `Identification of the Torque-Model Design Factors`
- **Line 316.** `I set the correction terms` → `I identified the design factors`
- **Line 322.** `For a candidate correction vector` → `For a candidate vector of design factors`
- **Line 342.** `the zero-correction baseline value of the same component, the same kinematics and force model with every correction term disabled` → `the value for the baseline model, which uses the same kinematics and actuator force model with the design factors disabled`
- **Line 342.** `the zero-correction baseline on` → `this baseline on`
- **Line 354.** `without the experimentally identified compliance corrections` → `without the design factors for attachment compliance`
- **Line 382.** `contact activation, wrapping, and compliance correction` → `contact activation, wrapping, and deflection under load`
- **Line 445.** `\begin{tabularx}{\textwidth}{@{}p{0.27\textwidth}p{0.16\textwidth}X@{}}` → `\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.27\textwidth}p{0.16\textwidth}>{\raggedright\arraybackslash}X@{}}`
- **Line 478.** `The corrections identified on the pinned knee` → `The design factors identified on the pinned knee`
- **Line 518.** `before and after correction` → `with and without the design factors`
- **Line 528.** `before and after correction, using the identified terms` → `with and without the identified design factors`
- **Line 553.** `the correction terms disabled` → `the design factors disabled`
- **Line 559.** `The correction terms` → `The design factors`
- **Line 559.** `enters the usable-length correction in` → `scales the loss of usable length in`
- **Line 569.** `before and after correction, with stiffness shared with the flexor` → `with and without the design factors. The fitted model uses the flexor stiffness values`
- **Line 602.** `The extensor corrections` → `The extensor design factors`
- **Line 616.** `The shared correction terms` → `The shared design factors`
- **Line 649.** `Reasons for the Correction Terms` → `Physical Basis of the Design Factors`
- **Line 651.** `The correction terms` → `The design factors`
- **Line 657.** `the zero-correction baseline` → `the baseline model without the design factors`
- **Line 661.** `Adopted correction terms` → `Adopted design factors`
- **Line 673.** `require the compliance corrections` → `require the design factors for compliance`
- **Line 677.** `the correction terms in` → `the design factors in`

### Chapter 5

[LaTeX source](../../Final/chapters/50-database.tex)

- **Line 4.** `different species, preparations, and tasks` → `different species, experimental conditions, and tasks`
- **Line 12.** `the species and preparation` → `the species and experimental conditions`
- **Line 36.** `without leg movement, a preparation called fictive locomotion` → `without leg movement, a pattern of neural activity called fictive locomotion`
- **Line 38.** `the preparation and response` → `the experimental conditions and response`
- **Line 60.** `compares the midbrain-evoked and acute spinal preparations` → `compares these results with those from cats whose spinal cords had been transected`
- **Line 111.** `why the preparation matters` → `why the conditions of the experiment matter`
- **Line 121.** `the source, preparation, stimulus, and measured response` → `the source, experimental conditions, stimulus, and measured response`

### Chapter 6

[LaTeX source](../../Final/chapters/60-sns.tex)

- **Line 172.** `For rhythm and sensory timing, the search scripts use tree-structured Parzen estimation to propose new parameter sets from earlier scores.` → `The Gait2392 stepping search uses Optuna \citep{akiba_optuna_2019}, and the staged rhythm and sensory-timing searches use Rustuna, its Rust implementation \citep{optuna_rustuna_2026}. Both use a tree-structured Parzen estimator, which favors settings typical of the best-scoring trials and uncommon in the rest \citep{bergstra_algorithms_2011}.`

### Chapter 7

[LaTeX source](../../Final/chapters/70-overall.tex)

- **Line 8.** `Figure~\ref{fig:design_workflow} traces actuator sizing and route selection. Figure~\ref{fig:discussion_neuromechanical_workflow} shows how to test a change in a sensory pathway or an actuator route against the same baseline.` → `Figure~\ref{fig:design_workflow} shows the mechanical design process for a BPA-actuated joint, from a human torque target to the choice of actuator size and route. Figure~\ref{fig:discussion_neuromechanical_workflow} shows the neuromechanical modeling workflow. A trial with the original controller and actuator route provides the baseline. Each repeat changes one sensory pathway or one actuator route while keeping the initial pose, stimulus, and support the same.`
- **Line 29.** `A Festo patent shows two concentric layers of fibers wound in opposite directions, with flexible rubber between them` → `A Festo patent describes two concentric fiber layers wound in opposite directions and identifies rubber as a material that can separate them`
- **Line 29.** `The model can represent the wall as a nearly incompressible hyperelastic material and the reinforcement as two embedded families of tension-only truss elements following the measured fiber angles.` → `In the proposed model, nearly incompressible hyperelastic rubber fills the space between the two helical aramid layers. Each layer is represented by tension-only truss elements embedded in the rubber at its measured winding angle.`
- **Line 45.** `Sectioning a test specimen would establish whether to include the separate fabric covering shown here and would provide the layer thicknesses, fiber angles, and interfaces for the finite-element model.` removed from the caption.
- **Line 65.** `its stimulus, preparation, and locomotor phase` → `the stimulus, experimental conditions, and locomotor phase`

### Chapter 8

[LaTeX source](../../Final/chapters/80-conclusion.tex)

- **Line 8.** `shared correction terms` → `shared design factors`

### Appendix C

[LaTeX source](../../Final/chapters/94-AppendixC.tex)

- **Line 79.** `The correction terms` → `The design factors`
- **Line 159.** `precedes the offset, compliance, and wrapping corrections` → `is the length before applying the design factors for offset, compliance, and wrapping`
- **Line 240.** `physical correction terms` → `physical design factors`
- **Line 243.** `physical correction-term bounds` → `bounds on the physical design factors`
- **Line 261.** `Adopted Correction Terms` → `Adopted Design Factors`
- **Line 263.** `adopted correction terms` → `adopted design factors`
- **Line 265.** `before and after applying the adopted corrections` → `with and without the adopted design factors`
- **Line 265.** `with zero length and wrapping corrections and rigid attachments` → `with zero length offset, no wrapping loss, and rigid attachments`
- **Line 265.** `the adopted correction terms` → `the adopted design factors`
- **Line 271.** `before and after correction` → `with and without the design factors`
- **Line 290.** `before and after correction` → `with and without the design factors`
- **Line 311.** `The adopted corrections` → `The adopted design factors`
- **Line 311.** `larger errors after correction` → `larger errors with the design factors`

### Appendix E

[LaTeX source](../../Final/chapters/96-AppendixE.tex)

- **Line 38.** `why the preparation belongs alongside the label` → `why the experimental conditions belong alongside the label`
- **Line 38.** `In the acute spinal preparation treated with clonidine` → `In cats given clonidine after acute spinal transection`
- **Line 38.** `seen in the midbrain-evoked preparation` → `seen with midbrain stimulation`
- **Line 38.** `A model of the spinal preparation` → `A model of the experiment after transection`

### Bibliography

Three records were appended to the bibliography. All earlier bytes are preserved. [Bibliography](../../Final/thesis.bib)

| Citation key | Source supporting Chapter 6 |
| --- | --- |
| `bergstra_algorithms_2011` | Bergstra, Bardenet, Bengio, and Kégl. “Algorithms for Hyper-Parameter Optimization” (2011). [Original paper](https://proceedings.neurips.cc/paper_files/paper/2011/hash/86e8f7ab32cfd12577bc2619bc635690-Abstract.html) |
| `akiba_optuna_2019` | Akiba, Sano, Yanase, Ohta, and Koyama. “Optuna: A Next-generation Hyperparameter Optimization Framework” (2019). [Authors’ paper](https://arxiv.org/abs/1907.10902) |
| `optuna_rustuna_2026` | Optuna Developers. Rustuna, the Rust implementation used by the staged trainer. [Project repository](https://github.com/optuna/rustuna) |

The salamander paper compares its model with experiments on the spinal cord. The revised Appendix E wording distinguishes two conditions in McCrea’s cat experiments. [Ijspeert et al.](https://infoscience.epfl.ch/record/101375/files/science2007_preprint.pdf) · [McCrea et al.](https://pubmed.ncbi.nlm.nih.gov/8558481/)

## Editorial feedback and verification

The preceding local judgment on the Chapter 4 Results opening overlooked “corrections.” Its original score and wording are retained with the new judgment. Ben’s choice of “design factors” and his direction to name vivisection are recorded as author feedback. The Chapter 2 score distinguishes his chosen final sentence from the revised paragraph. The advisor’s tracked edits and comments retain their original scores. [Passages and judgments](PASSAGE_REVIEW.json) · [Author examples](../advisor_feedback_history/EDITORIAL_EXAMPLES.json)

The [Overleaf dissertation](https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb) compiled to 327 pages with zero errors, zero warnings, and 14 informational entries. The [local PDF](../../Bolen_Dissertation.pdf) has 327 pages. All citations and references resolved. The local log contains 17 unrelated underfull horizontal boxes and zero overfull boxes. The [source ZIP](../../Bolen_Dissertation.zip) passed its integrity check. Eleven of its 165 members changed from the preceding release. The other 154 match byte for byte. The eleven are the nine chapter and appendix sources named above, the bibliography, and the bracket diagram. The numerical results, equation labels, and Figure 4.9 data are preserved.

The historical-feedback check matched all 7 saved passage judgments against the revised text. It compared 1,151 passages, captions, and table entries with 444 documented advisor comments and tracked edits. The existing Results heading in the published Chapter 3 opens directly into a subsection.

Full PDF SHA-256 `665a165ce4962494a39e2f657ea5bd21b3920f246557414d4205f48986476bd5`. Source ZIP SHA-256 `6239ce18e1d43e38964bb7050e631ef9e881c0e2cde68ab5e77b5905ee7d8f04`.
