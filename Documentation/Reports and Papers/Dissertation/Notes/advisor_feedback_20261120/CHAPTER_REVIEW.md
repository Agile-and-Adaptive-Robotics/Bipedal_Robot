# November 20 advisor review of Chapters 4, 5, 7, and 8 and Appendices C and E

The force-comparison finding for Chapter 7 below records the earlier draft. Chapter 7 now uses the matched foot-clearance experiment in [Muscle forces and pneumatic actuation](ACTUATOR_REVIEW.md).

The full Overleaf editor text matched the working copies of Chapter 4, Chapter 7, Chapter 8, Appendix C, and Appendix E before this pass. Chapter 5 differed in five tracked changes by Alexander Hunt. I used the complete editor text as the baseline, so all five appear in the revision, including two grammatical corrections to the inserted wording. Ben’s November 20 Chapter 4 deletion of “Pang, Scharzenberger,” was already present in the live and local copies. The author contributions were preserved. [Live review transcript](<D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_feedback_20261120/LIVE_REVIEW.txt>) · [Selections recovered from the editor](<D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_feedback_20261120/LIVE_ANCHORS.json>)

## Each Chapter 5 review item

Times are Pacific. The records identify Hunt as `ajh26` and Ben as `bbolen83`. The two 8:36 edits have distinct positions in the source sentence.

| Time | Source and attribution | Disposition in `50-database.tex` |
| --- | --- | --- |
| 8:36 | Hunt replaced “The” with “A” before the Neuroscience 2023 abstract and poster. | Retained at line 6. |
| 8:36 | Hunt made the corresponding replacement before Neuroscience 2024. | Retained at line 6. |
| 8:40 | Hunt inserted “AI” before “tools’ summaries”. | Retained at line 16 to identify which outputs Ben evaluated against his own reading. |
| 8:41 | Hunt asked why the McCrea study appeared. Bolen replied that Leo had requested a worked example. | The closing sentence of Section 5.1 explains that the example shows how to preserve the different experimental conditions behind two responses to one input. Section 5.1.1 then identifies McCrea and colleagues, the stimulus, and the measured response. |
| 8:44 | Hunt selected “To repeat the lookup,” and asked what was being repeated. Ben suggested “to follow the example” or “to launch the app” at 12:15. | Line 40 uses Ben’s “To follow this example” and gives the public app link, the pathway filter, and the record to select. |
| 8:46 | Hunt selected the second introduction of McCrea in Results and asked why the reader returned to this paper. Ben replied at 10:32 that Leo wanted an example of curation. | Line 50 links back to the Methods example and introduces a different question for the controller designer. What sign should the input have at rest and during extension? Figure 5.2 connects the record to the paper and Figure 5.3 shows the two net responses. |
| 8:46 | Hunt inserted “who wishes to implement this pathways”. | Retained his description of the modeler’s purpose and corrected agreement to “this pathway” at line 60. |
| 8:48 | Hunt added “a sensing organ in insects” after plural “campaniform sensilla”. | Line 80 defines them as sensory organs that detect strain in the insect cuticle. This supplies the requested definition and preserves plural agreement. [Primary recordings of load responses in cockroach legs](https://pubmed.ncbi.nlm.nih.gov/10082909/) |
| 8:49 | Hunt questioned “including Grillner.” The passage changed from the McCrea worked example to the separate paper visible in the Bubble map screenshot without explaining why. | Chapter 5 identifies Grillner’s review as the selected paper in Figure 5.5. Following its citation to the Forssberg experiment gives the reader a reason to open a source record. Appendix E gives the map counts and how the app counts unique papers. [Grillner review](https://doi.org/10.1126/science.3975635) |

The substantive request behind each thread is addressed in the local source. Clearing the four Overleaf comment threads and accepting individual tracked changes are separate editor operations. This pass made local manuscript edits only.

## Other chapters

| Passage | Decision and evidence |
| --- | --- |
| Chapter 4, §4.1.4, line 344 | Replaced a general statement about “checking” candidates with the exact comparison the reader can find in Tables 4.6 and 4.7. Each reports fitted and baseline errors on the training tests and the configurations held out of that fit. |
| Chapter 4, §4.3.2, line 657 | Replaced the opener about which procedure the figures used with its result. Prediction error decreased on the pinned-knee tests held out of the fit. The original training and validation percentages and all table entries remain. The existing §4.2 lead already states the physical result and §4.2.4 already explains the single-flexor shortfall. |
| Chapter 7, §7.2.4, line 79 | Identified the force-comparison observations that motivate the proposed future experiment. At the same pose and activation, MuJoCo produced a much larger rectus femoris force at one point than OpenSim. The MuJoCo mean was larger for the biceps femoris short head. Either muscle’s estimated force bears on the torque attributed to that muscle. Replaced “force law” with “BPA force model”, consistent with Ben’s 10:55 reply. Appendix D.1.2 supplies the quantitative comparison. |
| Chapter 8, line 12 | Replaced a list of checks with their meaning. The arm example demonstrates corrective feedback against a load. Under the same drive and support, the W2L model retained alternating heel strikes after the Ib-only intervention and lost that timing after the combined contact, pooled-force, and associated pathway intervention. Section D.4.2 gives heel-strike phases of 0.498 for the intact circuit, 0.494 after Ib-only deletion, and 0.223 after the combined deletion. The conclusion gives the implication without repeating those values. |
| Appendix C | Retained. Its tables document the design factors and the training, validation, and unused tests. Its geometric checks are conditions for accepting or releasing contact, with explicit consequences. These details do a reader’s reproducibility work. |
| Appendix E, §E.1, lines 8 and 36 | Explained what a researcher gains from searching the complete app before reporting the inventories. Distinguished links between publications from manually assigned pathway tags before presenting the citation totals. All counts, category definitions, and table values remain. |

Compared with the six complete live source files, all number strings and their multiplicities, table data rows, cited works, graphic includes, LaTeX environments, and pre-existing labels match. The only new label connects the curation example in Methods to its use in Results. `git diff --check` passed for all assigned files. The existing historical checker examined 507 units against 444 previously indexed advisor events. The fresh November 20 findings above come from Hunt’s actual review cards and highlighted text. No Chapter 5 comment depends on an unidentified selection.

The local changes are in `30-torque.tex`, `50-database.tex`, `70-overall.tex`, `80-conclusion.tex`, and `96-AppendixE.tex`. The unchanged `94-AppendixC.tex` matches the live capture. The parent session coordinates the full dissertation build and live synchronization.
