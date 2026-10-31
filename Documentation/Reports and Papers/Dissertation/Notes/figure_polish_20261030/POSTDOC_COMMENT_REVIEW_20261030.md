# Leo Micklem dissertation review - proposed responses

Review window: 30 October 2026, 6:12-7:23 pm Pacific. All 43 open threads and the additional 7:20 reply were read live: 44 messages across the abstract, introduction, and SADb chapter. This is a read-only review, not approval to implement or resolve comments.

Ben's explicit decisions and the advisor's requirements take precedence. Recommendations are evaluated on scientific correctness, reader clarity, duplication, and revision scope. The current abstract/introduction match local ProofFinal; live Chapter5 has its presentation credits in the opener and lacks the local OpenAlex gloss, so its line numbers differ by two. Locations below name both where necessary.

## Recommended priorities

1. Clarify the abstract's purpose and causal sequence while retaining one page, one paragraph, BPA-before-torque order, no SADb record count, and qualified virtual-harness outcomes.
2. Separate motivation from achieved contributions in the Introduction, using the existing four-aim organization. Retain all three broader research motivations and the master's-project chronology.
3. Add one genuine SADb worked example and make the Discussion explain what the database supplies versus what the researcher concludes. Clarify the candidate-pathway-to-controller-test handoff.
4. Make small terminology, ownership, and sentence-length repairs. Treat typography, a new curation flowchart, and a second abstract paragraph as optional.
5. Keep the AI-triage example pending unless a genuine documented example is recovered; no new formal assessment or invented evidence.

No manuscript, figure, bibliography, downloaded ZIP, Zotero record, review comment, or tracked change was changed. No live upload, compile, commit, or push was performed.

## A01 - 30 October 2026, 6:12 pm Pacific

**Abstract; 02-abstract.tex:7 (ProofFinal).**

> There's no real set up for this sentence. You mention biomimetic robotics in general and then jump to human robots

**Agree | High.** Name the musculoskeletal humanoid and its actuation/control problem in the opening, then state the specific gaps. The current general-biomimetic opener does not prepare the shift to a human-proportioned BPA robot.

## A02 - 30 October 2026, 6:14 pm Pacific

**Abstract; 02-abstract.tex:7 (ProofFinal).**

> It would be nice if there was a logical flow to these points or a story you're telling with the project rather than just a list

**Agree | High.** Replace the needs-list followed by a results-list with a causal sequence: force characterization enables torque prediction; joint losses motivate route redesign; sensory evidence and modeling tools enable controller tests. Preserve the four aims and their order.

## A03 - 30 October 2026, 6:15 pm Pacific

**Abstract; 02-abstract.tex:7 (ProofFinal).**

> This sentence would be easier to understand if it was in the first person. 'I produced a x model by...'

**Agree | Medium.** Use first person for the contribution: "I characterized Festo BPAs across resting lengths, contractions, and pressures and developed..." Keep BPA characterization before joint torque, as Ben requested.

## A04 - 30 October 2026, 6:19 pm Pacific

**Abstract; 02-abstract.tex:7 (ProofFinal).**

> I think each of these are effectively chapter conclusion summaries but they need some more set up. A thesis abstract can be much longer than a paper abstract. I would consider splitting this into two paragraphs and setting up each takeaway a little more from the needs prior. I find the two lists a little disconnected.

**Modify | High.** Accept the need to connect gaps and takeaways. Keep Ben's one-page, one-paragraph choice for now; a second paragraph is optional, not required to fix the argument. Do not expand merely because a thesis can be longer.

## A05 - 30 October 2026, 6:16 pm Pacific

**Abstract; 02-abstract.tex:7 (ProofFinal).**

> I think you need something like this at the start of the abstract to help the reader understand where you're going with it.

**Agree | High.** Move the purpose now stated in the closing sentence into the opening: experimentally informed robot design and controllable neuromechanical experiments. Keep prosthetic/orthotic applications as potential downstream use, not an evaluated outcome.

## I01 - 30 October 2026, 6:21 pm Pacific

**1.1 Motivation; 10-introduction.tex:6 (ProofFinal).**

> For me this is missing a high level intro of why we should be interested in biomimetic humanoid robots

**Agree | High.** Lead with the unstructured-terrain problem and the value of controlled biological experiments; place the deliverables after the motivation. Do not claim energy savings or terrain autonomy were demonstrated.

## I02 - 30 October 2026, 6:24 pm Pacific

**1.1 Motivation; 10-introduction.tex:8 (ProofFinal).**

> This paragraph is less a goal and more a motivation which I think should proceed the previous paragraph

**Modify | High.** Put both broader motivations before the contribution list. Keep Ben's priorities: better robots first, biological understanding second, and prosthetic/orthotic potential third.

## I03 - 30 October 2026, 6:25 pm Pacific

**1.1 Motivation; 10-introduction.tex:8 (ProofFinal).**

> Not sure this is a goal specifically and more motivation for what this work could also be beneficial for. You don't do anything related to this in your work

**Modify | Medium.** Retain prosthetic/orthotic design as a longer-term application: "These methods may also inform future prosthetic and orthotic devices." Do not describe it as a tested dissertation deliverable.

## I04 - 30 October 2026, 6:26 pm Pacific

**1.1 Motivation; 10-introduction.tex:10 (ProofFinal).**

> It's not clear to me what class of robot you mean here

**Agree | High.** Replace "this class of robot" with "musculoskeletal humanoid robots actuated by artificial muscles." General biomimetic robots need not use muscles.

## I05 - 30 October 2026, 6:27 pm Pacific

**1.1 Motivation; 10-introduction.tex:10 (ProofFinal).**

> Good

**Retain | None.** Retain the explanation of BPA-versus-skeletal-muscle differences. No corrective edit is needed.

## I06 - 30 October 2026, 6:37 pm Pacific

**1.1 Motivation; 10-introduction.tex:10 (ProofFinal).**

> This doesn't belong in the motivation. Abstract or contributions

**Agree | High.** Move "I therefore characterized...developed a force model" to the contribution/objective summary. Let this motivation paragraph end with the unmet prediction need.

## I07 - 30 October 2026, 6:41 pm Pacific

**1.1 Motivation; 10-introduction.tex:12 (ProofFinal).**

> I think you can say this without the previous sentence. You're just trying to motivate the work here.

**Modify | Medium.** An optional shorter transition is reasonable. First move the preceding completed-work sentence; keep a clear distinction between actuator force and transmitted joint torque. This is not a reason to delete the transmission-loss rationale.

## I08 - 30 October 2026, 6:42 pm Pacific

**1.1 Motivation; 10-introduction.tex:12 (ProofFinal).**

> Good. Maybe say 'a robot knee' so it's a little more generalised

**Modify | Low.** Use "In a BPA-actuated robot knee..." for the general mechanism, then "the tested knee" for the study-specific evidence. Do not generalize the measurements to all robot knees.

## I09 - 30 October 2026, 6:48 pm Pacific

**1.1 Motivation; 10-introduction.tex:12 (ProofFinal).**

> Muscle routes? Tendon routes?

**Agree | High.** Define "actuator-tendon routes" as the combined BPA and artificial-tendon paths through attachment and wrapping points, then use that term consistently.

## I10 - 30 October 2026, 6:43 pm Pacific

**1.1 Motivation; 10-introduction.tex:12 (ProofFinal).**

> You're citing a dissertation that isn't 'this dissertation' being your one

**Modify | Medium.** Clarify collaborative work done during Ben's doctorate. The bibliography records Morrow, Bolen, and Hunt (2020), a coauthored paper, not a different dissertation. Keep the citation and distinguish it from the 2019 master's work.

## I11 - 30 October 2026, 6:44 pm Pacific

**1.1 Motivation; 10-introduction.tex:14 (ProofFinal).**

> I would make it more clear that this is just the first problem

**Agree | High.** Explicitly label the first neural-modeling gap as access to physiological evidence and the second as controller transfer to different mechanics. Each motivates a separate aim.

## I12 - 30 October 2026, 6:45 pm Pacific

**1.1 Motivation; 10-introduction.tex:14 (ProofFinal).**

> Again, not part of motivation

**Agree | High.** Move "I addressed that evidence-access problem with the Sensory Afferent Database" to objective/contribution 3; keep the preceding evidence-access gap.

## I13 - 30 October 2026, 6:46 pm Pacific

**1.1 Motivation; 10-introduction.tex:14 (ProofFinal).**

> Good

**Retain | None.** Retain the explanation that circuit behavior depends on mechanics, sensory input, support, and numerical implementation.

## I14 - 30 October 2026, 6:47 pm Pacific

**1.1 Motivation; 10-introduction.tex:14 (ProofFinal).**

> I think this doesn't belong in the motivation either. Abstract or have a contributions section.

**Agree | High.** Move "I developed coupled..." and the representative-test outcomes into the fourth contribution. Keep the motivation focused on the need to test the coupled model.

## I15 - 30 October 2026, 6:48 pm Pacific

**1.1 Motivation; 10-introduction.tex:16 (ProofFinal).**

> This is good but you could strengthen it by moving those contributions that are dotted in the motivation to here.

**Agree | High.** Gather the relocated contribution sentences in the existing four-aim summary or Research Objectives section. Avoid adding another nearly identical contributions list.

## I16 - 30 October 2026, 6:50 pm Pacific

**1.2 The Biomimetic Humanoid Robot Project; 10-introduction.tex:20 (ProofFinal).**

> This is missing the 'why' of the work for me

**Modify | Medium.** Add one purpose-led link: the platform permits controlled tests of actuation and neural-control hypotheses. Do not repeat the full Motivation section.

## I17 - 30 October 2026, 6:49 pm Pacific

**1.2 The Biomimetic Humanoid Robot Project; 10-introduction.tex:20 (ProofFinal).**

> Not sure where the trunk begins and ends.

**Agree | Medium.** Replace "from the trunk down" with an anatomically explicit description, matched to the illustrated assembly: lumbar structure, pelvis, and legs. Do not imply a complete upper torso.

## I18 - 30 October 2026, 6:51 pm Pacific

**1.2 The Biomimetic Humanoid Robot Project; 10-introduction.tex:22 (ProofFinal).**

> I think this is the same as what you mean by routes previously but that could be clarified

**Agree | Medium.** Use the same actuator-tendon-route terminology for earlier attachment/wrapping placement. Preserve the master's-to-doctoral chronology and distinguish rigid placement searches from later design-factor redesign.

## I19 - 30 October 2026, 6:52 pm Pacific

**1.2 The Biomimetic Humanoid Robot Project; 10-introduction.tex:22 (ProofFinal).**

> This is one very long sentence and you lost me in what you're trying to say

**Agree | High.** Split the long platform-history sentence into separate knee-linkage and Muscle Mutt/SNS statements, keeping the source beside each claim. No scientific content needs to be lost.

## I20 - 30 October 2026, 6:54 pm Pacific

**1.3 Research Objectives; 10-introduction.tex:24 (ProofFinal).**

> It should be clear to me from the motivation paragraphs why you have chosen each of these objectives

**Agree | High.** Make the four motivating needs explicit before their matching objectives: resting-length force prediction, transmission-aware torque prediction, accessible sensory evidence, and controller-transfer tools.

## I21 - 30 October 2026, 6:52 pm Pacific

**1.3 Research Objectives; 10-introduction.tex:26 (ProofFinal).**

> Good

**Retain | None.** Retain the objective introduction and four-aim organization; absorb moved contributions without creating another full list.

## I22 - 30 October 2026, 6:55 pm Pacific

**1.3 Research Objectives; 10-introduction.tex:30 (ProofFinal).**

> I'm not sure this is described earlier?

**Modify | Medium.** Wrapping-induced loss of usable length already appears in line12. Make the terminology and connection to objective2 easier to follow; this is signposting, not missing science or a need for new equations.

## I23 - 30 October 2026, 7:02 pm Pacific

**1.3 Research Objectives; 10-introduction.tex:32 (ProofFinal).**

> This only really gets one motivation sentence - I had to go back and check why you were doing this

**Agree | High.** Give controller transfer a short separate rationale: hand tuning after transfer to different mechanics is inefficient and does not establish walking. Connect graphical editing and automated evaluation to Gait2392 and BPA robot models.

## I24 - 30 October 2026, 7:03 pm Pacific

**1.4 Dissertation Organization; 10-introduction.tex:35 (ProofFinal).**

> This and the figure are good

**Retain | None.** Retain the dissertation roadmap and figure. Update pointers only if later approved changes require it.

## S01 - 30 October 2026, 7:04 pm Pacific

**Chapter 5 opening; 50-database.tex:4 (ProofFinal).**

> I think you need to make it clear that this is your work in the first sentence

**Agree | High.** Start "I developed the Sensory Afferent Database (SADb), a curated literature resource..." Keep the AARL/C3NS context and both 2023/2024 presentation credits.

## S02 - 30 October 2026, 7:05 pm Pacific

**Chapter 5 opening; 50-database.tex:8 (ProofFinal); live line6.**

> Is this someone who is building a model?

**Agree | Medium.** Use "A researcher building a neural or neuromechanical model" at first mention instead of the unexplained "A modeler."

## S03 - 30 October 2026, 7:07 pm Pacific

**Chapter 5 opening; 50-database.tex:19 (ProofFinal); live line17.**

> I think it would be helpful to the reader if you put 'clickable' items on the site in italics

**Modify | Low.** Optional: italicize actual tab/control labels consistently, such as Table and Search (online). Keep "tabs" and the existing callouts. Do not italicize every clickable concept or redesign the app.

## S04 - 30 October 2026, 7:08 pm Pacific

**5.1 Curation Methods; 50-database.tex:23 (ProofFinal); live line21.**

> Was?

**Modify | Low.** The live anchor is "This was a curated research collection." Prefer "SADb is a curated research collection..." or "I curated..."; retain past tense for completed import and screening actions.

## S05 - 30 October 2026, 7:08 pm Pacific

**5.1 Curation Methods; 50-database.tex:23 (ProofFinal); live line21.**

> What keywords?

**Modify | Medium.** Clarify representative search terms, but distinguish historical manual searching from the current monthly watch. Its config verifies terms such as Ia/Ib afferent, muscle spindle, proprioceptive feedback, and central pattern generator; it does not prove the original manual queries.

## S06 - 30 October 2026, 7:11 pm Pacific

**5.1 Curation Methods; 50-database.tex:25 (ProofFinal); live line23.**

> I feel some kind of flow diagram for the order of operations could be helpful

**Modify; optional figure | Medium.** Make the sequence explicit: collect, deduplicate, screen, classify/annotate, and transfer to Explorer; monthly searches feed manual curation. A small flow diagram is optional and should replace repetitive prose, not add a second account.

## S07 - 30 October 2026, 7:16 pm Pacific

**5.2 Results: Application and Corpus Coverage; 50-database.tex:29 (ProofFinal); live line27.**

> I think it would be really helpful to pick an example that you can walk through of how you'd use this

**Agree | High.** Add one short real use case using the existing Grillner screenshot: select a paper, inspect notes, follow a cited/citing study, and inspect its evidence. Verify any named linked record; do not assume every paper has a pathway link. Replace repeated menu descriptions.

## S08 - 30 October 2026, 7:12 pm Pacific

**5.2 Results: Application and Corpus Coverage; 50-database.tex:33 (ProofFinal); live line31.**

> I'm assuming this is mentioned in the background section

**Retain local explanation | Low.** The anchor is Louvain. It is not explained in Background, but the present sentence already states its function: grouping citation-connected papers. Keep that brief local explanation; a new background subsection is unnecessary.

## S09 - 30 October 2026, 7:16 pm Pacific

**5.2.1 AI-Assisted Literature Triage; 50-database.tex:67 (ProofFinal); live line65.**

> This feels like more of an abstract for the subsection rather than the section itself. Where is the substance of what was done, examples of things not working etc?

**Modify; example pending | Medium.** The live anchor is the AI-triage heading, not corpus coverage. Acknowledge the lack of a concrete example. Consider one genuine example from the 2024 poster if available; otherwise leave this pending. Preserve Ben's subjective judgment against his own speed/quality; no invented protocol, timings, failures, or new rating exercise.

## S10 - 30 October 2026, 7:18 pm Pacific

**5.3.1 Uses of the Database; 50-database.tex:78 (ProofFinal); live line76.**

> Does the database directly tell you these differences or connect the papers? I think it's good to try to use the example but I'm not sure of the point.

**Agree | High.** Lead with what SADb supplies: curated records, notes, and source connections. The researcher compares study conditions and assumptions; do not imply the app automatically infers the biological distinctions.

## S11 - 30 October 2026, 7:19 pm Pacific

**5.3.1 Uses of the Database; 50-database.tex:80 (ProofFinal); live line78.**

> I'm not sure how this paragraph relates to the database?

**Agree | High.** Start the Geyer/Deng paragraph with its database-use point, then retain only the study detail needed to show how comparing records informs controller choice.

## S11r - 30 October 2026, 7:20 pm Pacific

**5.3.1 Uses of the Database; 50-database.tex:80 (ProofFinal); live line78.**

> These paragraphs should probably start with the point you're trying to make and then show evidence

**Agree | High.** Use point, source example, and modeling implication in that order. Preserve informative primary-study examples Ben requested, without claiming an unverified link between papers.

## S12 - 30 October 2026, 7:21 pm Pacific

**5.3.1 Uses of the Database; 50-database.tex:80 (ProofFinal); live line78.**

> Muscle pathway? Neural pathway?

**Agree | Medium.** Say "sensory-feedback pathway" here, to distinguish neural connections from the mechanical actuator-tendon routes in Chapter4.

## S13 - 30 October 2026, 7:22 pm Pacific

**5.3.1 Uses of the Database; 50-database.tex:82 (ProofFinal); live line80.**

> This should probably be in the motivation section? I'm not sure what relevance it has to the database

**Agree selectively | Medium.** The Ijspeert salamander example already supplies this broader motivation in Introduction line8. Keep that fuller example there; shorten its repeat in Chapter5 to a database-specific connection, or remove the duplication.

## S14 - 30 October 2026, 7:23 pm Pacific

**5.3.1 Uses of the Database; 50-database.tex:86 (ProofFinal); live line84.**

> Is the take home point supposed to be that you can use the database to select the optimal pathway to use for your robot? and then you go do the model development stuff?

**Modify | High.** Make the handoff explicit: SADb helps identify candidate sensory-feedback pathways and supporting evidence; Chapter6 implements and tests their encoding, connections, and gains. Do not claim the database selects an optimal pathway.

