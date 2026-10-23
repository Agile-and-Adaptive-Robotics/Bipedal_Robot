# chatgpt_report

Advisor comment inventory and SNS revision proposal

Reviewed 22 October 2026, after your comment cleanup. Source: the live Overleaf project Bolen-Dissertation. Advisor account: ajh26.

## Scope and reading conventions

The current overview contains **36 open advisor comment threads (37 advisor messages) and 8 pending advisor tracked edits**, across four files. SNS has 11 comment threads and 4 tracked edits; the database chapter has 24 threads and 2 tracked edits; the introduction has 1 thread; the background has 2 tracked edits. Resolved/deleted comments and your own review entries are excluded. Your replies are not presented as advisor comments.

Line numbers are LaTeX source lines, not PDF line numbers. A paragraph often occupies one source line. They describe the source at review time and will shift after edits. Highlighted text is transcribed from the editor's comment markings, with adjacent fragments joined across spellchecking markup. Some anchors are only one word even when the comment discusses a paragraph. LaTeX commands are preserved. Tracked deletions are reported at their current insertion/deletion point; deleted text has no remaining highlight.

No chapter prose, suggestions, or comment status was changed for this report. The inventory comes first; the background cross-check and proposed changes follow it.

## 1. SNS chapter — chapters/60-sns.tex

### SNS-01 · line 6 · comment

**Highlighted text:**

```latex
Two definitions cover most of the neuroscience vocabulary used below. First, the locomotor controller is a two-level central pattern generator: a rhythm generator sets the timing of the step cycle, and a pattern-formation layer shapes that timing into phase-appropriate activation of antagonist muscle groups \citep{rybak_modelling_2006}. Each leg's rhythm generator is a pair of mutually inhibiting antagonist populations, one extensor and one flexor, so that only one is strongly active at a time. Second, in the non-spiking (graded) formulation used for most of this work, each neuron is electrically a nonlinear RC circuit: synaptic current charges a membrane capacitance through a leak conductance, and the resulting membrane voltage, tens of millivolts, is both the neuron's state and its output \citep{nourse_sns_2023}. Graded chemical synapses couple these units: a depolarized presynaptic cell opens a conductance in the postsynaptic membrane whose sign makes the connection excitatory or inhibitory. The sensory channels that close the loop are the muscle afferents, Ia signaling muscle length and velocity, Ib signaling muscle force, and group II signaling length, together with heel and toe contact sensors.
```

**Advisor — 21 October, 11:48 am:**

> So many details and I don't understand what the goal is. What is the gap in knowledge? What did you do to fill it?

### SNS-02 · line 8 · comment

**Highlighted text:**

```latex
the reflex-based and learned alternatives \citep{geyer_musclereflex_2010, song_deep_2021} are treated there.
```

**Advisor — 21 October, 12:45 pm:**

> huh?

### SNS-03 · line 10 · comment

**Highlighted text:**

```latex
tooling
```

**Advisor — 21 October, 12:47 pm:**

> tooling tends to be a manufacturing term. Software Tool Development is more like what you are looking for here.

### SNS-04 · line 10 · tracked edit

**Advisor's tracked edit:**

```text
Added: virtual
```

### SNS-05 · line 10 · tracked edit

**Advisor's tracked edit:**

```text
Changed: that
to this
```

### SNS-06 · line 10 · comment

**Highlighted text:**

```latex
built in any of several simulation engines
```

**Advisor — 21 October, 12:49 pm:**

> why would we need several?

### SNS-07 · line 10 · comment

**Highlighted text:**

```latex
answer
```

**Advisor — 21 October, 12:53 pm:**

> This is all super awkward. You need to identify the gap in knowledge and how you fill that gap in knowledge. Much of what you wrote here should be saved for discussion, where you talk about the implications of the work.

The anchor is the word “answer” at the end of the long introductory paragraph on the virtual laboratory; the comment explicitly addresses that broader framing.

### SNS-08 · line 14 · tracked edit

**Advisor's tracked edit:**

```text
Deleted: Preliminary Results:
```

### SNS-09 · line 20 · comment

**Highlighted text:**

```latex
verification-first approach
```

**Advisor — 21 October, 12:59 pm:**

> I don't understand the connection.

### SNS-10 · line 24 · comment

**Highlighted text:**

```latex
build order
```

**Advisor — 21 October, 1:02 pm:**

> huh?

### SNS-11 · line 26 · comment

**Highlighted text:**

```latex
\begin{figure}[htbp]
```

**Advisor — 21 October, 1:00 pm:**

> Overlapping text in figure

Figure anchor: `pipeline.pdf`, caption “Simulation-to-robot pipeline,” label `fig:pipeline` (source lines 26–31).

### SNS-12 · line 35 · tracked edit

**Advisor's tracked edit:**

```text
Deleted: The laboratory stands on tools this program did not build.
```

### SNS-13 · line 37 · comment

**Highlighted text:**

```latex
The SADb of Chapter~\ref{ch:database} is the data resource among the existing tools \citep{bolen_sensory_2023}. Its segments are the five-table Airtable base (Papers, Feedback, Models, Review Papers, and Recorders), the laboratory's Zotero libraries that feed it, citation mining in the Research Rabbit style that expanded it, and the laboratory's own packaged explorer application that distributes it. The SADb framed this chapter's benchmarks and validation measurements: its pathway records fix what each afferent class encodes and what each connection rule does, so the sensory channels, the reflex topologies, and the acceptance phenomena against which the controllers are validated (the deafferentation and stimulation results cited in Section~\ref{sec:walker_variant_results}) are looked up rather than re-derived from the primary literature.
```

**Advisor — 21 October, 1:06 pm:**

> This whole paragraph seems unnecessary.

### SNS-14 · line 39 · comment

**Highlighted text:**

```latex
\subsection{New Tools}\label{sec:sns_new_tools}
```

**Advisor — 21 October, 1:08 pm:**

> Wait, that is all you are saying about the other tools? It is hard for the reader to understand the limitations and gaps without a deeper description of them.

### SNS-15 · line 41 · comment

**Highlighted text:**

```latex
Simulink blocks
```

**Advisor — 21 October, 1:11 pm:**

> Did you compare and discuss against the libraries myself and Fletcher made here?


## 2. Database chapter — chapters/50-database.tex

### DB-01 · line 6 · comment

**Highlighted text:**

```latex
and design needs parameters, not only principles
```

**Advisor — 21 October, 11:13 am:**

> huh?

### DB-02 · line 6 · comment

**Highlighted text:**

```latex
For each sensory channel, the designer has to know what the receptor encodes, over what ranges of length, velocity, and force it operates and with what receptor dynamics,
```

**Advisor — 21 October, 11:13 am:**

> not sure what this has to do with the database

### DB-03 · line 6 · comment

**Highlighted text:**

```latex
whether the target is ipsilateral or contralateral
```

**Advisor — 21 October, 11:14 am:**

> there are more options

### DB-04 · line 6 · comment

**Highlighted text:**

```latex
Every element of that specification was established experimentally
```

**Advisor — 21 October, 11:14 am:**

> There are still lots of unknowns

### DB-05 · line 6 · comment

**Highlighted text:**

```latex
literature
```

**Advisor — 21 October, 11:16 am:**

> After this sentence, define the gap in knowledge and what you did to fill it.

### DB-06 · line 6 · tracked edit

**Advisor's tracked edit:**

```text
Deleted:  rather than collected anywhere
```

### DB-07 · line 6 · comment

**Highlighted text:**

```latex
A controller designer is not the only reader stranded by that scatter. An animal researcher planning the next experiment wants to know which preparations established a pathway property, and how the property compares across species, before committing a preparation to it. A roboticist wants to know which receptor encoding a sensor should imitate. Those readers ask different questions of the same literature, and the literature answers none of them without being gathered, structured, and searchable first.
```

**Advisor — 21 October, 11:15 am:**

> This should probably be the beginning of the paragraph

### DB-08 · line 8 · comment

**Highlighted text:**

```latex
Reviews
```

**Advisor — 21 October, 11:19 am:**

> The tone of this paragraph is as an opening paragraph, and yet it is the 2nd. Probably doesn't need to be here. Discussion of the reviews as a starting point probably goes in your 'methods' section of this chapter

### DB-09 · line 10 · comment

**Highlighted text:**

```latex
Neuroscience 2023 \citep{bolen_sensory_2023}
```

**Advisor — 21 October, 11:20 am:**

> I thought it was presented 2x, once without the AI help, and the other with.

### DB-10 · line 10 · comment

**Highlighted text:**

```latex
distilled mechanism notes
```

**Advisor — 21 October, 11:20 am:**

> huh?

### DB-11 · line 10 · comment

**Highlighted text:**

```latex
first-class record
```

**Advisor — 21 October, 11:21 am:**

> huh?

### DB-12 · line 10 · comment

**Highlighted text:**

```latex
pathway records
```

**Advisor — 21 October, 11:21 am:**

> huh?

### DB-13 · line 10 · comment

**Highlighted text:**

```latex
The database lives in a five-table Airtable base (Papers, Feedback, Models, Review Papers, and Recorders)
```

**Advisor — 21 October, 11:25 am:**

> I thought you migrated it to a different system

### DB-14 · line 10 · comment

**Highlighted text:**

```latex
canonical
```

**Advisor — 21 October, 11:25 am:**

> database?

### DB-15 · line 12 · comment

**Highlighted text:**

```latex
This chapter and the controller chapter are therefore built as one program.
```

**Advisor — 21 October, 11:26 am:**

> huh?

### DB-16 · line 12 · comment

**Highlighted text:**

```latex
The database informs the sensory channels of Chapter~\ref{ch:sns}: the afferent classes the controller monitors, the signal encoding each channel carries, and the connection rules its circuit implements all trace to records held here, and Section~\ref{sec:sadb_use} documents that mapping.
```

**Advisor — 21 October, 11:28 am:**

> This is a discussion point for the end of the chapter.

### DB-17 · line 13 · tracked edit

**Advisor's tracked edit:**

```text
Added:  Methods
```

The insertion point is on line 13; the inserted heading text “Methods” appears on line 14 in the current copied source.

### DB-18 · line 14 · comment

**Highlighted text:**

```latex
Methods
```

**Advisor — 21 October, 11:33 am:**

> You don't have to follow this exactly, but you should generally put this in the Chochrane form. https://www.cochrane.org/sites/default/files/uploads/PDFs/guide_to_the_contents_of_a_cochrane_methodology_protocol_and_review.pdf

### DB-19 · line 16 · comment

**Highlighted text:**

```latex
audit
```

**Advisor — 21 October, 11:29 am:**

> ?

### DB-20 · line 16 · comment

**Highlighted text:**

```latex
The author's personal Zotero collection held 538 items, the AARL group library's Biology subtree held 911 items across 27 collections, and the Airtable base held 99 papers. Matching record by record on DOI and normalized title showed that the personal collection was essentially a subset of the group library (all 483 personal DOIs present) and that the Airtable base was a small curated subset rather than a mirror: 95 of its 99 papers resolved to records in the group library's Biology subtree, while 443 personal items had never been entered at all. The gap was imported in verified batches, with every written record re-read and checked after writing, bringing the Papers table to 532 records. The corpus then kept growing while curation caught up: citation mining in the style of Research Rabbit contributed a Crossref-verified digest of recent papers, newly saved references flowed in from the laboratory's Zotero libraries, and a bulk pass imported 372 papers auto-classified from the laboratory's staged collection of PDFs.
```

**Advisor — 21 October, 11:37 am:**

> This is super confusing. We don't need all these details. The airtable base is the database...it can't 'start' with anything if you 'start from the beginning'. Also, starting from the beginning does not mean you tell a chronological story of everything that happened for the database curation. Start with the finished product, then tell us what we need to know about how it could be re-created.

### DB-21 · line 18 · comment

**Highlighted text:**

```latex
note
```

**Advisor — 21 October, 11:37 am:**

> what are notes?

### DB-22 · line 18 · comment

**Highlighted text:**

```latex
A note is written from the paper's own text, from the abstract at minimum and from the full text when available, never from the curator's recollection of the paper;
```

**Advisor — 21 October, 11:39 am:**

> You are giving details before giving us the big picture. What does that database look like? What is in it? How is it organized?

### DB-23 · line 22 · comment

**Highlighted text:**

```latex
approved
```

**Advisor — 21 October, 11:40 am:**

> approved by who?

### DB-24 · line 22 · comment

**Highlighted text:**

```latex
left as off-topic or orphaned
```

**Advisor — 21 October, 11:40 am:**

> how was this decided?

**Advisor — 21 October, 11:41 am:**

> how were they identified?

### DB-25 · line 24 · comment

**Highlighted text:**

```latex
the laboratory
```

**Advisor — 21 October, 11:41 am:**

> don't you mean you? Take credit in your dissertation.

### DB-26 · line 24 · comment

**Highlighted text:**

```latex
curation tool, not a reader
```

**Advisor — 21 October, 11:42 am:**

> huh?


## 3. Introduction — chapters/10-introduction.tex

### INTRO-01 · line 8 · comment

**Highlighted text:**

```latex
Artificial
```

**Advisor — 14 October, 10:59 am:**

> This paragraph does a great job identifying one of the gaps of knowledge you worked to fill, but there are a few others you can worked on too.

This is the first word of the BPA motivation paragraph; the advisor's comment concerns the paragraph's gap statement, not the word itself.


## 4. Background — chapters/15-background.tex

### BG-01 · line 48 · tracked edit

**Advisor's tracked edit:**

```text
Changed: requirement
to improvement
```

### BG-02 · line 303 · tracked edit

**Advisor's tracked edit:**

```text
Deleted: these tools were
```

This deletion is anchored at the end of the current file (line 303), rather than recoverable in its original sentence. Treat it as a dangling tracked deletion; no original sentence is inferred.


## 5. Does 15-background.tex already answer the framing comments?

**Partly.** The background already supplies the broad scientific context and a proposed gap, but it does not complete the chapter-specific argument about software capability, prior libraries, and demonstrated scope. The judgments below compare the two current chapter texts; they are not an independent verification of the cited literature.

| Advisor concern | Existing background material | Assessment and use in SNS |
|---|---|---|
| What kind of controller is this, and why use it? (SNS-01/02) | Lines 180 and 194: two-level CPG organization and sensory channels. Lines 253–257: executable neural circuits and their verification phenomena. Lines 261–263: SNS lineage and SNS-Toolbox. | Already explained substantially. Cut the long SNS definitions paragraph (line 6) to essential terms and cross-reference these sections. |
| How does this relate to reflex and learned controllers? (SNS-02) | Lines 269–273 compare reflex-based, learned, and explicit SNS approaches. | Already answered at background level. Replace “are treated there” with a specific cross-reference and one sentence identifying the chosen approach and reason. |
| What gap does this dissertation address? (SNS-01/07) | Line 273 explicitly identifies the attempted combination of an anatomical two-level controller, validated pneumatic mechanics, and an eventual robot. SNS line 8 already repeats it. | The ingredients are present. Put a bounded, concrete gap statement at the beginning of SNS, followed immediately by what you implemented and what the tests establish. “Unclaimed” and “always” in the background are broad novelty claims; narrow them to the work reviewed unless a systematic comparison supports them. |
| Why several software environments? (SNS-06) | Lines 277–279 describe OpenSim, SCONE, AnimatLab, MuJoCo, MyoConverter, and SNS-Toolbox and motivate simulation throughput. | Partly answered. This describes functions more than limitations. It does not explain Simulink/Simscape or why a single environment does not meet the project's requirements. Add a short capability comparison, not another software inventory. |
| How are the new tools distinct from existing tools? (SNS-14) | Background establishes software and controller lineage, but no detailed comparison of prior and new capabilities. | Still needs work. State each missing interface or capability, its consequence for this project, your addition, and the demonstration that tests it. |
| How does the Simulink library compare with Hunt/Fletcher's? (SNS-15) | No Fletcher comparison or description of these libraries appears in the background. Hunt citations discuss other contributions. | Unanswered. Inspect the actual libraries and their documentation before claiming new component functionality. A comparison table should distinguish reused equations, implementation differences, integration features, and genuinely added capability. |
| Why does this work matter for experiments? | Lines 283 and 301–303 give the ethical and model-validation motivation. | Already covered broadly. SNS needs only a brief link to the specific interventions implemented here. Do not repeat the full ethical case. |

The background itself also needs claim alignment: line 301 says the controller is verified against published stimulation and deletion phenomena, whereas SNS line 234 explicitly says the pattern-formation and motoneuron tests have not been run. Specify which phenomena have been tested and retain the others as planned work.

## 6. Proposed changes — for review, not applied to chapter prose

### A. Make the chapter answer a small set of questions

The advisor is asking for a research argument even for a tools contribution. The useful chain is:

**Existing capability → specific limitation → consequence for your research → your addition → evidence → remaining limitation.**

A missing capability can motivate software development without inventing a biological hypothesis. The interfaces and verification procedures can be the contribution. Avoid treating a list of software combinations as sufficient proof of a field-wide knowledge gap.

Use your five combinations from the comment reply as candidate contributions, then consolidate them:

| Candidate gap, bounded to this project | Your addition | Evidence to retain in the main chapter |
|---|---|---|
| Existing model/controller files do not by themselves provide a verified closed loop on both the human and BPA robot representations. | Model conversion/repair, routing checks, BPA force interface, SNS coupling, and explicit identification of which body each controller drives. | One pipeline figure, geometry/force checks, a clear closed-loop example, and the supported-motion limitations. |
| A graphical circuit-building workflow is needed alongside the script-based simulation workflow and CAD-derived mechanics. | Simulink component implementation, controller exports, mechanical imports, and the MuJoCo bridge. | Hunt/Fletcher comparison, one readable library/circuit figure, reflex response, bridge comparison, and the incomplete imports identified honestly. |
| Reimplementation and coding changes can alter circuit behavior, so a controller port cannot be assumed equivalent. | Shared verification protocols across engines and the graded/spiking comparisons. | One agreement table and one informative failure example. Distinguish matched externally driven behavior from autonomous oscillation. |

The human-model CPG coupling, BPA robot coupling, CAD-based mechanism, and Li/Deng/contact combination must each be labeled **implemented**, **demonstrated**, **partially assembled**, or **planned** from the underlying records. The current text mixes models with similar names; a single model–controller–environment–status table would prevent readers from mistaking one demonstration for another.

### B. Replace the opening with a short gap-led introduction

Proposed draft (adapt after checking the prior-library comparison and the exact model identities):

> The preceding chapters establish how braided pneumatic actuator force and knee torque can be predicted. To study how neural control interacts with these mechanics, I also needed a simulation environment that connects anatomically organized spinal circuits to the robot's muscle routing and BPA force model. The existing tools reviewed in the background provide the neural, biomechanical, and physics components separately; using them for this robot requires implemented interfaces and tests of the behavior preserved through conversion.
>
> In this chapter, I develop and evaluate those interfaces in MuJoCo with SNS-Toolbox and in a graphical Simulink workflow. I compare controller implementations across environments and examine whether selected behaviors persist when graded neural components are replaced with spiking components. The contribution is the integration and its verification, with the distinction from existing Simulink libraries established by a direct comparison.
>
> The demonstrations establish specific component and closed-loop capabilities. They also identify failures of behavior transfer and limits of the assembled models. Supported stepping and progression are distinguished from unsupported walking, which has not been demonstrated.

This replaces SNS lines 4–12 with a focused opening, rather than adding another introduction ahead of them. Cite the relevant background sections and supporting sources in the LaTeX version. The draft intentionally avoids claiming first-ever software functionality before the prior-library comparison.

### C. Cut the main narrative substantially

Aim initially to remove or relocate **about 40–50% of the main narrative**, preserving the evidence needed to assess the contributions. This is an editorial target, not a measured reduction already performed. The current source contains about 12,600 whitespace-separated tokens/words including LaTeX, tables, and captions; that is not a clean prose word count.

| Current SNS material | Proposed action |
|---|---|
| Opening, lines 4–12 | Replace with the short gap/contribution/scope opening. Delete the extended roadmap and move basic definitions to a cross-reference. |
| Prior AnimatLab walker, lines 14–20 | Reduce to one paragraph: what existed, its relevant limitation, and the resulting requirement. Keep detailed architecture/history only if needed to interpret a retained comparison. |
| Existing/new tools, lines 33–43 | Replace the overlapping inventories with one capability-gap-contribution table. Explain the environments through their roles and limitations. |
| SADb paragraph, line 37 (SNS-13) | Delete the database composition and curation recap. If needed, retain one sentence pointing to the precise records used for sensory encoding and pathway selection. |
| AnimatLab organization and GUI views, lines 45–86 | Keep the clearest circuit/body architecture view. Move redundant GUI screenshots, instrumented-run provenance, and assembly chronology to an appendix or repository documentation. |
| Conversion, component library, body/controller ports, lines 88–142 | Keep interfaces, equations or references needed for reproduction, audit criteria, and key defects that affect interpretation. Move software counts, script mechanics, environment switches, and import troubleshooting to an appendix. |
| Tuning curriculum, currently lines 410–421 in Discussion | Move the procedure and scoring to Methods. Summarize the stages once. Discussion should interpret the performance and failure modes. |
| Results, lines 166–267 | Group by the questions above. Keep quantitative results and failures together. Move neural/body logging details and repeated “fundamental demonstration” justifications out of Results. |
| Spiking methods/results/discussion, lines 144–156, 269–408 | Condense to motivation, calibration/design, principal outcome, and limitation. Retain a compact comparison table and one revealing figure; relocate exhaustive gate records, block-level residuals, conversion variants, and repeated baseline confirmations. |
| BPA integration, lines 423–425 | Put the force-interface method in Methods; put quantitative demonstration outcomes in Results. Keep only the implications of quasi-static actuation in Discussion. |
| Robot and walking status, lines 427–437 | Consolidate into one results/status table. State the support protocol once per distinct experiment. Keep a short, explicit unsupported-walking limitation. |
| Closing summary, lines 439–443 | Replace repeated inventories with one paragraph stating which gaps were addressed, what the evidence supports, and what remains open. |

Do not achieve the cuts by removing the limitations. Remove repeated inventories, development chronology, software provenance, and repeated defenses of the work's value. Some new comparison material is necessary, but it should replace existing exposition.

A lean structure would be: **Motivation and contribution; Methods; Results organized by contribution; Discussion and limitations.** Fold the spiking comparison into that structure unless its scientific importance warrants a separate major section.

### D. Resolve the localized SNS comments

| Inventory item | Proposed response |
|---|---|
| SNS-01/02 | Put the gap first; cross-reference background for CPG, SNS, and controller alternatives. |
| SNS-03 | Use “software tool development” or a concrete phrase naming the interface/library. Avoid repeatedly calling the work “tooling.” |
| SNS-04/05 | The advisor proposes “virtual” and “this.” These are local wording edits already visible as suggestions; leave acceptance to your review. |
| SNS-06 | Explain the complementary roles of fast scripted closed-loop simulation, graphical circuit construction, and CAD-derived mechanics. Support any comparative speed, accuracy, or tendon-model claim with an actual benchmark or documentation; your comment reply alone is not that evidence. |
| SNS-07 | Move statements about portability, trust, and implications to Discussion; retain the actual purpose and testable contribution in the introduction. |
| SNS-08 | Use a neutral prior-work heading; keep the distinction between prior preliminary work and current results in the prose. |
| SNS-09 | Replace “verification-first approach” with the concrete connection: manual tuning left sensitivity unmeasured, so the current procedure records performance under parameter and input perturbations. |
| SNS-10 | Replace “build order” with an outline of method components. Delete construction-order narration. |
| SNS-11 | Redraw or simplify pipeline.pdf. Separate conversion, closed-loop simulation, and planned hardware deployment; reduce text inside boxes. This report records the comment but does not repair the figure. |
| SNS-12/13 | Remove the generic disclaimer about pre-existing tools and the database recap. Cite prior tools where used and distinguish your additions directly. |
| SNS-14/15 | Add the substantive limitations/prior-library comparison before claiming novelty for the new tools. |

### E. Correct internal inconsistencies before polishing

These are additional reading findings, not advisor quotations or independently audited simulation results:

1. **Scoring direction:** line 421 says the no-cycle fallback sits below every genuine gait so a frozen model cannot win. Lines 212 and 236–239 identify a −25 fallback, while genuine gait scores in the tables are around −144 to −279 and higher is better. Those statements conflict. Explain an independent validity gate if one exists; otherwise the scoring rule needs correction. Do not merely change the wording to hide the conflict.
2. **Replication versus validation:** cross-engine agreement verifies an implementation under the tested protocol; it does not by itself validate biological fidelity or unsupported walking. Use those terms consistently.
3. **Spiking scope:** line 146 says the controllers above all compute with graded neurons, while line 142 describes spiking rhythm cells in the W2L core. Name the network and layer being converted instead of treating all families as uniformly graded.
4. **Prediction versus observation:** line 404 calls a period-matched spiking controller “a short tuning effort away.” Replace that prediction with the measured mismatch and the specific tuning still required.
5. **Method/results placement:** the curriculum belongs in Methods, and the BPA walker outcomes belong in Results. Their current placement makes the chapter's central achievement hard to locate.

### F. Other-chapter follow-through

For the database chapter, the same editorial principle applies: show the finished database and its structure first, then give reproducible selection/curation methods and results. Replace the audit/import/pruning chronology with sources, scope, inclusion/exclusion criteria, extraction fields, quality checks, responsibility, and final counts. Define “notes,” remove software-data jargon such as “first-class record,” distinguish the Airtable working base from the distributed explorer, confirm the presentation history, and identify who approved the curation decisions. Use “I” where you performed the work. The advisor's Cochrane link is preserved in DB-18; this report does not claim to have evaluated that external guide.

For the introduction, retain the strong BPA motivation and add concise corresponding gap statements for the joint mechanics, database, and SNS contributions. For the background, retain the conceptual material but avoid repeating it in SNS; align the verification claims with the tests actually reported.

## 7. Recommended revision order

1. Verify the Hunt/Fletcher comparison and model–controller–environment identities.
2. Replace the SNS opening and establish the contribution/status table.
3. Reorganize Methods and Results around those contributions.
4. Move supporting technical detail to the appendix and remove repeated material.
5. Resolve the scientific/score inconsistencies, simplify the pipeline figure, and then edit sentence-level wording.

This report proposes the revisions; it does not rewrite the dissertation chapter or accept/resolve the advisor's review items.
