# Dissertation audit — pending only, v2

FINAL: current source anchors verified; approved figure corrections closed after integrated r2 checks

25 pending groups: 17 dissertation/cross-document groups and 8 companion-paper groups. These are proposed actions or verification questions, not new approvals. Existing numerical results remain frozen.

The original 83-page audit is preserved unchanged. Resolved, author-closed, disproven, and excluded items are omitted here; their crosswalk is retained separately. Chapter references below use current canonical one-based source lines. Bright yellow marks each pending heading and action in the final PDF.

## Dissertation and shared material

### D01 — Document the historical gait score behind each result family

**Status:** Existing-record verification. **Original IDs:** SNS-S01.

Appendix D still describes zero as perfect agreement and calls the tuning errors mean-removed, but does not identify the exact historical objective, its terms, or its version. Earlier code review identified additional mean-angle and bilateral-bonus terms; the link between that revision and every saved table remains unestablished.

**Pending action/decision:** Match the saved result families to their archived scorer. Give a short equation or term/weight table and distinguish the no-cycle sentinel from a genuine gait score.

- §D.3 Tuning and Supported Evaluations — [95-AppendixD.tex:127](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex:127>) — exact anchor: `The cycle-shape metric assigns zero to perfect agreement.`.
- §D.3.2 Direct Kinematics Overlay — [95-AppendixD.tex:157](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex:157>) — exact anchor: `These errors retain posture offsets and are distinct from the mean-removed cycle-shape errors used in tuning.`.

*Boundary:* Use existing records; preserve reported scores. No tuning or simulation is requested.

### D02 — Define duty and knee-minimum denominators for the zero-cycle case

**Status:** Existing-definition verification. **Original IDs:** SNS-S02.

The supported-ground table says knee minima and contact duty are computed over detected cycles, yet its non-spiking left leg has zero detected cycles and a contact duty of 0.985.

**Pending action/decision:** Check the saved diagnostic definition. If duty is a time fraction over the scored window, label it that way; state how the reported knee minimum is selected.

- §D.4.1 MuJoCo Network and Ground Evaluation — [95-AppendixD.tex:222](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex:222>) — exact anchor: `Knee minima and contact duty are computed over detected cycles.`.
- §D.4.1 MuJoCo Network and Ground Evaluation — [95-AppendixD.tex:229](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex:229>) — exact anchor: `Gait cycles, right / left`.
- §D.4.1 MuJoCo Network and Ground Evaluation — [95-AppendixD.tex:231](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex:231>) — exact anchor: `Contact duty, right / left`.

*Boundary:* Clarify the caption and definition without recalculating the run.

### D03 — Identify the earlier spiking knee comparison as a startup configuration

**Status:** Clarification proposal. **Original IDs:** SNS-S03, O-A08.

The five-second spiking comparison retains “matched operating point” and “settle” terminology. The main chapter now explicitly describes a different, broadened-ramp graded configuration that settles before its pulse. The earlier traces should not imply replication of that later experiment.

**Pending action/decision:** Label the retained comparison as the earlier startup test, use “final-window” where steady settling was not demonstrated, and briefly state its angle convention and relation to the later transient configuration.

- §D.4.2 Simulink Demonstrations — [95-AppendixD.tex:251](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex:251>) — exact anchor: `Knee-reflex demonstration, graded baseline versus spiking twin at the matched operating point`.
- §6.4.3 SNS-Simscape Block Library and Editor — [60-sns.tex:141](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/60-sns.tex:141>) — exact anchor: `The displayed configuration broadens the MN and Ia synaptic saturation ramps;`.

*Boundary:* Keep the earlier negative/transient evidence and all recorded values; do not request a new spiking experiment.

### D04 — Provide one compact model/controller/support map

**Status:** Presentation decision. **Original IDs:** SNS-S08, LATER-37.

The chapter distinguishes the model families in prose, but Table 6.1 is a tool map rather than a result-to-configuration map. Appendix D still spans human, W2L, spiking, and CAD configurations with different actuator laws and support conditions.

**Pending action/decision:** Consider one appendix table mapping model/mechanics, controller ID, actuator law, assistance, evaluated outcome, saved configuration, and relevant figure/table. Replace some inventory prose with that table; identify “production configuration” precisely.

- §6.4.2 Gait2392-Based Models — [60-sns.tex:122](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/60-sns.tex:122>) — exact anchor: `Removing Ia feedback eliminated detected cycles in the tested production configuration`.
- §D.3.1 Reference-Library Comparison — [95-AppendixD.tex:133](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex:133>) — exact anchor: `The production network (`.
- §D.5.4 CAD-Derived Robot Integration — [95-AppendixD.tex:358](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex:358>) — exact anchor: `The 262-neuron controller receives 230 inputs`.

*Boundary:* Do not turn different support configurations into a matched support-weaning claim or remove failed outcomes.

### D05 — Enlarge the inspectable SNS diagram details

**Status:** Figure-polish decision. **Original IDs:** SNS-S05.

The retained neuron implementation and the knee inset use much smaller labels than the surrounding explanatory content. These are the SNS library/knee assets, separate from the now-enlarged Chapter 2 comparison and circuit figures.

**Pending action/decision:** Improve internal area allocation: give the neuron signal-flow/ports a readable crop and enlarge the knee inset within its existing composite. Retain the twelve-block overview and all six feedback connections. Also change the existing knee-demo panel citation from `Figure~\ref{fig:sns_reflex_current}C` to `Figure~\ref{fig:sns_reflex_current}\textbf{(C)}`, preserving its figure label.

- §6.4.3 SNS-Simscape Block Library and Editor — [60-sns.tex:149](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/60-sns.tex:149>) — exact anchor: `{sns_library_implementation_current.pdf}`.
- §6.4.3 SNS-Simscape Block Library and Editor — [60-sns.tex:156](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/60-sns.tex:156>) — exact anchor: `{knee_transient_current_layout.pdf}`.

- §6.4.3 SNS-Simscape Block Library and Editor — [60-sns.tex:141](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/60-sns.tex:141>) — exact anchor: `Figure~\ref{fig:sns_reflex_current}C`.

*Boundary:* This is an unapproved presentation proposal; no model or connection changes are proposed.

### D06 — Show the contact-phase result in the W2L comparison itself

**Status:** Figure-annotation decision. **Original IDs:** SNS-S10.

The text supplies the central 0.498 versus 0.223 contact phases and 0.35–0.65 criterion. The retained figure caption declares the coordination failure, while the graphic relies on dense traces and differing force-axis limits.

**Pending action/decision:** Annotate the existing columns with the saved phase statistic and acceptance interval. Use matching comparison-axis limits where appropriate; an inset is optional.

- §6.4.1 Li--Deng/W2L Walker and Controller — [60-sns.tex:83](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/60-sns.tex:83>) — exact anchor: `disabling all sensory feedback shifted median bilateral contact phase from 0.498 to 0.223 cycles.`.
- §6.4.1 Li--Deng/W2L Walker and Controller — [60-sns.tex:95](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/60-sns.tex:95>) — exact anchor: `{w2l_feedback_verified.pdf}`.

*Boundary:* Preserve the saved traces and metrics; no retuning or new simulation.

### D07 — Separate optional software acceptance bookkeeping from appendix evidence

**Status:** Optional presentation decision. **Original IDs:** APP-A5.

Appendix B still includes the outstanding graphical-acceptance ruling and several full 17-digit identity scores. The headless-test limits and identity result matter; the project-task status and long strings need not occupy the scientific account.

**Pending action/decision:** If desired, retain the headless-only scope and pass/fail outcomes here, and move GUI task status/full identity strings to the linked verification report. Do not silently round or delete recorded numbers.

- §B.4.3 Verification Gates of the MuJoCo Mirror — [93-AppendixB.tex:152](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/93-AppendixB.tex:152>) — exact anchor: `Non-spiking bit-identity & PASS`.
- §B.4.4 AnimatLab Conversion Variants — [93-AppendixB.tex:200](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/93-AppendixB.tex:200>) — exact anchor: `graphical verification (or an equivalent acceptance ruling) remains outstanding`.

### D08 — Use reader-facing status terms in Appendix D/B

**Status:** Wording proposal. **Original IDs:** SNS-S11.

The surviving tables/prose still use “dead,” “no / up,” “zero hand drive,” and implementation aliases. These are less precise than the actual outcomes already recorded.

**Pending action/decision:** Use explicit outcomes such as “no sustained alternation” or “zero external neural drive”; separate the yes/no non-finite/fall fields. Keep alias bookkeeping in technical notes unless needed to interpret the model census.

- §D.4.1 MuJoCo Network and Ground Evaluation — [95-AppendixD.tex:233](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex:233>) — exact anchor: `Non-finite states / fell`.
- §D.5 Supported BPA-Robot Evaluations — [95-AppendixD.tex:304](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex:304>) — exact anchor: `retained bilateral stepping down to zero hand drive.`.
- §D.5.4 CAD-Derived Robot Integration — [95-AppendixD.tex:360](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex:360>) — exact anchor: `four named ankle indices alias knee-PF cells`.
- §B.4.4 AnimatLab Conversion Variants — [93-AppendixB.tex:177](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/93-AppendixB.tex:177>) — exact anchor: `dead (2 spikes / 5 s)`.

*Boundary:* Preserve negative results and the correct neuron counts.

### D09 — Improve the remaining route-geometry and transmission figure details

**Status:** Figure/readability decision. **Original IDs:** T14, APP-A3.

The main A/H cam-pose caption still does not state coordinate axes/units. The Appendix C eight-pose route view, small central 3-D bone view, and first transmission plots retain the original dense area allocation and labels. These are separate from the current 4.11/4.14 plot fixes.

**Pending action/decision:** Supply the main pose frame/axes/units, uncrowd contact labels, and enlarge meaningful detail at final size. Crop the bone view to the knee or omit it if it adds no unique geometry; reduce redundant poses before adding figures.

- §4.3.3 Redesigned Routes for the Biomimetic Robot Model — [30-torque.tex:715](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex:715>) — exact anchor: `Cam-extensor geometry at`.
- §C.6 Redesigned Extensor Route Geometry — [94-AppendixC.tex:325](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/94-AppendixC.tex:325>) — exact anchor: `{route_geometry_ext20mm_20260925.pdf}`.
- §C.6 Redesigned Extensor Route Geometry — [94-AppendixC.tex:333](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/94-AppendixC.tex:333>) — exact anchor: `{bone_mesh_static_extensor20mm.png}`.
- §C.8 Supporting Route Searches and Transmission Sweeps — [94-AppendixC.tex:406](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/94-AppendixC.tex:406>) — exact anchor: `{pulley_torque_resolved_publication.pdf}`.
- §C.8 Supporting Route Searches and Transmission Sweeps — [94-AppendixC.tex:413](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/94-AppendixC.tex:413>) — exact anchor: `{pulley_margin_resolved.pdf}`.

*Boundary:* Keep all route coordinates and plotted results unchanged. Direct-route and cam views are distinct configurations; do not assume they are duplicates.

### D10 — Remove repeated full-corpus inventories from Appendix E

**Status:** Selective compression decision. **Original IDs:** DB-S5, SADB_AUDIT—Appendix relevance (unnumbered P3).

Appendix E repeats Chapter 5’s full 943-paper field inventory and citation-graph counts before giving the distinct retained-set counts. Its III/IV signal cell also contains a long pathway-convergence explanation.

**Pending action/decision:** Use short Chapter 5 references in place of the repeated full-set lists; retain the unique 561-paper maintenance counts and denominator distinctions. A two-column coverage table is optional. Consider moving the III/IV convergence explanation to the pathway-role column.

- §5.2 Results: Application and Corpus Coverage — [50-database.tex:31](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/50-database.tex:31>) — exact anchor: `Within the application's Papers category, containing 943 entries,`.
- §E.1 Record Sets and Coverage — [96-AppendixE.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/96-AppendixE.tex:10>) — exact anchor: `Coverage of the application's Papers category is reported in Chapter`.
- §E.1 Record Sets and Coverage — [96-AppendixE.tex:14](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/96-AppendixE.tex:14>) — exact anchor: `The citation graph covers the 943-paper audit:`.
- §E.2 Afferent Classes and Linked Pathways — [96-AppendixE.tex:33](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/96-AppendixE.tex:33>) — exact anchor: `III/IV & free nerve endings`.

*Boundary:* The earlier “five text pages before the first image” diagnosis is closed; the current chapter shows the app earlier.

### D11 — Identify the retrievable SADb snapshot

**Status:** Availability/provenance decision. **Original IDs:** DB-S6.

The chapter describes the embedded 1,393-record app and monthly watch but still gives no retrievable application release or snapshot identity. A working local browser session is different from an artifact a dissertation reader can obtain.

**Pending action/decision:** Add one verified availability statement for the existing app/repository and the snapshot behind the reported counts/captures. Include an existing version/date if recoverable; do not invent a release or new hosting.

- Chapter 5 opening — [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/50-database.tex:8>) — exact anchor: `The primary deliverable is SADb Explorer,`.
- §5.1 Curation Methods — [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/50-database.tex:25>) — exact anchor: `I implemented a monthly literature watch`.
- Appendix E opening — [96-AppendixE.tex:4](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/96-AppendixE.tex:4>) — exact anchor: `This appendix records the curation snapshot`.

### D12 — Add only the useful local terminology reminders

**Status:** Optional gloss/title decision. **Original IDs:** LATER-34, LATER-35, LATER-36, LATER-42.

Remaining small reader aids are the local meaning of l620, identification/citation of the Hunt force-model baseline, the “Mono- and Bi-Pulley” heading, and OpenAlex’s role as a bibliographic index.

**Pending action/decision:** Consider a brief reminder/cross-reference at l620, a local Hunt citation, and a mechanical reeving heading or definition that distinguishes BPA count, cable parts, and sheave count. Give OpenAlex a short first-use gloss if helpful.

- §4.1.3 Hybrid Torque Calculation Method — [30-torque.tex:83](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex:83>) — exact anchor: `$l_{620}$ and $l_\text{m}$`.
- §4.2 Results — [30-torque.tex:472](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex:472>) — exact anchor: `The pinned-flexor results first compare the Hunt actuator model`.
- §4.3.5 Mono- and Bi-Pulley Tendon Routing — [30-torque.tex:764](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex:764>) — exact anchor: `\subsection{Mono- and Bi-Pulley Tendon Routing}`.
- Chapter 5 opening — [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/50-database.tex:8>) — exact anchor: `Search (online) finds candidate additions through OpenAlex.`.

*Boundary:* The pulley equations already define n and G correctly; no new transmission calculation is requested. Published Chapter 3 stays protected.

### D13 — Check the existing reuse records for three reproduced images

**Status:** Source-specific verification. **Original IDs:** C2-O15, C2-O16, C2-O17, LATER-56, LATER-57, LATER-58.

Kengoro, Rybak 2006 Figure 1, and Shevtsova’s Figure 5.8 have different source/credit contexts. The previous primary-source review did not establish a common CC BY grant or determine whether existing permissions were already held.

**Pending action/decision:** Locate each actual image’s source and existing reuse/permission basis, then use the required credit. Keep the correct 2006 Rybak source and the Springer chapter’s Figure 5.8 attribution.

- §2.1 Biomimetic and Musculoskeletal Robots — [15-background.tex:19](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/15-background.tex:19>) — exact anchor: `Kengoro human-mimetic musculoskeletal humanoid in neutral`.
- §2.5 Neural Control of Locomotion: Central Pattern Generators — [15-background.tex:169](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/15-background.tex:169>) — exact anchor: `Reproduced from \citet[Fig.~1]{rybak_modelling_2006}`.
- §2.5 Neural Control of Locomotion: Central Pattern Generators — [15-background.tex:181](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/15-background.tex:181>) — exact anchor: `Reproduced from \citet[Fig.~5.8]{shevtsova_modeling_2016}`.

*Boundary:* This is not a finding of unlawful reuse. Do not apply the 2015 review’s license to the 2006 article or add blanket CC BY wording.

Sources: [Kengoro primary record](https://pubmed.ncbi.nlm.nih.gov/33157878/), [Rybak 2006 source](https://pmc.ncbi.nlm.nih.gov/articles/PMC1890439/), [Shevtsova chapter](https://link.springer.com/chapter/10.1007/978-1-4939-3267-2_5).

### D14 — Confirm the torque manuscript author list and title

**Status:** Human metadata decision. **Original IDs:** T11.

The dissertation’s publication note/bibliographic record names seven torque co-authors and an older title. The companion paper names four authors and uses “Design Factors for Isometric Robotic Knee Torque Using Braided Pneumatic Actuators.”

**Pending action/decision:** Have the authors confirm the intended author list, contribution statement, title, and manuscript status, then align the dissertation note, bibliography, and paper metadata.

- Chapter 4 opening — [30-torque.tex:5](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex:5>) — exact anchor: `Publication note: This chapter is the subject of an unpublished manuscript`.
- Paper front matter — [Main.tex:48](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:48>) — exact anchor: `\def\Authors{`.
- Paper front matter — [Main.tex:69](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:69>) — exact anchor: `\title[Design Factors`.
- References — [thesis.bib:1113](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/thesis.bib:1113>) — exact anchor: `@unpublished{bolen_isometric_2025`.

*Boundary:* No agent should select or remove co-authors.

### D15 — Choose the scope of prior-publication/presentation notes

**Status:** Author decision; no insertion authorized. **Original IDs:** PUBLICATION_NOTE_RECOMMENDATION.

Chapter 5 already identifies the 2023/2024 presentations. Chapter 4 cites the earlier Morrow–Bolen–Hunt placement work and identifies the new unpublished torque manuscript. A citation alone does not tell whether a chapter actually reproduces/adapts prior co-authored material. The 2023 poster record also pairs an SfN note with a NeuroNex URL.

**Pending action/decision:** Choose one concise presentation credit for Chapter 5 rather than duplicating its opener. Confirm the 2023 record/link. For Chapter 4, add a scoped published-material/contribution note only if prior co-authored material is incorporated; merely citing related work does not by itself establish that case.

- Chapter 5 opening — [50-database.tex:4](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/50-database.tex:4>) — exact anchor: `The resource and the evaluation of AI-assisted literature curation were presented`.
- §4.1.5 Actuator Placement and Route Redesign — [30-torque.tex:342](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex:342>) — exact anchor: `\citep{bolen_determination_2019,morrow_optimization_2020}`.
- References — [thesis.bib:1132](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/thesis.bib:1132>) — exact anchor: `@misc{bolen_sensory_2023`.

*Boundary:* PSU’s requirement is conditional on incorporated material; a poster does not automatically make the entire later chapter a published paper. The detailed source/Zotero recommendation is preserved separately.

Sources: [PSU ETD handbook, March 2026](https://www.pdx.edu/gradschool/sites/gradschool.web.wdt.pdx.edu/files/2026-03/ETD_Formatting_Requirements_and_Format_Review_Handbook_March_2026.pdf).

### D16 — Choose limited compression of the mechanical methods

**Status:** Optional compression decision. **Original IDs:** T13, LATER-46.

Chapter 4 and the paper retain successive scalar/component substitutions in the compliance development. Chapter 4 also retains the legacy cube-search steps and numerical weights before the distinct later mesh objective.

**Pending action/decision:** If shortening is useful, inline/move only repetitive substitutions and legacy search implementation details. Keep the advisor-requested bracket mechanics, stiffness projection, series relation, implicit force balance, master’s/Morrow attribution, and distinction between the published cost and the later mesh objective.

- §4.1.4 / Length Offset and System Stiffness (unnumbered subheading) — [30-torque.tex:91](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex:91>) — exact anchor: `\subsubsection{Length Offset and System Stiffness}`.
- Paper §2.6 Improved Torque Model: Length Offset and Bracket Compliance — [Main.tex:197](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:197>) — exact anchor: `\subsection{Improved Torque Model: Length Offset and Bracket Compliance}`.
- §4.1.5 Actuator Placement and Route Redesign — [30-torque.tex:344](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex:344>) — exact anchor: `At each scale the search perturbs both attachments`.
- §4.1.5 Actuator Placement and Route Redesign — [30-torque.tex:346](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex:346>) — exact anchor: `The original cost combined the summed absolute torque mismatch`.

*Boundary:* Do not delete the author-retained route-builder chronology or change equations, weights, or results without an explicit decision.

### D17 — Separate the two meanings of rho in the contact-release description

**Status:** Notation-fix proposal. **Original IDs:** OP07, LATER-49.

The same paragraph uses rho(theta) for offset-curve radius of curvature and rho = 3 mm for collision-envelope relaxation. This collision survives in both Appendix C and the paper supplement.

**Pending action/decision:** Use a distinct symbol for the fixed clearance relaxation, retaining rho(theta) for curvature. Change the prose/notation consistently in the two documents and preserve the 3 mm value.

- §C.7 Contact and Release Angles of the Wrapped Route — [94-AppendixC.tex:377](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/94-AppendixC.tex:377>) — exact anchor: `each envelope radius may be relaxed by $\rho = \qty{3}{\mm}$`.
- Supplement §1.5 Contact and Release Angles of the Wrapped Route — [frontiers_SupplementaryMaterial.tex:251](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/frontiers_SupplementaryMaterial.tex:251>) — exact anchor: `each envelope radius may be relaxed by $\rho = \qty{3}{\mm}$`.

## Companion paper and supplement

### P01 — Update the published force-model reference

**Status:** Bibliographic correction proposal. **Original IDs:** T05.

The companion bibliography still defines bolen_2026 as an unpublished manuscript submitted for publication, although the force paper is published in Actuators 15(5):230.

**Pending action/decision:** Replace the stale record’s type/status with the verified journal metadata and DOI 10.3390/act15050230, keeping the cited key and byte-safe bibliography handling.

- References — [BolenFrontiers22.bib:1085](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/BolenFrontiers22.bib:1085>) — exact anchor: `@unpublished{bolen_2026`.
- Paper §2.2 Actuator Force Model — [Main.tex:112](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:112>) — exact anchor: `The actuator force model is fully described`.

*Boundary:* This changes reference metadata, not the protected published chapter or the force model.

Sources: [Published force paper](https://www.mdpi.com/2076-0825/15/5/230).

### P02 — Check the remaining small pinned/supplementary plot labels

**Status:** Figure-polish decision. **Original IDs:** T07.

The outstanding scope is the six-panel pinned-extensor figure and supplementary direct-route comparison. The prior audit measured legends around 6.3 and 5.3 pt, respectively. The repaired physical pair and extensor legend are separate, completed changes.

**Pending action/decision:** Re-export or allocate more usable area so the surviving plots meet the standing 10-pt final-size target. Preserve every trace and result; do not treat the repaired physical-pair styles as clearing these other assets.

- Paper §3.2 Extensor on Pin Joint — [Main.tex:443](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:443>) — exact anchor: `\label{fig:ExtPin_group}`.
- Supplement §1.3 Figures — [frontiers_SupplementaryMaterial.tex:173](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/frontiers_SupplementaryMaterial.tex:173>) — exact anchor: `\label{fig:RouteVariants}`.

*Boundary:* Confirm final inclusion scale when exporting; the inherited measured values are not a claim that every paper label is undersized.

### P03 — Decide whether the journal should place hardware evidence in Results

**Status:** Optional paper-structure decision. **Original IDs:** T08.

The redesigned-route capacity and physical comparisons remain within Discussion. The cost details have already moved to the supplement, and the capacity/physical passages already lead with their important results.

**Pending action/decision:** For journal organization, consider moving the new hardware capacity evidence into Results and leaving interpretation/trade-offs in Discussion. Preserve the currently approved result-first order and pulley placement after that evidence.

- Paper §4.3.2 Physical Validation Against Human Torque — [Main.tex:600](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:600>) — exact anchor: `\subsubsection{Physical Validation Against Human Torque}`.
- Paper §4.4 Whole-Leg Capability Against Human Torque Surfaces — [Main.tex:623](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:623>) — exact anchor: `\subsection{Whole-Leg Capability Against Human Torque Surfaces}`.

*Boundary:* This is not an instruction to move the specifically named dissertation §4.3.4 or repeat completed U12/U13 work.

### P04 — Replace the remaining identification settings wall with a compact specification

**Status:** Compression proposal. **Original IDs:** T09.

The Identification Procedure still places both stages’ memberships, populations, generations, tolerances, normalization, bounds, and initialization into one long source paragraph. The redesigned-route tables are already completed and are not this issue.

**Pending action/decision:** Consider a two-row stage/settings table, with one short explanatory paragraph. Keep every existing training/validation assignment, bound, solver setting, and adopted stiffness-pair constraint; remove duplicated settings prose rather than adding a decorative figure.

- Paper §2.8 Identification Procedure — [Main.tex:382](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:382>) — exact anchor: `The identification proceeded in two stages.`.

*Boundary:* No new optimization or changes to the FVU analysis.

### P05 — Represent both hardware routes in the paper abstract

**Status:** Framing/optional compression proposal. **Original IDs:** T12, OP05, LATER-54.

The abstract promotes hardware measurements of the redesigned extensor only; the current physical section reports that both redesigned routes met/exceeded their selected references at the measured angles. The abstract, final Introduction paragraph, and Methods Overview also repeat the work sequence.

**Pending action/decision:** Include the existing flexor as well as extensor capacity result in a concise abstract statement, scoped to the tested points/configurations. Optionally trim repeated contribution/process lists between Introduction and Overview.

- Paper front matter — [Main.tex:81](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:81>) — exact anchor: `together with the first hardware torque measurements of the redesigned extensor route.`.
- Paper §4.3.2 Physical Validation Against Human Torque — [Main.tex:602](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:602>) — exact anchor: `Both redesigned routes met or exceeded their selected human-muscle torque references`.
- Paper §2.1 Overview — [Main.tex:108](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:108>) — exact anchor: `A hybrid calculation isolated the error sources,`.

*Boundary:* The scientific segue and 27-site contribution have already been corrected; do not reopen them.

### P06 — Add the remaining pinned-bracket inventory callout

**Status:** Small callout proposal. **Original IDs:** T15, OP14.

The retained pinned-bracket table has an explanatory lead but still no numbered reference to tab:RoutePinnedBrackets. The other retained coordinate inventories and cam figures now have purposeful references; the old incorrect supplement-inventory sentence is fixed.

**Pending action/decision:** Add one explicit Supplementary Table reference to the existing pinned-bracket introduction or relevant main Methods sentence. Do not add redundant calls for items that are already introduced.

- Supplement §1.2 BPA routing and via points — [frontiers_SupplementaryMaterial.tex:80](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/frontiers_SupplementaryMaterial.tex:80>) — exact anchor: `Bracket reference points of the pinned-knee BPA routes, in mm.`.
- Supplement §1.2 BPA routing and via points — [frontiers_SupplementaryMaterial.tex:84](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/frontiers_SupplementaryMaterial.tex:84>) — exact anchor: `\label{tab:RoutePinnedBrackets}`.

### P07 — Make first table mentions follow table order

**Status:** Callout-order proposal. **Original IDs:** OP13, LATER-52.

The adopted-parameter table is still called in the first flexor/extensor Results paragraphs before the goodness-of-fit and biomimetic tables. Its natural introduction remains in Identified Design Parameters.

**Pending action/decision:** Defer/remove the two premature adopted-table references while retaining the parameter values. Verify the resulting first-mention sequence by labels after the build; do not manually renumber tables.

- Paper §3.1 Flexor on Pin Joint — [Main.tex:411](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:411>) — exact anchor: `Its adopted design parameters are a length offset`.
- Paper §3.2 Extensor on Pin Joint — [Main.tex:437](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:437>) — exact anchor: `The extensor identification of record locks the stiffness pair`.
- Paper §4.2 Identified Design Parameters — [Main.tex:551](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex:551>) — exact anchor: `The adopted correction terms are listed in Table`.

*Boundary:* The former premature hybrid-figure callout is already removed.

Sources: [Frontiers figure/table guidance](https://www.frontiersin.org/guidelines/author-guidelines).

### P08 — Optional bibliography source cleanup

**Status:** Optional source-hygiene decision. **Original IDs:** OP11, OP15, LATER-51, LATER-53.

The uncited supplement retains an unused bibliography command. The main bibliography also retains two records for the same Hunt 2017 DOI, but only the non-suffixed key is cited; this does not produce a duplicate printed reference.

**Pending action/decision:** If preparing a submission-only source package, remove the unused supplement bibliography invocation while it has no citations, and consolidate verified Hunt metadata while preserving the cited key. If citations are later added to the supplement, supply its style and an end reference list.

- Supplement §1.3 Figures — [frontiers_SupplementaryMaterial.tex:177](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/frontiers_SupplementaryMaterial.tex:177>) — exact anchor: `\bibliography{BolenFrontiers22}`.
- References — [BolenFrontiers22.bib:369](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/BolenFrontiers22.bib:369>) — exact anchor: `@inproceedings{hunt_modeling_2017,`.
- References — [BolenFrontiers22.bib:1020](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/BolenFrontiers22.bib:1020>) — exact anchor: `@inproceedings{hunt_modeling_2017-1,`.

*Boundary:* Low priority. There is no current visible mid-document bibliography or proven duplicate rendered citation.

