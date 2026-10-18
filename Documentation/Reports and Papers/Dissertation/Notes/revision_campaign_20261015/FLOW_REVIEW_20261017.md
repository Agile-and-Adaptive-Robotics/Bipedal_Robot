# FLOW Review — Full-Dissertation Read, 2026-10-17

Independent committee-style read of every chapter of `ProofFinal/chapters/` (10, 15, 20, 30, 50, 60, 70, 80, 92–94) plus `main.tex` and the abstract, in order, judged on FLOW: chapter openers, section handoffs, duplication, orphaned references, argument order, and residual campaign/administrative language. Line numbers refer to the `.tex` files as they stood at review time (30-torque was modified 05:41 on 2026-10-17 and is under concurrent overhaul — see the overlap note in its section). This review is read-only; no dissertation file was touched.

Mechanical sweep results up front (details inline):

| Flag | 10-intro | 15-bg | 20-mus | 30-torq | 50-db | 60-sns | 70-disc | 80-concl | 92-A | 93-B | 94-C |
|---|---|---|---|---|---|---|---|---|---|---|---|
| "campaign" in prose | 0 | 0 | 0 | 9 | 3 | 12 | 2 | 1 | 0 | 5 | 3 |
| "of record" in prose | 0 | 0 | 0 | 11 | 0 | 3 | 0 | 0 | 2 | 4 | 11 |
| em-dash `---` in prose | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 pair (L194; rest are table cells) | 1 (table cell) |
| "plant" for robot/model | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Em-dash overuse and "plant" are **clean** — no action needed. The live problems are campaign/of-record language, dates and file names in prose, and a set of state mismatches between chapters.

---

## main.tex

- **L87** — `\submitdate{July, 2026}` with submission deadline 2026-10-21. Why it matters: the title page will carry a date that predates submission. Fix: set the actual submission term before the ETD upload.
- **L145–146, L154–159, L162–165, L179–182** — scaffolding comment blocks (motivation chapter, list-of-appendices, structure notes). Invisible in the PDF; harmless, but clean at final pass if the source is shared.

## 10-introduction.tex

- **L18** — "(SNN)" is introduced for "synthetic neural network" and never used again; the rest of the dissertation standardizes on SNS (synthetic nervous system). Fix: drop the abbreviation or use SNS here.
- **L20 vs 15-background L238** — Intro names the quadruped "Muscle Mutt, formerly called DoggyDeux"; Background calls it "the DoggyDeux quadruped." Fix: pick one name and use it in both chapters.
- **L39** — Chapter-4 roadmap says the chapter "then reports the measured torque..., the follow-up identification, the redesigned routes..., before closing on the physical sources of prediction error and the transmission-design implications of mono- and bi-pulley routing." As written this is accurate, but the concurrent 30-torque overhaul (follow-up identification moving to Discussion) will make this sentence stale. Fix: re-check this roadmap after the overhaul lands and move "the follow-up identification" to match its new location.
- **L41** — Chapter-5 roadmap describes the corpus, the packaging, and the triage evaluation, but omits §5.4 "Toward Robot and Simulation Use," which is the chapter's payoff and its handoff to Chapter 6. Fix: add half a sentence ("and closes by mapping the corpus onto the controller's sensory channels").
- **L43** — The Chapter-6 tool list ("the OpenSim-to-MuJoCo conversion..., the Simulink SNS library..., and the cross-platform verification protocol") omits the spiking-neuron mirrors from the tools sentence; they appear only in the following sentence as "demonstrations." Minor: promote the mirrors to the tool list.

## 15-background.tex

- **L4** — Chapter opener is good: states scope, and explicitly parks dissertation work (prior walker, SADb, pipeline) in later chapters. The chapter then violates its own policy once (L256, below).
- **L98–107** — `fig:bpa_graphical_abstract`, "Overview of the BPA force-characterization study," is the dissertation's own graphical abstract placed in the Background chapter. Why it breaks flow: Background ends its prior-models subsection by exhibiting the contribution; Chapter 3 never shows this figure. Fix: move to Chapter 3 or cut.
- **L115** — "the SimTK repository that distributes OpenSim also hosts multiple projects dedicated specifically to the human knee" — vague, uncited claim doing load-bearing work in an argument about modeling choice. Fix: cite one or two named knee projects or delete the sentence.
- **L194** — "the tests planned in Chapter~\ref{ch:futurework}" renders as "Chapter 7," i.e., the Discussion; reads as if a dedicated future-work chapter exists (label survives the old structure). Acceptable, but prefer referencing the specific subsection.
- **L238** — "the DoggyDeux quadruped" (see Intro L20) and "is presented as prior dissertation work in Section~\ref{sec:prior_walker}" — awkward phrase; "prior work of this dissertation" or "as prior work in..." reads cleaner.
- **L256** — Three problems in one sentence: (a) "I used these tools to convert the laboratory's modified OpenSim leg model to MuJoCo, couple it to a synthetic nervous system..., and run the full neuromechanical loop faster than real time" is dissertation results inside the Background chapter, contradicting the opener's own policy (L4); (b) it is the only first-person "I" in the chapter; (c) "run the full neuromechanical loop faster than real time" is contradicted by 60-sns L73 ("Full-body closed-loop co-simulation is not yet running") and is at best stale given Chapter 6's ground evaluations. Fix: reduce to a forward pointer ("these tools are assembled into the pipeline of Chapter 6") and delete the real-time claim.
- **L274** — "A 50-minute classroom recording ... first shows an intact cat playing with its handler ..., and then performs the decerebration for the class" — subject-verb confusion: the recording does not perform surgery. Unclear on first read. Fix: "...shows ... and then records/shares the decerebration performed for the class."
- **L278** — "the neuromechanical controller can be verified against published tonic stimulation and deletion phenomena ... before any locomotion tuning is attempted" — Chapter 6 reports locomotion tuning already performed (four tuned production configurations). Fix: align tense with Chapter 6's actual state ("...phenomena, the verification step that precedes locomotion tuning in this program").

## 20-muscles.tex

- **L1** — Triple alias `\label{ch:muscles}\label{ch:methods}\label{ch:discussion}`. `ch:methods` is used by 60-sns/80-conclusion (renders correctly as Chapter 3 but is a maintenance landmine); `ch:discussion` is referenced nowhere. Fix: keep `ch:muscles` + `ch:methods` at most; delete `ch:discussion`.
- **L7** — Chapter opener accurately describes the chapter; matches the Intro roadmap. Good.
- **L295** — Text: "adjusted R²=0.9998, a RMSE=0.004537, and a maximum absolute residual of 10.6%." Table `tab:Fstar` (10 mm, Model column) reports RMSE 0.005118 and max error 10.3% (10.6% is the Validation column). Fix: reconcile text with the table (likely RMSE and max-error values were mixed between Model and Validation columns).
- **L312** — Two single-value `\qtylist{498.6}{\N}` / `\qtylist{350.9}{\N}` (should be `\qty`); renders oddly with a list separator context.
- **L330** — "We also anticipate the equation will also work" — doubled "also."
- **L405** — `\cite{chou_measurement_1996}` uses `\cite` in a natbib document where everything else is `\citep`/`\citet`. Cosmetic inconsistency in the rendered reference style.
- Voice: "we" throughout, justified by the publication note; no action.

## 30-torque.tex — reviewed as-is (concurrent overhaul in flight)

**Overlap note:** the in-flight overhaul (follow-up identification → Discussion, de-campaigning, added tables) targets exactly §4.5 (L421–462) and the campaign language catalogued below. Line numbers here describe the current file; whoever merges should re-anchor them.

- **L5/L7** — The chapter has **no opening paragraph**: publication note, then straight into `\section{Methods}`. It is the only main chapter without a statement of what it does and why, and it is the longest chapter. Fix: add a 6–10 line opener (actuator model → joint torque → correction terms → route redesign → what is validated against what), matching the Intro's L39 roadmap. (If the overhaul adds one, disregard.)
- **L5** — Publication note: sections "...are **also** the subject of an unpublished manuscript" — "also" is orphaned (nothing prior relates this chapter to a manuscript; inherited from Chapter 3's note shape). Fix: drop "also."
- **L33** — "the flexor is actuated with two BPAs **in with** pathways mirrored about the Z axis" — typo; and the sentence leaves unclear which knee type uses 10 mm single vs 20 mm paired BPAs. Fix the typo and split the configuration statement into one clean sentence per knee.
- **L35** — "Two vectors enter that moment, and each is a quantity with a magnitude and a direction, not a direction alone" — reads as a rebuttal to an objection the reader has not made; unclear purpose on first read. Fix: delete the second clause; the equation does the work.
- **L82** — "Data collected during the first flexor experiments (Chapter~\ref{ch:results})" — `ch:results` is this chapter's own alias; a self-reference that reads as an external pointer (orphan of the old methods/results split). Fix: "(Section~\ref{sec:results_torque})" or just "the first flexor experiments."
- **L254 vs L279 vs L423–425** — The identification is narrated in three stages (single-bracket fit in Methods, two-stage extension, adopted values in §4.5), and the reader absorbs two incompatible value sets: Table `tab:OptCoeff` (χ0 = 11.9 mm, χ1 = 1.02e5, χ2 = 1.06e4) at L338, then the adopted set (+3.94 mm, 3.998e4, 1.473e4) at L425. The disclaimer arrives only at L423. Fix: flag at L338 that these are the initial single-bracket values, superseded by §4.5's adopted set (one sentence).
- **L279** — Wall-of-text paragraph; contains a mid-paragraph `\par Figure~\ref{fig:xiProjection}...` splice. Fix: split at the figure introduction.
- **L318, L431, L433, L435, L437, L440** — File names in main-text prose: `Opt_run.m`, `Opt_run_Ext.m`, `objective_KneeExt20mm`, `objective_KneeFlexor20mm`, `Opt_run_pulley.m`, `OpenSim_Vasti_Results.txt`, `OpenSim_Bifem_Results.txt`, `buildKneeExtContext20mm.m`, `buildKneeFlexorContext20mm.m`, four `.mat` result files, `TorqueZ`/`TorqueInsZ`, `mjc/gait2392_robot`, `add_mag3`. Appendix C already documents these. Fix: strip from the chapter; keep one pointer to Appendix C.
- **L354** — "the new stiffness and length **offet** model" — typo.
- **L363** — "almost equal in magnitude, but opposite sign, **as** the flexor configuration" — grammar ("opposite in sign to").
- **L394/L400 vs Table `tab:TorqueFit` L416** — Text and figure caption say the 20 mm flexor is $l_0 = 42.0$ cm; the table row says 41.5 cm. Fix the table.
- **L421–462 (§4.5)** — Campaign-language epicenter, all in a subsection filed under Results that reads as a status chronicle ("has since been completed," "is likewise far advanced," "has since been redesigned once more"): "2026-09-20 campaign" (L425, 427, 433, 437), "design of record" (L427 ×2, 433, 435), "September design of record" (L433), "front row 77" (L425), "angle index 77 of the one-based MATLAB grid" (L433), "the preflight value of 0.6209 that appears once in the 2026-09-20 campaign log is superseded" (L437), "method of record" (L425), `.mat` filenames (L427, 433, 435, 440), "converted 2026-09-24" / "2026-09-26 pulley batch" (L460). Fix: this is the overhaul's target; strip dates, filenames, and index arithmetic; keep adopted values and outcomes.
- **L425/L427/L433/L435** — The same extensor margin is quoted three times with an improving sequence (+5.043% → +5.02% → cam +30.83%) across four paragraphs. The chronology forces the reader to recompute which number is current. Fix: state the adopted (cam) result once, with one sentence of history.
- **L440–442 vs L488 (§4.6.2)** — The Mono/Bi-pulley results are stated **twice nearly verbatim**: objective 0.08465, exit flag 1, 81.23→85.74 N·m, +4.51 N·m, off-axis ≤0.803 N·m, clearances 4.26–19.14 mm vs 25 mm, sentinel 101, the same five sites. The Discussion version (L488) is the better-written one. Fix: keep §4.6.2; reduce §4.5.1's pulley paragraph to two sentences and a pointer.
- **L460** — "The MuJoCo conversion tree ... predates both the 2026-09-26 pulley batch and the adductor-magnus restoration, so it carries one stale path point" — honest, but it contradicts L427's claim that "These redesign solutions are the model definitions carried into the neuromechanical simulation work." Fix: qualify L427 ("carried into the simulation work, whose converted model tree is being refreshed to match" or equivalent).
- **L462** — The 27-actuator whole-body paragraph reports concatenation errors and crash sets in Results. This is engineering status, not results. Fix: compress to 3–4 sentences and move to Discussion (or the future-work part of ch7), keeping the honest 8-of-15/6-of-8 accounting.
- **L479–484 (end of §4.6.1)** — The dangling oscillator derivation: "In equation $m\ddot\varepsilon^* + c(\dot\varepsilon^*)\dot\varepsilon^* + k(\varepsilon^*)\varepsilon^* = F$ given ... then if we set mass m and external force F to zero we get $c(\dot\varepsilon^*) = (k(\varepsilon^*)\varepsilon^*)/\dot\varepsilon^*$." No lead-in, no follow-up, no stated purpose, ending the Sources-of-Error subsection. A committee member will ask what it is for. Fix: introduce it (damping extraction from the fitted force surface for the dynamic model of Elzein et al.) or cut it.
- **L490** — "measured knee torque of 1.163 N·m at 10°, 7.849 N·m at 11°, and 23.437 N·m at 10°" — first and third at the same angle; the sentence explains neither why two points sit at ~10° nor that each point pairs a different pressure and hung load until the reader unpacks the lists. Fix: present as a 3-row mini table (angle, load, gauge pressure, torque).
- **Terminology** — "the Hunt model" (L394), "the original model," "improved BPA force model," "improved BPA characterization" are used interchangeably for two distinct baselines (Hunt 2017 prediction; pre-correction rigid prediction). Fix: define both once in §4.1 and use two fixed names.

## 50-database.tex

- **L4–L20** — Opening and assembly narrative are strong; the chapter states what it is for and hands off to Chapter 6 explicitly (L12). Best-flow chapter of the set.
- **L14** — "The 2026 audit and curation **campaign** began from a three-way reconciliation" — date + campaign in the first sentence of the narrative. Fix: "The audit and curation effort reported here began..."
- **L24** — "at the close of the **campaign**"; also "43 records below the laboratory's 1,000-record operating cap, with room to grow" — internal lab administration detail. Fix: drop the operating-cap clause.
- **L39** — "The 2026 **campaign** was built as the operational answer to that finding" — same fix; "the curation workflow reported here."
- **L83–85** — Closing reach to `sec:fw_senses` (Chapter 7) is a good forward handoff.
- **Cross-chapter (see below)** — Chapter 5's corpus statistics (561 retained; 523 notes; 86 afferent-classified; 258 translation notes) are re-quoted in Chapter 7 with different numbers and a different denominator. Fix belongs in 70-overall.

## 60-sns.tex

- **L4–L12** — Five-paragraph preamble; thorough but heavy. L6 contains a ~70-word sentence ("The chapter has a single stated goal: a verified neuromechanical laboratory ... meaning a set of tools in which a controller can be built in any of several simulation engines, checked against its own re-implementations ... and then interrogated with the stimulus and deletion experiments ..."). Fix: split; compress preamble to three paragraphs.
- **THE HANDOFF GAP — zero mentions of the Sensory Afferent Database.** Chapter 5 ends by promising that "the afferent classes the controller monitors, the signal encoding each channel carries, and the connection rules its circuit implements all trace to records held here," and §5.4 is titled "Toward Robot and Simulation Use" — but 60-sns contains no occurrence of "Sensory Afferent Database," `ch:database`, or `bolen_sensory_2023` (grep-verified). The dissertation's stated rationale for the database is never cashed in where the controller is built. Fix: one paragraph in the preamble or §6.2 wiring the Ia/Ib/II/contact channels to the SADb pathway records (Table 5.2 exists precisely for this), plus a citation.
- **L73** — "Full-body closed-loop co-simulation is not yet running: leg-joint states develop non-finite values during longer simulations, with the foot-contact meshes the primary suspect under investigation." This is contradicted inside the same chapter: §6.4.4–6.4.6 report 16 s ground evaluations of four production configurations with gait cycles, contact duty, and pelvis tilt, and §6.5 reports a spiking ground gait that "stays up." The tool description understates the state of the work its own results section reports. Fix: update L73 to the resolved/conditional truth (e.g., "early full-body runs developed non-finite leg-joint states, traced to the foot-contact meshes and resolved; the production evaluations of Section 6.4 run full-body on this coupling").
- **L169** — "the four production walking configurations ... (\texttt{s3k}) and its pruned derivative (\texttt{s3kpruned}), the six-synergy variant (\texttt{syn6}) ... and the W2L-layout variant (\texttt{w2lvar})" plus "the production pipeline's cycle-shape kinematic score" — configuration codenames and a scoring pipeline arrive with **no prior description of the tuning workflow that produced them**. Objective 4 of the Introduction promises "a staged tuning and verification workflow"; the chapter documents the verification half and never describes the tuning half (how the gains were tuned, against what reference, with what optimizer). Results-before-methods at chapter scale. Fix: add a short "tuning protocol" paragraph in §6.2 or §6.4.4 (even 6 lines: reference cycle, staged curriculum, optimizer, score definition) or soften the configurations' provenance ("tuned in prior laboratory iterations" with a pointer).
- **L169** — "\texttt{Case\_40\_motion}, the predicted-walking cycle of the Falisse predictive-simulation dataset" — file-name ID used as a proper noun throughout the table and text. Fix: give it a descriptive name at first use and use the filename only in the table.
- **L232** — "re-established on the **machine of record**" and the 17-significant-digit score "$-161.56754173676563$" in prose. Fix: "on the recorded software environment" (or drop) and round to $-161.57$ (the exact value already lives in Appendix B's gate table).
- **"campaign" ×12** — L6 ×2 ("a walking campaign to justify itself," "the heavier campaigns"), L131, L165, L228, L246, L325, L351, L355 ×2 ("a short tuning campaign away"), L359. Most mean "evaluation series" or "study." Fix: replace; keep at most one deliberate use.
- **L121, L325** — "calibration of record," "its period of record" — of-record idiom in main text; "the calibration" / "its documented period" reads better.
- **L6** — meta-commentary about the dissertation's own reporting policy ("reported for the questions they answer, not as evidence that the tools exist") is fine once; avoid repeating the move at L131 ("Failures are reported with the same weight as successes throughout") and L131's echo at L187 ("it is reported as such") — three declarations of honesty in one section reads defensive.
- **L363** — Chapter closing paragraph is effective and hands off to Chapter 7 correctly.

## 70-overall.tex

- **No chapter preamble** — straight from title into §7.1. One short paragraph stating the Discussion's job (synthesize the chain, then scope what remains) would fix it.
- **L6** — Database statistics re-quoted with different numbers and an unexplained different denominator: "of the 943 papers in the audited corpus, 833 carry structured curation notes, 95 carry an afferent classification, and 363 carry notes translating each pathway toward robot or simulation use." Chapter 5 reports, over its 561 **retained** papers: 523 notes, 86 afferent-classified, 258 translation notes (the 943/833 figures match Chapter 5's archived-inclusive corpus, but 95 and 363 match neither set). A committee member reading both chapters sees two different databases. Fix: quote one basis, matching Chapter 5 exactly, and say which corpus (retained vs audited-with-archived) each count uses.
- **L20** — "the SNS-Toolbox coupling runs at a **0.1 ms** neural step with force feedback" — 60-sns L73 says the coupling runs at a **2 ms** neural step matching the physics timestep (0.5 ms only as spiking sub-steps, L119). 0.1 ms appears nowhere in Chapter 6. Fix: correct to 2 ms (or whatever the delivered value is).
- **L20** — "The remaining blocker is non-finite leg-joint states in full-body co-simulation ... once it is resolved, the open work shifts from assembly to tuning and validation" — Chapter 6 reports four **tuned** production configurations with full-body ground evaluations and scores. The Discussion describes a pre-tuning state of the simulation work that Chapter 6 has moved past. This is the largest state mismatch in the document (together with 60-sns L73). Fix: rewrite §7.2.1's opening to the current state (assembly done; tuning begun with four production configurations; remaining work = closed-loop refinement, balance layers, bilateral gait).
- **L22–24** — §7.2.1 re-narrates the completed conversion in the future tense ("The model is kept in OpenSim 4.x format and converted to MuJoCo using MyoConverter ... The conversion is checked by inspecting the generated model ...") — duplicating 60-sns §6.2.2, which describes the same steps as complete and audited. Future-work section describing finished work confuses the reader about what remains. Fix: replace with a one-sentence pointer to §6.2 and keep only the genuinely open parts.
- **L26–28** — Fourth retelling of the prior AnimatLab walker (after ch6 preamble L12, ch6 §6.1, and L12 here), and it introduces a new number — "stepping-like motion at approximately 4 Hz" — that appears nowhere in Chapter 6's own account of that walker. Fix: keep the full story in §6.1 only; compress here to two sentences; either add the 4 Hz figure to §6.1 or drop it.
- **L41** — "The first of these tests is complete at the rhythm-generator level" — consistent with ch6; good.
- **L45–51** — Sensorium/hardware/kickflip closing is well-pitched and lands the arc.

## 80-conclusion.tex

- **L6** — "the analysis in Chapters~\ref{ch:methods} and~\ref{ch:results} produced novel equations for calculating isometric force in Festo BPAs" — renders "Chapters 3 and 4," but the actuator-level force equations are Chapter 3 only; Chapter 4 is joint-level. The alias labels (`ch:methods`, `ch:results`) survive the old structure and now mis-scope the claim. Fix: reference `ch:muscles` alone for the actuator-level sentence.
- **L10** — A single ~350-word paragraph covering methodological contribution, full future-work scope, cross-engine verification numbers, honest negatives, the ethics payoff, and whole-body maneuvers. It also duplicates Chapter 7's future-work specifics (embedded stack, tonic-stimulation-first verification) nearly verbatim. Fix: split into 2–3 paragraphs and reduce the future-work recap to one sentence pointing at Chapter 7.

## 92-AppendixA.tex

- **L4** — Good opener; scope matches the Intro's appendix description.
- **L224** — "As executed by the evaluator of record" — the of-record idiom is denser here than the appendix needs; acceptable as house provenance style, but two per page (L217 "two-bracket fit of record," L224) can be halved.
- No orphaned references found; every section traces to a main-text anchor. The GoF definitions (L176–177) are correctly back-referenced from Chapter 3 (20-muscles L119).

## 93-AppendixB.tex

- **L70** — "a **406-neuron** network" — Chapter 6 (L73) and this appendix's own L198 say **410** neurons. Stale count. Fix: 410 (and say once how the count is composed: 201 per leg + 8 shared).
- **L74** — "The campaign's artifacts live under \texttt{Code/MuJoCo_SNS/spinal/reports_spiking_20261002/}" — campaign word plus a date-bearing repo path; appendix provenance tolerates paths, but the date adds nothing. Trim.
- **L165** — "Connexions" — typo for "Connections" (table caption).
- **L194** — The document's only prose em-dash pair ("graphical verification --- or an equivalent acceptance ruling --- remains outstanding"). Fine as is; parenthesize if style demands.
- **L198** — "(\texttt{fitted\_walk\_params.json}, the IK/NNLS back-solve refit ...)" — the IK/NNLS back-solve is never introduced in any main chapter; appendix-only shorthand for a method the reader has not met. One clause of explanation (or a main-text sentence in ch6's tuning paragraph) would close it.

## 94-AppendixC.tex

- Generally well-organized provenance appendix; the stiffness-array and bounds tables do real work for Chapter 4.
- **L236** — The `removedText` paragraph: eight lines of forensics explaining that a `-p8` string in a results file is "a plotting leftover rather than a route statement," including the line numbers of the plotting block that produced it. This is workspace debugging, not dissertation content. Fix: cut to one sentence ("the results file's `removedText` field is a plotting artifact; the route of record includes p8 as scheduled in Table C.x") or move to a code comment.
- **"of record" ×11** (L34, 126, 157, 192 ×3, 195, 209 ×2, 236 ×2) — as appendix provenance style this is defensible, but "run of record," "calls of record," "fit of record" in three consecutive sentences (L209) reads as boilerplate. Halve the density.
- **L192, L209** — Date-bearing log/mat names (`frommat_flx_20260920.log`, `reeval_gof_20260926.m`). If main-text dates are stripped by the overhaul, decide one appendix policy (filenames yes, bare dates no) and apply it here too.
- **L173, L192** — "flexor front row 77" arithmetic mirrors 30-torque L425; if the main text is de-campaigned, re-point these to the adopted values table rather than Pareto-row indices.

---

## Cross-chapter findings

### Duplications (same content explained twice)

1. **Mono/Bi-pulley results twice in Chapter 4** — §4.5.1 (L440–442) vs §4.6.2 (L488): same numbers, same five sites, same sentinel discussion. Keep the Discussion version.
2. **Prior AnimatLab walker told four times** — ch6 preamble (L12), ch6 §6.1 (full), ch7 L12, ch7 L26–28. Keep §6.1; compress the other three to pointers.
3. **MyoConverter conversion described twice in different tenses** — ch6 §6.2.2 (complete, audited) vs ch7 L22–24 (future tense). The reader cannot tell whether it is done.
4. **Database statistics quoted twice with different numbers** — ch5 L24 vs ch7 L6 (different denominators, two values matching neither set).
5. **Conclusion recapitulates Chapter 7's future work** — 80-conclusion L10 vs 70-overall L43–51 (embedded stack, verification order, vestibular/ocular/cerebellar sequence).
6. *(Acceptable)* ch3 opener ≈ ch1 Chapter-3 roadmap — standard dissertation practice; no action.

### Handoff breaks

1. **ch5 → ch6 is one-way.** Chapter 5 promises its records inform Chapter 6's sensory channels; Chapter 6 never mentions the database (0 hits). The gap the Introduction motivates in objective 3 is never closed in the text. Highest-value single fix in the document after the ch6/ch7 state mismatch.
2. **ch4 → ch6 model provenance.** L427 claims the redesigns "are the model definitions carried into the neuromechanical simulation work"; L460 admits the converted tree predates them (stale `add_mag3`). Reconcile in one sentence.
3. **ch6 → ch7 state mismatch.** Chapter 6 reports tuned, ground-evaluated production configurations; Chapter 7 (L20) still describes a NaN blocker that precedes all tuning, a 0.1 ms step that matches nothing, and "the open work shifts from assembly to tuning." One of the two chapters describes an outdated pipeline.
4. **Background → ch6.** 15-background L256 claims the full loop ran "faster than real time"; Chapter 6 L73 says full-body co-simulation is not yet running (itself stale, per above). The Background claim should be deleted regardless.

### Numeric/terminology inconsistencies

| Where | Conflict |
|---|---|
| 70 L20 vs 60 L73 | neural step 0.1 ms vs 2 ms |
| 93 L70 vs 60 L73 / 93 L198 | 406 vs 410 neurons |
| 70 L6 vs 50 L24 | 943/833/95/363 vs 561/523/86/258 (different bases, two values match neither) |
| 20 L295 vs Table `tab:Fstar` | RMSE 0.004537 vs 0.005118; 10.6% is the validation max error, not model |
| 30 L394/L400 vs Table `tab:TorqueFit` L416 | 20 mm flexor 42.0 cm vs 41.5 cm |
| 30 L425 vs L437 | χ0 = +3.94 vs +3.937 mm (cosmetic rounding) |
| ch1 L20 vs ch2 L238 | Muscle Mutt vs DoggyDeux |
| ch1 L18 vs passim | SNN (once) vs SNS (everywhere) |
| 30 L394 et seq. | "Hunt model" / "original model" / "improved BPA characterization" used interchangeably |

### Orphaned references (things that no longer exist as described)

1. 30-torque L82 — `Chapter~\ref{ch:results}` resolves to the chapter it sits in.
2. 80-conclusion L6 — `ch:methods`/`ch:results` aliases mis-scope the actuator-level claim to two chapters.
3. 20-muscles L1 — `ch:discussion` alias referenced nowhere.
4. 15-background L131 — placement optimization "discussed in Chapter~\ref{ch:torque}"; Chapter 4 cites it once in passing (L50); the actual discussion is ch7 §7.1.
5. 15-background L278 / 70 L20 — "before any locomotion tuning is attempted" / "shifts from assembly to tuning" describe a state Chapter 6 no longer matches.

---

## TOP 10 fixes, ranked by reading-experience improvement before 2026-10-21

1. **Reconcile Chapter 6 with Chapter 7 on the state of the simulation work** (60-sns L73 "not yet running" vs §6.4.4+ tuned ground evaluations; 70 L20 NaN blocker + "assembly to tuning" + the 0.1 ms step). A committee member reading both will ask which document is true; nothing else in the review undermines confidence as much.
2. **Close the ch5→ch6 handoff** — add the Sensory Afferent Database paragraph to Chapter 6 (afferent channels trace to SADb records; cite `bolen_sensory_2023`). One paragraph converts Chapter 5 from an isolated study into the controller chapter's sensory foundation, which is the dissertation's declared logic.
3. **De-campaign §4.5 of Chapter 4** (30-torque L421–462): strip dates, `.mat`/`.m` filenames, "design of record," "front row 77," "angle index 77," collapse the three-margin chronology to the adopted cam result. (Concurrent overhaul target — this review confirms scope and adds the filename list.)
4. **Delete the in-chapter pulley duplication** — keep 30-torque §4.6.2, reduce §4.5.1's pulley/bi-pulley paragraph (L440–442) to a pointer. Saves a full page of repeated numbers.
5. **Add a chapter-opening paragraph to 30-torque** — the only main chapter without one, and the longest; 6–10 lines from actuator model to route redesign, matching the Intro roadmap (and updating that roadmap if the overhaul moves §4.5 to Discussion).
6. **Fix the Chapter 7 database statistics** (70 L6) to match Chapter 5 exactly and name the corpus basis — or quote only Chapter 5's retained-corpus numbers.
7. **Resolve or cut the dangling damping derivation at the end of 30-torque §4.6.1 (L479–484)** — currently an unintroduced three-line equation block ending the error-sources discussion.
8. **Split the Conclusion's 350-word final paragraph and fix its chapter reference** (80 L10; L6 `ch:methods`+`ch:results` → Chapter 3 only).
9. **Patch the concrete mismatches:** 20 L295 RMSE/max-error vs Table Fstar; 30 Table `tab:TorqueFit` 41.5→42.0 cm; 93 L70 406→410 neurons; 15 L256 "faster than real time" claim; typos 30 L33 ("in with"), L354 ("offet"), 93 L165 ("Connexions"), 20 L330 (doubled "also").
10. **De-duplicate the prior-walker narrative and the conversion re-description** — keep the full story in 60-sns §6.1; compress 60 L12, 70 L12, and 70 L26–28 to pointers; delete 70 L22–24's future-tense re-narration of the completed MyoConverter conversion (and its orphaned 4 Hz figure).

*Checked and clean, for the record: em-dash overuse (zero prose em-dashes in all main chapters; one pair in Appendix B) and "plant" used for the robot/model (zero occurrences — only "transplanting" and "plantarflexion").*
