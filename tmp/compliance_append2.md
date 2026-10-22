

## Supervisor review round (diss-supervisor agent, four gates) — fixes applied

The supervisor returned BLOCK with three factual failures and proofreading flags; all
mandated fixes applied and verified in the rebuild (285 pp, 0 undefined, 0 missing,
overfull multiset still 0 added vs. baseline):

| Finding | Fix applied |
|---|---|
| 60-sns opener: "graded neuron and synapse models … defined in Sections 2.5–2.8" false (only CPG + afferents defined there; equations live in App. B); also "is defined" with list subject; hardcoded numbers | Reworded: vocabulary/CPG/afferents → `Sections~\ref{sec:bg_neural}--\ref{sec:bg_sns}`; neuron and synapse equations → `Appendix~\ref{app:neuron_equations}`; "are defined" |
| 60-sns tuning ¶: "detailed in Section~\ref{sec:gait_library_results}" dangled (composition deleted) | Restored the baseline's verified composition sentence at the target subsection (32 running / 8 subjects / 4 speeds; ten walking 0.5–2.0 m/s + 3 self-selected; Falisse cycle) |
| 10-intro: "evaluated on the converted robot-body model" mislabel (Ch. 4 evaluates on the OpenSim model; conversion predates redesigns) | "converted" dropped |
| 10-intro: "every piece of it is hardware this dissertation builds on" overclaim (Muscle Mutt is a demonstration) | "and this dissertation builds on its hardware and demonstrations" |
| 10-intro: "spiking-neuron variants" vs ch. 6 "mirrors" term mismatch | "spiking-neuron mirrors"; closing list "the spiking comparison" |
| 50-database: "the afferent classes that controller monitors" | "that the controller monitors" |
| 60-sns:32 "documents the tools in build order" (advisor SNS-10 "huh?") | "documents the tools and what each contributes" |
| 30-torque:588 verbatim manuscript frames (Gate 1 FLAG: 11/18/12/10-word runs vs bolen_isometric_2025 Main.tex) | Prose frames reworded (opening sentence, adjoint-transform clause, load-cell clause, ceiling clause); numbers and data recitation untouched |
| 30-torque:588 "at the low end of the 269–507 kPa band" (256 < 269) | "just below the low end" |
| 30-torque:588 mixed `\newton\cdot\meter` / `\newton\meter` in one sentence | unified to `\newton\meter` (chapter-dominant form) |
| 30-torque:7 "follow the manuscript" understates the newer sweep | added "; the physical-validation paragraph of Section~\ref{sec:disc_pulley} additionally reports the newer nine-point sweep results, which postdate the manuscript draft" |
| INTRO-31 tracked edit (outcome-first Chapter-7 sentence) | advisor's wording applied verbatim at 10-introduction.tex:45 |

Escalated to Ben (NOT decided here): Gate-1 disposition is now moot (frames reworded),
but the flag is recorded; BG-24 "Delete?" ruling; SNS-14 depth-vs-reduce ruling;
first-person sweep scope (chs. 1–2 still "the author") + style-contract.md §2
supersession note; BG-03 and SNS-11 figure-level fixes.

## Final measurements (post-supervisor-fix)

| Chapter | Baseline narrative | Proposal narrative | Δ |
|---|---:|---:|---:|
| 10-introduction | 2,093 | 2,229 | +136 (repairs + INTRO-31) |
| 15-background | 5,217 | 5,392 | +175 (repairs; −1 figure) |
| 30-torque | 8,993 | 9,154 | +161 (verified sweep + reworded frames) |
| 50-database | 2,588 | 2,202 | −386 (−14.9%) |
| 60-sns | 10,851 | 7,910 | −2,941 (−27.1%) |
| Total (5 chapters) | 29,742 | 27,887 | −1,855 (−6.2%) |

PDF: 294 → 285 pages. Verification build of record: `tmp/build_live5.log`
(patches regenerated 2026-10-22 against the same fresh live baseline).
