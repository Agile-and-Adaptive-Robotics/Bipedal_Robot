# SUPERVISOR VERDICT - batch 20261026_deng-original-replication (easteregg2_2)

Reviewed 2026-10-26 on easteregg2 (hostname checked). Scope: the two insertion
fragments (patches/P01_deng_original_replication.tex,
patches/P02_bib_deng_biarticular_2023.tex), the proposed caption/alt text, and
REVIEW.md, against the four gates of
Notes/revision_campaign_20261015/supervisor-rubric.md under the 2026-10-22
style-contract supersessions (first person allowed; self-plagiarism gating
removed; own-paper citation still required - no Ben paper exists for this
content, it is his new simulation result). Pure insertions; anchors confirmed
at ProofFinal/chapters/60-sns.tex:202 (paragraph ending "...reproducing the
published behavior of the original model in its original engine.") and :204
("The replication series has a defined extension, ...").

---

GATE 1 (plagiarism / self-plagiarism): PASS
- External 7-gram overlap of P01, caption, and alt text against the Deng et al.
  2023 full text (D:\GitHub\_zcode_backups\lit_20261026\kaiyu2023_fulltext.txt):
  0 hits in all three. 10-gram overlap against the existing
  ProofFinal/chapters/*.tex: 0 hits.
- The paper's predecessor sentence ("perturbing the pattern formation network of
  a single joint did not affect the neuron activations of the rhythm generator")
  is paraphrased in P01 with citation, not quoted; no verbatim run.
- Ben's two papers are not implicated by subject; per the batch gate note no
  own-paper citation is due. The publication notes (20-muscles.tex:5,
  30-torque.tex:5) are untouched by this batch and unaffected.
- P02's bibliographic facts check against the paper header: authors
  "Kaiyu Deng, Alexander J. Hunt, Hillel J. Chiel, Roger D. Quinn" and title
  "Biarticular muscles improve the stability of a neuromechanical model of the
  rat hindlimb" match the fulltext; volume/pages/DOI intentionally blank for
  Ben's Zotero record (do not invent them).

GATE 2 (AI-slop): PASS
- Tier A regexes: 0 hits in P01/caption/alt. Tier B: 0 hits. No
  "not only ... but also". Em-dash count (---): 0 in all three. En-dash ranges
  (2.0--2.5 s, \SIrange macros) correct and exempt.
- No reviewer/comment text or meta-commentary in the manuscript prose; P01's
  leading % lines are patch scaffolding outside the fragment body.
- First person ("I exported", "I then applied") is allowed by the 2026-10-22
  supersession (style-contract.md section 2). Paragraph rhythm: 8 sentences,
  varied openings and lengths (4/7 "The/This + noun" - below the uniformity
  thresholds); register matches the surrounding section (60-sns.tex:202).
- Tense discipline correct: past for procedures, present for model behavior.

GATE 3 (factual drift): PASS on numbers (all traced); FLAG on two claim-level
items (F1, F2 below) that need dispositions before application.

Numbers table (every numeric token in P01 + caption + alt):

| Value | Where quoted | Source verified | Check |
|---|---|---|---|
| 10.1 s window | P01 | DENG_ORIGINAL_20261026_REPORT.md ("All runs 10.1 s") | OK |
| ModelDB 264844 | P01 | REPORT; REVIEW.md (author-confirmed mapping, Ben's Order B1) | OK |
| 0.23 s period | P01 | summary_stats_20261026.json B0: LH_RG_E 0.2331, LH_RG_F 0.2316, RH_RG_E 0.2268, RH_RG_F 0.2257 s | OK (all round to 0.23) |
| two legs in phase | P01 | summary_stats_20261026.json B0 L_vs_R_RGflx_lag 0.0073 s | OK |
| 2 nA, 2.0-2.1 s | P01 | paper fulltext ("a magnitude of 2 nA to the flexor neuron of the rhythm generator from 2 s to 2.1 s"); REPORT P1b/P2b | OK |
| extra burst 2.000 s | P01, cap (A), alt | run_P1b_rgflexor_2p1/analysis.json LH_RG_F onsets include 2.0 | OK |
| baseline 2.433 s | P01, cap (A), alt | run_B0_baseline/analysis.json LH_RG_F onsets include 2.433 | OK |
| early extensor 2.124 s / 2.589 s | cap (B), alt | P1b LH_RG_E onset 2.124; B0 LH_RG_E onset 2.589 | OK |
| post 0.224 s vs 0.221 s | P01 | REPORT table (P1b post 0.2238 s vs B0 0.2214 s); recomputed from onsets (P1b post diffs settle 0.222-0.223; B0 0.210-0.231) | OK |
| -2.3 +/- 2.4 percent | P01, cap (C) | phase_table_20261026.txt:34 (P2b LH_RG_E post -2.3, sd 2.4); definition in compare_phase_20261026.py (offset vs nearest baseline onset, wrapped to +/- half period) | OK |
| in-window mean 22.3 mV below baseline | cap (D) | summary_stats_20261026.json: P2 LH_KneeMN_E in-window -86.74 mV vs B0 -64.45 mV = 22.29 mV | OK |
| clipped at 0 mV; +current depolarizing | cap, alt | v_max = 0.0 in every column; REPORT sign check (P1 F -58.26 vs -60.79 mV) | OK |
| spans 1.6-3.6 s, -100..0 mV, band 2.0-2.5, dashed 2.1, blue/red | alt/cap | make_figure_review_20261026.py (xlim 1.6-3.6, ylim -105..5, ticks -100..0, axvspan 2.0-2.5, axvline 2.1, #4477AA/#EE6677, Arial, margin letters) | OK |

Method note: the figure was verified from its generator script and the saved
chart/analysis files; no pixel-level visual inspection was performed in this
session. The paper's 15.44/15 and 5 percent numbers are correctly NOT quoted in
P01 (REVIEW.md so claims; confirmed against the fulltext).

GATE 3 FLAGS (claims, not numbers - dispositions required):

- F1 - attribution mis-scope, P01 line 7: "The original produces the response
  that paper reports for its antagonist-only predecessor: an additional short
  flexor burst begins at stimulus onset ...". The 2023 paper reports the listed
  response (rapid short step, rhythm continuation) for its own biarticular
  model (Fig. 4); its only antagonist-only-predecessor perturbation statement
  is the pattern-formation-level no-effect claim, which P01's NEXT sentence
  cites correctly ("consistent with the antagonist-only predecessor behavior
  reported there"). The evidence report keeps the two apart ("reproduces the
  paper's PERTURBATION FAMILY ... the paper's own account of this 2019-era
  model is confirmed"; "Any claim of 'reproducing Deng 2023's Figure 4/5' must
  be scoped to: original-project response matches the paper's description of
  its own perturbation protocol family") - the patch sentence, as parsed most
  naturally, collapses them. Disposition: reword so the Fig-4-family response
  is attributed to the paper's protocol (e.g., "produces the rapid-short-step
  response family that paper reports under this stimulation protocol"),
  keeping "antagonist-only predecessor" only where the paper actually reports
  predecessor behavior; or Ben accepts the phrasing explicitly.
- F2 - ports clause vs the banked twin run, P01 line 7 (closing sentence): "the
  program's own ports of the same organization have not yet been given the
  published stimulation protocol." Literally true only under the paper-exact
  reading (2.0-2.1 s). The Simulink twin - a port of the same organization -
  received 2 nA on the RG flexor over 2.0-2.5 s on 2026-10-26
  (Code/Matlab/SNS_Simscape/dev/run_deng_perturbation_20261026.m:1-3, header:
  "Order B2: Deng 2023 perturbation protocol on the Simulink twin. 2 nA to the
  RG FLEXOR, t = 2.0-2.5 s"), and that script itself calls the 2.0-2.5 s
  stimulus "Deng 2023 perturbation protocol". The program's own artifacts
  therefore disagree about what "the published stimulation protocol" means,
  and REVIEW.md does not disposition this tension (its twin discussion
  concerns suppression-vs-phase-delay, not this clause). Needs Ben's ruling,
  or tightened wording (e.g., "have not yet been given the paper-exact
  stimulus window of 2.0-2.1 s").

GATE 4 (comment-intent): PARTIALLY ADDRESSED

Comment (verbatim): "Rybak's identified implementation is Python. Scope each
published-model claim to its actual implementation and demonstrated behavior."
(REVIEW_INBOX/README.md:41, "Current editorial rulings"; REVIEW.md states no
open advisor thread names this passage - the addition is justified by Ben's
Order B1 newly checked result.)

Claimed edits: P01 inserted between 60-sns.tex:202 and :204; P02 bib insertion
after deng_biomechanical_2022 (thesis.bib:1750); proposed new figure
fig:deng_original_perturbation with caption/alt.

Judgment: scoping is mostly exemplary - original-engine-only closing sentence;
period-level recovery numbers only; the 2023 quantitative phase magnitudes
explicitly disclaimed to the biarticular model ("a different network from the
archived project"); run identity stated so the 0.23 s period cannot be
conflated with the ports' 1.938 s at 60-sns.tex:181. Two clauses fall short of
the ruling's letter: F1 (the Fig-4-family response attributed to the
predecessor) and F2 (ports clause contradicted-at-the-edge by the B2 twin
run's own header).

Verdict: PARTIALLY ADDRESSED - remainder: reword F1; resolve F2 (Ben's ruling
or tightened wording).

MINOR FLAGS (mechanical; fix at application, no ruling needed):
- M1: caption.txt uses bare "2.000 s", "-2.3 +/- 2.4 percent" (percent as a
  prose word), "22.3 mV". At application, put units through siunitx (\qty) per
  style-contract section 5.
- M2: caption's global sentence "the dashed line marks 2.1 s, the end of the
  ... stimulus in the protocol" is untrue for panel (D) (0.5 s hold, stated
  per-panel). Add "paper-exact" to pre-empt the misread.
- M3: P02 header says "alphabetical position: after the deng_biomechanical_2022
  entry" - 'biarticular' sorts BEFORE 'biomechanical' (b-i-a < b-i-o), and the
  bib is not strictly alphabetical anyway (deng_neuromechanical_2019 at
  thesis.bib:1653 sits about 100 lines from deng_biomechanical_2022 at :1750).
  Say "near the other Deng entries" instead.
- M4: P02 title in Title Case vs house sentence case and the paper's own title
  casing; harmless to the bst render but inconsistent with neighbors (cf.
  deng_neuromechanical_2019).
- M5: no figure environment/asset-placement fragment is included (label and
  placement are suggestions; the asset must be copied into the Figures tree
  and P01 needs a Figure~\ref{} if the figure is referenced). Disclosed in
  REVIEW.md ("The figure is proposed, not placed"), so dispositioned - carry
  it as an application-time step.

OVERALL: BLOCK

No gate failures; the number trail is complete and the slop/overlap scans are
clean. Blocked on two undispositioned claim-level flags: F1 (reword the
predecessor-attribution clause in P01) and F2 (ports clause vs the banked B2
twin run). F2 is a scope question I do not resolve.

Escalations for Ben:
1. F2 ruling: is "the program's own ports ... have not yet been given the
   published stimulation protocol" acceptable as written (paper-exact reading),
   given run_deng_perturbation_20261026.m applied 2 nA to the twin's RG flexor
   over 2.0-2.5 s and its own header calls that "Deng 2023 perturbation
   protocol"? Or should P01 name the exact window / acknowledge the twin run?
2. F1 approval of the reword (or explicit acceptance of the current phrasing).
3. P02: supply volume/pages/DOI from Zotero before application (by design).

- diss-supervisor, 2026-10-26

=======================================================================
ROUND 2 (REV 2 fragments, reviewed 2026-10-26)
=======================================================================

Scope: revised P01 (F1/F2 fixes), P02 (M3/M4), new OPTIONAL P03 (M5), revised
caption.txt (M1/M2), REVIEW.md revision history. Round-1 number table remains
valid: the revised fragments introduce no new numeric values except
\qty{0.5}{\second} (the extended window 2.0-2.5 s, = the figure generator's
axvspan(2.0, 2.5)); all other tokens unchanged from round 1 and already
traced.

GATE 1 (ROUND 2): PASS
- 7-gram overlap of P01r2 / P02r2 / P03 / captionR2 vs the Deng-2023 fulltext:
  0 hits each. 10-gram vs existing chapters: 0 hits each.
- P02 title now sentence case, matching both the paper's own title and house
  style; authors/venue facts unchanged from the round-1 check.

GATE 2 (ROUND 2): PASS
- Tier A/B regexes over all four revised files: 0 hits. Em-dash count: 0.
  No reviewer-text leak. P01 rev 2 is 9 sentences with varied openings;
  register unchanged.

GATE 3 (ROUND 2): PASS (round-1 flags F1/F2 resolved; no new drift)
- F1 RESOLVED - verified against the paper fulltext. The revised sentence
  ("Under that stimulus the original shows the same phenomenology that paper
  reports for its biarticular model under the rhythm-generator stimulus: ...")
  attributes the rapid-short-step/rhythm-continuation phenomenology to the
  paper's Fig.-4 subject (its biarticular model), and the knee-PF sentence
  ("...which is the no-effect response that paper reports for its
  antagonist-only predecessor, the model class this archived project belongs
  to") attributes exactly what the paper states for the predecessor
  ("perturbing the pattern formation network of a single joint did not affect
  the neuron activations of the rhythm generator"). The magnitudes disclaimer
  is retained and strengthened ("so only the response family is claimed
  here"). Attribution now matches the source text clause for clause.
- F2 RESOLVED - verified against Code/Matlab/SNS_Simscape/dev/
  run_deng_perturbation_20261026.m:1-3 (2 nA to the RG flexor, 2.0-2.5 s).
  The revised closer ("the program's own port of the same organization has
  received this current only at an extended stimulus window, not at the
  paper-exact \SIrange{2}{2.1}{\second} timing") is exactly true of that run,
  makes no claim about the twin's response behavior, and no port is claimed to
  have received the paper-exact timing. The round-1 conflation risk is gone.
- M1/M2 RESOLVED: caption.txt is fully siunitx, and the dashed-line sentence
  is now correct for all panels ("the stimulus in panels (A)--(C) ends there,
  while the local knee-level stimulus in panel (D) is held for the full
  band"), consistent with the generator (P1b/P2b windows 2.0-2.1 s; P2 2.0-2.5
  s). P03's caption uses the shorter "the paper-exact stimulus end" form,
  which is also correct given panel (D)'s own full-window statement - see N1.
- M3/M4 RESOLVED: P02 anchor now "near the existing deng_* entries (for
  example immediately after deng_neuromechanical_2019)" with the correct
  non-alphabetical caveat; title sentence case.
- M5 RESOLVED: P03 supplies the full figure environment with siunitx caption,
  short caption, and label fig:deng_original_perturbation (no collision in the
  chapters). All numbers in P03 match the round-1 traced set.

GATE 4 (ROUND 2): ADDRESSED
- Comment (verbatim, unchanged): "Rybak's identified implementation is Python.
  Scope each published-model claim to its actual implementation and
  demonstrated behavior." (REVIEW_INBOX/README.md:41)
- Every published-model claim in P01 rev 2 is now scoped to an actual
  implementation (archived 2019 project, original Bullet engine, ModelDB
  264844) and demonstrated behavior (the logged runs): the Fig.-4-family
  phenomenology attributed to the paper's biarticular model, the no-effect
  knee-PF response attributed to the paper's predecessor statement, the
  quantitative magnitudes disclaimed to the biarticular model, and the ports
  statement matching the program's records (extended-window run acknowledged,
  paper-exact run correctly stated as not yet done).
- REVIEW.md's revision history records both round-1 flags and their fixes
  accurately; nothing was silently dropped; the two items left open for Ben
  are carried explicitly.

ROUND-2 NOTES (minor, editorial; dispositioned here, no ruling needed):
- N1: P03's caption omits caption.txt's clarifying parenthetical (panels A-C
  end at 2.1 s; panel D holds the full band). P03's "paper-exact stimulus
  end" is adequate since panel (D) states its full 0.5 s hold, but harmonize
  the two versions at application (prefer the explicit parenthetical).
- N2: P03 uses \begin{figure}[t] and width=0.72\textwidth; the chapter's
  existing figures use [htbp] with \textwidth-scale widths. P03's own comments
  flag this as adjustable - match the chapter convention at application.
- N3 (optional wording): "the program's own port" (singular) reads loosely
  given multiple Deng-lineage ports exist (60-sns.tex:437); "the program's
  Simulink port" or the plural would be tighter. The sentence is accurate as
  scoped either way; not blocking.

OVERALL (ROUND 2): GO

All four gates pass; round-1 flags F1, F2 and minors M1-M5 are fixed and
verified against their sources; round-2 notes N1-N3 are dispositioned above as
application-time editorial steps on the OPTIONAL P03 fragment and a wording
preference. Conditions carried into application (already open in REVIEW.md,
for Ben): (a) P02 volume/pages/DOI from Zotero before the bib entry is
applied - do not invent them; (b) the twin suppression-vs-phase-delay
reinterpretation remains a question for Ben/ChatGPT (REFERENCE_ONLY figure)
and no manuscript claim depends on it; (c) re-anchor P01/P03 against the live
Overleaf text if the anchor paragraphs have moved (REVIEW.md's standing
caveat); (d) copy the figure asset into the chapter's Figures tree and add a
Figure~\ref{fig:deng_original_perturbation} in P01 if the figure is accepted.

FINAL (ROUND 2): GO

- diss-supervisor, 2026-10-26 (round 2)
