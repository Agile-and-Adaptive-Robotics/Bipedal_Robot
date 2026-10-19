# EVIDENCE_SUMMARY — Bolen_Dissertation_Defense.pptx (whole-deck rollup)

**2026-10-18 A-E pass (after the EB475WS4 result fixes were pulled in):** slides 17, 25, 27,
33, 34, 38, 40, 41 were reworked to the assertion-evidence skill (call-outs instead of prose
cards, figures enlarged). Every number removed from those slide faces moved verbatim into the
speaker notes, so the manifests below remain the numbers of record. Slide 15's stand figure
was replaced with the chapter's real three-panel `testJigs_reviewed.pdf` composite (the old
asset carried the isometric force jig).

**Total slides: 54** (verified by python-pptx read of `Defense\Bolen_Dissertation_Defense.pptx`; renders
`deck_builder\deck_render\slide_01..54.png` post-date the pptx by 3–7 s, same build).

Manifest coverage (no gaps, no overlaps; matches `deck_ranges.json` and the pptx slide count):
`evidence_a.md` = slides 1–14 · `evidence_b.md` = 15–27 · `evidence_c.md` = 28–43 ·
`evidence_d.md` = 44–54 (rows lettered A–K).

Chapter numbering key (ProofFinal include order, verified against `build_r5/main.toc/.lot/.lof`):
10-introduction = Ch. 1 · 15-background = Ch. 2 · 20-muscles = Ch. 3 · 30-torque = Ch. 4 ·
50-database = Ch. 5 · 60-sns = Ch. 6 · 70-overall = Ch. 7 · 80-conclusion = Ch. 8 ·
appendices 92/93/94 = A/B/C. Chapter paths below are
`Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\`.

## Caveats carried from the manifests + final-gate findings

1. **Slides 35 and 36 carried off-by-one table cites (final-gate finding; FIXED 2026-10-18 after
   the run).** The cite lines had read "Table 6.2"/"Table 6.3" where the compiled dissertation
   numbers them **Table 6.1** (gait-validation sweep, `tab:gait_sweep`, 60-sns.tex:171–173) and
   **Table 6.2** (pathway removal, `tab:walker_ablations`, 60-sns.tex:195–197). Slide cites,
   speaker note, and evidence_c rows 35/36 were corrected and the deck rebuilt; the slide
   *numbers* (−144.5/−5808 … −9961; deletion scores) verified correct throughout. Slide 53's
   cites (Table 6.4, Figs. 6.5/6.11) verify correct against the same .lot/.lof.
2. **Slide 43 speaker notes used stale appendix letters (FIXED 2026-10-18 after the run).**
   Notes had routed "bracket geometry and adopted χ values" to "Appendix B" and "spiking
   calibration" to "Appendix A tables"; they now route to Appendices A–C (force model, GoF),
   H (bracket data, slide 51), F (protocol/provenance), and J (implementation, slide 53).
3. Slide 1's "November 2026" is a CONTENT_SPEC directive (defense 2026-11-06), not a
   chapter-sourced number. Verified present in the pptx text and in the rendered footer band.
4. Slide 19: the source figure's own panel titles read "48.0 cm" where the chapter caption/text
   say 48.5 cm (flagged in the slide's speaker notes). Slide 23: chapter text says l₀ = 42.0 cm
   for the ⌀20 mm test while Table TorqueFit lists 41.5 cm; the slide follows the text
   (evidence_b caveat).
5. Slide 34 carries the deck's only **since submission** content (spiking-twin 62.7°/61.9° and
   curriculum 266.9 vs 185.9), labeled on-slide and cited to artifact paths per CONTENT_SPEC.
6. Hardware extensor numbers on slide 26 use the chapter's 3-dp values (1.163/7.849/23.437 N·m,
   30-torque.tex:544), not the 4-dp artifact variants.
7. Page numbers: present and equal to position on slides 2–42 and 44–54; absent on the title
   (1) and the thank-you slide (43). Speaker notes: present on all 54; thinnest is slide 42
   (19 words). Kickers: chapter-based throughout; no "Aim 1/2/3", no "September 2026",
   no "Ch. 5 (Discussion)" anywhere in the deck (full-text scan of the pptx).
8. Final-gate spot re-verification (this session, against the chapter files): manifest rows
   a2, a12, b20, b21, c28, c29, c33, c39, dH — all numbers matched their cited lines exactly.
   Rows marked ✔ below were re-verified this session.

## Per-slide evidence table

| # | Assertion | Key number(s) | Source |
|---|-----------|---------------|--------|
| 1 | Title: dissertation defense, name, advisor, lab | November 2026 | CONTENT_SPEC s01 date directive (defense 2026-11-06); not chapter-sourced |
| 2 | Four objectives connect one actuator to a walking, neurally controlled robot ✔ | 4 cards; 537 tests; 943-paper corpus | 10-introduction.tex:26-31; 20-muscles.tex:290; 30-torque.tex:504; 50-database.tex:14; 60-sns.tex:6,12 |
| 3 | A humanoid robot from the trunk down | 92 muscle paths | 10-introduction.tex:18,20; 60-sns.tex:73 |
| 4 | Robots run neuromechanical experiments impossible in living subjects | — | 10-introduction.tex:6 |
| 5 | BPAs: muscle-like strengths, three modeling headaches | ε₆₂₀ | 10-introduction.tex:8; 20-muscles.tex:22-28,312 |
| 6 | FOUR gaps between existing design and neural walking | 4 gap cards | 10-introduction.tex:14 (+8,10,12,24) |
| 7 | Custom jig: isometric force across l₀, P, contraction | 10/20 mm DMSP; 112–518 / 300–509 mm; 620 kPa; 537 = 321+216 | 20-muscles.tex:71,74-78,290 |
| 8 | Max isometric force depends on resting length (unreported) | — | 20-muscles.tex:96 |
| 9 | Two coefficients per diameter capture the l₀ effect | 470 vs 490 N; 1460 vs 1570 N; 10%; 7.5 / 13 mm | 20-muscles.tex:98,99,57,103,310,312 |
| 10 | Festo tool mis-predicts short BPAs by 42% | 4.5% vs 42.1%; 350.9/335.3/498.6 N; l₀ = 112 mm | 20-muscles.tex:312,316 |
| 11 | All 537 measurements collapse onto one surface | ε* = ε/ε₆₂₀; P* = P/620 kPa | 20-muscles.tex:290,33,52,264 |
| 12 | Three coefficients per diameter; OOS to a few percent ✔ | 0.9998/0.9994/10.3%/10.6%; 0.992/0.9943/6.8%/5.7%; RMSE 0.005/0.023 | 20-muscles.tex:221-248,290 |
| 13 | Lowest error of compared models everywhere | RMSE grid 22.9…862.1 N | 20-muscles.tex:369-372,391 |
| 14 | Force model = two cheap maps for real-time control | — | 20-muscles.tex:399,62-64 |
| 15 | Two test stands: pinned + biomimetic four-bar | ⌀10/⌀20 mm; 1.83 m donor | 30-torque.tex:33,22,55,24 |
| 16 | Hybrid method: the error is kinematic, not force | M = r × F, measured Lm, rk | 30-torque.tex:63-78,338,346 |
| 17 | Three effects close the gap: offset, compliance, wrap loss | ζ = 1 − ε*; bracketed fzero | 30-torque.tex:88-93,97-99,141-151,211-230,278-291,246-252,126-132,254-270 |
| 18 | Terms identified by multiobjective optimization, held-out judged | 4 tests, 2 train / 2 hold per fold; 47.9 cm not used; bounds; pop 150/600, 50/150 | 30-torque.tex:301,305,373-377,319; 94-AppendixC.tex:18,175-186 |
| 19 | Corrected prediction tracks measurement incl. the dip | RMSE 1.616 / 2.1214 N·m | 30-torque.tex:346,348,469-470 |
| 20 | Flexor training-mean RMSE 8.66 → 1.62; validation 7.59 → 1.49 ✔ | per-test 2.165/1.084/1.368/1.611; FVU ≤ 0.18; 47.9 post hoc | 30-torque.tex:358-360,373-380,502; 94-AppendixC.tex:212-236 |
| 21 | Extensor wraps; χ₃ recovers lost length ✔ | χ₀ −6.38 mm; χ₃ 0.158; training 2.50 → 1.00, validation 2.55 → 1.10 N·m | 30-torque.tex:502,409-417,398,358,387,396 |
| 22 | One stiffness pair of record spans configurations ✔ | +3.94 mm; 3.998×10⁴; 1.473×10⁴ N/m; rows 77/32 | 94-AppendixC.tex:196-203,207,192; 30-torque.tex:502 |
| 23 | Biomimetic knee meets/exceeds human flexor torque | shortfall 60–100°; RMSE 1.80 / 0.30 N·m | 30-torque.tex:451,457,472-473 |
| 24 | Bio-extensor: corrected still above measurement | 51.8 cm; 0.89 / 1.17 N·m; FVU 0.67 | 30-torque.tex:430,471,502; 94-AppendixC.tex:209 |
| 25 | χ values drive a redesign ladder ✔ | +5.01% @ −107°; +5.02%; cam +30.83%/+39.77%; 85.74 N·m | 30-torque.tex:504,510,512,514,516,518,508,330,328 |
| 26 | Cam extensor = design of record | 269–507 kPa; 96.9 mm; 43.1 mm; 72.4 mm mean arm; 1.163/7.849/23.437 N·m | 30-torque.tex:512,57,544; fig Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_4.png |
| 27 | Pulleys buy reach; geometry binds | 81.23 → 85.74 N·m; ≤0.803 N·m; ×3.51/×3.70; 4.26–19.14 vs 25 mm | 30-torque.tex:542,516,527-531,518 |
| 28 | One audited, offline-searchable afferent resource ✔ | 943; 561; 382; 730/11,130/18; 957 = 561+39+244+109+4 | 50-database.tex:14,20,26,28,33 |
| 29 | Machine drafts at scale; judgments stay human ✔ | 374 / 15 batches; 164 flagged; 436 removed; 351 protected; 523/426/258… | 50-database.tex:14,18,20,26,41 |
| 30 | Pathway record → circuit edge | pathway counts; 39 records; 258 notes | 50-database.tex:45,54-58,71-78,83 (Tables 5.1–5.2; numbering verified in main.lot) |
| 31 | Direct ancestor: hand-tuned AnimatLab walker | 9 bodies; ≈4 Hz vs ≈1 Hz | 60-sns.tex:18,20,135; 70-overall.tex:26,28 |
| 32 | A toolchain that did not exist | 92 actuators; −35°; 7/12 blocks; 2.00° vs 7.85°; 410 neurons; 2 ms | 60-sns.tex:71,73,85,91,93 |
| 33 | Rhythm core in two engines: exact agreement ✔ | 1.000/1.000 s; r −0.525; r ≥ 0.99996; RMSE 0.099–0.111 mV; 11/11 | 60-sns.tex:113,146,148,159 |
| 34 | Characterized BPA model drives a MuJoCo walker | 1.542 s; 63.5°; 880/1372 N; 372 kPa; +0.001 m ground | 60-sns.tex:163,165; since-submission: bpa_spk_best_20261003.npz + CURRICULUM_RUN_20261003.log:136-138 |
| 35 | One-reference tuning transfers poorly | −144.5/−5808 … −200.4/−9961; 43 = 32+10+1 | 60-sns.tex:169,179-182,187 (Table 6.1) |
| 36 | First lesion and stimulus experiments | Ia −25; all-aff −267.8; Renshaw −267.3; II −255.2; sway 0.123 m | 60-sns.tex:193,202-218,224,191 (Table 6.2); connectome JSON 110n/216e |
| 37 | Spiking mirrors judged by gates | 0.399 vs 0.621 s (+55.6%); 1,186/7,908/7,762; 12/12; −236.35 vs −161.57 | 60-sns.tex:234,236,246,255-259,232 |
| 38 | Spiking changes timing, not the interface | 15.79°; 5.85°/5.83°; 2.00/7.85°; 442–450 Hz; 12/12 | 60-sns.tex:277,303-305,355,357,327,341-342,123,275,353,359 |
| 39 | Each stage feeds the next ✔ | 3 coefficients; meets/exceeds; 943; r > 0.9999; 60% ceiling | 70-overall.tex:10,12,14; 50-database.tex:14; 80-conclusion.tex:8 |
| 40 | What these results do not claim | ε₆₂₀ a priori; family of pairs; 4.49→0.67; 3.03→1.17 | 20-muscles.tex:318,397; 30-torque.tex:502,488; 60-sns.tex:73,165,187,222,363; 50-database.tex:87 |
| 41 | Program the tools define | staged loop order; Teensy/Jetson; kickflip | 70-overall.tex:20,24,45,47,49,51 |
| 42 | Acknowledgments | names only | unchanged from 20260910 deck per CONTENT_SPEC |
| 43 | Thank you — questions | appendix routing | unchanged; **notes carry stale appendix letters (A/B/C vs deck A–K)** |
| 44 | Appendix A: complete force model | Eqs. (3.1), (3.4)–(3.6) | 20-muscles.tex:17,39-50,54-58,60-64; 92-AppendixA.tex:190-202 |
| 45 | Appendix B: coefficients + 95% CIs | 303.5 (300,308); 922.4; c₀ 0.5682… | 20-muscles.tex:167-199,221-248; 92-AppendixA.tex:197 |
| 46 | Appendix C: GoF vs published models | full RMSE/FVU/max grid | 20-muscles.tex:369-372 |
| 47 | Appendix D: compliance chain | Eqs. (4.18), (4.20), (4.22); fzero bracket | 30-torque.tex:215-220,228-232,244-250; 92-AppendixA.tex:233-238 |
| 48 | Appendix E: wrap loss per arc | R₁ 12 mm 27°→−23°; R₂ 40 mm; Eqs. (4.24)–(4.26) | 30-torque.tex:279-297; 92-AppendixA.tex:242; 94-AppendixC.tex:173 |
| 49 | Appendix F: protocol + data decisions | +5.3° encoder; flexor 2/2 folds, 47.9 not used; extensor 3/3/3; rows 77/32; +5.01/+5.02 | 94-AppendixC.tex:18,72,176-188,192,201-207; 30-torque.tex:301-305,319,502-504 |
| 50 | Appendix G: frame invariance | iff k_x = k_z | 30-torque.tex:139,303; 94-AppendixC.tex:126,156-171 |
| 51 | Appendix H: bracket points + adopted terms ✔ | six points (Pbri…); χ +3.94/−6.38 mm; 3.998e4/1.473e4; 0.158 | 94-AppendixC.tex:128-143,194-205,207; 30-torque.tex:502 |
| 52 | Appendix I: moment-arm conventions | perpendicular foot; ε < −0.02 | 30-torque.tex:37-46,48,326; 92-AppendixA.tex:213,32-36 |
| 53 | Appendix J: SNS library + spiking twin | 7 → 12 → 15 blocks; 15.79°; 6.2×10⁻⁷ µS; 62/62 | 60-sns.tex:91,123,275,277,286-288,93; Figs. 6.5/6.11 (verified in main.lof) |
| 54 | Appendix K: walker numbers of record | sweep + deletion tables verbatim | 60-sns.tex:179-182,187,193,202-217,224 (Tables 6.1/6.2; numbering verified) |
