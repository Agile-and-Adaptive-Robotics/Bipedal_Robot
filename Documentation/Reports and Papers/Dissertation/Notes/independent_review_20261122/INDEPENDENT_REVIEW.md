# Independent dissertation text and figure review

**Disposition for Ben.** The live dissertation makes four identifiable contributions. Its strongest evidence connects measured actuator force to knee design and a specific sensory experiment to a model rule. The current text needs targeted corrections to the provenance of a fitted parameter set, the conditions under which the flexor met its target, and several source credits. Those corrections affect how a reader would reproduce or cite the work. They do not erase the measured knee-torque result. The neuromechanics results support controller construction, controlled interventions, and specific behaviors under external support. Two prominent claims about gait need the analysis window and phase statistic that made them true. A wholesale rewrite would lose effective explanations already present.

**Authority and independence.** The primary reading copy was the 353-page live Overleaf export on EB475WS4, read in printed order and checked against its extracted source ZIP. The finalized 353-page handoff copy and the reconstructed November 20 candidate were compared separately. Locations below use printed page numbers. Add 20 to an Arabic page number to find the PDF leaf. The title and copyright precede roman p. i. Source lines are relative to the live ZIP unless a different version is named. No manuscript or approved author wording was edited. This assessment was prepared before consulting a separate easteregg2 review. The companion version comparison records the file origins, hashes, and identifiable differences.

## Four objectives traced to their evidence

| Objective in §1.3, pp. 6–7 | Contribution and strongest evidence | What the conclusion can support |
|---|---|---|
| Predict isometric BPA force across actuator sizes and resting lengths | Chapter 3, pp. 42–67. The resting-length fit at 620 kPa and the normalized pressure–contraction surface address different parts of the prediction. Tables 3.2–3.3 give fit and comparison results. Appendix A, pp. 217–219, describes the implementation. | Chapter 8, p. 184, correctly confines the contribution to tested isometric conditions. Check the short-actuator worked example in the designated published chapter. |
| Predict joint torque and use it to redesign the knee | Chapter 4, pp. 68–122. Pinned-knee experiments identify length and compliance corrections. Biomimetic-knee curves test transfer, and the two redesigned routes are physically tested against the selected muscle torques. Figures 4.13–4.15, pp. 117–120, and Appendix C give the geometry and evidence. | At the measured angles and pressures both routes exceeded their respective reference torques. Keep the modeled full-angle design envelope, physical samples, and parameter-selection evidence distinct. |
| Make sensory-afferent evidence usable for model design | Chapter 5, pp. 123–135, and Appendix E, pp. 329–333. A reader can find the McCrea paper, inspect the preparation and perturbation, find both labels in the record, and use the conditional sign of the response to constrain a model rule. The five record categories account for all 1,393 records. | The application makes physiological evidence and its qualifications retrievable. Mixed group-I stimulation constrains a conditional response without individually identifying the strength of an Ib connection. The 943 Papers records are one category. |
| Build and test tools for connecting neural control to mechanical models | Chapter 6, pp. 136–167. Figures 6.2–6.8 connect model preparation, editable neural blocks, and controlled single-joint tests. Figures 6.9–6.12 and Appendices B and D separate circuit behavior, body response, and mechanical support. | Chapters 7–8 identify original tools, controlled demonstrations, and remaining transfer and balance questions. The 0.94° knee extension followed an injected current. The pneumatic arm descended 2.00° with feedback and 7.85° with it disabled under the same 0.625 kg load. Neither result requires a walking claim. |

## Decisions with the largest effect on interpretation

| Priority and location in the live edition | What a reader could conclude | Evidence and smallest repair |
|---|---|---|
| 1. Chapter 4 and Appendix C, pinned-extensor training and validation | A third trial contributed to the fitted row and its reported training improvement. | The adopted row 32 in the reference-commit MAT stores five eligible indices and holds out three. Averaging the two remaining trials reproduces all three saved objectives. Adding test 7 fails that check. The two-test training mean RMSE drops from 2.658 to 1.251 N·m (52.9%). The three-test validation improvement remains 57.0%. Correct the training labels and the associated summary. A real prediction for test 7 can be reported with its actual status. See the mechanical evidence below for exact locations, identifiers, and the independent-selection qualification. |
| 2. Abstract p. i, §4.3.4 p. 116, Fig. 4.14 p. 118, conclusion p. 184 | Readers may carry the modeled 620 kPa ceiling into the measured full-range result. | All 15 sampled flexor poses met the human target at the recorded pressures. Five readings were 628–638 kPa. The +0.42 N·m margin at −105° occurred at 638 kPa. Retain the finalized abstract's “at measured configurations.” Add “at the measured pressures” in the result and give the over-ceiling count once. |
| 3. Chapter 2 Fig. 2.4 p. 22, bibliography pp. 186–206 | A reader following the credit cannot find the two-panel illustration. Several bibliography entries lead to different experiments, investigators, or organisms. | Both panels match Fig. 6 of Steele, Etoundi, and Hunt's later ASME paper (DOI 10.1115/1.4054441). The earlier conference figure supplies only the mechanism view. Update the caption citation and retain the appropriate credit. Correct the high-confidence reference identities listed below using publishers' records and original PDFs. |
| 4. Chapter 6 §6.2.1 p. 142, Appendix D §D.3.5 pp. 302–305 | A reader may take stance within three percentage points of the reference as a steady-drive walking result. | The 62.3% and 61.4% estimates include the transition toward standing. With drive held constant the values were 77.5% and 75.4%, compared with the 61% reference. State the window alongside the first favorable comparison. |
| 5. Chapters 6 and 8, pp. 160–162 and 185, Appendix D pp. 309–311 | “Alternating heel strikes” can imply each left contact was followed by one right contact. | The intact record counts 34 left and 22 right strikes. Median left phase within a right cycle was 0.498 intact, 0.494 after Ib-only removal, and 0.223 after combined removal. Say which median phase was preserved or shifted. The three-second traces in Fig. 6.11 are illustrations. Give the scored interval and affected pathways next to the comparison. |
| 6. Chapter 3 §3.3.1 p. 58, 20-muscles.tex:349 | A reader substituting the inputs cannot reproduce the reported 112 N or 8% error. | Direct evaluation of the printed equations at 54 mm and 300 kPa gives 106.38 N and 3.28% against 103 N. The published text uses intermediate values 0.49 and 230 N. Flag this for Ben and the coauthor as a possible author correction or dissertation note. Preserve the published chapter. |

## Additional checks at the point of use

- **Modeled target and optimization, pp. 94–96, 214.** At 30-torque.tex:387,414, “no shortfall at a single angle can be traded away for a better average” does not follow from positive finite weights of 100,000 and 1,000. If a separate constraint enforces pointwise coverage, name that constraint. Otherwise say the weights give priority to reducing the worst sampled shortfall, and report the selected design's actual minimum margin. At 30-torque.tex:457, “maximum isometric force” can be read as the constant muscle parameter. At 92-AppendixA.tex:161 the reference is the isometric tensile force at each pose. Use the Appendix A definition at first use and retain the names of the exported curves.
- **Effect of the muscle equations, p. 141.** At 60-sns.tex:81, specify what “mattered most” means. Rectus femoris and biceps femoris short head showed the largest gains in R² (0.2783 and 0.1706). Their absolute reductions in force RMSE were 24.23 and 37.45 N. Soleus and medial gastrocnemius had larger absolute reductions of 142.76 and 110.0 N (Table D2, p. 272). A two-word metric qualifier resolves the apparent contradiction.
- **Ethics, pp. 38–41.** The section explains why reversible interventions on a validated model can answer questions that are ethically difficult to investigate in living subjects. Preserve Ben's position and wording. At 15-background.tex:290 the reference for the “loop” goes to §1.2, the robot project. The bidirectional loop is introduced in §1.1. Retarget the reference. Treat reduction of future experiments as the scientific opportunity the author describes, conditional on model validation for the specific question.
- **Dedication, p. ii.** Quotation marks already identify borrowed wording. The phrase “Los corazones obreros” at 04-dedication.tex:20–23 matches the same verse of *Hijos del pueblo*. A short title credit below the stanza would let readers identify it while preserving the author's dedication. The documentary account discusses historical uncertainty about the personal attribution. [Published source comparison](https://www.penguinlibros.com/es/historia/11502-libro-que-me-estas-cantando-9788499921341/fragmento)
- **Physical clearance, Appendix C, p. 257.** At 94-AppendixC.tex:377 contraction is the stated reason for relaxing the clearance by 3 mm at release. Axial contraction of a braided actuator generally increases its diameter. Pressure, local bending, and partial contraction can still give a diameter below the assumed maximum. Cite the diameter at those release poses or identify the relaxation as a numerical tolerance and state where physical clearance was checked. The passage alone cannot establish an interference or invalidate the routing result.

## Wording and figures worth preserving

- The abstract on p. i identifies four contributions in one paragraph and qualifies the knee result by measured configurations. The acknowledgments on p. iii retain individual and collaborative credit. Both were read in full. No wording change to either is proposed. The detailed conditions on force and locomotion belong with the evidence and are noted above.
- Chapter 5 makes the McCrea example earn its space. It gives the reader a paper, an actual app operation, the original evidence, its limitation, and an engineering decision. Preserve the distinction between mixed group-I stimulation and isolated Ib stimulation. The database's citation map complements the evidence without substituting citation frequency for physiological support.
- Preserve the side-by-side knee and arm experiments (Figs. 6.6–6.7, pp. 150–151). The captions and axes distinguish injected neural current, muscle force, a physical load, arm displacement, and commanded pressure. The knee plots align at pulse onset and use different time units for onset and recovery. Keep the mechanical stretch test referenced separately.
- Preserve the model identifiers and support conditions next to Fig. 6.10 (p. 161) and Fig. 6.12 (p. 164). The former plots a feet-off-ground human-model trial aligned to maximum right-knee flexion. The stance bars belong to the walking references. The latter depicts the CAD robot with 54 ideal tensile actuators under substantial pelvic support. The pressure-driven robot counterpart is a third experiment. Table 6.4 on p. 164 reports execution time and separates the entire run from neural updates.
- Figure 7.2, p. 172, distinguishes fittings from the end effect in the adjacent actuator wall. Its three loading conditions make the proposed test concrete. Preserve the distinction between a proposed measurement and a measured mechanism.

## Evidence for the mechanical findings

### Accuracy and interpretation

1. **The archived extensor fit used two training tests. A third test is mislabeled as training.** Live pp. 88–90, 105–107, 112, 240, 242, 250–251. `30-torque.tex` lines 316, 565, 569–595, 657; `94-AppendixC.tex` lines 79, 94, 267, 304. The words “trains on the 43.6 cm, 43.2 cm-tendon, and 45.7 cm tests,” the panel C “Training” tag, and Tables 4.7, C4, and C12 classify test 7 as a fit. The specified source, `Testing_Data/2022_02_Festo/minimizeExt10mmX3_results.mat`, has `allBPA=[1,2,5,6,8]`, `numBPA=5`, and all ten three-test fold combinations of those five. Selected row 32 holds out {1,2,8}. Its three saved normalized training objectives are (0.46542348, 0.22521556, 0.46972816). Recomputed from `f` and `a0`, these match the mean ratios for **{5,6}**. Including 7 produces (0.38521591, 0.16698883, 0.39396581). The stored prediction for test 7 is real and can be reported as a retrospective test. The two-test training mean RMSE is 2.657605 → 1.251093 N·m, a 52.92% reduction. Adding the retrospective case yields the published three-case 60.07%. The three held-out cases still give 2.554118 → 1.097665 N·m (57.02%). Relabel test 7 in Fig. 4.9 and both chapters, recompute the training mean row from cases 5 and 6, and update the 60% claim. The extensor MAT file is byte-identical to its Git blob `4e1cdfff8290eee1e0d749de93d9ec48066dc410` at the stated 477ae532 checkpoint (SHA-256 97e3d61e6aa19e5d13b7b49cb8951432bbcf97048551775a0629d8a043bfe420). If another run is intended, its primary archive and matching parameter/score row are required.

2. **Separate parameter fitting, candidate selection, and evidence of transfer.** Live pp. 88–91, 99, 112–113, 249–251. `30-torque.tex` lines 316, 342–344, 478, 657, 673; `94-AppendixC.tex` lines 263–267, 311. The account does say candidates were filtered, but “reduced torque-prediction error even on the pinned-knee tests held out of each fit” does not explain that those same observations helped screen candidates. Both drivers calculate held-out scores for every candidate, sort, and reject baseline failures before a member is selected (`minimizeFlxPin10mm_2brk.m` lines 116, 135, 168–197; `minimizeExt10mmX3.m` lines 124, 140, 174–211). Moreover, the extensor MAT’s `LOCKSRC` credits the locked stiffness choice in part to “bio-flexor score”; the documentary record on later external design requirements should be reconciled before describing the biomimetic flexor as independent. The unused 48 cm extensor test worsens on all metrics (C12). Keep the measured reductions. Add a single sentence explaining that held-out configurations were excluded from optimization and used to filter and select among its solutions. Reserve independent-transfer language for configurations shown to be absent from the entire selection process. Identify the retrospective result where this distinction matters.

3. **Report the actual flexor test pressures beside the 620 kPa design requirement.** Live pp. 113, 116, 118; `30-torque.tex` lines 677, 740–744, 756. “With margins from 0.42 N·m at −105°” is a discrete physical-test result. The extensor gives its pressure range, but the adjacent flexor paragraph and caption omit it. In `Testing_Data/2026_06_Festo/Results_table_20mm.xlsx`, sheet `FlxTest20mm_56cm`, C7:Q8, the 15 measured pressures range from 226.5 to 638 kPa. Five exceed 620 kPa: −119° (633), −114° (634), −105° (638), −87° (635), and −83° (628). The stored values agree exactly with the array plotted by `Figures/30-torque/Components/publication_20261028/publication_flexor_pressure.m` and its companion JSON. Keep the observed claim and add the range/count next to it. The +0.42 N·m minimum was observed at 638 kPa. Those five measurements alone do not demonstrate physical attainment at or below 620 kPa. The abstract’s “at measured configurations” is appropriately scoped. No pressure correction or rerun is necessary for this wording.

4. **Joint-torque reference requires pose-dependent muscle force.** Live p. 95, `30-torque.tex` line 457; compare p. 214, `92-AppendixA.tex` lines 160–165. “Moment arm is multiplied by the muscle’s maximum isometric force” can be implemented using Gait2392’s constant optimal-fiber force parameter. The archived `OpenSim_{Vasti,Bifem}_{Results,MomentArm}.txt` files show otherwise. Torque divided by moment arm for vastus medialis is 1,231.39 N at −120° and 643.15 N at +10°. Short-head biceps femoris gives 553.94 and 738.70 N respectively. Specify “the muscle’s isometric tensile force at that pose, as defined in Appendix A.” The definition in A is sound.

5. **The final paired flexor’s proximal attachments are described two ways.** Live pp. 73 and 99, `30-torque.tex` lines 45 and 470. The first occurrence says “pathways are mirrored about the Z axis”; the second explains the final staggered proximal pair. The figure generator calls `flexorBpa2Endpoints20mm`; lines 18–19 give p1Bz=p1Az−2p2Az and p2Bz=−p2Az. One origin does not generally reflect the other. Identify the distal attachments as reflected across the sagittal plane and the proximal attachments as staggered. If the early sentence describes separate physical hardware, name that configuration there and identify where its tests appear.

6. **Clarify the two generations of the flexor calculation.** Live pp. 109–110, 113, 116–118, 243–244; `30-torque.tex` lines 616, 677, 740, 756; `94-AppendixC.tex` lines 107–115, 135–159. The design optimization uses maximum usable contraction 0.255. The plotted final 620 kPa and pressure-matched predictions explicitly use 0.235 in `publication_flexor_pressure.m` (lines 67–71), consistent with the 0.2335 and 0.2357 measured fractions in worksheet cells I3:J3. Give those two model settings alongside the existing measured-versus-shared-length explanation in Table C7, or add a short parenthesis at first use of Fig. 4.14. In the single-BPA case, distinguish the physical 21 mm flexor tendon (p. 109) from the effective 11 mm used in C5. These need clear labels, not a parameter substitution. The worksheet’s raw tendon entries are 0.28/0.30 and 0.92/0.85, whereas C7 reports 28/30 mm and 92/85 mm. Confirm the scale of these four cells or the original caliper/lab record if the table is to be called a direct metric transcription; the bare decimals do not establish the units.

7. **Protected Chapter 3 contains a reproducible worked-example discrepancy.** Live p. 58, `20-muscles.tex` line 349, Eqs. (3.4), (3.6), (3.7), pp. 44–45, 49. The text reports “0.49 · 230 N = 112 N” and 8% error. At the stated 54 mm length and zero contraction, F620 = 303.5 atan[19.03(0.054−0.0075)] = 219.8575 N, F* = 300/620, and the prediction is 106.3826 N. Relative to 103 N the error is 3.28%. The conclusion that this predicts the measurement better than the 270 N manufacturer estimate survives. The published paper also prints the 112 N example on pp. 10–11 (Bolen et al., Actuators 15, 230, published 22 April 2026, DOI 10.3390/act15050230; original publisher PDF in local Zotero attachment F5GKFJ8Q, SHA-256 7d148b0f3793323ccba7d477d1e7021bbb43ea3379f69b107e1adcb92c7eb541). Preserve the designated paper and its credit. If Ben wishes to document an erratum, give the arithmetic in a dissertation-authored note outside the reproduced text. Verify whether the original 54 mm example was a separate measurement, since the main experimental range on p. 46 starts at 112 mm.

8. **A finite weighted objective prioritizes worst-case error; it does not mathematically forbid a trade.** Live pp. 94, 96, `30-torque.tex` lines 387, 414; Eqs. (4.33–4.36). “No shortfall at a single angle can be traded away for a better average” claims a lexicographic guarantee from a finite sum. An arbitrarily small worsening of the maximum can accompany a larger improvement across other poses. Say that the weight strongly prioritizes the worst angle, then cite the independently evaluated minimum-margin row of Table 4.10 for the returned design. Enforcing a genuine hard torque constraint would be a separate method change.

9. **Give a reference for the human band in the single-flexor result.** Live pp. 109–110, `30-torque.tex` lines 618, 624; Fig. 4.11(B). The paragraph calls a ±20% band “about the human mean”; the panel identifies only “Human model (max).” No population mean, sample, or band appears there. If it is a sensitivity band, identify the underlying single-muscle curve and the threshold used to obtain 78°–103°. If it represents human variation, supply the underlying evidence and show/identify the comparison. “Maximum human isometric flexion torque” should identify the selected single muscle. Direct readers to §4.1.6 for this constrained redesign.

10. **Repair the unit label on the pulse-readout calibration.** Live p. 229, `93-AppendixB.tex` lines 119, 123, 131. Table B3 labels 5.00 mV as a “Residual,” which would imply an error as large as the full-scale output. The primary file `Code/MuJoCo_SNS/spinal/spiking_calibration.json` records `g_inc_readout_check_mV=4.9975183311936044`; B.4.2 gives the full-scale target. The signed readout error is −0.00248 mV relative to 5 mV. Label the last column “Residual or check” and write “4.998 mV output (5 mV target),” or report the signed residual there.

11. **Direct the reader to the actual contact-release angles.** Live pp. 240–241, 251–252, `94-AppendixC.tex` lines 39, 42, 317. The phrase “each release angle in Table C2” sends the reader to a table of seed coordinates with no angles. The six first-inactive angles, from −89.8° for p5 to +7.4° for p8, appear in the prose on p. 240. Fig. C2 deliberately plots the *preceding* active poses. Point to that paragraph or put the six release angles in C2; keep the two pose definitions.

12. **Name the denominator for “force gain 1.”** Live pp. 119–120, 260–261, `30-torque.tex` lines 761–770; `94-AppendixC.tex` lines 403–408. The two-part tackle gives insertion force equal to half the two-BPA input, or equal to one BPA’s force. “Force gain 1, travel gain 2” beside Fig. C3 can be misread against Eq. (4.39). Say “output force equals one BPA’s force (half the pair’s total), with twice the travel of one BPA.”

### Figures and layout

- Figs. 4.7 and 4.8, pp. 101 and 103, `30-torque.tex` lines 492 and 524. The printed vector text for the legends is 5.2 pt, compared with 9.4 pt in Fig. 4.9 and 10.3 pt in Figs. 4.13–4.14. Legend entries occupy one corner of panel D, so a reader must identify the shared line styles through small gray and colored samples. A common legend across the width of each four-panel figure would permit larger lettering without covering measurements. On p. 102, line 518, the “loss of torque near full extension” chiefly describes the modest downturn in panel A. Panels B–D do not all include measurements near extension. Reference panel A if that is the observation intended.
- The live PDF does contain a nearly empty printed p. 122. The entire body text is “Chapter 6 draws on.” (`30-torque.tex` line 799); the line ends a paragraph on p. 121 just before Chapter 5. Shorten that bridge to one sentence, e.g. “Chapter 5 turns from actuator mechanics to the sensory evidence used in controller design.” This is a layout-only issue. Table 3.3 on p. 67 similarly floats eight pages beyond its first citation on p. 59; keep the protected table and paper intact.
- Preserve Figs. 4.13–4.14’s signed torque axes, measured versus pressure-specific markers, separate 620 kPa curve, and human reference label. Preserve Fig. 4.10’s apparatus/trace pairing. Fig. B1 labels its three conditions, the compared signals, both axes, and the 0.5–5 s window. The control with retained voltage-dependent synapses makes the comparison interpretable. B2/B3 distinguish connectivity from synaptic strength. Fig. C2’s last-active pose and Fig. C5’s distinct cam release sequence should retain their separate captions.
- C.1 introduces Table C1 twice (p. 239, `94-AppendixC.tex` lines 15, 21), interleaving pinned-flexor, pinned-extensor, and biomimetic route records. Combine the two flexor introductions and retain the specific note about the +5.3° encoder correction. In C.6 (p. 252, lines 317–321), file sizes and “plotting-title metadata” can be omitted from narrative without sacrificing geometry, event definitions, or reproducibility.

### What the evidence establishes

The published force work and the torque chapter address different physical questions. Normalizing force with independently measured free contraction provides a tractable actuator model across the tested pressures and lengths. The hybrid calculation then exposes the effect of errors in *geometry* at the joint. The coupled force balance, separate human-reference normalization in Eq. (4.38), and pressure-specific comparisons are worth preserving. The best demonstration of mechanical redesign is the sequence from the single-BPA flexor shortfall to the paired physical tests. The whole-leg results are explicitly predictions. Table C14’s 27 site entries give arithmetic means of approximately 37.35% and 65.81% torque coverage, subject to grid and configuration limits. They should remain separate from physical knee hardware validation. The calibration/ablation comparisons in B exemplify how to explain what stayed fixed and where a change in signaling failed to preserve motion.

The MAT files and workbook used here match their Git blobs at the finalized working checkpoint. No direct ProQuest artifact or heavy simulation was required for any listed check. Bench readings are interpreted as recorded; no sensor-error estimates or inferred corrected measurements were introduced. The compilers’ bookmark destinations can precede the headings after floats. Page ranges above follow the headings and figures on the printed pages.

## Neuromechanics and SADb evidence

### Findings in priority order

| Priority and location | Reader’s question and independent check | Smallest proposed repair |
|---|---|---|
| Accuracy. Ch. 6, p. 160, 60-sns.tex:265, “preserved that alternation”; Appendix D, pp. 309–311, 95-AppendixD.tex:381,387,393 | Which measure of alternation survived? Median left-strike phase within a right cycle was 0.498 intact, 0.494 with Ib deleted, 0.223 after the combined deletion. The intact recording had 34 left and 22 right detections. In a single sequence that count requires at least 11 same-side adjacencies, even when every right event splits two left runs. The historical gate requires at least six strikes per side and a median phase in [0.35,0.65]. It does not require one-for-one ordering. | Say “Deleting Ib input alone preserved the median heel-strike phase near half a cycle. The combined deletion shifted it to 0.223 cycles.” For Fig. 6.11 (p. 162, line 270), put the scoring-window medians in the caption, and identify the combined intervention over the right column. The pictured 12–15 s is a sample of the 4.5–30 s scoring interval. |
| Accuracy. Ch. 6, p. 142, 60-sns.tex:87, “within three percentage points of the human recording in the short trial” | Did this fit reproduce stance during sustained stepping? In Appendix D, p. 304, 95-AppendixD.tex:334, the comparison combines walking, reduction of drive, and return toward standing over 2–16 s. Right/left loads exceed 50 N for 62.3/61.4% of that window, but for 77.5/75.4% of the constant-drive window. Human reference is 61% per side. These are fractions of the analysis window. | Replace “in the short trial” with “over an interval that included the return to standing.” Add “During constant drive, both feet remained loaded for about three quarters of the interval.” Line 87 already mentions the longer-run failure. |
| Accuracy. Ch. 6, p. 141, 60-sns.tex:81, “mattered most for the rectus femoris and the short head of biceps femoris” | Most by which quantity? The paragraph reports RMSE (N). Table D2, p. 272, 95-AppendixD.tex:62–69 gives larger absolute error reductions for soleus (154.7 to 11.94 N) and medial gastrocnemius (125.4 to 15.40 N) than for the named pair (27.8 to 3.57 N and 39.4 to 1.95 N). The two named muscles have the largest gains in R² (0.2783, 0.1706, respectively). | “The two largest gains in R² occurred in rectus femoris and the short head of biceps femoris (Table D2).” Keep the absolute errors as supporting evidence. |
| Interpretation. Ch. 6, pp. 142–143, 60-sns.tex:91,93,98, Fig. 6.4 | Was the 0.25 N change in vastus medialis a pneumatic-pressure result? The prior paragraph identifies native biological actuators for an airborne model. The next experiment changes pressure at 54 sites. The standalone Fig. 6.4 caption identifies the model by its routes and bridge, but leaves the actuation in its pictured 25 ms run unspecified. Table 6.1 combines the two versions. | Identify the actuator law in this particular Simulink bridge trial in the caption or line 93, after checking the exact robotbody.xml used. Reserve the pressure-driven label for the verified pressure case in D.1.5. No numerical conclusion depends on an inference about the screenshot. |
| Interpretation. Ch. 6, p. 156, 60-sns.tex:228, “all 18 held-out tests”; Appendix D p. 306, 95-AppendixD.tex:364 | What was actually new after tuning? A retained ankle validation driver, Code/MuJoCo_SNS/spinal/postural_pilot_20261104/robust_adapter.py:9–16, uses two mirrored tilt/push directions × three added delays (0, 30, 70 ms) × three new noise seeds, with 0.05 nA noise. Training already used mirrored 16 N·m pushes and some of the delays at lower noise. These are 18 combined conditions, not 18 new mechanics or push magnitudes. | In Ch. 6 say “after tuning, the fixed controller recovered in all 18 reserved disturbance trials.” In D.3.6 add one clause with the factorial design after confirming the archived result uses this driver. |
| Figure interpretation. Ch. 6, pp. 150–151, Figs. 6.6A/6.7A, 60-sns.tex:161,170; 95-AppendixD.tex:101,103, pp. 274–275 | The inhibitory connection runs directly from the modeled afferent to the antagonist motor neuron. The diagrams contain four sensory cells and two MNs. Anatomical Ia reciprocal inhibition uses a spinal inhibitory interneuron as the dissertation itself explains on p. 134 (50-database.tex:114). In 6.6A “solid: excitation / dashed: inhibition” sits above a second legend where solid/dashed mean extensor/flexor for the traces. | State once that the signed links represent net effects through spinal circuitry. Move the extensor/flexor trace legend below the schematic or title it “Traces B–E.” Keep the plotted signals and neural topology. |
| Figure interpretation. Ch. 6, p. 161, Fig. 6.10, 60-sns.tex:259 | The reference bars are keyed by source, but the gray and dark segments lack stance/swing labels. Either bar-state assignment is consistent with the caption in isolation. The thin per-cycle model traces have no legend entry. Phase zero and unloaded model condition are already explicit. | Add “Dark bars show reference stance, gray bars swing; thin blue curves show single model cycles” after confirming the plotting mask. |
| Figure reference. App. D pp. 295–296, Fig. D3, 95-AppendixD.tex:265,270–273 | A reader citing “Figure D3(A)” could mean either applied current on p. 295 or hip angle on p. 296. The continued float restarts A–C, and the caption heading on the second page still reads “Figure D3”. The plotted windows agree. | Give the second page panels E–G and identify its caption as a continuation. |
| Reading flow. Ch. 4 to 5 pp. 121–123, 30-torque.tex:799 | The sole body text on printed p. 122 is “Chapter 6 draws on.” Ch. 5 begins on p. 123. An entire nearly blank page interrupts the mechanical-to-sensory handoff. | Shorten the tail of the transition so it fits on p. 121, then compile and inspect. For example, “Chapter 5 connects the controller’s sensory pathways to measured afferent responses.” |
| Optional. Ch. 5 pp. 124–125, 50-database.tex:34,38; Ch. 5 p. 133, line 107 | “The following example shows …” says what the example is for, but a heading separates it from the experiment. The AI ratings “adequate” and “good” can be mistaken for score thresholds. The 2024 author poster labels its three levels a qualitative evaluation. | Preserve the sentence’s purpose, optionally move it just after the subsection heading. Prefix “rating” on p. 133 with “qualitative” and leave the decision rule with the curator. |

### Positive findings and checks

* On pp. 125–129 (50-database.tex:38–67), the recorded experiment and proposed switch-deletion test are distinct. The local SADb Explorer paper record recBMsJS4yKuHNuo5, DOI 10.1113/jphysiol.1995.sp020897, contains both signs, the plantaris stimulus, medial gastrocnemius recording, MLR versus spinal conditions, the Ia tendon-stretch result, and the limits on inferring Ib alone. The original McCrea et al. (1995), *J Physiol* 487, 527–539, in Ben’s Zotero attachment 7TBPMM9A, confirms the at-rest/extension comparison (Fig. 3), selective Ia evidence (Fig. 4), and proposed circuit (Fig. 8, discussion p. 536). Fig. 5.3 represents a net response. Fig. 5.2D explicitly labels its text as summaries.
* The actual HTML data match the reported counts on pp. 130–133 (50-database.tex:73,95,101–103). Counting the app JSON yields 943 papers, 247 models, 157 reviews, 42 pathways, and four recorders. Of the papers, 833 have a curation note and 119 a named pathway. The citation payload yields 11,130 directed paper pairs and 816 endpoints. The Cat × Ib intersection is 21. Grillner record recFA1PeRSPTv9POu lists Forssberg reclhHfzxQojBL8Jh among its cited papers. The direct record and original-paper cross-check justify preserving the example and its body-space allocation.
* On pp. 134–135 (50-database.tex:116–122), Ting and Macpherson’s 2005 J Neurophysiol paper (local CJ3PCZPI) relates postural muscle synergies to endpoint force during support-surface translation. Deng et al. 2022 (N3MKG97M) study feedback and biomechanical entrainment in their rat model. Both illustrate the measured variable researchers would use when selecting a paper.
* Chapter 6 pp. 136–139, 60-sns.tex:4–61, names the original Li–Deng course experiment and Morrow’s collaboration, the tool limitations, and what the author contributed. Table 6.1 separates fixed-support, pressure-driven, native-muscle, and ideal tensile cases. On pp. 147–151, the biological knee pulse and the pneumatic arm load ramp are distinct experiments. Fig. 6.6B–F uses zero at pulse onset in both millisecond and second views. The pulse lasts 30 ms, and the knee reaches roughly 0.94° relative to its control. Fig. 6.7 shows 0.625 kg at 10 s and elevations consistent with 2.00° and 7.85° sag at the same time. The D.2.2 discussion guards against inferring measured human reflex latency, especially from a current injected downstream of the spindle.
* On pp. 153–158, the matched NumPy/Simulink input comparison and its two-millisecond relative timing distinguish a driven rhythm from spontaneous oscillation. Chapter 6 credits Hunt et al.’s staged method and Szczecinski et al.’s functional-subnetwork method. The local original Hunt et al. (2AH8F3JU, §2.4) gives the intrinsic-period and sensory-entrainment steps. Szczecinski et al. (2UF9XUFB) design circuit subnetworks for specified functions. D.1.4’s numerical example is checkable: 0.5 μS × 0.2 × 5 mV = 0.5 nA. Rybak’s deletions paper (attachment ZQCQ3K39, DOI 10.1113/jphysiol.2006.118703, esp. Fig. 8) supports the phase-preserving PF intervention behind D.2.3.
* On pp. 157–158, 291–293, the text discloses that both sensory return and sustained drive changed. The shorter periods cannot isolate a feedback effect. Figure 6.9 captions the moving, unloaded legs, their pelvic spring, and disabled ground contact. The ground examples report 50.5% body-weight harness lift for the crouched pressure-driven Gait2392 counterpart and 77.5% for the CAD model with ideal tension actuators. The latter achieved foot exchange with a small net backward displacement. At a specified CAD posture, requested force, filtered ideal force, and candidate BPA force are distinguished (p. 314). These limits are scientifically important.
* 60-sns.tex:116’s 7,868 connected neuron pairs is supported by the **full recipe**, Code/MuJoCo_SNS/spinal/fsa_design_20261101/DISSERTATION_WIRING_DECISIONS_20261102.md:13. The 7,866 default in an older regression check is a different configuration. Preserve the full-recipe value.
* Table 6.4, pp. 163–164, correctly scopes 18.31× to isolated neural updates and 6.82× to one entire matched process. (199.624 − 29.284) × 100 / 3600 = 4.73 h. D.6 qualifies the separate 1.016× coupled benchmark with its different timer and controller.
* Keep D.1.5’s distinction between 54 pressure commands and 44 routes that developed force in the loaded nominal trial (pp. 279–281). D.3.5 explicitly discloses the double counting of contact in an earlier recorder, the corrected threshold, and the distinction between the new scoring and signals that still drive the controller (pp. 302–304). The D.5 spiking knee has pneumatic actuation and pre-stimulus oscillation, unlike the nonspiking biological knee. D.6.1 limits the hardware claim to interfaces awaiting sensor-to-valve timing measurements.

### Limits and dependencies

The original acute-spinal McCrea result, label membership and app counts were independently checked. The Grillner/Forssberg link was tested as a directed graph edge; transfer of the detailed Forssberg stance result was not independently rederived from its recording. The subset of primary papers checked is described above. The screenshots and all figure callouts in these sections were inspected in the rendered live PDF, including zoomed views of the shared time axes and labels. Identical raw arrays, fit reproducibility, hardware response and new heavy simulations are separate verification tasks. In particular the difference between the p. 309 intact baseline (28 L / 21 R) and p. 311 ablation-series control (34 L / 22 R) is not an established numerical error. Give each its configuration or recording identifier if they are intentionally distinct; the remedy does not require harmonizing the counts. These priority judgments are the reader’s own, not observed advisor grades.

### Historical calibration

Hunt’s captured 14 October comment C15 (ajh26, 4:09 p.m., HISTORICAL_SIGNALS.json record OCT-C-39bd4fbbd654af42, source 15_report line 296) asked for “a clear storyline of what is being built and why up front” on the former standalone prior-walker heading. In the live pp. 136–139, lines 4–61, a reader can now name the scientific need, the failed ground release motivating it, Bolen/Morrow’s original experiment, and the new graphical coupling and tuning tools. This is an independent assessment of the current answer to the historical question. The historical event’s −2 belongs only to that recorded event. Ben’s November 20 note “Harness-assisted seems to be more specific” (EDITORIAL_EXAMPLES.json, BOLEN-NOV20-APPARATUS) is author feedback. The distinct spring, rail, contact, and drive descriptions satisfy that request in the cited cases. Neither present judgment implies subsequent advisor praise.

The local Zotero database lists Forssberg DOI 10.1152/jn.1979.42.4.936 under MUTCUC9Z and HQHQALMG without an attached primary PDF. The verified claim about Grillner citing Forssberg refers to the local app’s directed link. Rybak et al.’s cited deletions paper (not their same-year companion on afferent stimulation) describes RG-E resetting in Fig. 6 and PF-E non-resetting in Fig. 8. The pulse in the local knee demonstration tests a circuit implementation, so physiology and measured-latency claims remain bounded as in D.2.2.

## Advisor-calibrated reading and observed feedback

### Findings for Ben

1. **Accuracy, pp.20–21 and 214. OpenSim moment arms.** In 15-background.tex:96, “OpenSim uses length sensitivities” follows the tendon-excursion explanation. In 92-AppendixA.tex:158–167 the signed derivative is then introduced for each human path, including paths with moving points. A researcher attempting to reproduce the reference curves could finite-difference musculotendon lengths and silently obtain a different answer. Sherman, Seth, and Delp distinguish the effective joint torque divided by tendon tension from the tendon-excursion expression. The expressions agree under the workless-constraint assumptions; prescribed moving points can violate them. In the documented OpenSim geometry-path solver, unit tendon tension is applied, the forces are mapped to generalized forces, and the coordinate coupling is included. **Small repair:** Specify the pathway and coordinate used for the exported reference, state the workless condition for the length derivative, and identify whether the archived curves came from the OpenSim moment-arm output or an explicit length-difference postprocess. This review establishes a wording defect, not a changed torque target. Primary evidence: Sherman et al., *What Is a Moment Arm?*, local Zotero attachment 3ZECQB39, PDF pp.3–4 and 7, and [OpenSim MomentArmSolver](https://github.com/opensim-org/opensim-core/blob/main/OpenSim/Simulation/MomentArmSolver.cpp).
2. **Accuracy, pp.18 and 67. Fit quality versus transfer to another specimen.** 15-background.tex:60 calls Martens and Boblan “even more accurate” immediately after quoting two coefficients of determination from Sárosi. The examples refer to different actuators and error measures. On the dissertation’s common 20 mm, 300 mm specimen, maximum absolute errors are 137.8 N (Martens and Boblan) and 313.4 N (Sárosi). At 450 mm the ordering reverses, 1321.1 N and 294.8 N, respectively (20-muscles.tex:378–382). **Small repair:** Delete the cross-study ranking in Chapter 2 or qualify it to the original actuators and name the metric and comparison that supports it. Point the reader to Table 3.3 for the matched comparison. The published chapter requires no edit.
3. **Navigation, contents iv–x.** On the live PDF, 25 of 125 chapter, section, and subsection links matched to their printed contents entries point to a different physical page from the printed page number. Examples independently verified against actual headings: §1.2 links to p.4, heading on p.5 (10-introduction.tex:18); §6.2.1 links to p.138, heading on p.140 (60-sns.tex:66); §6.4.3 links to p.149, heading on p.152 (60-sns.tex:174); §7.2.2 links to p.173, heading on p.174. The clickable TOC entries and outline bookmarks share the wrong named destinations. PyMuPDF and pypdf agreed on these samples. Same scale of defect in the working and rebuilt ETD PDFs, 25/125 each; individual affected entries vary. **Small repair:** Correct heading-anchor creation in the shared preamble/class, then verify destinations against printed headings in a disposable build. Inspect hyperref/titlesec interaction (main.tex:6,25; psutex.cls:440–463) before touching individual chapter headings. Cause remains untested. The complete cross-version navigation count is in the accompanying version comparison.
4. **Reference accuracy, p.41.** 15-background.tex:290 directs the reader to §1.2 for the “loop described.” Section 1.2 introduces the laboratory platform. Section 1.1, pp.1–2, explicitly defines the two halves. **Small repair for Ben's approval:** change the cross-reference to §1.1 (add a label at 10-introduction.tex:6). The technical check calls for no change to Ben’s ethical position or wording.
5. **Graphical interpretation, pp.23 and 25.** The text says the torque surfaces “differ substantially” (15-background.tex:134). Figure 2.7 E and F (caption line 149) reverse the apparent front-axis order for knee and hip angle, and each uses its own color scale. The human hip panel includes positive values; the robot hip panel’s color scale tops out at zero. Matching yellow patches therefore conceal a mechanical difference. **Small repair:** tell readers in the caption that axes/view and color scales differ and specify the difference being illustrated in the preceding sentence. Retain credit to Bolen (2019), Fig.8.
6. **Source credit, dedication ii.** 04-dedication.tex:20–23 marks the four Spanish lines as a quotation but displays no source. The opening in the dedication, “Los corazones obreros que laten”, matches the line under “HIJOS DEL PUEBLO” in CNT’s centenary dossier, followed by the same stanza. **Small repair:** if Ben wants to identify the quotation, add the song title and the source or version he used beneath it. The dossier and other editions vary on historical details, so the title is the reliable minimum. [CNT centenary dossier, section on the historic songs](https://www.cnt.es/wp-content/uploads/anterior/DossierPremsaCentenari.pdf).
7. **Reader continuity, pp.31–34.** The last sentence of p.31 stops after “recordings or with”, and two full pages of preceding-section figures come before “responses to stimulation and pathway removal” on p.34 (15-background.tex, §2.7 opening). **Small repair:** retain the §2.6 figures but place them before the §2.7 heading, or use a float barrier before §2.7 in a test build. Check pagination. Figures 2.11–2.12 themselves have substantial explanatory content.
8. **Reader scope, abstract i and p.185.** “Virtual-harness-assisted balance and walking behaviors” (02-abstract.tex:7) could invite a reader to identify the balance and walking tasks as the same experiment. The final Chapter 8 paragraph on controllers identifies the supported W2L body and retains independent balance as open; Appendix D.3.6 describes separate fixed-base pendulum balance. The abstract already qualifies its claim. Preserve the approved wording. If an examiner asks, distinguish the selected balance task, supported stepping, and independent balance.
9. **Layout, x and xvii.** Actual page images confirm “Harness317” and “Parameters319” in TOC rows D.5.1 and D.5.3 and almost no clearance between “ground” and “158” on Fig.6.9’s list entry. Minimal fix is shorter optional contents titles at 95-AppendixD.tex:438,473 and a shorter optional caption at 60-sns.tex:245. Full on-page headings and captions can remain. This is independent of the clickable-destination defect.
10. **Detail and purpose, p.173.** The proposed pressure-ramp tests end with a reference to biological muscle formulations, after the paragraph's conclusion on wall shape versus measured held-length force (70-overall.tex:39). Appendix A cannot answer the Festo deformation question being discussed. If shortening Chapter 7, cut that sentence here; the human-model formulations are already signposted in the earlier technical chapters.
11. **Evidence-sensitive wording, p.2.** The motivation states “subsequent robotic experiments generated predictions” (10-introduction.tex:8). The cited Ijspeert paper attributes the prediction of different oscillator frequencies to the spinal model; the robot tests the circuit’s locomotor behavior; biological data support the frequency prediction. A precision edit, if desired, is “the model predicted a frequency difference that the authors tested in the salamander spinal cord.” Preserve the explicitly praised two-way motivation and the role of the robot. [Ijspeert et al., *Science* (2007), pp.1416–1420](https://pubmed.ncbi.nlm.nih.gov/17347441/). Verified against local Zotero WQAH7HJL, PDF p.5.

### Advisor evidence and calibration

The scratch classifier output used the handoff’s own review.py and all 18 explicitly enumerated *.tex chapter, preliminary, bibliography, and appendix files from the live ZIP. With --check-results-openings it returned exit code 0, 1,220 prose/caption/table units, and 67 opening records across six files, including Chapter 6 and Appendix D. Zero Results-opening flags were raised. No prior passage-decision file was supplied; its zero stale or reopened counts cannot certify acceptance. These tools rank similarities to earlier critiques. Findings 1–7 required reading, figure inspection, source comparison, or inspecting PDF destinations.

The historical ledger includes 388 pre-November attributable events and 98 unique November events. July: 4, −7; September: 64, −98; October: 320, −404; November: 98, −126. Overall historical sum −635 for 486 events. The November 18 capture includes November 17 and contains 41 events, −42; November 19 adds 15, −20; November 20 adds 42, −64 including the positive event. Repeated captures are counted once. Distinct tracked changes cost one each on the earlier wording. Substantive and interpretive comments cost two or three as recorded. Explicitly identifiable praise follows the applicable ledger rule; silence is unlabeled. None of these totals grades Ben’s live dissertation.

- **Advisor, November 18, INTRO-T36.** Original reads “supports the selection of feedback pathways and the definition of biologically informed interactions in neuromechanical models.” The tracked replacement says the chapter gives an overview of development, organization, and features for sorting and understanding material. The live closing sentence at 10-introduction.tex:14 retains this distinction. Chapter 7, p.175 (70-overall.tex:57), supplies the concrete application: the same stimulus led to different signs under distinct conditions. The reader can distinguish searching evidence from assigning a controller weight.
- **Advisor, November 19 and 20, Introduction.** The November 19 mixed comment said the paragraph was significantly better but that the highlighted section jumped around. The −2 applies to that draft. On November 20 the identifiable revised motivation received “Much better.” (+1). Preserve the two-way loop. A source-specific correction to the salamander example would answer a narrower question.
- **Advisor, November 19–20, SNS.** “Who is ‘our’” sought ownership; the two filler comments sought less narration. On November 20 “Is this part of myoconverter or something you did extra?” and “what is wrong and how you fixed it” requested concrete ownership and consequences. Current 10-introduction.tex:16 and 80-conclusion.tex:12 explain what the tools let a researcher hold fixed and what the W2L deletion actually changed. First person at the point of contribution is appropriate.
- **Advisor, November 20, NOV20-05 and -07.** Asked why one paper appeared and why the narrative returned to it. Section 7.2.3, p.175, names the McCrea stimulus/sign contrast and explains the choice it informs. Preserve that example. **Ben, November 20:** “Harness-assisted seems to be more specific” and the animal-experiment examples belong to the author ledger; no advisor points attach. Section 7.2.4 explains ground contact with upward harness force and the retained shared drive when supplementary drive reached zero, directly answering that concern.
- **Earlier local matrix.** Archive/Notes/advisor_revision_20261023/Revision_E/review/self_scoring_matrix.json has 12 local draft scores of 2 and 12 acceptance scores of 0. Its own definition of 2 is local revision with evidence, and 0 acceptance means unaccepted. The separate October 0–4 readiness rubric reserves 3 for a context/source-verified revision and 4 for author/advisor confirmation. The September B− concerned the coherence of the combined papers. These scales have different units, dates, and owners. No numeric current approval is inferred.

### Evidence scope

Used C:/Users/Ben Bolen/Dissertation_Handoff/Advisor_Scoring/Notes, the local date-stamped November 18–20 classifier captures, 32 editorial examples, and the historical signals; compared their actual comment and replacement strings. Checked repo Notes/advisor_feedback_20261120, its live anchor files, historical archive matrix, the local Codex sessions index and available November rollouts, and local skill directories. Local dissertation skill and plugin diss-style-contract were available under Documents/GitHub/ZCode_Skills. Zotero originals for Sherman, Selionov, Brown, Asano, and Ijspeert were readable. The Selionov study used suspended human legs, confirming the scope in §2.10, and Brown describes rhythm in low spinal preparations with deafferented muscles. Russell and Burch’s original account supports the Three Rs classifications. The 2025 classroom recording’s uploader/title is independently identifiable, but the film-specific procedure sequence and the comparison with the horse require a time-coded transcript or viewing of that specific footage to verify. Treat those details as unverified rather than erroneous. Figure 2.13 credits its source, labels restoration and machine colorization, and confines its role to documentary context. [Russell and Burch, original text](https://caat.jhsph.edu/wp-content/uploads/2024/08/17.-The-Removal-of-Inhumanity.pdf); [Selionov et al. (2009)](https://pubmed.ncbi.nlm.nih.gov/19339461/); [UTHealth Houston School of Dentistry recording](https://www.youtube.com/watch?v=1nkczBM5YiM).

The scientific narrative at p.170 distinguishes tube end-effects from cap/attachment compliance and assigns each explanation a separate measurement. Figure 7.2 A/B and continued C support this distinction. The proposed 500 kPa protocol is prospective. Keep this paragraph and the Chapter 8 measured-angle qualifier. No simulator was run for this reading.

**Machine-local chat and skill search.** The Codex sessions index on EB475WS4 lists earlier dissertation review threads. Reading original user-role messages in the October 23/26 “Compile advisor comments” lineage recovered Ben’s instruction to include his reply to the advisor (October 24 01:01 UTC) and his direction to present tools and relevant verification (October 26 13:08 UTC). The October 29 “Polish dissertation Chapters 4 and 5” thread records Ben’s requested heading “Exceeding Human Torque Magnitude.” Those are author rulings. Four locally accessible November rollout files were searched for dissertation/advisor messages; the relevant user messages concerned defense slides. No newer attributable advisor mark was recovered from those rollouts. The latest dated observed manuscript marks in the supplied captures were November 20. Searched the local ZCode skills checkout and plugin marketplace/cache, Codex plugin cache, user skill locations, the Academic_Writing checkout, and Bipedal_Robot for a separately named dissertation-advisor-review or technical-writing SKILL.md. Neither was found in those locations. The dissertation skill, diss-style-contract, actual attributed ledger, and local October matrix were available.

## Attribution and bibliography evidence

### Figure 2.4 credits the wrong publication

**Accuracy and reuse credit.** On p. 22, Fig. 2.4, the caption credits both images to Steele et al. (2018), Fig. 6, under CC BY 4.0. See `chapters/15-background.tex:106–110`, the embedded asset `figs/15-background/steeleknee.pdf`, and `thesis.bib:993`. The cited 2018 Fig. 6 contains the single mechanical-stop view. The two images in the dissertation match the two panels of Fig. 6 on p. 014503-4 of Steele, Etoundi, and Hunt, *Experimental Verification of Kinematics and Kinetics in a Biomimetic Bipedal Robot*. Panel A shows the sagittal mechanism at full extension; panel B shows the stop and its blue link. I compared the three graphics locally. The ASME article appeared online in 2022 and in the February 2023 issue. Its first page identifies the CC BY 4.0 license.

**Smallest repair.** Cite the ASME article and its Fig. 6 in this caption, keep the license credit, and update the source comment. The separate 2018 conference reference supports the surrounding design results, including the 0.079 sliding-to-rolling ratio and 151° range. [ASME article](https://doi.org/10.1115/1.4054441) · [Published PDF with figure and license](https://uwe-repository.worktribe.com/OutputFile/9689241).

### Two titles change the scientific meaning

| Live location | Reader consequence and evidence | Smallest repair |
|---|---|---|
| p. 198; `thesis.bib:1375–1382`; cited p. 34, `chapters/15-background.tex:238` | Nourse et al. (2019) is described as a study “in Stick Insect.” The original studies *Drosophila*. The reference credits Hunt and omits three contributors. The DOI leads to a chapter authored by William Nourse, Nicholas Szczecinski, Moritz Haustein, Till Bockemühl, Ansgar Büschges, and Roger Quinn. The body’s “insect-inspired inter-leg coordination” remains valid. | Import the original title and six authors into the existing entry. [Chapter DOI](https://doi.org/10.1007/978-3-030-24741-6_34) · [Coauthor’s publication list](https://zoologie.uni-koeln.de/arbeitsgruppen/ag-bueschges/inhalt/mitarbeiterinnen/dr-moritz-haustein). |
| p. 194; `thesis.bib:1396`; p. 34, `chapters/15-background.tex:238`; p. 179, `chapters/70-overall.tex:89` | The reference title ends “Similar to Human Test Subjects.” The published title ends “Similar to an Analogous Classical Controller.” The original compares a neural balance controller with its classical counterpart. A reader following the listing could infer a human-subject validation. | Restore the published title. On p. 34, “models of human balance control” would clarify what the paper contributes. At the hardware claim on p. 179, attribute each demonstration to the paper or poster that performed it. [Author/institutional article record](https://pdxscholar.library.pdx.edu/mengin_fac/334/). |

The same survey sentence on p. 34 calls the work of Schilling et al. “adaptive hindlimb walking controllers.” The cited study develops a model for adaptive hexapod walking. The current wording can suggest a mammalian hindlimb study. Use “adaptive control of hexapod walking.” [Author/institutional publication record](https://www.hsbi.de/publikationsserver/record/3049).

### Additional confirmed reference corrections

Each row identifies the field to change in the existing entry. These corrections affect attribution, searchability, and retrieval. They do not require rewriting the protected Chapter 3.

| Live printed page | Source line | Current entry and confirmed repair |
|---|---|---|
| 186 | `thesis.bib:1274` | Akay et al. (2014) credits “William G. Tourtellotte.” The published byline gives **Warren G. Tourtellotte**. [Original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC4250167/). |
| 186 | `thesis.bib:1342` | Akay and Murray (2021) ends “An Overview.” The subtitle is **Opportunities for Discovery in the Age of Molecular Science**. [Publisher](https://www.mdpi.com/1422-0067/22/3/1467). |
| 196 | `thesis.bib:1646` | Li et al. (2017) credits “Wansheng Li.” The journal article credits **Wei Li**. A similarly titled conference paper has a different byline. Use the journal DOI in this entry. [Journal record](https://pubmed.ncbi.nlm.nih.gov/28748830/). |
| 196 | `thesis.bib:1816` | Li et al. (2024) credits “Yinlong Li, Yinlai Jiang, Kou Hosoda.” The published names are **Yiqi Li, Yelin Jiang, and Koh Hosoda**. [Original article](https://doi.org/10.3389/frobt.2024.1449721). |
| 198 | `thesis.bib:1370` | Nourse et al. (2023), the SNS-Toolbox paper, credits “Axel Jackson.” The publisher and the local article’s title page give **Clayton Jackson**. [Publisher](https://www.mdpi.com/2313-7673/8/2/247). |
| 199 | `thesis.bib:1201` | Rybak et al. (2015) credits “Kevin J. Dougherty.” The article credits **Kimberly J. Dougherty**. [Original article](https://www.eneuro.org/content/2/5/ENEURO.0069-15.2015). |
| 199 | `thesis.bib:1181` | Rybak et al. (2006), the fictive-locomotion deletions paper, starts at 627 in the bibliography. Correct the range to **617–639**. [Original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC1890439/). |
| 204 | `thesis.bib:1435` | Wang et al. (2022), MyoSim, lists pp. 10410–10417. The ICRA paper runs **8104–8111**. [Author/institutional record](https://research.utwente.nl/en/publications/myosim-fast-and-physiologically-realistic-mujoco-models-for-muscu/). |
| 206 | `thesis.bib:1732` | Zhagiparova et al. (2025) prints “Altynai,” “Viktor,” and “Daekyu.” The publisher credits **Aliya Zhagiparova, Vladimir Golubev, and Daewon Kim**. [Original article](https://www.mdpi.com/2076-0825/14/12/582). |
| 186 | `thesis.bib:1796` | Armstead et al. (2025) credits “Bryson Armstead” and “Christian Jackson.” The conference’s official award for this paper names **Benjamin Armstead and Clayton Jackson**. [Official conference record](https://livingmachinesconference.eu/2025/awards/). |
| 190 | `thesis.bib:1675` | The AnimatLab DOI ending **2010.02.005** links to an unrelated study of mirror-chamber and elevated-plus-maze tests in mice. Change that month component to **2010.01.005**. The other journal coordinates (187, 280–288) already match AnimatLab. [Correct paper](https://pubmed.ncbi.nlm.nih.gov/20074588/) · [Current DOI’s destination](https://pubmed.ncbi.nlm.nih.gov/20149823/). |
| 192 | `thesis.bib:1885` | Fukuoka et al. (2021) lists “Yasuyuki Habu” and “Ning Inoue,” and omits two authors. Restore **Yasuhiro Fukuoka, Yasushi Habu, Kouta Inoue, Satoshi Ogura, and Yoshikazu Mori**. [Published PDF, p. 1](https://journals.sagepub.com/doi/pdf/10.1177/17298814211044936). |
| 193 | `thesis.bib:1784` | Grummer et al. (2024) credits “Christian B. Jackson.” The conference program names **Clayton Jackson**, published as **Clayton B. Jackson**. [Official conference program](https://livingmachinesconference.eu/2024/posters/) · [Chapter](https://doi.org/10.1007/978-3-031-72597-5_26). |
| 193 | `thesis.bib:1865` | Haeufle et al. (2018) credits “Benjamin Schmortte” and “Andre Seyfarth.” The original authors are **Daniel F. B. Haeufle, Birgit Schmortte, Hartmut Geyer, Roy Müller, and Syn Schmitt**. Two actual contributors have been omitted. [Original article](https://doi.org/10.3389/fncom.2018.00080). |
| 194 | `thesis.bib:1875` | Ichimura and Yamazaki (2019) prints “Tetsuro Yamazaki.” Correct to **Tadashi Yamazaki**. [Original byline](https://pubmed.ncbi.nlm.nih.gov/31616276/). |

### Checks that support the existing manuscript

- Chapter 5’s McCrea example (pp. 125–129) distinguishes the group I stimulation actually delivered from the evidence about individual afferent classes. The primary paper compares resting and locomotor conditions; preserve this context when tightening the passage. [Original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC1156590/).
- On p. 32, Fig. 2.11, `chapters/15-background.tex:216` credits the source diagram to Di Russo et al., Fig. 5. The article and source figure agree. On p. 33 the Shinohara Fig. 2.12 credit identifies a preprint; the local author PDF lists the credited contributors and the CC BY 4.0 license. [Di Russo et al., published PDF](https://mitkof6.gitlab.io/personal-site/publications/jne2023/dirusso-jne-2023.pdf).
- Domínguez-Rodríguez et al. (2020), p. 191, `thesis.bib:2097–2104`, agrees with the accepted manuscript’s title, byline, and pp. 96–112. The reference describes **candidate** interneurons and identifies cat as the species. This is a useful limit on how strongly a model connection can be called experimentally established.
- Nourse et al. (2023) already demonstrated SNS-Toolbox coupled to MuJoCo with sensory feedback. The coupling belongs to the existing tools. Credit Ben for his specific models, circuits, and integration choices. [Original article, Fig. 6](https://www.mdpi.com/2313-7673/8/2/247).

### Evidence scope

I screened the 129 distinct explicit DOI entries cited by the live chapter and appendix sources against the local Zotero database. The local library had 102 matching DOI strings. A match enabled comparison of fields; it did not itself authenticate the record. I confirmed the corrections above against publisher, author, and bibliographic records and the available original PDFs. Local Zotero attachment keys for reinspection include SNXQVWXW and 4MQUGRSR (Steele), 8RYYKWP4 (Hilts), 5IJZ32B2 (SNS-Toolbox), HXP3SJV9 (Domínguez-Rodríguez), D6W75UDD (MyoSim), and 4PXKC5UR (Shinohara). This was a targeted primary-source audit. It did not authenticate all bibliography fields or compare the entire dissertation against all cited texts for distinctive wording. I found no verified pair warranting an unmarked-quotation finding in the inspected passages.

### Findings requiring correction

The highlighted bibliography errors change authorship, the identity of an experiment, or the information a reader uses to locate it. They are independent of the scientific merits of the underlying papers. The title fragments below identify the printed entries. Repair the shared bibliography from the intended paper’s publication record. These corrections can leave the protected Chapter 3 prose intact.

| Live page | Source in live thesis.bib | Printed evidence and reader impact | Verified small repair and source |
|---|---|---|---|
| 196 | 1507–1515, martens_framework_2018 | “A Sensor-Less Torque Control Interface” credits Martin Martens, Alexander Passon, and Ingo Boblan. Its DOI leads to a different title and authors. Chapter 3 cites the key at 20-muscles.tex:407. | If the 2018 journal work is intended, use its framework title and byline Mirco Martens, Thomas Seel, Johannes Zawatzki, and Ivo Boblan. If specifically citing the 2017 conference sensor-less paper, use the distinct record in reference 26 of that journal article. [Publisher, title and references](https://www.mdpi.com/2076-0825/7/4/82) |
| 204 | 1451–1459, ugurlu2019stable | “Pneumatic Muscle Actuator Joint” credits Dallali, Caldwell, and Tsagarakis and gives DOI 10.1109/TII.2018.2849330, 15(2), 1240–1251. This prevents retrieval of the work invoked at 20-muscles.tex:407. | The torque/stiffness and position/stiffness study is by Barkan Ugurlu, Paolo Forni, Corinne Doppmann, Emre Sariyildiz, and Jun Morimoto. Use the paper’s title, its DOI 10.1109/TII.2019.2916228, and 15(12), 6270–6279. [Original paper at authors’ laboratory](https://bicr.atr.jp/bri/wp-content/uploads/2020/01/Stable-control-of-force_Barkan.pdf) |
| 198 | 1375–1382, nourse_analyzing_2019 | “in Stick Insect” changes the experimental species and the printed byline adds Alexander Hunt while omitting three contributors. Ch. 2 uses the general description “insect-inspired inter-leg coordination” at 15-background.tex:238. | The experiment is in *Drosophila*. Restore the byline Nourse, Szczecinski, Haustein, Bockemühl, Büschges, and Quinn. The chapter’s general description can stand. [Published volume contents](https://opac.admin.ch/toc/toc1902092622.pdf) · [Full bibliographic record](https://dblp.org/rec/conf/lm/NourseSHBBQ19) |
| 193 | 1864–1871, haeufle_benefit_2018 | The step-down paper credits “Benjamin Schmortte” and “Andre Seyfarth” while omitting two authors. | Use Daniel F. B. Haeufle, Birgit Schmortte, Hartmut Geyer, Roy Müller, and Syn Schmitt. [Publisher](https://doi.org/10.3389/fncom.2018.00080) |
| 186 and 198 | 1795–1796, armstead_design_2025, and 1364–1370, nourse_sns_2023 | The robot-cat paper names “Bryson Armstead, Christian Jackson.” The SNS-Toolbox paper names “Axel Jackson.” Readers cannot reliably identify the contributors to the tools described in Chapter 2. | Conference lists Benjamin Armstead and Clayton Jackson for the cat. Clayton Jackson also coauthored SNS-Toolbox. Import each actual byline. [Conference](https://livingmachinesconference.eu/2025/awards/) · [SNS-Toolbox publisher](https://www.mdpi.com/2313-7673/8/2/247) |
| 206 | 1731–1732, zhagiparova_recent_2025 | “Altynai Zhagiparova, Viktor Golubev, and Daekyu Kim.” The first names of all three authors differ from the paper. | Use Aliya Zhagiparova, Vladimir Golubev, and Daewon Kim. [Publisher](https://www.mdpi.com/2076-0825/14/12/582) |
| 200 and 204 | 532 onward, seth_opensim_2011, and 1429–1435, wang_myosim_2022 | The 2011 OpenSim reference ends at the title and year; the MyoSim conference pages read “10410–10417.” | Include *Procedia IUTAM* 2, 212–232, DOI 10.1016/j.piutam.2011.04.021 for Seth et al. Correct Wang et al. to 8104–8111. [OpenSim article](https://pmc.ncbi.nlm.nih.gov/articles/PMC4397580/) · [MyoSim authors’ university record](https://research.utwente.nl/en/publications/myosim-fast-and-physiologically-realistic-mujoco-models-for-muscu/) |

Inspect the other active entries against the original bylines before treating the bibliography as verified. Citation resolution and a successful LaTeX build cannot detect substitutions within an existing reference.

**Frey et al., p. 192, thesis.bib:1713–1720.** The live entry lists Roi Mattamala and Michael A. Osborne. The original byline names Matías Mattamala and Maike Osborne. The article identifier is eaee0787. Keep the DOI suffix aee0787. [PubMed record](https://pubmed.ncbi.nlm.nih.gov/42525724/)

## Figure readability, navigation, and optional layout repairs

### Figure readability and navigation

| Live page and location | Evidence in the printed PDF | Minimal repair |
|---|---|---|
| 295–296, Figure D3, 95-AppendixD.tex:265 and :273 | Both pages say “Figure D3.” A–D on p. 295 are added current and neural voltages. A–C on p. 296 are joint angles. D3B thus names either RG voltage or knee angle. Time windows and 100-ms shaded pulse are aligned on both pages. | Begin the second caption “Continued,” and relabel second-page panels E–G. Figure 7.2, pp. 172–173, already demonstrates how to continue to C across a page break. Preserve actual recorded-time axes. |
| 101 and 103, Figures 4.7 and 4.8, 30-torque.tex:492 and :524 | In the embedded plots “Measured,” “Unoptimized,” and the other legend entries are **5.18 pt**; major tick labels are approximately 7.8 pt. Readers must distinguish the observed, hybrid, baseline, and fitted responses from these entries. | Keep panel arrangement. Move one shared legend below the plot grid and enlarge the legend and tick labels at final placement. |
| 117–118, Figures 4.13 and 4.14, 30-torque.tex:749 and :756 | Measured and predicted torques both use filled circles. At some angles one partly covers the other, including near −14° in Figure 4.14. Color distinguishes these two series. | Use open contrasting diamonds for predictions and filled circles for observations. Keep pressure-scaled marker areas where applicable and the flexor’s negative torque convention. |
| 140–141, Figure 6.3, 60-sns.tex:70–75 | The reader has just learned that the rectus femoris path was corrected at the knee. The full-body models are about 106 PDF points (37 mm) tall; image widths are 21 and 36 points. They identify platforms, but an eye cannot follow that path around the patella at this size. | Keep the orientation views and add a labelled knee crop if the figure is meant to show the named repair. Otherwise refer to this figure explicitly as the initial pose only. |
| 121–122, Chapter 4 close, 30-torque.tex:799 | “Chapter 6 draws on.” is the entire body of printed p. 122. This leaves the reader turning a page to finish a four-word dependent clause. | Shorten the transition by one composed line, e.g., “Chapter 5 identifies sensory responses that can guide the controller tests in Chapter 6.” Recheck the resulting page break. |
| xvii, List of Figures, 60-sns.tex:245 | Entry 6.9 actually prints “ground158” with no visible break between the last word and linked page number. | Shorten its optional list caption, e.g., “BiW2Lvar+ coordination with both feet off the ground.” |
| 218, Appendix A, 92-AppendixA.tex:232 | Three coefficients for an unused 40 mm actuator interrupt the replication instructions for the studied force fits. Readers could mistake the adjacent coefficients for the operative fit. | Retain the Table 3.2 reference and code locator; omit the unused coefficients in the review draft if Ben agrees. |

### What survives inspection

- Figure 1.1 (10) has labelled Background-to-actuation and Background-to-SADb paths and a dashed future integration route. The arrowheads are present in the live export.
- Figures 5.2–5.3 (128–129) lead from the actual McCrea article and pathway fields to opposite net actions at rest and during extensor activity. The artwork identifies the observed input as group I, and the text describes how a curator’s Ib tag retrieves the broader group-I study. The source credit and proposed status are explicit. Figure 5.5 retains its labelled map callouts.
- Figures 6.6–6.7 (150–151) use a 30-ms Ia step at zero for the knee and a loading ramp in seconds for the arm. The knee’s two axes explicitly switch from ms to s and define zero in the caption. There is no evidence of an axis-origin error in these figures. Figure 6.9 (158) places four signals on 7–11 s. Figure 6.11 (162) distinguishes its 12–15 s illustration from the longer metric window. Figure D6 (313) shows native flexion-negative joint angles over 16–28 s, includes reference lines and soft-contact explanation, and identifies ideal actuators.
- The human OpenSim/MuJoCo pair (6.3), the rerouted robot model and actuator interface (6.4), and the torque-recruited ideal-actuator CAD case (6.12 and D6) have different caption descriptions. The CAD result explicitly states the harness support. Keep these descriptions with the artwork.
- Appendix A (211–216) gives a concrete consequence of a shape factor, states the sign of shortening, and evaluates a half-speed example. Its active MuJoCo shortening expression gives (1−0.5)² = 0.25, as reported. Appendix A’s separate description of minimizeFlxPin.m at line 260 matches the current executable transforms in that MATLAB file at lines 289 and 306 despite conflicting historical comments inside the MATLAB file. No change is supported on those line anchors.
- Appendix E (329–333) distinguishes records from unique studies. 943+247+157+42+4 = 1,393, 561+244+109+39+4 = 957, and the difference is 436. The 382 archived Papers plus 54 other records reconcile. The overlapping Ia/Ib example explains why adding Pivot rows double-counts a paper. Table E3 preserves the distinction between a recorded pathway name and the group-I fibers actually stimulated.
- Appendix D’s split Figure D3 and the two-page Figure 7.2 are deliberate in the live source. The issue specific to D3 is its repeated panel letters and absent continuation label. The standalone appendix divider (207), sideways tables, and appendix-script longtable with repeated headers are intentional. Caption text, roman pagination and the 1–333 sequence all render. Chapter float pages with generous margins do not, by themselves, establish a fault.

**Verification boundary.** The bibliography was read end to end and selected high-consequence or suspicious records were checked individually against the sources listed. Those checks are not an exhaustive bibliographic or close-paraphrase audit. Diagram timing and labels were checked against the printed work and matching sources; the other readers checked selected underlying test records and image credits.

**Electronic access.** The PDF catalogs have no structure tree. None of the three editions contains an indirect object with an /Alt field. The live ZIP includes descriptive sidecars, for example figs/60-sns/gait2392_human_models.alt.txt, but 60-sns.tex:74–75 does not link that description to a tagged figure in the distributed PDF. If the electronic edition is required to expose those descriptions to screen-reader users, test one tagged figure and its reading order before changing the full build. The supplied artifact cannot establish the institution's submission requirements.

## Reading coverage

I read the live PDF in printed order, including the title, copyright, abstract (i), dedication (ii), acknowledgments (iii), contents (iv–x), list of tables (xi–xiii), list of figures (xiv–xviii), Chapters 1–8 (1–185), bibliography (186–206), appendix divider (207), and Appendices A–E (208–333). The review included captions, table notes, callouts, and rotated pages. Six independent readers provided specialist and cross-version checks. The artwork inventory records all 74 caption occurrences for 72 numbered figures. Coverage is the material inspected. It does not mean every plotted point was regenerated, every paper was available, or each source was independently authenticated.

The printed contents supplied the starts below. For a section that ends as another begins, the coverage interval extends to the next heading's printed page, so shared boundary pages can appear twice. A chapter ends at its actual last page. The finer subheadings in the published paper, equations, long tables, and figure continuations are also covered by the technical reader notes and inventory.

| Section | Printed heading | Pages covered |
|---|---|---|
| 1 | Introduction | 1–10 |
| 1.1 | Motivation | 1–5 |
| 1.2 | The Biomimetic Humanoid Robot Project | 5–6 |
| 1.3 | Research Objectives | 6–7 |
| 1.4 | Dissertation Organization | 7–10 |
| 2 | Background | 11–41 |
| 2.1 | Biomimetic and Musculoskeletal Robots | 11–12 |
| 2.2 | Robot Actuation Technologies | 12–20 |
| 2.2.1 | Overview | 12–15 |
| 2.2.2 | Braided Pneumatic Actuators | 15–18 |
| 2.2.3 | BPA Force–Length–Pressure Models | 18–20 |
| 2.3 | Human Lower-Limb Biomechanics | 20–21 |
| 2.4 | Biomimetic Knee and Bipedal Robot Model | 21–26 |
| 2.5 | Neural Control of Locomotion: Central Pattern Generators | 26–28 |
| 2.6 | Sensory Feedback in Locomotion | 28–31 |
| 2.7 | Computational Neuroscience | 31–34 |
| 2.8 | Synthetic Nervous Systems | 34–35 |
| 2.9 | Neuromechanical Simulation | 35–38 |
| 2.9.1 | Simulation Tools | 36–38 |
| 2.10 | Testing Neuromechanical Hypotheses on Robots: The Ethical Case | 38–41 |
| 3 | Isometric Force Characterization of Braided Pneumatic Actuators | 42–67 |
| 3.1 | Methods | 43–48 |
| 3.1.1 | Actuator Force Model | 43–45 |
| 3.1.2 | BPA Isometric Force Characterization Experiment | 45–48 |
| 3.2 | Results | 48–55 |
| 3.2.1 | BPA Force Characterization | 48–55 |
| 3.3 | Discussion | 55–67 |
| 3.3.1 | Effects of Resting Length | 55–58 |
| 3.3.2 | Nondimensionalized Force Prediction | 58–59 |
| 3.3.3 | Comparison with Other Published BPA Models | 59–63 |
| 3.3.4 | Implications for Real-Time Control and Simulation | 63–64 |
| 3.3.5 | Biomimetic Artificial Muscle Analogy | 64–67 |
| 4 | Isometric Knee Torque | 68–122 |
| 4.1 | Methods | 69–99 |
| 4.1.1 | Robot Leg Model and Knee Test Stands | 69–74 |
| 4.1.2 | Knee Torque Data Collection and Equipment | 74–76 |
| 4.1.3 | Hybrid Torque Calculation Method | 76–77 |
| 4.1.4 | Improved Torque Model | 77–92 |
| 4.1.5 | Actuator Placement and Route Redesign | 92–93 |
| 4.1.6 | Constrained Knee-Route Redesign | 93–99 |
| 4.2 | Results | 99–111 |
| 4.2.1 | Flexor on Pin Joint | 100–105 |
| 4.2.2 | Extensor on Pin Joint | 105–107 |
| 4.2.3 | Extensor on Biomimetic Knee | 107–109 |
| 4.2.4 | Flexor on Biomimetic Knee | 109–111 |
| 4.3 | Discussion | 111–122 |
| 4.3.1 | Physical Basis of the Design Factors | 111–112 |
| 4.3.2 | Design-Factor Optimization | 112–113 |
| 4.3.3 | Redesigned Routes for the Biomimetic Robot Model | 113–116 |
| 4.3.4 | Exceeding Human Torque Magnitude | 116 |
| 4.3.5 | Tendon Reeving and Parallel BPAs | 116–119 |
| 4.3.6 | Whole-Leg Actuator and Tackle Capability | 119–122 |
| 5 | The Sensory Afferent Database | 123–135 |
| 5.1 | Curation Methods | 123–126 |
| 5.1.1 | Classifying a Feedback Experiment | 125–126 |
| 5.2 | Results | 126–134 |
| 5.2.1 | Evidence for a Phase-Dependent Feedback Rule | 127 |
| 5.2.2 | Comparing Groups of Papers | 127–129 |
| 5.2.3 | Following Citations in the Bubble Map | 129–131 |
| 5.2.4 | Coverage of the Curated Collection | 131–133 |
| 5.2.5 | AI-Assisted Literature Triage | 133–134 |
| 5.3 | Discussion | 134–135 |
| 5.3.1 | Uses of the Database | 134–135 |
| 6 | Tools for Constructing and Testing Neuromechanical Models | 136–167 |
| 6.1 | The Limits of Hand Tuning in AnimatLab | 136–137 |
| 6.2 | Moving the Neural and Mechanical Models into MuJoCo | 137–144 |
| 6.2.1 | Converting the Human Model and Comparing Muscle Force | 140–142 |
| 6.2.2 | Relocating the Actuators on the Robot Counterpart | 142–144 |
| 6.3 | Building the Mechanical Model from Robot CAD | 144–145 |
| 6.4 | SNS-Simscape Block Programming for Neural Controllers | 145–154 |
| 6.4.1 | Designing Neural Circuits in the SNS-Simscape Graphical Interface | 145–147 |
| 6.4.2 | Sensory Feedback at a Single Joint | 147–152 |
| 6.4.3 | Retesting a Controller after a CAD Change | 152–154 |
| 6.5 | Optimizing Controller Parameters in Stages | 154–165 |
| 6.5.1 | Methods for Choosing, Scoring, and Replaying Trials | 154–157 |
| 6.5.2 | Testing Bilateral Coordination with the Feet Off the Ground | 157–159 |
| 6.5.3 | Tracing a Timing Error in the Six-Synergy Controller | 159–160 |
| 6.5.4 | Contact Feedback and Locomotion | 160–163 |
| 6.5.5 | Speeding Up Repeated Controller Trials | 163–165 |
| 6.6 | Discussion | 165–167 |
| 6.6.1 | Transferring Neural Controllers Between Mechanical Models | 165–166 |
| 6.6.2 | Distinguishing Neural Drive from Mechanical Support | 166–167 |
| 6.7 | Conclusion | 167 |
| 7 | Discussion | 168–183 |
| 7.1 | Overview | 168–170 |
| 7.2 | Research Contributions and Future Work | 170–182 |
| 7.2.1 | BPA Actuator Characterization | 170–174 |
| 7.2.2 | Joint Design Factors | 174–175 |
| 7.2.3 | Sensory Afferent Database | 175–176 |
| 7.2.4 | SNS Tools and Neuromechanical Experiments | 176–182 |
| 7.3 | From Design Tools to Neuromechanical Experiments | 182–183 |
| 8 | Conclusion | 184–185 |
| A | Appendix A: Muscle Force and Joint Torque Calculations | 208–222 |
| A.1 | Human Muscle Force Calculations | 208–215 |
| A.1.1 | OpenSim Thelen2003 and Millard2012 Formulations | 208–214 |
| A.1.2 | OpenSim Moment Arms and Human Torque Targets | 214–215 |
| A.2 | Converted MuJoCo Muscle Formulation | 215–217 |
| A.3 | BPA Force Model Implementation | 217–219 |
| A.3.1 | Strain Measures | 217–218 |
| A.3.2 | Normalized Force Surface | 218–219 |
| A.3.3 | Maximum Force | 219 |
| A.4 | Muscle Path Length and Moment Arm | 219–220 |
| A.5 | Compliance Model Implementation | 220–221 |
| A.6 | Wrapping-Loss Implementation | 221–222 |
| B | Appendix B: Neuron Model Equations | 223–237 |
| B.1 | Spiking Neuron Model | 223–224 |
| B.2 | Voltage-Dependent Chemical Synapses | 224–225 |
| B.3 | Relation to the SNS-Toolbox Implementation | 225–226 |
| B.4 | Spiking-Mirror Implementation, Calibration, and Verification | 226–234 |
| B.4.1 | Adapting Integrate-and-Fire Cells and Hybrid Synapses | 226–228 |
| B.4.2 | Calibrating Transmission Between Neurons | 228 |
| B.4.3 | Rhythm and Motion With the Spiking Controller | 228–231 |
| B.4.4 | AnimatLab Conversion Variants | 231–234 |
| B.5 | Connectome of the Verified W2L Model | 234–237 |
| C | Appendix C: Model Geometry, Bracket Frames, and Compliance Parameters | 238–267 |
| C.1 | Static BPA Pathway Data | 238–245 |
| C.2 | Joint Transformation Matrices | 245–246 |
| C.3 | Bracket Reference Points and Frames | 246–247 |
| C.4 | Stiffness Arrays and Compliance | 247–249 |
| C.5 | Adopted Design Factors and Complete Pinned-Knee Results | 249–251 |
| C.6 | Redesigned Extensor Route Geometry | 251–252 |
| C.7 | Contact and Release Angles of the Wrapped Route | 252–258 |
| C.8 | Supporting Route Searches and Transmission Sweeps | 258–264 |
| C.9 | Master’s Robot Model: Joint Coordinates and Actuator Sizes | 264–267 |
| C.9.1 | Torque and Active-Force Coverage | 264–267 |
| D | Appendix D: Supplementary Neuromechanical Verification | 268–328 |
| D.1 | Adapting Mechanical Models and Connecting Signals | 268–282 |
| D.1.1 | Repairs to the Converted Models | 269–271 |
| D.1.2 | Matched-Pose Force Comparison | 271–273 |
| D.1.3 | Signals Between Neurons, Sensors, and Actuators | 273–276 |
| D.1.4 | Calculating Synaptic Conductances for a Connectome | 276–278 |
| D.1.5 | Pressure-Driven Robot Counterpart | 278–282 |
| D.2 | Circuit Implementations and Responses to Stimulation | 282–289 |
| D.2.1 | Matched Rhythm-Core Implementations | 282–284 |
| D.2.2 | Biological Knee Reflex and Mechanical Stretch | 284–286 |
| D.2.3 | A Published Neural Mechanism as a Verification Target | 286 |
| D.2.4 | Layer-Specific Controls in the Reduced Deng Circuit | 286–288 |
| D.2.5 | Li’s Walker and the Li–Deng Adaptation in AnimatLab | 288–289 |
| D.3 | Staged Calibration and Coupled-Model Experiments | 289–307 |
| D.3.1 | BiW2Lvar+ Coordination with the Feet Off the Ground | 291–293 |
| D.3.2 | Fixed-Weight Neural Noise and Phase-Specific Stimuli | 293–294 |
| D.3.3 | Locating and Testing a Timing Error in syn6 | 294–301 |
| D.3.4 | Slowing the W2L108 Walker Toward a Self-Supported Waddle | 301–302 |
| D.3.5 | Gait2392 Motion with a Pelvic Harness Compared with Human Joint Angles | 302–306 |
| D.3.6 | Balance on Fixed-Base Inverted Pendulums | 306–307 |
| D.4 | Pelvic Harnesses, Applied Loads, and Pathway Interventions | 307–316 |
| D.4.1 | Effect of Muscle Mechanics on Foot Clearance | 307–308 |
| D.4.2 | W2L Heel Strikes and Load Carried by the Pelvic Harness | 308–309 |
| D.4.3 | Pathway Deletions and Pushes during W2L Stepping with a Pelvic Harness | 309–311 |
| D.4.4 | CAD Robot Foot Clearance with a Virtual Pelvic Harness | 311–315 |
| D.4.5 | Partial Transfer of the Di Russo Controller to Gait2392 | 315–316 |
| D.5 | Spiking-Neuron Comparisons | 316–320 |
| D.5.1 | MuJoCo Neural Timing and Gait2392 Behavior with a Pelvic Harness | 317–319 |
| D.5.2 | Simulink Demonstrations | 319 |
| D.5.3 | Sensitivity of the AnimatLab Rhythm to Spiking Synapse Parameters | 319–320 |
| D.6 | Execution, Provenance, and Reproduction Tools | 320–328 |
| D.6.1 | Hardware Interface Boundary | 328 |
| E | Appendix E: Supporting SADb Records | 329–333 |
| E.1 | Record Sets and Coverage | 329–331 |
| E.2 | Afferent Classes and Linked Pathways | 331–333 |

## Printed figure inventory

The live source locations below point to the full captions.

| Figure | Printed page | Source and caption line |
|---|---:|---|
| 1.1 | 10 | 10-introduction.tex:66 |
| 2.1 | 13 | 15-background.tex:19 |
| 2.2 | 15 | 15-background.tex:39 |
| 2.3 | 17 | 15-background.tex:52 |
| 2.4 | 22 | 15-background.tex:110 |
| 2.5 | 23 | 15-background.tex:130 |
| 2.6 | 24 | 15-background.tex:139 |
| 2.7 | 25 | 15-background.tex:149 |
| 2.8 | 27 | 15-background.tex:169 |
| 2.9 | 29 | 15-background.tex:181 |
| 2.10 | 31 | 15-background.tex:204 |
| 2.11 | 32 | 15-background.tex:216 |
| 2.12 | 33 | 15-background.tex:228 |
| 2.13 | 39 | 15-background.tex:280 |
| 3.1 | 47 | 20-muscles.tex:101 |
| 3.2 | 50 | 20-muscles.tex:132 |
| 3.3 | 53 | 20-muscles.tex:298 |
| 3.4 | 54 | 20-muscles.tex:313 |
| 3.5 | 61 | 20-muscles.tex:398 |
| 4.1 | 71 | 30-torque.tex:29 |
| 4.2 | 72 | 30-torque.tex:42 |
| 4.3 | 79 | 30-torque.tex:116 |
| 4.4 | 81 | 30-torque.tex:138 |
| 4.5 | 87 | 30-torque.tex:277 |
| 4.6 | 89 | 30-torque.tex:310 |
| 4.7 | 101 | 30-torque.tex:492 |
| 4.8 | 103 | 30-torque.tex:524 |
| 4.9 | 106 | 30-torque.tex:565 |
| 4.10 | 108 | 30-torque.tex:607 |
| 4.11 | 110 | 30-torque.tex:624 |
| 4.12 | 115 | 30-torque.tex:711 |
| 4.13 | 117 | 30-torque.tex:749 |
| 4.14 | 118 | 30-torque.tex:756 |
| 4.15 | 120 | 30-torque.tex:775 |
| 5.1 | 125 | 50-database.tex:30 |
| 5.2 | 128 | 50-database.tex:58 |
| 5.3 | 129 | 50-database.tex:67 |
| 5.4 | 130 | 50-database.tex:78 |
| 5.5 | 132 | 50-database.tex:95 |
| 6.1 | 137 | 60-sns.tex:20 |
| 6.2 | 140 | 60-sns.tex:43 |
| 6.3 | 141 | 60-sns.tex:75 |
| 6.4 | 143 | 60-sns.tex:98 |
| 6.5 | 148 | 60-sns.tex:146 |
| 6.6 | 150 | 60-sns.tex:161 |
| 6.7 | 151 | 60-sns.tex:170 |
| 6.8 | 153 | 60-sns.tex:189 |
| 6.9 | 158 | 60-sns.tex:245 |
| 6.10 | 161 | 60-sns.tex:259 |
| 6.11 | 162 | 60-sns.tex:270 |
| 6.12 | 164 | 60-sns.tex:285 |
| 7.1 | 169 | 70-overall.tex:13 |
| 7.2 | 172 | 70-overall.tex:35 |
| 7.2 | 173 | 70-overall.tex:45 |
| 7.3 | 180 | 70-overall.tex:96 |
| 7.4 | 181 | 70-overall.tex:107 |
| 7.5 | 183 | 70-overall.tex:126 |
| B1 | 233 | 93-AppendixB.tex:196 |
| B2 | 236 | 93-AppendixB.tex:214 |
| B3 | 237 | 93-AppendixB.tex:223 |
| C1 | 253 | 94-AppendixC.tex:327 |
| C2 | 254 | 94-AppendixC.tex:335 |
| C3 | 261 | 94-AppendixC.tex:408 |
| C4 | 262 | 94-AppendixC.tex:418 |
| C5 | 263 | 94-AppendixC.tex:425 |
| D1 | 283 | 95-AppendixD.tex:149 |
| D2 | 290 | 95-AppendixD.tex:220 |
| D3 | 295 | 95-AppendixD.tex:265 |
| D3 | 296 | 95-AppendixD.tex:273 |
| D4 | 298 | 95-AppendixD.tex:285 |
| D5 | 305 | 95-AppendixD.tex:356 |
| D6 | 313 | 95-AppendixD.tex:408 |
| D7 | 318 | 95-AppendixD.tex:463 |
| D8 | 321 | 95-AppendixD.tex:503 |

## Comparison after the independent assessment

**Order of evidence.** The independent text and figure assessment above was saved before either ZCode file was opened. Its first 89,661 bytes have SHA-256 `9891cfaa1ee50ccae5710794cdbda08d53eae380069e8421d7cd4a5db8cb0862`. This addendum uses the November 22 easteregg2 computational report and run manifest delivered to the EB475WS4 handoff. Their SHA-256 values are `aede2e2e4386257d470696841c105ff7d459ab0e19e78fb486d864ee844d0f1f` and `e7c8ddbafbab4082b1b1af6bc70995477e5d6bb041a3589e93be0def827706fe`. ZCode performed the heavy computations. The observations below identify which conclusions EB475WS4 could verify from independent local records.

| Question | Evidence after comparison | Consequence for the live dissertation |
|---|---|---|
| Which edition was reviewed? | Both machines found 353 corresponding pages with identical live/working text and the same 17 text-difference pages against the November 20 reconstruction. The EB475WS4 comparison additionally inspected raster differences and embedded art. | The earlier stale-source warning does not apply to the finalized abstract or acknowledgments in the reviewed handoff. |
| Why do two W2L baselines give different strike counts and progression? | **Resolved from saved jobs.** In `Code/MuJoCo_SNS/spinal/reports_20261027_codex/trial458_d30.json` and `panel_control.json`, 69 of 70 saved settings agree. The remaining parameter, `wsup`, was 0.70 in §D.4.2 and 0.25 in §D.4.3. This is prescribed feedforward lift as a fraction of body weight. The net harness force also includes the spring and damper. Both runs used the +40 N forward pulse from 8 to 8.15 s (`w2l_mujoco/test_w2l_ground.py`, lines 148 and 636–640). The profiles give 28 L/21 R, phase 0.496, 2.246 m versus 34 L/22 R, phase 0.498, 2.239 m. The available panel NPZ also reproduces the combined-deletion result. | At pp. 308–311, 95-AppendixD.tex:381–395, add one sentence at the boundary: “The pathway-intervention series prescribed 25% body-weight lift, compared with 70% for the preceding replay. Both included the same forward pulse.” Preserve each run's numbers. Describe the scoring criterion as median phase, as in the independent finding. |
| Are the two timing measurements in Table 6.4 traceable? | **Resolved for the reported values.** The original benchmark bundle was unavailable to ZCode. On EB475WS4, `Documentation/Reports and Papers/Dissertation/Notes/ch6_authoritative_20261116/TIMING_RESULTS.json` (SHA-256 `d500e760092dce03439d95a6cb5be2bb6b6999979f1e64ed1a77a210d8eada45`) gives external wall times 199.6235726 and 29.2842800 s (lines 73, 348) and the corresponding neural-update times 179.0791121 and 9.7795164 s (lines 11, 286). It records matching compiled model and aligned-trace hashes. The November 4 defense evidence inventory independently points to the producing campaign. | Keep Table 6.4, p. 164, 60-sns.tex:295–307. The 6.82× complete-process and 18.31× neural speedups apply to one matched 1,025-cell trial. The additional 1,033-cell benchmark in the ZCode manifest appears separately in the live Appendix D. The evidence supports the measurements as recorded and permits tracing the original workload. No new wall-clock timing was collected here. |
| Did the selected human controller reach both stance and joint-angle targets? | ZCode recomputed 531 physical evaluations, of which 428 were eligible. None met the joint-angle target and both stance bounds together. All nine recorded finalist-by-duration cells reproduced. Of eight perturbations of the stance-selected case, eight retained at least six loading cycles, one met both stance bounds, and zero met both stance and angle targets. Earlier stored and physically rescored studies also differed in effective load threshold. | This strengthens the independent finding on pp. 142 and 302–305. Keep the selected interval, transition back toward standing, and 50 N loading criterion visible wherever the short-trial stance result is interpreted. |
| Do the pressure-driven robot and the human-model controller results reproduce? | ZCode's selected pressure-driven trial gives +1.05 cm **net** fore–aft displacement over 16 s, despite a 9.30 cm within-trial excursion. Eight of nine fresh pressure-curriculum replays reached full duration. The baseline stage-3 run exited at 4.1 s when a sleeve left the model's force domain. The rank-six compliant recruitment fit reproduced at 84.77% right and 84.31% left centered variance accounted for. | Preserve “about 1.1 cm” as net progression for the pressure-driven counterpart (pp. 160–163, 278–281), and retain the limit imposed by the actuator model. The fit in §D.1.4, 95-AppendixD.tex:115, describes the compliant-recruitment configuration. The separate native-muscle fit has different ranks and percentages. These simulation pressure limits are distinct from the physical Chapter 4 flexor pressure readings. |
| Do Tables 4.12 and C14 support the whole-leg design summary? | ZCode's angular-domain quadrature reproduced the saved record bit for bit. Mean torque/active-force coverage was 37.3527%/92.9977% at baseline and 65.8091%/99.3056% in the selected design. All 27 site rows agreed at printed precision. | Keep the published rounded values and specify modeled coverage over each coordinate domain. These tables make a different claim from the 15 measured flexor torques. |

**Extent of the comparison.** ZCode identified the November 21 knee/arm, muscle-force, four-muscle retuning, CAD, and selected robot replays as prior evidence carried forward. Those are not additional November 22 replications. The delivered computational summary expressly withholds the separate full text review. The new W2L and timing findings above resolve its two stated record questions using local evidence.
