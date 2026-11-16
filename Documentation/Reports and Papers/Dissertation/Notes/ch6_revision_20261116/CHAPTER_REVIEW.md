# Chapter 6 revision and evidence record

Source edited: Final/chapters/60-sns.tex. Baseline: latest downloaded Bolen_Dissertation.zip, extracted without modification in tmp/ch6_20261116/baseline. Other Final chapters and bibliography were not edited by this task. User-authorized scope is the complete Chapter 6/associated Appendix D rewrite, removal of irrelevant material, current figures/results, and independent supervision.

## Structure and directive compliance

| Directive | Implementation evidence in 60-sns.tex |
|---|---|
| Align with rewritten introduction objective 4 | Title and introduction, lines 1–26; contribution explicitly attributed to implementations, integrations, and experiment procedures |
| Coherent dissertation chapter | Introduction at line 4, Methods at 28, Results at 148, Discussion at 253, Conclusion at 279 |
| Distinguish design, calibration, and physical curriculum | Methods subsection line 108 and physical-curriculum table; one Hunt-method explanation replaces duplicated explanations |
| Establish model identities | Model table and explanatory paragraphs, lines 43–61; no controller name is used as a mechanical-body identity |
| Connect construction with testing | Workflow and capability table, lines 30–79; mechanical checks, signal contracts, protocols, and corresponding Results |
| Update current results | BiW2Lvar+ air/noise/stimulus at lines 180–195; PF failure at 208–212; newer W2L108 and CAD at 214–229; runtime at 231–251 |
| Interpret evidence and remaining requirements | Three Discussion subsections, lines 255–277 |
| Preserve technical details without inventory bloat | Detailed configurations/scripts and secondary tests referenced to specific Appendix D anchors |
| Protected material | No published muscle chapter, torque sources, torque results, or unrelated chapter prose changed |
| Source verification | Current numerical evidence and changes documented below; independent supervisor/full build handled by root |

Line numbers describe the first complete draft and may move during supervisor corrections.

## Numerical additions, replacements, and interpretation changes

| Baseline material | Revised material | Evidence and scope |
|---|---|---|
| Unspecified trained human candidate had passed air stages and assisted standing | Current revised **BiW2Lvar+** passes only same-weights air stages; no standing/ground inheritance | Defense/Session1_Handoff_20261104/experiments_20261104/EXPERIMENT_UPDATE.md and GAIT_AIR_VALIDATION.json. Two 16 s runs, expression 0.1. Principal cycles R/L 4/4 and 7/7; periods 2.466/1.300 s; antiphase errors 0.0337/0.0933 cycles. Input/command conditions may differ, so feedback alone is not assigned causal responsibility for the period change. |
| No current body-noise example | Six continuous RG-current noise trials pass air-stage movement/knee/separation checks | Same folder, BiW2LvarX_RG_NOISE_VERIFICATION.json and EXPERIMENT_UPDATE.md. Three seeds at 0.01/0.025 nA stationary SD, correlation 20 ms; independent repeated seed exactly matches. Not fixed-weight perturbation or ground-balance evidence. |
| Older large-current reduced Deng layer tests dominated main text | Current 2 nA, 100 ms phase-specific BiW2Lvar+ perturbations summarized; reduced tests retained in Appendix | Same folder, BiW2LvarX_STIMULUS_PAIRED_METRICS.json and EXPERIMENT_UPDATE.md. Active PF changes up to 1.14 mV and iliacus MN 0.16 mV; silent PF pulse remains below transmission threshold. |
| No PF design failure reported | Syn6 0/12 timing-channel pass, actual centered-variance components 3/side versus intended 6; s3k 2/side and undeclared timing targets | Same folder, circuit_trace_20261105/PF_DESIGN_AUDIT_20261105.md. Bounded repair remained unsuccessful; no new accepted controller. Prior output-fit passes cannot establish correct PF timing. |
| Newer W2L displayed 9L/8R swings and 40.5% whole-trial net assistance | Newer bounded-fit W2L108 described with 8L/7R swings, 29.519% upward assistance after 1 s and 35.069% whole-trial signed net assistance | Same folder, LiDeng_W2L108_RESULT.md. Separate condition, not a numerical correction to the older trial. 108 cells/216 connections; 10 s; knee flexion minima 1.946/2.087 degrees. |
| No slow-shuffle failure | Slow candidate has 7/7 swings and approximately 10.1% assistance over 2–16 s, but failed excursion/double-support gates | Same folder, SHORT_SHUFFLE_RECORD_20261104.md. Not promoted as accepted short shuffle. |
| Old CAD abductor-control still, 80% BW, 67%-to-zero reversed ankle ordering, provisional welded foot mount | Current corrected-foot/contact-handover condition: 77.5% measured assistance during 16–28 s, 7R/6L qualifying swings, 5.10-degree knee margin, 5/8 local checks | Defense/Session1_Handoff_20261104/records/CAD_CONTACT_HANDOVER_EXPERIMENT_RECORD.md and CAD_RESULTS_AND_FILES.md. Selected run 28 s; SNS timing plus engineering allocator and 54 ideal reference-tension actuators. Old comparison is not silently reinterpreted as current model behavior. |
| Runtime table led 9.21× balance/17.75× suspended workload and assorted 56–584× kernels | Matched Gait full-process 199.624/29.284 s =6.82×; neural 179.079/9.780 s =18.31×; later sparse/compact kernel .895287/.555062 s=1.613×, whole diagnostic 92.944/91.516 s=1.0156× | Defense/Session1_Handoff_20261104/docs/TOOLCHAIN_EVIDENCE_INVENTORY.md §7 and records/NEURAL_ACCELERATION_20261104.md. Different workloads/scopes identified. No new benchmarking. 62 arrays/model match exactly in compact pair; full diagnostic timing is one pair, warm kernel median of three. |
| Cup sag 2.00 versus 7.85 degrees without endpoint detail | Same numbers, explicitly 10 s ramp endpoint reaching 0.625 kg, not equilibrium | Graphics owner's source-builder audit. No numerical response recalculation by prose task. |
| Knee described as generic antagonist BPA example | Same 1.09-degree peak, 111 ms, recovery 6.355 s; explicitly simplified linear BPA-force block, not full Festo model | Baseline and Code/Matlab/SNS_Simscape/results/knee_transient_reflex_20261028; root/graphics source audit |

Unchanged main results retained: matched-pose force exceptions and eight-of-twelve scope, W2L rhythm-core 1.000 s/ten bursts/correlation ≥0.99996, historical fixed-controller feedback phases .498/.223 and displacements 2.239/1.996 m over 4.5–30 s, bridge readback and stated exported-network test. None is merged with the new controller's stage results.

Historical original AnimatLab approximately 2.2 Hz figure is summarized by the suspended/fall contrast; aligned recording and source context remain in Appendix D. Isolated compliant-bank, Di Russo equation-kernel, and other secondary timing claims are omitted from the main argument rather than reinterpreted as whole-body performance.

## Figure dispositions

No scientific source graphic was deleted from disk by this task.

| Baseline asset/label | Disposition |
|---|---|
| OurWalker.png / LI_DENG_FALL_FRAME.png | Retained together in short motivating figure |
| connor_ben_shared_rg_original.png / fig:animatlab_shared_rg | Root asked Appendix owner to retain in Appendix D |
| Gait2392_human_reference_oblique.png / fig:gait2392_human_geometry | Static model inventory view removed from main; model identity explained in table and current coupled example |
| w2l_selected_left_subnetwork.pdf / fig:w2l108_subnetworks | Main inventory schematic omitted; Appendix contains relevant pathway/intervention evidence |
| DENG_RAT_NATIVE_SUBNETWORKS.pdf / fig:rat_native_subnetworks | Main native-page inventory omitted; unrelated to central verified experiment sequence |
| W2LNEW_TREADMILL_STILL.png / fig:walker_fall_to_supported | Old supported trial still omitted; current W2L108 result described with exact condition and Appendix source |
| simulation_workflow.pdf | Replaced by root's simulation_workflow_20261116.pdf |
| sns_library_clean.pdf | Replaced by sns_library_20261116.pdf; current catalog only, no invented internal panel |
| knee_transient_current_layout.pdf | Replaced by sns_knee_reflex_20261116.pdf, actual six paths and current/angle/MN signals |
| beverage_current_layout.pdf | Replaced by sns_cup_reflex_20261116.pdf, actual circuit/geometry/load/elevation/command panels |
| ah001_supported_control.png | Replaced by current cad_supported_20261116.png; label retained |
| w2l_feedback_short_20261031.pdf | Retained without numerical alteration |
| New biw2lvar_air_stage2_20261116.pdf | Root's one-page selected neural/mechanical panel; caption states retained positive-flexion display and actual source sign |

## Checks and remaining verification

- 18 begin/end environment pairs balanced in the completed source.
- All cited keys exist in the latest ZIP's thesis.bib.
- No bibliography rewrite, simulation, optimization, commit, push, or Overleaf interaction.
- No dates or machine/campaign narration added to chapter prose; Oxford commas and model terminology retained.
- Root owns full compilation, cross-reference/asset checks, rendered-page inspection, and independent supervisor corrections.
- The seven-step GLIF procedure was not presented as a universal neural-design method. The baseline bibliography lacks the specific Szczecinski 2020 DOI 10.3389/fnbot.2020.577804; root was notified if a named procedure/reference is wanted.

