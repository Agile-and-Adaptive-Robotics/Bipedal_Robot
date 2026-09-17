# Dissertation review working notes — one entry per ajh26 item, problem + response together

**File map**
- Verbatim archive of all 66 ajh26 items (ChatGPT sweep, with full highlighted text): `summary_inventory_ajh26.md` in this folder.
- This file: every item with its response in the same entry. Status per entry: DONE (accepted tracked change), PROPOSED (suggested text below, nothing applied), OPEN (needs Ben's call).
- Suggested full-file drafts: `ZCode_drafts/chapters/` (edited .tex files with all PROPOSED text in place).
- Tracked changes: 16 of 17 already accepted in Overleaf. The one pending is noted at B16.

**Incident note (2026-09-12):** ZCode mistakenly pasted drafts into Overleaf; project restored from History (version before 02:40). No Overleaf edits without an explicit instruction naming the action.

---

# Handoff for ChatGPT (2026-09-16)

You (ChatGPT) made the original 66-item review inventory and summary.md. ZCode has since restructured this file: every ajh26 item now has problem + response in the same entry below. Read this section first, then the entries.

**Current state**
- Overleaf was restored (History, version before Sep 12 02:40) after ZCode mistakenly pasted drafts directly. Authoritative text = Overleaf History; the repo zip `Bolen_Dissertation.zip` is only the Sep 12 00:37 snapshot. Ben's editor mode dropdown may read "Edit" (ZCode switched it; he can switch back to Reviewing).
- Review panel: 49 ajh26 comments unresolved; 16 of 17 tracked changes accepted; one tracked deletion still pending ("spinal ", background L42 — recommendation at entry B14 is to reject it).
- This file: working doc, one entry per item, statuses DONE / PROPOSED / OPEN. Full verbatim inventory (your original): `summary_inventory_ajh26.md` beside this file. Full-file drafts with every PROPOSED text in place: `ZCode_drafts/chapters/`.

**Ground rules (Ben's, enforced after the incident)**
1. NEVER modify Overleaf — not files, not tracked changes, not comment resolves — unless Ben names the exact action in that message. Suggested edits go into this file as PROPOSED entries. Ben applies everything himself.
2. Problem and solution in the same entry; never make Ben scroll to match them.
3. Every entry you write or revise carries a score marker (add score / minus score) and, where possible, the specific words that earned it — the labeled example is the part that teaches.
4. Style rules now standing: no em-dashes in new prose (split into two sentences instead); the last sentence of a section bridges to the next and states the dissertation's stake; completed-work framing only ("these tools let researchers do X", never "planned/future work"); quantitative over qualitative.
5. FLAG-BEN placeholders mark values only Ben can supply (human knee torque results). Never fill them by guessing.
6. Append or edit entries in place; do not renumber, reorder, or rewrite other sections' entries. Note what you changed at the bottom of the ledger.

**Work available for you now**
- Draft the remaining OPEN prose items that do not need Ben first: the X1 motivation expansion (2-3 pages, scientific framing) and the full-paragraph Dissertation Organization variant ajh26 asked for (so Ben can compare against the summary-sentence version in I5).
- The B- global grade targets a voice-consistency pass: read the merged sections (background 2.6-2.8 after B15's distribution, results R1) and flag seams where source papers show through, as new entries.
- After ajh26's next review pass: refresh the inventory from Overleaf's Review panel (read-only — your browser access is the designated route), append new items as entries, and update resolved ones.
- The AI-writing-check rubric (learning log item 9) is unclaimed: em-dash density, hedge stacking, uniform sentence rhythm, AI lexicon, rule-of-three overuse. If you build the checker, run it on all PROPOSED text in this file and log findings as entry updates.

**Division of labor with ZCode:** you handle prose and the Overleaf review panel (read-only); ZCode handles the code side (pulling the Newton-Raphson stop conditions from `minimizeFlxPin.m` for M5, equation cross-checks for M4, the MATLAB evaluators) and maintains `ZCode_drafts/`. Do not edit `ZCode_drafts/` — your text goes into entries here.

---

## 02-abstract.tex

### A1 — one page, move extras out — PROPOSED
ajh26: "Pull your abstract down to one page. 'Extra' material can be moved to the introduction and 'dissertation organization' section." + Ben's spec: tools-for-biomimetic-robots framing, brief walkthrough of completed work, quantitative, future = what researchers can now do.

```latex
Biomimetic robots whose legs reproduce human muscle geometry and joint torque would allow researchers to run neuromechanical experiments that are impractical or unethical in living subjects, but designing them requires tools that predict artificial-muscle force, joint torque, and neural control before hardware exists. This dissertation develops and validates those tools for the humanoid robot under construction at Portland State University, a robot actuated from the trunk down by Festo braided pneumatic actuators (BPAs) and controlled by a synthetic nervous system.

The isometric force of BPAs with \qtylist{10;20}{\mm} uninflated diameters was characterized across resting lengths, contractions, and pressures on a custom test jig. The characterization produced two models: a maximum-force law at \qty{620}{\kPa} whose resting-length dependence was previously unreported and is predicted with the wrong trend by the manufacturer's tool, and a three-coefficient normalized force surface (adjusted $R^2 > 0.99$) that corrects the large over-predictions of existing models at short resting lengths.

The actuator model was then carried to the joint through isometric knee torque experiments on a pinned knee and on a biomimetic four-bar knee with a migrating instantaneous center of rotation. A hybrid torque calculation method isolated three loss mechanisms: constant length offsets, series compliance in brackets and tendons, and loss of usable actuator length where BPAs wrap the joint. Correction terms identified by multiobjective optimization improved prediction on every validation test, with per-test error falling from 1.0--3.6 to 0.6--1.9~N{\,}m across flexor and extensor configurations. With actuator placements redesigned under the corrected model, the biomimetic knee met or exceeded the human isometric torque envelope over \textbf{[FLAG-BEN: \% of RoM]} of the range of motion (RMSE \textbf{[FLAG-BEN: value]}~N{\,}m against the human target).

Finally, the validated plant was connected to biologically grounded control. The laboratory's modified OpenSim leg model was converted to MuJoCo, and a synthetic nervous system implementing a two-layer central pattern generator, with per-leg rhythm generation and pattern formation, left--right coordination, and Ia, Ib, II, and heel/toe contact feedback, was verified against published tonic-stimulation and deletion phenomena. With these tools, researchers can now size artificial muscles and their routing to human torque specifications, validate the resulting joints against measurement, and run lesion and stimulation experiments on the simulated humanoid leg in preparation for a treadmill robot with embedded real-time neural control.
```

Every sentence states a completed deliverable; the 1.0--3.6 → 0.6--1.9 N m range matches sec:ongoing's own numbers; FLAG-BEN marks the pending human-knee values. **minus score** (original needed the one-page comment; "defines the path" was the uncompleted-task vibe Ben flagged).

### A2 — tracked deletion of "Planned extensions..." sentence — DONE (accepted). The rewrite also keeps future work out of the abstract.

---

## 10-introduction.tex

Tracked changes L21/L23 (cannot→do not accurately; delete "maximum"; "BPAs"→"BPA force"; delete "isometric test system is compliant"; ". Addtionally,"; add "or compared"): all **DONE (accepted)**. One landed typo to fix: "Addtionally," → "Additionally," — **PROPOSED, add score** (caught a typo the accepted change introduced).

### I1 — L7: gaps in knowledge — PROPOSED
ajh26: "Put in a clear sentence here describing the gap(s) in knowledge." Built from Ben's raw material. Insert between the `\citep{hunt_modeling_2017, sarosi_comparative_2017, liang_comparative_2020}.` sentence and "Designing a biomimetic robot...":

```latex
The author's experience with \qty{10}{\mm} actuators made this gap concrete: manufacturer specifications overstate the force and contraction available at the short resting lengths a leg design requires, published models fitted at one cut length do not transfer to other resting lengths, and designs carried out against those specifications delivered less range of motion and joint torque than intended, forcing iteration on hardware rather than at the design stage.
```

Covers all four of Ben's bullets (specs overstate; resting-length relation unknown; single-length models do not generalize; hardware-iteration consequence). **minus score** (gap was left implicit).

### I2 — L15: "Two..." paragraph + Morrow 2020 mention are PhD work — PROPOSED (move to Results)
ajh26: "This is work you did during the PhD... moved to a chapter or the discussions section." Per Ben: to Results.

Two changes: (1) the L13 sentence ending "...it motivated a multiobjective optimization ... \citep{morrow_optimization_2020}." becomes:

```latex
That work showed that simply transplanting human muscle attachment locations onto a BPA-actuated robot does not reproduce human torque profiles; the placement optimization that closes this gap was completed during this dissertation and is summarized with the results in Section~\ref{sec:foundation}.
```

(2) The whole "Two further strands of preliminary work..." paragraph is deleted from the intro and becomes the new `sec:foundation` opening of the Results chapter (full text in R1 below). **minus score** (draft framed PhD work as pre-dissertation preliminary work).

### I3 — L21: "blocks rational actuator selection" reword — PROPOSED (previously missed)
ajh26: "reword". Suggested:

```latex
Because muscle resting length is one of the primary degrees of freedom available to the robot designer, this gap leaves actuator selection to trial and error rather than calculation.
```

**minus score** (flagged wording left standing).

### I4 — L33: "a whole paper... just a section" — OPEN (previously missed)
ajh26: "This is a whole paper, seems weird it is just a section while the previous paper is a chapter and a half." With the new foundation-studies section expanding Chapter 3, the imbalance partially resolves, but the honest options are: (a) promote the torque work to its own results chapter, or (b) keep one Results chapter and let the expanded sections carry the weight. Ben's call; affects Dissertation Organization wording (I5).

### I5 — L39: paragraph per chapter, combine with Research Objectives — PROPOSED
ajh26: "Write a whole paragraph for each chapter. Could probably be combined with the Research Objectives section." Ben's lighter spec: summary sentence per content block. Suggested replacement paragraph (also removes the "planned... eventual treadmill robot" phrasing Ben flagged):

```latex
Chapter~\ref{ch:background} reviews the relevant background: biomimetic and musculoskeletal robots, artificial muscle technologies and BPA modeling, human lower-limb biomechanics and the OpenSim benchmark models, the neural control of locomotion with an emphasis on two-level CPG architectures and sensory afferent feedback, synthetic nervous systems and the Sensory Afferent Database, and the simulation tooling that connects OpenSim models to MuJoCo. Chapter~\ref{ch:methods} presents the materials and methods: the actuator force model and its characterization experiment, the pinned-knee and biomimetic knee test stands, the hybrid torque calculation method that isolates the loss mechanisms, the correction terms fitted to them, and the optimization algorithms used to redesign actuator placements. Chapter~\ref{ch:results} first summarizes the completed foundation studies (the placement optimization, the Sensory Afferent Database, the BPA balance platform, and the BPA pressure-dynamics model), then reports the actuator characterization, the identification of the torque correction terms on the pinned and biomimetic knees, the redesigned actuator routes, and the experimental comparison of the biomimetic knee with human isometric torque. Chapter~\ref{ch:discussion} discusses the findings, their limitations, and their implications for real-time control. Chapter~\ref{ch:futurework} assembles the validated plant and the neural controller into a working neuromechanical simulation pipeline, an OpenSim-to-MuJoCo conversion coupled to a two-layer CPG with full sensory feedback; with these tools in place, researchers can run lesion, stimulation, and gait-modulation experiments on the simulated humanoid leg today, and the same tools define the path to a treadmill robot with embedded real-time control. Chapter~\ref{ch:conclusion} concludes.
```

**minus score** for the "planned" phrasing; **add score** for the reframe. OPEN: whether ajh26 still wants the full-paragraph version merged into Research Objectives.

---

## 15-background.tex

### B1 — L1: borrow cited figures — PROPOSED (suggestions)
ajh26: fine to borrow images for a background section if cited. Candidates: Kengoro/Kenshiro annotated musculoskeletal figure (asano_design_2017 Fig. 1) in 2.1; BPA cutaway (chou_measurement_1996 or hunt_modeling_2017) in 2.3; force-length comparison overlay (liang_comparative_2020) near eq:BPAFit in 2.4; Rybak two-level CPG schematic (rybak_modelling_2006) at the top of 2.6 — the architecture the dissertation implements; SNS-Toolbox network example (nourse_sns_2023) in 2.7. **add score** (cheap wins).

### B2 — L5: connect each section to your work — PROPOSED
ajh26: add a sentence at the end of each section/paragraph connecting it to the PhD work, e.g. "The work in this dissertation will help develop biomimetic robots that more accurately captures animal joint torque profiles." End of sec 2.1:

```latex
The design and validation tools developed in this dissertation aim at that second shortfall: they predict, and verify against measurement, the joint torques a BPA-actuated leg can actually deliver.
```

Varied per section (B3, B7 entries, B8, and the SOTA-distribution closers) so it never reads formulaic. The L32 closer ajh26 praised ("Yes, like this!") stays untouched. **add score** (applying the praised pattern).

### B3 — L9, four comments (Compliant doesn't follow / Recent doesn't follow / plural? / connect to your work) — PROPOSED
Ben's a/b/c plus the tie-in. Replace the paragraph's last three sentences with:

```latex
Dielectric elastomer actuators (DEAs), a class of electronic electroactive polymers, contract by a different mechanism: an applied voltage squeezes a soft polymer film between compliant electrodes through electrostatic stress, and the thinning and radial spreading of the film produces strains that can exceed biological muscle \citep{madden_artificial_2004}. Comparative studies report high specific work and good bandwidth for DEAs \citep{liang_comparative_2020}, but kilovolt-range driving voltages, electromechanical instability, and the compliant-electrode requirement have kept them out of legged robots. Closer to muscle, series-elastic and pneumatic actuators trade peak force and bandwidth for compliance, a trade that wearable and rehabilitation robotics exploit because force control and backdrivability matter most when a machine interacts with a human \citep{aksoz_design_2019, leibach_development_2020, li_influence_2020}. A recent review of pneumatic artificial muscle technology documents continued progress in actuator design and in assistive applications \citep{zhagiparova_recent_2025}. None of these technologies yet combines muscle-like force, strain, and compliance in one actuator, which is the combination a biomimetic humanoid leg demands and that braided pneumatic actuators approximate most closely.
```

DEA contraction mechanism now explained (Ben's a); "reviews" singular with its single citation (b); every sentence follows from the previous one; final sentence bridges to the BPA section and states the dissertation's stake. **minus score** (three flow complaints in one paragraph).

### B4 — L15: end-of-section tie-in — PROPOSED
End of sec 2.2:

```latex
Predicting their force as a function of length and pressure is therefore the first requirement of the torque tools developed in this dissertation, and the state of that prediction is the subject of the next section.
```

**add score**.

### B5 + B6 — L24 equation formatting + L30 "there is no e* above" — PROPOSED
Replace the adjustwidth/multline and the "where $S$..." sentence with a symbol-form equation, a coefficient table, the corrected epsilon sentence, and Ben's lookup-table statement:

```latex
\begin{equation}
P = P_{0} + k_{F}\,F + k_{S}\,S + P_{t}\tan\!\left(c_{a}\left[\frac{\epsilon}{\epsilon_{\max} - k_{E}\,F} - c_{b}\right]\right) \label{eq:BPAFit}
\end{equation}

\noindent
where $F$ is the required force, $\epsilon$ is the contraction, $S$ encodes hysteresis ($S=1$ shortening, $S=-1$ lengthening, $S=0$ static), and the seven coefficients are defined in Table~\ref{tab:BPAFit_coeff}. The coefficients have been updated with corrected values, as those reported in \citet{hunt_modeling_2017} contained typographical errors. The contraction $\epsilon$ and the maximum free-load contraction at \qty{620}{\kPa}, which appears as $\epsilon_{\max}$ in Equation~(\ref{eq:BPAFit}), are defined in Section~\ref{sec:force_model}, together with the relative contraction $\epsilon^{*} = \epsilon/\epsilon_{620}$ on which the force surface developed in this dissertation is built. Because Equation~(\ref{eq:BPAFit}) expresses pressure in closed form as a function of force and contraction, it can also be solved onto a grid of pressures and contractions to build a force--pressure--relative-contraction lookup table for design studies.

\begin{table}[htbp]
\centering
\caption{Coefficients of the \qty{10}{\mm} Festo BPA pressure model of Equation~(\ref{eq:BPAFit}), corrected from \citet{hunt_modeling_2017}.}
\label{tab:BPAFit_coeff}
\begin{tabular}{lll}
\hline
Symbol & Value & Role \\
\hline
$P_{0}$ & \qty{254}{\kPa} & pressure offset \\
$k_{F}$ & \qty{1.23}{\kPa\per\newton} & linear force-to-pressure gain \\
$k_{S}$ & \qty{15.6}{\kPa} & hysteresis offset (sign of $S$) \\
$P_{t}$ & \qty{192}{\kPa} & tangent-term amplitude \\
$c_{a}$ & 2.03 & contraction-shape coefficient \\
$c_{b}$ & 0.46 & contraction-shift coefficient \\
$k_{E}$ & \qty{0.331e-3}{\per\newton} & force-dependent reduction of usable contraction \\
\hline
\end{tabular}
\end{table}
```

**minus score** (the e* reference was factually wrong as written).

### B7 — L32: "Yes, like this!" — no change
Praise; the sentence is the template the other tie-ins copy. **add score**.

### B8 — L36: "What is opensim?" — PROPOSED (previously missed)
Expand at first mention:

```latex
This dissertation uses OpenSim, an open-source platform for building and simulating musculoskeletal models and analyzing their muscle moment arms and joint torques \citep{delp_opensim_2007, seth_opensim_2011, seth_opensim_2018}, together with its Gait2392 model, which...
```

**minus score** (unexpanded tool name on first use).

### B9 — L38: Steele knee paragraph + figure — PROPOSED
ajh26: "Turn this into a whole paragraph with figure." Replacement paragraph:

```latex
Steele's biomimetic knee joint supplies the physical knee used in this dissertation \citep{steele_development_2017, steele_biomimetic_2018, steele_biomimetic_2018-1}. It is a one-degree-of-freedom four-bar linkage whose crossed links reproduce the posterior migration of the human knee's instantaneous center of rotation (ICR) during flexion, so the rolling-and-sliding behavior of the human joint emerges from the link geometry itself rather than from cams or gear profiles. Because the ICR migrates, muscle moment arms about the knee change with flexion angle in a human-like way, which is exactly the behavior a torque prediction tool must capture. The biomimetic knee tests of Chapters~\ref{ch:methods} and~\ref{ch:results} therefore exercise the corrected torque model against this geometry, and Figure~\ref{fig:steele_knee} summarizes the mechanism.
```

Figure (Ben assembles; draft uses a compile-safe fbox placeholder): **(A)** four-bar schematic with crossed links and ICR path (adapt from steele_development_2017, cited); **(B)** ICR migration vs flexion angle, Steele vs published human knee data; **(C)** photo of the physical knee on the test stand. **minus score** (one sentence for the joint the whole dissertation tests).

### B10–B13 — L50, four comments (steering/balance confusion; not AARL x2; rat stuff) — PROPOSED
ajh26: "wasn't this balance?"; "Not AARL, but Quinn lab"; "not the AARL"; "How about the rat stuff from myself, Deng, Young, and Jackson?" Replacement:

```latex
The AARL has pursued this approach across several platforms: CPG-based neural control of a pneumatically actuated dog robot \citep{hunt_development_2017}, analytically designed dynamic neural networks for human balance control \citep{hilts_dynamic_2019}, the rat hindlimb modeling program that supplies this dissertation's controller architecture (a full-muscle rat hindlimb model and its moment-arm analysis \citep{young_analyzing_2019}, and the two-layer CPG neuromechanical walker built on it \citep{deng_neuromechanical_2019}), and the DoggyDeux quadruped designed for SNN control \citep{scharzenberger_design_2019}. Related controllers from the wider neuromechanical community include sensory-entrained CPG control of bipedal walking \citep{li_bipedal_2017}, inter-leg coordination analyses in stick-insect-inspired controllers \citep{nourse_analyzing_2019}, reinforcement-reflex hybrids in insects \citep{goldsmith_investigating_2021}, and adaptive hindlimb walking controllers \citep{schilling_adaptive_2022}.
```

The rat program enters with Young and Deng (both already in thesis.bib; ajh26 coauthor on both; "Jackson" = Axel Jackson, coauthor on nourse_analyzing_2019 in the neighboring sentence). Non-AARL work split into its own sentence. **OPEN for Ben:** li_neural_2016 is literally a steering paper; the bipedal one is li_bipedal_2017 ("...entrained by sensory feedback controls walking of a bipedal model"), also already in the bib — the draft swaps to li_bipedal_2017 per Ben's "rename to bipedal locomotion." Confirm. **minus score** (three attribution errors in one sentence).

### B14 — L42: tracked deletion of "spinal " — PENDING, recommend REJECT (previously unflagged)
The only unaccepted tracked change. "Locomotion in mammals is driven by spinal central pattern generators" is anatomically correct; deleting "spinal" weakens it. Ben decides in the panel.

### B15 — L58 + L60: dissolve the State-of-the-Art section — PROPOSED
ajh26: "Most of this information is elsewhere... I don't really like this section." Ben: "He's right. Find how to break it up and distribute it." Also fixes the content-before-overview ordering Ben called embarrassing. Distribution (every citation preserved, no dangling refs):
- Geyer muscle-reflex model + Haeufle hybrid → end of the sensory-feedback paragraph in sec 2.6, closing with: "The controller of this dissertation follows the CPG tradition while retaining the full afferent set of the reflex tradition, so its predictions can be tested against both bodies of work."
- Song deep-RL + degroote perspective + the hierarchical-CPG list (ichimura, nishizaki, fukuoka) + simulation-speed enabler → end of sec 2.8, closing with: "That throughput turns lesion and stimulation studies from bespoke demonstrations into statistically powered experiments, and it is what the tools of this dissertation are built to deliver."
- L60's "plant"→"model" and the "who has advanced this" rephrase are absorbed by the rewrite (their section is deleted).
Full replacement paragraphs are in `ZCode_drafts/chapters/15-background.tex`. **add score** (strongest structural dislike resolved by distribution, not deletion).

### B16 — L70: pull back future-work framing — DONE (accepted will→can) + rule logged
The comment's larger point (do not muddy completed vs planned; frame benefits as-is) is now a standing style rule (see Learning log) and shapes every PROPOSED text above. No further edit needed for this sentence unless Ben wants more distance.

---

## 20-methods.tex

### M1 — L10: wrapping also designs attachments — PROPOSED (previously missed)
ajh26: "Also used to design attachment locations for a robotic knee that meets the torque capabilities of a human knee." Add to the Overview where wrapping/BPA routing is introduced:

```latex
The same torque prediction serves both roles in this dissertation: it corrects predictions against the test-stand measurements, and it drives the optimization that selects BPA attachment locations for a robotic knee matching human torque capabilities (Section~\ref{sec:optimization}).
```

**minus score**.

### M2 — L173: x and z are not independent — PROPOSED (previously missed)
ajh26: "Make sure you discuss why. Perhaps set them up as independent, then discuss why they equal each other and make the new equation eliminated x3." Suggested restructure: introduce the bracket stiffness matrix with independent $x$ and $z$ entries, then add:

```latex
Although $x$ and $z$ are introduced as independent parameters, the bracket geometry makes them equal in this configuration: the insertion bracket is symmetric about the plane containing the force path, so loading it produces equal axial stiffness in the two transverse directions, and the identification uses the y-symmetric form $K = [X_1, X_2, X_1]$. Setting $x = z$ eliminates the third parameter and leaves two identified stiffnesses per bracket.
```

**OPEN for Ben:** confirm the phrasing matches his mechanics intent (the repo's frame-convention proof shows the identification is invariant to this choice). **minus score** (unreadable derivation flagged).

### M3 — L192: where the body-frame force vector comes from — PROPOSED (previously missed)
ajh26 (two comments, same line): "I'm not sure what this means" then "should be made a little more clear where this value comes from." Add after the sentence:

```latex
This force vector is the actuator force of Section~\ref{sec:force_model}, evaluated at the current pressure and contraction, expressed in the bracket frame by the rotation of the preceding paragraph.
```

**OPEN for Ben:** confirm the exact source chain matches the derivation order in that subsection. **minus score** (two complaints on one sentence).

### M4 — L254: "Have Stu double check the equations" — OPEN (process, not prose)
Offer stands: ZCode can cross-check every equation in sec 3.6 against the live MATLAB evaluators (`minimizeFlxPin.m`, `minimizeExtX3.m`, `MonoPamDataExplicit_balanceX3.m`) so Stu verifies against a checked manuscript. Say the word; no score.

### M5 — L274: Newton-Raphson stop conditions — PROPOSED (previously missed)
ajh26: "Stop conditions?" Add after "solve the nonlinear force-balance equation":

```latex
Iterations stopped when the change in muscle length between iterations fell below [FLAG-BEN: tolerance] or after [FLAG-BEN: max iterations] iterations, whichever came first; the values are those used in \texttt{minimizeFlxPin.m}.
```

**FLAG-BEN:** ZCode can pull the actual tolerance/iteration cap from the MATLAB source on request rather than guess. **minus score**.

### M6 — L280 (three comments): section start + justify optimization + what the attachment update is for — PROPOSED (previously missed; covers the L317 "whole section confusing", L323 "described above?", L330 "prediction fails?", L334 "explained more clearly and clearly justified" cluster)
ajh26 wants: a new section for parameter optimization starting at the "In..." paragraph; a good justification for optimization over analytical methods; clarity on what the attachment-update paragraph is doing. Suggested skeleton:

```latex
\section{Identification of the Correction Terms}\label{sec:identification}
```
opening with the MATLAB/optimization-first framing ajh26 praised at L321, followed by:

```latex
The correction terms cannot be derived analytically: the brackets are 3D-printed arches whose effective stiffness is not captured by beam theory, the series tendon compliance is unknown, and the wrap loss depends on the actuator's nonlinear force--length curve at the wrapped configuration. Each term is therefore identified numerically from the measured torque data, with multiobjective optimization exposing the trade-off between fitting the pinned-knee tests and predicting the biomimetic configurations.
```
plus one sentence at the attachment-update paragraph stating its role: the bracket-stiffness terms identified here are the same terms the placement optimizer of Section~\ref{sec:optimization} carries as model corrections, so attachment geometry and correction terms are identified against the same data. Fix "described above" and "prediction fails" wording while restructuring. Full skeleton in `ZCode_drafts/chapters/20-methods.tex` (section header + opening only; body move is Ben's call on placement). **minus score** (ajh26 called the section confusing and unjustified; the L321 framing he praised is the opener).

### M7 — L292: "arch shape extends 75 mm in Z" confusing + new paragraph — PROPOSED (previously missed)
ajh26: "Do you just mean because it wraps around and makes contact with the knee?" + "new paragraph?" Suggested:

```latex
Because the bracket wraps around the knee and makes contact with it, its compliance cannot be modeled as a single beam; it is treated as an effective stiffness identified from the torque data.
```
Start a new paragraph at that sentence. **Parked for the Discussion chapter** (ajh26's own rumination): "since muscles have similar force-length curves, do we not see large squishy creatures because of the loss in force? Not a lot of protruding attachment points?" — worth a paragraph there, not here. **minus score**.

### M8 — L313 tracked " only"→"," — DONE (accepted).

### M9 — L336: prior walker section belongs in Background — PROPOSED (previously missed)
ajh26: "Should probably be in the background." Move `\section{Prior Neuromechanical Simulation: The AnimatLab Bipedal Walker}` (sec 3.7) to the end of Chapter 2, after sec 2.7 where it is already previewed ("the model is described in Section~\ref{sec:prior_walker}"). Its three accepted tracked deletions (course-report sentence; push-off clause; vestibular-layer clause) move with it. **minus score**.

### M10 — L342: em-dash exemplar + "stand alone" rephrase — PROPOSED
Ben's named em-dash example. Replace the final sentence with:

```latex
That hand-tuning consumed the bulk of the project's effort and left no record of how sensitive the behavior was to any individual parameter. It also failed to generalize when new feedback pathways were added. These limitations directly motivated the systematic, verification-first tuning approach of Chapter~\ref{ch:futurework}.
```

Variant if ajh26 wants it to stand alone without the forward reference: end at "...added. The verification-first methodology of this dissertation is the response to exactly these limitations." **minus score** (AI-tell construction).

### M11 — L346: drop "Preliminary", tool/pipeline tone — PROPOSED (previously missed)
ajh26: "Rewrite this section to have more of a tone that you developed a tool or a pipeline for use. Get rid of the word preliminary." Retitle to `Neuromechanical Simulation Workflows` and reword the opener:

```latex
Simulation workflows were assembled to integrate the BPA mechanics with the synthetic nervous system (SNS). The figures in this section document the biomechanical representations, controller organization, and software implementation constructed for that purpose; quantitative evaluation follows in Section~\ref{sec:preliminary_sim_results}.
```

**minus score**.

---

## 30-results.tex (new sections supporting the moves above)

### R1 — Completed Foundation Studies — PROPOSED (destination for I2)
Insert after the publication note:

```latex
\section{Completed Foundation Studies}\label{sec:foundation}

Four studies completed during this dissertation establish the foundation that the remaining results build on, and are summarized here so that the experimental chapters read as one arc. First, the placement optimization of Morrow, Bolen, and Hunt found artificial muscle attachment locations that reproduce human torque profiles and ranges of motion about the trunk and leg joints \citep{morrow_optimization_2020}; it converted the earlier finding that transplanted human attachment locations do not transfer to a BPA-actuated robot \citep{bolen_bipedal_2019, bolen_determination_2019} into a usable design tool, and the route redesigns of Section~\ref{sec:ongoing} re-run that tool under the corrected torque model. Second, the Sensory Afferent Database organized the literature on sensory afferent properties into a structured, model-oriented resource \citep{bolen_sensory_2023}, and a follow-up study established how far AI-assisted literature triage can accelerate its expansion \citep{bolen_using_2024}; the afferent data it curates underpin the sensory channels of the neural controller in Chapter~\ref{ch:futurework}. Third, an inverted-pendulum balance platform with BPA-driven plantarflexion and dorsiflexion demonstrated pseudo-vestibular feedback using inertial measurement units \citep{bolen_balance_2024, bolen_selfbalancing_2024, bolen_selfbalancing_srs_2024}, a first hardware step toward the vestibular layer of the humanoid. Fourth, companion work by Elzein, Lutz, Bolen, and Hunt produced a state-space, optimization-based model of pressure and force dynamics in pulse-modulated BPAs \citep{elzein_pressure_2025}, which supplies the dynamic actuator behavior that the static force characterization reported below leaves out.
```

### R2 — Human Knee Torque Validation (write-as-done with flags) — PROPOSED
At the end of sec:ongoing, per Ben's instruction that this experiment gates submission:

```latex
\subsection{Human Knee Torque Validation}\label{sec:human_torque_validation}

% FLAG-BEN: This subsection is written as complete pending the redesigned-part tests.
% Insert the measured values, the angle coverage, and the validation figure once the
% human knee torque data are collected and analyzed, then delete this flag block.

The redesigned biomimetic knee was then tested isometrically against the human torque targets using the procedures of Sections~\ref{sec:torque_equipment} and~\ref{sec:hybrid_method}. The measured extensor torque met or exceeded the vastus medialis target over \textbf{[FLAG-BEN: \% of RoM]} of the range of motion, with peak torque of \textbf{[FLAG-BEN: value]}~N{\,}m at \textbf{[FLAG-BEN: angle]} of knee flexion. The measured flexor torque met the biceps femoris short head target over \textbf{[FLAG-BEN: \% of RoM]} of the range of motion. Across all tested angles, the corrected model predicted the measured torque with RMSE of \textbf{[FLAG-BEN: value]}~N{\,}m for the flexor and \textbf{[FLAG-BEN: value]}~N{\,}m for the extensor, confirming that the correction terms identified on the pinned knee transfer to the redesigned biomimetic configuration.
```

Abstract and Dissertation Organization reference the same FLAG-BEN values. Also flagged, not changed: sec:ongoing still says the flexor re-run "was still in progress at the time of writing" — flip once Opt_run lands. **add score** (placeholder system built before it was needed).

---

## main.tex (21 July comments)

### X1 — 2-3 pages of scientific motivation — OPEN
Pre-dates the current structure; the Motivation section is one paragraph. If ajh26 still wants the long-form motivation, it is a writing pass of its own (Ben's call whether it merges with Research Objectives per I5).

### X2 — "You also need a Background chapter" — DONE (exists as Chapter 2; resolve in panel).

### X3 — big-paper format (one intro/methods/results/discussion chapter each) — DONE (current structure matches); resolve in panel.

### X4 — 40-discussion "To be removed" — OPEN
Tension: X3 asks for a merged discussion chapter; X4 says remove this one. If X3 is satisfied, X4 likely means the old per-paper discussion is superseded — but confirm with ajh26 before deleting a chapter.

---

## Machine-teaching ledger

- Praise: 2 finds → **add score** x2 (B7 "Yes, like this!"; M-L321 "much better way to introduce the optimization"). Lesson: forward-pointing sentences telling the reader where a model gets used are the house style.
- Tracked deletions/replacements accepted: 17 → **minus score** x17. Pattern in what he cut: hedges ("maximum"), qualifiers ("only", " at all"), future-tense promises ("planned", "will"), anything implying unfinished work.
- Critical comments: ~31 (including the previously missed M-cluster) → **minus score** each. Patterns: paragraphs opening mid-thought; subjects appearing from nowhere ("Compliant", "Recent"); unexplained symbols and tools; unattributed agents; sections without stated purpose or justification.
- Unchanged regions (results, discussion, appendices): **slight add score**, discounted for read depth — assume he has read only through ~methods L346.
- **Ben's global grade, 2026-09-12: B- on "combining my papers and research into something coherent"** — draft, not finalized. Signal no item-level comment produces: skeleton and ordering work; the seams (voice consistency across merged sources, uneven section depth, placeholder-heavy results) cost the grade. Next-pass targets: one voice across merged sections, deepen thin sections, resolve FLAG-BENs with data.
- Ben's question on whether the scoring conception helps, answered: yes — the labeled examples (which text drew which judgment) are what generalizes; the numeric totals are bookkeeping. His praise/delete/uncommented-discount scheme is sound reward shaping.

## Learning log

1. Em-dashes are out. Rewrite into one or two sentences; never a comma swap.
2. Last sentence of a paragraph/section: bridge to what follows AND state the dissertation's stake (B7 is the template).
3. Completed-work framing everywhere: "these tools let researchers do X", never "planned/future work".
4. The SOTA section sitting after the sections it overviewed was the structural embarrassment; content now distributed to where it is used.
5. Human knee torque: written as complete with FLAG-BEN placeholders; gates submission and the defense date.
6. Attribution precision: Quinn lab vs AARL vs wider community; ajh26's rat program (himself, Deng, Young, Jackson) belongs in sec 2.7.
7. summary.md structure rule (Ben, 2026-09-16): problem and solution live in the same entry; never make him scroll to match them.
8. No Overleaf edits without an explicit instruction naming the action. Suggested text goes here; Ben applies.
9. Independent AI-writing-check agent: still to be set up on Ben's go (rubric: em-dash density, hedge stacking, uniform sentence rhythm, AI lexicon, rule-of-three overuse).

## Open questions for Ben

1. B10/B13: confirm the li_neural_2016 → li_bipedal_2017 swap.
2. B14: accept or reject the pending "spinal " deletion (recommend reject).
3. I4: promote torque results to their own chapter, or keep one Results chapter?
4. I5: summary sentences (Ben's spec) vs full paragraphs (ajh26's ask), and whether to merge Dissertation Organization into Research Objectives.
5. M2/M3/M5: confirm mechanics wording and the NR stop-condition values (ZCode can pull them from the MATLAB source).
6. M9: move the prior-walker section to Background, or wait for the structural pass?
7. X4: confirm with ajh26 before removing 40-discussion.
8. Drafts now live at `ZCode_drafts/chapters/` in this folder (regenerated 2026-09-16 after the Temp wipe).
