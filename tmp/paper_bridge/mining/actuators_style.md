# Style card — the published Actuators force paper (bolen_2026)

Calibration reference for section writers of the Frontiers knee-torque paper.
Every number below was measured from the final proofread manuscript source
(not eyeballed), as described in "Provenance" at the end.

## 1. At a glance

| Metric | Value |
|---|---|
| Body prose (excl. floats, equations, captions) | ~4,430 words |
| Abstract | 153 words, 5 sentences, 1 paragraph, zero citations |
| Sections | Introduction, Background, Materials and Methods, Results, Discussion, Conclusions |
| Section share of prose | Discussion 52%, Intro 13%, Results 14%, Background 8%, Methods 7%, Conclusions 6% |
| Figures | 5 (1 hardware photo in Methods; 4 data figures in Results) |
| Tables | 3 (all coefficient/goodness-of-fit tables; no descriptive-equipment tables) |
| Labeled equations | 9 (4 prior-work/definition eqs in Background; 4 new fits in Results; 1 composition eq in Discussion) |
| References | 30 unique keys, 55 in-text calls (~1.2 calls per 100 body words) |
| Highlights | 3 bullets, one clause each |
| Compiled length | 17 pp in MDPI accept-mode proof; ~10 pp typeset |

The headline structural fact: **the Discussion is half the paper** (~2,300 words
— more than Methods + Results + Conclusions combined). Results subsections are
short (100–200 words each) because each one parks its evidence in a figure or
table and moves on; interpretation, mechanism, limitations, comparisons, and the
worked example all live in the Discussion. The BRIEF's "LONG Discussion" is not
an embellishment of this paper's style — it is the paper's style.

## 2. Section skeleton with measured lengths

1. **Introduction** — 581 words, 5 paragraphs, 20 citation calls. Funnel:
   (P1) biomimetic robots as experimental platforms for motor-control questions;
   (P2) tour of actuation alternatives (motors, hydraulics, DEAs) each dismissed
   in one clause with citations, landing on BPAs; (P3) three numbered
   BPA-vs-muscle differences ("First… Second… Third…"); (P4) model-by-model
   review with `\citet` narrative citations, ending on the gap ("none of the
   existing models effectively take into account differences in maximum force
   due to initial actuator length"); (P5) "The present work develops…" — one
   paragraph stating what was collected, what was derived, and what it supports.
2. **Background** — 353 words, ~4 paragraphs, 4 equations, 4 citation calls.
   Quotes prior-work equations verbatim (with the correction note: "the
   coefficients have been updated with the correct values as the values
   reported in [Hunt] contained typographical errors"), then chains the strain
   definitions ε, ε₆₂₀, ε* used later. This is where all borrowed math lives.
3. **Materials and Methods** — 291 words, 2 subsections ("Overview"; "BPA Force
   Characterization Experiment"), 1 figure (labeled test-stand photo), 2
   citation calls, **zero equations**. Hardware named with brand and series
   (80/20 1515 extruded aluminum, S-shaped load cell, Freescale MPX5700),
   measurement conventions defined against datasheet definitions with
   citations, test matrix given as siunitx ranges.
4. **Results** — 618 words, 4 subsections, each anchored to one figure/table:
   "Maximum Force at 620 kPa"; "Force as a Function of Pressure and Resting
   Length"; "Maximum Contraction at 620 kPa"; "Force as a Function of Pressure
   and Contraction". Contains the four new fit equations.
5. **Discussion** — 2,297 words, 5 subsections: "Effects of Resting Length";
   "Nondimensionalization"; "Comparison with Other Published BPA Models";
   "Implications for Real-Time Control and Simulation"; "Biomimetic Artificial
   Muscle Analogy". 26 citation calls — heavier cited than the Introduction.
   Includes the composition equation F = F*·F₆₂₀ (eq:force) and a fully worked
   numerical example ("given l_rest = 54 mm and P = 300 kPa, the measured force
   was … This is an error of 8%…").
6. **Conclusions** — 287 words, **one paragraph**, 0 citations. "The analysis
   in this study has created novel equations for…" → recap of each deliverable
   with its equation/table reference → the design-workflow payoff → expected
   downstream use. No new information, no hedging.
7. Back matter: CRediT author-contributions block, funding, data availability
   ("available from the authors upon reasonable request"), acknowledgments,
   conflicts, one-row abbreviations table, references.

## 3. Abstract shape (153 words, 5 sentences)

The rhetorical template, sentence by sentence:

1. *Motivation + asset*: "Artificial muscles such as braided pneumatic
   actuators (BPAs) offer many advantages for robotic systems, including…"
2. *Gap*: "However, their use in robotic systems is still extremely limited,
   in part due to their poor force, length, and pressure characterization."
3. *What we did*: "In this work, a test setup is created to compare force
   produced by Festo fluidic BPAs with leading models."
4. *Enumerated contributions*: "Our analysis of the data has resulted in (1)
   the development of new equations to calculate force … and (2) a novel
   equation for the maximum force … as a function of their resting length."
5. *Implication*: "This will lead to faster design processes and the
   development of new systems such as biomimetic robots that are able to more
   accurately reproduce…"

No citations, no numbers, no abbreviations beyond BPA/PAM, no methods detail.
The "(1) … and (2) …" enumeration inside a single sentence is the paper's
signature move for contributions — it recurs verbatim in spirit in the
Conclusions.

## 4. How equations are introduced

- **Pattern A — fit equations (Results):** a sentence names the observed shape
  and the fitting method, ends with a colon, then the equation:
  "Using the Nonlinear Least Squares method and a Least Absolute Residual
  robustness, we fit an arctan curve to the data to get the maximum force at
  620 kPa given a resting length as:" → equation. Immediately after:
  coefficients with units in running text ("where a₁ = 0.4895 N/kPa and
  a₂ = 0.03068 kPa⁻¹m⁻¹ for the 10 mm actuator…"), then "Goodness-of-fit
  measures are given in Table~\ref{tab:Fmax}."
- **Pattern B — prior-work equations (Background):** "\citet{hunt_modeling_2017}
  looked at six resting lengths… the scalar pressure P … can be determined by
  solving the equation:" → equation → "where S is the artificial muscle
  hysteresis factor such that S = 1 indicates…" Every symbol gets a "where…"
  definition sentence; units ride inside the equation via `\qty` on
  coefficients (e.g. `\qty{254}{\kPa}`, `\qty{1.23}{\kPa\per\newton}`).
- **Pattern C — definition equations (ε, ε₆₂₀, ε*):** minimal one-line
  intro ("The amount of contraction, ε, is calculated as") + equation; later
  equations reference earlier ones ("using the results from
  Equations~(\ref{eq:strain}) and (\ref{eq:max_strain})").
- Equations are referenced as "Equation~(\ref{eq:x})" (always parenthesized)
  or "Equations~(\ref{eq:a}) and (\ref{eq:b})"; 22 such references in the
  body. The compounding equation appears in the Discussion, not Results —
  interpretation-level math is allowed to live in the Discussion.

## 5. Results voice

- **Data as subject, present tense**: "Data of the maximum force … show a
  dependency on resting length", "The data show a force response resembling an
  arctan curve", "it is clear that the Festo tool over-predicts the expected
  force".
- **"We" for actions, past tense for procedures**: "we fit an arctan curve",
  "We collected extensive force data", "The upper end was attached to an
  S-shaped load cell."
- **Novelty flagged plainly, without adjectives**: "This is a previously
  unreported characteristic of these artificial muscles."
- **Negative results reported without apology**: the ε₆₂₀-vs-l_rest subsection
  is one paragraph concluding "there is no direct, predictable relationship",
  with the practical instruction "this value should be recorded in each muscle
  used on a robot".
- Statistics quoted inline to 3–4 significant digits ("adjusted R² = 0.9998, a
  RMSE = 0.004537, and a maximum absolute residual of 10.6%"), with the full
  set in the table.
- Almost no citations in Results (exactly 1 call) — Results reports own data;
  comparison and credit happen in the Discussion.

## 6. Figures and tables: caption and reference conventions

**Figures (5):**
- Referenced as capitalized `Figure~\ref{fig:x}`, usually parenthetical:
  "(Figure~\ref{fig:Fmax_fit})"; 10 in-text figure references.
- Multi-panel figures use `\subfloat` with panel letters referenced as
  "(Figure~\ref{fig:F-P_contraction}\textbf{b})" — bold lowercase letter
  appended outside the ref.
- Captions are long and self-contained (60–110 words): what is plotted, then
  panel-by-panel "(\textbf{a}) … (\textbf{b}) …" descriptions, then the visual
  encoding legend IN WORDS ("Dashed line is the F₆₂₀ data from Festo. Solid
  line is the fit from Equation (\ref{eq:Fmax10})."), plus data provenance
  ("In total, 537 data points were collected … 80% of these data were used for
  training, and 20% were used for validation."). A caption that explains the
  error bars is present where relevant ("The error bars show the effect of
  being ±1 mm in measuring l_rest and l₆₂₀.").
- The Methods figure is a labeled photo of the test stand.

**Tables (3):**
- All three are coefficient tables with 95% CIs and goodness-of-fit columns
  (Adj. R², RMSE, Max. Error; one adds FVU), one row-block per equation or per
  diameter, booktabs rules (`\toprule/\midrule/\bottomrule`), `tabularx` to
  `\fulllength` in an `adjustwidth` for wide ones.
- Captions sit ABOVE the table (the `\caption` precedes the tabular), are full
  sentences, define every abbreviation used ("FVU denotes the fraction of
  variance unexplained"), and state units ("RMSE and maximum error are
  reported in newtons"). 6 in-text "Table~\ref" references, always
  capitalized, never abbreviated.
- Coefficients carry units inside the table via `\qty` (e.g.
  `\qty{0.4895}{\N\per\kPa}`), CIs as "(0.4822, 0.4968)".

## 7. Citation density and style

- natbib author-year: `\citet` (16 uses) for narrative model-by-model review;
  `\citep` (35 uses) for clustered support — 2–3 keys per parenthetical is the
  norm at claim-heavy moments.
- Per-section calls: Introduction 20, Background 4, Methods 2, Results 1,
  Discussion 26, Conclusions 0. Target ≈ 1.2 calls per 100 words overall, with
  the Introduction and Discussion carrying nearly all of it.
- Datasheets are first-class citations: three Festo datasheets cited with
  page-specific pointers — "\citet{festo_2026} (p.~12) states that force will
  be reduced by approximately 10%…", "\citep{festo_2018} (pp.~10–11)". This is
  the house precedent for citing manufacturer data with page numbers.
- The Discussion credits and compares rather than attacking: "The
  \citet{martens_modeling_2017} model is derived from a geometric and
  material-based framework and provides a more physically grounded
  representation… However, its accuracy depends strongly on…"

## 8. Register notes (voice, tense, person)

- First-person plural throughout for the team's actions; passive/past for
  equipment procedure. Present tense for what data and models do.
- Hedged speculation is explicit and personal: "We still suspect that it is a
  function of l_rest", "It could be a reversed sign coefficient in their
  internal curve-fit model", "It remains to be seen how ε₆₂₀ can be known a
  priori."
- Direct practical address to the reader: "Researchers using BPA resting
  lengths under 300 mm should take note of these results."
- Paragraph rhythm is irregular (3–10 sentences), topic-sentence-led; sections
  end on implication, not summary.
- The word "plant" never appears (grep-verified); test hardware is "test jig",
  "test stand", "test setup".
- Occasional genuine human roughness survived copyediting ("alumnium";
  "the describes the passive force–length curve has been described") — proof
  of a non-template register, not a feature to imitate.

## 9. Do-not-copy quirks (deviations the new contract overrides)

- **Diameter symbol**: the Actuators paper uses `$\phi$\qty{10}{\mm}`. The
  Frontiers contract mandates the engineering symbol `$\diameter$` (mathabx) —
  a deliberate deviation from this exemplar.
- The paper has a standalone **Background** section for prior-work equations;
  the Frontiers structure in the BRIEF is Introduction → Methods and Materials
  → Results → long Discussion → Conclusion. Borrowed math (rigid-body torque,
  BPA force model) should fold into Introduction/Methods instead unless Ben
  rules otherwise.
- One `\label{eq:BPAFit}` sits AFTER `\end{equation}` instead of inside the
  environment (works, but fragile) — put labels inside.
- Inconsistent unit markup in places ("10 mm and 20 mm muscles" bare, next to
  `\qty` everywhere else); the new contract bans bare units outright.
- Author-facing artifacts (MDPI `\hl` highlights, proof comments) obviously do
  not transfer.

## 10. Calibration targets for the Frontiers knee-torque paper

Scaled to the BRIEF's scope (Methods-heavy experimental paper with additional
system-identification and route-redesign results folded into the Discussion):

- Body prose budget ~5,500–7,000 words (the added Xi/two-bracket and
  mesh-optimization content justifies ~30–50% more than the 4,430-word
  exemplar, still front-loading the Discussion).
- Ratios to aim for: Discussion ≥ 2× Results; Results ≥ 1.5× Methods words
  (the torque rig and protocol are more involved than the force jig, so
  Methods grows relative to the exemplar's 7% share — toward 12–18%).
- 6–8 display items: 1 test-stand figure in Methods, the rest in Results/
  Discussion; every results figure paired with a coefficient/GoF table where a
  fit exists. Caption style per §6 (encoding-in-words, panel letters, counts).
- Equations: pattern per §4 — method sentence ending in a colon, `\qty`-units
  on coefficients, "where…" definitions, GoF pointer to a table. Prior-work
  and definition equations before the new fit equations in reading order.
- Citations: ~1–1.5 calls per 100 words, concentrated in Introduction and
  Discussion, ~1 in Results, 0 in Conclusions; datasheet cites with page
  pointers where manufacturer data is invoked.
- Abstract: ~150 words, 5 sentences, the exact §3 template with a "(1) … (2)
  … (3) …" enumerated-contribution sentence covering test-rig results,
  correction factors, and the redesign outcome.

## Provenance — what was read and measured

- Source of record: `Documentation/Reports and Papers/Isometric_Force_Test/
  Final/Final_proofread/manuscript.v5/actuators-4207821-edited.tex`
  (final proofread v5, compiled PDF in the same folder; read in full,
  lines 1–931). This is the edited accepted manuscript (MDPI article
  actuators-4207821, \pubyear{2026}, accepted 17 April 2026 = Actuators
  15(5):230 per the BRIEF).
- Counts were computed by script, not estimated:
  `tmp/paper_bridge/mining/_analyze_actuators_style.py` (first pass),
  `_analyze_fix.py` (corrected comment stripping — escaped `\%` had eaten
  trailing line content, hiding one figure label),
  `_analyze_final.py` (final per-section word counts + equation locations),
  all run with `D:/Anaconda/python.exe` from the repo root this session.
  Intermediate first-pass numbers (5 figs vs 4, 53 vs 55 citation calls)
  are superseded by the corrected run reported above.
- Page count: `grep "Output written" actuators-4207821-edited.log` →
  "17 pages, 511547 bytes". The sibling typeset PDF
  `IsometricForceCharacterizationofBPAs.pdf` at the Isometric_Force_Test root
  counts ~10 /Page objects (crude binary count; no pypdf available) — treat
  the 17-page proof compile as the length reference.
- "plant" absence: `grep -in "plant" actuators-4207821-edited.tex` → no
  matches (exit 1).
- Review-round artifacts (`Review_round_1/`, response-to-reviewers) were not
  mined for style — only the final text, per the ask.
