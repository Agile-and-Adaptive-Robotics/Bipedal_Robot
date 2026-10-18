# SNS + DISCUSSION chapter redesign — Ben's directive, 2026-10-17 (captured near-verbatim)

Source: Ben's message of 2026-10-17 (while machining final torque parts). This file is the
execution spec for the next chapter passes. Submission deadline 2026-10-21; defense 11-06.

## A. Global rule

The WHOLE torque chapter, SNS chapter, and the Discussion chapter (before conclusion) are
to be written to journal standard, or as close as dissertation-appropriate.

## B. SNS chapter (60-sns.tex) redesign

Framing: the chapter is about TOOLING — how existing tools were combined into NEW tools.
We did not invent SNS-Toolbox, SCONE/Hyfydy, OpenSim, or Simscape. What is NEW:
- OpenSim + MyoConverter + SNS-Toolbox together on the FULL lower-body model.
- Bridging to Simulink/Simscape: our own SNS-Simscape derivative library.
- Deng and Li models converted into gait2392_simbody (and robotbody) with different
  neural controllers in different programs; the SolidWorks robot assembly into SNS-Simscape.
"If we just used AnimatLab and the arm-position we would have nothing to write about
(that's the existing AnimatLab paper of record). If we ran the Li and Deng models in
AnimatLab, nothing. But recreating them in SNS-Simscape and/or MuJoCo SNS-Toolbox —
then we have something interesting." Engine differences are expected and part of the story.

METHODS subsections:
1. Existing tools USED: AnimatLab, MuJoCo, MyoConverter, SNS-Toolbox, Simscape Multibody,
   SCONE/Hyfydy, SADb (at this point in the dissertation) — the SADb framed the benchmarks
   and validation measurements; name its segments (Airtable, Zotero, Research Rabbit, our app).
2. Tools CREATED (list above).
3. Training methodology: perturbations and slight scaling differences for robustness
   (per Ben's standing instructions); IK results from MANY subjects of different sizes and
   speeds (the reference library Ben provided — stored at Code\MuJoCo_SNS\spinal\gait_refs\
   *.npz, 43 refs: Falisse predicted + Arnold running + Ong walking; staging history in
   D:\temp\gait_lib_staging; see AGENTS.md spinal notes + digests).

RESULTS:
- Basic reflex demos for SNS-Simscape (it is new).
- Moving models from their original home (AnimatLab Deng and Li) into other programs and
  comparing against their papers: deafferented stepping; afferented feedback WITH STIMULUS
  to RG, PF, and MN layers.
- If original models are lacking: replication fidelity in the other programs via the SAME
  neural-pathway deletions and stimulus protocols. Use image analysis / circuit-vision
  tooling where it helps. If very long: a few examples in chapter, rest to appendices.
- ALSO (Ben previously instructed, status UNCERTAIN): replicate the Shinohara, Di Russo,
  and Shevtsova models in SNS-Toolbox and SNS-Simscape from their papers. Current status
  on this machine: circuit-rule JSONs + first-pass drafts exist in
  Code\MuJoCo_SNS\spinal\replication\ but were NOT built/verified as running models.
  Check recent GitHub commits from the other two machines for progress before claiming.

DISCUSSION subsections:
- The curriculum (5 stages) and the scoring.
- Adding the BPA force actuator model into the programs.
- gait2392_robotbody and the SolidWorks-derived models; scoring results of the robot
  physical models under different neural controllers. Known status: all walk deafferented
  in air, afferented in air, self-balance, and gross walking-like movement under harness;
  harness-free walking NOT confirmed — state honestly after checking other machines' commits.
- The BPA-actuated robot SNS-Toolbox-MuJoCo work belongs in the DISCUSSION (Ben's explicit
  placement; worth further discussion with him).
- LAST subsection: briefly summarize the new tools and segue into future work other
  researchers can do with them (leans heavily into the next chapter).

## C. Discussion chapter (70-overall.tex) redesign

7.1 = overview and summary of everything written so far — NOT SADb-statistics-heavy
(currently half of it is paper counts, ignoring the ~400 deleted). Lead high-level:

"Improved first-pass robot design is now achievable" — name the design factors (Xi0-Xi3,
arguably Xi4) FIRST, then the enablers, in no particular order (not chronological):
* Improved BPA characterization.
* Design factors removing rigid-body assumptions — hold up not only for BPAs on a pinned
  knee but for the more biomimetic knee with reinforced frame.
* Attachment and routing optimization, tiered:
  - the simpler 2020 placement optimization (Connor + Ben, morrow_optimization_2020);
  - Opt_run / Opt_run_Ext route selection (depending on BPA bending): 1+ BPAs meeting or
    exceeding human torque magnitude over the full knee RoM;
  - further generalizations: 1+ BPAs replacing 1+ human muscles over 1+ joints for two
    important DoFs, starting from muscles close to the bone toward the outside, accounting
    for spatial constraints of these and already-mounted BPAs;
  - if only route and tendon length can change: the cam-shape solver (trades biomimicry
    for torque-profile optimization and reduced component forces).
* Sensory-afferent database for locomotion control — including its detailed table, pivot
  table, bubble chart, and the monthly literature skimmer/curator workflow.
* Novel neuromechanical simulation tools.

Subsequent subsections (future-work framed, for other researchers):
1. BPA characterization expansion: drop in the max-contraction-vs-resting-length function
   (20 mm, single product batch — state the batch limitation). Others: test multiple
   batches of 10/20/40 mm for general fit surfaces. Active vs passive vs dynamic force
   (tie to the Hill muscle model appendix; ADD more Thelen and Millard equations there —
   Gaussian shape factors, tau, etc.). Isometric was our work; researchers should extend
   the framework to isokinetic and isotonic/isobaric. Instron stretching with the Onyx or
   aluminum end caps, unpressurized and pressurized. FEA program Ben sketched from the
   Festo patents: symmetry half-model (half circle fixed at outer+inner diameter one end,
   diameter expansion at the other, optional lateral force/fixity), thin-shell/membrane
   layers with inter-layer contact, two separate handed aramid weaves as beam/truss elements.
2. Routing expansion: e.g., design a patella.
3. SADb expansion: AI prompt passes over the database (discover notes/rules, add entries);
   ideas for accuracy and economy; human "real intelligence" evaluation still required.
4. Simulation tools: hypothesis testing, stimulus/deletion/pathway-ADDITION studies;
   ocular, vestibular, cerebellar, basal-ganglia inhibition/excitation; cross-program
   workflow example: full-body OpenSim model (arms+trunk, on MyoSuite) → port to
   SNS-Simscape (GUI if wanted) → import back to SNS-Toolbox MuJoCo for speed.
5. Real-time hardware how-to (step-by-step; overflow to appendices): run SNS-Toolbox on
   the bipedal humanoid in real time on the treadmill setup; incorporate feedback from
   the lab's existing encoders, pressure sensors, liquid-wire strain gauges, pressure
   manifolds, and IMUs; compute Ia, II, and Ib from them; Jetson Orin setup; the code and
   steps; computer-in-the-loop tuning; wiring schematics; neural-pathway schematics.
6. Close: these are options for future researchers → segue to the conclusion.

## D. Related rulings from the same message

- "Faster than real time" claim stands (literature-supported).
- Rybak/Shevtsova diagrams are in (2006 Fig 1, 2016 Fig 5.8); ALSO include the recent
  (~1 yr) Rybak-group paper for the evolution narrative — Shinohara 2025 is already
  included (fig:shinohara_connectome); check for the Shevtsova 2026 eLife (RP107480) and
  add if suitable, with a sentence tracing model evolution 2006 → 2016 → 2023/2025.
- Pulley duplication is WITHIN the torque chapter (Results vs Discussion) — keep the
  Discussion copy, cut the Results duplicate (flow-review item 4).
- Tables 4.2/4.3 stay two-model: the Hunt-model comparison lives in the FlxPin figures
  (slight improvement over original is the point there); the rest of the plots compare
  improved-BPA-rigid vs optimized. Optional: update results plots from Testing_Data per
  Ben's figure-spec guidelines (may go to easteregg2 session).
- CPG_airstepping_figs: Ben removes from Overleaf; locally move to Deprecated.
- The 382 pruned papers stay archived in the app (verify archive/ contents).
