# Outline — Background Chapter 2 additions

Status: drafted, polished version in `Notes/overleaf/ch2_additions.tex`.
Not yet in Overleaf. Numbering will auto-adjust (currently after §2.8).

## §2.9 State of the Art in Neuromechanical Simulation

Goal: position the dissertation among the three existing lines of
neuromechanical simulation, then claim the gap this program fills.

1. **Line 1 — reflex-physiology controllers.**
   Geyer & Herr 2010 muscle-reflex walking (positive force/length/velocity
   feedback, no CPG needed). Connects to the cat neuromechanical models already
   cited in §2.6 (Markin 2010, Shevtsova 2016).
2. **Line 2 — learned control (deep RL).**
   Song et al. 2021: RL policies on OpenSim musculoskeletal models, reproduce
   kinematics + activations. Weakness for science use: opaque controller —
   "neurons" are not physiological circuits, so lesion/deletion/stimulation
   experiments lose meaning. De Groote & Falisse 2021 perspective frames plant
   maturity vs interpretability as the open problem.
3. **Line 3 — synthetic nervous systems (where this dissertation sits).**
   Explicit conductance-based neurons mapped to identified circuits; real-time
   (SNS-Toolbox 2023). Unique combination claimed by this program: grounded
   controller + validated pneumatic plant + eventual hardware.
4. **Speed paragraph.** MuJoCo faster than real time; SNS-Toolbox real-time down
   to embedded boards. Overnight = thousands of parameter sweeps / hundreds of
   gait cycles per condition → statistically powered lesion studies, not demos.
   This is the throughput no human or animal protocol can match.

References needed (9 total, all verified — see candidate_references.bib):
geyer_musclereflex_2010, song_deep_2021, degroote_perspective_2021,
todorov_mujoco_2012 (in bib), caggiano_myosuite_2022 (in bib),
nourse_sns_2023 (in bib).

## §2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case

Goal: the ethics argument Ben asked for. Four moves:

1. **Impossible in humans.** Muscle deletion, pathway blocking, neural lesioning,
   tonic spinal stimulation — invasive, irreversible, limited to clinical
   populations, no within-subject control.
2. **Costly in animals → robot as Replacement.** Three Rs (Russell & Burch 1959;
   Tannenbaum & Bennett 2015). Robot = replacement method; duplicated subjects =
   reduction; failures damage hardware not tissue = refinement.
3. **Model standing must be earned (Webb standard).** Webb 2001 (BBS), Webb 2002
   (Nature): robot is a good model when mechanism is tested at the same level of
   organization as the target phenomenon. Map onto our validation chain:
   actuator (bench data) → torque (isometric + human envelopes) → controller
   (tonic-stimulation and deletion phenomena: Selionov 2009, Rybak 2006).
4. **Bidirectional close.** Unethical-impractical-unreachable experiments become
   routine; insights feed biology, next robot iteration, and
   prosthetics/orthotics (same problem: controlling muscle-like actuators for
   human users). Ethical case = scientific case.

## Open decisions for Ben

- Does §2.10 (ethics) belong in Ch. 2 Background, or partially folded into
  Ch. 1 §1.1 Motivation (1 short paragraph there + full version here)?
- Section title wording for §2.10.
- Whether the speed/throughput paragraph in §2.9 should instead move to Ch. 6
  §6.1.2 (it cites forward to the pipeline either way).
