# 27-actuator capability table (2026-10-22)

Each cell is the row average: the share of the 25-point range-of-motion grid
at the actuator's primary coordinate where the condition holds (every other
joint the muscle crosses pinned at 0 deg). Human targets: single same-named
muscle torque curves, Gait2392 via OpenSim 4.4.1
(`human_torque_curves_20261022.csv`).

| Row | % of RoM torque meets or exceeds human | % of RoM with active force |
|---|---|---|
| Hand-tuned (27 of 27 graded) | 72.89% | 98.81% |
| Optimized (16 of 27 graded; 23 built) | 52.00% | 82.25% |

- **Hand-tuned row:** averages all 27 covered actuators (0 excluded).
- **Optimized row:** Ben asked for all 27; **23 of 27 were built**, and the
  average runs over the **16 built muscles that could be graded**. The row
  label states the count actually averaged; the missing are stated, never
  faked:
  - Could not build (4): Erector Spinae, Internal Oblique, External Oblique,
    Vastus Intermedius - the master muscle path crosses no joint (all points
    on one body), so no joint actuator can be formed.
  - Built but grade-unable (7): Rectus Femoris, Tibialis Posterior, Peroneus
    Brevis, Peroneus Longus, Peroneus Tertius, Adductor Magnus, Psoas X -
    built from the seed route, but the regrade spec's hinge chain does not
    model their primary coordinate (knee_angle_r for Rectus Femoris;
    subtalar_angle_r for the three peroneals and Tibialis Posterior;
    hip_adduction_r for Adductor Magnus; hip_flexion_r for Psoas X), so no
    regraded percentage exists. Their opt grade cells are blank in the
    per-muscle CSV.
  - Full reasons, the nine transcription suspects, the six artifact-flagged
    human curves, and the budget disclosure: `CAVEATS_REPORT_20261022.md`.
