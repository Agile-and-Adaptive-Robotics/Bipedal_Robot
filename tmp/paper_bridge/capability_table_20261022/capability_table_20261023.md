# 27-actuator capability table (2026-10-23, v2 coupled-DOF surfaces, 19 graded)

Each cell is the row average of the per-actuator AREA SHARE on the 5x5
coupled-DOF grading grid (`human_surfaces_5x5_grading_20261023.csv`): the
share of the actuator's graded domain where the condition holds. "Meets or
exceeds human" = the robot torque has the human torque's sign AND
|tau_robot| >= |tau_human| at the cell. "Active force" = positive BPA strain
AND taut, positive torque. Hand-tuned row: old-code routing, old-code BPA
sizes, 620 kPa. Optimized row: fresh BiPulley survey2 solves, same old-code
BPA sizes. Both rows graded against the same single-muscle human surfaces
(stock gait2392_simbody routing, OpenSim 4.4.1). The hip adduction/
abduction envelope (linear coupling, Ben 2026-10-23) bounds the graded hip
domains.

| Row | % of RoM torque meets or exceeds human | % of RoM with active force |
|---|---|---|
| Hand-tuned (19 of 27 graded) | 66.32% | 97.05% |
| Optimized (19 of 27 graded) | 75.30% | 95.23% |

**Excluded from both rows (8 of 27, Ben's ruling 2026-10-23):** the two knee
extensors (Rectus Femoris, Vastus Intermedius) — the legacy hand-code patella
moment arms were computed with an incorrect method (the dL/ds derivative form
per the literature is the correct one) and the optimized side uses the rigid
fixed-tibia bypass, so neither row is trustworthy there — and six actuators
whose optimized-row solver records are internally inconsistent (torque
recorded with zero positive-strain cells; resting lengths that contradict the
recorded strain): Erector Spinae, Semimembranosus, Medial Gastrocnemius,
Lateral Gastrocnemius, Gluteus Maximus, Adductor Magnus. Excluded rows keep
their recorded values in `per_muscle_capability_20261023.csv` with the
exclusion reason; nothing was deleted.

Over the 19 graded actuators, the optimized routing meets or exceeds the
human torque over a strictly larger share of the RoM on 11
actuators (3 ties, 5 worse). Full caveats:
`CAVEATS_REPORT_20261023.md`.
