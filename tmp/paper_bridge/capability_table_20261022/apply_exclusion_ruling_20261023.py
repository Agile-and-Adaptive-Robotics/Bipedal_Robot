import csv
from datetime import date

CAMPAIGN = "tmp/paper_bridge/capability_table_20261022"

EXCLUDED = {
    "Rectus Femoris": "knee extensor; legacy hand-code patella moment arms used an incorrect method (should be the dL/ds derivative form per the literature, Ben 2026-10-23) and the optimized side uses the rigid fixed-tibia bypass - neither row trustworthy",
    "Vastus Intermedius": "knee extensor; legacy hand-code patella moment arms used an incorrect method (should be the dL/ds derivative form per the literature, Ben 2026-10-23) and the optimized side uses the rigid fixed-tibia bypass - neither row trustworthy",
    "Erector Spinae": "optimized-row solver record internally inconsistent (torque present with zero positive-strain cells; rest/strain/torque records contradict), unresolved - excluded per Ben 2026-10-23",
    "Semimembranosus": "optimized-row solver record internally inconsistent (torque present with zero positive-strain cells), unresolved - excluded per Ben 2026-10-23",
    "Medial Gastroocnemius": "optimized-row solver record internally inconsistent (torque present with zero positive-strain cells; recorded rest implausibly short for the route), unresolved - excluded per Ben 2026-10-23",
    "Lateral Gastrocenemius": "optimized-row solver record internally inconsistent (torque present with zero positive-strain cells), unresolved - excluded per Ben 2026-10-23",
    "Gluteus Maximus": "optimized-row solver record internally inconsistent (torque present with zero positive-strain cells), unresolved - excluded per Ben 2026-10-23",
    "Adductor Magnus": "optimized-row solver record internally inconsistent (torque present with zero positive-strain cells), unresolved - excluded per Ben 2026-10-23",
}

rows = list(csv.DictReader(open(f"{CAMPAIGN}/per_muscle_capability_20261023.csv", newline="", encoding="utf-8")))
kept = [r for r in rows if r["muscle"] not in EXCLUDED]
excl = [r for r in rows if r["muscle"] in EXCLUDED]
assert len(kept) == 19 and len(excl) == 8, (len(kept), len(excl))

def avg(field):
    return sum(float(r[field]) for r in kept) / len(kept)

hand_meet, hand_act = avg("hand_meet_pct"), avg("hand_active_pct")
opt_meet, opt_act = avg("opt_meet_pct"), avg("opt_active_pct")
try:
    hand_meet_u, hand_act_u = avg("hand_meet_pct_unflagged"), avg("hand_active_pct_unflagged")
    opt_meet_u, opt_act_u = avg("opt_meet_pct_unflagged"), avg("opt_active_pct_unflagged")
except KeyError:
    hand_meet_u = hand_act_u = opt_meet_u = opt_act_u = None

better = [r["muscle"] for r in kept if float(r["opt_meet_pct"]) > float(r["hand_meet_pct"])]
worse = [r["muscle"] for r in kept if float(r["opt_meet_pct"]) < float(r["hand_meet_pct"])]
ties = [r["muscle"] for r in kept if float(r["opt_meet_pct"]) == float(r["hand_meet_pct"])]
print(f"kept 19 | hand {hand_meet:.2f}/{hand_act:.2f} | opt {opt_meet:.2f}/{opt_act:.2f}")
print(f"better {len(better)}, tie {len(ties)}, worse {len(worse)}")

# --- per-muscle CSV: add exclusion columns, keep all 27 rows
fields = list(rows[0].keys()) + ["excluded", "exclusion_reason"]
with open(f"{CAMPAIGN}/per_muscle_capability_20261023.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for r in rows:
        r["excluded"] = "yes" if r["muscle"] in EXCLUDED else "no"
        r["exclusion_reason"] = EXCLUDED.get(r["muscle"], "")
        w.writerow(r)

# --- table markdown
md = f"""# 27-actuator capability table (2026-10-23, v2 coupled-DOF surfaces, 19 graded)

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
| Hand-tuned (19 of 27 graded) | {hand_meet:.2f}% | {hand_act:.2f}% |
| Optimized (19 of 27 graded) | {opt_meet:.2f}% | {opt_act:.2f}% |

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
human torque over a strictly larger share of the RoM on {len(better)}
actuators ({len(ties)} ties, {len(worse)} worse). Full caveats:
`CAVEATS_REPORT_20261023.md`.
"""
open(f"{CAMPAIGN}/capability_table_20261023.md", "w", encoding="utf-8").write(md)

# --- table LaTeX
u_note = ""
if opt_meet_u is not None:
    u_note = (
        f"% without-flagged-cell variants (caveats report): hand {hand_meet_u:.2f}/{hand_act_u:.2f}, opt {opt_meet_u:.2f}/{opt_act_u:.2f}\n"
    )
tex = (
    "% 27-actuator capability table (2026-10-23, v2 coupled-DOF surfaces, 19 graded of 27)\n"
    + u_note
    + "% Suggested caption: Coupled-DOF capability of the hand-tuned and BiPulley-optimized\n"
    + "% BPA routings against single-muscle human torque surfaces (Gait2392, OpenSim 4.4.1),\n"
    + "% mean per-actuator area share over the 19 actuators graded on both rows.\n"
    + "\\begin{tabular}{lcc}\n\\toprule\n"
    + " & \\% of RoM torque $\\geq$ human & \\% of RoM with active force \\\\\n\\midrule\n"
    + f"Hand-tuned (19 of 27) & {hand_meet:.1f}\\% & {hand_act:.1f}\\% \\\\\n"
    + f"Optimized (19 of 27) & {opt_meet:.1f}\\% & {opt_act:.1f}\\% \\\\\n"
    + "\\bottomrule\n\\end{tabular}\n"
)
open(f"{CAMPAIGN}/capability_table_20261023.tex", "w", encoding="utf-8").write(tex)

# --- xlsx
from openpyxl import Workbook, load_workbook
wb = Workbook()
ws = wb.active
ws.title = "capability 2x2"
ws.append(["Row", "% of RoM torque meets or exceeds human", "% of RoM with active force"])
ws.append([f"Hand-tuned (19 of 27 graded)", round(hand_meet, 2), round(hand_act, 2)])
ws.append([f"Optimized (19 of 27 graded)", round(opt_meet, 2), round(opt_act, 2)])
ws2 = wb.create_sheet("per muscle")
with open(f"{CAMPAIGN}/per_muscle_capability_20261023.csv", newline="", encoding="utf-8") as f:
    for row in csv.reader(f):
        ws2.append(row)
wb.save(f"{CAMPAIGN}/capability_table_20261023.xlsx")

# --- caveats append
with open(f"{CAMPAIGN}/CAVEATS_REPORT_20261023.md", "a", encoding="utf-8") as f:
    f.write(f"""

## Exclusion ruling (Ben, 2026-10-23) — the table averages 19 of 27

Eight actuators are EXCLUDED from both rows of the final table (values kept
in the per-muscle CSV with reasons; nothing deleted):

1. **Knee extensors (Rectus Femoris, Vastus Intermedius):** Ben ruled that
   the legacy hand-code patella location moment arms were calculated
   incorrectly — the literature's ds/dl (dL/d-theta) derivative form is the
   correct moment-arm calculation — and the optimized row uses the rigid
   fixed-tibia bypass, so neither row is trustworthy for the quads.
2. **Six solver-record-inconsistent rows (Erector Spinae, Semimembranosus,
   Medial Gastrocnemius, Lateral Gastrocnemius, Gluteus Maximus, Adductor
   Magnus):** optimized winners record torque in cells with zero
   positive-strain cells anywhere, with resting lengths that contradict the
   recorded strain (e.g., medial gastrocnemius rest 0.139 m on a ~0.4 m+
   route). Root cause unresolved by ruling of no-time; the re-grade pass was
   stopped and the rows excluded rather than fixed (Ben: "We don't have time
   to fix everything").

Resulting averages over the 19 graded actuators: hand {hand_meet:.2f}% meet /
{hand_act:.2f}% active; optimized {opt_meet:.2f}% / {opt_act:.2f}%. Strictly
better {len(better)} / ties {len(ties)} / worse {len(worse)} on meet share.
""")
print("done: table md/tex/xlsx + per-muscle CSV + caveats updated")
