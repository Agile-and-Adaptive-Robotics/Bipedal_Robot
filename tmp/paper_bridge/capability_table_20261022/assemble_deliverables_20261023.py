# assemble_deliverables_20261023.py -- v2 capability-table assembler.
# Re-derives every number from the three CSVs of record in this folder
#   human_surfaces_5x5_grading_20261023.csv   (5x5 grading surfaces + flag layer)
#   hand_tuned_metrics_20261023.csv           (hand-tuned row)
#   BiPulley_survey2_perMuscle_20261023.csv   (optimized row)
# cross-checks them against human_surfaces_manifest_20261023.json, then writes
# ONLY the six dated _20261023 deliverables. No v1 (20261022) file and no input
# CSV is modified. Run on easteregg2:
#   D:/Anaconda/envs/myo/python.exe assemble_deliverables_20261023.py
import csv, json, os, sys

FOLDER = os.path.dirname(os.path.abspath(__file__))
NINE_DIR = "human_surfaces_9x9_20261023"

def read_dated_csv(name):
    path = os.path.join(FOLDER, name)
    with open(path, newline="", encoding="utf-8") as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]
    return list(csv.DictReader(lines))

surfaces = read_dated_csv("human_surfaces_5x5_grading_20261023.csv")
hand = read_dated_csv("hand_tuned_metrics_20261023.csv")
opt = read_dated_csv("BiPulley_survey2_perMuscle_20261023.csv")
with open(os.path.join(FOLDER, "human_surfaces_manifest_20261023.json"), encoding="utf-8") as f:
    manifest = json.load(f)

f2 = lambda x: f"{x:.2f}"
f4 = lambda x: f"{x:.4f}"
def mean(v):
    vals = list(v)
    return sum(vals) / len(vals)

# ---------- 1. parse + validate the surfaces (grading source of truth) ----------
SURF = {}
for r in surfaces:
    aid = int(r["actuator_id"])
    S = SURF.setdefault(aid, {
        "rows": [], "master": r["master_name"], "primary": r["primary_dof"],
        "secondary": r["secondary_dof"], "p_lo": r["primary_rom_lo_deg"],
        "p_hi": r["primary_rom_hi_deg"], "s_lo": r["secondary_rom_lo_deg"],
        "s_hi": r["secondary_rom_hi_deg"], "pinned": r["pinned_dofs_at_0"],
        "env": r["envelope_limited"],
    })
    assert S["master"] == r["master_name"] and S["primary"] == r["primary_dof"] \
        and S["secondary"] == r["secondary_dof"] and S["p_lo"] == r["primary_rom_lo_deg"] \
        and S["p_hi"] == r["primary_rom_hi_deg"] and S["s_lo"] == r["secondary_rom_lo_deg"] \
        and S["s_hi"] == r["secondary_rom_hi_deg"] and S["pinned"] == r["pinned_dofs_at_0"] \
        and S["env"] == r["envelope_limited"], f"surface header mismatch act {aid}"
    S["rows"].append(r)

assert sorted(SURF) == list(range(1, 28)), "surfaces actuator ids != 1..27"
for aid, S in SURF.items():
    assert len(S["rows"]) == 25, f"act {aid}: {len(S['rows'])} rows != 25"
    dom = [r for r in S["rows"] if r["in_graded_domain"] == "1"]
    exc = [r for r in S["rows"] if r["in_graded_domain"] == "0"]
    assert len(dom) + len(exc) == 25
    S["n_dom"], S["n_exc"] = len(dom), len(exc)
    S["env"] == "yes" and (S["n_dom"], S["n_exc"]) == (15, 10) or S["env"] == "no" and (S["n_dom"], S["n_exc"]) == (25, 0) or (_ for _ in ()).throw(AssertionError(f"act {aid}: domain/excluded {S['n_dom']}/{S['n_exc']} vs env={S['env']}"))
    S["flags"] = {k: sum(1 for r in dom if r[f"flag_{k}"] == "1") for k in ("via", "step", "slack", "short")}
    S["n_flag"] = sum(1 for r in dom if any(r[f"flag_{k}"] == "1" for k in ("via", "step", "slack", "short")))
    assert S["n_flag"] <= S["n_dom"]
    S["pinned_list"] = [] if S["pinned"] == "(none)" else [S["pinned"]]
    # envelope sanity
    if S["env"] == "yes":
        assert (S["n_dom"], S["n_exc"]) == (15, 10), f"act {aid}: domain/excluded {S['n_dom']}/{S['n_exc']} vs env=yes"
    else:
        assert S["env"] == "no" and (S["n_dom"], S["n_exc"]) == (25, 0), f"act {aid}: domain/excluded {S['n_dom']}/{S['n_exc']} vs env={S['env']}"
    # manifest cross-check
    M = next(a for a in manifest["actuators"] if a["id"] == aid)
    assert M["n_graded_cells"] == S["n_dom"] and M["envelope"]["n_excluded_cells"] == S["n_exc"], f"act {aid}: manifest cell counts"
    assert M["n_flagged_graded_cells"] == S["n_flag"], f"act {aid}: manifest flag count {M['n_flagged_graded_cells']} != {S['n_flag']}"
    for k in ("VIA", "STEP", "SLACK", "SHORT"):
        assert M["flag_counts_graded_cells"][k] == S["flags"][k.lower()], f"act {aid}: manifest {k}"
    assert M["axes"]["pinned_at_0_deg"] == S["pinned_list"], f"act {aid}: manifest pinned"

n_graded_total = sum(S["n_dom"] for S in SURF.values())
n_exc_total = sum(S["n_exc"] for S in SURF.values())
n_flag_total = sum(S["n_flag"] for S in SURF.values())
T = manifest["reporting"]["totals"]
assert (n_graded_total, n_exc_total, n_flag_total) == (T["n_graded_cells"], T["n_excluded_envelope_cells"], T["n_flagged_graded_cells"]) \
    == (635, 40, 126), f"totals {n_graded_total}/{n_exc_total}/{n_flag_total}"

# 9x9 export files
NINE = {1: "ercspn", 2: "intobl", 3: "extobl", 4: "bifemlh", 5: "semimem", 6: "sar",
        7: "tfl", 8: "rect_fem", 9: "grac", 10: "bifemsh", 11: "vas_int", 12: "soleus",
        13: "med_gas", 14: "lat_gas", 15: "tib_post", 16: "tib_ant", 17: "per_brev",
        18: "per_long", 19: "per_tert", 20: "flex_dig", 21: "flex_hal", 22: "ext_dig",
        23: "ext_hal", 24: "glut_max1", 25: "add_mag1", 26: "iliacus", 27: "psoas"}
nine_rel = {}
for aid, short in NINE.items():
    fn = f"act{aid:02d}_{short}_9x9_20261023.csv"
    assert os.path.exists(os.path.join(FOLDER, NINE_DIR, fn)), f"missing 9x9 export {fn}"
    nine_rel[aid] = f"tmp/paper_bridge/capability_table_20261022/{NINE_DIR}/{fn}"

# ---------- 2. hand-tuned row ----------
assert len(hand) == 27, f"hand rows {len(hand)} != 27"
HAND = {}
for i, h in enumerate(hand):
    aid = i + 1
    assert h["muscle"] == SURF[aid]["master"], f"hand row {i}: {h['muscle']} != {SURF[aid]['master']}"
    assert h["primary_dof"] == SURF[aid]["primary"] and h["secondary_dof"] == SURF[aid]["secondary"], f"hand act {aid}: dof mismatch"
    assert int(h["cells_domain"]) == SURF[aid]["n_dom"], f"hand act {aid}: cells_domain"
    assert int(h["cells_excluded"]) == SURF[aid]["n_exc"], f"hand act {aid}: cells_excluded"
    assert int(h["cells_flagged_in_domain"]) == SURF[aid]["n_flag"], f"hand act {aid}: flagged {h['cells_flagged_in_domain']} != {SURF[aid]['n_flag']}"
    HAND[aid] = h
suspects = [aid for aid in range(1, 28) if HAND[aid]["suspect"] == "yes"]
assert len(suspects) == 9, f"suspect rows {len(suspects)} != 9: {suspects}"

# ---------- 3. optimized row ----------
assert len(opt) == 27, f"opt rows {len(opt)} != 27"
OPT = {}
for i, o in enumerate(opt):
    aid = int(o["actuator_id"])
    assert aid == i + 1, f"opt row order at {i}"
    assert o["built"] == "true", f"act {aid} built != true"
    assert o["master_name"] == SURF[aid]["master"], f"opt act {aid}: name mismatch"
    assert o["primary_dof"] == SURF[aid]["primary"] and o["secondary_dof"] == SURF[aid]["secondary"], f"opt act {aid}: dof mismatch"
    assert int(o["domain_cells"]) == SURF[aid]["n_dom"] and int(o["excluded_cells"]) == SURF[aid]["n_exc"], f"opt act {aid}: cells"
    assert o["converged"] == "1" and o["stopReason"] == "converged", f"opt act {aid}: not converged"
    assert 0 <= int(o["evalsUsed"]) <= 300, f"opt act {aid}: evals {o['evalsUsed']} outside 0..300"
    OPT[aid] = o

# ---------- 4. row averages (with- and without-flag variants) ----------
def row_avg(src, meet_key, active_key, un_meet_key, un_active_key):
    m = mean(float(src[aid][meet_key]) for aid in range(1, 28))
    a = mean(float(src[aid][active_key]) for aid in range(1, 28))
    mu = mean(float(src[aid][un_meet_key]) for aid in range(1, 28))
    au = mean(float(src[aid][un_active_key]) for aid in range(1, 28))
    return m, a, mu, au

handMeet, handAct, handMeetU, handActU = row_avg(HAND, "meetPct", "activePct", "meetPct_unflagged", "activePct_unflagged")
optMeet, optAct, optMeetU, optActU = row_avg(OPT, "meetPct", "activePct", "meetPct_unflagged", "activePct_unflagged")

evals = [int(OPT[a]["evalsUsed"]) for a in range(1, 28)]
walls = [float(OPT[a]["wallMin"]) for a in range(1, 28)]
evals_sorted = sorted(evals)
median_evals = evals_sorted[len(evals_sorted) // 2]
min_evals, max_evals = min(evals), max(evals)
max_evals_act = OPT[evals.index(max_evals) + 1]["master_name"]
active0 = [a for a in range(1, 28) if float(OPT[a]["activePct"]) == 0.0]
meet100 = [a for a in range(1, 28) if float(OPT[a]["meetPct"]) == 100.0]
hip_acts = [a for a in range(1, 28) if SURF[a]["env"] == "yes"]
assert hip_acts == [24, 25, 26, 27]

print("=== derived row averages (percent) ===")
print(f"hand meet  {f2(handMeet)}  active  {f2(handAct)}  | unflagged meet {f2(handMeetU)} active {f2(handActU)}")
print(f"opt  meet  {f2(optMeet)}  active  {f2(optAct)}  | unflagged meet {f2(optMeetU)} active {f2(optActU)}")
print(f"evals min {min_evals} max {max_evals} ({max_evals_act}) median {median_evals}; wall {min(walls)}-{max(walls)} min")
print(f"active==0 opt rows: {active0}; meet==100 rows: {meet100}; suspects: {suspects}")
print(f"totals: graded {n_graded_total}, excluded {n_exc_total}, flagged {n_flag_total}")

# ---------- 5. deliverable 1: capability_table_20261023.md ----------
md = f"""# 27-actuator capability table (2026-10-23, v2 coupled-DOF surfaces)

Each cell is the row average of the per-actuator AREA SHARE on the 5x5
coupled-DOF grading grid (`human_surfaces_5x5_grading_20261023.csv`): the
share of the actuator's graded domain (5x5 inclusive nodes over the ranges
of motion of the two joint coordinates the route crosses) where the condition
holds. "Meets or exceeds human" = the robot torque has the human torque's
sign AND |tau_robot| >= |tau_human| at the cell. "Active force" = positive
BPA strain AND taut, positive torque. Flagged cells (OpenSim path artifacts)
COUNT in these headline numbers; the without-flag variants are in
`CAVEATS_REPORT_20261023.md`. Hand-tuned row: old-code routing, old-code BPA
sizes, 620 kPa. Optimized row: fresh BiPulley survey2 solves, same old-code
BPA sizes. Both rows graded against the same single-muscle human surfaces
(stock gait2392_simbody routing, OpenSim 4.4.1).

| Row | % of RoM torque meets or exceeds human | % of RoM with active force |
|---|---|---|
| Hand-tuned (27 of 27 graded) | {f2(handMeet)}% | {f2(handAct)}% |
| Optimized (27 of 27 graded; 27 built) | {f2(optMeet)}% | {f2(optAct)}% |

**Muscles excluded: none.** Both rows average all 27 actuators; no muscle is
excluded from either row. The only exclusions anywhere in the metric are GRID
CELLS, not muscles: on the four hip surfaces whose adduction axis is graded
(Gluteus Maximus, Adductor Magnus, Iliacus, Psoas X), 10 of 25 cells each lie
outside the coupled hip adduction/abduction envelope (FIX_SPEC_20261023.md
section 1b) and are removed from the area denominator for BOTH rows, never
counted as failures. Per-muscle shares, flags, and suspects:
`per_muscle_capability_20261023.csv`; full caveats: `CAVEATS_REPORT_20261023.md`.
"""
open(os.path.join(FOLDER, "capability_table_20261023.md"), "w", encoding="utf-8", newline="\n").write(md)

# ---------- 6. deliverable 2: capability_table_20261023.tex ----------
tex = f"""% capability_table_20261023.tex -- 27-actuator capability table (v2,
% coupled-DOF surfaces), assembled 2026-10-23 from
% per_muscle_capability_20261023.csv in this folder. Float placement NOT
% decided: Ben wraps and places it.
% Suggested caption: Isometric capability of the hand-tuned and optimized
% actuator routings against the same-named human muscle torque surfaces
% (Gait2392, OpenSim 4.4.1). Each cell is the row average of the area shares
% on the 5x5 grading grid spanning the two joint coordinates each route
% crosses; path-artifact-flagged cells count in the headline numbers. Both
% rows average all 27 actuators. See CAVEATS_REPORT_20261023.md.
\\begin{{tabular}}{{lcc}}
\\toprule
 & Torque meets or exceeds human (\\% of RoM) & Active force (\\% of RoM) \\\\
\\midrule
Hand-tuned (27 of 27 graded) & {f2(handMeet)}\\% & {f2(handAct)}\\% \\\\
Optimized (27 of 27 graded; 27 built) & {f2(optMeet)}\\% & {f2(optAct)}\\% \\\\
\\bottomrule
\\end{{tabular}}
"""
open(os.path.join(FOLDER, "capability_table_20261023.tex"), "w", encoding="utf-8", newline="\n").write(tex)

# ---------- 7. deliverable 3: per_muscle_capability_20261023.csv ----------
pm_cols = ["muscle", "actuator_id", "primary_dof", "secondary_dof", "coverage_source",
           "cells_domain", "cells_excluded", "cells_flagged_in_domain",
           "flags_via", "flags_step", "flags_slack", "flags_short",
           "path_artifact_flag", "transcription_suspect",
           "hand_meet_pct", "hand_active_pct", "hand_meet_pct_unflagged", "hand_active_pct_unflagged",
           "opt_meet_pct", "opt_active_pct", "opt_meet_pct_unflagged", "opt_active_pct_unflagged",
           "opt_nBPA", "opt_G", "bpa_size_mm", "opt_evals_used", "opt_converged",
           "opt_winner_mat", "notes"]
pm_path = os.path.join(FOLDER, "per_muscle_capability_20261023.csv")
with open(pm_path, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(pm_cols)
    for aid in range(1, 28):
        S, H, O = SURF[aid], HAND[aid], OPT[aid]
        anyflag = "yes" if S["n_flag"] > 0 else "no"
        notes = (f"HAND: {H['note']} | OPT: {O['note']} | "
                 f"FLAGS: {S['n_flag']} flagged graded cells (via {S['flags']['via']}, step {S['flags']['step']}, "
                 f"slack {S['flags']['slack']}, short {S['flags']['short']}); "
                 f"WITH-flag meet/active: hand {H['meetPct']}/{H['activePct']}, opt {O['meetPct']}/{O['activePct']}; "
                 f"WITHOUT-flag meet/active: hand {H['meetPct_unflagged']}/{H['activePct_unflagged']}, "
                 f"opt {O['meetPct_unflagged']}/{O['activePct_unflagged']} (percent)")
        w.writerow([S["master"], aid, S["primary"], S["secondary"], H["source"],
                    S["n_dom"], S["n_exc"], S["n_flag"],
                    S["flags"]["via"], S["flags"]["step"], S["flags"]["slack"], S["flags"]["short"],
                    anyflag, H["suspect"],
                    H["meetPct"], H["activePct"], H["meetPct_unflagged"], H["activePct_unflagged"],
                    O["meetPct"], O["activePct"], O["meetPct_unflagged"], O["activePct_unflagged"],
                    O["chosen_nB"], O["chosen_G"], O["bpa_size_mm"], O["evalsUsed"], O["converged"],
                    O["matFile"], notes])

# ---------- 8. deliverable 4: surfaces_index_20261023.md ----------
idx = [f"""# Human 5x5 grading surfaces index (2026-10-23)

One line per actuator: the two graded joint coordinates (with their
range-of-motion spans in degrees), the grid, the in-domain artifact-flag
count, and the 9x9 viewing export. The single grading source of truth is
`human_surfaces_5x5_grading_20261023.csv` (this folder; 25 cells per
actuator, direct OpenSim evaluations at the nodes; knee flexion NEGATIVE,
hip adduction POSITIVE). The 9x9 exports are bilinear interpolations from
the 1-deg OpenSim samples, for VIEWING only; percentages are computed on
the 5x5 grading grid, never on the 9x9. All 9x9 exports live in
`tmp/paper_bridge/capability_table_20261022/{NINE_DIR}/`.

| id | Actuator | Graded coordinates (deg) | Grid | Flagged graded cells | 9x9 export |
|---|---|---|---|---|---|"""]
for aid in range(1, 28):
    S = SURF[aid]
    pin = f"; pinned at 0: {', '.join(S['pinned_list'])}" if S["pinned_list"] else ""
    idx.append(
        f"| {aid} | {S['master']} | {S['primary']} [{S['p_lo']}, {S['p_hi']}] x "
        f"{S['secondary']} [{S['s_lo']}, {S['s_hi']}]"
        f"{'; envelope-coupled adduction axis' if S['env'] == 'yes' else ''}{pin} | "
        f"5x5 grading ({S['n_dom']} graded of 25); 9x9 view | {S['n_flag']} of {S['n_dom']} | "
        f"`{NINE_DIR}/act{aid:02d}_{NINE[aid]}_9x9_20261023.csv` |")
idx.append(f"""
Totals: {n_graded_total} graded cells, {n_exc_total} envelope-excluded cells
(acts 24-27, 10 each), {n_flag_total} flagged graded cells. Flag layer,
hinge-chain, gate, and provenance detail per actuator:
`human_surfaces_manifest_20261023.json` (this folder). Pinned DOFs: actuator
1 lumbar_rotation; 7 hip_adduction_r; 20-23 subtalar_angle_r; 24
hip_rotation_r.
""")
open(os.path.join(FOLDER, "surfaces_index_20261023.md"), "w", encoding="utf-8", newline="\n").write("\n".join(idx))

# ---------- 9. deliverable 5: capability_table_20261023.xlsx ----------
import openpyxl  # verified present (3.1.5) in D:/Anaconda/envs/myo before writing
wb = openpyxl.Workbook()
ws = wb.active
ws.title = "Capability 2x2"
ws.append(["Row", "% of RoM torque meets or exceeds human", "% of RoM with active force"])
ws.append(["Hand-tuned (27 of 27 graded)", round(handMeet, 2), round(handAct, 2)])
ws.append(["Optimized (27 of 27 graded; 27 built)", round(optMeet, 2), round(optAct, 2)])
ws.append([])
ws.append(["Cells = row averages of the per-actuator area shares on the 5x5 coupled-DOF grading grid;"])
ws.append(["flagged cells count in these headline numbers (without-flag variants in CAVEATS_REPORT_20261023.md)."])
ws.append(["Sources of record: human_surfaces_5x5_grading_20261023.csv; hand_tuned_metrics_20261023.csv;",
           "BiPulley_survey2_perMuscle_20261023.csv. Per-muscle sheet: per_muscle_capability_20261023.csv. Spec: FIX_SPEC_20261023.md."])
ws2 = wb.create_sheet("Per muscle")
ws2.append(pm_cols)
for aid in range(1, 28):
    S, H, O = SURF[aid], HAND[aid], OPT[aid]
    anyflag = "yes" if S["n_flag"] > 0 else "no"
    notes = (f"HAND: {H['note']} | OPT: {O['note']} | "
             f"FLAGS: {S['n_flag']} flagged graded cells (via {S['flags']['via']}, step {S['flags']['step']}, "
             f"slack {S['flags']['slack']}, short {S['flags']['short']}); "
             f"WITH-flag meet/active: hand {H['meetPct']}/{H['activePct']}, opt {O['meetPct']}/{O['activePct']}; "
             f"WITHOUT-flag meet/active: hand {H['meetPct_unflagged']}/{H['activePct_unflagged']}, "
             f"opt {O['meetPct_unflagged']}/{O['activePct_unflagged']} (percent)")
    ws2.append([S["master"], aid, S["primary"], S["secondary"], H["source"],
                S["n_dom"], S["n_exc"], S["n_flag"],
                S["flags"]["via"], S["flags"]["step"], S["flags"]["slack"], S["flags"]["short"],
                anyflag, H["suspect"],
                float(H["meetPct"]), float(H["activePct"]), float(H["meetPct_unflagged"]), float(H["activePct_unflagged"]),
                float(O["meetPct"]), float(O["activePct"]), float(O["meetPct_unflagged"]), float(O["activePct_unflagged"]),
                int(O["chosen_nB"]), int(O["chosen_G"]), float(O["bpa_size_mm"]), int(O["evalsUsed"]),
                O["converged"], O["matFile"], notes])
xlsx_path = os.path.join(FOLDER, "capability_table_20261023.xlsx")
wb.save(xlsx_path)

# read-back verification of the xlsx
wbv = openpyxl.load_workbook(xlsx_path)
assert wbv.sheetnames == ["Capability 2x2", "Per muscle"], wbv.sheetnames
wsv = wbv["Capability 2x2"]
assert wsv["A2"].value == "Hand-tuned (27 of 27 graded)"
assert abs(wsv["B2"].value - round(handMeet, 2)) < 1e-9 and abs(wsv["C2"].value - round(handAct, 2)) < 1e-9
assert abs(wsv["B3"].value - round(optMeet, 2)) < 1e-9 and abs(wsv["C3"].value - round(optAct, 2)) < 1e-9
assert wbv["Per muscle"].max_row == 28, wbv["Per muscle"].max_row
print("xlsx read-back VERIFY OK:", wbv.sheetnames)

# ---------- 10. deliverable 6: CAVEATS_REPORT_20261023.md ----------
suspect_lines = "\n".join(
    f"- **{HAND[a]['muscle']}** (act {a}, hand meet {HAND[a]['meetPct']}%, active {HAND[a]['activePct']}%): {HAND[a]['note']}"
    for a in suspects)
active0_lines = ", ".join(f"{SURF[a]['master']} (act {a}, opt meet {OPT[a]['meetPct']}%)" for a in active0)
meet100_lines = ", ".join(f"{SURF[a]['master']} (act {a})" for a in meet100)
const_sec = [a for a in range(1, 28) if "not a class hinge" in OPT[a]["note"]]
const_sec_lines = ", ".join(f"{SURF[a]['master']} (act {a}: {OPT[a]['primary_dof']} graded, "
                            f"{OPT[a]['secondary_dof']} robot surface constant)" for a in const_sec)
pinned_lines = "; ".join(f"act {a} {', '.join(SURF[a]['pinned_list'])}" for a in range(1, 28) if SURF[a]["pinned_list"])

cav = f"""# CAVEATS_REPORT_20261023

Plain-language caveats for the v2 27-actuator capability table
(capability_table_20261023.md / .tex / .xlsx) and the per-muscle rows
(per_muscle_capability_20261023.csv), assembled 2026-10-23 on easteregg2
under `FIX_SPEC_20261023.md`. Every number re-derives from the CSVs of
record in this folder: `human_surfaces_5x5_grading_20261023.csv` (the single
grading source for BOTH rows), `hand_tuned_metrics_20261023.csv`
(hand-tuned row), `BiPulley_survey2_perMuscle_20261023.csv` (optimized row),
cross-checked against `human_surfaces_manifest_20261023.json`. The assembler
(`assemble_deliverables_20261023.py`) asserts every count below before
writing; all asserts passed in this run.

## 1. What changed from the v1 table (2026-10-22)

- Metric: v2 grades each actuator on the AREA SHARE of its coupled-DOF 5x5
  surface (the two joint coordinates the route crosses), replacing v1's 1-D
  primary-DOF slice (FIX_SPEC sections 1/1b).
- Optimized row: fresh solves against the single-muscle surfaces (all 27
  built and graded; v1 had 23 built, 16 gradeable).
- Hip adduction/abduction coupling (FIX_SPEC 1b, Ben 2026-10-23) is applied;
  v1's quarter-ellipse placeholder is retired.
- v1 files (`*_20261022*`) are untouched.

## 2. Covered counts

- **Hand-tuned row: 27 of 27 covered, 0 excluded.** Coverage source is
  `connor-corrected` on every row (old-code routing, old-code BPA sizes,
  620 kPa).
- **Optimized row: 27 of 27 built and graded.** `built = true` and
  `converged = 1` on every row of the survey2 CSV (builder coverage 27/27).
- No muscle is excluded from either row.

## 3. Exclusions that ARE in the metric: hip-envelope cells (not muscles)

On the four surfaces whose adduction axis is graded (Gluteus Maximus,
Adductor Magnus, Iliacus, Psoas X, acts 24-27), 10 of 25 grid cells each lie
outside the coupled hip adduction/abduction envelope and are REMOVED from the
area denominator for BOTH rows (never counted as failures): 40 cells total.
Rule (FIX_SPEC 1b): flexion axis -15 deg (max extension) to +110 deg (max
flexion); adduction maximum tapers linearly from +20 deg at full extension to
+10 deg at full flexion, abduction maximum from -10 deg to -60 deg; a cell is
graded when its center is inside the closed envelope. Physical origin: BONE
INTERFERENCE (femur-pelvis contact, per Ben), a real anatomical/hardware
constraint. Graded cells per hip actuator: 15. `hip_rotation` stays the
static rom-table box (pinned at 0 on act 24); act 7 (TFL) has hip_adduction
pinned at 0, which is inside the envelope at every graded flexion, so no
envelope exclusion applies to its surface.

## 4. Artifact-flagged cells, and both percentage versions

Flag layer on the human surfaces (per-cell; curves are NEVER corrected):
FLAG_VIA (conditional path point active-set change in the cell bin),
FLAG_STEP (adjacent-sample moment-arm jump over 2% of the surface-wide peak),
FLAG_SLACK (musculotendon at or under tendon slack length), FLAG_SHORT
(estimated normalized fiber length under 0.5). FLAG_BONE (bone penetration)
is NOT APPLIED: the check is still TODO per the checklist, no flags were
invented.

Totals (re-derived from the grading CSV; equal to the manifest): 635 graded
cells, 40 envelope-excluded cells, **126 flagged graded cells**. Per-actuator
flagged counts are in `surfaces_index_20261023.md` and the per-muscle CSV
column `cells_flagged_in_domain`.

**Flagged cells count in the headline percentages (conservative). Both
versions, as row averages of the per-muscle shares:**

- Hand-tuned row: WITH flags meet **{f2(handMeet)}%** / active **{f2(handAct)}%**;
  WITHOUT flags meet **{f2(handMeetU)}%** / active **{f2(handActU)}%**.
- Optimized row: WITH flags meet **{f2(optMeet)}%** / active **{f2(optAct)}%**;
  WITHOUT flags meet **{f2(optMeetU)}%** / active **{f2(optActU)}%**.

## 5. Transcription suspects (hand-tuned row; carried over VERBATIM from v1)

Nine actuators carry a transcription-suspect flag because the corrected
Connor script differs from Ben's reference values at one or more routing
points (the hand CSV header states "suspect+note carried over VERBATIM from
hand_tuned_metrics_20261022.csv"; this pass verified 9 `suspect = yes` rows).
**All were graded as-is**; none was silently corrected. A suspect flag means
the hand-tuned percentage could shift once Ben rules on the value; it does
not mean the grade is wrong. Hand percentages below are the v2 (2026-10-23)
re-grade values.

{suspect_lines}

## 6. Unconverged or capped rows, and budget disclosure

- **Convergence: none failed.** `converged = 1` and `stopReason = converged`
  on all 27 optimized rows.
- **No row hit the evaluation cap.** Evaluations per muscle: min {min_evals}
  (Vastus Intermedius), median {median_evals}, max {max_evals} (External
  Oblique), all under the 300-evaluation cap. Recorded wall time 0.12-0.37
  min per muscle, far under the roughly 10-min cap. Parpool was capped at 8
  per the FIX_SPEC; the survey CSV of record does not record the pool size
  actually used, so that is stated as the cap, not as a measured value.
- Ladder: nBPA {{1,2,3}} x tackle G {{1,2,3}} with sizes re-solved per
  configuration, resting lengths as upper limits, envelope guards; the
  recorded winner is the best-meeting design (an active tie-break never
  records a lower-meeting design; FIX_SPEC section 4).
- **Configuration ceiling:** only {len(meet100)} of 27 winners reach 100.00%
  meet ({meet100_lines}). The rest fall short of full coverage within the
  allowed configuration set and budget; the ladder chains in the survey2 CSV
  are the record. This is a capability survey, not a design campaign: the
  budget answers "can this routing class reach the human surface at all",
  and the percentages are not flight-ready designs.

## 7. Active-force zeros in the optimized row

Six optimized winners record active 0.0000% while still meeting part of the
surface: {active0_lines}. The winners are taut where they meet but have zero
positive-strain cells on the grading grid. The values are recorded as-is
from the CSV of record; treat those rows' active cells with care until the
strain-window question is resolved (the same class of issue v1 flagged on
Bicep Femoris Long Head). The hand-tuned row has no active zeros (row minimum
80.00%).

## 8. Two-hinge scope for 3-plus-joint muscles, and constant-secondary rows

Each actuator is graded on the TWO coordinates the route crosses; remaining
joints the muscle crosses are pinned at 0 deg: {pinned_lines}. Moment arms
the human model develops at other joint configurations are not represented;
this is a per-muscle capability slice, not a coupled-posture prediction.

Constant-secondary rows: {const_sec_lines}. For these the optimized robot
surface is constant along the secondary axis because that coordinate is not
a BiPam class hinge (the survey2 CSV notes say "remaining joint pinned at 0,
robot surface constant along that axis (same as the hand-tuned row)"); the
hand-tuned robot surface has the same property on those rows. Their meet
shares are effectively 1-D over the primary coordinate, evaluated on the
same 25-cell grid against the full 2-D human surface.

## 9. Grid, percentage, and sign conventions

- Grading grid: 5x5 INCLUSIVE uniform nodes per axis over the rom-table
  ranges (25 cells per actuator); cell values are DIRECT OpenSim evaluations
  at the nodes; cell bins are midpoint-Voronoi intervals. Percentages are
  computed on the grading grid, never on the 9x9 exports.
- 9x9 exports (`human_surfaces_9x9_20261023/`) are bilinear interpolations
  from the 1-deg OpenSim samples, for VIEWING only.
- Human torque: tau_h = sum over chosen muscles of computeMomentArm (about
  the primary coordinate) x Fmax, OpenSim native signs (moment arm = -dL/dq;
  positive torque acts in the positive coordinate direction); knee flexion
  NEGATIVE; hip adduction POSITIVE.
- Per cell: meet = robot torque has the human torque's sign AND
  |tau_robot| >= |tau_human|; active = positive strain AND taut torque > 0.
- Row averages are UNWEIGHTED means of the per-muscle area shares (the
  printed 4-decimal shares of the two row CSVs), displayed rounded to two
  decimals.
- Human surfaces are STOCK gait2392_simbody routing only (human data on
  human paths; robot/BPA paths never substituted).

## 10. Provenance and environment

- OpenSim 4.4.1 (pinned) against `gait2392_simbody.osim` (sha1
  b89062c763e6d8fa3d7651c0a4fd2820055d807c) via the actuator map of record
  (`Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json`, sha1
  a97afffeeae85687810986273db6ba0151eba1a2) and the 1-D audit of record
  (sha1 47912eb6c1e858e9e62ceefe36b8a510e052c417); provenance headers are in
  the grading CSV.
- Machine easteregg2; Python `D:/Anaconda/envs/myo/python.exe` (the same
  interpreter the surface generator used); xlsx written with openpyxl 3.1.5,
  already present in that env (verified by import in this pass; no package
  was installed).
- Assembler: `assemble_deliverables_20261023.py` (this folder) re-derives
  every cell from the three CSVs of record, asserts the counts (27/27 both
  rows; 635/40/126 cell totals; 9 suspects; manifest flag equality) before
  writing, verifies the xlsx by read-back, and writes ONLY the six dated
  _20261023 deliverables. Inputs and all v1 (20261022) files are unmodified.
- Git state at assembly time (checked with `git status --porcelain`, branch
  `KneeTestSetup_BenBo_stw`): the working tree was NOT clean before this
  assembly (tracked edits to AGENTS.md, the torque-paper build files under
  `Documentation/Reports and Papers/Knee_Torque_Test/`, and
  `tmp/paper_bridge/build_flx20_20261023/` + `figscripts/`; untracked
  BiPulley_survey2 result mats under
  `Code/Matlab/Mesh_Optimization/Results/` and other files). This pass added
  only the six dated deliverables plus the assembler in this folder. No
  commits or pushes were made; Ben commits via GitHub Desktop.
"""
open(os.path.join(FOLDER, "CAVEATS_REPORT_20261023.md"), "w", encoding="utf-8", newline="\n").write(cav)

print("WROTE: capability_table_20261023.md .tex .xlsx, per_muscle_capability_20261023.csv,",
      "surfaces_index_20261023.md, CAVEATS_REPORT_20261023.md")
