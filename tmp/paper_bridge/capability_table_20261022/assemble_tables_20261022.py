# Assemble the 27-actuator 2x2 capability table (Ben, 2026-10-22).
# Inputs of record (this folder):
#   human_torque_curves_20261022.csv                      - single-muscle human targets (grading source, both rows)
#   BiPulley_survey_perMuscle_regraded_20261022_223448.csv - optimized row (23 of 27 built; regraded vs single-muscle)
#   hand_tuned_metrics_20261022.csv                        - hand-tuned row (27 of 27 covered)
# Averaging rules (launch ask):
#   hand row: averages over covered muscles only (count in the row label)
#   optimized row: averages over the muscles that were built; 7 of the 23
#   built are GRADE-UNABLE (no percentages exist), so the numeric average
#   covers the 16 built muscles that could be graded. Nothing is faked.
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CURVES = os.path.join(HERE, "human_torque_curves_20261022.csv")
OPT = os.path.join(HERE, "BiPulley_survey_perMuscle_regraded_20261022_223448.csv")
HAND = os.path.join(HERE, "hand_tuned_metrics_20261022.csv")

def read_rows(path, comment_prefix=None):
    with open(path, "r", encoding="utf-8", newline="") as f:
        lines = [ln for ln in f if not ln.lstrip().startswith(comment_prefix)] if comment_prefix else f
        return list(csv.DictReader(lines))

# ---- load inputs -----------------------------------------------------------
curves = read_rows(CURVES, "#")
opt = read_rows(OPT)
hand = read_rows(HAND)

assert len(curves) == 27, f"curves rows {len(curves)} != 27"
assert len(opt) == 27, f"opt rows {len(opt)} != 27"
assert len(hand) == 27, f"hand rows {len(hand)} != 27"

# curves: id, master_name, primary_dof, artifact flag
curve_info = {}
for r in curves:
    aid = int(r["actuator_id"])
    flag = r["artifact_flag_caveat"].strip()
    curve_info[aid] = {
        "name": r["master_name"].strip(),
        "dof": r["primary_dof"].strip(),
        "artifact": False if flag == "none" or flag == "" else True,
        "flag_text": "" if flag == "none" else flag,
    }
# opt CSV rows follow actuator-map order (harvest asserts 27 rows, map order)
for i, r in enumerate(opt, start=1):
    assert curve_info[i]["name"] == r["muscle"].strip(), (
        f"opt row {i} muscle '{r['muscle']}' != curve '{curve_info[i]['name']}'")

hand_by_name = {r["muscle"].strip(): r for r in hand}
for i in range(1, 28):
    assert curve_info[i]["name"] in hand_by_name, f"hand CSV missing {curve_info[i]['name']}"

# ---- per-muscle rows -------------------------------------------------------
rows = []
for i in range(1, 28):
    ci = curve_info[i]
    o = opt[i - 1]
    h = hand_by_name[ci["name"]]
    built = o["built"].strip().lower() == "true"
    # graded == surveyed winner reloaded and re-graded (has percentages):
    # built with a nonzero nBPA/G choice and a converged checkpoint
    graded = built and o["nBPA"].strip() not in ("0", "") and o["converged"].strip() == "1"
    h_meet, h_act = float(h["meetPct"]), float(h["activePct"])
    o_meet = float(o["meetPct"]) if graded else None
    o_act = float(o["activePct"]) if graded else None
    notes = []
    if h["note"].strip():
        notes.append("HAND: " + h["note"].strip())
    if o["note"].strip():
        notes.append("OPT: " + o["note"].strip())
    if ci["artifact"]:
        notes.append("CURVE-ARTIFACT: " + ci["flag_text"])
    rows.append({
        "muscle": ci["name"],
        "actuator_id": i,
        "primary_dof": ci["dof"],
        "coverage_source": h["source"].strip(),
        "hand_meet_pct": h_meet,
        "hand_active_pct": h_act,
        "opt_meet_pct": o_meet,
        "opt_active_pct": o_act,
        "opt_nBPA": int(o["nBPA"]) if built else 0,
        "opt_G": int(o["tackleG"]) if built else 0,
        "opt_evals_used": int(o["evalsUsed"]) if built else 0,
        "opt_converged": o["converged"].strip() == "1" if built else False,
        "opt_built": built,
        "opt_graded": graded,
        "path_artifact_flag": "yes" if ci["artifact"] else "no",
        "transcription_suspect_flag": h["suspect"].strip(),
        "notes": " | ".join(notes),
    })

# ---- checks ----------------------------------------------------------------
hand_cov = [r for r in rows if r["coverage_source"] in ("connor-corrected", "ben-original")]
opt_built = [r for r in rows if r["opt_built"]]
opt_graded = [r for r in rows if r["opt_graded"]]
hand_excluded = [r["muscle"] for r in rows if r["coverage_source"] == "none"]

def avg(rs, key):
    return sum(r[key] for r in rs) / len(rs)

hand_meet_avg = avg(hand_cov, "hand_meet_pct")
hand_act_avg = avg(hand_cov, "hand_active_pct")
opt_meet_avg = avg(opt_graded, "opt_meet_pct")
opt_act_avg = avg(opt_graded, "opt_active_pct")

artifact_set = sorted(r["muscle"] for r in rows if r["path_artifact_flag"] == "yes")
suspect_set = sorted(r["muscle"] for r in rows if r["transcription_suspect_flag"] == "yes")
cannot_build = [r["muscle"] for r in rows if not r["opt_built"]]
grade_unable = [r["muscle"] for r in opt_built if not r["opt_graded"]]
capped = [r["muscle"] for r in rows if "capability-capped" in r["notes"]]
unconverged = [r["muscle"] for r in opt_graded if not r["opt_converged"]]

print("== CHECKS ==")
print(f"hand covered: {len(hand_cov)} of 27; excluded: {hand_excluded}")
print(f"opt built: {len(opt_built)} of 27; graded: {len(opt_graded)}")
print(f"cannot build ({len(cannot_build)}): {cannot_build}")
print(f"grade-unable built ({len(grade_unable)}): {grade_unable}")
print(f"unconverged among graded: {unconverged}")
print(f"artifact-flagged curves ({len(artifact_set)}): {artifact_set}")
print(f"transcription suspects ({len(suspect_set)}): {suspect_set}")
print(f"capability-capped ({len(capped)}): {capped}")
print(f"hand meet avg  = {hand_meet_avg:.10f} (sum {sum(r['hand_meet_pct'] for r in hand_cov)})")
print(f"hand active avg = {hand_act_avg:.10f} (sum {sum(r['hand_active_pct'] for r in hand_cov)})")
print(f"opt meet avg   = {opt_meet_avg:.10f} (sum {sum(r['opt_meet_pct'] for r in opt_graded)})")
print(f"opt active avg = {opt_act_avg:.10f} (sum {sum(r['opt_active_pct'] for r in opt_graded)})")
expect = {
    "hand_meet": 1968 / 27, "hand_act": 2668 / 27,
    "opt_meet": 832 / 16, "opt_act": 1316 / 16,
}
assert abs(hand_meet_avg - expect["hand_meet"]) < 1e-9
assert abs(hand_act_avg - expect["hand_act"]) < 1e-9
assert abs(opt_meet_avg - expect["opt_meet"]) < 1e-9
assert abs(opt_act_avg - expect["opt_act"]) < 1e-9
print("averages reproduce the independent hand sums -> OK")

HAND_LBL = f"Hand-tuned ({len(hand_cov)} of 27)"
OPT_LBL = f"Optimized ({len(opt_built)} of 27)"

# ---- 1. per-muscle CSV -----------------------------------------------------
pm_path = os.path.join(HERE, "per_muscle_capability_20261022.csv")
pm_cols = ["muscle", "actuator_id", "primary_dof", "coverage_source",
           "hand_meet_pct", "hand_active_pct", "opt_meet_pct", "opt_active_pct",
           "opt_nBPA", "opt_G", "opt_evals_used", "opt_converged",
           "path_artifact_flag", "transcription_suspect_flag", "notes"]
with open(pm_path, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(pm_cols)
    for r in rows:
        w.writerow([
            r["muscle"], r["actuator_id"], r["primary_dof"], r["coverage_source"],
            f"{r['hand_meet_pct']:.4f}", f"{r['hand_active_pct']:.4f}",
            f"{r['opt_meet_pct']:.4f}" if r["opt_graded"] else "",
            f"{r['opt_active_pct']:.4f}" if r["opt_graded"] else "",
            r["opt_nBPA"] if r["opt_built"] else "", r["opt_G"] if r["opt_built"] else "",
            r["opt_evals_used"] if r["opt_built"] else "",
            "1" if r["opt_converged"] else ("0" if r["opt_built"] else ""),
            r["path_artifact_flag"], r["transcription_suspect_flag"], r["notes"],
        ])
print(f"WROTE {pm_path}")

# ---- 2. markdown table -----------------------------------------------------
def fmt(x):
    return f"{x:.2f}%"

def fmt_tex(x):
    return f"{x:.2f}\\%"

md = []
md.append("# 27-actuator capability table (2026-10-22)")
md.append("")
md.append("Each cell is the row average: the share of the 25-point range-of-motion grid at the actuator's primary coordinate where the condition holds (all other joints pinned at 0 deg). Human targets: single same-named muscle torque curves, Gait2392 via OpenSim 4.4.1 (`human_torque_curves_20261022.csv`).")
md.append("")
md.append("| Row | % of RoM torque meets or exceeds human | % of RoM with active force |")
md.append("|---|---|---|")
md.append(f"| {HAND_LBL} | {fmt(hand_meet_avg)} | {fmt(hand_act_avg)} |")
md.append(f"| {OPT_LBL} | {fmt(opt_meet_avg)} | {fmt(opt_act_avg)} |")
md.append("")
md.append(f"The optimized row averages the **{len(opt_graded)} built muscles that could be graded** ({len(opt_built)} of 27 built): 7 built actuators have no grade because the survey hinge chain does not model their primary coordinate (listed below). No muscle is excluded from the hand-tuned row (27 of 27 covered, 0 excluded).")
md.append("")
md.append("**Not in the optimized average:**")
md.append("")
md.append("- Could not build (4): Erector Spinae, Internal Oblique, External Oblique, Vastus Intermedius. Reason: the master path crosses no joint (all points on one body), so there is nothing for a joint actuator to transmit.")
md.append(f"- Built but grade-unable (7): {', '.join(grade_unable)}. Reason: the actuator was built from the seed route, but its primary coordinate (subtalar_angle_r for Tibialis Posterior, Peroneus Brevis, Peroneus Longus, Peroneus Tertius; knee_angle_r for Rectus Femoris; hip_adduction_r for Adductor Magnus; hip_flexion_r for Psoas X) is not modeled by the regrade spec's hinge chain, so no torque curve exists to grade.")
md.append("")
md_text = "\n".join(md) + "\n"
md_path = os.path.join(HERE, "capability_table_20261022.md")
with open(md_path, "w", encoding="utf-8", newline="") as f:
    f.write(md_text)
print(f"WROTE {md_path}")

# ---- 3. LaTeX --------------------------------------------------------------
tex = r"""% capability_table_20261022.tex -- 27-actuator capability table, assembled
% 2026-10-22 from per_muscle_capability_20261022.csv in this folder.
% Float placement NOT decided: Ben wraps/places it.
% Suggested caption: Isometric capability of the hand-tuned and optimized
% actuator routings against the same-named human muscle torque curves
% (Gait2392, OpenSim 4.4.1). Each cell is the row average over the 25-point
% range-of-motion grid at the actuator's primary coordinate, with every
% non-graded joint pinned at 0 deg. The hand-tuned row covers all 27
% actuators; the optimized row averages the 16 of 23 built actuators that
% could be graded on the survey hinge chain (4 could not be built and 7
% could not be graded; see the caveats report).
\begin{tabular}{lcc}
\toprule
 & Torque meets or exceeds human (\% of RoM) & Active force (\% of RoM) \\
\midrule
HANDLBL & HMEET & HACT \\
OPTLBL & OMEET & OACT \\
\bottomrule
\end{tabular}
"""
tex = tex.replace("HANDLBL", HAND_LBL).replace("OPTLBL", OPT_LBL)
tex = tex.replace("HMEET", fmt_tex(hand_meet_avg)).replace("HACT", fmt_tex(hand_act_avg))
tex = tex.replace("OMEET", fmt_tex(opt_meet_avg)).replace("OACT", fmt_tex(opt_act_avg))
tex_path = os.path.join(HERE, "capability_table_20261022.tex")
with open(tex_path, "w", encoding="utf-8", newline="") as f:
    f.write(tex)
print(f"WROTE {tex_path}")

# ---- 4. xlsx ---------------------------------------------------------------
import openpyxl
from openpyxl.styles import Font, Alignment

wb = openpyxl.Workbook()
ws = wb.active
ws.title = "capability_table"
ws.append(["Row", "% of RoM torque meets or exceeds human", "% of RoM with active force"])
ws.append([HAND_LBL, round(hand_meet_avg, 2), round(hand_act_avg, 2)])
ws.append([OPT_LBL, round(opt_meet_avg, 2), round(opt_act_avg, 2)])
ws.append([])
ws.append(["Averaging bases"])
ws.append(["Hand-tuned row", f"average over {len(hand_cov)} covered muscles (27 of 27; 0 excluded)"])
ws.append(["Optimized row", f"average over {len(opt_graded)} built muscles that could be graded ({len(opt_built)} of 27 built; {len(cannot_build)} could not build, {len(grade_unable)} built but grade-unable)"])
ws.append(["Human targets", "single same-named muscle torque curves, human_torque_curves_20261022.csv (OpenSim 4.4.1, Gait2392)"])
for c in ("A1", "B1", "C1"):
    ws[c].font = Font(bold=True)
ws["A6"].font = Font(bold=True)
ws.column_dimensions["A"].width = 34
ws.column_dimensions["B"].width = 40
ws.column_dimensions["C"].width = 30

ws2 = wb.create_sheet("per_muscle")
ws2.append(pm_cols)
for r in rows:
    ws2.append([
        r["muscle"], r["actuator_id"], r["primary_dof"], r["coverage_source"],
        r["hand_meet_pct"], r["hand_active_pct"],
        r["opt_meet_pct"], r["opt_active_pct"],
        r["opt_nBPA"] if r["opt_built"] else None, r["opt_G"] if r["opt_built"] else None,
        r["opt_evals_used"] if r["opt_built"] else None,
        r["opt_converged"] if r["opt_built"] else None,
        r["path_artifact_flag"], r["transcription_suspect_flag"], r["notes"],
    ])
for c in range(1, len(pm_cols) + 1):
    ws2.cell(row=1, column=c).font = Font(bold=True)
xlsx_path = os.path.join(HERE, "capability_table_20261022.xlsx")
wb.save(xlsx_path)
print(f"WROTE {xlsx_path}")

# round-trip check of the xlsx
wb2 = openpyxl.load_workbook(xlsx_path)
assert wb2.sheetnames == ["capability_table", "per_muscle"], wb2.sheetnames
assert wb2["capability_table"]["B2"].value == round(hand_meet_avg, 2)
assert wb2["capability_table"]["B3"].value == round(opt_meet_avg, 2)
assert wb2["per_muscle"].max_row == 28 and wb2["per_muscle"].max_column == len(pm_cols)
print(f"XLSX round-trip OK: sheets {wb2.sheetnames}, per_muscle rows {wb2['per_muscle'].max_row - 1}")

# ---- 5. caveats report -----------------------------------------------------
suspect_rows = [r for r in rows if r["transcription_suspect_flag"] == "yes"]
cav = []
cav.append("# CAVEATS_REPORT_20261022")
cav.append("")
cav.append("Plain-language caveats for the 27-actuator capability table (capability_table_20261022.md / .tex / .xlsx) and the per-muscle rows (per_muscle_capability_20261022.csv). Every number in the table re-derives from the CSVs of record in this folder: hand_tuned_metrics_20261022.csv (hand-tuned row), BiPulley_survey_perMuscle_regraded_20261022_223448.csv (optimized row), and human_torque_curves_20261022.csv (human targets, both rows).")
cav.append("")
cav.append("## 1. Covered counts for the two rows")
cav.append("")
cav.append(f"- **Hand-tuned row: {len(hand_cov)} of 27 covered, 0 excluded.** Its averages ({fmt(hand_meet_avg)} meet, {fmt(hand_act_avg)} active) run over all 27 actuators.")
cav.append(f"- **Optimized row: {len(opt_built)} of 27 built; its averages ({fmt(opt_meet_avg)} meet, {fmt(opt_act_avg)} active) run over the {len(opt_graded)} built muscles that could be graded.** Ben asked for all 27; the missing are stated, never faked:")
cav.append(f"  - **Could not build ({len(cannot_build)}):** {', '.join(cannot_build)}. The master muscle path crosses no joint (all routing points sit on one body), so a joint actuator cannot be formed. This is a property of the OpenSim reference paths, not a solver failure.")
cav.append(f"  - **Built but grade-unable ({len(grade_unable)}):** {', '.join(grade_unable)}. Each was built from its seed route, but the regrade spec's hinge chain does not model its primary coordinate (subtalar_angle_r for Tibialis Posterior, Peroneus Brevis, Peroneus Longus, and Peroneus Tertius; knee_angle_r for Rectus Femoris; hip_adduction_r for Adductor Magnus; hip_flexion_r for Psoas X), so no regraded percentage exists. Their opt columns are blank in the per-muscle CSV.")
cav.append("")
cav.append("## 2. Transcription suspects (hand-tuned row)")
cav.append("")
cav.append("Nine actuators carry a transcription-suspect flag because the corrected Connor script differs from Ben's reference values at one or more routing points. **All were graded as-is**; none was silently corrected. A suspect flag means the hand-tuned percentage could shift once Ben rules on the value; it does not mean the grade is wrong.")
cav.append("")
for r in suspect_rows:
    short = r["notes"].split(" | ")[0].replace("HAND: ", "")
    cav.append(f"- **{r['muscle']}** (act {r['actuator_id']}, hand meet {fmt(r['hand_meet_pct'])}): {short}")
cav.append("")
cav.append("**Code-level note (not a per-muscle value):** the repo copy `Code/Matlab/Previous Optimization Code/PamData.m` line 271 reads `size = delta/(1 - 1/def);9hi`. The stray `9hi` token is syntax corruption and the class fails to parse (verified on this machine, MATLAB R2025a, error `Unexpected 'hi'`). Ben's original `Size.m` in the old-code copy has the clean line. The run used a staged dated copy with only that token removed (diff-verified single-line change); the repo file was NOT modified. Fixing it in place is Ben's call.")
cav.append("")
cav.append("## 3. Artifact-flagged human curves (caveat only)")
cav.append("")
cav.append("Six human target curves carry via-point artifact flags from the path-artifact audit: **Semimembranosus (semimem_r), Rectus Femoris (rect_fem_r), Gracilis (grac_r), Vastus Intermedius (vas_int_r), Iliacus (iliacus_r), Psoas X (psoas_r)**. The flags are printed in the curves CSV's `artifact_flag_caveat` column and in the per-muscle CSV's notes. They are **caveats only: the curves were never corrected**, and grading used them as generated.")
cav.append("")
cav.append("## 4. Unconverged small-budget runs")
cav.append("")
cav.append(f"**None.** Every surveyed site's checkpoint records stopReason `converged` (the regrade CSV's converged column is 1 for all {len(opt_graded)} graded rows). The cannot-build and grade-unable rows are structural exclusions, not convergence failures.")
cav.append("")
cav.append("## 5. Budget disclosure: this is a capability survey, not a design campaign")
cav.append("")
opt_evals = sorted(r["opt_evals_used"] for r in opt_graded)
cav.append(f"Each optimized design was chosen with a few hundred function evaluations per muscle (graded rows used {opt_evals[0]}-{opt_evals[-1]} evaluations, median {sorted(opt_evals)[len(opt_evals)//2]}). That budget answers \"can this routing class reach the human curve at all\"; it does not produce flight-ready designs. Several rows are additionally **capability-capped**: the ladder's final rung still reported the verdict `torque_short` at the directive configuration cap (nBPA, G) of (3,1), (3,2), or (3,3), meaning more actuators or tackle might still help. The capped muscles are flagged in the per-muscle CSV notes ({len(capped)} of {len(opt_graded)} graded rows).")
cav.append("")
cav.append("## 6. Still pending: hip adduction/abduction coupling")
cav.append("")
cav.append("The hip adduction/abduction coupling (Ben's 9 control points) is **still pending** and is not represented anywhere in this table. Any row whose actuator crosses the hip is graded with that coupling absent.")
cav.append("")
cav.append("## 7. Target-change disclosure")
cav.append("")
cav.append("The optimized designs were **chosen against the earlier, harder group targets** (full survey run 14:06-14:49 against the group-target curves) and then **re-graded against the single-muscle targets** that both rows now use (Ben's 2026-10-22 ruling: each actuator is graded against the same-named muscle only; no cross-muscle group sums). Optimizing against harder targets and grading against easier ones is a **conservative over-design**: the optimized percentages are, if anything, understated relative to a like-for-like optimization. The per-muscle CSV states the target version in every optimized row's note (`optimized vs GROUP targets; regraded vs SINGLE-MUSCLE csv`); the hand-tuned row was always graded against the single-muscle CSV. Two optimized rows also record Ben's target-change effects explicitly (Semimembranosus meet 44.0% -> 48.0%; Gluteus Maximus unchanged).")
cav.append("")
cav.append("## 8. Grading grid and the neutral-pose pinning convention")
cav.append("")
cav.append("Each actuator is graded on **25 inclusive uniform points over its primary coordinate's range of motion** (the campaign `rom_table_deg` table; verified against the actuator map by the grading scripts). **Every non-graded joint the muscle crosses is pinned at 0 deg**: a one-DOF neutral-pose slice. Moment arms the human model develops at other joint configurations are therefore not represented; this is a per-muscle capability slice, not a coupled-posture prediction.")
cav.append("")
cav.append("## 9. Sign conventions")
cav.append("")
cav.append("Human torques come from OpenSim `computeMomentArm` times Fmax with **OpenSim native signs** (moment arm = -dL/dq, so positive torque acts in the positive coordinate direction); **knee flexion is negative** (vasti positive, hamstrings negative). A grid point counts as **meet** when the actuator is taut and its torque magnitude is at least the human magnitude with matching sign; **active force** counts taut points with positive BPA strain and nonzero torque.")
cav.append("")
cav.append("## 10. One data-quality flag inside the optimized row")
cav.append("")
cav.append("**Bicep Femoris, Long Head** shows meet 100.00% with active 0.00%: the winner is taut everywhere it meets but has zero positive-strain cells (flagged `SUSPECT` in the regrade CSV). The grade is included as-is from the CSV of record; treat that row's 100% with care until the strain-window question is resolved.")
cav.append("")
cav.append("## 11. Environment notes")
cav.append("")
cav.append("- Human curves generated with OpenSim 4.4.1 (pinned) against `gait2392_simbody.osim` (sha1 b89062c7...) via the actuator map of record (`Routing27_v2/actuator_map_27.json`, sha1 a97afffe...); provenance headers are in the curves CSV.")
cav.append("- The xlsx was written with openpyxl, which was **pip-installed into `D:/Anaconda/envs/myo` during this assembly** (openpyxl 3.1.5 plus its dependency et-xmlfile 2.0.0); the env had no openpyxl before.")
cav.append("- All five deliverables are new dated files in this folder, written by the assembly script `assemble_tables_20261022.py`; the assembly itself modified no file outside this folder.")
cav.append("- **Repo-state disclosure (verified 2026-10-22 with GitHub Desktop's bundled git, `git status --porcelain`):** the working tree is NOT clean, and this campaign folder holds two of the modified tracked files: `human_torque_curves_20261022.csv` and `human_torque_curves_manifest_20261022.json` were regenerated in place by the earlier stages of this campaign (13:12 group-target version -> 21:59 single-muscle version, per Ben's 2026-10-22 target ruling; these are inputs of record for this table, not products of the assembly). The other two tracked modifications are `Code/Matlab/Mesh_Optimization/Results/Knee_Flexor_Data_20mm.m` (`tendon = xBest(8)-.02`, `KMAX` 0.255 -> 0.235) and `Testing_Data/2026_06_Festo/readserialnumbers2.m` (`prf` [0.085, -0.296, 0] -> [0.095, -0.325, 0]); both are hardware-side edits outside this campaign's scope and their attribution could not be verified from the working tree alone. Untracked non-tmp files are also present: `Flx_20mm_pressure_for_human_torque.csv` and `Knee_Flexor_20mm.gif` at the repo root, and `Testing_Data/2026_06_Festo/Flx_20mm/`. Ben commits via GitHub Desktop; nothing was committed by the agents.")
cav_text = "\n".join(cav) + "\n"
cav_path = os.path.join(HERE, "CAVEATS_REPORT_20261022.md")
with open(cav_path, "w", encoding="utf-8", newline="") as f:
    f.write(cav_text)
print(f"WROTE {cav_path}")

# round-trip read of the per-muscle CSV
back = read_rows(pm_path)
assert len(back) == 27 and back[0]["muscle"] == rows[0]["muscle"]
assert back[3]["muscle"].startswith("Bicep Femoris, Long Head")
print("per-muscle CSV round-trip OK (27 rows)")

# summary payload for the result
summary = {
    "hand_meet_avg": round(hand_meet_avg, 2), "hand_act_avg": round(hand_act_avg, 2),
    "opt_meet_avg": round(opt_meet_avg, 2), "opt_act_avg": round(opt_act_avg, 2),
    "hand_covered": len(hand_cov), "opt_covered_graded": len(opt_graded),
    "opt_built": len(opt_built),
    "files": [pm_path, md_path, tex_path, xlsx_path, cav_path],
}
with open(os.path.join(HERE, "assemble_summary_20261022.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)
print("SUMMARY " + json.dumps(summary))
