# _verify_v2_deliverables_20261023.py -- independent verification of the v2
# capability-table deliverables against the three source CSVs.
# Read-only w.r.t. all inputs; writes only stdout.
import csv, json, os, re, sys
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))

def load_csv(fname, comment="#"):
    with open(os.path.join(BASE, fname), newline="", encoding="utf-8") as f:
        lines = [ln for ln in f if not ln.startswith(comment)]
    return list(csv.DictReader(lines))

surfaces = load_csv("human_surfaces_5x5_grading_20261023.csv")
hand = load_csv("hand_tuned_metrics_20261023.csv")
opt = load_csv("BiPulley_survey2_perMuscle_20261023.csv")
pm = load_csv("per_muscle_capability_20261023.csv")

problems = []   # (severity, message)
notes = []      # informational

def bad(msg): problems.append(("DISCREPANCY", msg))
def miss(msg): problems.append(("CAVEAT-MISSING", msg))
def info(msg): notes.append(msg)

f4 = lambda s: float(s)

# ---------- Check A: source row counts, labels, coverage ----------
if len(hand) != 27: bad(f"hand CSV has {len(hand)} rows, expected 27")
if len(opt) != 27: bad(f"opt CSV has {len(opt)} rows, expected 27")
if len(pm) != 27: bad(f"per-muscle deliverable has {len(pm)} rows, expected 27")
if len(surfaces) != 675:
    bad(f"surfaces CSV has {len(surfaces)} data rows, expected 27*25=675")

hand_sources = {r["source"] for r in hand}
if hand_sources != {"connor-corrected"}:
    bad(f"hand CSV coverage_source set = {hand_sources}, expected only connor-corrected")
opt_built = {r["built"] for r in opt}
opt_conv = {r["converged"] for r in opt}
opt_stop = {r["stopReason"] for r in opt}
if opt_built != {"true"}: bad(f"opt CSV built values = {opt_built}, expected all true")
if opt_conv != {"1"}: bad(f"opt CSV converged values = {opt_conv}, expected all 1")
if opt_stop != {"converged"}: bad(f"opt CSV stopReason values = {opt_stop}, expected all converged")

# actuator ids aligned 1..27; hand CSV carries no actuator_id column -> id by row order
for name, rows in (("opt", opt), ("pm", pm)):
    ids = [r["actuator_id"] for r in rows]
    if ids != [str(i) for i in range(1, 28)]:
        bad(f"{name} actuator_id sequence = {ids[:5]}..., expected 1..27 in order")
for idx, r in enumerate(hand, start=1):
    r["actuator_id"] = str(idx)  # row-order assignment; names cross-checked below
srf_ids = sorted({int(r["actuator_id"]) for r in surfaces})
if srf_ids != list(range(1, 28)):
    bad(f"surfaces actuator_id set mismatch: {srf_ids}")

# ---------- Check B: recompute the two row averages ----------
def mean(vals): return sum(vals) / len(vals)

hand_meet  = [f4(r["meetPct"]) for r in hand]
hand_act   = [f4(r["activePct"]) for r in hand]
hand_meetU = [f4(r["meetPct_unflagged"]) for r in hand]
hand_actU  = [f4(r["activePct_unflagged"]) for r in hand]
opt_meet   = [f4(r["meetPct"]) for r in opt]
opt_act    = [f4(r["activePct"]) for r in opt]
opt_meetU  = [f4(r["meetPct_unflagged"]) for r in opt]
opt_actU   = [f4(r["activePct_unflagged"]) for r in opt]

recomputed = {
    ("hand", "meet"): mean(hand_meet),   ("hand", "active"): mean(hand_act),
    ("hand", "meetU"): mean(hand_meetU), ("hand", "activeU"): mean(hand_actU),
    ("opt",  "meet"):  mean(opt_meet),   ("opt",  "active"):  mean(opt_act),
    ("opt",  "meetU"): mean(opt_meetU),  ("opt",  "activeU"): mean(opt_actU),
}
print("== Recomputed row averages (mean of printed 4-dp shares) ==")
for k, v in recomputed.items():
    print(f"  {k}: {v:.6f}  -> displayed {round(v, 2):.2f}")

# headline values claimed in md/tex/xlsx/caveats
claims = {
    ("hand", "meet"): 72.49, ("hand", "active"): 96.64,
    ("opt",  "meet"): 65.19, ("opt",  "active"): 74.27,
    ("hand", "meetU"): 73.21, ("hand", "activeU"): 95.23,
    ("opt",  "meetU"): 67.00, ("opt",  "activeU"): 72.75,
}
for k, claimed in claims.items():
    got = round(recomputed[k], 2)
    if abs(got - claimed) > 1e-9:
        bad(f"row average {k}: recomputed {got:.2f}, claimed {claimed:.2f} (raw {recomputed[k]:.6f})")
    else:
        info(f"row average {k} = {got:.2f} matches claim {claimed:.2f}")

# also verify the deliverable per-muscle CSV averages (should equal the same)
pm_hand_meet = [f4(r["hand_meet_pct"]) for r in pm]
pm_hand_act  = [f4(r["hand_active_pct"]) for r in pm]
pm_opt_meet  = [f4(r["opt_meet_pct"]) for r in pm]
pm_opt_act   = [f4(r["opt_active_pct"]) for r in pm]
for label, a, b in (
    ("pm hand meet vs hand CSV", pm_hand_meet, hand_meet),
    ("pm hand active vs hand CSV", pm_hand_act, hand_act),
    ("pm opt meet vs opt CSV", pm_opt_meet, opt_meet),
    ("pm opt active vs opt CSV", pm_opt_act, opt_act),
):
    if a != b:
        bad(f"{label}: lists differ")
        for i, (x, y) in enumerate(zip(a, b)):
            if abs(x - y) > 1e-9:
                bad(f"  act {i+1}: pm {x} vs source {y}")

# ---------- Check C: per-muscle deliverable vs all three sources, every cell ----------
hand_by_id = {r["actuator_id"]: r for r in hand}
opt_by_id = {r["actuator_id"]: r for r in opt}
pm_by_id = {r["actuator_id"]: r for r in pm}

# surfaces per-actuator counts
srf = defaultdict(lambda: {"domain": 0, "excluded": 0, "via": 0, "step": 0, "slack": 0, "short": 0, "flagged": 0})
for r in surfaces:
    a = srf[r["actuator_id"]]
    if r["in_graded_domain"] == "1":
        a["domain"] += 1
        v, s, sl, sh = (int(r["flag_via"]), int(r["flag_step"]), int(r["flag_slack"]), int(r["flag_short"]))
        a["via"] += v; a["step"] += s; a["slack"] += sl; a["short"] += sh
        if v or s or sl or sh: a["flagged"] += 1
    else:
        a["excluded"] += 1

for i in range(1, 28):
    sid = str(i)
    p, h, o = pm_by_id[sid], hand_by_id[sid], opt_by_id[sid]
    # counts vs surfaces
    for col, key in (("cells_domain", "domain"), ("cells_excluded", "excluded"),
                     ("cells_flagged_in_domain", "flagged"),
                     ("flags_via", "via"), ("flags_step", "step"),
                     ("flags_slack", "slack"), ("flags_short", "short")):
        if int(p[col]) != srf[sid][key]:
            bad(f"act {i} {col}: deliverable {p[col]} vs surfaces {srf[sid][key]}")
    # path_artifact_flag yes/no consistency
    want_pa = "yes" if srf[sid]["flagged"] > 0 else "no"
    if p["path_artifact_flag"] != want_pa:
        bad(f"act {i} path_artifact_flag: deliverable {p['path_artifact_flag']} vs surfaces-derived {want_pa}")
    # hand row vs hand CSV
    for pc, hc in (("hand_meet_pct", "meetPct"), ("hand_active_pct", "activePct"),
                   ("hand_meet_pct_unflagged", "meetPct_unflagged"), ("hand_active_pct_unflagged", "activePct_unflagged")):
        if abs(f4(p[pc]) - f4(h[hc])) > 1e-9:
            bad(f"act {i} {pc}: deliverable {p[pc]} vs hand CSV {h[hc]}")
    if p["transcription_suspect"] != h["suspect"]:
        bad(f"act {i} transcription_suspect: deliverable {p['transcription_suspect']} vs hand CSV {h['suspect']}")
    # opt row vs opt CSV
    for pc, oc in (("opt_meet_pct", "meetPct"), ("opt_active_pct", "activePct"),
                   ("opt_meet_pct_unflagged", "meetPct_unflagged"), ("opt_active_pct_unflagged", "activePct_unflagged"),
                   ("opt_nBPA", "chosen_nB"), ("opt_G", "chosen_G"),
                   ("bpa_size_mm", "bpa_size_mm"), ("opt_evals_used", "evalsUsed"),
                   ("opt_converged", "converged")):
        pv = p[pc]; ov = o[oc]
        same = (abs(f4(pv) - f4(ov)) < 1e-9) if re.fullmatch(r"-?\d+\.\d+", pv) and re.fullmatch(r"-?\d+\.\d+", ov) else (pv == ov)
        if not same:
            bad(f"act {i} {pc}: deliverable {pv!r} vs opt CSV {ov!r}")
    if p["opt_winner_mat"] != o["matFile"]:
        bad(f"act {i} opt_winner_mat: deliverable {p['opt_winner_mat']!r} vs opt CSV {o['matFile']!r}")
# dofs and muscle-name labels across all sources
for i in range(1, 28):
    sid = str(i)
    p, h, o = pm_by_id[sid], hand_by_id[sid], opt_by_id[sid]
    sp = {r["primary_dof"] for r in surfaces if r["actuator_id"] == sid}
    ss = {r["secondary_dof"] for r in surfaces if r["actuator_id"] == sid}
    if len(sp) != 1 or len(ss) != 1:
        bad(f"act {i}: surfaces dof columns not constant: {sp} {ss}")
        continue
    sp, ss = sp.pop(), ss.pop()
    if not (p["primary_dof"] == h["primary_dof"] == o["primary_dof"] == sp):
        bad(f"act {i} primary_dof mismatch: pm {p['primary_dof']} hand {h['primary_dof']} opt {o['primary_dof']} surfaces {sp}")
    if not (p["secondary_dof"] == h["secondary_dof"] == o["secondary_dof"] == ss):
        bad(f"act {i} secondary_dof mismatch: pm {p['secondary_dof']} hand {h['secondary_dof']} opt {o['secondary_dof']} surfaces {ss}")
    # muscle-name labels vs sources (act 16 duplicate-name handling)
    hm, om, pname = h["muscle"], o["master_name"], p["muscle"]
    if hm != om and i != 16:
        bad(f"act {i} muscle name: hand CSV {hm!r} vs opt CSV {om!r}")
    if pname != hm:
        bad(f"act {i} deliverable muscle label {pname!r} vs hand CSV {hm!r}")
    if i == 16:
        info(f"act 16 labels: hand {hm!r} / opt {om!r} / deliverable {pname!r} (duplicate-name bug documented in deliverable note; row CSVs disambiguate consistently)")

# ---------- Check D: totals claimed in caveats (635/40/126, hip 15/10) ----------
tot_dom = sum(v["domain"] for v in srf.values())
tot_exc = sum(v["excluded"] for v in srf.values())
tot_flag = sum(v["flagged"] for v in srf.values())
print(f"\n== Surfaces totals: graded {tot_dom}, excluded {tot_exc}, flagged {tot_flag} ==")
for got, want, what in ((tot_dom, 635, "graded cells"), (tot_exc, 40, "excluded cells"), (tot_flag, 126, "flagged graded cells")):
    if got != want: bad(f"surfaces {what}: recomputed {got}, caveats claim {want}")
    else: info(f"surfaces {what} = {got} matches caveats")

for i in (24, 25, 26, 27):
    if srf[str(i)]["domain"] != 15 or srf[str(i)]["excluded"] != 10:
        bad(f"act {i} hip surface: domain {srf[str(i)]['domain']}/excluded {srf[str(i)]['excluded']}, caveats claim 15/10")
for i in range(1, 24):
    if srf[str(i)]["excluded"] != 0:
        bad(f"act {i} non-hip surface has {srf[str(i)]['excluded']} excluded cells, expected 0")

# manifest cross-check (per-act n_graded_cells / envelope n_excluded_cells / n_flagged_graded_cells)
try:
    with open(os.path.join(BASE, "human_surfaces_manifest_20261023.json"), encoding="utf-8") as f:
        man = json.load(f)
    macts = {a["id"]: a for a in man["actuators"]}
    m_dom = sum(a["n_graded_cells"] for a in macts.values())
    m_exc = sum(a["envelope"].get("n_excluded_cells", 0) for a in macts.values())
    m_flag = sum(a["n_flagged_graded_cells"] for a in macts.values())
    print(f"manifest totals: graded={m_dom} excluded={m_exc} flagged={m_flag}")
    for mv, got, what in ((m_dom, tot_dom, "graded"), (m_exc, tot_exc, "excluded"), (m_flag, tot_flag, "flagged")):
        if mv != got:
            bad(f"manifest {what} total {mv} != recomputed-from-CSV {got}")
        else:
            info(f"manifest {what} total {mv} == recomputed-from-CSV {got} (caveats claim 'equal to the manifest' holds)")
    # per-act manifest vs surfaces CSV
    for i in range(1, 28):
        if macts[i]["n_graded_cells"] != srf[str(i)]["domain"] or \
           macts[i]["envelope"].get("n_excluded_cells", 0) != srf[str(i)]["excluded"] or \
           macts[i]["n_flagged_graded_cells"] != srf[str(i)]["flagged"]:
            bad(f"act {i}: manifest per-act counts differ from surfaces CSV counts")
    # sha1 pins quoted in caveats section 10
    for k, want in (("model_sha1", "b89062c763e6d8fa3d7651c0a4fd2820055d807c"),
                    ("map_sha1", "a97afffeeae85687810986273db6ba0151eba1a2"),
                    ("audit_csv_sha1", "47912eb6c1e858e9e62ceefe36b8a510e052c417")):
        if man[k] != want:
            bad(f"manifest {k} = {man[k]}, caveats quote {want}")
    if man["opensim_version"] != "4.4.1":
        bad(f"manifest opensim_version {man['opensim_version']!r}, caveats claim 4.4.1")
except Exception as e:
    bad(f"manifest cross-check failed to run: {e}")

# ---------- Check E: caveats section 5 (suspects) ----------
suspects_src = [(r["actuator_id"], r["muscle"], f4(r["meetPct"]), f4(r["activePct"])) for r in hand if r["suspect"] == "yes"]
print(f"\n== hand CSV suspect rows: {len(suspects_src)} ==")
if len(suspects_src) != 9:
    bad(f"hand CSV suspect count = {len(suspects_src)}, caveats claim 9")
expected_suspects = {
    "6": ("Sartorius", 36.0, 96.0), "9": ("Gracilis", 20.0, 100.0),
    "10": ("Bicep Femoris, Short Head", 100.0, 100.0), "15": ("Tibialis Posterior", 88.0, 96.0),
    "18": ("Peroneus Longus", 36.0, 96.0), "19": ("Peroneus Tertius", 80.0, 96.0),
    "22": ("Extensor Digitorum Longus", 40.0, 96.0), "25": ("Adductor Magnus", 93.3333, 93.3333),
    "27": ("Psoas X", 66.6667, 100.0),
}
for aid, name, m, a in suspects_src:
    exp = expected_suspects.get(aid)
    if exp is None:
        bad(f"suspect act {aid} ({name}) not in caveats section-5 list")
    else:
        en, em, ea = exp
        if name != en: bad(f"suspect act {aid} name: CSV {name!r} vs caveats {en!r}")
        if abs(m - em) > 0.00005: bad(f"suspect act {aid} ({name}) meet: CSV {m} vs caveats {em}")
        if abs(a - ea) > 0.00005: bad(f"suspect act {aid} ({name}) active: CSV {a} vs caveats {ea}")
for aid in expected_suspects:
    if aid not in {s[0] for s in suspects_src}:
        bad(f"caveats section-5 lists act {aid} ({expected_suspects[aid][0]}) as suspect but hand CSV says no")
# deliverable transcription_suspect column agrees
pm_sus = {r["actuator_id"] for r in pm if r["transcription_suspect"] == "yes"}
if pm_sus != set(expected_suspects.keys()):
    bad(f"deliverable suspect set {sorted(pm_sus)} != caveats suspect set {sorted(expected_suspects)}")

# ---------- Check F: caveats section 6 (budget / ceiling) ----------
evals = [int(r["evalsUsed"]) for r in opt]
srt = sorted(evals); n = len(srt)
median = (srt[n//2] if n % 2 else (srt[n//2 - 1] + srt[n//2]) / 2)
print(f"\n== opt evals: min {min(evals)} (act {opt[evals.index(min(evals))]['master_name']}), "
      f"median {median}, max {max(evals)} (act {opt[evals.index(max(evals))]['master_name']}), cap claim 300 ==")
if min(evals) != 68: bad(f"evals min {min(evals)}, caveats claim 68 (Vastus Intermedius)")
if opt[evals.index(68)]["master_name"] != "Vastus Intermedius": bad("evals-min muscle is not Vastus Intermedius")
if median != 163: bad(f"evals median {median}, caveats claim 163")
if max(evals) != 252: bad(f"evals max {max(evals)}, caveats claim 252 (External Oblique)")
if opt[evals.index(252)]["master_name"] != "External Oblique": bad("evals-max muscle is not External Oblique")
if max(evals) >= 300: bad(f"eval cap 300 was hit (max {max(evals)})")
walls = [f4(r["wallMin"]) for r in opt]
if min(walls) < 0.12 - 1e-9 or max(walls) > 0.37 + 1e-9:
    bad(f"wallMin range {min(walls)}-{max(walls)}, caveats claim 0.12-0.37")

full_meet = [(r["actuator_id"], r["master_name"]) for r in opt if f4(r["meetPct"]) == 100.0]
print(f"opt winners at 100.00% meet: {full_meet}")
exp_full = {"1": "Erector Spinae", "10": "Bicep Femoris, Short Head", "12": "Soleus", "16": None}
got_ids = {a for a, _ in full_meet}
if got_ids != {"1", "10", "12", "16"}:
    bad(f"opt 100% meet acts {sorted(got_ids)}, caveats claim [1, 10, 12, 16]")

# ---------- Check G: caveats section 7 (active zeros) ----------
act_zero = [(r["actuator_id"], r["master_name"], r["meetPct"]) for r in opt if f4(r["activePct"]) == 0.0]
print(f"\n== opt active==0 rows: {act_zero} ==")
exp_zero = {"1": "100.0000", "5": "32.0000", "13": "20.0000", "14": "20.0000", "24": "53.3333", "25": "60.0000"}
if {a for a, _, _ in act_zero} != set(exp_zero.keys()):
    bad(f"opt active-zero act set {sorted(a for a,_,_ in act_zero)} != caveats {sorted(exp_zero)}")
for a, _, m in act_zero:
    if m != exp_zero.get(a): bad(f"act {a} active-zero meet {m}, caveats say {exp_zero.get(a)}")
hand_min_active = min(f4(r["activePct"]) for r in hand)
if abs(hand_min_active - 80.0) > 1e-9:
    bad(f"hand row min active {hand_min_active}, caveats claim 80.00")
else:
    info(f"hand row min active = {hand_min_active:.2f} matches caveats claim 80.00")

# ---------- Check H: caveats section 8 (pinned dofs, constant-secondary) ----------
pinned_by_act = {}
for r in surfaces:
    pinned_by_act.setdefault(r["actuator_id"], set()).add(r["pinned_dofs_at_0"])
for i in range(1, 28):
    if len(pinned_by_act[str(i)]) != 1:
        bad(f"act {i}: pinned_dofs_at_0 not constant across surface rows: {pinned_by_act[str(i)]}")
exp_pin = {"1": "lumbar_rotation", "7": "hip_adduction_r", "20": "subtalar_angle_r",
           "21": "subtalar_angle_r", "22": "subtalar_angle_r", "23": "subtalar_angle_r",
           "24": "hip_rotation_r"}
for aid, want in exp_pin.items():
    got = next(iter(pinned_by_act.get(aid, {"<missing>"})))
    if got != want:
        bad(f"act {aid} pinned dof: surfaces {got!r}, caveats claim {want!r}")
for i in range(1, 28):
    if str(i) not in exp_pin:
        got = next(iter(pinned_by_act[str(i)]))
        if got not in ("none", "(none)", ""):
            bad(f"act {i} pinned dof {got!r} not in caveats section-8 list")

const_sec = {"10", "11", "12", "27"}
for aid in const_sec:
    o = opt_by_id[aid]
    if "robot surface constant along that axis" not in o["note"]:
        bad(f"act {aid} opt note lacks constant-secondary wording: {o['note']!r}")
for aid in map(str, range(1, 28)):
    if aid not in const_sec and "robot surface constant along that axis" in opt_by_id[aid]["note"]:
        bad(f"act {aid} opt note has constant-secondary wording but is not in caveats list")

# ---------- Check I: md header/prose claims vs sources ----------
# 620 kPa hand pressure: grep the hand grading script
ht = open(os.path.join(BASE, "ht04_grade_hand_tuned_2d_20261023.m"), encoding="utf-8", errors="replace").read()
if "620" in ht:
    m620 = [ln.strip() for ln in ht.splitlines() if "620" in ln][:4]
    info(f"620 kPa found in ht04 hand grading script: {m620[:2]}")
else:
    miss("could not find 620 kPa in ht04_grade_hand_tuned_2d_20261023.m (md claims hand row at 620 kPa)")
# act 16 hand note documents duplicate-name
if "duplicate" not in hand_by_id["16"]["note"]:
    bad("act 16 hand note does not document the duplicate-name bug")

# ---------- Check J: mat files exist ----------
missing_mats = [r["matFile"] for r in opt if not os.path.exists(r["matFile"])]
if missing_mats:
    for mm in missing_mats: bad(f"opt winner mat missing on disk: {mm}")
else:
    info("all 27 opt winner mat files exist on disk")

# ---------- Report ----------
print("\n===== PROBLEMS =====" if problems else "\n===== NO PROBLEMS =====")
for sev, msg in problems:
    print(f"[{sev}] {msg}")
print("\n===== NOTES / PASSED CHECKS =====")
for m in notes:
    print(f"[ok] {m}")
print(f"\nTotals: {len(problems)} problems, {len(notes)} notes")
