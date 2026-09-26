"""AUDIT (goal4 variant curricula, 2026-09-25) - READ-ONLY db census.

Opens optuna_w2lvar.db / optuna_syn6.db via sqlite mode=ro (no writes
possible), and optuna_walk.db for the untouched check. For each claimed
study: trial counts by state, full value list, argmax among COMPLETE
trials, winner params from trial_params; then cross-checks against the
curriculum_*_stage*.json winner files (score / trial / params).
Sentinel census: NaN<=-199.5, frozen -320, fall -150, clip -315.
"""
import json
import os
import sqlite3
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
os.chdir(SPINAL)

EXPECTED = {
    "optuna_w2lvar.db": {
        "curr_w2lvar_s1_air_deaff": dict(n=20, wtrial=17, wval=80.835587),
        "curr_w2lvar_s2_air_aff": dict(n=18, wtrial=13, wval=60.230493),
        "curr_w2lvar_s3_balance": dict(n=18, wtrial=0, wval=25.077367),
        "curr_w2lvar_s4_walk_nocontact": dict(n=23, wtrial=17, wval=-225.890865),
        "curr_w2lvar_s5_walk_contact": dict(n=22, wtrial=4, wval=-165.040506),
    },
    "optuna_syn6.db": {
        "curr_syn6_s1_air_deaff": dict(n=20, wtrial=17, wval=95.6388),
        "curr_syn6_s2_air_aff": dict(n=18, wtrial=12, wval=86.7508),
        "curr_syn6_s3_balance": dict(n=18, wtrial=11, wval=37.4630),
        "curr_syn6_s4_walk_nocontact": dict(n=23, wtrial=14, wval=-237.2534),
        "curr_syn6_s5_walk_contact": dict(n=22, wtrial=0, wval=-197.2146),
    },
}
JSON_STEM = {
    "curr_w2lvar_s1_air_deaff": "curriculum_w2lvar_stage1.json",
    "curr_w2lvar_s2_air_aff": "curriculum_w2lvar_stage2.json",
    "curr_w2lvar_s3_balance": "curriculum_w2lvar_stage3.json",
    "curr_w2lvar_s4_walk_nocontact": "curriculum_w2lvar_stage4.json",
    "curr_w2lvar_s5_walk_contact": "curriculum_w2lvar_stage5.json",
    "curr_syn6_s1_air_deaff": "curriculum_syn6_stage1.json",
    "curr_syn6_s2_air_aff": "curriculum_syn6_stage2.json",
    "curr_syn6_s3_balance": "curriculum_syn6_stage3.json",
    "curr_syn6_s4_walk_nocontact": "curriculum_syn6_stage4.json",
    "curr_syn6_s5_walk_contact": "curriculum_syn6_stage5.json",
}


def q(con, sql, args=()):
    return con.execute(sql, args).fetchall()


def census(dbfile, expected):
    print("=" * 78)
    print(f"DB (read-only): {dbfile}")
    con = sqlite3.connect(f"file:{os.path.join(SPINAL, dbfile)}?mode=ro",
                          uri=True)
    studies = q(con, "SELECT study_id, study_name FROM studies "
                     "ORDER BY study_id")
    print(f"studies ({len(studies)}):")
    for sid, sname in studies:
        print(f"  [{sid}] {sname}")
    problems = []
    for sname, exp in expected.items():
        row = q(con, "SELECT study_id FROM studies WHERE study_name=?",
                (sname,))
        if not row:
            print(f"!! {sname}: STUDY MISSING")
            problems.append(f"{sname}: missing")
            continue
        sid = row[0][0]
        rows = q(con,
                 "SELECT t.number, t.state, tv.value "
                 "FROM trials t LEFT JOIN trial_values tv "
                 "ON tv.trial_id = t.trial_id AND tv.objective = 0 "
                 "WHERE t.study_id = ? ORDER BY t.number", (sid,))
        n_total = len(rows)
        states = {}
        for _, st, _v in rows:
            states[st] = states.get(st, 0) + 1
        comp = [(num, v) for num, st, v in rows
                if st == "COMPLETE" and v is not None]
        n_comp = len(comp)
        best_num, best_val = max(comp, key=lambda r: r[1])
        vals = [v for _, v in comp]
        n_nan = sum(1 for v in vals if v <= -199.5)
        n_frozen = sum(1 for v in vals if -320.5 < v < -319.5)
        n_fall = sum(1 for v in vals if -150.5 < v < -149.5)
        n_clip = sum(1 for v in vals if -315.5 < v < -314.5)
        print(f"\n--- {sname}")
        print(f"  trials total={n_total} states={states} "
              f"COMPLETE={n_comp} (expected {exp['n']})")
        print(f"  values: {['%.4f' % v for v in vals]}")
        print(f"  sentinels: nan<={-199.5}: {n_nan}  frozen -320: "
              f"{n_frozen}  fall -150: {n_fall}  clip -315: {n_clip}")
        print(f"  argmax(COMPLETE): trial {best_num} value {best_val!r}")
        ok = True
        if n_total != exp["n"]:
            problems.append(f"{sname}: total {n_total} != claimed {exp['n']}")
            ok = False
        if n_comp != exp["n"]:
            problems.append(f"{sname}: COMPLETE {n_comp} != claimed "
                            f"{exp['n']}")
            ok = False
        if best_num != exp["wtrial"] or abs(best_val - exp["wval"]) > 5e-4:
            problems.append(f"{sname}: db best (trial {best_num}, "
                            f"{best_val!r}) != claimed (trial "
                            f"{exp['wtrial']}, {exp['wval']!r})")
            ok = False
        # winner params from db
        dbp = {name: val for name, val in
               q(con, "SELECT tp.param_name, tp.param_value "
                      "FROM trial_params tp JOIN trials t "
                      "ON t.trial_id = tp.trial_id "
                      "WHERE t.study_id = ? AND t.number = ?",
                 (sid, best_num))}
        # winner json cross-check
        jf = JSON_STEM[sname]
        jp = os.path.join(SPINAL, jf)
        if not os.path.exists(jp):
            print(f"  !! winner json {jf} MISSING")
            problems.append(f"{sname}: json missing")
            continue
        J = json.loads(open(jp, encoding="utf-8").read())
        print(f"  json {jf}: score={J.get('score')!r} trial={J.get('trial')}"
              f" study={J.get('study')} db={J.get('db')}")
        print(f"  json params ({len(J.get('params', {}))} keys): "
              f"{J.get('params')}")
        if J.get("trial") != best_num:
            problems.append(f"{sname}: json trial {J.get('trial')} != db "
                            f"best {best_num}")
        if abs(float(J.get("score", 9e9)) - best_val) > 5e-4:
            problems.append(f"{sname}: json score {J.get('score')} != db "
                            f"{best_val}")
        missing_in_json = set(dbp) - set(J.get("params", {}))
        if missing_in_json:
            problems.append(f"{sname}: json params missing db keys "
                            f"{sorted(missing_in_json)}")
        for k, v in dbp.items():
            if k in J.get("params", {}) and abs(
                    float(J["params"][k]) - float(v)) > 1e-9:
                problems.append(f"{sname}: param {k} json "
                                f"{J['params'][k]} != db {v}")
        if ok:
            print("  trial-count / winner-vs-claimed: OK")
    con.close()
    return problems


allp = []
for dbf, exp in EXPECTED.items():
    allp += census(dbf, exp)

print("\n" + "=" * 78)
print("optuna_walk.db (protected) check:")
con = sqlite3.connect(
    f"file:{os.path.join(SPINAL, 'optuna_walk.db')}?mode=ro", uri=True)
studies = q(con, "SELECT study_name FROM studies")
print(f"  {len(studies)} studies: {[s for (s,) in studies]}")
leak = [s for (s,) in studies if "w2lvar" in s or "syn6" in s]
print(f"  variant-name leak: {leak if leak else 'none'}")
con.close()

print("\n" + "=" * 78)
print("PROBLEMS FOUND:" if allp else "NO PROBLEMS: db census matches all "
      "claims (counts, winners, jsons).")
for p in allp:
    print(f"  - {p}")
