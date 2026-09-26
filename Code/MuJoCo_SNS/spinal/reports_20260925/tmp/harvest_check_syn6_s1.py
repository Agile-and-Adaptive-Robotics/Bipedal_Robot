"""Harvest state check for curr_syn6_s1_air_deaff (stage 1, syn6 chain).
Reads the study from optuna_syn6.db and reports per-state trial counts,
every trial's objective value (trial_values join), and a sentinel-floor
census.
"""
import json
import os
import sqlite3

SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
db = os.path.join(SP, "optuna_syn6.db")
con = sqlite3.connect(db)
cur = con.cursor()
cur.execute("SELECT study_id FROM studies WHERE study_name=?",
            ("curr_syn6_s1_air_deaff",))
sid = cur.fetchone()[0]
cur.execute("SELECT state, COUNT(*) FROM trials WHERE study_id=? "
            "GROUP BY state", (sid,))
print("trial states:", cur.fetchall())
cur.execute("SELECT t.number, t.state, tv.value FROM trials t "
            "LEFT JOIN trial_values tv ON tv.trial_id = t.trial_id "
            "AND tv.objective = 0 WHERE t.study_id=? ORDER BY t.number",
            (sid,))
vals = cur.fetchall()
for n, st, v in vals:
    print(f"trial {n:3d} {st:9s} value={v}")
c_nan = c_frozen = c_gate = c_real = 0
for n, st, v in vals:
    if v is None:
        continue
    if abs(v + 200.0) < 1e-9:
        c_nan += 1
    elif abs(v + 320.0) < 1e-9:
        c_frozen += 1
    elif v <= 25.0:
        c_gate += 1
    else:
        c_real += 1
print(f"floor census: -200(NaN/unphysical/flutter)={c_nan} "
      f"-320(frozen)={c_frozen} <=25(gate/static)={c_gate} "
      f">25(real rhythm)={c_real}")
con.close()
j = json.load(open(os.path.join(SP, "curriculum_syn6_stage1.json"),
                   encoding="utf-8"))
print("winner json: score", j["score"], "trial", j["trial"],
      "params", json.dumps(j["params"]))
